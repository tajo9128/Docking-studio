"""
Backend Connection Manager
Handles communication with Docker backend API
"""

import requests
import logging
from typing import Optional, Dict, Any, Callable
from dataclasses import dataclass
from enum import Enum

logger = logging.getLogger(__name__)

BACKEND_URL = "http://localhost:8000"


class DockingWorker:
    """Worker class for streaming docking progress via SSE"""
    
    def __init__(self, job_id: str, progress_callback: Callable[[int, str, str], None] = None):
        self.job_id = job_id
        self.progress_callback = progress_callback
        self.running = False
    
    def start_docking(self, total_ligands: int = 10) -> bool:
        """Start docking job and stream progress"""
        import uuid
        
        if not self.job_id:
            self.job_id = str(uuid.uuid4())
        
        try:
            response = requests.post(
                f"{BACKEND_URL}/dock/start",
                data={"job_id": self.job_id, "total_ligands": str(total_ligands)},
                timeout=10
            )
            if response.status_code != 200:
                logger.error(f"Failed to start docking: {response.text}")
                return False
            
            self.running = True
            self._stream_progress()
            return True
        except Exception as e:
            logger.error(f"Failed to start docking job: {e}")
            return False
    
    def _stream_progress(self):
        """Stream progress from SSE endpoint"""
        try:
            response = requests.get(
                f"{BACKEND_URL}/dock/{self.job_id}/stream",
                stream=True,
                timeout=30
            )
            response.raise_for_status()
            
            for line in response.iter_lines():
                if not self.running:
                    break
                
                if line:
                    decoded = line.decode('utf-8')
                    if decoded.startswith('data: '):
                        import json
                        try:
                            data = json.loads(decoded[6:])
                            progress = data.get('progress', 0)
                            status = data.get('status', '')
                            message = data.get('message', '')
                            
                            if self.progress_callback:
                                self.progress_callback(progress, status, message)
                            
                            if status in ['completed', 'cancelled', 'failed']:
                                break
                        except json.JSONDecodeError:
                            pass
        except requests.exceptions.RequestException as e:
            logger.error(f"SSE stream error: {e}")
            if self.progress_callback:
                self.progress_callback(0, 'error', str(e))
    
    def cancel(self):
        """Cancel the running docking job"""
        self.running = False
        try:
            requests.post(f"{BACKEND_URL}/dock/{self.job_id}/cancel", timeout=5)
        except Exception as e:
            logger.warning(f"Failed to cancel job: {e}")


class BackendStatus(Enum):
    """Backend connection status"""
    CONNECTED = "connected"
    DISCONNECTED = "disconnected"
    UNHEALTHY = "unhealthy"
    UNKNOWN = "unknown"


@dataclass
class BackendInfo:
    """Backend information"""
    status: BackendStatus
    version: Optional[str] = None
    ollama_available: bool = False
    ollama_provider: Optional[str] = None
    security_status: Optional[Dict] = None
    error_message: Optional[str] = None


class BackendConnection:
    """
    Manages connection to Docker backend
    Handles health checks, API calls, and error handling
    """
    
    def __init__(self, base_url: str = BACKEND_URL, timeout: int = 10):
        self.base_url = base_url.rstrip('/')
        self.timeout = timeout
        self._cached_info: Optional[BackendInfo] = None
    
    def is_connected(self) -> bool:
        """Check if backend is connected"""
        try:
            response = requests.get(
                f"{self.base_url}/health",
                timeout=self.timeout
            )
            return response.status_code == 200
        except requests.exceptions.ConnectionError:
            logger.debug("Backend not connected")
            return False
        except requests.exceptions.Timeout:
            logger.warning("Backend connection timeout")
            return False
        except Exception as e:
            logger.error(f"Backend connection error: {e}")
            return False
    
    def check_health(self) -> BackendInfo:
        """Get full backend health status"""
        info = BackendInfo(status=BackendStatus.UNKNOWN)
        
        try:
            # Check basic health
            response = requests.get(
                f"{self.base_url}/health",
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                info.status = BackendStatus.CONNECTED
                data = response.json()
                info.version = data.get("version", "1.0.0")
            else:
                info.status = BackendStatus.UNHEALTHY
                info.error_message = f"Status code: {response.status_code}"
                return info
            
            # Check chat/AI status
            try:
                chat_response = requests.get(
                    f"{self.base_url}/chat/status",
                    timeout=self.timeout
                )
                if chat_response.status_code == 200:
                    chat_data = chat_response.json()
                    info.ollama_available = chat_data.get("ollama_available", False)
                    info.ollama_provider = chat_data.get("provider", "offline")
            except:
                pass
            
            # Check security status
            try:
                sec_response = requests.get(
                    f"{self.base_url}/security/status",
                    timeout=self.timeout
                )
                if sec_response.status_code == 200:
                    info.security_status = sec_response.json()
            except:
                pass
            
            self._cached_info = info
            logger.info(f"Backend connected: v{info.version}, AI: {info.ollama_provider}")
            
        except requests.exceptions.ConnectionError:
            info.status = BackendStatus.DISCONNECTED
            info.error_message = "Cannot connect to backend"
            logger.warning("Backend disconnected")
        except requests.exceptions.Timeout:
            info.status = BackendStatus.UNHEALTHY
            info.error_message = "Connection timeout"
            logger.error("Backend timeout")
        except Exception as e:
            info.status = BackendStatus.UNKNOWN
            info.error_message = str(e)
            logger.error(f"Backend check error: {e}")
        
        return info
    
    def get_cached_info(self) -> Optional[BackendInfo]:
        """Get cached backend info"""
        return self._cached_info
    
    def post(self, endpoint: str, json_data: Dict = None, timeout: int = None) -> Optional[Dict]:
        """POST request to backend"""
        try:
            response = requests.post(
                f"{self.base_url}{endpoint}",
                json=json_data,
                timeout=timeout or self.timeout
            )
            response.raise_for_status()
            return response.json()
        except requests.exceptions.ConnectionError:
            logger.error(f"Connection error on POST {endpoint}")
            return None
        except requests.exceptions.Timeout:
            logger.error(f"Timeout on POST {endpoint}")
            return None
        except Exception as e:
            logger.error(f"POST error on {endpoint}: {e}")
            return None
    
    def get(self, endpoint: str, timeout: int = None) -> Optional[Dict]:
        """GET request to backend"""
        try:
            response = requests.get(
                f"{self.base_url}{endpoint}",
                timeout=timeout or self.timeout
            )
            response.raise_for_status()
            return response.json()
        except requests.exceptions.ConnectionError:
            logger.error(f"Connection error on GET {endpoint}")
            return None
        except requests.exceptions.Timeout:
            logger.error(f"Timeout on GET {endpoint}")
            return None
        except Exception as e:
            logger.error(f"GET error on {endpoint}: {e}")
            return None
    
    def send_chat(self, message: str) -> str:
        """Send chat message to AI"""
        result = self.post("/chat", json_data={"message": message}, timeout=60)
        if result:
            return result.get("response", "No response")
        return "Backend unavailable"
    
    def get_jobs(self) -> list:
        """Get all jobs"""
        result = self.get("/jobs")
        if result:
            return result.get("jobs", [])
        return []
    
    def create_job(self, job_name: str, receptor: str, ligand: str, engine: str = "vina") -> Optional[str]:
        """Create new docking job"""
        result = self.post("/jobs", json_data={
            "job_name": job_name,
            "receptor_path": receptor,
            "ligand_path": ligand,
            "engine": engine
        })
        if result:
            return result.get("job_uuid")
        return None
    
    def analyze_pose(self, receptor: str, ligand: str) -> Optional[Dict]:
        """Analyze pose interactions"""
        return self.post("/analyze", json_data={
            "receptor": receptor,
            "ligand": ligand
        })
    
    def calculate_rmsd(self, pdb1: str, pdb2: str) -> float:
        """Calculate RMSD between two poses"""
        result = self.post("/rmsd", json_data={
            "pdb1": pdb1,
            "pdb2": pdb2
        })
        if result:
            return result.get("rmsd", -1.0)
        return -1.0
    
    def run_security_scan(self) -> Optional[Dict]:
        """Run security scan"""
        return self.post("/security/scan")
    
    def get_gpu_status(self) -> Optional[Dict]:
        """Get GPU status"""
        return self.get("/gpu/status")
    
    def start_docking_job(self, job_id: str, total_ligands: int = 10) -> Optional[Dict]:
        """Start a docking job"""
        try:
            response = requests.post(
                f"{self.base_url}/dock/start",
                data={"job_id": job_id, "total_ligands": str(total_ligands)},
                timeout=10
            )
            response.raise_for_status()
            return response.json()
        except Exception as e:
            logger.error(f"Failed to start docking job: {e}")
            return None
    
    def get_docking_status(self, job_id: str) -> Optional[Dict]:
        """Get docking job status"""
        return self.get(f"/dock/{job_id}/status")


# Global connection instance
_connection: Optional[BackendConnection] = None


def get_connection() -> BackendConnection:
    """Get global backend connection"""
    global _connection
    if _connection is None:
        _connection = BackendConnection()
    return _connection


def check_backend() -> bool:
    """Quick check if backend is available"""
    return get_connection().is_connected()


if __name__ == "__main__":
    conn = BackendConnection()
    
    print("Checking backend connection...")
    
    if conn.is_connected():
        print("✓ Backend is connected")
        
        info = conn.check_health()
        print(f"  Version: {info.version}")
        print(f"  Status: {info.status.value}")
        print(f"  AI Provider: {info.ollama_provider}")
        print(f"  Ollama Available: {info.ollama_available}")
    else:
        print("✗ Backend is not connected")
        print("  Run: docker compose up -d")
