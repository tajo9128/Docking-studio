"""
Docker Monitor - Real-time Container Monitoring
Uses Docker API for container status and log monitoring
"""

import os
import time
import threading
from typing import Dict, List, Optional, Callable, Any
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
import logging

logger = logging.getLogger(__name__)

DOCKER_AVAILABLE = False
try:
    import docker
    from docker.errors import NotFound, APIError
    DOCKER_AVAILABLE = True
except ImportError:
    logger.warning("Docker SDK not available. Install with: pip install docker")


class ContainerStatus(Enum):
    """Container status states"""
    RUNNING = "running"
    EXITED = "exited"
    RESTARTING = "restarting"
    DEAD = "dead"
    CREATED = "created"
    PAUSED = "paused"


@dataclass
class ContainerInfo:
    """Container information"""
    container_id: str
    name: str
    image: str
    status: ContainerStatus
    created: str
    started_at: Optional[str] = None
    finished_at: Optional[str] = None
    exit_code: Optional[int] = None
    gpu_used: bool = False


@dataclass
class ResourceUsage:
    """Container resource usage"""
    cpu_percent: float = 0.0
    memory_usage_mb: float = 0.0
    memory_limit_mb: float = 0.0
    memory_percent: float = 0.0
    gpu_memory_mb: float = 0.0
    gpu_percent: float = 0.0


@dataclass
class LogEntry:
    """A log entry"""
    timestamp: str
    message: str
    stream: str  # stdout/stderr


class DockerMonitor:
    """
    Docker container monitor using Docker API.
    Provides real-time status, resource usage, and log monitoring.
    """
    
    def __init__(self):
        """Initialize Docker monitor"""
        self.client = None
        self._connect()
        
        self.monitors: Dict[str, Any] = {}
        self.log_handlers: Dict[str, Callable] = {}
        
        logger.info("DockerMonitor initialized")
    
    def _connect(self):
        """Connect to Docker daemon"""
        if not DOCKER_AVAILABLE:
            logger.warning("Docker SDK not available")
            return
        
        try:
            # Try default socket
            self.client = docker.from_env()
            # Test connection
            self.client.ping()
            logger.info("Connected to Docker daemon")
        except Exception as e:
            logger.warning(f"Could not connect to Docker: {e}")
            self.client = None
    
    def is_available(self) -> bool:
        """Check if Docker is available"""
        return self.client is not None
    
    def list_containers(self, all_containers: bool = True) -> List[ContainerInfo]:
        """List all containers"""
        if not self.client:
            return []
        
        try:
            containers = self.client.containers.list(all=all_containers)
            
            result = []
            for c in containers:
                info = ContainerInfo(
                    container_id=c.id[:12],
                    name=c.name,
                    image=c.image.tags[0] if c.image.tags else c.image.short_id,
                    status=ContainerStatus(c.status),
                    created=c.attrs.get('Created', ''),
                    started_at=c.attrs.get('State', {}).get('StartedAt'),
                    finished_at=c.attrs.get('State', {}).get('FinishedAt'),
                    exit_code=c.attrs.get('State', {}).get('ExitCode'),
                )
                result.append(info)
            
            return result
            
        except Exception as e:
            logger.error(f"Failed to list containers: {e}")
            return []
    
    def get_container(self, container_id: str) -> Optional[ContainerInfo]:
        """Get container info by ID or name"""
        if not self.client:
            return None
        
        try:
            c = self.client.containers.get(container_id)
            
            return ContainerInfo(
                container_id=c.id[:12],
                name=c.name,
                image=c.image.tags[0] if c.image.tags else c.image.short_id,
                status=ContainerStatus(c.status),
                created=c.attrs.get('Created', ''),
                started_at=c.attrs.get('State', {}).get('StartedAt'),
                finished_at=c.attrs.get('State', {}).get('FinishedAt'),
                exit_code=c.attrs.get('State', {}).get('ExitCode'),
            )
            
        except NotFound:
            logger.warning(f"Container not found: {container_id}")
            return None
        except Exception as e:
            logger.error(f"Failed to get container: {e}")
            return None
    
    def get_container_status(self, container_id: str) -> Optional[ContainerStatus]:
        """Get container status"""
        info = self.get_container(container_id)
        return info.status if info else None
    
    def get_resource_usage(self, container_id: str) -> Optional[ResourceUsage]:
        """Get container resource usage"""
        if not self.client:
            return None
        
        try:
            c = self.client.containers.get(container_id)
            
            # Get stats
            stats = c.stats(stream=False)
            
            # CPU calculation
            cpu_delta = stats['cpu_stats']['cpu_usage']['total_usage'] - \
                       stats['precpu_stats']['cpu_usage']['total_usage']
            system_delta = stats['cpu_stats']['system_cpu_usage'] - \
                          stats['precpu_stats']['system_cpu_usage']
            
            cpu_percent = 0.0
            if system_delta > 0:
                cpu_percent = (cpu_delta / system_delta) * \
                            len(stats['cpu_stats']['cpu_usage'].get('percpu_usage', [0])) * 100
            
            # Memory
            memory_usage = stats['memory_stats'].get('usage', 0)
            memory_limit = stats['memory_stats'].get('limit', 1)
            
            memory_usage_mb = memory_usage / (1024 * 1024)
            memory_limit_mb = memory_limit / (1024 * 1024)
            memory_percent = (memory_usage / memory_limit) * 100
            
            # GPU (if available via nvidia-smi)
            gpu_memory_mb = 0
            gpu_percent = 0
            try:
                import subprocess
                result = subprocess.run(
                    ['nvidia-smi', '--query-compute-apps=used_memory', 
                     '--format=csv,noheader,nounits', '-c', '1'],
                    capture_output=True, text=True, timeout=2
                )
                if result.returncode == 0:
                    gpu_memory_mb = float(result.stdout.strip())
            except:
                pass
            
            return ResourceUsage(
                cpu_percent=cpu_percent,
                memory_usage_mb=memory_usage_mb,
                memory_limit_mb=memory_limit_mb,
                memory_percent=memory_percent,
                gpu_memory_mb=gpu_memory_mb,
                gpu_percent=gpu_percent,
            )
            
        except Exception as e:
            logger.error(f"Failed to get resource usage: {e}")
            return None
    
    def start_log_monitor(
        self,
        container_id: str,
        callback: Callable[[LogEntry], None],
        follow: bool = True,
        tail: int = 100
    ):
        """
        Start monitoring container logs.
        
        Args:
            container_id: Container ID or name
            callback: Function to call with each log entry
            follow: Follow log stream
            tail: Number of lines to read from end
        """
        if not self.client:
            return
        
        def monitor_logs():
            try:
                c = self.client.containers.get(container_id)
                
                for line in c.logs(stream=follow, tail=tail):
                    if line:
                        # Parse log line
                        text = line.decode('utf-8', errors='ignore').rstrip()
                        
                        # Determine stream type
                        stream = 'stdout'
                        if hasattr(line, 'stream'):
                            stream = line.stream
                        
                        entry = LogEntry(
                            timestamp=datetime.now().isoformat(),
                            message=text,
                            stream=stream,
                        )
                        
                        callback(entry)
                        
            except Exception as e:
                logger.error(f"Log monitor error: {e}")
        
        # Start in background thread
        thread = threading.Thread(target=monitor_logs, daemon=True)
        thread.start()
        
        self.log_handlers[container_id] = thread
    
    def stop_log_monitor(self, container_id: str):
        """Stop monitoring logs for a container"""
        if container_id in self.monitors:
            # Thread will stop on its own when container dies
            del self.monitors[container_id]
    
    def get_container_logs(
        self,
        container_id: str,
        tail: int = 100,
        timestamps: bool = True
    ) -> List[LogEntry]:
        """Get container logs (non-streaming)"""
        if not self.client:
            return []
        
        try:
            c = self.client.containers.get(container_id)
            logs = c.logs(tail=tail, timestamps=timestamps).decode('utf-8', errors='ignore')
            
            entries = []
            for line in logs.split('\n'):
                if line.strip():
                    entries.append(LogEntry(
                        timestamp=datetime.now().isoformat(),
                        message=line,
                        stream='stdout',
                    ))
            
            return entries
            
        except Exception as e:
            logger.error(f"Failed to get logs: {e}")
            return []
    
    def restart_container(self, container_id: str) -> bool:
        """Restart a container"""
        if not self.client:
            return False
        
        try:
            c = self.client.containers.get(container_id)
            c.restart()
            logger.info(f"Container restarted: {container_id}")
            return True
        except Exception as e:
            logger.error(f"Failed to restart container: {e}")
            return False
    
    def stop_container(self, container_id: str, timeout: int = 10) -> bool:
        """Stop a container"""
        if not self.client:
            return False
        
        try:
            c = self.client.containers.get(container_id)
            c.stop(timeout=timeout)
            logger.info(f"Container stopped: {container_id}")
            return True
        except Exception as e:
            logger.error(f"Failed to stop container: {e}")
            return False
    
    def remove_container(self, container_id: str, force: bool = False) -> bool:
        """Remove a container"""
        if not self.client:
            return False
        
        try:
            c = self.client.containers.get(container_id)
            c.remove(force=force)
            logger.info(f"Container removed: {container_id}")
            return True
        except Exception as e:
            logger.error(f"Failed to remove container: {e}")
            return False


# Global monitor instance
_docker_monitor: Optional[DockerMonitor] = None


def get_docker_monitor() -> DockerMonitor:
    """Get or create Docker monitor instance"""
    global _docker_monitor
    if _docker_monitor is None:
        _docker_monitor = DockerMonitor()
    return _docker_monitor
