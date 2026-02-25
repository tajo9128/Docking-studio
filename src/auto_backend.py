"""
Auto Backend Starter
Automatically starts Docker backend when needed
"""

import subprocess
import logging
import time
import sys

logger = logging.getLogger(__name__)


class AutoBackendStarter:
    """
    Automatically starts Docker backend if not running
    """
    
    def __init__(self, docker_compose_dir: str = "."):
        self.docker_compose_dir = docker_compose_dir
    
    def is_backend_running(self) -> bool:
        """Check if backend is running"""
        try:
            import requests
            r = requests.get("http://localhost:8000/health", timeout=2)
            return r.status_code == 200
        except:
            return False
    
    def start_backend(self) -> bool:
        """Start Docker backend"""
        logger.info("Starting Docker backend...")
        
        try:
            # Run docker compose up -d
            result = subprocess.run(
                ["docker", "compose", "up", "-d"],
                cwd=self.docker_compose_dir,
                capture_output=True,
                text=True,
                timeout=120
            )
            
            if result.returncode != 0:
                logger.error(f"Failed to start backend: {result.stderr}")
                return False
            
            # Wait for backend to be ready
            logger.info("Waiting for backend to be ready...")
            for i in range(30):
                time.sleep(2)
                if self.is_backend_running():
                    logger.info("Backend is ready!")
                    return True
            
            logger.error("Backend did not become ready in time")
            return False
            
        except subprocess.TimeoutExpired:
            logger.error("Backend start timed out")
            return False
        except Exception as e:
            logger.error(f"Error starting backend: {e}")
            return False
    
    def ensure_backend(self) -> bool:
        """Ensure backend is running, start if not"""
        if self.is_backend_running():
            logger.info("Backend is already running")
            return True
        
        logger.warning("Backend not running, attempting to start...")
        return self.start_backend()


def auto_start_backend() -> bool:
    """
    Convenience function to auto-start backend
    Returns True if backend is running
    """
    starter = AutoBackendStarter()
    return starter.ensure_backend()


if __name__ == "__main__":
    print("Checking backend...")
    
    starter = AutoBackendStarter()
    
    if starter.is_backend_running():
        print("✓ Backend is running")
    else:
        print("✗ Backend not running")
        print("Starting backend...")
        
        if starter.start_backend():
            print("✓ Backend started successfully")
        else:
            print("✗ Failed to start backend")
            sys.exit(1)
