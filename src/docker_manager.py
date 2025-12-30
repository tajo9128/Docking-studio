"""
BioDockify Docking Studio - Docker Manager
Handles Docker container operations and lifecycle
"""

import docker
import logging
from typing import Optional, Dict, Any, List, Tuple
import os

from src.utils.log_utils import get_logger, log_exception

logger = get_logger(__name__)

class DockerManager:
    """Manages Docker containers for BioDockify"""
    
    def __init__(self):
        self.client = None
        self.container = None
        self.container_name = "biodockify_worker"
        self.image_name = "biodockify/biodockify:latest"
        self._connect()
        
    def _connect(self):
        """Connect to Docker Desktop"""
        try:
            self.client = docker.from_env()
            logger.info("Connected to Docker Desktop")
        except Exception as e:
            logger.error(f"Failed to connect to Docker: {e}")
            self.client = None

    def is_docker_available(self) -> bool:
        """Check if Docker is available"""
        return self.client is not None

    def is_docker_running(self) -> bool:
        """Check if Docker daemon is running"""
        if not self.client:
            return False
        try:
            self.client.ping()
            return True
        except Exception:
            return False

    def start_container(self, receptor_file: str, ligand_file: str, 
                       parameters: Dict[str, Any], job_id: str) -> bool:
        """Start docking container"""
        if not self.is_docker_running():
            logger.error("Docker is not running")
            return False
            
        try:
            # Cleanup existing container if any
            self.cleanup_container()
            
            # Prepare volumes
            # Assuming files are absolute paths
            receptor_dir = os.path.dirname(os.path.abspath(receptor_file))
            ligand_dir = os.path.dirname(os.path.abspath(ligand_file))
            output_dir = os.path.abspath(os.path.join(os.getcwd(), "data"))
            
            volumes = {
                receptor_dir: {'bind': '/data/receptor', 'mode': 'ro'},
                ligand_dir: {'bind': '/data/ligand', 'mode': 'ro'},
                output_dir: {'bind': '/data/output', 'mode': 'rw'}
            }
            
            # Prepare environment variables
            env = {
                "JOB_ID": job_id,
                # Add other params...
            }
            env.update(parameters)
            
            # Pull image if missing (Optimization #10)
            try:
                self.client.images.get(self.image_name)
            except docker.errors.ImageNotFound:
                logger.info(f"Image {self.image_name} not found locally. Pulling...")
                self.client.images.pull(self.image_name)
            
            # Use unique container name for this job (Refinement #1)
            # This prevents collisions with zombie containers from previous runs
            unique_container_name = f"biodockify-job-{job_id}"
            
            # Run container
            logger.info(f"Starting container {unique_container_name} for job {job_id}")
            self.container = self.client.containers.run(
                self.image_name,
                name=unique_container_name,
                detach=True,
                volumes=volumes,
                environment=env,
                # auto_remove=False # We handle removal manually now
            )
            
            # Update current container name so stop_container knows what to kill
            self.container_name = unique_container_name
            
            return True
            
        except Exception as e:
            log_exception(logger, e, "Failed to start container")
            return False

    def stop_container(self, container_name: Optional[str] = None) -> bool:
        """Stop and remove Docker container"""
        target_name = container_name or self.container_name
        logger.info(f"Stopping container: {target_name}")
        
        if not self.client:
            return False
            
        try:
            # Find and stop/remove container (Fix #1)
            found = False
            for container in self.client.containers.list(all=True):
                if container.name == target_name:
                    if container.status == 'running':
                        container.stop(timeout=10)
                    container.remove(force=True)
                    logger.info(f"Container {target_name} stopped and removed.")
                    found = True
                    break
            
            if not found:
                logger.warning(f"Container {target_name} not found to stop.")
            
            self.container = None
            return True
            
        except Exception as e:
            log_exception(logger, e, f"Failed to stop container {target_name}")
            return False

    def get_container_logs(self, container_name: Optional[str] = None, tail: int = 100) -> str:
        """Get logs from container"""
        target_name = container_name or self.container_name
        
        if not self.client:
            return ""
            
        try:
            container = self.client.containers.get(target_name)
            logs = container.logs(tail=tail).decode('utf-8')
            return logs
        except Exception as e:
            logger.error(f"Failed to get logs for {target_name}: {e}")
            return ""

    def get_container_status(self, container_name: Optional[str] = None) -> str:
        """Get container status"""
        target_name = container_name or self.container_name
        
        if not self.client:
            return "unknown"
            
        try:
            container = self.client.containers.get(target_name)
            return container.status
        except Exception:
            return "not_found"

    def cleanup_container(self, container_name: Optional[str] = None) -> bool:
        """Example cleanup alias"""
        return self.stop_container(container_name)
