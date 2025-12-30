"""
BioDockify Docking Studio - Docker API Router
Handles Docker operations (container management, logs)
"""

from fastapi import APIRouter, HTTPException, status
from fastapi.responses import JSONResponse, StreamingResponse
from src.docker_manager import DockerManager
from src.config import Config
import logging
import json
from typing import Optional, Dict, Any

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/docker")

docker_manager = DockerManager()

@router.get("/status", response_model=dict)
async def get_docker_status():
    """Get Docker status"""
    
    is_available = docker_manager.is_docker_available()
    is_running = docker_manager.is_docker_running()
    
    status_info = {
        "available": is_available,
        "running": is_running,
        "container_status": docker_manager.get_container_status() if is_running else None
    }
    
    return JSONResponse(
        status_code=status.HTTP_200_OK,
        content=status_info
    )

@router.post("/container/start/{job_id}", response_model=dict)
async def start_container(job_id: str, receptor_file: str, ligand_file: str, parameters: dict):
    """Start Docker container for job"""
    
    # Build parameters dict with defaults
    # Note: config is not instantiated in global scope here, assuming it's available or should be imported/instantiated
    # Fixing potential NameError by instantiating Config
    config = Config()
    
    job_params = {
        "exhaustiveness": parameters.get("exhaustiveness", config.get("exhaustiveness", 8)),
        "num_modes": parameters.get("num_modes", config.get("num_modes", 9)),
        "box_size": parameters.get("box_size", config.get("box_size", 20.0))
    }
    
    success = docker_manager.start_container(
        receptor_file=receptor_file,
        ligand_file=ligand_file,
        parameters=job_params,
        job_id=job_id
    )
    
    if success:
        return JSONResponse(
            status_code=status.HTTP_200_OK,
            content={
                "job_id": job_id,
                "status": "container_started"
            }
        )
    else:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to start Docker container"
        )

@router.post("/container/stop/{container_name}", response_model=dict)
async def stop_container(container_name: Optional[str] = None):
    """Stop Docker container"""
    
    success = docker_manager.stop_container(container_name)
    
    if success:
        return JSONResponse(
            status_code=status.HTTP_200_OK,
            content={
                "status": "container_stopped"
            }
        )
    else:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to stop Docker container"
        )

@router.get("/container/logs/{container_name}", response_model=dict)
async def get_container_logs(container_name: Optional[str] = None, tail: int = 100):
    """Get Docker container logs"""
    
    logs = docker_manager.get_container_logs(container_name, tail)
    
    if logs:
        return JSONResponse(
            status_code=status.HTTP_200_OK,
            content={
                "container_name": container_name,
                "logs": logs,
                "tail": tail
            }
        )
    else:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Container logs not found"
        )

@router.post("/container/cleanup/{container_name}", response_model=dict)
async def cleanup_container(container_name: Optional[str] = None):
    """Clean up Docker container"""
    
    success = docker_manager.cleanup_container(container_name)
    
    if success:
        return JSONResponse(
            status_code=status.HTTP_200_OK,
            content={
                "status": "container_cleaned"
            }
        )
    else:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to clean up Docker container"
        )
