"""
BioDockify Docking Studio - Dependencies API Router
Handles system dependency checks
"""

from fastapi import APIRouter, HTTPException, status
from fastapi.responses import JSONResponse
from src.utils.docker_utils import check_docker_availability as check_docker, get_docker_info  # Use alias
from src.config import Config
import logging

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/dependencies")

@router.get("/check", response_model=dict)
async def check_dependencies():
    """Check all system dependencies"""
    
    deps_status = {
        "docker_available": False,
        "docker_running": False,
        "docker_info": None,
        "status": "unknown"
    }
    
    # Check Docker
    docker_available, docker_message = check_docker()
    
    if docker_available:
        deps_status["docker_available"] = True
        docker_info = get_docker_info()
        
        if docker_info:
            deps_status["docker_info"] = docker_info
            deps_status["docker_running"] = True
            deps_status["status"] = "ready"
        else:
            deps_status["status"] = "error"
    else:
        deps_status["status"] = "unavailable"
        deps_status["message"] = docker_message
    
    return JSONResponse(
        status_code=status.HTTP_200_OK,
        content=deps_status
    )

@router.get("/docker-info", response_model=dict)
async def get_docker_information():
    """Get detailed Docker information"""
    
    docker_info = get_docker_info()
    
    if docker_info:
        return JSONResponse(
            status_code=status.HTTP_200_OK,
            content=docker_info
        )
    else:
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Docker Desktop is not running"
        )
