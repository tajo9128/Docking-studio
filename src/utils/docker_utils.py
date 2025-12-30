"""
BioDockify Docking Studio - Docker Utilities
Helper functions for Docker operations
"""

import subprocess
import logging
import shutil
from typing import Tuple, Optional, Dict

logger = logging.getLogger(__name__)

def check_docker_availability() -> bool:
    """
    Check if Docker is available and running.
    Uses direct client ping for speed with fallback to subprocess.
    """
    # 1. Try using docker client first (faster/more reliable)
    try:
        import docker
        client = docker.from_env()
        client.ping()
        logger.info("Docker is available (verified via client ping)")
        return True
    except ImportError:
        logger.debug("docker python package not installed")
    except Exception as e:
        logger.debug(f"Docker client ping failed: {e}")
        
    # 2. Fallback to subprocess with timeout
    try:
        # Popen is non-blocking until communicate, so we can timeout
        proc = subprocess.Popen(
            ["docker", "info"], 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE
        )
        try:
            outs, errs = proc.communicate(timeout=2)
            if proc.returncode == 0:
                logger.info("Docker is available (verified via subprocess)")
                return True
            else:
                # Only log detailed warning if return code is non-zero (meaning command ran but failed)
                # If command simply doesn't exist, FileNotFoundError is caught below
                if outs:
                    logger.debug(f"Docker info stdout: {outs.decode().strip()}")
                if errs:
                    logger.warning(f"Docker check failed: {errs.decode().strip()}")
        except subprocess.TimeoutExpired:
            proc.kill()
            logger.warning("Docker check timed out")
            
        return False
        
    except FileNotFoundError:
        logger.error("Docker executable not found")
        return False
    except Exception as e:
        logger.error(f"Failed to check Docker availability: {e}")
        return False

def get_docker_info() -> Tuple[bool, Optional[str], Optional[str]]:
    """
    Get detailed Docker information.
    Returns: (is_available, error_message, version_info)
    """
    if not shutil.which("docker"):
        return False, "Docker executable not found on PATH.", None
        
    if not check_docker_availability():
        return False, "Docker daemon is not running.", None
        
    try:
        result = subprocess.run(
            ["docker", "--version"], 
            capture_output=True, 
            text=True, 
            timeout=5
        )
        return True, None, result.stdout.strip()
    except Exception as e:
        return True, None, "Unknown version"

def format_docker_command(container_name: str, command: str) -> list:
    """Format command for execution in container"""
    return ["docker", "exec", container_name, "sh", "-c", command]
