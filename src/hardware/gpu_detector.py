"""
Hardware Detection Module
Detects NVIDIA GPU and provides hardware information for intelligent engine selection.
"""

import subprocess
import logging
from typing import Optional, Dict, Any
from enum import Enum

logger = logging.getLogger(__name__)

VRAM_DISABLE_GPU = 4096
VRAM_REDUCE_BATCH = 6144
VRAM_ENABLE_MMGBSA = 6144


class VRAMLevel(Enum):
    """VRAM capability levels"""
    LOW = "low"
    MEDIUM = "medium" 
    HIGH = "high"


def has_gpu() -> bool:
    """
    Quick check if NVIDIA GPU is available.
    
    Returns:
        bool: True if GPU detected, False otherwise
    """
    try:
        result = subprocess.run(
            ["nvidia-smi"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=5
        )
        return result.returncode == 0
    except FileNotFoundError:
        logger.debug("nvidia-smi not found - no NVIDIA GPU")
        return False
    except subprocess.TimeoutExpired:
        logger.warning("nvidia-smi command timed out")
        return False
    except Exception as e:
        logger.warning(f"GPU detection error: {e}")
        return False


def get_gpu_info() -> Optional[Dict[str, Any]]:
    """
    Get detailed GPU information.
    
    Returns:
        dict: GPU info {name, vram_mb, driver_version, cuda_version} or None
    """
    if not has_gpu():
        return None
    
    try:
        result = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=name,memory.total,driver_version,cuda.version",
                "--format=csv,noheader,nounits"
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=10,
            text=True
        )
        
        if result.returncode != 0:
            return None
        
        lines = result.stdout.strip().split('\n')
        if not lines:
            return None
        
        parts = lines[0].split(',')
        if len(parts) < 4:
            return None
        
        return {
            "name": parts[0].strip(),
            "vram_mb": int(parts[1].strip()),
            "driver_version": parts[2].strip(),
            "cuda_version": parts[3].strip()
        }
        
    except Exception as e:
        logger.error(f"Failed to get GPU info: {e}")
        return None


def get_optimal_engine_config() -> Dict[str, Any]:
    """
    Get optimal engine configuration based on hardware.
    
    Returns:
        dict: {engine_type, reason, recommended_settings}
    """
    gpu_info = get_gpu_info()
    
    if gpu_info is None:
        return {
            "engine_type": "vina-cpu",
            "gpu_available": False,
            "reason": "No GPU detected",
            "recommended_settings": {
                "exhaustiveness": 16,
                "num_modes": 9,
                "cpu": 0  # Use all available
            }
        }
    
    vram_mb = gpu_info["vram_mb"]
    gpu_name = gpu_info["name"]
    
    if vram_mb <= 4096:
        return {
            "engine_type": "vina-gpu",
            "gpu_available": True,
            "gpu_name": gpu_name,
            "vram_mb": vram_mb,
            "reason": f"GPU detected ({gpu_name}, {vram_mb}MB) - using Vina-GPU",
            "recommended_settings": {
                "exhaustiveness": 8,
                "num_modes": 9,
                "cpu": 0
            }
        }
    
    if vram_mb <= 8192:
        return {
            "engine_type": "gnina",
            "gpu_available": True,
            "gpu_name": gpu_name,
            "vram_mb": vram_mb,
            "reason": f"GPU detected ({gpu_name}, {vram_mb}MB) - using GNINA with CNN scoring",
            "recommended_settings": {
                "exhaustiveness": 8,
                "num_modes": 9,
                "cnn_scoring": "rescore"
            }
        }
    
    return {
        "engine_type": "gnina",
        "gpu_available": True,
        "gpu_name": gpu_name,
        "vram_mb": vram_mb,
        "reason": f"Powerful GPU detected ({gpu_name}, {vram_mb}MB) - using GNINA with CNN scoring",
        "recommended_settings": {
            "exhaustiveness": 16,
            "num_modes": 9,
            "cnn_scoring": "rescore"
        }
    }


def get_vram_usage() -> Optional[int]:
    """
    Get current VRAM usage in MB.
    
    Returns:
        int: VRAM usage in MB or None if unavailable
    """
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=5,
            text=True
        )
        
        if result.returncode == 0:
            return int(result.stdout.strip())
        
    except Exception:
        pass
    
    return None


def check_cuda_available() -> bool:
    """
    Check if CUDA is available for GPU computation.
    
    Returns:
        bool: True if CUDA available, False otherwise
    """
    try:
        import torch
        return torch.cuda.is_available()
    except ImportError:
        pass
    
    try:
        result = subprocess.run(
            ["nvidia-smi"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=5
        )
        return result.returncode == 0
    except Exception:
        return False
