"""
Resource Guard - VRAM Monitoring and GPU Throttling
Protects 4GB GPUs from OOM crashes and optimizes performance.
"""

import subprocess
import logging
from typing import Dict, Any, Optional, Tuple
from dataclasses import dataclass
from enum import Enum

logger = logging.getLogger(__name__)


class VRAMLevel(Enum):
    """VRAM capability levels"""
    LOW = "low"           # < 4GB - CPU only
    MEDIUM = "medium"    # 4-6GB - Limited GPU
    HIGH = "high"        # 6GB+ - Full GPU


@dataclass
class GPUCapabilities:
    """GPU capabilities based on VRAM"""
    vram_level: VRAMLevel
    total_vram_mb: int
    free_vram_mb: int
    supports_vina_gpu: bool
    supports_gnina: bool
    supports_mmgbsa: bool
    recommended_cnn_batch: int
    recommended_exhaustiveness: int
    max_parallel_jobs: int


VRAM_DISABLE_GPU = 4096
VRAM_REDUCE_BATCH = 6144
VRAM_ENABLE_MMGBSA = 6144


class ResourceGuard:
    """
    Monitors GPU resources and determines safe operational parameters.
    Essential for GTX 1650 (4GB) and similar GPUs.
    """
    
    def __init__(self):
        self._cached_info: Optional[Dict[str, Any]] = None
        self._cache_time: float = 0
    
    def check_gpu_available(self) -> bool:
        """Quick check if GPU is available"""
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
    
    def get_vram_info(self) -> Tuple[int, int, int]:
        """
        Get VRAM information.
        
        Returns:
            tuple: (total_mb, free_mb, used_mb)
        """
        try:
            result = subprocess.run(
                ["nvidia-smi", "--query-gpu=memory.total,memory.free,memory.used",
                 "--format=csv,noheader,nounits"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=5,
                text=True
            )
            
            if result.returncode == 0:
                lines = result.stdout.strip().split('\n')
                if lines:
                    parts = lines[0].split(',')
                    if len(parts) >= 3:
                        total = int(parts[0].strip())
                        free = int(parts[1].strip())
                        used = int(parts[2].strip())
                        return total, free, used
            
        except Exception as e:
            logger.warning(f"Failed to get VRAM info: {e}")
        
        return 0, 0, 0
    
    def get_gpu_capabilities(self) -> GPUCapabilities:
        """
        Determine GPU capabilities based on VRAM.
        
        Returns:
            GPUCapabilities: Safe operational parameters
        """
        if not self.check_gpu_available():
            return GPUCapabilities(
                vram_level=VRAMLevel.LOW,
                total_vram_mb=0,
                free_vram_mb=0,
                supports_vina_gpu=False,
                supports_gnina=False,
                supports_mmgbsa=False,
                recommended_cnn_batch=0,
                recommended_exhaustiveness=8,
                max_parallel_jobs=0
            )
        
        total, free, _ = self.get_vram_info()
        
        if total < VRAM_DISABLE_GPU:
            return GPUCapabilities(
                vram_level=VRAMLevel.LOW,
                total_vram_mb=total,
                free_vram_mb=free,
                supports_vina_gpu=True,
                supports_gnina=False,
                supports_mmgbsa=False,
                recommended_cnn_batch=4,
                recommended_exhaustiveness=4,
                max_parallel_jobs=1
            )
        
        if total < VRAM_REDUCE_BATCH:
            return GPUCapabilities(
                vram_level=VRAMLevel.MEDIUM,
                total_vram_mb=total,
                free_vram_mb=free,
                supports_vina_gpu=True,
                supports_gnina=True,
                supports_mmgbsa=False,
                recommended_cnn_batch=8,
                recommended_exhaustiveness=8,
                max_parallel_jobs=1
            )
        
        return GPUCapabilities(
            vram_level=VRAMLevel.HIGH,
            total_vram_mb=total,
            free_vram_mb=free,
            supports_vina_gpu=True,
            supports_gnina=True,
            supports_mmgbsa=True,
            recommended_cnn_batch=16,
            recommended_exhaustiveness=16,
            max_parallel_jobs=2 if total >= 12000 else 1
        )
    
    def can_run_gpu_job(self, required_vram_mb: int = 1500) -> Tuple[bool, str]:
        """
        Check if GPU job can run safely.
        
        Args:
            required_vram_mb: Minimum VRAM required (default 1.5GB)
            
        Returns:
            tuple: (can_run, reason)
        """
        if not self.check_gpu_available():
            return False, "No GPU detected"
        
        _, free, _ = self.get_vram_info()
        
        if free < required_vram_mb:
            return False, f"Insufficient VRAM: {free}MB free, need {required_vram_mb}MB"
        
        return True, f"OK: {free}MB VRAM available"
    
    def get_optimal_parameters(self, job_type: str) -> Dict[str, Any]:
        """
        Get optimal parameters for job type based on GPU.
        
        Args:
            job_type: Type of job (vina, gnina, mmgbsa)
            
        Returns:
            dict: Optimal parameters
        """
        caps = self.get_gpu_capabilities()
        
        params = {
            "use_gpu": caps.supports_vina_gpu,
            "exhaustiveness": caps.recommended_exhaustiveness,
            "max_parallel": caps.max_parallel_jobs
        }
        
        if job_type == "vina":
            params["cpu"] = 0 if caps.supports_vina_gpu else 8
        
        elif job_type == "gnina":
            params["cnn_batch_size"] = caps.recommended_cnn_batch if caps.supports_gnina else 0
            params["cnn_scoring"] = "rescore" if caps.supports_gnina else "none"
        
        elif job_type == "mmgbsa":
            params["enabled"] = caps.supports_mmgbsa
            params["top_poses"] = 5 if caps.supports_mmgbsa else 0
        
        return params
    
    def wait_for_vram(self, required_mb: int = 1500, max_wait_seconds: int = 60) -> bool:
        """
        Wait for sufficient VRAM to become available.
        
        Args:
            required_mb: Required free VRAM in MB
            max_wait_seconds: Maximum time to wait
            
        Returns:
            bool: True if VRAM available, False if timeout
        """
        import time
        
        start = time.time()
        
        while time.time() - start < max_wait_seconds:
            _, free, _ = self.get_vram_info()
            
            if free >= required_mb:
                logger.info(f"VRAM available: {free}MB")
                return True
            
            logger.info(f"Waiting for VRAM: {free}MB / {required_mb}MB")
            time.sleep(5)
        
        logger.warning(f"Timeout waiting for VRAM: needed {required_mb}MB")
        return False
    
    def get_status_summary(self) -> Dict[str, Any]:
        """Get GPU status summary for UI display"""
        caps = self.get_gpu_capabilities()
        
        return {
            "gpu_available": caps.supports_vina_gpu,
            "vram_level": caps.vram_level.value,
            "total_vram_mb": caps.total_vram_mb,
            "free_vram_mb": caps.free_vram_mb,
            "supports_gnina": caps.supports_gnina,
            "supports_mmgbsa": caps.supports_mmgbsa,
            "recommended_settings": {
                "cnn_batch": caps.recommended_cnn_batch,
                "exhaustiveness": caps.recommended_exhaustiveness
            }
        }
