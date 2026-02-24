"""
Agent Zero v4 - Compute Selector Module
Safe, autonomous GPU/CPU selection based on scoring system.
No forced compute - always chooses best available.
"""

import logging
import subprocess
import psutil
from typing import Optional, Dict, Any, Tuple
from dataclasses import dataclass
from enum import Enum

logger = logging.getLogger(__name__)


class ComputeMode(Enum):
    """Compute mode selection policy"""
    AUTO = "auto"  # Autonomous adaptive selection
    GPU_FORCE = "gpu_force"  # Force GPU (not recommended)
    CPU_FORCE = "cpu_force"  # Force CPU (not recommended)


@dataclass
class ComputeScore:
    """Scoring result for a compute resource"""
    score: float
    reason: str
    available: bool
    details: Dict[str, Any]


@dataclass
class ComputeDecision:
    """Final compute decision"""
    mode: str  # "GPU" or "CPU"
    engine: str  # Engine variant to use
    confidence: float  # 0-1
    reason: str
    fallback_available: bool


class ComputeSelector:
    """
    Autonomous compute selection engine.
    Scores GPU and CPU resources and chooses best available.
    """
    
    MIN_VRAM_FOR_GPU = 1200  # MB - minimum for basic GPU docking
    HIGH_VRAM_THRESHOLD = 6000  # MB - for full features
    
    def __init__(self):
        self._gpu_info: Optional[Dict[str, Any]] = None
        self._cache_score_until = 0  # Timestamp for caching
    
    def score_gpu(self) -> ComputeScore:
        """
        Score GPU availability and health.
        
        Scoring:
        - No GPU: 0 points
        - GPU not accessible by Docker: 0 points
        - Free VRAM < 1.2GB: 10 points (limited)
        - Otherwise: 100 - utilization (higher is better)
        
        Returns:
            ComputeScore with score, reason, availability
        """
        # Check GPU exists
        if not self._has_gpu():
            return ComputeScore(
                score=0,
                reason="No NVIDIA GPU detected",
                available=False,
                details={}
            )
        
        # Check Docker GPU access
        if not self._docker_gpu_accessible():
            return ComputeScore(
                score=0,
                reason="Docker GPU access not available (nvidia-docker not configured)",
                available=False,
                details={"docker_gpu": False}
            )
        
        # Get VRAM info
        free_vram = self._get_free_vram()
        gpu_util = self._get_gpu_utilization()
        
        if free_vram is None:
            return ComputeScore(
                score=50,
                reason="GPU detected but VRAM info unavailable",
                available=True,
                details={"vram": "unknown"}
            )
        
        # Score based on VRAM
        if free_vram < self.MIN_VRAM_FOR_GPU:
            return ComputeScore(
                score=10,
                reason=f"Insufficient free VRAM ({free_vram}MB < {self.MIN_VRAM_FOR_GPU}MB)",
                available=True,
                details={"free_vram_mb": free_vram, "limited": True}
            )
        
        # Full score based on utilization
        # Lower utilization = higher score
        score = max(10, 100 - gpu_util)
        
        reason = f"GPU available ({free_vram}MB free, {gpu_util}% utilized)"
        
        return ComputeScore(
            score=score,
            reason=reason,
            available=True,
            details={
                "free_vram_mb": free_vram,
                "gpu_utilization": gpu_util,
                "limited": False
            }
        )
    
    def score_cpu(self) -> ComputeScore:
        """
        Score CPU availability and health.
        
        Scoring:
        - CPU load > 95%: 10 points (almost saturated)
        - Otherwise: 80 - load (lower is better)
        
        Returns:
            ComputeScore with score, reason, availability
        """
        try:
            cpu_load = psutil.cpu_percent(interval=0.5)
        except Exception as e:
            logger.warning(f"Failed to get CPU load: {e}")
            cpu_load = 50  # Assume moderate load
        
        cpu_count = psutil.cpu_count()
        
        if cpu_load > 95:
            return ComputeScore(
                score=10,
                reason=f"CPU heavily loaded ({cpu_load}%)",
                available=True,
                details={"cpu_load": cpu_load, "cpu_count": cpu_count}
            )
        
        # Score: 80 - load (so lower load = higher score)
        score = max(10, 80 - cpu_load)
        
        return ComputeScore(
            score=score,
            reason=f"CPU available ({cpu_load}% load, {cpu_count} cores)",
            available=True,
            details={
                "cpu_load": cpu_load,
                "cpu_count": cpu_count
            }
        )
    
    def choose_best(self, mode: ComputeMode = ComputeMode.AUTO) -> ComputeDecision:
        """
        Choose best compute resource based on scoring.
        
        Args:
            mode: Selection policy (AUTO, GPU_FORCE, CPU_FORCE)
            
        Returns:
            ComputeDecision with mode, engine, confidence, reason
        """
        gpu_score = self.score_gpu()
        cpu_score = self.score_cpu()
        
        logger.info(f"Compute scoring - GPU: {gpu_score.score:.1f}, CPU: {cpu_score.score:.1f}")
        
        # Handle forced modes (for advanced override)
        if mode == ComputeMode.GPU_FORCE:
            if gpu_score.available:
                return ComputeDecision(
                    mode="GPU",
                    engine="gnina-gpu",
                    confidence=0.9,
                    reason=f"GPU forced: {gpu_score.reason}",
                    fallback_available=True
                )
            else:
                logger.warning("GPU forced but unavailable, falling back to CPU")
                return ComputeDecision(
                    mode="CPU",
                    engine="vina-cpu",
                    confidence=0.5,
                    reason=f"GPU forced but unavailable: {gpu_score.reason}",
                    fallback_available=False
                )
        
        if mode == ComputeMode.CPU_FORCE:
            return ComputeDecision(
                mode="CPU",
                engine="vina-cpu",
                confidence=0.95,
                reason=f"CPU forced: {cpu_score.reason}",
                fallback_available=gpu_score.available
            )
        
        # AUTO mode: Compare scores
        if gpu_score.score > cpu_score.score:
            # GPU is better
            confidence = min(0.95, gpu_score.score / 100 + 0.5)
            return ComputeDecision(
                mode="GPU",
                engine="gnina-gpu",
                confidence=confidence,
                reason=f"GPU scores higher ({gpu_score.score:.1f} vs {cpu_score.score:.1f}): {gpu_score.reason}",
                fallback_available=True
            )
        elif cpu_score.score > gpu_score.score:
            # CPU is better
            confidence = min(0.95, cpu_score.score / 100 + 0.5)
            return ComputeDecision(
                mode="CPU",
                engine="vina-cpu",
                confidence=confidence,
                reason=f"CPU scores higher ({cpu_score.score:.1f} vs {gpu_score.score:.1f}): {cpu_score.reason}",
                fallback_available=True
            )
        else:
            # Equal scores - prefer CPU (more stable)
            return ComputeDecision(
                mode="CPU",
                engine="vina-cpu",
                confidence=0.6,
                reason=f"Equal scores ({gpu_score.score:.1f}), defaulting to CPU",
                fallback_available=gpu_score.available
            )
    
    def _has_gpu(self) -> bool:
        """Check if NVIDIA GPU exists"""
        try:
            result = subprocess.run(
                ["nvidia-smi"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=5
            )
            return result.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False
    
    def _docker_gpu_accessible(self) -> bool:
        """Check if Docker can access GPU"""
        try:
            result = subprocess.run(
                ["docker", "run", "--rm", "--gpus", "all", "nvidia/cuda:12.0-base", "nvidia-smi"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=30
            )
            return result.returncode == 0
        except Exception:
            return False
    
    def _get_free_vram(self) -> Optional[int]:
        """Get free VRAM in MB"""
        try:
            result = subprocess.run(
                ["nvidia-smi", "--query-gpu=memory.free", "--format=csv,noheader,nounits"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=5,
                text=True
            )
            if result.returncode == 0:
                return int(result.stdout.strip().split('\n')[0])
        except Exception:
            pass
        return None
    
    def _get_gpu_utilization(self) -> int:
        """Get GPU utilization percentage (0-100)"""
        try:
            result = subprocess.run(
                ["nvidia-smi", "--query-gpu=utilization.gpu", "--format=csv,noheader,nounits"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=5,
                text=True
            )
            if result.returncode == 0:
                return int(result.stdout.strip().split('\n')[0])
        except Exception:
            pass
        return 0
    
    def get_status_summary(self) -> Dict[str, Any]:
        """Get comprehensive status summary for UI display"""
        gpu = self.score_gpu()
        cpu = self.score_cpu()
        decision = self.choose_best()
        
        return {
            "gpu": {
                "available": gpu.available,
                "score": gpu.score,
                "reason": gpu.reason,
                "details": gpu.details
            },
            "cpu": {
                "available": cpu.available,
                "score": cpu.score,
                "reason": cpu.reason,
                "details": cpu.details
            },
            "selected": {
                "mode": decision.mode,
                "engine": decision.engine,
                "confidence": decision.confidence,
                "reason": decision.reason
            }
        }


# Singleton instance
_selector: Optional[ComputeSelector] = None

def get_compute_selector() -> ComputeSelector:
    """Get singleton compute selector instance"""
    global _selector
    if _selector is None:
        _selector = ComputeSelector()
    return _selector


def auto_select_compute() -> ComputeDecision:
    """Convenience function for auto-compute selection"""
    return get_compute_selector().choose_best(ComputeMode.AUTO)
