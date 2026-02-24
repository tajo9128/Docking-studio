"""
Retry Manager - Auto-Retry and Self-Healing
Implements automatic retry logic with parameter adjustment
"""

import time
from typing import Dict, Optional, Callable, List, Tuple, Any
from dataclasses import dataclass, field
from enum import Enum
import logging

logger = logging.getLogger(__name__)

from .failure_detector import ErrorType, Failure


class RetryStrategy(Enum):
    """Retry strategies for different error types"""
    REDUCE_RESOURCES = "reduce_resources"  # Reduce exhaustiveness, batch size
    SWITCH_TO_CPU = "switch_to_cpu"  # Switch from GPU to CPU
    RESTART_CONTAINER = "restart_container"  # Restart the container
    WAIT_AND_RETRY = "wait_and_retry"  # Wait before retrying
    REDUCE_TIMEOUT = "reduce_timeout"  # Reduce expected runtime
    USE_FALLBACK = "use_fallback"  # Use fallback engine


@dataclass
class RetryConfig:
    """Retry configuration"""
    max_retries: int = 3
    base_delay_seconds: int = 30
    max_delay_seconds: int = 300
    exponential_backoff: bool = True
    reduce_exhaustiveness: bool = True
    switch_to_cpu_on_oom: bool = True


@dataclass
class RetryAction:
    """A retry action to take"""
    strategy: RetryStrategy
    description: str
    parameters: Dict = field(default_factory=dict)


class RetryManager:
    """
    Manages retry logic for failed docking jobs.
    Implements self-healing strategies based on error type.
    """
    
    # Strategy mapping for each error type
    STRATEGY_MAP = {
        # GPU Errors
        ErrorType.GPU_OOM: RetryStrategy.REDUCE_RESOURCES,
        ErrorType.CUDA_ERROR: RetryStrategy.SWITCH_TO_CPU,
        ErrorType.CUDA_OUT_OF_MEMORY: RetryStrategy.REDUCE_RESOURCES,
        ErrorType.GPU_NOT_FOUND: RetryStrategy.SWITCH_TO_CPU,
        
        # Container/Docker Errors
        ErrorType.CONTAINER_CRASH: RetryStrategy.RESTART_CONTAINER,
        ErrorType.CONTAINER_EXITED: RetryStrategy.RESTART_CONTAINER,
        ErrorType.DOCKER_ERROR: RetryStrategy.RESTART_CONTAINER,
        
        # System Errors
        ErrorType.SEGMENTATION_FAULT: RetryStrategy.SWITCH_TO_CPU,
        ErrorType.SYSTEM_KILL: RetryStrategy.WAIT_AND_RETRY,
        ErrorType.OUT_OF_MEMORY: RetryStrategy.REDUCE_RESOURCES,
        ErrorType.TIMEOUT: RetryStrategy.REDUCE_TIMEOUT,
        
        # Docking Errors
        ErrorType.RECEPTOR_ERROR: RetryStrategy.USE_FALLBACK,
        ErrorType.LIGAND_ERROR: RetryStrategy.USE_FALLBACK,
        ErrorType.GRID_ERROR: RetryStrategy.USE_FALLBACK,
        ErrorType.INVALID_OUTPUT: RetryStrategy.USE_FALLBACK,
        
        # IO Errors
        ErrorType.FILE_NOT_FOUND: RetryStrategy.USE_FALLBACK,
        ErrorType.PERMISSION_DENIED: RetryStrategy.USE_FALLBACK,
        ErrorType.DISK_FULL: RetryStrategy.WAIT_AND_RETRY,
        
        # Unknown
        ErrorType.UNKNOWN: RetryStrategy.WAIT_AND_RETRY,
    }
    
    def __init__(self, config: Optional[RetryConfig] = None):
        """Initialize retry manager"""
        self.config = config or RetryConfig()
        self.retry_history: Dict[str, List[RetryAction]] = {}
        
        logger.info("RetryManager initialized")
    
    def should_retry(
        self,
        job_id: str,
        failure: Failure,
        current_retry_count: int
    ) -> bool:
        """
        Determine if a job should be retried.
        
        Args:
            job_id: Job identifier
            failure: The failure that occurred
            current_retry_count: Current retry count
            
        Returns:
            True if should retry
        """
        # Check max retries
        if current_retry_count >= self.config.max_retries:
            logger.info(f"Job {job_id}: Max retries ({self.config.max_retries}) reached")
            return False
        
        # Check if recoverable
        if not failure.recoverable:
            logger.info(f"Job {job_id}: Error not recoverable")
            return False
        
        return True
    
    def get_retry_strategy(
        self,
        job_id: str,
        failure: Failure,
        current_retry_count: int
    ) -> Optional[RetryAction]:
        """
        Get the retry strategy for a failure.
        
        Args:
            job_id: Job identifier
            failure: The failure that occurred
            current_retry_count: Current retry count
            
        Returns:
            RetryAction to take, or None if should not retry
        """
        if not self.should_retry(job_id, failure, current_retry_count):
            return None
        
        # Get base strategy from error type
        strategy = self.STRATEGY_MAP.get(
            failure.error_type,
            RetryStrategy.WAIT_AND_RETRY
        )
        
        # Build action
        action = self._build_action(strategy, failure, current_retry_count)
        
        # Record in history
        if job_id not in self.retry_history:
            self.retry_history[job_id] = []
        self.retry_history[job_id].append(action)
        
        logger.info(
            f"Job {job_id}: Retry {current_retry_count + 1} - "
            f"Strategy: {strategy.value} - {action.description}"
        )
        
        return action
    
    def _build_action(
        self,
        strategy: RetryStrategy,
        failure: Failure,
        retry_count: int
    ) -> RetryAction:
        """Build a retry action"""
        
        if strategy == RetryStrategy.REDUCE_RESOURCES:
            # Reduce exhaustiveness and batch size
            exhaustiveness = max(4, 8 // (retry_count + 1))
            batch_size = max(1, 5 // (retry_count + 1))
            
            return RetryAction(
                strategy=strategy,
                description=f"Reduce resources: exhaustiveness={exhaustiveness}, batch_size={batch_size}",
                parameters={
                    "exhaustiveness": exhaustiveness,
                    "batch_size": batch_size,
                }
            )
        
        elif strategy == RetryStrategy.SWITCH_TO_CPU:
            return RetryAction(
                strategy=strategy,
                description="Switch from GPU to CPU engine",
                parameters={"engine": "vina_cpu"}
            )
        
        elif strategy == RetryStrategy.RESTART_CONTAINER:
            return RetryAction(
                strategy=strategy,
                description="Restart container before retry",
                parameters={"restart": True}
            )
        
        elif strategy == RetryStrategy.WAIT_AND_RETRY:
            # Calculate delay with exponential backoff
            delay = self._calculate_delay(retry_count)
            
            return RetryAction(
                strategy=strategy,
                description=f"Wait {delay}s then retry",
                parameters={"delay_seconds": delay}
            )
        
        elif strategy == RetryStrategy.REDUCE_TIMEOUT:
            return RetryAction(
                strategy=strategy,
                description="Reduce timeout and retry",
                parameters={"timeout_multiplier": 0.8}
            )
        
        elif strategy == RetryStrategy.USE_FALLBACK:
            return RetryAction(
                strategy=strategy,
                description="Use fallback configuration",
                parameters={"use_fallback": True}
            )
        
        # Default
        delay = self._calculate_delay(retry_count)
        return RetryAction(
            strategy=RetryStrategy.WAIT_AND_RETRY,
            description=f"Wait {delay}s then retry",
            parameters={"delay_seconds": delay}
        )
    
    def _calculate_delay(self, retry_count: int) -> int:
        """Calculate delay before retry"""
        if self.config.exponential_backoff:
            delay = self.config.base_delay_seconds * (2 ** retry_count)
        else:
            delay = self.config.base_delay_seconds
        
        return min(delay, self.config.max_delay_seconds)
    
    def apply_retry(
        self,
        job_id: str,
        action: RetryAction,
        job_config: Dict
    ) -> Dict:
        """
        Apply retry action to job configuration.
        
        Args:
            job_id: Job identifier
            action: Retry action to apply
            job_config: Current job configuration
            
        Returns:
            Updated job configuration
        """
        updated_config = job_config.copy()
        
        if action.strategy == RetryStrategy.REDUCE_RESOURCES:
            if "exhaustiveness" in action.parameters:
                updated_config["exhaustiveness"] = action.parameters["exhaustiveness"]
            if "batch_size" in action.parameters:
                updated_config["batch_size"] = action.parameters["batch_size"]
        
        elif action.strategy == RetryStrategy.SWITCH_TO_CPU:
            updated_config["engine"] = action.parameters.get("engine", "vina_cpu")
        
        elif action.strategy == RetryStrategy.REDUCE_TIMEOUT:
            multiplier = action.parameters.get("timeout_multiplier", 1.0)
            if "timeout" in updated_config:
                updated_config["timeout"] = int(updated_config["timeout"] * multiplier)
        
        logger.info(f"Applied retry action to job {job_id}: {action.description}")
        
        return updated_config
    
    def get_retry_summary(self, job_id: str) -> Dict:
        """Get retry summary for a job"""
        history = self.retry_history.get(job_id, [])
        
        return {
            "retry_count": len(history),
            "actions": [
                {
                    "strategy": a.strategy.value,
                    "description": a.description,
                    "parameters": a.parameters,
                }
                for a in history
            ]
        }
    
    def clear_history(self, job_id: Optional[str] = None):
        """Clear retry history"""
        if job_id:
            self.retry_history.pop(job_id, None)
        else:
            self.retry_history.clear()


class ResourceGuard:
    """
    Guards against resource exhaustion.
    Checks resources before launching jobs.
    """
    
    # Minimum required resources
    MIN_VRAM_MB = 1500  # Minimum 1.5GB free VRAM
    MIN_RAM_MB = 2000   # Minimum 2GB free RAM
    MIN_DISK_GB = 5     # Minimum 5GB free disk
    
    def __init__(self):
        """Initialize resource guard"""
        self.last_check: Dict = {}
        logger.info("ResourceGuard initialized")
    
    def check_resources(self) -> Dict[str, bool]:
        """
        Check if resources are available.
        
        Returns:
            Dict of resource checks
        """
        import subprocess
        
        checks = {
            "gpu_available": False,
            "vram_sufficient": False,
            "ram_sufficient": False,
            "disk_sufficient": False,
        }
        
        # Check GPU/VRAM
        try:
            result = subprocess.run(
                ['nvidia-smi', '--query-gpu=memory.free', 
                 '--format=csv,noheader,nounits'],
                capture_output=True, text=True, timeout=5
            )
            if result.returncode == 0:
                vram_free = int(result.stdout.strip().split('\n')[0])
                checks["vram_sufficient"] = vram_free >= self.MIN_VRAM_MB
                checks["gpu_available"] = True
                self.last_check["vram_free_mb"] = vram_free
        except:
            checks["gpu_available"] = False
        
        # Check RAM
        try:
            import psutil
            mem = psutil.virtual_memory()
            ram_free_mb = mem.available / (1024 * 1024)
            checks["ram_sufficient"] = ram_free_mb >= self.MIN_RAM_MB
            self.last_check["ram_free_mb"] = ram_free_mb
        except:
            pass
        
        # Check disk
        try:
            import shutil
            disk_free_gb = shutil.disk_usage('.').free / (1024**3)
            checks["disk_sufficient"] = disk_free_gb >= self.MIN_DISK_GB
            self.last_check["disk_free_gb"] = disk_free_gb
        except:
            pass
        
        self.last_check["timestamp"] = str(subprocess.check_output(['date']))
        
        all_ok = all(checks.values())
        logger.info(f"Resource check: {checks} - {'OK' if all_ok else 'INSUFFICIENT'}")
        
        return checks
    
    def can_launch_gpu_job(self) -> Tuple[bool, str]:
        """Check if can launch GPU job"""
        checks = self.check_resources()
        
        if not checks["gpu_available"]:
            return False, "No GPU available"
        
        if not checks["vram_sufficient"]:
            vram = self.last_check.get("vram_free_mb", 0)
            return False, f"Insufficient VRAM ({vram}MB < {self.MIN_VRAM_MB}MB)"
        
        if not checks["ram_sufficient"]:
            return False, "Insufficient system RAM"
        
        if not checks["disk_sufficient"]:
            return False, "Insufficient disk space"
        
        return True, "OK"
    
    def get_status_summary(self) -> str:
        """Get status summary for UI"""
        if not self.last_check:
            return "Not checked"
        
        parts = []
        if "vram_free_mb" in self.last_check:
            parts.append(f"VRAM: {self.last_check['vram_free_mb']}MB")
        if "ram_free_mb" in self.last_check:
            parts.append(f"RAM: {self.last_check['ram_free_mb']:.0f}MB")
        if "disk_free_gb" in self.last_check:
            parts.append(f"Disk: {self.last_check['disk_free_gb']:.1f}GB")
        
        return " | ".join(parts) if parts else "Unknown"
