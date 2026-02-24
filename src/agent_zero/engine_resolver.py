"""
Agent Zero v4 - Safe Engine Resolution
Auto-switches engines when GPU fails, with predefined safe actions.
"""

import logging
from typing import Optional, Dict, Any, List, Callable
from dataclasses import dataclass
from enum import Enum
import time

logger = logging.getLogger(__name__)


class FailureType(Enum):
    """Types of execution failures"""
    GPU_OOM = "gpu_oom"
    CUDA_ERROR = "cuda_error"
    CONTAINER_CRASH = "container_crash"
    TIMEOUT = "timeout"
    SEGFAULT = "segfault"
    UNKNOWN = "unknown"


class EngineType(Enum):
    """Available engine types"""
    VINA_CPU = "vina-cpu"
    VINA_GPU = "vina-gpu"
    GNINA_GPU = "gnina-gpu"


@dataclass
class SafeAction:
    """A predefined safe action that can be taken"""
    name: str
    description: str
    apply: Callable[[Dict[str, Any]], Dict[str, Any]]


class SafeEngineResolver:
    """
    Resolves engine based on compute selection and handles failures.
    Implements auto-downgrade and predefined safe actions.
    """
    
    # Engine variants for each compute mode
    ENGINE_MAP = {
        "GPU": {
            "primary": "gnina-gpu",
            "fallback": "vina-cpu"
        },
        "CPU": {
            "primary": "vina-cpu",
            "fallback": None  # No further fallback
        }
    }
    
    def __init__(self):
        self._has_downgraded = False
        self._downgrade_count = 0
        self._max_downgrades = 2
    
    def resolve_engine(self, compute_mode: str, job_config: Dict[str, Any]) -> str:
        """
        Resolve the appropriate engine based on compute mode.
        
        Args:
            compute_mode: "GPU" or "CPU"
            job_config: Job configuration
            
        Returns:
            Engine type string
        """
        if compute_mode == "GPU":
            return self.ENGINE_MAP["GPU"]["primary"]
        else:
            return self.ENGINE_MAP["CPU"]["primary"]
    
    def handle_failure(
        self, 
        failure_type: FailureType, 
        current_engine: str,
        job_config: Dict[str, Any]
    ) -> Optional[Dict[str, Any]]:
        """
        Handle runtime failure with auto-downgrade.
        
        Args:
            failure_type: Type of failure that occurred
            current_engine: Engine that failed
            job_config: Current job configuration
            
        Returns:
            New configuration with adjusted parameters, or None if should stop
        """
        self._downgrade_count += 1
        
        logger.warning(f"Handling failure: {failure_type.value} with engine {current_engine}")
        
        # If already downgraded too many times, stop
        if self._downgrade_count > self._max_downgrades:
            logger.error("Maximum downgrades exceeded, stopping job")
            return None
        
        # Determine if we should switch to CPU
        should_switch_to_cpu = self._should_switch_to_cpu(failure_type)
        
        if should_switch_to_cpu and current_engine != "vina-cpu":
            logger.info("Switching to CPU engine due to GPU failure")
            self._has_downgraded = True
            return self._create_cpu_config(job_config)
        
        # Try to adjust parameters without switching engine
        adjusted = self._apply_safe_adjustments(failure_type, job_config)
        
        if adjusted:
            return adjusted
        
        # If adjustments fail, try CPU fallback
        if current_engine != "vina-cpu":
            return self._create_cpu_config(job_config)
        
        # CPU failed too - no more options
        return None
    
    def _should_switch_to_cpu(self, failure_type: FailureType) -> bool:
        """Determine if failure warrants switching to CPU"""
        gpu_failure_types = {
            FailureType.GPU_OOM,
            FailureType.CUDA_ERROR,
        }
        return failure_type in gpu_failure_types
    
    def _create_cpu_config(self, job_config: Dict[str, Any]) -> Dict[str, Any]:
        """Create CPU-compatible configuration"""
        config = job_config.copy()
        config["engine"] = "vina-cpu"
        config["compute_mode"] = "CPU"
        config["exhaustiveness"] = min(
            job_config.get("exhaustiveness", 8) + 4,  # Increase exhaustiveness for CPU
            32  # Cap at 32
        )
        return config
    
    def _apply_safe_adjustments(
        self, 
        failure_type: FailureType, 
        job_config: Dict[str, Any]
    ) -> Optional[Dict[str, Any]]:
        """
        Apply predefined safe adjustments based on failure type.
        
        Returns adjusted config or None if no adjustment possible.
        """
        config = job_config.copy()
        
        if failure_type == FailureType.GPU_OOM:
            # Reduce batch size and exhaustiveness
            if "cnn_batch_size" in config:
                config["cnn_batch_size"] = max(1, config.get("cnn_batch_size", 8) // 2)
            if "exhaustiveness" in config:
                config["exhaustiveness"] = max(4, config.get("exhaustiveness", 8) // 2)
            config["reduced_for_oom"] = True
            logger.info("Applied OOM adjustments: reduced batch and exhaustiveness")
            return config
            
        elif failure_type == FailureType.CUDA_ERROR:
            # Clear CUDA cache hint, reduce parameters
            if "cnn_batch_size" in config:
                config["cnn_batch_size"] = max(1, config.get("cnn_batch_size", 8) // 2)
            config["clear_cuda_cache"] = True
            logger.info("Applied CUDA error adjustments")
            return config
            
        elif failure_type == FailureType.TIMEOUT:
            # Reduce exhaustiveness to speed up
            if "exhaustiveness" in config:
                config["exhaustiveness"] = max(4, config.get("exhaustiveness", 8) - 2)
            config["num_modes"] = min(config.get("num_modes", 9), 3)
            logger.info("Applied timeout adjustments")
            return config
        
        return None
    
    def get_safe_actions(self) -> List[Dict[str, str]]:
        """
        Get list of predefined safe actions (for documentation/UI).
        
        Returns:
            List of safe action descriptions
        """
        return [
            {
                "name": "reduce_exhaustiveness",
                "description": "Reduce exhaustiveness to lower computation time"
            },
            {
                "name": "reduce_cnn_batch", 
                "description": "Reduce CNN batch size to lower VRAM usage"
            },
            {
                "name": "switch_gpu_to_cpu",
                "description": "Switch from GPU engine to CPU engine"
            },
            {
                "name": "restart_container",
                "description": "Restart the Docker container to clear state"
            },
            {
                "name": "delay_job",
                "description": "Delay job execution to wait for resources"
            },
            {
                "name": "lower_parallelism",
                "description": "Reduce number of parallel jobs"
            },
            {
                "name": "clear_cuda_cache",
                "description": "Clear CUDA cache before retry"
            }
        ]
    
    def reset(self):
        """Reset downgrade state for new job"""
        self._has_downgraded = False
        self._downgrade_count = 0
    
    @property
    def has_downgraded(self) -> bool:
        return self._has_downgraded
    
    @property
    def downgrade_count(self) -> int:
        return self._downgrade_count


def parse_failure_type(error_message: str) -> FailureType:
    """
    Parse failure type from error message.
    
    Args:
        error_message: Error message from logs
        
    Returns:
        FailureType enum value
    """
    error_lower = error_message.lower()
    
    if "out of memory" in error_lower or "oom" in error_lower:
        return FailureType.GPU_OOM
    elif "cuda" in error_lower and "error" in error_lower:
        return FailureType.CUDA_ERROR
    elif "segfault" in error_lower or "segmentation" in error_lower:
        return FailureType.SEGFAULT
    elif "timeout" in error_lower:
        return FailureType.TIMEOUT
    elif "container" in error_lower and "exit" in error_lower:
        return FailureType.CONTAINER_CRASH
    else:
        return FailureType.UNKNOWN


# Singleton
_resolver: Optional[SafeEngineResolver] = None

def get_engine_resolver() -> SafeEngineResolver:
    """Get singleton engine resolver"""
    global _resolver
    if _resolver is None:
        _resolver = SafeEngineResolver()
    return _resolver
