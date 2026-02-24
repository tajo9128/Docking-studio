"""
Agent Zero v4 - Trusted Autonomous Orchestrator
Integrates all v4 components into unified execution flow.
"""

import logging
from typing import Optional, Dict, Any
from dataclasses import dataclass

from src.agent_zero.compute_selector import (
    ComputeSelector, ComputeMode, get_compute_selector, auto_select_compute
)
from src.agent_zero.resource_guard import ResourceGuard
from src.agent_zero.engine_resolver import (
    SafeEngineResolver, FailureType, parse_failure_type, get_engine_resolver
)
from src.agent_zero.internet_diagnostics import (
    InternetDiagnostics, DiagnosticsMode, get_internet_diagnostics
)

logger = logging.getLogger(__name__)


@dataclass
class ExecutionPlan:
    """Complete execution plan from v4 orchestration"""
    compute_mode: str  # "GPU" or "CPU"
    engine: str
    config: Dict[str, Any]
    can_proceed: bool
    resource_check_passed: bool
    warnings: list
    recommended_actions: list


class AgentZeroV4:
    """
    Trusted Autonomous Orchestrator for Agent Zero v4.
    
    Integrates:
    - Compute Selector (auto GPU/CPU selection)
    - Resource Guard (pre-flight checks)
    - Safe Engine Resolver (auto-downgrade on failure)
    - Internet Diagnostics (advisory-only)
    """
    
    def __init__(self):
        self.compute_selector = get_compute_selector()
        self.resource_guard = ResourceGuard()
        self.engine_resolver = get_engine_resolver()
        self.internet_diagnostics = get_internet_diagnostics()
        
        self._compute_mode = ComputeMode.AUTO
        self._last_failure: Optional[FailureType] = None
        
        logger.info("Agent Zero v4 initialized")
    
    def set_compute_mode(self, mode: str):
        """Set compute mode: 'auto', 'gpu_force', 'cpu_force'"""
        if mode == "auto":
            self._compute_mode = ComputeMode.AUTO
        elif mode == "gpu_force":
            self._compute_mode = ComputeMode.GPU_FORCE
        elif mode == "cpu_force":
            self._compute_mode = ComputeMode.CPU_FORCE
        logger.info(f"Compute mode set to: {mode}")
    
    def get_compute_mode(self) -> str:
        """Get current compute mode"""
        return self._compute_mode.value
    
    def plan_execution(self, job_config: Dict[str, Any]) -> ExecutionPlan:
        """
        Create execution plan using autonomous selection.
        
        Args:
            job_config: Job configuration dictionary
            
        Returns:
            ExecutionPlan with compute decision, engine, and safety checks
        """
        # Step 1: Autonomous compute selection
        compute_decision = self.compute_selector.choose_best(self._compute_mode)
        
        # Step 2: Resource guard pre-flight checks
        caps = self.resource_guard.get_gpu_capabilities()
        
        # Build a simple resource check result
        class ResourceCheckResult:
            def __init__(self, caps):
                # CPU mode is always available
                self.can_launch = True
                self.warnings = []
                self.recommended_actions = []
                if not caps.supports_vina_gpu:
                    self.warnings.append("GPU not available - using CPU mode")
                if caps.vram_level.value == "low" and caps.total_vram_mb > 0:
                    self.warnings.append(f"Low VRAM ({caps.free_vram_mb}MB)")
        
        resource_check = ResourceCheckResult(caps)
        
        # Step 3: Resolve engine
        engine = self.engine_resolver.resolve_engine(
            compute_decision.mode, 
            job_config
        )
        
        # Build execution config
        config = job_config.copy()
        config["compute_mode"] = compute_decision.mode
        config["engine"] = engine
        config["confidence"] = compute_decision.confidence
        
        warnings = []
        recommended_actions = []
        
        # Add compute selection warning if confidence low
        if compute_decision.confidence < 0.7:
            warnings.append(f"Low confidence ({compute_decision.confidence:.0%}): {compute_decision.reason}")
        
        # Add resource warnings
        warnings.extend(resource_check.warnings)
        recommended_actions.extend(resource_check.recommended_actions)
        
        # Determine if can proceed
        can_proceed = resource_check.can_launch
        
        # If GPU selected but no Docker GPU access, switch to CPU
        if compute_decision.mode == "GPU" and not self._check_docker_gpu():
            logger.warning("Docker GPU not accessible, switching to CPU")
            compute_decision.mode = "CPU"
            engine = "vina-cpu"
            config["compute_mode"] = "CPU"
            config["engine"] = "vina-cpu"
            warnings.append("Docker GPU not accessible - using CPU")
        
        return ExecutionPlan(
            compute_mode=compute_decision.mode,
            engine=engine,
            config=config,
            can_proceed=can_proceed,
            resource_check_passed=resource_check.can_launch,
            warnings=warnings,
            recommended_actions=recommended_actions
        )
    
    def handle_failure(
        self, 
        error_message: str, 
        current_config: Dict[str, Any]
    ) -> Optional[Dict[str, Any]]:
        """
        Handle execution failure with auto-recovery.
        
        Args:
            error_message: Error message from execution
            current_config: Current job configuration
            
        Returns:
            New configuration for retry, or None if should stop
        """
        # Parse failure type
        failure_type = parse_failure_type(error_message)
        self._last_failure = failure_type
        
        logger.info(f"Handling failure: {failure_type.value}")
        
        # Get advisory suggestions
        advisory = self.internet_diagnostics.lookup(error_message)
        if advisory.found:
            logger.info(f"Advisory suggestions: {advisory.suggestions}")
        
        # Handle through engine resolver
        new_config = self.engine_resolver.handle_failure(
            failure_type,
            current_config.get("engine", "unknown"),
            current_config
        )
        
        if new_config:
            logger.info("Retrying with adjusted configuration")
        else:
            logger.error("No more recovery options - job failed")
        
        return new_config
    
    def get_status(self) -> Dict[str, Any]:
        """Get comprehensive v4 system status"""
        compute_status = self.compute_selector.get_status_summary()
        resource_status = self.resource_guard.get_status_summary()
        
        return {
            "compute_mode": self._compute_mode.value,
            "compute": compute_status,
            "resources": resource_status,
            "failure_history": self._last_failure.value if self._last_failure else None,
            "safe_actions": self.engine_resolver.get_safe_actions()
        }
    
    def reset(self):
        """Reset state for new job"""
        self.engine_resolver.reset()
        self._last_failure = None
    
    def _check_docker_gpu(self) -> bool:
        """Check if Docker can access GPU"""
        import subprocess
        try:
            result = subprocess.run(
                ["docker", "run", "--rm", "--gpus", "all", 
                 "nvidia/cuda:12.0-base", "nvidia-smi"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=30
            )
            return result.returncode == 0
        except Exception:
            return False


# Singleton
_v4_orchestrator: Optional[AgentZeroV4] = None

def get_v4_orchestrator() -> AgentZeroV4:
    """Get singleton v4 orchestrator"""
    global _v4_orchestrator
    if _v4_orchestrator is None:
        _v4_orchestrator = AgentZeroV4()
    return _v4_orchestrator
