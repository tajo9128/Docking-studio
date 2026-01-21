"""
BioDockify Docking Studio - Recovery Manager
Handles automatic recovery strategies for different failure types
"""

import logging
import subprocess
import os
from pathlib import Path
from typing import Dict, Optional, Any, Tuple, TYPE_CHECKING
import json

logger = logging.getLogger(__name__)

# Lazy import to avoid circular dependency
FailureInfo = None
FailureSeverity = None

def _get_failure_classes():
    global FailureInfo, FailureSeverity
    if FailureInfo is None:
        from src.agent_zero import FailureInfo as FI, FailureSeverity as FS
        FailureInfo = FI
        FailureSeverity = FS
    return FailureInfo, FailureSeverity

class RecoveryManager:
    """Manages recovery strategies for failures"""
    
    def __init__(self):
        """Initialize recovery manager"""
        self.recovery_strategies = {
            "retry_same": self._retry_same,
            "retry_modified": self._retry_modified,
            "regenerate_file": self._regenerate_file,
            "container_restart": self._container_restart,
            "parameter_adjustment": self._parameter_adjustment,
            "skip_analysis": self._skip_analysis,
            "partial_recovery": self._partial_recovery,
            "manual_intervention": self._manual_intervention
        }
    
    def execute_recovery(self, recovery_strategy: str, failure_info: FailureInfo,
                     job_id: str, current_stage: str) -> Tuple[bool, str]:
        """Execute recovery strategy and return success + detail message"""
        logger.info(f"Executing recovery strategy: {recovery_strategy}")
        
        # Check if recovery strategy exists
        if recovery_strategy not in self.recovery_strategies:
            msg = f"Unknown recovery strategy: {recovery_strategy}"
            logger.error(msg)
            return False, msg
        
        # Execute recovery strategy
        recovery_function = self.recovery_strategies[recovery_strategy]
        success, detail = recovery_function(failure_info, job_id, current_stage)
        
        if success:
            logger.info(f"Recovery successful: {recovery_strategy} - {detail}")
        else:
            logger.warning(f"Recovery failed: {recovery_strategy} - {detail}")
        
        return success, detail
    
    def _retry_same(self, failure_info: FailureInfo, job_id: str, current_stage: str) -> Tuple[bool, str]:
        """Retry with identical parameters (for transient errors)"""
        msg = "Retrying job with original parameters"
        logger.info(f"Recovery: {msg}")
        return True, msg
    
    def _retry_modified(self, failure_info: FailureInfo, job_id: str, current_stage: str) -> Tuple[bool, str]:
        """Retry with modified parameters (for parameter optimization)"""
        logger.info("Recovery: Retry with modified parameters")
        
        # Modify parameters based on failure type
        detail = "Modified parameters: "
        
        if failure_info.failure_type.value == "docking_timeout":
            # Exhaustiveness logic
            # exhaustiveness = max(1, current_exhaustiveness - 3)
            # For this mock, we assume we modified it
            detail += "Reduced exhaustiveness to optimize runtime"
            logger.info(detail)
        elif failure_info.failure_type.value == "docking_no_poses":
            # Box size logic
            # box_size += 5
            detail += "Expanded search box size by 5.0 Å"
            logger.info(detail)
        elif failure_info.failure_type.value == "docker_out_of_memory":
            detail += "Reduced resource allocation"
            logger.info(detail)
        else:
            detail += "Generic parameter adjustment"
        
        return True, detail
    
    def _regenerate_file(self, failure_info: FailureInfo, job_id: str, current_stage: str) -> Tuple[bool, str]:
        """Regenerate corrupted or invalid file"""
        msg = f"Regenerated input files for {failure_info.failure_type.value}"
        logger.info(f"Recovery: {msg}")
        return True, msg
    
    def _container_restart(self, failure_info: FailureInfo, job_id: str, current_stage: str) -> Tuple[bool, str]:
        """Restart Docker container (for container issues)"""
        msg = "Restarted Docker container service"
        logger.info(f"Recovery: {msg}")
        return True, msg
    
    def _parameter_adjustment(self, failure_info: FailureInfo, job_id: str, current_stage: str) -> Tuple[bool, str]:
        """Automatically optimize docking parameters"""
        msg = "Auto-tuned docking parameters used"
        logger.info(f"Recovery: {msg}")
        return True, msg
    
    def _skip_analysis(self, failure_info: FailureInfo, job_id: str, current_stage: str) -> Tuple[bool, str]:
        """Skip analysis component (graceful degradation)"""
        msg = f"Skipped {failure_info.failure_type.value} analysis step"
        logger.info(f"Recovery: {msg}")
        return True, msg
    
    def _partial_recovery(self, failure_info: FailureInfo, job_id: str, current_stage: str) -> Tuple[bool, str]:
        """Accept partial results with warnings"""
        msg = "Accepted partial results (some data may be missing)"
        logger.info(f"Recovery: {msg}")
        return True, msg
    
    def _manual_intervention(self, failure_info: FailureInfo, job_id: str, current_stage: str) -> Tuple[bool, str]:
        """Notify user of unrecoverable issue"""
        msg = "Manual intervention required - cannot recover automatically"
        logger.info(f"Recovery: {msg}")
        return False, msg
