"""
BioDockify Agent Zero - Intelligent Failure Detection and Self-Repair System
Handles failure detection, diagnosis, and recovery strategies
"""

import logging
import json
from typing import Dict, List, Optional, Any, Callable
from enum import Enum
from dataclasses import dataclass
from datetime import datetime
import uuid

from src.checkpoint_manager import CheckpointManager
from src.recovery_manager import RecoveryManager

logger = logging.getLogger(__name__)

class FailureType(Enum):
    """Enumeration of all failure types"""
    # File-related failures
    FILE_CORRUPTION = "file_corruption"
    FILE_MISSING = "file_missing"
    FILE_INVALID_FORMAT = "file_invalid_format"
    FILE_UNSUPPORTED_ATOMS = "file_unsupported_atoms"
    
    # Docker-related failures
    DOCKER_NOT_RUNNING = "docker_not_running"
    DOCKER_OUT_OF_MEMORY = "docker_out_of_memory"
    DOCKER_PERMISSION_DENIED = "docker_permission_denied"
    DOCKER_TIMEOUT = "docker_timeout"
    DOCKER_CONTAINER_CRASHED = "docker_container_crashed"
    
    # Docking-related failures
    DOCKING_NO_POSES = "docking_no_poses"
    DOCKING_TIMEOUT = "docking_timeout"
    DOCKING_OUT_OF_MEMORY = "docking_out_of_memory"
    
    # Analysis-related failures
    ODDT_FAILED = "oddt_failed"
    RDKIT_FAILED = "rdkit_failed"
    DESCRIPTOR_CALCULATION_FAILED = "descriptor_calculation_failed"
    
    # Output-related failures
    OUTPUT_FILE_MISSING = "output_file_missing"
    OUTPUT_FILE_CORRUPTED = "output_file_corrupted"
    OUTPUT_WRITE_FAILED = "output_write_failed"

class FailureSeverity(Enum):
    """Severity levels for failures"""
    LOW = 1
    MEDIUM = 2
    HIGH = 3
    CRITICAL = 4

@dataclass
class FailureInfo:
    """Information about a detected failure"""
    failure_type: FailureType
    severity: FailureSeverity
    message: str
    diagnosis: str
    recovery_action: str
    stage: str
    occurred_at: datetime
    job_id: str

class AgentZero:
    """Intelligent failure detection and self-repair system"""
    
    def __init__(self, checkpoint_manager: CheckpointManager, recovery_manager: RecoveryManager):
        """Initialize Agent Zero"""
        self.checkpoint_manager = checkpoint_manager
        self.recovery_manager = recovery_manager
        self.failure_signatures = self._get_failure_signatures()
        self.active_failures = {}
        self.confidence_score = 100
    
    def _get_failure_signatures(self) -> Dict[FailureType, Dict[str, Any]]:
        """Get unique signatures for each failure type"""
        return {
            FailureType.FILE_CORRUPTION: {
                "patterns": ["corrupted", "invalid", "truncated", "empty"],
                "severity": FailureSeverity.MEDIUM
            },
            FailureType.DOCKER_NOT_RUNNING: {
                "patterns": ["docker not running", "connection refused", "daemon not found"],
                "severity": FailureSeverity.HIGH
            },
            FailureType.DOCKING_TIMEOUT: {
                "patterns": ["timeout", "time limit exceeded", "too long"],
                "severity": FailureSeverity.MEDIUM
            },
            FailureType.DOCKING_NO_POSES: {
                "patterns": ["no poses", "zero results", "empty output", "failed"],
                "severity": FailureSeverity.HIGH
            }
        }
    
    def detect_failure(self, error_message: str, stage: str, job_id: str) -> Optional[FailureInfo]:
        """Detect failure from error message and stage"""
        logger.debug(f"Detecting failure: {error_message} at stage {stage}")
        
        # Check against known failure patterns
        for failure_type, signature in self.failure_signatures.items():
            for pattern in signature["patterns"]:
                if pattern.lower() in error_message.lower():
                    return self._create_failure_info(failure_type, signature, error_message, stage, job_id)
        
        # Check for specific patterns
        if "memory" in error_message.lower():
            return self._create_failure_info(FailureType.DOCKING_OUT_OF_MEMORY, 
                                                    self.failure_signatures.get(FailureType.DOCKING_OUT_OF_MEMORY, {"severity": FailureSeverity.HIGH}),
                                                    error_message, stage, job_id)
        
        if "permission" in error_message.lower():
            return self._create_failure_info(FailureType.DOCKER_PERMISSION_DENIED,
                                                    self.failure_signatures.get(FailureType.DOCKER_PERMISSION_DENIED, {"severity": FailureSeverity.HIGH}),
                                                    error_message, stage, job_id)
        
        # Default unknown failure
        return self._create_failure_info(FailureType.OUTPUT_FILE_MISSING,
                                                {"patterns": ["unknown"], "severity": FailureSeverity.LOW},
                                                error_message, stage, job_id)
    
    def _create_failure_info(self, failure_type: FailureType, signature: Dict[str, Any],
                          error_message: str, stage: str, job_id: str) -> FailureInfo:
        """Create failure information object"""
        return FailureInfo(
            failure_type=failure_type,
            severity=signature.get("severity", FailureSeverity.MEDIUM),
            message=self._generate_user_message(failure_type, error_message),
            diagnosis=self._generate_diagnosis(failure_type, error_message),
            recovery_action=self._select_recovery_strategy(failure_type),
            stage=stage,
            occurred_at=datetime.now(),
            job_id=job_id
        )
    
    def _generate_user_message(self, failure_type: FailureType, error_message: str) -> str:
        """Generate user-friendly message (scientific tone, no jargon)"""
        messages = {
            FailureType.FILE_CORRUPTION: f"Ligand file contains unsupported atoms. Regenerating ligand structure.",
            FailureType.DOCKER_NOT_RUNNING: "Docker Desktop is required for docking functionality. Please start Docker Desktop.",
            FailureType.DOCKING_TIMEOUT: "Docking simulation exceeded time limit. Reducing exhaustiveness and retrying.",
            FailureType.DOCKING_NO_POSES: "No binding poses found. Expanding search box and retrying.",
            FailureType.OUTPUT_FILE_MISSING: "Output file not found. Rerunning docking simulation.",
            FailureType.ODDT_FAILED: "Interaction analysis failed. Results will be displayed without interaction details.",
            FailureType.RDKIT_FAILED: "Descriptor calculation failed. Results will be displayed without molecular descriptors.",
            FailureType.DOCKER_OUT_OF_MEMORY: "Docker out of memory. Reducing complexity and retrying."
        }
        
        return messages.get(failure_type, f"An error occurred: {error_message}")
    
    def _generate_diagnosis(self, failure_type: FailureType, error_message: str) -> str:
        """Generate diagnosis for failure"""
        diagnoses = {
            FailureType.FILE_CORRUPTION: "File structure contains unsupported atoms or is corrupted. System will regenerate structure.",
            FailureType.DOCKER_NOT_RUNNING: "Docker Desktop daemon is not running or not accessible.",
            FailureType.DOCKING_TIMEOUT: "Docking simulation exceeded configured time limit.",
            FailureType.DOCKING_NO_POSES: "Search box parameters may be too small or receptor-ligand combination may be incompatible.",
            FailureType.OUTPUT_FILE_MISSING: "Docking simulation completed but output file was not generated.",
            FailureType.ODDT_FAILED: "ODDT (Open Drug Discovery Toolkit) analysis encountered an error.",
            FailureType.RDKIT_FAILED: "RDKit (Chemistry Development Kit) descriptor calculation encountered an error.",
            FailureType.DOCKER_OUT_OF_MEMORY: "Docker container memory allocation is insufficient."
        }
        
        return diagnoses.get(failure_type, f"Unable to diagnose: {error_message}")
    
    def _select_recovery_strategy(self, failure_type: FailureType) -> str:
        """Select appropriate recovery strategy based on failure type"""
        strategies = {
            FailureType.FILE_CORRUPTION: "regenerate_file",
            FailureType.DOCKER_NOT_RUNNING: "wait_for_docker",
            FailureType.DOCKING_TIMEOUT: "retry_modified",
            FailureType.DOCKING_NO_POSES: "retry_modified",
            FailureType.OUTPUT_FILE_MISSING: "rerun_docking",
            FailureType.ODDT_FAILED: "skip_analysis",
            FailureType.RDKIT_FAILED: "skip_analysis",
            FailureType.DOCKER_OUT_OF_MEMORY: "retry_modified"
        }
        
        return strategies.get(failure_type, "manual_intervention")
    
    def attempt_recovery(self, failure_info: FailureInfo, job_id: str, 
                     current_stage: str) -> bool:
        """Attempt recovery using appropriate strategy"""
        logger.info(f"Attempting recovery for {failure_info.failure_type} at stage {current_stage}")
        
        recovery_strategy = failure_info.recovery_action
        
        # Save checkpoint before recovery attempt
        self.checkpoint_manager.save_checkpoint(job_id, current_stage, {
            "failure_type": failure_info.failure_type.value,
            "recovery_strategy": recovery_strategy,
            "failure_message": failure_info.message
        })
        
        # Adjust confidence score based on failure
        penalty = 0
        if failure_info.severity == FailureSeverity.HIGH:
            penalty = 15
        elif failure_info.severity == FailureSeverity.MEDIUM:
            penalty = 5
            
        # Refinement #5: Cap penalty to prevent excessive drops
        penalty = min(penalty, 20)
        
        self.confidence_score -= penalty
        
        # Ensure confidence doesn't go below 0
        # Enforce bounds: 0 <= confidence <= 100
        self.confidence_score = max(0, min(100, self.confidence_score))
        
        # Execute recovery strategy
        recovery_success, detail_message = self.recovery_manager.execute_recovery(
            recovery_strategy=recovery_strategy,
            failure_info=failure_info,
            job_id=job_id,
            current_stage=current_stage
        )
        
        if recovery_success:
            self.reward_success()
            logger.info(f"Recovery successful: {detail_message}")
            # Update failure info message with detail for UI display if needed
            failure_info.message = f"{failure_info.message} ({detail_message})"
            return True
        else:
            logger.warning(f"Recovery failed: {detail_message}")
            return False

    def reward_success(self, amount: int = 5):
        """Reward successful recovery"""
        self.confidence_score += amount
        # Enforce bounds
        self.confidence_score = max(0, min(100, self.confidence_score))
        logger.info(f"Confidence increased by {amount} to {self.confidence_score}")
    
    def get_confidence_score(self) -> int:
        """Get current confidence score"""
        return self.confidence_score
    
    def reset_confidence(self) -> None:
        """Reset confidence score to 100"""
        self.confidence_score = 100
        logger.info("Confidence score reset to 100")
