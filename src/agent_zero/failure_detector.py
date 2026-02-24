"""
Failure Detector - Bug Detection for Docking Jobs
Rule-based failure scanner for log analysis
"""

import re
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
from enum import Enum
import logging

logger = logging.getLogger(__name__)


class ErrorType(Enum):
    """Types of errors that can occur during docking"""
    # GPU Errors
    GPU_OOM = "GPU_OOM"
    CUDA_ERROR = "CUDA_ERROR"
    CUDA_OUT_OF_MEMORY = "CUDA_OUT_OF_MEMORY"
    GPU_NOT_FOUND = "GPU_NOT_FOUND"
    
    # Container/Docker Errors
    CONTAINER_CRASH = "CONTAINER_CRASH"
    CONTAINER_EXITED = "CONTAINER_EXITED"
    DOCKER_ERROR = "DOCKER_ERROR"
    
    # System Errors
    SEGMENTATION_FAULT = "SEGMENTATION_FAULT"
    SYSTEM_KILL = "SYSTEM_KILL"
    OUT_OF_MEMORY = "OUT_OF_MEMORY"
    TIMEOUT = "TIMEOUT"
    
    # Docking Errors
    RECEPTOR_ERROR = "RECEPTOR_ERROR"
    LIGAND_ERROR = "LIGAND_ERROR"
    GRID_ERROR = "GRID_ERROR"
    INVALID_OUTPUT = "INVALID_OUTPUT"
    
    # Network/IO Errors
    FILE_NOT_FOUND = "FILE_NOT_FOUND"
    PERMISSION_DENIED = "PERMISSION_DENIED"
    DISK_FULL = "DISK_FULL"
    
    # Unknown
    UNKNOWN = "UNKNOWN"


@dataclass
class Failure:
    """A detected failure"""
    error_type: ErrorType
    message: str
    severity: str  # critical, error, warning
    line: Optional[str] = None
    timestamp: Optional[str] = None
    recoverable: bool = True
    suggested_fix: str = ""


class FailureDetector:
    """
    Rule-based failure detection for docking jobs.
    Analyzes log output to identify issues.
    """
    
    # Error patterns to search for
    ERROR_PATTERNS = {
        # GPU Errors
        ErrorType.GPU_OOM: [
            r"CUDA out of memory",
            r"out of memory",
            r"OutOfMemoryError",
            r"GPU memory.*exceeded",
            r"memory allocation.*failed",
        ],
        ErrorType.CUDA_ERROR: [
            r"CUDA error",
            r"cuda.*error",
            r"NVidia driver.*not",
            r"CUDA_VISIBLE_DEVICES",
        ],
        
        # Container/System Errors
        ErrorType.SEGMENTATION_FAULT: [
            r"Segmentation fault",
            r"segfault",
            r"SIGSEGV",
        ],
        ErrorType.SYSTEM_KILL: [
            r"killed",
            r"Killed",
            r"SIGKILL",
            r"Memory cgroup out of memory",
            r"oom-killer",
        ],
        ErrorType.CONTAINER_CRASH: [
            r"container.*exited",
            r"exit code.*1",
            r"exit code.*137",
            r"exit code.*139",
        ],
        
        # Docking Specific Errors
        ErrorType.RECEPTOR_ERROR: [
            r"cannot read receptor",
            r"invalid receptor",
            r"receptor.*error",
            r"failed to load receptor",
        ],
        ErrorType.LIGAND_ERROR: [
            r"cannot read ligand",
            r"invalid ligand",
            r"ligand.*error",
            r"failed to load ligand",
            r"no ligands found",
        ],
        ErrorType.GRID_ERROR: [
            r"invalid grid",
            r"grid.*error",
            r"center.*invalid",
            r"size.*invalid",
        ],
        
        # File Errors
        ErrorType.FILE_NOT_FOUND: [
            r"No such file",
            r"file not found",
            r"cannot open",
            r"ENOENT",
        ],
        ErrorType.PERMISSION_DENIED: [
            r"permission denied",
            r"cannot access",
            r"EACCES",
        ],
        ErrorType.DISK_FULL: [
            r"No space left",
            r"disk full",
            r"ENOSPC",
        ],
        
        # Timeout
        ErrorType.TIMEOUT: [
            r"timeout",
            r"timed out",
            r"deadline exceeded",
        ],
        
        # Additional GPU Errors
        ErrorType.CUDA_OUT_OF_MEMORY: [
            r"CUDA out of memory",
            r"cuda.*out of memory",
            r"CUDA.*OOM",
        ],
        ErrorType.GPU_NOT_FOUND: [
            r"GPU.*not found",
            r"no.*GPU",
            r"nvidia-smi.*not found",
            r"cannot access.*GPU",
        ],
        
        # Container Exit
        ErrorType.CONTAINER_EXITED: [
            r"container.*exited",
            r"exited with code",
        ],
        
        # Docker Error
        ErrorType.DOCKER_ERROR: [
            r"docker.*error",
            r"docker.*failed",
            r"docker.*exception",
        ],
        
        # Invalid Output
        ErrorType.INVALID_OUTPUT: [
            r"invalid output",
            r"output.*error",
            r"cannot write output",
            r"failed to write",
        ],
    }
    
    # Severity mapping
    SEVERITY_MAP = {
        # GPU Errors
        ErrorType.GPU_OOM: "critical",
        ErrorType.CUDA_ERROR: "critical",
        ErrorType.CUDA_OUT_OF_MEMORY: "critical",
        ErrorType.GPU_NOT_FOUND: "error",
        
        # Container/Docker Errors
        ErrorType.CONTAINER_CRASH: "error",
        ErrorType.CONTAINER_EXITED: "error",
        ErrorType.DOCKER_ERROR: "error",
        
        # System Errors
        ErrorType.SEGMENTATION_FAULT: "critical",
        ErrorType.SYSTEM_KILL: "critical",
        ErrorType.OUT_OF_MEMORY: "critical",
        ErrorType.TIMEOUT: "warning",
        
        # Docking Errors
        ErrorType.RECEPTOR_ERROR: "error",
        ErrorType.LIGAND_ERROR: "error",
        ErrorType.GRID_ERROR: "error",
        ErrorType.INVALID_OUTPUT: "error",
        
        # IO Errors
        ErrorType.FILE_NOT_FOUND: "error",
        ErrorType.PERMISSION_DENIED: "error",
        ErrorType.DISK_FULL: "critical",
    }
    
    # Recovery suggestions
    FIX_SUGGESTIONS = {
        # GPU Errors
        ErrorType.GPU_OOM: "Reduce exhaustiveness, reduce batch size, or switch to CPU",
        ErrorType.CUDA_ERROR: "Check GPU drivers, try CPU mode instead",
        ErrorType.CUDA_OUT_OF_MEMORY: "Reduce CNN batch size, clear CUDA cache, or use CPU",
        ErrorType.GPU_NOT_FOUND: "Check GPU availability, use CPU mode",
        
        # Container/Docker Errors
        ErrorType.CONTAINER_CRASH: "Check container logs, try restarting container",
        ErrorType.CONTAINER_EXITED: "Container exited unexpectedly, check logs and restart",
        ErrorType.DOCKER_ERROR: "Verify Docker is running, check Docker logs",
        
        # System Errors
        ErrorType.SEGMENTATION_FAULT: "Container crashed - restart container or use CPU",
        ErrorType.SYSTEM_KILL: "System resource exhaustion - wait and retry with lower settings",
        ErrorType.OUT_OF_MEMORY: "Reduce memory usage, switch to CPU, or wait for resources",
        ErrorType.TIMEOUT: "Increase timeout, reduce exhaustiveness, or simplify parameters",
        
        # Docking Errors
        ErrorType.RECEPTOR_ERROR: "Validate receptor file format (PDB/PDBQT)",
        ErrorType.LIGAND_ERROR: "Validate ligand file format (SDF/MOL2/SMILES)",
        ErrorType.GRID_ERROR: "Check grid parameters (center, size) are valid",
        ErrorType.INVALID_OUTPUT: "Check docking output, validate input files",
        
        # IO Errors
        ErrorType.FILE_NOT_FOUND: "Check file paths are correct and files exist",
        ErrorType.PERMISSION_DENIED: "Check file/directory permissions",
        ErrorType.DISK_FULL: "Free up disk space before running",
    }
    
    def __init__(self):
        """Initialize failure detector"""
        # Compile regex patterns
        self._compiled_patterns: Dict[ErrorType, List[re.Pattern]] = {}
        for error_type, patterns in self.ERROR_PATTERNS.items():
            self._compiled_patterns[error_type] = [
                re.compile(p, re.IGNORECASE) for p in patterns
            ]
        
        logger.info("FailureDetector initialized")
    
    def analyze_line(self, line: str) -> Optional[Failure]:
        """
        Analyze a single log line for errors.
        
        Args:
            line: Log line to analyze
            
        Returns:
            Failure if error detected, None otherwise
        """
        line = line.strip()
        if not line:
            return None
        
        for error_type, patterns in self._compiled_patterns.items():
            for pattern in patterns:
                if pattern.search(line):
                    return Failure(
                        error_type=error_type,
                        message=line,
                        severity=self.SEVERITY_MAP.get(error_type, "error"),
                        line=line,
                        recoverable=self._is_recoverable(error_type),
                        suggested_fix=self.FIX_SUGGESTIONS.get(error_type, "Check logs for details"),
                    )
        
        return None
    
    def analyze_log(self, log_lines: List[str]) -> List[Failure]:
        """
        Analyze multiple log lines for errors.
        
        Args:
            log_lines: List of log lines
            
        Returns:
            List of failures detected
        """
        failures = []
        
        for line in log_lines:
            failure = self.analyze_line(line)
            if failure:
                failures.append(failure)
        
        # Log summary
        if failures:
            error_types = set(f.error_type for f in failures)
            logger.warning(f"Detected {len(failures)} failures: {[e.value for e in error_types]}")
        
        return failures
    
    def analyze_exit_code(self, exit_code: Optional[int]) -> Optional[Failure]:
        """
        Analyze container exit code.
        
        Args:
            exit_code: Container exit code
            
        Returns:
            Failure if error detected
        """
        if exit_code is None:
            return None
        
        if exit_code == 0:
            return None  # Success
        
        if exit_code == 137:
            # SIGKILL (likely OOM)
            return Failure(
                error_type=ErrorType.SYSTEM_KILL,
                message=f"Container killed with exit code {exit_code}",
                severity="critical",
                recoverable=True,
                suggested_fix=self.FIX_SUGGESTIONS[ErrorType.SYSTEM_KILL],
            )
        
        if exit_code == 139:
            # SIGSEGV (segmentation fault)
            return Failure(
                error_type=ErrorType.SEGMENTATION_FAULT,
                message=f"Segmentation fault (exit code {exit_code})",
                severity="critical",
                recoverable=False,
                suggested_fix=self.FIX_SUGGESTIONS[ErrorType.SEGMENTATION_FAULT],
            )
        
        if exit_code == 1:
            return Failure(
                error_type=ErrorType.DOCKER_ERROR,
                message=f"Container exited with code {exit_code}",
                severity="error",
                recoverable=True,
                suggested_fix="Check container logs for details",
            )
        
        return Failure(
            error_type=ErrorType.UNKNOWN,
            message=f"Unknown exit code: {exit_code}",
            severity="error",
            recoverable=False,
            suggested_fix="Check logs for details",
        )
    
    def _is_recoverable(self, error_type: ErrorType) -> bool:
        """Check if error is recoverable"""
        non_recoverable = [
            ErrorType.SEGMENTATION_FAULT,
            ErrorType.UNKNOWN,
            ErrorType.FILE_NOT_FOUND,
            ErrorType.PERMISSION_DENIED,
        ]
        return error_type not in non_recoverable
    
    def get_error_summary(self, failures: List[Failure]) -> Dict[str, Any]:
        """Get summary of failures"""
        if not failures:
            return {"total": 0, "critical": 0, "errors": 0, "warnings": 0}
        
        critical = sum(1 for f in failures if f.severity == "critical")
        errors = sum(1 for f in failures if f.severity == "error")
        warnings = sum(1 for f in failures if f.severity == "warning")
        
        error_types = {}
        for f in failures:
            etype = f.error_type.value
            error_types[etype] = error_types.get(etype, 0) + 1
        
        return {
            "total": len(failures),
            "critical": critical,
            "errors": errors,
            "warnings": warnings,
            "by_type": error_types,
            "recoverable": sum(1 for f in failures if f.recoverable),
        }
