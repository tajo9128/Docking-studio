"""
Agent Zero v4 - Internet Diagnostics Module
Safe, advisory-only internet access for troubleshooting.
NEVER executes remote code - read-only advisory only.
"""

import logging
import requests
from typing import Optional, Dict, Any, List
from dataclasses import dataclass
from enum import Enum
import re

logger = logging.getLogger(__name__)


class DiagnosticsMode(Enum):
    """Internet diagnostics mode"""
    DISABLED = "disabled"
    ENABLED = "enabled"


@dataclass
class AdvisoryResult:
    """Result from advisory lookup"""
    found: bool
    topic: str
    suggestions: List[str]
    source: str
    error: Optional[str]


class InternetDiagnostics:
    """
    Advisory-only internet diagnostics.
    Fetches documentation and troubleshooting info - never executes.
    """
    
    # Knowledge base of common errors and solutions (local fallback)
    LOCAL_KB = {
        "gnina cuda oom": [
            "Reduce CNN batch size: --cnn_batch_size 4",
            "Lower exhaustiveness: --exhaustiveness 8",
            "Switch to CPU mode: use vina-cpu engine",
            "Clear GPU memory: restart Docker container"
        ],
        "vina gpu error": [
            "Check Docker GPU access: docker run --gpus all nvidia/cuda:12.0-base nvidia-smi",
            "Verify NVIDIA driver: nvidia-smi",
            "Try CPU fallback: use vina-cpu engine"
        ],
        "container crash": [
            "Check Docker is running: docker info",
            "Verify image exists: docker images",
            "Check logs: docker logs <container_id>",
            "Try pulling latest image"
        ],
        "docker gpu access": [
            "Install nvidia-docker2: distribution=$(. /etc/os-release;echo $ID$VERSION_ID)",
            "Update Docker config: sudo nvidia-ctk runtime configure --runtime=docker",
            "Restart Docker: sudo systemctl restart docker",
            "Test: docker run --rm --gpus all nvidia/cuda:12.0-base nvidia-smi"
        ],
        "out of memory": [
            "Close other GPU applications",
            "Reduce batch size",
            "Use CPU engine instead",
            "Wait for memory to free up"
        ],
        "segmentation fault": [
            "Container may have crashed - check Docker logs",
            "Try running with --privileged flag",
            "Check system memory",
            "Try CPU fallback"
        ]
    }
    
    def __init__(self, mode: DiagnosticsMode = DiagnosticsMode.DISABLED):
        self.mode = mode
        self._session: Optional[requests.Session] = None
    
    def set_mode(self, mode: DiagnosticsMode):
        """Enable or disable internet diagnostics"""
        self.mode = mode
        logger.info(f"Internet diagnostics mode: {mode.value}")
    
    def lookup(self, error_message: str) -> AdvisoryResult:
        """
        Look up troubleshooting for an error.
        
        Args:
            error_message: Error message to look up
            
        Returns:
            AdvisoryResult with suggestions
        """
        # First try local knowledge base
        local_result = self._lookup_local(error_message)
        if local_result.found:
            return local_result
        
        # If enabled, try internet lookup
        if self.mode == DiagnosticsMode.ENABLED:
            return self._lookup_internet(error_message)
        
        return AdvisoryResult(
            found=False,
            topic=error_message[:50],
            suggestions=[],
            source="local",
            error="Diagnostics disabled"
        )
    
    def _lookup_local(self, error_message: str) -> AdvisoryResult:
        """Search local knowledge base"""
        error_lower = error_message.lower()
        
        for kb_key, solutions in self.LOCAL_KB.items():
            if kb_key in error_lower:
                return AdvisoryResult(
                    found=True,
                    topic=kb_key,
                    suggestions=solutions,
                    source="local_kb",
                    error=None
                )
        
        return AdvisoryResult(
            found=False,
            topic=error_message[:50],
            suggestions=[],
            source="local",
            error=None
        )
    
    def _lookup_internet(self, error_message: str) -> AdvisoryResult:
        """
        Advisory-only internet lookup.
        Fetches search results - does NOT execute code.
        """
        # This is a safe implementation - it would fetch from a known
        # documentation source, but for security we keep it local-only
        # by default. Internet access is read-only advisory.
        
        logger.info(f"Internet diagnostics lookup for: {error_message[:50]}")
        
        # For security, we don't implement actual HTTP fetches
        # Instead, we return local suggestions
        return AdvisoryResult(
            found=False,
            topic=error_message[:50],
            suggestions=[
                "Enable internet diagnostics in settings for online lookup",
                "Check documentation at: https://github.com/biodock/docking-studio",
                "Search common error: " + error_message[:30]
            ],
            source="internet",
            error="Internet lookup not implemented (security)"
        )
    
    def get_suggested_actions(self, failure_type: str) -> List[str]:
        """
        Get predefined safe actions based on failure type.
        
        Args:
            failure_type: Type of failure (gpu_oom, cuda_error, etc.)
            
        Returns:
            List of safe action descriptions
        """
        suggestions = {
            "gpu_oom": [
                "Reduce CNN batch size",
                "Switch to CPU engine",
                "Wait for VRAM to free",
                "Close other GPU applications"
            ],
            "cuda_error": [
                "Restart Docker container",
                "Check NVIDIA driver",
                "Switch to CPU engine",
                "Verify CUDA installation"
            ],
            "container_crash": [
                "Check Docker logs",
                "Restart Docker daemon",
                "Verify image is up to date",
                "Try with --privileged flag"
            ],
            "timeout": [
                "Reduce exhaustiveness",
                "Reduce number of poses",
                "Use faster engine"
            ],
            "segfault": [
                "Check system memory",
                "Try CPU engine",
                "Restart container"
            ]
        }
        
        return suggestions.get(failure_type, [
            "Review error logs",
            "Try restarting the job",
            "Switch compute mode"
        ])


# Singleton
_diagnostics: Optional[InternetDiagnostics] = None

def get_internet_diagnostics() -> InternetDiagnostics:
    """Get singleton diagnostics instance"""
    global _diagnostics
    if _diagnostics is None:
        _diagnostics = InternetDiagnostics()
    return _diagnostics
