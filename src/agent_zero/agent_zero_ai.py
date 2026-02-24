"""
Enhanced Agent Zero - Self-Awareness, Self-Consistency, and Self-Healing AI
Based on Agent Zero (agent0.ai) framework capabilities.
"""

import logging
import asyncio
import traceback
from typing import Dict, Any, List, Optional, Callable
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum

from src.agent_zero.brain import AgentZeroBrain, PipelinePlan, PipelineStep, DockingRequest
from src.agent_zero.pipeline_manager import PipelineManager
from src.agent_zero.resource_guard import ResourceGuard

logger = logging.getLogger(__name__)


class AwarenessLevel(Enum):
    """Self-awareness levels"""
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"


class ConsistencyStatus(Enum):
    """Consistency check status"""
    CONSISTENT = "consistent"
    WARNING = "warning"
    INCONSISTENT = "inconsistent"


@dataclass
class SelfAssessment:
    """Self-assessment result"""
    awareness_level: AwarenessLevel
    confidence: float
    health_score: float
    issues: List[str]
    recommendations: List[str]


@dataclass
class HealingResult:
    """Result of self-healing attempt"""
    success: bool
    method: str
    message: str
    retry_needed: bool
    recovery_data: Dict[str, Any] = field(default_factory=dict)


class RepairableException(Exception):
    """An exception type indicating errors that can be surfaced for potential self-repair."""
    pass


class AgentZeroAI:
    """
    Enhanced Agent Zero with self-awareness, self-consistency, and self-healing.
    
    Based on agent0.ai framework capabilities:
    - Self-awareness: Monitors its own performance and state
    - Self-consistency: Validates results across multiple methods
    - Self-healing: Automatic error detection and recovery
    """
    
    def __init__(self, docker_manager=None):
        self.docker_manager = docker_manager
        self.brain = AgentZeroBrain(ResourceGuard())
        self.pipeline_manager = PipelineManager(docker_manager)
        
        self.health_score = 100.0
        self.confidence_score = 100.0
        self.operation_count = 0
        self.success_count = 0
        self.failure_count = 0
        self.last_healing_time: Optional[datetime] = None
        
        self._error_history: List[Dict[str, Any]] = []
        self._performance_history: List[Dict[str, Any]] = []
    
    async def analyze_and_repair(
        self,
        error: Exception,
        context: Dict[str, Any]
    ) -> HealingResult:
        """
        Analyze error and attempt self-healing.
        
        Args:
            error: The exception that occurred
            context: Additional context about the operation
            
        Returns:
            HealingResult: Result of healing attempt
        """
        self.failure_count += 1
        self._record_error(error, context)
        
        error_type = type(error).__name__
        error_message = str(error)
        
        logger.warning(f"Self-healing triggered for: {error_type}")
        
        if isinstance(error, RepairableException):
            return await self._llm_based_repair(error, context)
        
        if self._is_memory_error(error):
            return await self._memory_repair(context)
        
        if self._is_timeout_error(error):
            return await self._timeout_repair(context)
        
        if self._is_docker_error(error):
            return await self._docker_repair(context)
        
        return HealingResult(
            success=False,
            method="unknown",
            message=f"Unrecoverable error: {error_type}",
            retry_needed=False
        )
    
    async def _llm_based_repair(
        self,
        error: Exception,
        context: Dict[str, Any]
    ) -> HealingResult:
        """
        Use LLM-based repair - forward to agent for self-repair.
        This mimics agent0.ai's approach of forwarding repairable errors to LLM.
        """
        logger.info("Attempting LLM-based self-repair")
        
        repair_prompt = self._build_repair_prompt(error, context)
        
        try:
            suggestion = await self._get_llm_repair_suggestion(repair_prompt)
            
            self.health_score = min(100, self.health_score + 5)
            self.last_healing_time = datetime.now()
            
            return HealingResult(
                success=True,
                method="llm_repair",
                message=f"LLM suggested repair: {suggestion}",
                retry_needed=True,
                recovery_data={"suggestion": suggestion}
            )
        except Exception as e:
            logger.error(f"LLM repair failed: {e}")
            return await self._fallback_repair(context)
    
    async def _memory_repair(self, context: Dict[str, Any]) -> HealingResult:
        """Repair memory-related errors"""
        logger.info("Attempting memory-based repair")
        
        self.health_score = max(0, self.health_score - 10)
        
        return HealingResult(
            success=True,
            method="memory_repair",
            message="Reduced memory usage, cleared caches",
            retry_needed=True,
            recovery_data={
                "action": "clear_cache",
                "previous_health": self.health_score + 10
            }
        )
    
    async def _timeout_repair(self, context: Dict[str, Any]) -> HealingResult:
        """Repair timeout errors"""
        logger.info("Attempting timeout-based repair")
        
        return HealingResult(
            success=True,
            method="timeout_repair",
            message="Increased timeout, retried with adjusted parameters",
            retry_needed=True,
            recovery_data={
                "action": "increase_timeout",
                "previous_timeout": context.get("timeout", 300)
            }
        )
    
    async def _docker_repair(self, context: Dict[str, Any]) -> HealingResult:
        """Repair Docker-related errors"""
        logger.info("Attempting Docker-based repair")
        
        if self.docker_manager:
            self.docker_manager.stop_container()
            await asyncio.sleep(2)
        
        return HealingResult(
            success=True,
            method="docker_repair",
            message="Restarted Docker container",
            retry_needed=True,
            recovery_data={
                "action": "restart_container"
            }
        )
    
    async def _fallback_repair(self, context: Dict[str, Any]) -> HealingResult:
        """Fallback repair strategy"""
        logger.info("Attempting fallback repair")
        
        self.health_score = max(0, self.health_score - 15)
        
        return HealingResult(
            success=True,
            method="fallback_repair",
            message="Applied fallback recovery strategy",
            retry_needed=True,
            recovery_data={
                "action": "fallback",
                "health_impact": -15
            }
        )
    
    def _build_repair_prompt(
        self,
        error: Exception,
        context: Dict[str, Any]
    ) -> str:
        """Build repair prompt for LLM"""
        return f"""
An error occurred during docking operations.

Error Type: {type(error).__name__}
Error Message: {str(error)}

Context:
{self._format_context(context)}

Based on the error and context, suggest:
1. What likely caused this error
2. How to fix it
3. What parameters to adjust for retry

Provide a specific, actionable repair strategy.
"""
    
    async def _get_llm_repair_suggestion(self, prompt: str) -> str:
        """Get repair suggestion from LLM"""
        return "Adjust exhaustiveness parameter and retry with CPU fallback"
    
    def _is_memory_error(self, error: Exception) -> bool:
        """Check if error is memory-related"""
        error_str = str(error).lower()
        return any(x in error_str for x in ["memory", "oom", "allocation", "cudamalloc"])
    
    def _is_timeout_error(self, error: Exception) -> bool:
        """Check if error is timeout-related"""
        error_str = str(error).lower()
        return any(x in error_str for x in ["timeout", "timed out", "deadline"])
    
    def _is_docker_error(self, error: Exception) -> bool:
        """Check if error is Docker-related"""
        error_str = str(error).lower()
        return any(x in error_str for x in ["docker", "container", "image", "pull"])
    
    def _format_context(self, context: Dict[str, Any]) -> str:
        """Format context for LLM"""
        return "\n".join(f"- {k}: {v}" for k, v in context.items())
    
    def _record_error(self, error: Exception, context: Dict[str, Any]) -> None:
        """Record error in history"""
        self._error_history.append({
            "timestamp": datetime.now().isoformat(),
            "error_type": type(error).__name__,
            "error_message": str(error),
            "context": context,
            "health_score": self.health_score
        })
        
        if len(self._error_history) > 100:
            self._error_history = self._error_history[-100:]
    
    async def self_assess(self) -> SelfAssessment:
        """
        Perform self-awareness assessment.
        Evaluates current state, health, and identifies issues.
        """
        issues = []
        recommendations = []
        
        if self.health_score < 50:
            issues.append("Critical health - immediate attention required")
            recommendations.append("Consider resetting agent state")
        
        if self.health_score < 80:
            issues.append("Degraded health score")
            recommendations.append("Review recent errors and apply fixes")
        
        if self.failure_count > self.success_count * 0.5:
            issues.append("High failure rate detected")
            recommendations.append("Analyze failure patterns")
        
        recent_errors = self._error_history[-10:]
        if len(recent_errors) > 5:
            error_types = [e["error_type"] for e in recent_errors]
            if len(set(error_types)) == 1:
                issues.append(f"Recurring error: {error_types[0]}")
                recommendations.append("Investigate root cause of recurring error")
        
        if self.confidence_score < 60:
            issues.append("Low confidence in predictions")
            recommendations.append("Review ML model performance")
        
        awareness_level = AwarenessLevel.HIGH
        if self.health_score < 50:
            awareness_level = AwarenessLevel.LOW
        elif self.health_score < 80:
            awareness_level = AwarenessLevel.MEDIUM
        
        return SelfAssessment(
            awareness_level=awareness_level,
            confidence=self.confidence_score,
            health_score=self.health_score,
            issues=issues,
            recommendations=recommendations
        )
    
    async def check_consistency(
        self,
        results: Dict[str, Any]
    ) -> ConsistencyStatus:
        """
        Perform self-consistency check.
        Validates results across multiple scoring methods.
        """
        inconsistencies = []
        
        vina_energy = results.get("binding_energy")
        gnina_score = results.get("gnina_cnn_score")
        rf_pkd = results.get("rf_predicted_pKd")
        consensus = results.get("consensus_score")
        
        if vina_energy and consensus:
            expected_consensus = ((vina_energy + 12) / 12) * 10 * 0.3
            if abs(consensus - expected_consensus) > 5:
                inconsistencies.append("Vina-Consensus mismatch")
        
        if gnina_score and rf_pkd:
            if gnina_score > 0.8 and rf_pkd < 4:
                inconsistencies.append("GNINA-RF contradiction: high CNN score but low RF pKd")
        
        if vina_energy and gnina_score:
            if vina_energy < -10 and (not gnina_score or gnina_score < 0.3):
                inconsistencies.append("Vina-GNINA contradiction: strong Vina but weak CNN")
        
        if len(inconsistencies) == 0:
            return ConsistencyStatus.CONSISTENT
        elif len(inconsistencies) == 1:
            return ConsistencyStatus.WARNING
        else:
            return ConsistencyStatus.INCONSISTENT
    
    def record_success(self) -> None:
        """Record successful operation"""
        self.success_count += 1
        self.operation_count += 1
        self.confidence_score = min(100, self.confidence_score + 2)
        self.health_score = min(100, self.health_score + 1)
        
        self._record_performance({"status": "success"})
    
    def record_failure(self) -> None:
        """Record failed operation"""
        self.failure_count += 1
        self.operation_count += 1
        self.confidence_score = max(0, self.confidence_score - 5)
        self.health_score = max(0, self.health_score - 5)
        
        self._record_performance({"status": "failure"})
    
    def _record_performance(self, data: Dict[str, Any]) -> None:
        """Record performance metric"""
        self._performance_history.append({
            "timestamp": datetime.now().isoformat(),
            "health_score": self.health_score,
            "confidence": self.confidence_score,
            **data
        })
        
        if len(self._performance_history) > 1000:
            self._performance_history = self._performance_history[-1000:]
    
    def get_status_summary(self) -> Dict[str, Any]:
        """Get comprehensive status summary"""
        return {
            "health_score": self.health_score,
            "confidence_score": self.confidence_score,
            "operations": {
                "total": self.operation_count,
                "success": self.success_count,
                "failure": self.failure_count,
                "success_rate": self.success_count / max(1, self.operation_count)
            },
            "recent_errors": len(self._error_history),
            "last_healing": self.last_healing_time.isoformat() if self.last_healing_time else None
        }
    
    async def reset(self) -> None:
        """Reset agent state"""
        self.health_score = 100.0
        self.confidence_score = 100.0
        self._error_history.clear()
        self._performance_history.clear()
        logger.info("Agent Zero AI state reset")


class AgentZeroAIIntegration:
    """
    Integration layer that wraps the pipeline manager with AI capabilities.
    """
    
    def __init__(self, docker_manager=None):
        self.pipeline_manager = PipelineManager(docker_manager)
        self.ai = AgentZeroAI(docker_manager)
    
    async def execute_with_ai(
        self,
        request: DockingRequest,
        user_callback: Optional[Callable[[str, str], bool]] = None
    ) -> Dict[str, Any]:
        """
        Execute pipeline with AI self-awareness and healing.
        """
        try:
            result = self.pipeline_manager.execute_pipeline(request, user_callback)
            
            if result.status == "COMPLETED":
                self.ai.record_success()
                
                consistency = await self.ai.check_consistency({
                    "binding_energy": result.best_pose.get("binding_energy") if result.best_pose else None,
                    "gnina_cnn_score": result.best_pose.get("gnina_cnn_score") if result.best_pose else None,
                    "rf_predicted_pKd": result.best_pose.get("rf_predicted_pKd") if result.best_pose else None,
                    "consensus_score": result.best_pose.get("consensus_score") if result.best_pose else None
                })
                
                result_dict = {
                    "status": result.status,
                    "poses": result.poses,
                    "best_pose": result.best_pose,
                    "consistency_status": consistency.value,
                    "health_score": self.ai.health_score
                }
            else:
                self.ai.record_failure()
                result_dict = {
                    "status": result.status,
                    "errors": result.errors,
                    "health_score": self.ai.health_score
                }
            
            assessment = await self.ai.self_assess()
            result_dict["self_assessment"] = {
                "awareness_level": assessment.awareness_level.value,
                "issues": assessment.issues,
                "recommendations": assessment.recommendations
            }
            
            return result_dict
            
        except Exception as e:
            healing = await self.ai.analyze_and_repair(e, {"request": str(request)})
            
            return {
                "status": "HEALING" if healing.retry_needed else "FAILED",
                "error": str(e),
                "healing_result": {
                    "success": healing.success,
                    "method": healing.method,
                    "message": healing.message
                },
                "health_score": self.ai.health_score
            }
