"""
BioDockify Docking Studio - Agent Zero API Router
Handles Agent Zero failure detection and recovery
"""

from fastapi import APIRouter, HTTPException, status
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from typing import Dict, Any, List, Optional
import logging
from datetime import datetime

from src.agent_zero import AgentZero, FailureType, FailureInfo
from src.checkpoint_manager import CheckpointManager
from src.recovery_manager import RecoveryManager

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/agent-zero")

checkpoint_manager = CheckpointManager()
recovery_manager = RecoveryManager()
agent_zero = AgentZero(checkpoint_manager, recovery_manager)

class FailureReport(BaseModel):
    """Failure report model"""
    failure_type: str
    severity: str
    message: str
    diagnosis: str
    recovery_action: str
    stage: str
    occurred_at: str

class ConfidenceScore(BaseModel):
    """Confidence score model"""
    job_id: str
    score: int
    level: str

class RecoveryReport(BaseModel):
    """Recovery report model"""
    job_id: str
    recovery_action: str
    success: bool
    timestamp: str

@router.post("/detect-failure", response_model=FailureReport)
async def detect_failure(error_message: str, stage: str, job_id: str):
    """Detect failure from error message"""
    
    failure_info = agent_zero.detect_failure(error_message, stage, job_id)
    
    return FailureReport(
        failure_type=failure_info.failure_type.value,
        severity=failure_info.severity.name,
        message=failure_info.message,
        diagnosis=failure_info.diagnosis,
        recovery_action=failure_info.recovery_action,
        stage=failure_info.stage,
        occurred_at=failure_info.occurred_at.isoformat()
    )

@router.post("/attempt-recovery/{job_id}", response_model=RecoveryReport)
async def attempt_recovery(job_id: str, current_stage: str, error_message: str):
    """Attempt recovery for failure"""
    
    failure_info = agent_zero.detect_failure(error_message, current_stage, job_id)
    recovery_success = agent_zero.attempt_recovery(failure_info, job_id, current_stage)
    
    return RecoveryReport(
        job_id=job_id,
        recovery_action=failure_info.recovery_action,
        success=recovery_success,
        timestamp=datetime.now().isoformat()
    )

@router.get("/confidence/{job_id}", response_model=ConfidenceScore)
async def get_confidence_score(job_id: str):
    """Get confidence score for job"""
    
    score = agent_zero.get_confidence_score()
    
    if score >= 80:
        level = "HIGH"
    elif score >= 60:
        level = "MEDIUM"
    else:
        level = "LOW"
    
    return ConfidenceScore(
        job_id=job_id,
        score=score,
        level=level
    )

@router.post("/reset-confidence", response_model=dict)
async def reset_confidence_score():
    """Reset confidence score to 100"""
    
    agent_zero.reset_confidence()
    
    return JSONResponse(
        status_code=status.HTTP_200_OK,
        content={
            "status": "confidence_reset",
            "score": 100,
            "level": "HIGH"
        }
    )

@router.get("/safe-checkpoints", response_model=dict)
async def get_safe_checkpoints():
    """Get list of safe checkpoints"""
    
    return JSONResponse(
        status_code=status.HTTP_200_OK,
        content={
            "safe_checkpoints": CheckpointManager.SAFE_CHECKPOINTS
        }
    )
