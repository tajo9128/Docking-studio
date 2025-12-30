"""
BioDockify Docking Studio - Checkpoint API Router
Handles checkpoint management for job resumption
"""

from fastapi import APIRouter, HTTPException, status
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from typing import Dict, Any, Optional
import logging
from datetime import datetime

from src.checkpoint_manager import CheckpointManager
from src.database import Database

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/checkpoints")

checkpoint_manager = CheckpointManager()

class CheckpointSave(BaseModel):
    """Checkpoint save request model"""
    job_id: str
    stage: str
    data: Dict[str, Any]

class CheckpointGet(BaseModel):
    """Checkpoint get request model"""
    job_id: str
    stage: str

@router.post("/save", response_model=dict)
async def save_checkpoint(checkpoint: CheckpointSave):
    """Save checkpoint"""
    
    success = checkpoint_manager.save_checkpoint(
        job_id=checkpoint.job_id,
        stage=checkpoint.stage,
        data=checkpoint.data
    )
    
    if success:
        return JSONResponse(
            status_code=status.HTTP_200_OK,
            content={
                "job_id": checkpoint.job_id,
                "stage": checkpoint.stage,
                "timestamp": datetime.now().isoformat(),
                "status": "saved"
            }
        )
    else:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to save checkpoint"
        )

@router.get("/{job_id}/{stage}", response_model=dict)
async def get_checkpoint(job_id: str, stage: str):
    """Get checkpoint by job ID and stage"""
    
    checkpoint = checkpoint_manager.get_checkpoint(job_id, stage)
    
    if not checkpoint:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Checkpoint not found"
        )
    
    return JSONResponse(
        status_code=status.HTTP_200_OK,
        content=checkpoint
    )

@router.get("/{job_id}/previous/{stage}", response_model=dict)
async def get_previous_checkpoint(job_id: str, stage: str):
    """Get previous safe checkpoint"""
    
    previous_checkpoint = checkpoint_manager.get_previous_checkpoint(job_id, stage)
    
    if not previous_checkpoint:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="No previous checkpoint found"
        )
    
    return JSONResponse(
        status_code=status.HTTP_200_OK,
        content=previous_checkpoint
    )

@router.delete("/{job_id}", response_model=dict)
async def clear_checkpoints(job_id: str):
    """Clear all checkpoints for job"""
    
    checkpoint_manager.clear_checkpoints(job_id)
    
    return JSONResponse(
        status_code=status.HTTP_200_OK,
        content={
            "job_id": job_id,
            "status": "checkpoints_cleared",
            "timestamp": datetime.now().isoformat()
        }
    )

@router.get("/safe-stages", response_model=dict)
async def get_safe_stages():
    """Get list of safe checkpoint stages"""
    
    return JSONResponse(
        status_code=status.HTTP_200_OK,
        content={
            "safe_stages": CheckpointManager.SAFE_CHECKPOINTS
        }
    )
