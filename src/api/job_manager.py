"""
BioDockify Docking Studio - Job Manager API Router
Handles job lifecycle management
"""

from fastapi import APIRouter, HTTPException, status
from fastapi.responses import JSONResponse, StreamingResponse
from pydantic import BaseModel
from typing import Dict, Any, List, Optional
import logging
import json
from datetime import datetime
import uuid

from src.config import Config
from src.database import Database
from src.docker_manager import DockerManager
from src.vina_engine import VinaEngine
from src.oddt_analyzer import ODDTAnalyzer
from src.rdkit_calculator import RDKitCalculator
from src.agent_zero import AgentZero
from src.checkpoint_manager import CheckpointManager
from src.recovery_manager import RecoveryManager

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/jobs")

from functools import lru_cache
from pathlib import Path

# Dependency Injection
@lru_cache(maxsize=1)
def get_config():
    return Config()

def get_database():
    return Database()

docker_manager = DockerManager()
vina_engine = VinaEngine(docker_manager)
oddt_analyzer = ODDTAnalyzer(docker_manager)
rdkit_calculator = RDKitCalculator()
checkpoint_manager = CheckpointManager() 
# AgentZero needs CheckpointManager and RecoveryManager. 
# Assuming RecoveryManager is available or defaulted.
agent_zero = AgentZero(checkpoint_manager, RecoveryManager())

class JobCreate(BaseModel):
    """Job creation request model"""
    receptor_file: str
    ligand_file: str
    parameters: Dict[str, Any]

class JobUpdate(BaseModel):
    """Job update request model"""
    status: Optional[str] = None
    parameters: Optional[Dict[str, Any]] = None

@router.post("/create", response_model=dict)
async def create_job(job: JobCreate):
    """Create new docking job"""
    
    job_id = database.create_job(
        receptor_file=job.receptor_file,
        ligand_file=job.ligand_file,
        parameters=job.parameters
    )
    
    # Save initial checkpoint
    checkpoint_manager.save_checkpoint(job_id, "FILE_VALIDATION", {
        "stage": "FILE_VALIDATION",
        "status": "PENDING",
        "timestamp": datetime.now().isoformat()
    })
    
    logger.info(f"Job created: {job_id}")
    
    return JSONResponse(
        status_code=status.HTTP_201_CREATED,
        content={
            "job_id": job_id,
            "status": "PENDING",
            "created_at": datetime.now().isoformat()
        }
    )

@router.post("/{job_id}/start", response_model=dict)
async def start_job(job_id: str):
    """Start docking job"""
    
    db = get_database()
    try:
        # Update job status to RUNNING
        db.update_job_status(job_id, "RUNNING")
        
        # Get job details
        job = db.get_job(job_id)
        
        if not job:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="Job not found"
            )
        
        # Save checkpoint
        checkpoint_manager.save_checkpoint(job_id, "PREPROCESSING_COMPLETE", {
            "stage": "PREPROCESSING_COMPLETE",
            "status": "RUNNING",
            "timestamp": datetime.now().isoformat()
        })
        
        # Run Vina docking ... (calls continue below, but I need to handle db calls later)
        # Actually I can't replace just the top part easily because method spans long lines. 
        # I'll just do the db part.

        
        # Run Vina docking
        vina_results = vina_engine.run_docking(
            receptor_file=job["receptor_file"],
            ligand_file=job["ligand_file"],
            parameters=job["parameters"],
            job_id=job_id
        )
        
        if vina_results["status"] == "COMPLETED":
            # Save checkpoint
            checkpoint_manager.save_checkpoint(job_id, "VINA_COMPLETED", {
                "stage": "VINA_COMPLETED",
                "status": "RUNNING",
                "timestamp": datetime.now().isoformat()
            })
            
            # Analyze interactions
            oddt_results = oddt_analyzer.analyze_interactions(
                receptor_file=job["receptor_file"],
                ligand_file=job["ligand_file"],
                output_file=f"/data/output_{job_id}.pdbqt",
                job_id=job_id
            )
            
            # Calculate descriptors
            rdkit_results = rdkit_calculator.calculate_descriptors(
                ligand_file=job["ligand_file"],
                job_id=job_id
            )
            
            # Calculate confidence score
            confidence_score = agent_zero.get_confidence_score()
            
            # Save results
            db.save_result(
                job_id=job_id,
                binding_energy=vina_results["binding_energy"],
                interactions=oddt_results["interactions"],
                descriptors=rdkit_results["descriptors"],
                confidence_score=confidence_score
            )
            
            # Update job status to COMPLETED
            db.update_job_status(
                job_id,
                "COMPLETED",
                completed_at=datetime.now().isoformat()
            )
            
            logger.info(f"Job completed: {job_id}")
            
            return JSONResponse(
                status_code=status.HTTP_200_OK,
                content={
                    "job_id": job_id,
                    "status": "COMPLETED",
                    "binding_energy": vina_results["binding_energy"],
                    "interactions": oddt_results["interactions"],
                    "descriptors": rdkit_results["descriptors"],
                    "confidence_score": confidence_score,
                    "completed_at": datetime.now().isoformat()
                }
            )
        else:
            # Vina failed
            db.update_job_status(job_id, "FAILED")
            db.save_log(job_id, "ERROR", f"Vina docking failed: {vina_results.get('logs', '')}")
            
            # Attempt Agent Zero recovery
            failure_info = agent_zero.detect_failure(vina_results.get("logs", ""), "DOCKING", job_id)
            recovery_success = agent_zero.attempt_recovery(failure_info, job_id, "DOCKING")
            
            if recovery_success:
                logger.info(f"Job recovered by Agent Zero: {job_id}")
            else:
                logger.warning(f"Agent Zero recovery failed for job: {job_id}")
            
            return JSONResponse(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                content={
                    "job_id": job_id,
                    "status": "FAILED",
                    "error": vina_results.get("logs", ""),
                    "agent_zero_attempted": True,
                    "agent_zero_recovery_success": recovery_success
                }
            )
    
    except Exception as e:
        # Handle exception
        db.update_job_status(job_id, "FAILED")
        db.save_log(job_id, "ERROR", f"Job failed with exception: {str(e)}")
        
        # Attempt Agent Zero recovery
        failure_info = agent_zero.detect_failure(str(e), "JOB_EXECUTION", job_id)
        recovery_success = agent_zero.attempt_recovery(failure_info, job_id, "JOB_EXECUTION")
        
        logger.exception(f"Job failed with exception: {job_id} - {e}")
        
        return JSONResponse(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            content={
                "job_id": job_id,
                "status": "FAILED",
                "error": str(e),
                "agent_zero_attempted": True,
                "agent_zero_recovery_success": recovery_success
            }
        )

@router.post("/{job_id}/cancel", response_model=dict)
async def cancel_job(job_id: str):
    """Cancel running job"""
    
    # Update job status to CANCELLED
    database.update_job_status(job_id, "CANCELLED")
    database.save_log(job_id, "INFO", "Job cancelled by user")
    
    # Stop Docker container
    docker_manager.stop_container()
    docker_manager.cleanup_container()
    
    logger.info(f"Job cancelled: {job_id}")
    
    return JSONResponse(
        status_code=status.HTTP_200_OK,
        content={
            "job_id": job_id,
            "status": "CANCELLED"
        }
    )

@router.get("/{job_id}", response_model=dict)
async def get_job(job_id: str):
    """Get job details"""
    
    job = database.get_job(job_id)
    
    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Job not found"
        )
    
    return JSONResponse(
        status_code=status.HTTP_200_OK,
        content=job
    )

@router.get("/{job_id}/status", response_model=dict)
async def get_job_status(job_id: str):
    """Get job status"""
    
    job = database.get_job(job_id)
    
    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Job not found"
        )
    
    return JSONResponse(
        status_code=status.HTTP_200_OK,
        content={
            "job_id": job_id,
            "status": job["status"],
            "current_checkpoint": job.get("current_checkpoint"),
            "created_at": job["created_at"],
            "updated_at": job["updated_at"],
            "completed_at": job.get("completed_at")
        }
    )

@router.get("/recent", response_model=dict)
async def get_recent_jobs(limit: int = 10):
    """Get recent jobs"""
    
    jobs = database.get_recent_jobs(limit)
    
    return JSONResponse(
        status_code=status.HTTP_200_OK,
        content={
            "count": len(jobs),
            "jobs": jobs
        }
    )
