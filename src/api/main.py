"""
BioDockify Docking Studio - FastAPI Main Application
Handles HTTP API requests for docking operations
"""

from fastapi import FastAPI, HTTPException, status, Request
from fastapi.responses import JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import logging
from typing import Optional, Dict, Any, List
from datetime import datetime
import uuid

from src.database import Database
from src.docker_manager import DockerManager
from src.vina_engine import VinaEngine
from src.oddt_analyzer import ODDTAnalyzer
from src.rdkit_calculator import RDKitCalculator
from src.agent_zero import AgentZero
from src.recovery_manager import RecoveryManager
from src.checkpoint_manager import CheckpointManager
from src.models.docking_job import DockingJob as JobModel
from src.models.result import Result as ResultModel
from src.config import Config
from src.utils.log_utils import setup_logging

# Setup logging
logger = setup_logging()

# Create FastAPI application
app = FastAPI(
    title="BioDockify API",
    description="API for molecular docking operations",
    version="1.0.0"
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Allow all origins for local use
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Initialize components
config = Config()
database = Database()
docker_manager = DockerManager()
vina_engine = VinaEngine(docker_manager)
oddt_analyzer = ODDTAnalyzer(docker_manager)
rdkit_calculator = RDKitCalculator()
checkpoint_manager = CheckpointManager()
recovery_manager = RecoveryManager()
agent_zero = AgentZero(checkpoint_manager, recovery_manager)

# Pydantic models
class HealthCheck(BaseModel):
    """Health check response model"""
    status: str
    docker_available: bool
    docker_running: bool
    timestamp: str

class DockingJobCreate(BaseModel):
    """Create docking job request model"""
    receptor_file: str
    ligand_file: str
    parameters: Dict[str, Any]

class DockingJobStatus(BaseModel):
    """Docking job status response model"""
    job_id: str
    status: str
    created_at: str
    updated_at: str
    completed_at: Optional[str] = None
    progress_percentage: Optional[int] = 0
    current_stage: Optional[str] = None
    error_message: Optional[str] = None

class DockingJobCancel(BaseModel):
    """Cancel docking job request model"""
    job_id: str

class AgentZeroStatus(BaseModel):
    """Agent Zero status response model"""
    job_id: str
    confidence_score: int
    failures_detected: Optional[Dict[str, Any]] = {}
    repairs_attempted: Optional[Dict[str, Any]] = {}

class ResultExport(BaseModel):
    """Export results request model"""
    job_id: str
    format: str  # "csv", "json", "pymol"

# Health check endpoint
@app.get("/health", response_model=HealthCheck)
async def health_check():
    """Check API health and Docker status"""
    docker_available = docker_manager.is_docker_available()
    docker_running = docker_manager.is_docker_running()
    
    return HealthCheck(
        status="healthy" if docker_available else "unhealthy",
        docker_available=docker_available,
        docker_running=docker_running,
        timestamp=datetime.now().isoformat()
    )

# Create docking job endpoint
@app.post("/api/jobs", response_model=DockingJobStatus)
async def create_docking_job(job: DockingJobCreate):
    """Create new docking job"""
    job_id = database.create_job(
        receptor_file=job.receptor_file,
        ligand_file=job.ligand_file,
        parameters=job.parameters
    )
    
    logger.info(f"Docking job created: {job_id}")
    
    return DockingJobStatus(
        job_id=job_id,
        status="PENDING",
        created_at=datetime.now().isoformat(),
        updated_at=datetime.now().isoformat()
    )

# Start docking job endpoint
@app.post("/api/jobs/{job_id}/start", response_model=DockingJobStatus)
async def start_docking_job(job_id: str):
    """Start docking job"""
    
    # Update job status to RUNNING
    database.update_job_status(job_id, "RUNNING")
    
    try:
        # Start Docker container
        success = docker_manager.start_container(
            receptor_file=job.receptor_file,
            ligand_file=job.ligand_file,
            parameters=job.parameters,
            job_id=job_id
        )
        
        if not success:
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail="Failed to start Docker container"
            )
        
        # Run Vina docking
        vina_results = vina_engine.run_docking(
            receptor_file=job.receptor_file,
            ligand_file=job.ligand_file,
            parameters=job.parameters,
            job_id=job_id
        )
        
        # Analyze interactions with ODDT
        oddt_results = oddt_analyzer.analyze_interactions(
            receptor_file=job.receptor_file,
            ligand_file=job.ligand_file,
            output_file=f"/data/output_{job_id}.pdbqt",
            job_id=job_id
        )
        
        # Calculate descriptors with RDKit
        rdkit_results = rdkit_calculator.calculate_descriptors(
            ligand_file=job.ligand_file,
            job_id=job_id
        )
        
        # Calculate confidence score
        confidence_score = agent_zero.get_confidence_score()
        
        # Save results
        success = database.save_result(
            job_id=job_id,
            binding_energy=vina_results["binding_energy"],
            interactions=oddt_results["interactions"],
            descriptors=rdkit_results["descriptors"],
            confidence_score=confidence_score
        )
        
        # Update job status to COMPLETED
        database.update_job_status(
            job_id,
            "COMPLETED",
            completed_at=datetime.now().isoformat()
        )
        
        logger.info(f"Docking job completed: {job_id}")
        
        return DockingJobStatus(
            job_id=job_id,
            status="COMPLETED",
            created_at=datetime.now().isoformat(),
            updated_at=datetime.now().isoformat(),
            completed_at=datetime.now().isoformat(),
            progress_percentage=100,
            current_stage="COMPLETED"
        )
    
    except HTTPException:
        # Re-raise HTTP exceptions
        raise
    
    except Exception as e:
        # Log exception
        logger.exception(f"Failed to start docking job {job_id}: {e}")
        
        # Update job status to FAILED
        database.update_job_status(job_id, "FAILED")
        database.save_log(job_id, "CRITICAL", f"Docking failed: {str(e)}")
        
        # Attempt Agent Zero recovery
        failure_info = agent_zero.detect_failure(str(e), "DOCKING", job_id)
        recovery_success = agent_zero.attempt_recovery(failure_info, job_id, "DOCKING")
        
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Docking failed: {str(e)}"
        )

# Get docking job status endpoint
@app.get("/api/jobs/{job_id}", response_model=DockingJobStatus)
async def get_docking_job_status(job_id: str):
    """Get docking job status"""
    
    job = database.get_job(job_id)
    
    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Job not found"
        )
    
    return DockingJobStatus(
        job_id=job["id"],
        status=job["status"],
        created_at=job["created_at"],
        updated_at=job["updated_at"],
        completed_at=job.get("completed_at"),
        progress_percentage=0,  # Would calculate from database
        current_stage=job["status"],
        error_message=None
    )

# Cancel docking job endpoint
@app.post("/api/jobs/{job_id}/cancel", response_model=DockingJobStatus)
async def cancel_docking_job(job_id: str, cancel_request: DockingJobCancel):
    """Cancel docking job"""
    
    # Update job status to CANCELLED
    success = database.update_job_status(job_id, "CANCELLED")
    
    if not success:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to cancel job"
        )
    
    # Stop Docker container
    docker_manager.stop_container()
    
    # Clean up container
    docker_manager.cleanup_container()
    
    logger.info(f"Docking job cancelled: {job_id}")
    
    return DockingJobStatus(
        job_id=job_id,
        status="CANCELLED",
        created_at=datetime.now().isoformat(),
        updated_at=datetime.now().isoformat(),
        completed_at=datetime.now().isoformat()
    )

# Get Agent Zero status endpoint
@app.get("/api/jobs/{job_id}/agent-zero", response_model=AgentZeroStatus)
async def get_agent_zero_status(job_id: str):
    """Get Agent Zero status for job"""
    
    confidence_score = agent_zero.get_confidence_score()
    
    # In real implementation, would retrieve failures and repairs from database
    failures_detected = {}  # database.get_job_failures(job_id)
    repairs_attempted = {}  # database.get_job_repairs(job_id)
    
    return AgentZeroStatus(
        job_id=job_id,
        confidence_score=confidence_score,
        failures_detected=failures_detected,
        repairs_attempted=repairs_attempted
    )

# Export results endpoint
@app.get("/api/jobs/{job_id}/export", response_model=Dict[str, Any])
async def export_results(job_id: str, export_format: str):
    """Export docking results"""
    
    # Get results from database
    job = database.get_job(job_id)
    
    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Job not found"
        )
    
    # In real implementation, would export results in requested format
    if export_format == "csv":
        content = "Binding Energy,Interactions,Descriptors\n"
        content += f"{job.get('binding_energy', '')},Mock Interactions,Mock Descriptors\n"
    elif export_format == "json":
        content = {
            "job_id": job_id,
            "binding_energy": job.get("binding_energy"),
            "interactions": job.get("interactions"),
            "descriptors": job.get("descriptors"),
            "confidence_score": agent_zero.get_confidence_score()
        }
    elif export_format == "pymol":
        content = "load your_ligand.pdb\nzoom 100\n"
    else:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Invalid export format"
        )
    
    logger.info(f"Results exported for job {job_id} as {export_format}")
    
    return JSONResponse(
        status="success",
        format=export_format,
        content=content,
        timestamp=datetime.now().isoformat()
    )

# Get recent jobs endpoint
@app.get("/api/jobs/recent", response_model=List[DockingJobStatus])
async def get_recent_jobs(limit: int = 10):
    """Get recent docking jobs"""
    
    jobs = database.get_recent_jobs(limit)
    
    job_statuses = [
        DockingJobStatus(
            job_id=job["id"],
            status=job["status"],
            created_at=job["created_at"],
            updated_at=job["updated_at"],
            completed_at=job.get("completed_at"),
            progress_percentage=0,
            current_stage=job["status"],
            error_message=None
        )
        for job in jobs
    ]
    
    return job_statuses

# Exception handler
@app.exception_handler(Exception)
async def global_exception_handler(request: Request, exc: Exception):
    """Global exception handler"""
    logger.exception(f"Unhandled exception: {exc}")
    
    return JSONResponse(
        status="error",
        detail=str(exc),
        timestamp=datetime.now().isoformat()
    )

# Startup event
@app.on_event("startup")
async def startup_event():
    """Run on application startup"""
    logger.info("BioDockify API starting up")
    logger.info("Agent Zero self-repair system initialized")
    logger.info("Docker manager initialized")
    logger.info("Vina engine initialized")
    logger.info("ODDT analyzer initialized")
    logger.info("RDKit calculator initialized")

# Shutdown event
@app.on_event("shutdown")
async def shutdown_event():
    """Run on application shutdown"""
    logger.info("BioDockify API shutting down")
    database.close()
    logger.info("Database connection closed")
