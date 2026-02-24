"""
Job Manager - Docking Job State Machine
Manages job states from QUEUED to COMPLETED/FAILED
"""

import os
import json
import time
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field, asdict
from datetime import datetime
from enum import Enum
import logging

logger = logging.getLogger(__name__)


class JobState(Enum):
    """Job state machine states"""
    QUEUED = "QUEUED"
    PREPARING = "PREPARING"
    RUNNING = "RUNNING"
    RESCORING = "RESCORING"
    REFINING = "REFINING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    CANCELLED = "CANCELLED"


class EngineType(Enum):
    """Docking engine types"""
    VINA_CPU = "vina_cpu"
    VINA_GPU = "vina_gpu"
    GNINA_CPU = "gnina_cpu"
    GNINA_GPU = "gnina_gpu"
    CONSENSUS = "consensus"


@dataclass
class JobConfig:
    """Job configuration"""
    receptor_file: str
    ligand_files: List[str]
    engine: EngineType = EngineType.VINA_CPU
    grid_center: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    grid_size: List[float] = field(default_factory=lambda: [20.0, 20.0, 20.0])
    exhaustiveness: int = 8
    num_modes: int = 9
    batch_size: int = 5


@dataclass
class JobMetrics:
    """Job performance metrics"""
    start_time: Optional[str] = None
    end_time: Optional[str] = None
    duration_seconds: float = 0.0
    retry_count: int = 0
    gpu_used: bool = False
    gpu_memory_mb: int = 0
    cpu_usage_percent: float = 0.0
    peak_memory_mb: int = 0


@dataclass
class Job:
    """Complete docking job"""
    job_id: str
    state: JobState = JobState.QUEUED
    config: JobConfig = field(default_factory=JobConfig)
    progress: int = 0
    message: str = ""
    error: Optional[str] = None
    metrics: JobMetrics = field(default_factory=JobMetrics)
    results: Dict[str, Any] = field(default_factory=dict)
    created_at: str = field(default_factory=lambda: datetime.now().isoformat())
    updated_at: str = field(default_factory=lambda: datetime.now().isoformat())


class JobManager:
    """
    Manages docking job lifecycle and state.
    """
    
    def __init__(self, data_dir: str = "data"):
        """
        Initialize job manager.
        
        Args:
            data_dir: Directory for job data storage
        """
        self.data_dir = data_dir
        self.jobs_dir = os.path.join(data_dir, "jobs")
        os.makedirs(self.jobs_dir, exist_ok=True)
        
        self.active_jobs: Dict[str, Job] = {}
        
        logger.info(f"JobManager initialized: {self.jobs_dir}")
    
    def create_job(
        self,
        job_id: str,
        receptor: str,
        ligands: List[str],
        engine: str = "vina_cpu",
        **kwargs
    ) -> Job:
        """
        Create a new docking job.
        
        Args:
            job_id: Unique job identifier
            receptor: Receptor file path
            ligands: List of ligand file paths
            engine: Engine type
            **kwargs: Additional configuration
            
        Returns:
            Job object
        """
        job = Job(
            job_id=job_id,
            config=JobConfig(
                receptor_file=receptor,
                ligand_files=ligands,
                engine=EngineType(engine),
                grid_center=kwargs.get('grid_center', [0.0, 0.0, 0.0]),
                grid_size=kwargs.get('grid_size', [20.0, 20.0, 20.0]),
                exhaustiveness=kwargs.get('exhaustiveness', 8),
                num_modes=kwargs.get('num_modes', 9),
                batch_size=kwargs.get('batch_size', 5),
            )
        )
        
        self.active_jobs[job_id] = job
        self._save_job(job)
        
        logger.info(f"Job created: {job_id}")
        
        return job
    
    def get_job(self, job_id: str) -> Optional[Job]:
        """Get job by ID"""
        if job_id in self.active_jobs:
            return self.active_jobs[job_id]
        
        # Try loading from disk
        job = self._load_job(job_id)
        if job:
            self.active_jobs[job_id] = job
        
        return job
    
    def update_state(
        self,
        job_id: str,
        state: JobState,
        progress: Optional[int] = None,
        message: Optional[str] = None,
        error: Optional[str] = None
    ) -> bool:
        """
        Update job state.
        
        Args:
            job_id: Job identifier
            state: New state
            progress: Optional progress percentage
            message: Optional status message
            error: Optional error message
            
        Returns:
            True if successful
        """
        job = self.get_job(job_id)
        if not job:
            logger.warning(f"Job not found: {job_id}")
            return False
        
        old_state = job.state
        job.state = state
        job.updated_at = datetime.now().isoformat()
        
        if progress is not None:
            job.progress = progress
        
        if message:
            job.message = message
        
        if error:
            job.error = error
            job.state = JobState.FAILED
        
        # Update timing
        if state == JobState.RUNNING and old_state != JobState.RUNNING:
            job.metrics.start_time = datetime.now().isoformat()
        
        if state in [JobState.COMPLETED, JobState.FAILED]:
            job.metrics.end_time = datetime.now().isoformat()
            if job.metrics.start_time:
                start = datetime.fromisoformat(job.metrics.start_time)
                end = datetime.fromisoformat(job.metrics.end_time)
                job.metrics.duration_seconds = (end - start).total_seconds()
        
        self._save_job(job)
        
        logger.info(f"Job {job_id}: {old_state.value} → {state.value}")
        
        return True
    
    def update_progress(self, job_id: str, progress: int, message: str = ""):
        """Update job progress"""
        job = self.get_job(job_id)
        if job:
            job.progress = min(100, max(0, progress))
            job.message = message
            job.updated_at = datetime.now().isoformat()
            self._save_job(job)
    
    def set_results(self, job_id: str, results: Dict[str, Any]):
        """Set job results"""
        job = self.get_job(job_id)
        if job:
            job.results = results
            job.updated_at = datetime.now().isoformat()
            self._save_job(job)
    
    def increment_retry(self, job_id: str):
        """Increment retry count"""
        job = self.get_job(job_id)
        if job:
            job.metrics.retry_count += 1
            job.updated_at = datetime.now().isoformat()
            self._save_job(job)
            logger.info(f"Job {job_id} retry count: {job.metrics.retry_count}")
    
    def get_all_jobs(self) -> List[Job]:
        """Get all jobs"""
        jobs = []
        
        for filename in os.listdir(self.jobs_dir):
            if filename.endswith('.json'):
                job_id = filename.replace('.json', '')
                job = self.get_job(job_id)
                if job:
                    jobs.append(job)
        
        return sorted(jobs, key=lambda j: j.created_at, reverse=True)
    
    def get_active_jobs(self) -> List[Job]:
        """Get jobs that are not completed/failed/cancelled"""
        return [
            j for j in self.active_jobs.values()
            if j.state not in [JobState.COMPLETED, JobState.FAILED, JobState.CANCELLED]
        ]
    
    def get_job_status(self, job_id: str) -> Dict[str, Any]:
        """Get job status for GUI"""
        job = self.get_job(job_id)
        if not job:
            return {"status": "not_found"}
        
        return {
            "job_id": job.job_id,
            "state": job.state.value,
            "progress": job.progress,
            "message": job.message,
            "error": job.error,
            "created_at": job.created_at,
            "updated_at": job.updated_at,
            "duration_seconds": job.metrics.duration_seconds,
            "retry_count": job.metrics.retry_count,
            "engine": job.config.engine.value,
        }
    
    def cancel_job(self, job_id: str) -> bool:
        """Cancel a job"""
        return self.update_state(job_id, JobState.CANCELLED, message="Cancelled by user")
    
    def cleanup_completed(self, keep_days: int = 7):
        """Clean up old completed jobs"""
        cutoff = datetime.now().timestamp() - (keep_days * 86400)
        
        for job in self.get_all_jobs():
            if job.state in [JobState.COMPLETED, JobState.FAILED, JobState.CANCELLED]:
                completed_time = datetime.fromisoformat(job.updated_at).timestamp()
                if completed_time < cutoff:
                    # Archive or delete
                    logger.info(f"Cleaning up old job: {job.job_id}")
    
    def _job_file_path(self, job_id: str) -> str:
        """Get job file path"""
        return os.path.join(self.jobs_dir, f"{job_id}.json")
    
    def _save_job(self, job: Job):
        """Save job to disk"""
        filepath = self._job_file_path(job.job_id)
        
        # Convert to dict
        data = {
            "job_id": job.job_id,
            "state": job.state.value,
            "config": asdict(job.config),
            "progress": job.progress,
            "message": job.message,
            "error": job.error,
            "metrics": asdict(job.metrics),
            "results": job.results,
            "created_at": job.created_at,
            "updated_at": job.updated_at,
        }
        
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2)
    
    def _load_job(self, job_id: str) -> Optional[Job]:
        """Load job from disk"""
        filepath = self._job_file_path(job_id)
        
        if not os.path.exists(filepath):
            return None
        
        try:
            with open(filepath, 'r') as f:
                data = json.load(f)
            
            # Reconstruct Job object
            job = Job(
                job_id=data["job_id"],
                state=JobState(data["state"]),
                progress=data.get("progress", 0),
                message=data.get("message", ""),
                error=data.get("error"),
                results=data.get("results", {}),
                created_at=data.get("created_at", datetime.now().isoformat()),
                updated_at=data.get("updated_at", datetime.now().isoformat()),
            )
            
            # Reconstruct config
            if "config" in data:
                job.config = JobConfig(**data["config"])
                job.config.engine = EngineType(data["config"].get("engine", "vina_cpu"))
            
            # Reconstruct metrics
            if "metrics" in data:
                job.metrics = JobMetrics(**data["metrics"])
            
            return job
            
        except Exception as e:
            logger.error(f"Failed to load job {job_id}: {e}")
            return None


def create_job_id() -> str:
    """Generate a unique job ID"""
    return f"job_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
