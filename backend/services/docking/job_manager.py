import json
import os
from pathlib import Path
from typing import Optional

class DockingJobManager:
    def __init__(self, jobs_dir: Path):
        self.jobs_dir = jobs_dir
        self.jobs_dir.mkdir(parents=True, exist_ok=True)
    
    def create_job(self, job_id: str, params: dict) -> dict:
        job_info = {
            'job_id': job_id,
            'status': 'created',
            'params': params,
            'created_at': self._timestamp()
        }
        self._save_job(job_id, job_info)
        return job_info
    
    def update_status(self, job_id: str, status: str, results: Optional[dict] = None) -> dict:
        job_info = self._load_job(job_id)
        if job_info:
            job_info['status'] = status
            job_info['updated_at'] = self._timestamp()
            if results:
                job_info['results'] = results
            self._save_job(job_id, job_info)
        return job_info or {}
    
    def get_job(self, job_id: str) -> Optional[dict]:
        return self._load_job(job_id)
    
    def list_jobs(self, status: Optional[str] = None) -> list:
        jobs = []
        for job_file in self.jobs_dir.glob("*.json"):
            job = self._load_job(job_file.stem)
            if job:
                if status is None or job.get('status') == status:
                    jobs.append(job)
        return sorted(jobs, key=lambda x: x.get('created_at', ''), reverse=True)
    
    def _save_job(self, job_id: str, job_info: dict):
        job_file = self.jobs_dir / f"{job_id}.json"
        job_file.write_text(json.dumps(job_info, indent=2))
    
    def _load_job(self, job_id: str) -> Optional[dict]:
        job_file = self.jobs_dir / f"{job_id}.json"
        if job_file.exists():
            return json.loads(job_file.read_text())
        return None
    
    def _timestamp(self) -> str:
        from datetime import datetime
        return datetime.utcnow().isoformat()
