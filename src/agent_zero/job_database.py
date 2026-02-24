"""
Job Database Schema - SQLite
Job tracking and results storage
"""

import sqlite3
import json
from datetime import datetime
from typing import Optional, List, Dict, Any
from contextlib import contextmanager
import logging

logger = logging.getLogger(__name__)


class JobDatabase:
    """
    SQLite database for job management.
    """
    
    SCHEMA = """
    -- Jobs table
    CREATE TABLE IF NOT EXISTS jobs (
        job_id TEXT PRIMARY KEY,
        state TEXT NOT NULL DEFAULT 'QUEUED',
        created_at TEXT NOT NULL,
        updated_at TEXT NOT NULL,
        started_at TEXT,
        completed_at TEXT,
        
        -- Configuration
        receptor_file TEXT,
        ligand_files TEXT,  -- JSON array
        engine TEXT DEFAULT 'vina_cpu',
        
        -- Grid parameters
        center_x REAL DEFAULT 0,
        center_y REAL DEFAULT 0,
        center_z REAL DEFAULT 0,
        size_x REAL DEFAULT 20,
        size_y REAL DEFAULT 20,
        size_z REAL DEFAULT 20,
        
        -- Execution parameters
        exhaustiveness INTEGER DEFAULT 8,
        num_modes INTEGER DEFAULT 9,
        batch_size INTEGER DEFAULT 5,
        
        -- Progress
        progress INTEGER DEFAULT 0,
        message TEXT,
        error TEXT,
        
        -- Metrics
        retry_count INTEGER DEFAULT 0,
        runtime_seconds REAL DEFAULT 0,
        gpu_used INTEGER DEFAULT 0,
        
        -- Results (JSON)
        results_json TEXT
    );
    
    -- Job events for tracking
    CREATE TABLE IF NOT EXISTS job_events (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        job_id TEXT NOT NULL,
        event_type TEXT NOT NULL,
        event_data TEXT,  -- JSON
        timestamp TEXT NOT NULL,
        FOREIGN KEY (job_id) REFERENCES jobs(job_id)
    );
    
    -- Job logs
    CREATE TABLE IF NOT EXISTS job_logs (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        job_id TEXT NOT NULL,
        level TEXT NOT NULL,
        message TEXT NOT NULL,
        timestamp TEXT NOT NULL,
        FOREIGN KEY (job_id) REFERENCES jobs(job_id)
    );
    
    -- Docking poses
    CREATE TABLE IF NOT EXISTS poses (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        job_id TEXT NOT NULL,
        pose_index INTEGER NOT NULL,
        rank INTEGER NOT NULL,
        
        -- Scores
        vina_score REAL,
        gnina_score REAL,
        consensus_score REAL,
        
        -- File path
        pdb_file TEXT,
        
        -- Metadata
        created_at TEXT NOT NULL,
        
        FOREIGN KEY (job_id) REFERENCES jobs(job_id),
        UNIQUE(job_id, pose_index)
    );
    
    -- Failures for diagnostics
    CREATE TABLE IF NOT EXISTS failures (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        job_id TEXT NOT NULL,
        error_type TEXT NOT NULL,
        error_message TEXT,
        severity TEXT,
        retry_number INTEGER,
        fixed INTEGER DEFAULT 0,
        timestamp TEXT NOT NULL,
        FOREIGN KEY (job_id) REFERENCES jobs(job_id)
    );
    
    -- Indexes
    CREATE INDEX IF NOT EXISTS idx_jobs_state ON jobs(state);
    CREATE INDEX IF NOT EXISTS idx_jobs_created ON jobs(created_at);
    CREATE INDEX IF NOT EXISTS idx_poses_job ON poses(job_id);
    CREATE INDEX IF NOT EXISTS idx_logs_job ON job_logs(job_id);
    CREATE INDEX IF NOT EXISTS idx_failures_job ON failures(job_id);
    """
    
    def __init__(self, db_path: str = "data/docking.db"):
        """Initialize database"""
        self.db_path = db_path
        
        # Ensure directory exists
        import os
        os.makedirs(os.path.dirname(db_path), exist_ok=True)
        
        self._init_db()
        logger.info(f"JobDatabase initialized: {db_path}")
    
    def _init_db(self):
        """Initialize database schema"""
        with self._get_connection() as conn:
            conn.executescript(self.SCHEMA)
            conn.commit()
    
    @contextmanager
    def _get_connection(self):
        """Get database connection"""
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        try:
            yield conn
        finally:
            conn.close()
    
    # ==================== Job Operations ====================
    
    def create_job(
        self,
        job_id: str,
        receptor_file: str,
        ligand_files: List[str],
        engine: str = "vina_cpu",
        **config
    ) -> bool:
        """Create a new job"""
        now = datetime.now().isoformat()
        
        with self._get_connection() as conn:
            conn.execute("""
                INSERT INTO jobs (
                    job_id, state, created_at, updated_at,
                    receptor_file, ligand_files, engine,
                    center_x, center_y, center_z,
                    size_x, size_y, size_z,
                    exhaustiveness, num_modes, batch_size
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                job_id, 'QUEUED', now, now,
                receptor_file, json.dumps(ligand_files), engine,
                config.get('center_x', 0),
                config.get('center_y', 0),
                config.get('center_z', 0),
                config.get('size_x', 20),
                config.get('size_y', 20),
                config.get('size_z', 20),
                config.get('exhaustiveness', 8),
                config.get('num_modes', 9),
                config.get('batch_size', 5),
            ))
            conn.commit()
        
        logger.info(f"Job created: {job_id}")
        return True
    
    def update_job_state(
        self,
        job_id: str,
        state: str,
        progress: Optional[int] = None,
        message: Optional[str] = None,
        error: Optional[str] = None
    ) -> bool:
        """Update job state"""
        now = datetime.now().isoformat()
        
        updates = ["state = ?", "updated_at = ?"]
        values = [state, now]
        
        if progress is not None:
            updates.append("progress = ?")
            values.append(progress)
        
        if message is not None:
            updates.append("message = ?")
            values.append(message)
        
        if error is not None:
            updates.append("error = ?")
            values.append(error)
        
        if state == "RUNNING":
            updates.append("started_at = ?")
            values.append(now)
        
        if state in ["COMPLETED", "FAILED", "CANCELLED"]:
            updates.append("completed_at = ?")
            values.append(now)
        
        values.append(job_id)
        
        with self._get_connection() as conn:
            conn.execute(
                f"UPDATE jobs SET {', '.join(updates)} WHERE job_id = ?",
                values
            )
            conn.commit()
        
        return True
    
    def get_job(self, job_id: str) -> Optional[Dict]:
        """Get job by ID"""
        with self._get_connection() as conn:
            row = conn.execute(
                "SELECT * FROM jobs WHERE job_id = ?", (job_id,)
            ).fetchone()
        
        if row:
            return dict(row)
        return None
    
    def get_all_jobs(self, state: Optional[str] = None) -> List[Dict]:
        """Get all jobs, optionally filtered by state"""
        with self._get_connection() as conn:
            if state:
                rows = conn.execute(
                    "SELECT * FROM jobs WHERE state = ? ORDER BY created_at DESC",
                    (state,)
                ).fetchall()
            else:
                rows = conn.execute(
                    "SELECT * FROM jobs ORDER BY created_at DESC"
                ).fetchall()
        
        return [dict(row) for row in rows]
    
    def delete_job(self, job_id: str) -> bool:
        """Delete a job"""
        with self._get_connection() as conn:
            conn.execute("DELETE FROM job_logs WHERE job_id = ?", (job_id,))
            conn.execute("DELETE FROM failures WHERE job_id = ?", (job_id,))
            conn.execute("DELETE FROM poses WHERE job_id = ?", (job_id,))
            conn.execute("DELETE FROM job_events WHERE job_id = ?", (job_id,))
            conn.execute("DELETE FROM jobs WHERE job_id = ?", (job_id,))
            conn.commit()
        
        return True
    
    # ==================== Pose Operations ====================
    
    def add_pose(
        self,
        job_id: str,
        pose_index: int,
        rank: int,
        vina_score: Optional[float] = None,
        gnina_score: Optional[float] = None,
        consensus_score: Optional[float] = None,
        pdb_file: Optional[str] = None
    ) -> bool:
        """Add a docking pose"""
        now = datetime.now().isoformat()
        
        with self._get_connection() as conn:
            conn.execute("""
                INSERT OR REPLACE INTO poses (
                    job_id, pose_index, rank,
                    vina_score, gnina_score, consensus_score,
                    pdb_file, created_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                job_id, pose_index, rank,
                vina_score, gnina_score, consensus_score,
                pdb_file, now
            ))
            conn.commit()
        
        return True
    
    def get_poses(self, job_id: str) -> List[Dict]:
        """Get all poses for a job"""
        with self._get_connection() as conn:
            rows = conn.execute(
                "SELECT * FROM poses WHERE job_id = ? ORDER BY rank",
                (job_id,)
            ).fetchall()
        
        return [dict(row) for row in rows]
    
    # ==================== Event Operations ====================
    
    def log_event(
        self,
        job_id: str,
        event_type: str,
        event_data: Optional[Dict] = None
    ):
        """Log a job event"""
        now = datetime.now().isoformat()
        
        with self._get_connection() as conn:
            conn.execute("""
                INSERT INTO job_events (job_id, event_type, event_data, timestamp)
                VALUES (?, ?, ?, ?)
            """, (
                job_id, event_type,
                json.dumps(event_data) if event_data else None,
                now
            ))
            conn.commit()
    
    # ==================== Failure Operations ====================
    
    def log_failure(
        self,
        job_id: str,
        error_type: str,
        error_message: str,
        severity: str = "error",
        retry_number: int = 0
    ) -> bool:
        """Log a failure"""
        now = datetime.now().isoformat()
        
        with self._get_connection() as conn:
            conn.execute("""
                INSERT INTO failures (
                    job_id, error_type, error_message,
                    severity, retry_number, timestamp
                ) VALUES (?, ?, ?, ?, ?, ?)
            """, (
                job_id, error_type, error_message,
                severity, retry_number, now
            ))
            conn.commit()
        
        return True
    
    def get_failures(self, job_id: str) -> List[Dict]:
        """Get all failures for a job"""
        with self._get_connection() as conn:
            rows = conn.execute(
                "SELECT * FROM failures WHERE job_id = ? ORDER BY timestamp",
                (job_id,)
            ).fetchall()
        
        return [dict(row) for row in rows]
    
    # ==================== Log Operations ====================
    
    def log_message(
        self,
        job_id: str,
        level: str,
        message: str
    ):
        """Log a message"""
        now = datetime.now().isoformat()
        
        with self._get_connection() as conn:
            conn.execute("""
                INSERT INTO job_logs (job_id, level, message, timestamp)
                VALUES (?, ?, ?, ?)
            """, (job_id, level, message, now))
            conn.commit()
    
    def get_logs(self, job_id: str) -> List[Dict]:
        """Get all logs for a job"""
        with self._get_connection() as conn:
            rows = conn.execute(
                "SELECT * FROM job_logs WHERE job_id = ? ORDER BY timestamp",
                (job_id,)
            ).fetchall()
        
        return [dict(row) for row in rows]
    
    # ==================== Statistics ====================
    
    def get_stats(self) -> Dict:
        """Get database statistics"""
        with self._get_connection() as conn:
            total = conn.execute("SELECT COUNT(*) FROM jobs").fetchone()[0]
            queued = conn.execute("SELECT COUNT(*) FROM jobs WHERE state = 'QUEUED'").fetchone()[0]
            running = conn.execute("SELECT COUNT(*) FROM jobs WHERE state = 'RUNNING'").fetchone()[0]
            completed = conn.execute("SELECT COUNT(*) FROM jobs WHERE state = 'COMPLETED'").fetchone()[0]
            failed = conn.execute("SELECT COUNT(*) FROM jobs WHERE state = 'FAILED'").fetchone()[0]
            
            avg_runtime = conn.execute(
                "SELECT AVG(runtime_seconds) FROM jobs WHERE state = 'COMPLETED'"
            ).fetchone()[0] or 0
            
            return {
                "total_jobs": total,
                "queued": queued,
                "running": running,
                "completed": completed,
                "failed": failed,
                "avg_runtime_seconds": avg_runtime,
            }


def get_database(db_path: str = "data/docking.db") -> JobDatabase:
    """Get or create database instance"""
    return JobDatabase(db_path)
