"""
BioDockify Docking Studio - Database Management
Handles SQLite database operations for job storage and results persistence
"""

import sqlite3
import json
from pathlib import Path
from typing import Dict, List, Optional, Any
import logging
from datetime import datetime
import uuid
import os

logger = logging.getLogger(__name__)

class Database:
    """Database manager for BioDockify Docking Studio"""
    
    def __init__(self, db_path: Optional[Path] = None):
        """Initialize database connection"""
        if db_path is None:
            self.db_path = self._get_default_db_path()
        else:
            self.db_path = db_path
        
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self.connection = None
        self.cursor = None
        self._connect()
        self._initialize_tables()
    
    def _get_default_db_path(self) -> Path:
        """Get default database path"""
        if os.name == "nt":  # Windows
            db_dir = Path.home() / "AppData" / "Local" / "BioDockify" / "data"
        else:  # macOS/Linux
            db_dir = Path.home() / ".config" / "BioDockify" / "data"
        
        db_dir.mkdir(parents=True, exist_ok=True)
        return db_dir / "BioDockify.db"
    
    def _connect(self) -> None:
        """Connect to SQLite database"""
        try:
            self.connection = sqlite3.connect(self.db_path)
            # Enable foreign keys
            self.connection.execute("PRAGMA foreign_keys = ON")
            self.cursor = self.connection.cursor()
            logger.info(f"Connected to database: {self.db_path}")
        except sqlite3.Error as e:
            logger.error(f"Failed to connect to database: {e}")
            raise
    
    def _initialize_tables(self) -> None:
        """Initialize database tables"""
        try:
            with self.connection:
                self._create_jobs_table()
                self._create_checkpoints_table()
                self._create_results_table()
                self._create_logs_table()
            logger.info("Database tables initialized")
        except sqlite3.Error as e:
            logger.error(f"Failed to initialize tables: {e}")
            raise
    
    def _create_jobs_table(self) -> None:
        """Create jobs table"""
        self.connection.execute("""
            CREATE TABLE IF NOT EXISTS jobs (
                id TEXT PRIMARY KEY,
                status TEXT NOT NULL,
                receptor_file TEXT NOT NULL,
                ligand_file TEXT NOT NULL,
                parameters TEXT NOT NULL,
                created_at TIMESTAMP NOT NULL,
                updated_at TIMESTAMP NOT NULL,
                completed_at TIMESTAMP
            )
        """)
    
    def _create_checkpoints_table(self) -> None:
        """Create checkpoints table"""
        self.connection.execute("""
            CREATE TABLE IF NOT EXISTS checkpoints (
                id TEXT PRIMARY KEY,
                job_id TEXT NOT NULL,
                stage TEXT NOT NULL,
                data TEXT NOT NULL,
                created_at TIMESTAMP NOT NULL,
                FOREIGN KEY(job_id) REFERENCES jobs(id)
            )
        """)
    
    def _create_results_table(self) -> None:
        """Create results table"""
        self.connection.execute("""
            CREATE TABLE IF NOT EXISTS results (
                id TEXT PRIMARY KEY,
                job_id TEXT NOT NULL,
                binding_energy REAL NOT NULL,
                interactions TEXT NOT NULL,
                descriptors TEXT NOT NULL,
                confidence_score INTEGER NOT NULL,
                created_at TIMESTAMP NOT NULL,
                FOREIGN KEY(job_id) REFERENCES jobs(id)
            )
        """)
    
    def _create_logs_table(self) -> None:
        """Create logs table"""
        self.connection.execute("""
            CREATE TABLE IF NOT EXISTS logs (
                id TEXT PRIMARY KEY,
                job_id TEXT NOT NULL,
                level TEXT NOT NULL,
                message TEXT NOT NULL,
                details TEXT,
                created_at TIMESTAMP NOT NULL,
                FOREIGN KEY(job_id) REFERENCES jobs(id)
            )
        """)
    
    def create_job(self, receptor_file: str, ligand_file: str, 
                  parameters: Dict[str, Any]) -> str:
        """Create new docking job"""
        job_id = str(uuid.uuid4())
        created_at = datetime.now().isoformat()
        
        try:
            with self.connection:
                self.connection.execute("""
                    INSERT INTO jobs (id, status, receptor_file, ligand_file, parameters, created_at, updated_at)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                """, (job_id, "PENDING", receptor_file, ligand_file, 
                       json.dumps(parameters), created_at, created_at))
            logger.info(f"Job created: {job_id}")
            return job_id
        except sqlite3.Error as e:
            logger.error(f"Failed to create job: {e}")
            raise
    
    def update_job_status(self, job_id: str, status: str, 
                       completed_at: Optional[str] = None) -> bool:
        """Update job status"""
        updated_at = datetime.now().isoformat()
        
        try:
            with self.connection:
                if completed_at:
                    self.connection.execute("""
                        UPDATE jobs SET status = ?, updated_at = ?, completed_at = ?
                        WHERE id = ?
                    """, (status, updated_at, completed_at, job_id))
                else:
                    self.connection.execute("""
                        UPDATE jobs SET status = ?, updated_at = ?
                        WHERE id = ?
                    """, (status, updated_at, job_id))
            logger.info(f"Job {job_id} status updated to: {status}")
            return True
        except sqlite3.Error as e:
            logger.error(f"Failed to update job status: {e}")
            return False
    
    def get_job(self, job_id: str) -> Optional[Dict[str, Any]]:
        """Get job by ID"""
        try:
            cursor = self.connection.cursor()
            cursor.execute("SELECT * FROM jobs WHERE id = ?", (job_id,))
            row = cursor.fetchone()
            if row:
                return {
                    "id": row[0],
                    "status": row[1],
                    "receptor_file": row[2],
                    "ligand_file": row[3],
                    "parameters": json.loads(row[4]),
                    "created_at": row[5],
                    "updated_at": row[6],
                    "completed_at": row[7]
                }
            return None
        except sqlite3.Error as e:
            logger.error(f"Failed to get job {job_id}: {e}")
            return None
    
    def save_checkpoint(self, job_id: str, stage: str, data: str) -> bool:
        """Save checkpoint"""
        checkpoint_id = f"{job_id}_{stage}"
        created_at = datetime.now().isoformat()
        
        try:
            with self.connection:
                self.connection.execute("""
                    INSERT OR REPLACE INTO checkpoints (id, job_id, stage, data, created_at)
                    VALUES (?, ?, ?, ?, ?)
                """, (checkpoint_id, job_id, stage, json.dumps(data), created_at))
            logger.info(f"Checkpoint saved: {checkpoint_id}")
            return True
        except sqlite3.Error as e:
            logger.error(f"Failed to save checkpoint: {e}")
            return False
    
    def get_checkpoint(self, job_id: str, stage: str) -> Optional[Dict[str, Any]]:
        """Get checkpoint by job ID and stage"""
        checkpoint_id = f"{job_id}_{stage}"
        
        try:
            cursor = self.connection.cursor()
            cursor.execute("SELECT * FROM checkpoints WHERE id = ?", (checkpoint_id,))
            row = cursor.fetchone()
            if row:
                return {
                    "id": row[0],
                    "job_id": row[1],
                    "stage": row[2],
                    "data": json.loads(row[3]),
                    "created_at": row[4]
                }
            return None
        except sqlite3.Error as e:
            logger.error(f"Failed to get checkpoint {checkpoint_id}: {e}")
            return None
    
    def save_result(self, job_id: str, binding_energy: float, 
                 interactions: Dict[str, Any], descriptors: Dict[str, Any],
                 confidence_score: int) -> bool:
        """Save docking result"""
        result_id = str(uuid.uuid4())
        created_at = datetime.now().isoformat()
        
        try:
            with self.connection:
                self.connection.execute("""
                    INSERT INTO results (id, job_id, binding_energy, interactions, 
                                   descriptors, confidence_score, created_at)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                """, (result_id, job_id, binding_energy, json.dumps(interactions),
                       json.dumps(descriptors), confidence_score, created_at))
            logger.info(f"Result saved: {result_id}")
            return True
        except sqlite3.Error as e:
            logger.error(f"Failed to save result: {e}")
            return False
    
    def save_log(self, job_id: str, level: str, message: str, 
               details: Optional[str] = None) -> bool:
        """Save log entry"""
        log_id = str(uuid.uuid4())
        created_at = datetime.now().isoformat()
        
        try:
            with self.connection:
                self.connection.execute("""
                    INSERT INTO logs (id, job_id, level, message, details, created_at)
                    VALUES (?, ?, ?, ?, ?, ?)
                """, (log_id, job_id, level, message, details or json.dumps({}), created_at))
            logger.debug(f"Log saved: {log_id}")
            return True
        except sqlite3.Error as e:
            logger.error(f"Failed to save log: {e}")
            return False
    
    def get_recent_jobs(self, limit: int = 10) -> List[Dict[str, Any]]:
        """Get recent jobs"""
        try:
            cursor = self.connection.cursor()
            cursor.execute("""
                SELECT * FROM jobs 
                ORDER BY created_at DESC 
                LIMIT ?
            """, (limit,))
            rows = cursor.fetchall()
            
            jobs = []
            for row in rows:
                jobs.append({
                    "id": row[0],
                    "status": row[1],
                    "receptor_file": row[2],
                    "ligand_file": row[3],
                    "parameters": json.loads(row[4]),
                    "created_at": row[5],
                    "updated_at": row[6],
                    "completed_at": row[7]
                })
            
            return jobs
        except sqlite3.Error as e:
            logger.error(f"Failed to get recent jobs: {e}")
            return []
    
    def close(self) -> None:
        """Close database connection"""
        if self.connection:
            try:
                self.connection.close()
                self.connection = None
                logger.info("Database connection closed")
            except sqlite3.Error as e:
                logger.error(f"Error closing database: {e}")

    def __del__(self):
        """Clean up database connection on deletion"""
        if hasattr(self, 'connection') and self.connection:
            try:
                self.connection.close()
            except Exception:
                pass
            logger.info("Database connection closed via __del__")
