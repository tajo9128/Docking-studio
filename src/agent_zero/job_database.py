"""
BioDockify Job Database
Production-grade SQLite job tracking system.
"""

import sqlite3
import uuid
import os
from datetime import datetime
from typing import Optional, List, Dict, Any
from dataclasses import dataclass
from pathlib import Path


DB_PATH = "data/biodockify_jobs.db"


@dataclass
class User:
    """User record."""
    id: int
    email: str
    password_hash: Optional[str] = None
    created_at: Optional[str] = None


@dataclass
class Job:
    """Job record."""
    id: int
    user_id: int
    job_uuid: str
    job_name: Optional[str]
    status: str
    compute_mode: Optional[str]
    gpu_used: Optional[bool]
    gpu_name: Optional[str]
    vina_version: Optional[str]
    gnina_version: Optional[str]
    error_message: Optional[str]
    created_at: str
    completed_at: Optional[str]


@dataclass
class JobFile:
    """Job file record."""
    id: int
    job_id: int
    file_type: str
    file_path: str
    file_size: Optional[int]
    created_at: str


@dataclass
class DockingResult:
    """Docking result record."""
    id: int
    job_id: int
    ligand_name: str
    vina_score: Optional[float]
    gnina_score: Optional[float]
    rf_score: Optional[float]
    consensus_score: float
    rank: int


@dataclass
class GridMetadata:
    """Grid metadata record."""
    job_id: int
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float
    exhaustiveness: int
    seed: int


class JobDB:
    """
    SQLite job tracking database manager.
    """
    
    def __init__(self, db_path: str = None):
        self.db_path = db_path or DB_PATH
        self._ensure_dir()
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        self.create_tables()
    
    def _ensure_dir(self):
        """Ensure database directory exists."""
        Path(self.db_path).parent.mkdir(parents=True, exist_ok=True)
    
    def create_tables(self):
        """Create all database tables (backwards compatible)."""
        cursor = self.conn.cursor()
        
        # Main jobs table - keep existing schema + add new columns safely
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS jobs (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                user_id INTEGER NOT NULL,
                job_uuid TEXT UNIQUE NOT NULL,
                job_name TEXT,
                status TEXT CHECK(status IN ('pending', 'running', 'completed', 'failed', 'cancelled')),
                compute_mode TEXT,
                gpu_used BOOLEAN,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                completed_at DATETIME,
                FOREIGN KEY(user_id) REFERENCES users(id)
            )
        """)
        
        # Add new columns if not exist (migration)
        try:
            cursor.execute("ALTER TABLE jobs ADD COLUMN gpu_name TEXT")
        except sqlite3.OperationalError:
            pass  # Column exists
        
        try:
            cursor.execute("ALTER TABLE jobs ADD COLUMN vina_version TEXT")
        except sqlite3.OperationalError:
            pass
        
        try:
            cursor.execute("ALTER TABLE jobs ADD COLUMN gnina_version TEXT")
        except sqlite3.OperationalError:
            pass
        
        try:
            cursor.execute("ALTER TABLE jobs ADD COLUMN error_message TEXT")
        except sqlite3.OperationalError:
            pass
        
        try:
            cursor.execute("ALTER TABLE jobs ADD COLUMN progress INTEGER DEFAULT 0")
        except sqlite3.OperationalError:
            pass
        
        # Users table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS users (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                email TEXT UNIQUE NOT NULL,
                password_hash TEXT NOT NULL,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP
            )
        """)
        
        # Job files table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS job_files (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                job_id INTEGER NOT NULL,
                file_type TEXT,
                file_path TEXT NOT NULL,
                file_size INTEGER,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY(job_id) REFERENCES jobs(id)
            )
        """)
        
        # Docking results table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS docking_results (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                job_id INTEGER NOT NULL,
                ligand_name TEXT,
                vina_score REAL,
                gnina_score REAL,
                rf_score REAL,
                consensus_score REAL,
                rank INTEGER,
                FOREIGN KEY(job_id) REFERENCES jobs(id)
            )
        """)
        
        # Grid metadata table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS grid_metadata (
                job_id INTEGER PRIMARY KEY,
                center_x REAL,
                center_y REAL,
                center_z REAL,
                size_x REAL,
                size_y REAL,
                size_z REAL,
                exhaustiveness INTEGER,
                seed INTEGER,
                FOREIGN KEY(job_id) REFERENCES jobs(id)
            )
        """)
        
        # Create indexes
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_jobs_user_id ON jobs(user_id)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_jobs_uuid ON jobs(job_uuid)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_jobs_status ON jobs(status)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_docking_results_job ON docking_results(job_id)")
        
        self.conn.commit()
    
    # ==================== USER METHODS ====================
    
    def create_user(self, email: str, password_hash: str) -> int:
        """Create new user."""
        cursor = self.conn.cursor()
        cursor.execute(
            "INSERT INTO users (email, password_hash) VALUES (?, ?)",
            (email, password_hash)
        )
        self.conn.commit()
        return cursor.lastrowid
    
    def get_user_by_email(self, email: str) -> Optional[User]:
        """Get user by email."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT * FROM users WHERE email = ?", (email,))
        row = cursor.fetchone()
        return User(**dict(row)) if row else None
    
    def get_user(self, user_id: int) -> Optional[User]:
        """Get user by ID."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT * FROM users WHERE id = ?", (user_id,))
        row = cursor.fetchone()
        return User(**dict(row)) if row else None
    
    # ==================== JOB METHODS ====================
    
    def create_job(
        self,
        user_id: int,
        job_name: str,
        compute_mode: str = "auto",
        gpu_used: bool = False,
        gpu_name: str = None,
        vina_version: str = None,
        gnina_version: str = None
    ) -> str:
        """Create new job (backwards compatible)."""
        job_uuid = str(uuid.uuid4())
        cursor = self.conn.cursor()
        cursor.execute("""
            INSERT INTO jobs (
                user_id, job_uuid, job_name, status, compute_mode, gpu_used
            ) VALUES (?, ?, ?, ?, ?, ?)
        """, (user_id, job_uuid, job_name, "pending", compute_mode, gpu_used))
        
        # Add optional fields if they exist
        if gpu_name:
            try:
                cursor.execute("UPDATE jobs SET gpu_name = ? WHERE job_uuid = ?", (gpu_name, job_uuid))
            except:
                pass
        if vina_version:
            try:
                cursor.execute("UPDATE jobs SET vina_version = ? WHERE job_uuid = ?", (vina_version, job_uuid))
            except:
                pass
        if gnina_version:
            try:
                cursor.execute("UPDATE jobs SET gnina_version = ? WHERE job_uuid = ?", (gnina_version, job_uuid))
            except:
                pass
        
        self.conn.commit()
        return job_uuid
    
    def update_progress(self, job_uuid: str, progress: int) -> None:
        """Update job progress (0-100)."""
        cursor = self.conn.cursor()
        cursor.execute("UPDATE jobs SET progress = ? WHERE job_uuid = ?", (progress, job_uuid))
        self.conn.commit()
    
    def update_job_status(
        self,
        job_uuid: str,
        status: str,
        error_message: str = None
    ) -> None:
        """Update job status."""
        cursor = self.conn.cursor()
        completed_at = datetime.utcnow().isoformat() if status in ("completed", "failed", "cancelled") else None
        cursor.execute("""
            UPDATE jobs
            SET status = ?, completed_at = ?, error_message = ?
            WHERE job_uuid = ?
        """, (status, completed_at, error_message, job_uuid))
        self.conn.commit()
    
    def get_job(self, job_uuid: str) -> Optional[Job]:
        """Get job by UUID."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT * FROM jobs WHERE job_uuid = ?", (job_uuid,))
        row = cursor.fetchone()
        return Job(**dict(row)) if row else None
    
    def get_jobs(self, user_id: int, status: str = None, limit: int = 50) -> List[Job]:
        """Get user's jobs."""
        cursor = self.conn.cursor()
        if status:
            cursor.execute("""
                SELECT * FROM jobs
                WHERE user_id = ? AND status = ?
                ORDER BY created_at DESC
                LIMIT ?
            """, (user_id, status, limit))
        else:
            cursor.execute("""
                SELECT * FROM jobs
                WHERE user_id = ?
                ORDER BY created_at DESC
                LIMIT ?
            """, (user_id, limit))
        return [Job(**dict(row)) for row in cursor.fetchall()]
    
    def get_all_jobs(self, status: str = None, limit: int = 100) -> List[Job]:
        """Get all jobs (admin)."""
        cursor = self.conn.cursor()
        if status:
            cursor.execute("""
                SELECT * FROM jobs
                WHERE status = ?
                ORDER BY created_at DESC
                LIMIT ?
            """, (status, limit))
        else:
            cursor.execute("SELECT * FROM jobs ORDER BY created_at DESC LIMIT ?", (limit,))
        return [Job(**dict(row)) for row in cursor.fetchall()]
    
    def delete_job(self, job_uuid: str) -> bool:
        """Delete job and all related data."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT id FROM jobs WHERE job_uuid = ?", (job_uuid,))
        job = cursor.fetchone()
        if not job:
            return False
        
        job_id = job[0]
        cursor.execute("DELETE FROM docking_results WHERE job_id = ?", (job_id,))
        cursor.execute("DELETE FROM job_files WHERE job_id = ?", (job_id,))
        cursor.execute("DELETE FROM grid_metadata WHERE job_id = ?", (job_id,))
        cursor.execute("DELETE FROM jobs WHERE job_uuid = ?", (job_uuid,))
        self.conn.commit()
        return True
    
    # ==================== FILE METHODS ====================
    
    def add_job_file(
        self,
        job_uuid: str,
        file_type: str,
        file_path: str,
        file_size: int = None
    ) -> None:
        """Add file record."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT id FROM jobs WHERE job_uuid = ?", (job_uuid,))
        job = cursor.fetchone()
        if not job:
            raise ValueError(f"Job not found: {job_uuid}")
        
        cursor.execute("""
            INSERT INTO job_files (job_id, file_type, file_path, file_size)
            VALUES (?, ?, ?, ?)
        """, (job[0], file_type, file_path, file_size))
        self.conn.commit()
    
    def get_job_files(self, job_uuid: str) -> List[JobFile]:
        """Get all files for a job."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT id FROM jobs WHERE job_uuid = ?", (job_uuid,))
        job = cursor.fetchone()
        if not job:
            return []
        
        cursor.execute("""
            SELECT * FROM job_files WHERE job_id = ?
            ORDER BY created_at
        """, (job[0],))
        return [JobFile(**dict(row)) for row in cursor.fetchall()]
    
    # ==================== RESULT METHODS ====================
    
    def add_result(
        self,
        job_uuid: str,
        ligand_name: str,
        vina_score: float = None,
        gnina_score: float = None,
        rf_score: float = None,
        consensus_score: float = None,
        rank: int = None
    ) -> None:
        """Add docking result."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT id FROM jobs WHERE job_uuid = ?", (job_uuid,))
        job = cursor.fetchone()
        if not job:
            raise ValueError(f"Job not found: {job_uuid}")
        
        cursor.execute("""
            INSERT INTO docking_results (
                job_id, ligand_name, vina_score, gnina_score,
                rf_score, consensus_score, rank
            ) VALUES (?, ?, ?, ?, ?, ?, ?)
        """, (
            job[0], ligand_name, vina_score, gnina_score,
            rf_score, consensus_score, rank
        ))
        self.conn.commit()
    
    def add_results_batch(self, job_uuid: str, results: List[Dict]) -> None:
        """Add multiple results."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT id FROM jobs WHERE job_uuid = ?", (job_uuid,))
        job = cursor.fetchone()
        if not job:
            raise ValueError(f"Job not found: {job_uuid}")
        
        job_id = job[0]
        for r in results:
            cursor.execute("""
                INSERT INTO docking_results (
                    job_id, ligand_name, vina_score, gnina_score,
                    rf_score, consensus_score, rank
                ) VALUES (?, ?, ?, ?, ?, ?, ?)
            """, (
                job_id, r.get("ligand"), r.get("vina_affinity"),
                r.get("gnina_cnn_affinity"), r.get("rf_score"),
                r.get("consensus_score"), r.get("rank")
            ))
        self.conn.commit()
    
    def get_results(self, job_uuid: str) -> List[DockingResult]:
        """Get all results for a job."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT id FROM jobs WHERE job_uuid = ?", (job_uuid,))
        job = cursor.fetchone()
        if not job:
            return []
        
        cursor.execute("""
            SELECT * FROM docking_results
            WHERE job_id = ?
            ORDER BY rank ASC
        """, (job[0],))
        return [DockingResult(**dict(row)) for row in cursor.fetchall()]
    
    # ==================== GRID METHODS ====================
    
    def add_grid_metadata(
        self,
        job_uuid: str,
        center_x: float,
        center_y: float,
        center_z: float,
        size_x: float,
        size_y: float,
        size_z: float,
        exhaustiveness: int,
        seed: int
    ) -> None:
        """Add grid metadata."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT id FROM jobs WHERE job_uuid = ?", (job_uuid,))
        job = cursor.fetchone()
        if not job:
            raise ValueError(f"Job not found: {job_uuid}")
        
        cursor.execute("""
            INSERT OR REPLACE INTO grid_metadata VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (job[0], center_x, center_y, center_z, size_x, size_y, size_z, exhaustiveness, seed))
        self.conn.commit()
    
    def get_grid_metadata(self, job_uuid: str) -> Optional[GridMetadata]:
        """Get grid metadata for a job."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT id FROM jobs WHERE job_uuid = ?", (job_uuid,))
        job = cursor.fetchone()
        if not job:
            return None
        
        cursor.execute("SELECT * FROM grid_metadata WHERE job_id = ?", (job[0],))
        row = cursor.fetchone()
        return GridMetadata(**dict(row)) if row else None
    
    # ==================== STATS ====================
    
    def get_user_stats(self, user_id: int) -> Dict:
        """Get user statistics."""
        cursor = self.conn.cursor()
        
        cursor.execute("""
            SELECT status, COUNT(*) as count
            FROM jobs WHERE user_id = ?
            GROUP BY status
        """, (user_id,))
        status_counts = {row[0]: row[1] for row in cursor.fetchall()}
        
        cursor.execute("""
            SELECT COUNT(*) FROM jobs WHERE user_id = ?
        """, (user_id,))
        total = cursor.fetchone()[0]
        
        return {
            "total_jobs": total,
            "by_status": status_counts
        }
    
    def get_job_summary(self, job_uuid: str) -> Optional[Dict]:
        """Get full job summary."""
        job = self.get_job(job_uuid)
        if not job:
            return None
        
        results = self.get_results(job_uuid)
        files = self.get_job_files(job_uuid)
        grid = self.get_grid_metadata(job_uuid)
        
        return {
            "job": job,
            "results": results,
            "files": files,
            "grid": grid
        }
    
    def close(self):
        """Close database connection."""
        self.conn.close()


# ==================== HELPER FUNCTIONS ====================

def get_db(db_path: str = None) -> JobDB:
    """Get database instance."""
    return JobDB(db_path)


def init_demo_user(db: JobDB = None) -> int:
    """Create demo user for testing."""
    if db is None:
        db = get_db()
    
    demo = db.get_user_by_email("demo@biodockify.com")
    if demo:
        return demo.id
    
    import hashlib
    password_hash = hashlib.sha256(b"demo123").hexdigest()
    return db.create_user("demo@biodockify.com", password_hash)


if __name__ == "__main__":
    db = get_db()
    
    user_id = init_demo_user(db)
    print(f"Demo user ID: {user_id}")
    
    job_uuid = db.create_job(
        user_id=user_id,
        job_name="Test Docking Job",
        compute_mode="auto",
        gpu_used=True,
        gpu_name="NVIDIA RTX 3080"
    )
    print(f"Created job: {job_uuid}")
    
    db.update_job_status(job_uuid, "completed")
    
    jobs = db.get_jobs(user_id)
    print(f"User has {len(jobs)} jobs")
    
    db.close()
    print("Database test complete!")
