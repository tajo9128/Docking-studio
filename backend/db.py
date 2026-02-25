import sqlite3
import os
from datetime import datetime
from typing import Optional, List, Dict, Any

DB_PATH = os.path.join(os.path.dirname(__file__), "storage", "jobs.db")

def init_db():
    """Initialize SQLite database for job storage"""
    os.makedirs(os.path.dirname(DB_PATH), exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()

    cur.execute("""
        CREATE TABLE IF NOT EXISTS jobs (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            job_uuid TEXT UNIQUE,
            job_name TEXT,
            receptor_file TEXT,
            ligand_file TEXT,
            status TEXT DEFAULT 'pending',
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            completed_at DATETIME,
            binding_energy REAL,
            confidence_score REAL,
            engine TEXT
        )
    """)

    cur.execute("""
        CREATE TABLE IF NOT EXISTS docking_results (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            job_uuid TEXT,
            pose_id INTEGER,
            ligand_name TEXT,
            vina_score REAL,
            gnina_score REAL,
            rf_score REAL,
            consensus REAL,
            pdb_data TEXT,
            FOREIGN KEY (job_uuid) REFERENCES jobs(job_uuid)
        )
    """)

    cur.execute("""
        CREATE TABLE IF NOT EXISTS interactions (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            job_uuid TEXT,
            pose_id INTEGER,
            interaction_type TEXT,
            atom_a TEXT,
            atom_b TEXT,
            distance REAL,
            FOREIGN KEY (job_uuid) REFERENCES jobs(job_uuid)
        )
    """)

    conn.commit()
    conn.close()


def create_job(job_uuid: str, job_name: str, receptor_file: str, ligand_file: str, engine: str = "vina") -> bool:
    """Create a new job entry"""
    try:
        conn = sqlite3.connect(DB_PATH)
        cur = conn.cursor()
        cur.execute("""
            INSERT INTO jobs (job_uuid, job_name, receptor_file, ligand_file, engine, status)
            VALUES (?, ?, ?, ?, ?, 'pending')
        """, (job_uuid, job_name, receptor_file, ligand_file, engine))
        conn.commit()
        conn.close()
        return True
    except Exception as e:
        print(f"Error creating job: {e}")
        return False


def update_job_status(job_uuid: str, status: str, binding_energy: Optional[float] = None, 
                     confidence_score: Optional[float] = None) -> bool:
    """Update job status"""
    try:
        conn = sqlite3.connect(DB_PATH)
        cur = conn.cursor()
        
        completed_at = "CURRENT_TIMESTAMP" if status in ['completed', 'failed', 'cancelled'] else "NULL"
        
        query = f"""
            UPDATE jobs 
            SET status = ?, 
                completed_at = {completed_at}
        """
        params = [status]
        
        if binding_energy is not None:
            query += ", binding_energy = ?"
            params.append(binding_energy)
        
        if confidence_score is not None:
            query += ", confidence_score = ?"
            params.append(confidence_score)
        
        query += " WHERE job_uuid = ?"
        params.append(job_uuid)
        
        cur.execute(query, params)
        conn.commit()
        conn.close()
        return True
    except Exception as e:
        print(f"Error updating job: {e}")
        return False


def add_docking_result(job_uuid: str, pose_id: int, ligand_name: str, 
                      vina_score: Optional[float] = None, gnina_score: Optional[float] = None,
                      rf_score: Optional[float] = None, pdb_data: Optional[str] = None) -> bool:
    """Add docking result for a pose"""
    try:
        conn = sqlite3.connect(DB_PATH)
        cur = conn.cursor()
        
        consensus = None
        scores = [s for s in [vina_score, gnina_score, rf_score] if s is not None]
        if scores:
            consensus = sum(scores) / len(scores)
        
        cur.execute("""
            INSERT INTO docking_results (job_uuid, pose_id, ligand_name, vina_score, gnina_score, rf_score, consensus, pdb_data)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, (job_uuid, pose_id, ligand_name, vina_score, gnina_score, rf_score, consensus, pdb_data))
        
        conn.commit()
        conn.close()
        return True
    except Exception as e:
        print(f"Error adding docking result: {e}")
        return False


def add_interaction(job_uuid: str, pose_id: int, interaction_type: str, atom_a: str, 
                   atom_b: str, distance: float) -> bool:
    """Add interaction data"""
    try:
        conn = sqlite3.connect(DB_PATH)
        cur = conn.cursor()
        cur.execute("""
            INSERT INTO interactions (job_uuid, pose_id, interaction_type, atom_a, atom_b, distance)
            VALUES (?, ?, ?, ?, ?, ?)
        """, (job_uuid, pose_id, interaction_type, atom_a, atom_b, distance))
        conn.commit()
        conn.close()
        return True
    except Exception as e:
        print(f"Error adding interaction: {e}")
        return False


def get_job(job_uuid: str) -> Optional[Dict[str, Any]]:
    """Get job by UUID"""
    try:
        conn = sqlite3.connect(DB_PATH)
        cur = conn.cursor()
        cur.execute("SELECT * FROM jobs WHERE job_uuid = ?", (job_uuid,))
        row = cur.fetchone()
        conn.close()
        
        if row:
            columns = ['id', 'job_uuid', 'job_name', 'receptor_file', 'ligand_file', 'status', 
                      'created_at', 'completed_at', 'binding_energy', 'confidence_score', 'engine']
            return dict(zip(columns, row))
        return None
    except Exception as e:
        print(f"Error getting job: {e}")
        return None


def get_all_jobs(limit: int = 50) -> List[Dict[str, Any]]:
    """Get all jobs"""
    try:
        conn = sqlite3.connect(DB_PATH)
        cur = conn.cursor()
        cur.execute("SELECT * FROM jobs ORDER BY created_at DESC LIMIT ?", (limit,))
        rows = cur.fetchall()
        conn.close()
        
        columns = ['id', 'job_uuid', 'job_name', 'receptor_file', 'ligand_file', 'status',
                  'created_at', 'completed_at', 'binding_energy', 'confidence_score', 'engine']
        return [dict(zip(columns, row)) for row in rows]
    except Exception as e:
        print(f"Error getting jobs: {e}")
        return []


def get_docking_results(job_uuid: str) -> List[Dict[str, Any]]:
    """Get docking results for a job"""
    try:
        conn = sqlite3.connect(DB_PATH)
        cur = conn.cursor()
        cur.execute("SELECT * FROM docking_results WHERE job_uuid = ? ORDER BY pose_id", (job_uuid,))
        rows = cur.fetchall()
        conn.close()
        
        columns = ['id', 'job_uuid', 'pose_id', 'ligand_name', 'vina_score', 'gnina_score', 
                  'rf_score', 'consensus', 'pdb_data']
        return [dict(zip(columns, row)) for row in rows]
    except Exception as e:
        print(f"Error getting results: {e}")
        return []


def get_interactions(job_uuid: str, pose_id: Optional[int] = None) -> List[Dict[str, Any]]:
    """Get interactions for a job"""
    try:
        conn = sqlite3.connect(DB_PATH)
        cur = conn.cursor()
        
        if pose_id is not None:
            cur.execute("""
                SELECT * FROM interactions WHERE job_uuid = ? AND pose_id = ?
                ORDER BY distance
            """, (job_uuid, pose_id))
        else:
            cur.execute("SELECT * FROM interactions WHERE job_uuid = ? ORDER BY pose_id, distance", (job_uuid,))
        
        rows = cur.fetchall()
        conn.close()
        
        columns = ['id', 'job_uuid', 'pose_id', 'interaction_type', 'atom_a', 'atom_b', 'distance']
        return [dict(zip(columns, row)) for row in rows]
    except Exception as e:
        print(f"Error getting interactions: {e}")
        return []


def delete_job(job_uuid: str) -> bool:
    """Delete a job and its results"""
    try:
        conn = sqlite3.connect(DB_PATH)
        cur = conn.cursor()
        cur.execute("DELETE FROM interactions WHERE job_uuid = ?", (job_uuid,))
        cur.execute("DELETE FROM docking_results WHERE job_uuid = ?", (job_uuid,))
        cur.execute("DELETE FROM jobs WHERE job_uuid = ?", (job_uuid,))
        conn.commit()
        conn.close()
        return True
    except Exception as e:
        print(f"Error deleting job: {e}")
        return False


if __name__ == "__main__":
    init_db()
    print("Database initialized successfully!")
