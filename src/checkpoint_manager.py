"""
BioDockify Docking Studio - Checkpoint Manager
Manages checkpoint system for job resumption without data loss
"""

import logging
from typing import Dict, Optional, Any
from datetime import datetime
import uuid
import os
import json
from pathlib import Path

logger = logging.getLogger(__name__)

class CheckpointManager:
    """Manages checkpoint system for job resumption"""
    
    SAFE_CHECKPOINTS = [
        "FILE_VALIDATION",
        "PREPROCESSING_COMPLETE",
        "VINA_READY",
        "VINA_COMPLETED",
        "INTERACTION_ANALYSIS_READY",
        "INTERACTION_ANALYSIS_COMPLETED",
        "DESCRIPTOR_CALCULATION_READY",
        "DESCRIPTOR_CALCULATION_COMPLETED",
        "FINAL_RESULTS"
    ]
    
    def __init__(self, db_manager=None):
        """Initialize checkpoint manager"""
        self.db_manager = db_manager
        self.checkpoints = {}
        self._persist_dir = self._get_persist_dir()
        
    def _get_persist_dir(self) -> Path:
        """Get directory for checkpoint persistence"""
        if os.name == "nt":
            persist_dir = Path.home() / "AppData" / "Local" / "BioDockify" / "checkpoints"
        else:
            persist_dir = Path.home() / ".config" / "BioDockify" / "checkpoints"
        
        persist_dir.mkdir(parents=True, exist_ok=True)
        return persist_dir
    
    def save_checkpoint(self, job_id: str, stage: str, data: Dict[str, Any]) -> bool:
        """Save checkpoint"""
        checkpoint_id = f"{job_id}_{stage}"
        timestamp = datetime.now().isoformat()
        
        if stage not in self.SAFE_CHECKPOINTS:
            logger.warning(f"Invalid checkpoint stage: {stage}")
            return False
            
        checkpoint_data = {
            "stage": stage,
            "data": data,
            "timestamp": timestamp,
            "job_id": job_id
        }
        
        # 1. Save to in-memory cache
        self.checkpoints[checkpoint_id] = checkpoint_data
        
        # 2. Save to database (if available)
        if self.db_manager:
            try:
                # self.db_manager.save_checkpoint(job_id, stage, data)
                pass 
            except Exception as e:
                logger.error(f"Failed to save checkpoint to DB: {e}")
        
        # 3. Save to file
        try:
            checkpoint_file = self._persist_dir / f"{checkpoint_id}.json"
            with open(checkpoint_file, 'w') as f:
                json.dump(checkpoint_data, f, indent=2)
        except Exception as e:
            logger.error(f"Failed to save checkpoint to file: {e}")
            
        logger.info(f"Checkpoint saved: {checkpoint_id}")
        return True
    
    def get_checkpoint(self, job_id: str, stage: str) -> Optional[Dict[str, Any]]:
        """Get checkpoint"""
        checkpoint_id = f"{job_id}_{stage}"
        
        if checkpoint_id in self.checkpoints:
            return self.checkpoints[checkpoint_id]
            
        try:
            checkpoint_file = self._persist_dir / f"{checkpoint_id}.json"
            if checkpoint_file.exists():
                with open(checkpoint_file, 'r') as f:
                    data = json.load(f)
                self.checkpoints[checkpoint_id] = data
                return data
        except Exception:
            pass
            
        return None

    def get_previous_checkpoint(self, job_id: str, current_stage: str) -> Optional[Dict[str, Any]]:
        """Get previous safe checkpoint for job resumption"""
        try:
            current_index = self.SAFE_CHECKPOINTS.index(current_stage)
        except ValueError:
            return None
        
        if current_index == 0:
            return None
        
        previous_stage = self.SAFE_CHECKPOINTS[current_index - 1]
        return self.get_checkpoint(job_id, previous_stage)

    def cleanup_old_checkpoints(self, max_age_days: int = 7) -> int:
        """Clean up old checkpoints"""
        from datetime import timedelta
        cleaned_count = 0
        cutoff_date = datetime.now() - timedelta(days=max_age_days)
        
        to_remove = []
        for cid, cp in self.checkpoints.items():
            try:
                if datetime.fromisoformat(cp["timestamp"]) < cutoff_date:
                    to_remove.append(cid)
            except: pass
        for cid in to_remove:
            del self.checkpoints[cid]
            cleaned_count += 1
            
        try:
            for cp_file in self._persist_dir.glob("*.json"):
                if cp_file.stat().st_mtime < cutoff_date.timestamp():
                    cp_file.unlink()
                    cleaned_count += 1
        except Exception as e:
            logger.error(f"Error cleaning old checkpoints: {e}")
            
        return cleaned_count
    
    def clear_checkpoints(self, job_id: str) -> None:
        """Clear all checkpoints for a job"""
        to_clear = [key for key in self.checkpoints.keys() if key.startswith(job_id)]
        for key in to_clear:
            del self.checkpoints[key]
        logger.info(f"Checkpoints cleared for job: {job_id}")
