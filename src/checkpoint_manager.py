"""
BioDockify Docking Studio - Checkpoint Manager
Manages checkpoint system for job resumption without data loss
"""

import logging
from typing import Dict, Optional, Any
from datetime import datetime
import uuid

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
    
    def __init__(self):
        """Initialize checkpoint manager"""
        self.checkpoints = {}
    
    def save_checkpoint(self, job_id: str, stage: str, data: Dict[str, Any]) -> bool:
        """Save checkpoint"""
        checkpoint_id = f"{job_id}_{stage}"
        timestamp = datetime.now().isoformat()
        
        # Validate checkpoint stage
        if stage not in self.SAFE_CHECKPOINTS:
            logger.warning(f"Invalid checkpoint stage: {stage}")
            return False
        
        self.checkpoints[checkpoint_id] = {
            "stage": stage,
            "data": data,
            "timestamp": timestamp,
            "job_id": job_id
        }
        
        logger.info(f"Checkpoint saved: {checkpoint_id} at {timestamp}")
        return True
    
    def get_checkpoint(self, job_id: str, stage: str) -> Optional[Dict[str, Any]]:
        """Get checkpoint by job ID and stage"""
        checkpoint_id = f"{job_id}_{stage}"
        
        return self.checkpoints.get(checkpoint_id)
    
    def get_previous_checkpoint(self, job_id: str, current_stage: str) -> Optional[Dict[str, Any]]:
        """Get previous safe checkpoint for job resumption"""
        current_index = self.SAFE_CHECKPOINTS.index(current_stage)
        
        if current_index == 0:
            return None  # No previous checkpoint
        
        previous_stage = self.SAFE_CHECKPOINTS[current_index - 1]
        return self.get_checkpoint(job_id, previous_stage)
    
    def clear_checkpoints(self, job_id: str) -> None:
        """Clear all checkpoints for a job"""
        to_clear = [key for key in self.checkpoints.keys() if key.startswith(job_id)]
        
        for key in to_clear:
            del self.checkpoints[key]
        
        logger.info(f"Checkpoints cleared for job: {job_id}")
