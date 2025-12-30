"""
BioDockify Docking Studio - Checkpoint Model
Data model for checkpoint data
"""

from sqlalchemy import Column, String, JSON, DateTime
from sqlalchemy.ext.declarative import declarative_base
import json
from datetime import datetime

Base = declarative_base()

class Checkpoint(Base):
    """Checkpoint model"""
    __tablename__ = "checkpoints"
    
    id = Column(String, primary_key=True)
    job_id = Column(String, nullable=False, index=True)  # Index for efficient queries
    stage = Column(String, nullable=False)  # Checkpoint stage
    data = Column(JSON, nullable=False)  # Checkpoint data
    
    # Metadata
    created_at = Column(DateTime, default=datetime.now)
    
    def to_dict(self) -> dict:
        """Convert model to dictionary"""
        return {
            "id": self.id,
            "job_id": self.job_id,
            "stage": self.stage,
            "data": json.loads(self.data) if self.data else {},
            "created_at": self.created_at.isoformat()
        }
