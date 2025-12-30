"""
BioDockify Docking Studio - Result Model
Data model for docking results
"""

from sqlalchemy import Column, String, Float, Integer, JSON, DateTime
from sqlalchemy.ext.declarative import declarative_base
import json
from datetime import datetime

Base = declarative_base()

class Result(Base):
    """Result model"""
    __tablename__ = "results"
    
    id = Column(String, primary_key=True)
    job_id = Column(String, nullable=False, index=True)  # Index for efficient queries
    binding_energy = Column(Float, nullable=False)
    num_modes = Column(Integer, default=9)
    
    # Poses
    poses = Column(JSON, nullable=False)  # Store all poses
    
    # Analysis results
    interactions = Column(JSON, nullable=False)  # Store interaction data
    descriptors = Column(JSON, nullable=False)  # Store descriptor data
    confidence_score = Column(Integer, nullable=False)
    
    # Agent Zero data
    failures_detected = Column(JSON, nullable=True)
    repairs_attempted = Column(JSON, nullable=True)
    confidence_adjustment = Column(Integer, default=0)  # Penalty points
    
    # Quality metrics
    rmsd_ub = Column(Float, nullable=True)  # Root-mean-square deviation
    rmsd_lb = Column(Float, nullable=True)  # Best pose RMSD
    best_pose = Column(Integer, nullable=False)  # Index of best pose
    
    # Timestamps
    created_at = Column(DateTime, default=datetime.now)
    analyzed_at = Column(DateTime, nullable=True)
    
    def to_dict(self) -> dict:
        """Convert model to dictionary"""
        return {
            "id": self.id,
            "job_id": self.job_id,
            "binding_energy": self.binding_energy,
            "num_modes": self.num_modes,
            "poses": json.loads(self.poses) if self.poses else [],
            "interactions": json.loads(self.interactions) if self.interactions else {},
            "descriptors": json.loads(self.descriptors) if self.descriptors else {},
            "confidence_score": self.confidence_score,
            "failures_detected": json.loads(self.failures_detected) if self.failures_detected else [],
            "repairs_attempted": json.loads(self.repairs_attempted) if self.repairs_attempted else [],
            "confidence_adjustment": self.confidence_adjustment,
            "rmsd_ub": self.rmsd_ub,
            "rmsd_lb": self.rmsd_lb,
            "best_pose": self.best_pose,
            "created_at": self.created_at.isoformat(),
            "analyzed_at": self.analyzed_at.isoformat() if self.analyzed_at else None
        }
