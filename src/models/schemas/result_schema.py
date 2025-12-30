"""
BioDockify Docking Studio - Result Schema
SQLAlchemy schema definition for docking results
"""

from sqlalchemy import Column, String, Float, Integer, JSON, Boolean, DateTime
from sqlalchemy.orm import declarative_base
from models import Base

class Result(Base):
    """Result SQLAlchemy model"""
    __tablename__ = "results"
    __table_args__ = {"extend_existing": True}
    
    id = Column(String, primary_key=True)
    job_id = Column(String, nullable=False, index=True)
    binding_energy = Column(Float, nullable=False)
    num_modes = Column(Integer, default=9)
    poses = Column(JSON, nullable=False)
    interactions = Column(JSON, nullable=False)
    descriptors = Column(JSON, nullable=False)
    confidence_score = Column(Integer, nullable=False)
    failures_detected = Column(JSON, nullable=True)
    repairs_attempted = Column(JSON, nullable=True)
    confidence_adjustment = Column(Integer, default=0)
    rmsd_ub = Column(Float, nullable=True)
    rmsd_lb = Column(Float, nullable=True)
    best_pose = Column(Integer, nullable=False)
    created_at = Column(DateTime, default="CURRENT_TIMESTAMP")
    analyzed_at = Column(DateTime, nullable=True)
