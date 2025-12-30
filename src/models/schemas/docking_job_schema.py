"""
BioDockify Docking Studio - Docking Job Schema
SQLAlchemy schema definition for docking jobs
"""

from sqlalchemy import Column, Integer, String, Float, JSON, Boolean, DateTime, Text
from sqlalchemy.orm import declarative_base
from models import Base

class DockingJob(Base):
    """Docking job SQLAlchemy model"""
    __tablename__ = "docking_jobs"
    __table_args__ = {"extend_existing": True}
    
    id = Column(String, primary_key=True)
    status = Column(String, default="PENDING")
    receptor_file = Column(String, nullable=False)
    ligand_file = Column(String, nullable=False)
    parameters = Column(JSON, nullable=False)
    
    center_x = Column(Float, default=0.0)
    center_y = Column(Float, default=0.0)
    center_z = Column(Float, default=0.0)
    size_x = Column(Float, default=20.0)
    size_y = Column(Float, default=20.0)
    size_z = Column(Float, default=20.0)
    exhaustiveness = Column(Integer, default=8)
    num_modes = Column(Integer, default=9)
    
    binding_energy = Column(Float, nullable=True)
    num_modes = Column(Integer, nullable=True)
    vina_output = Column(Text, nullable=True)
    
    interactions = Column(JSON, nullable=True)
    descriptors = Column(JSON, nullable=True)
    confidence_score = Column(Integer, nullable=True)
    
    agent_zero_failures = Column(JSON, nullable=True)
    agent_zero_repairs = Column(JSON, nullable=True)
    current_checkpoint = Column(String, nullable=True)
    
    created_at = Column(DateTime, default="CURRENT_TIMESTAMP")
    updated_at = Column(DateTime, default="CURRENT_TIMESTAMP")
    completed_at = Column(DateTime, nullable=True)
