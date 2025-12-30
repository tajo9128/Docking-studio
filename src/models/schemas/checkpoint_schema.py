"""
BioDockify Docking Studio - Checkpoint Schema
SQLAlchemy schema definition for checkpoints
"""

from sqlalchemy import Column, String, JSON, DateTime
from sqlalchemy.orm import declarative_base
from models import Base

class Checkpoint(Base):
    """Checkpoint SQLAlchemy model"""
    __tablename__ = "checkpoints"
    __table_args__ = {"extend_existing": True}
    
    id = Column(String, primary_key=True)
    job_id = Column(String, nullable=False, index=True)
    stage = Column(String, nullable=False)
    data = Column(JSON, nullable=False)
    created_at = Column(DateTime, default="CURRENT_TIMESTAMP")
