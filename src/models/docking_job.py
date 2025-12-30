"""
BioDockify Docking Studio - Docking Job Model
Data model for docking jobs
"""

from sqlalchemy import Column, Integer, String, Float, Boolean, DateTime, Text, JSON
from sqlalchemy.ext.declarative import declarative_base
from datetime import datetime
import uuid
import json

Base = declarative_base()

class DockingJob(Base):
    """Docking job model"""
    __tablename__ = "docking_jobs"
    
    id = Column(String, primary_key=True, default=lambda: str(uuid.uuid4()))
    status = Column(String, default="PENDING")
    receptor_file = Column(String, nullable=False)
    ligand_file = Column(String, nullable=False)
    parameters = Column(JSON, nullable=False)
    
    # Vina parameters
    center_x = Column(Float, default=0.0)
    center_y = Column(Float, default=0.0)
    center_z = Column(Float, default=0.0)
    size_x = Column(Float, default=20.0)
    size_y = Column(Float, default=20.0)
    size_z = Column(Float, default=20.0)
    exhaustiveness = Column(Integer, default=8)
    num_modes = Column(Integer, default=9)
    
    # Results
    binding_energy = Column(Float, nullable=True)
    num_modes = Column(Integer, nullable=True)
    vina_output = Column(Text, nullable=True)
    
    # Analysis results
    interactions = Column(JSON, nullable=True)
    descriptors = Column(JSON, nullable=True)
    confidence_score = Column(Integer, nullable=True)
    
    # Agent Zero
    agent_zero_failures = Column(JSON, nullable=True)
    agent_zero_repairs = Column(JSON, nullable=True)
    
    # Checkpoints
    current_checkpoint = Column(String, nullable=True)
    
    # Timestamps
    created_at = Column(DateTime, default=datetime.now)
    updated_at = Column(DateTime, default=datetime.now)
    completed_at = Column(DateTime, nullable=True)
    
    def to_dict(self) -> dict:
        """Convert model to dictionary"""
        return {
            "id": self.id,
            "status": self.status,
            "receptor_file": self.receptor_file,
            "ligand_file": self.ligand_file,
            "parameters": json.loads(self.parameters) if self.parameters else {},
            "center_x": self.center_x,
            "center_y": self.center_y,
            "center_z": self.center_z,
            "size_x": self.size_x,
            "size_y": self.size_y,
            "size_z": self.size_z,
            "exhaustiveness": self.exhaustiveness,
            "num_modes": self.num_modes,
            "binding_energy": self.binding_energy,
            "num_modes": self.num_modes,
            "vina_output": self.vina_output,
            "interactions": json.loads(self.interactions) if self.interactions else {},
            "descriptors": json.loads(self.descriptors) if self.descriptors else {},
            "confidence_score": self.confidence_score,
            "agent_zero_failures": json.loads(self.agent_zero_failures) if self.agent_zero_failures else [],
            "agent_zero_repairs": json.loads(self.agent_zero_repairs) if self.agent_zero_repairs else [],
            "current_checkpoint": self.current_checkpoint,
            "created_at": self.created_at.isoformat(),
            "updated_at": self.updated_at.isoformat(),
            "completed_at": self.completed_at.isoformat() if self.completed_at else None
        }
