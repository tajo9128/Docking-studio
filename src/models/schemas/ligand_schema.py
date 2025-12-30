"""
BioDockify Docking Studio - Ligand Schema
SQLAlchemy schema definition for ligands
"""

from sqlalchemy import Column, String, Float, Integer, JSON, Boolean, DateTime
from sqlalchemy.orm import declarative_base
from models import Base

class Ligand(Base):
    """Ligand SQLAlchemy model"""
    __tablename__ = "ligands"
    __table_args__ = {"extend_existing": True}
    
    id = Column(String, primary_key=True)
    filename = Column(String, nullable=False, unique=True)
    file_path = Column(String, nullable=False)
    file_size = Column(Integer)
    file_hash = Column(String)
    num_atoms = Column(Integer, default=0)
    molecular_weight = Column(Float, default=0.0)
    smiles = Column(String, nullable=True)
    validation_status = Column(String, default="PENDING")
    validation_errors = Column(JSON, nullable=True)
    is_3d = Column(Boolean, default=True)
    num_rotatable_bonds = Column(Integer, default=0)
    has_metal = Column(Boolean, default=False)
    metal_atoms = Column(JSON, nullable=True)
    uploaded_at = Column(DateTime, default="CURRENT_TIMESTAMP")
    validated_at = Column(DateTime, nullable=True)
