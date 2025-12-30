"""
BioDockify Docking Studio - Ligand Model
Data model for ligand molecules
"""

from sqlalchemy import Column, String, Float, Integer, JSON, Boolean
from sqlalchemy.ext.declarative import declarative_base
import json

Base = declarative_base()

class Ligand(Base):
    """Ligand model"""
    __tablename__ = "ligands"
    
    id = Column(String, primary_key=True)
    filename = Column(String, nullable=False, unique=True)
    file_path = Column(String, nullable=False)
    file_size = Column(Integer)
    file_hash = Column(String)
    
    # Structure information
    num_atoms = Column(Integer, default=0)
    molecular_weight = Column(Float, default=0.0)
    smiles = Column(String, nullable=True)
    
    # Validation status
    validation_status = Column(String, default="PENDING")
    validation_errors = Column(JSON, nullable=True)
    
    # Analysis information
    is_3d = Column(Boolean, default=True)
    num_rotatable_bonds = Column(Integer, default=0)
    has_metal = Column(Boolean, default=False)
    metal_atoms = Column(JSON, nullable=True)
    
    # Timestamps
    uploaded_at = Column(String, nullable=False)
    validated_at = Column(String, nullable=True)
    
    def to_dict(self) -> dict:
        """Convert model to dictionary"""
        return {
            "id": self.id,
            "filename": self.filename,
            "file_path": self.file_path,
            "file_size": self.file_size,
            "file_hash": self.file_hash,
            "num_atoms": self.num_atoms,
            "molecular_weight": self.molecular_weight,
            "smiles": self.smiles,
            "validation_status": self.validation_status,
            "validation_errors": json.loads(self.validation_errors) if self.validation_errors else {},
            "is_3d": self.is_3d,
            "num_rotatable_bonds": self.num_rotatable_bonds,
            "has_metal": self.has_metal,
            "metal_atoms": json.loads(self.metal_atoms) if self.metal_atoms else {},
            "uploaded_at": self.uploaded_at,
            "validated_at": self.validated_at
        }
