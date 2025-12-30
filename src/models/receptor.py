"""
BioDockify Docking Studio - Receptor Model
Data model for receptor molecules
"""

from sqlalchemy import Column, String, Float, Integer, JSON
from sqlalchemy.ext.declarative import declarative_base
import json

Base = declarative_base()

class Receptor(Base):
    """Receptor model"""
    __tablename__ = "receptors"
    
    id = Column(String, primary_key=True)
    filename = Column(String, nullable=False, unique=True)
    file_path = Column(String, nullable=False)
    file_size = Column(Integer)
    file_hash = Column(String)
    
    # Structure information
    num_atoms = Column(Integer, default=0)
    num_residues = Column(Integer, default=0)
    protein_type = Column(String, default="unknown")
    molecular_weight = Column(Float, default=0.0)
    
    # Validation status
    validation_status = Column(String, default="PENDING")
    validation_errors = Column(JSON, nullable=True)
    
    # Analysis information
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
            "num_residues": self.num_residues,
            "protein_type": self.protein_type,
            "molecular_weight": self.molecular_weight,
            "validation_status": self.validation_status,
            "validation_errors": json.loads(self.validation_errors) if self.validation_errors else {},
            "has_metal": self.has_metal,
            "metal_atoms": json.loads(self.metal_atoms) if self.metal_atoms else {},
            "uploaded_at": self.uploaded_at,
            "validated_at": self.validated_at
        }
