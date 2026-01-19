"""Pydantic Schemas for Request/Response Validation"""

from pydantic import BaseModel, Field, field_validator
from typing import Optional, List, Dict, Any
from datetime import datetime
from enum import Enum

# Error Schemas
class ErrorType(str, Enum):
    validation = "validation_error"
    network = "network_error"
    system = "system_error"
    http = "http_error"

class ErrorResponse(BaseModel):
    """Error response schema"""
    
    type: str = Field(..., description="Type of error")
    code: str = Field(..., description="Error code")
    message: str = Field(..., description="Human-readable error message")
    details: Optional[str] = Field(None, description="Additional error details (development only)")
    correlation_id: Optional[str] = Field(None, description="Correlation ID for tracing")
    timestamp: str = Field(default_factory=lambda: datetime.now().isoformat(), description="Error timestamp")
    
    class Config:
        json_schema_extra = {
            "example": {
                "type": "validation_error",
                "code": "INVALID_FILE_TYPE",
                "message": "Unsupported file type",
                "details": "The file you uploaded is not supported.",
                "correlation_id": "123e4567-e89b-12d3-a456-426614174000",
                "timestamp": "2024-01-01T00:00:00Z",
            }
        }

class ValidationError(BaseModel):
    """Validation error"""
    field: str
    message: str

# Structure Schemas
class StructureUploadResponse(BaseModel):
    """Structure upload response"""
    
    structure_id: str = Field(..., description="Unique structure ID")
    file_name: str = Field(..., description="Uploaded file name")
    file_size: int = Field(..., description="File size in bytes")
    file_hash: str = Field(..., description="SHA-256 hash of file content")
    content_type: str = Field(..., description="Content type of file")
    stage: str = Field(default="upload", description="Current stage (upload/parse/analyze)")
    timestamp: str = Field(default_factory=lambda: datetime.now().isoformat(), description="Upload timestamp")

class AtomModel(BaseModel):
    """Atom model"""
    
    index: int = Field(..., description="0-indexed atom index")
    serial: int = Field(..., description="PDB serial number")
    name: str = Field(..., description="Atom name")
    alt_loc: str = Field(default="", description="Alternate location indicator")
    res_name: str = Field(..., description="Residue name")
    chain_id: str = Field(default="", description="Chain identifier")
    res_seq: int = Field(..., description="Residue sequence number")
    i_code: str = Field(default="", description="Insertion code")
    x: float = Field(..., description="X coordinate in Angstroms")
    y: float = Field(..., description="Y coordinate in Angstroms")
    z: float = Field(..., description="Z coordinate in Angstroms")
    occupancy: float = Field(default=1.0, description="Occupancy")
    temp_factor: float = Field(default=0.0, description="Temperature factor")
    element: str = Field(..., description="Element symbol")
    charge: float = Field(default=0.0, description="Formal charge")

class BondModel(BaseModel):
    """Bond model"""
    
    atom1_index: int = Field(..., description="Index of first atom")
    atom2_index: int = Field(..., description="Index of second atom")
    type: str = Field(..., description="Bond type (single/double/triple/aromatic)")
    order: int = Field(..., description="Bond order (1/2/3/1.5)")
    distance: float = Field(..., description="Bond distance in Angstroms")

class StructureMetadata(BaseModel):
    """Structure metadata"""
    
    file_name: str = Field(..., description="File name")
    file_size: int = Field(..., description="File size in bytes")
    atom_count: int = Field(..., description="Total atom count")
    bond_count: int = Field(..., description="Total bond count")
    chain_count: int = Field(..., description="Total chain count")
    model_count: int = Field(..., description="Total model count")
    title: Optional[str] = Field(None, description="Structure title (from PDB)")
    experimental_technique: Optional[str] = Field(None, description="Experimental technique (from PDB)")
    resolution: Optional[float] = Field(None, description="Resolution in Angstroms")
    warnings: List[str] = Field(default_factory=list, description="Parser warnings")

class StructureParseResponse(BaseModel):
    """Structure parse response"""
    
    structure_id: str = Field(..., description="Unique structure ID")
    metadata: StructureMetadata = Field(..., description="Structure metadata")
    atoms: List[AtomModel] = Field(..., description="List of atoms")
    bonds: Optional[List[BondModel]] = Field(None, description="List of bonds (if parsed)")
    stage: str = Field(default="parsed", description="Current stage (upload/parse/analyze)")
    timestamp: str = Field(default_factory=lambda: datetime.now().isoformat(), description="Parse timestamp")

# Analysis Schemas
class HydrogenBond(BaseModel):
    """Hydrogen bond model"""
    
    atom1_index: int = Field(..., description="Index of first atom")
    atom2_index: int = Field(..., description="Index of second atom")
    distance: float = Field(..., ge=1.5, le=2.5, description="Hydrogen bond distance in Angstroms (1.5-2.5)")
    angle: Optional[float] = Field(None, ge=120, le=180, description="D-H...A angle (120-180)")
    atom1_residue: str = Field(..., description="Residue of first atom")
    atom1_residue_seq: int = Field(..., description="Residue sequence number of first atom")
    atom2_residue: str = Field(..., description="Residue of second atom")
    atom2_residue_seq: int = Field(..., description="Residue sequence number of second atom")
    confidence: float = Field(..., ge=0.0, le=1.0, description="Confidence score (0-1)")
    is_predicted: bool = Field(default=False, description="Whether this is a predicted bond")

class VDWContact(BaseModel):
    """Van der Waals contact model"""
    
    atom1_index: int = Field(..., description="Index of first atom")
    atom2_index: int = Field(..., description="Index of second atom")
    distance: float = Field(..., ge=0.0, le=10.0, description="VdW contact distance in Angstroms")
    atom1_residue: str = Field(..., description="Residue of first atom")
    atom1_residue_seq: int = Field(..., description="Residue sequence number of first atom")
    atom2_residue: str = Field(..., description="Residue of second atom")
    atom2_residue_seq: int = Field(..., description="Residue sequence number of second atom")
    confidence: float = Field(..., ge=0.0, le=1.0, description="Confidence score (0-1)")
    is_predicted: bool = Field(default=False, description="Whether this is a predicted contact")

class SaltBridge(BaseModel):
    """Salt bridge model"""
    
    atom1_index: int = Field(..., description="Index of first atom")
    atom2_index: int = Field(..., description="Index of second atom")
    distance: float = Field(..., ge=0.0, le=4.0, description="Salt bridge distance in Angstroms (0-4.0)")
    atom1_residue: str = Field(..., description="Residue of first atom")
    atom1_residue_seq: int = Field(..., description="Residue sequence number of first atom")
    atom2_residue: str = Field(..., description="Residue of second atom")
    atom2_residue_seq: int = Field(..., description="Residue sequence number of second atom")
    confidence: float = Field(..., ge=0.0, le=1.0, description="Confidence score (0-1)")
    is_predicted: bool = Field(default=False, description="Whether this is a predicted bridge")

class AnalysisMetadata(BaseModel):
    """Analysis metadata"""
    
    processing_time_ms: float = Field(..., description="Processing time in milliseconds")
    atom_count: int = Field(..., description="Total atom count")
    bond_count: int = Field(..., description="Total bond count")
    algorithm: str = Field(..., description="Algorithm used for analysis")
    thresholds: Dict[str, Any] = Field(..., description="Thresholds used in analysis")

class AnalysisResponse(BaseModel):
    """Analysis response"""
    
    structure_id: str = Field(..., description="Unique structure ID")
    hydrogen_bonds: List[HydrogenBond] = Field(default_factory=list, description="List of hydrogen bonds")
    vdw_contacts: List[VDWContact] = Field(default_factory=list, description="List of VdW contacts")
    salt_bridges: List[SaltBridge] = Field(default_factory=list, description="List of salt bridges")
    total_interactions: int = Field(..., description="Total interaction count")
    metadata: AnalysisMetadata = Field(..., description="Analysis metadata")
    stage: str = Field(default="analyzed", description="Current stage (upload/parse/analyze)")
    timestamp: str = Field(default_factory=lambda: datetime.now().isoformat(), description="Analysis timestamp")

# Export Schemas
class SnapshotMetadata(BaseModel):
    """Snapshot metadata"""
    
    timestamp: str = Field(..., description="Snapshot timestamp")
    structure: Dict[str, Any] = Field(..., description="Structure metadata")
    visualization: Dict[str, Any] = Field(..., description="Visualization state")
    interactions: Dict[str, Any] = Field(..., description="Interaction state")
    export_format: str = Field(..., description="Export format (png/jpeg/webp)")
    resolution: int = Field(..., description="Export resolution")
    dpi: int = Field(..., description="Export DPI")
    background_color: Optional[str] = Field(None, description="Background color")
    visible_atoms: List[int] = Field(default_factory=list, description="Indices of visible atoms")
    hidden_atoms: List[int] = Field(default_factory=list, description="Indices of hidden atoms")

class ExportResponse(BaseModel):
    """Export response"""
    
    export_id: str = Field(..., description="Unique export ID")
    file_name: str = Field(..., description="Exported file name")
    file_size: int = Field(..., description="Exported file size in bytes")
    content_type: str = Field(..., description="Content type of exported file")
    metadata: SnapshotMetadata = Field(..., description="Snapshot metadata")
    stage: str = Field(default="exported", description="Current stage (upload/parse/analyze/export)")
    timestamp: str = Field(default_factory=lambda: datetime.now().isoformat(), description="Export timestamp")

# Upload Schemas
class AnalysisRequest(BaseModel):
    """Analysis request"""
    
    structure_id: str = Field(..., description="Unique structure ID")
    options: Optional[Dict[str, Any]] = Field(None, description="Analysis options")
