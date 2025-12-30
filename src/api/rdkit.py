"""
BioDockify Docking Studio - RDKit API Router
Handles RDKit descriptor calculation operations
"""

from fastapi import APIRouter, HTTPException, status
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from typing import Dict, Any
import logging

from src.rdkit_calculator import RDKitCalculator

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/rdkit")

rdkit_calculator = RDKitCalculator()

class DescriptorResult(BaseModel):
    """Descriptor result model"""
    molecular_weight: float
    log_p: float
    tpsa: float
    num_rotatable_bonds: int
    num_hbd: int
    num_hba: int
    num_aromatic_rings: int
    lipinski_violations: int
    drug_likeness: str
    veber_rule: str
    egan_rule: str

class FullDescriptors(BaseModel):
    """Full descriptors model"""
    all_descriptors: Dict[str, Any]
    lipinski_rule: bool
    veber_rule: bool
    egan_rule: bool
    drug_likeness: str

@router.post("/calculate", response_model=FullDescriptors)
async def calculate_descriptors(ligand_file: str, job_id: str):
    """Calculate RDKit molecular descriptors"""
    
    results = rdkit_calculator.calculate_descriptors(
        ligand_file=ligand_file,
        job_id=job_id
    )
    
    if results["status"] == "COMPLETED":
        descriptors = results["descriptors"]
        
        all_descriptors = {
            "molecular_weight": descriptors["molecular_weight"],
            "log_p": descriptors["log_p"],
            "tpsa": descriptors["tpsa"],
            "num_rotatable_bonds": descriptors["num_rotatable_bonds"],
            "num_hbd": descriptors["num_hbd"],
            "num_hba": descriptors["num_hba"],
            "num_aromatic_rings": descriptors["num_aromatic_rings"]
        }
        
        return FullDescriptors(
            all_descriptors=all_descriptors,
            lipinski_rule=descriptors["lipinski_violations"] == 0,
            veber_rule=descriptors["veber_rule"] == "Yes",
            egan_rule=descriptors["egan_rule"] == "Yes",
            drug_likeness=descriptors["drug_likeness"]
        )
    else:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=results.get("error", "RDKit calculation failed")
        )

@router.get("/lipinski/{job_id}", response_model=dict)
async def check_lipinski_rule(job_id: str):
    """Check Lipinski's Rule of Five"""
    
    # In real implementation, would retrieve from database
    descriptors = {}
    lipinski_violations = 0
    
    return JSONResponse(
        status_code=status.HTTP_200_OK,
        content={
            "job_id": job_id,
            "descriptors": descriptors,
            "violations": lipinski_violations,
            "compliant": lipinski_violations == 0
        }
    )
