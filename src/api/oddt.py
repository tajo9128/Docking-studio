"""
BioDockify Docking Studio - ODDT API Router
Handles ODDT interaction analysis operations
"""

from fastapi import APIRouter, HTTPException, status
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from typing import Dict, Any, List
import logging

from src.oddt_analyzer import ODDTAnalyzer
from src.docker_manager import DockerManager

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/oddt")

docker_manager = DockerManager()
oddt_analyzer = ODDTAnalyzer(docker_manager)

class InteractionType(BaseModel):
    """Interaction type model"""
    type: str
    count: int

class InteractionAnalysis(BaseModel):
    """Interaction analysis model"""
    total_count: int
    hydrogen_bonds: int
    hydrophobic_contacts: int
    pi_stacking: int
    halogen_bonds: int
    salt_bridges: int
    cation_pi: int
    metal_coordination: int

@router.post("/analyze", response_model=InteractionAnalysis)
async def analyze_interactions(receptor_file: str, ligand_file: str, output_file: str, job_id: str):
    """Analyze molecular interactions with ODDT"""
    
    results = oddt_analyzer.analyze_interactions(
        receptor_file=receptor_file,
        ligand_file=ligand_file,
        output_file=output_file,
        job_id=job_id
    )
    
    if results["status"] == "COMPLETED":
        interactions = results["interactions"]
        
        return InteractionAnalysis(
            total_count=interactions.get("total_count", 0),
            hydrogen_bonds=len(interactions.get("hydrogen_bonds", [])),
            hydrophobic_contacts=len(interactions.get("hydrophobic_contacts", [])),
            pi_stacking=len(interactions.get("pi_stacking", [])),
            halogen_bonds=len(interactions.get("halogen_bonds", [])),
            salt_bridges=len(interactions.get("salt_bridges", [])),
            cation_pi=len(interactions.get("cation_pi_interactions", [])),
            metal_coordination=len(interactions.get("metal_coordination", []))
        )
    else:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=results.get("error", "ODDT analysis failed")
        )

@router.get("/visualize/{job_id}", response_model=dict)
async def visualize_interactions(job_id: str):
    """Visualize interactions (2D diagram)"""
    
    # In real implementation, would generate 2D diagram
    diagram_data = {
        "job_id": job_id,
        "diagram": "Mock 2D interaction diagram"
    }
    
    return JSONResponse(
        status_code=status.HTTP_200_OK,
        content=diagram_data
    )
