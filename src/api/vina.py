"""
BioDockify Docking Studio - Vina API Router
Handles AutoDock Vina docking operations
"""

from fastapi import APIRouter, HTTPException, status
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from typing import Dict, Any, List
import logging

from src.vina_engine import VinaEngine
from src.docker_manager import DockerManager

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/vina")

docker_manager = DockerManager()
vina_engine = VinaEngine(docker_manager)

class VinaConfig(BaseModel):
    """Vina configuration model"""
    exhaustiveness: int = 8
    num_modes: int = 9
    box_size: float = 20.0

class VinaResult(BaseModel):
    """Vina result model"""
    binding_energy: float
    num_modes: int
    poses: List[Dict[str, Any]]

@router.post("/run", response_model=VinaResult)
async def run_vina_docking(receptor_file: str, ligand_file: str, config: VinaConfig):
    """Run AutoDock Vina docking"""
    
    parameters = {
        "center_x": 0.0,
        "center_y": 0.0,
        "center_z": 0.0,
        "size_x": config.box_size,
        "size_y": config.box_size,
        "size_z": config.box_size,
        "exhaustiveness": config.exhaustiveness,
        "num_modes": config.num_modes
    }
    
    import uuid
    job_id = str(uuid.uuid4())
    
    results = vina_engine.run_docking(
        receptor_file=receptor_file,
        ligand_file=ligand_file,
        parameters=parameters,
        job_id=job_id
    )
    
    if results["status"] == "COMPLETED":
        return VinaResult(
            binding_energy=results["binding_energy"],
            num_modes=results["num_modes"],
            poses=results["poses"]
        )
    else:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=results.get("logs", "Vina docking failed")
        )

@router.get("/output/{job_id}", response_model=dict)
async def get_vina_output(job_id: str):
    """Get Vina output for job"""
    
    # In real implementation, would read output file from database
    output_content = f"Mock Vina output for job {job_id}"
    
    return JSONResponse(
        status_code=status.HTTP_200_OK,
        content={
            "job_id": job_id,
            "output": output_content
        }
    )
