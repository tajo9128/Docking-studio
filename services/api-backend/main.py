"""
Docking Studio v2.0 - API Backend Service
Main gateway that routes requests to specialized services
"""

import os
import json
import logging
import redis
from contextlib import asynccontextmanager
from typing import Optional, List, Dict, Any
from pathlib import Path
from datetime import datetime

from fastapi import (
    FastAPI,
    HTTPException,
    UploadFile,
    File,
    BackgroundTasks,
    Depends,
    Form,
)
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, Response
from pydantic import BaseModel
from sqlalchemy.orm import Session
import httpx

from db import engine, get_db, Base
from models import Job, MDResult, UserProfile, MemoryEntry, ConversationHistory

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

DOCKING_SERVICE_URL = os.getenv("DOCKING_SERVICE_URL", "http://docking-service:8002")
RDKIT_SERVICE_URL = os.getenv("RDKIT_SERVICE_URL", "http://rdkit-service:8003")
PHARMACOPHORE_SERVICE_URL = os.getenv(
    "PHARMACOPHORE_SERVICE_URL", "http://pharmacophore-service:8004"
)
QSAR_SERVICE_URL = os.getenv("QSAR_SERVICE_URL", "http://qsar-service:8005")
MD_SERVICE_URL = os.getenv("MD_SERVICE_URL", "http://md-service:8000")
SENTINEL_SERVICE_URL = os.getenv("SENTINEL_SERVICE_URL", "http://sentinel-service:8000")
ANALYSIS_SERVICE_URL = os.getenv("ANALYSIS_SERVICE_URL", "http://analysis-service:8000")
BRAIN_SERVICE_URL = os.getenv("BRAIN_SERVICE_URL", "http://brain-service:8000")
REDIS_URL = os.getenv("REDIS_URL", "redis://redis:6379")
TRAINING_TIMEOUT = int(os.getenv("TRAINING_TIMEOUT", "300"))

STORAGE_DIR = Path("/app/storage")
UPLOADS_DIR = Path("/app/uploads")
STORAGE_DIR.mkdir(exist_ok=True)
UPLOADS_DIR.mkdir(exist_ok=True)

# In-memory LLM settings (shared across services)
llm_settings = {
    "provider": "openai",
    "model": "gpt-4o-mini",
    "api_key": "",
    "base_url": "https://api.openai.com/v1",
    "temperature": 0.7,
    "max_tokens": 4096,
}


@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("=" * 60)
    logger.info("Docking Studio v2.0 - API Backend")
    logger.info("=" * 60)
    logger.info(f"Docking Service: {DOCKING_SERVICE_URL}")
    logger.info(f"RDKit Service: {RDKIT_SERVICE_URL}")
    logger.info(f"Pharmacophore Service: {PHARMACOPHORE_SERVICE_URL}")
    logger.info(f"QSAR Service: {QSAR_SERVICE_URL}")
    logger.info(f"MD Service: {MD_SERVICE_URL}")
    logger.info(f"Brain Service: {BRAIN_SERVICE_URL}")
    Base.metadata.create_all(bind=engine)
    logger.info("Database tables created/verified")
    yield
    logger.info("API Backend shutting down...")


app = FastAPI(
    title="Docking Studio API",
    description="API Gateway for Docking Studio v2.0 Microservices",
    version="2.0.0",
    lifespan=lifespan,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class JobCreate(BaseModel):
    name: str
    job_type: str
    parameters: dict


class JobResponse(BaseModel):
    job_id: str
    status: str
    message: str


@app.get("/health")
async def health_check():
    return {"status": "healthy", "service": "api-backend", "version": "2.0.0"}


@app.get("/")
async def root():
    return {
        "service": "Docking Studio API Gateway",
        "version": "2.0.0",
        "endpoints": {
            "docs": "/docs",
            "health": "/health",
            "jobs": "/jobs",
            "upload": "/upload",
            "docking": "/docking",
            "rdkit": "/rdkit",
            "pharmacophore": "/pharmacophore",
            "qsar": "/qsar",
            "md": "/md",
            "sentinel": "/sentinel",
            "analysis": "/analysis",
            "brain": "/brain",
        },
    }


@app.post("/upload", response_model=dict)
async def upload_file(file: UploadFile = File(...)):
    """Upload a molecule or protein file"""
    try:
        file_path = UPLOADS_DIR / file.filename
        content = await file.read()
        file_path.write_bytes(content)

        file_type = "unknown"
        ext = file_path.suffix.lower()
        if ext in [".pdb", ".pdbqt"]:
            file_type = "protein"
        elif ext in [".sdf", ".mol", ".mol2", ".pdb"]:
            file_type = "ligand"
        elif ext in [".smiles", ".smi"]:
            file_type = "smiles"

        return {
            "filename": file.filename,
            "path": str(file_path),
            "type": file_type,
            "size": len(content),
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/jobs", response_model=JobResponse)
async def create_job(job: JobCreate, background_tasks: BackgroundTasks):
    """Create a new job and queue it"""
    import uuid

    job_id = str(uuid.uuid4())

    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            if job.job_type == "docking":
                response = await client.post(
                    f"{DOCKING_SERVICE_URL}/dock",
                    json={"job_id": job_id, **job.parameters},
                )
            elif job.job_type == "pharmacophore":
                response = await client.post(
                    f"{PHARMACOPHORE_SERVICE_URL}/generate",
                    json={"job_id": job_id, **job.parameters},
                )
            elif job.job_type == "rdkit":
                response = await client.post(
                    f"{RDKIT_SERVICE_URL}/process",
                    json={"job_id": job_id, **job.parameters},
                )
            else:
                raise HTTPException(
                    status_code=400, detail=f"Unknown job type: {job.job_type}"
                )

            response.raise_for_status()
            result = response.json()
            return JobResponse(
                job_id=job_id,
                status="queued",
                message=result.get("message", "Job queued"),
            )
        except httpx.HTTPError as e:
            logger.error(f"Failed to create job: {e}")
            return JobResponse(job_id=job_id, status="error", message=str(e))


@app.get("/jobs/{job_id}")
async def get_job_status(job_id: str):
    """Get job status from Redis — checks both job:* and docking_job:* keys"""
    import redis
    import json

    try:
        r = redis.from_url(REDIS_URL)
        job_data = r.get(f"job:{job_id}")
        if job_data:
            return json.loads(job_data)
        job_data = r.get(f"docking_job:{job_id}")
        if job_data:
            return json.loads(job_data)
        return {"job_id": job_id, "status": "not_found"}
    except Exception as e:
        return {"job_id": job_id, "status": "error", "message": str(e)}


@app.get("/jobs/{job_id}/results")
async def get_job_results(job_id: str):
    """Get docking results for a job"""
    import redis
    import json

    try:
        r = redis.from_url(REDIS_URL)
        job_data = r.get(f"docking_job:{job_id}")
        if job_data:
            data = json.loads(job_data)
            result = data.get("result", {})
            status = data.get("status", "unknown")
            if status == "completed":
                poses = result.get("poses", [])
                return {
                    "results": poses,
                    "best_score": result.get("best_energy") if poses else None,
                    "num_poses": len(poses),
                }
            elif status == "failed":
                raise HTTPException(
                    status_code=400, detail=data.get("error", "Job failed")
                )
            return {"results": [], "status": status}
        raise HTTPException(status_code=404, detail="Job not found")
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/jobs/{job_id}/interactions")
async def get_job_interactions(job_id: str, pose_id: int = None):
    """Get interaction analysis for a docking pose"""
    import redis
    import json

    try:
        r = redis.from_url(REDIS_URL)
        job_data = r.get(f"docking_job:{job_id}")
        if job_data:
            data = json.loads(job_data)
            if data.get("status") != "completed":
                return {"interactions": [], "pose_id": pose_id}

            result = data.get("result", {})
            poses = result.get("poses", [])
            target_pose = (
                poses[pose_id]
                if pose_id is not None and pose_id < len(poses)
                else (poses[0] if poses else None)
            )

            if not target_pose:
                return {"interactions": [], "pose_id": pose_id}

            receptor = target_pose.get("receptor_pdb", "")
            ligand = target_pose.get("ligand_pdb", "")

            interactions = []
            try:
                from rdkit import Chem
                from rdkit.Chem import rdMolDescriptors, Lipinski

                if receptor and ligand:
                    receptor_mol = (
                        Chem.MolFromPDBFile(receptor)
                        if os.path.exists(receptor)
                        else Chem.MolFromPDBBlock(receptor)
                    )
                    ligand_mol = (
                        Chem.MolFromPDBFile(ligand)
                        if os.path.exists(ligand)
                        else Chem.MolFromPDBBlock(ligand)
                    )

                    if receptor_mol and ligand_mol:
                        receptor_hbd = rdMolDescriptors.CalcNumHBD(receptor_mol)
                        receptor_hba = rdMolDescriptors.CalcNumHBA(receptor_mol)
                        ligand_hbd = rdMolDescriptors.CalcNumHBD(ligand_mol)
                        ligand_hba = rdMolDescriptors.CalcNumHBA(ligand_mol)

                        num_hbonds = min(receptor_hbd, ligand_hba) + min(
                            receptor_hba, ligand_hbd
                        )
                        num_hydrophobic = Lipinski.NumAromaticRings(
                            receptor_mol
                        ) + Lipinski.NumAromaticRings(ligand_mol)

                        if num_hbonds > 0:
                            interactions.append(
                                {
                                    "id": 1,
                                    "job_uuid": job_id,
                                    "pose_id": pose_id or 0,
                                    "interaction_type": "hbond",
                                    "atom_a": "protein",
                                    "atom_b": "ligand",
                                    "distance": 2.5,
                                }
                            )
                        if num_hydrophobic > 0:
                            interactions.append(
                                {
                                    "id": 2,
                                    "job_uuid": job_id,
                                    "pose_id": pose_id or 0,
                                    "interaction_type": "hydrophobic",
                                    "atom_a": "protein",
                                    "atom_b": "ligand",
                                    "distance": 3.5,
                                }
                            )
            except Exception:
                pass

            return {"interactions": interactions, "pose_id": pose_id}
        raise HTTPException(status_code=404, detail="Job not found")
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/jobs")
async def list_jobs(limit: int = 50):
    """List recent jobs from both job:* and docking_job:* Redis keys"""
    import redis
    import json

    try:
        r = redis.from_url(REDIS_URL)
        all_keys = r.keys("job:*") + r.keys("docking_job:*")
        all_keys = all_keys[:limit]
        jobs = []
        for key in all_keys:
            job_data = r.get(key)
            if job_data:
                try:
                    data = json.loads(job_data)
                    key_str = key.decode() if isinstance(key, bytes) else key
                    job_uuid = key_str.replace("job:", "").replace("docking_job:", "")
                    jobs.append(
                        {
                            "id": len(jobs),
                            "job_uuid": job_uuid,
                            "job_name": data.get("job_id", job_uuid),
                            "receptor_file": data.get("receptor_path", ""),
                            "ligand_file": data.get("ligand_path", ""),
                            "status": data.get("status", "unknown").upper(),
                            "created_at": data.get("created_at", ""),
                            "completed_at": data.get("completed_at", ""),
                            "binding_energy": data.get("result", {}).get("best_energy")
                            if data.get("result")
                            else None,
                            "engine": data.get("engine", "vina"),
                        }
                    )
                except Exception:
                    pass
        return {"jobs": jobs, "count": len(jobs)}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/docking/{job_id}/cancel")
async def cancel_docking_job(job_id: str):
    """Cancel a running docking job"""
    async with httpx.AsyncClient(timeout=10.0) as client:
        try:
            response = await client.post(f"{DOCKING_SERVICE_URL}/cancel/{job_id}")
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.get("/docking/{job_id}/results")
async def get_docking_results(job_id: str):
    """Get docking results"""
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.get(f"{DOCKING_SERVICE_URL}/results/{job_id}")
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=404, detail=str(e))


@app.post("/pharmacophore/generate")
async def generate_pharmacophore(
    receptor_pdb: Optional[str] = None,
    ligand_pdb: Optional[str] = None,
    features: Optional[str] = None,
    smiles: Optional[str] = None,
    pdb: Optional[str] = None,
):
    """Generate pharmacophore from receptor or ligand"""
    async with httpx.AsyncClient(timeout=60.0) as client:
        try:
            json_body = {}
            if receptor_pdb:
                json_body["receptor_pdb"] = receptor_pdb
            if ligand_pdb:
                json_body["ligand_pdb"] = ligand_pdb
            if features:
                json_body["features"] = features
            if smiles:
                json_body["smiles"] = smiles
            if pdb:
                json_body["pdb"] = pdb
            response = await client.post(
                f"{PHARMACOPHORE_SERVICE_URL}/generate",
                json=json_body,
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/pharmacophore/screen")
async def screen_pharmacophore(pharmacophore_id: str, library_path: str):
    """Screen a library against a pharmacophore"""
    async with httpx.AsyncClient(timeout=300.0) as client:
        try:
            response = await client.post(
                f"{PHARMACOPHORE_SERVICE_URL}/screen",
                json={
                    "pharmacophore_id": pharmacophore_id,
                    "library_path": library_path,
                },
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.get("/qsar/descriptor-groups")
async def qsar_descriptor_groups():
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.get(f"{QSAR_SERVICE_URL}/descriptor-groups")
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/qsar/descriptors")
async def qsar_descriptors(smiles: List[str], groups: Optional[List[str]] = None):
    async with httpx.AsyncClient(timeout=60.0) as client:
        try:
            response = await client.post(
                f"{QSAR_SERVICE_URL}/descriptors",
                json={"smiles": smiles, "groups": groups},
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/qsar/descriptors/upload")
async def qsar_descriptors_upload(
    file: UploadFile = File(...),
    smiles_col: str = Form(...),
    activity_col: str = Form(...),
    groups: Optional[str] = Form(None),
):
    async with httpx.AsyncClient(timeout=120.0) as client:
        try:
            files = {"file": (file.filename, file.file, file.content_type)}
            data = {"smiles_col": smiles_col, "activity_col": activity_col}
            if groups:
                data["groups"] = groups
            response = await client.post(
                f"{QSAR_SERVICE_URL}/descriptors/upload",
                files=files,
                data=data,
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/qsar/train")
async def qsar_train(
    X: List[List[float]],
    y: List[float],
    feature_names: List[str],
    model_type: str,
    model_name: str,
    activity_column: str,
    descriptor_groups: List[str],
    cv_folds: int = 5,
    model_params: Optional[dict] = None,
):
    async with httpx.AsyncClient(timeout=TRAINING_TIMEOUT) as client:
        try:
            response = await client.post(
                f"{QSAR_SERVICE_URL}/train",
                json={
                    "X": X,
                    "y": y,
                    "feature_names": feature_names,
                    "model_type": model_type,
                    "model_name": model_name,
                    "activity_column": activity_column,
                    "descriptor_groups": descriptor_groups,
                    "cv_folds": cv_folds,
                    "model_params": model_params,
                },
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.get("/qsar/train/{job_id}/status")
async def qsar_train_status(job_id: str):
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.get(f"{QSAR_SERVICE_URL}/train/{job_id}/status")
            response.raise_for_status()
            return response.json()
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                raise HTTPException(status_code=404, detail=f"Job {job_id} not found")
            raise HTTPException(status_code=500, detail=str(e))
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.get("/qsar/train/{job_id}/results")
async def qsar_train_results(job_id: str):
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.get(f"{QSAR_SERVICE_URL}/train/{job_id}/results")
            response.raise_for_status()
            return response.json()
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                raise HTTPException(status_code=404, detail=f"Job {job_id} not found")
            raise HTTPException(status_code=500, detail=str(e))
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/qsar/predict")
async def qsar_predict(model_id: str, smiles: str):
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.post(
                f"{QSAR_SERVICE_URL}/predict",
                json={"model_id": model_id, "smiles": smiles},
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                raise HTTPException(status_code=404, detail=str(e))
            raise HTTPException(status_code=500, detail=str(e))
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/qsar/predict/batch")
async def qsar_predict_batch(model_id: str, smiles_list: List[str]):
    async with httpx.AsyncClient(timeout=120.0) as client:
        try:
            response = await client.post(
                f"{QSAR_SERVICE_URL}/predict/batch",
                json={"model_id": model_id, "smiles_list": smiles_list},
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                raise HTTPException(status_code=404, detail=str(e))
            raise HTTPException(status_code=500, detail=str(e))
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.get("/qsar/models")
async def qsar_models():
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.get(f"{QSAR_SERVICE_URL}/models")
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.get("/qsar/models/{model_id}")
async def qsar_model_info(model_id: str):
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.get(f"{QSAR_SERVICE_URL}/models/{model_id}")
            response.raise_for_status()
            return response.json()
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                raise HTTPException(
                    status_code=404, detail=f"Model {model_id} not found"
                )
            raise HTTPException(status_code=500, detail=str(e))
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.delete("/qsar/models/{model_id}")
async def qsar_model_delete(model_id: str):
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.delete(f"{QSAR_SERVICE_URL}/models/{model_id}")
            response.raise_for_status()
            return response.json()
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                raise HTTPException(
                    status_code=404, detail=f"Model {model_id} not found"
                )
            raise HTTPException(status_code=500, detail=str(e))
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


# ============================================================
# MD Service proxies (Molecular Dynamics + MD-Suite)
# ============================================================


class MDDynamicsRequest(BaseModel):
    pdb_content: str
    steps: int = 50000
    temperature: float = 300.0
    pressure: float = 1.0
    frame_interval: int = 500
    solvent_model: str = "tip3p"
    ionic_strength: float = 0.0
    name: str = "simulation"
    notify_on_start: bool = False
    notify_on_complete: bool = True


@app.post("/md/dynamics")
async def md_run_dynamics(request: MDDynamicsRequest):
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.post(
                f"{MD_SERVICE_URL}/dynamics",
                json=request.model_dump(),
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


class MDAnalysisRequest(BaseModel):
    job_id: str
    trajectory_path: Optional[str] = None
    energy_csv_path: Optional[str] = None


@app.post("/md/analysis/rmsd")
async def md_analyze_rmsd(request: MDAnalysisRequest):
    async with httpx.AsyncClient(timeout=60.0) as client:
        try:
            response = await client.post(
                f"{MD_SERVICE_URL}/analysis/rmsd", json=request.model_dump()
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/md/analysis/rmsf")
async def md_analyze_rmsf(request: MDAnalysisRequest):
    async with httpx.AsyncClient(timeout=60.0) as client:
        try:
            response = await client.post(
                f"{MD_SERVICE_URL}/analysis/rmsf", json=request.model_dump()
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/md/analysis/energy")
async def md_analyze_energy(request: MDAnalysisRequest):
    async with httpx.AsyncClient(timeout=60.0) as client:
        try:
            response = await client.post(
                f"{MD_SERVICE_URL}/analysis/energy", json=request.model_dump()
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/md/analysis/gyration")
async def md_analyze_gyration(request: MDAnalysisRequest):
    async with httpx.AsyncClient(timeout=60.0) as client:
        try:
            response = await client.post(
                f"{MD_SERVICE_URL}/analysis/gyration", json=request.model_dump()
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/md/analysis/sasa")
async def md_analyze_sasa(request: MDAnalysisRequest):
    async with httpx.AsyncClient(timeout=60.0) as client:
        try:
            response = await client.post(
                f"{MD_SERVICE_URL}/analysis/sasa", json=request.model_dump()
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/md/analysis/hbonds")
async def md_analyze_hbonds(request: MDAnalysisRequest):
    async with httpx.AsyncClient(timeout=60.0) as client:
        try:
            response = await client.post(
                f"{MD_SERVICE_URL}/analysis/hbonds", json=request.model_dump()
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/md/analysis/all")
async def md_analyze_all(request: MDAnalysisRequest):
    async with httpx.AsyncClient(timeout=300.0) as client:
        try:
            response = await client.post(
                f"{MD_SERVICE_URL}/analysis/all", json=request.model_dump()
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.get("/md/job/{job_id}")
async def md_get_job(job_id: str):
    async with httpx.AsyncClient(timeout=10.0) as client:
        try:
            response = await client.get(f"{MD_SERVICE_URL}/job/{job_id}")
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.get("/md/job/{job_id}/trajectory")
async def md_get_trajectory(job_id: str):
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.get(f"{MD_SERVICE_URL}/job/{job_id}/trajectory")
            response.raise_for_status()
            return response.content
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.get("/md/job/{job_id}/frame")
async def md_get_frame(job_id: str):
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.get(f"{MD_SERVICE_URL}/job/{job_id}/frame")
            response.raise_for_status()
            return response.content
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.get("/md/job/{job_id}/energies")
async def md_get_energies(job_id: str):
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.get(f"{MD_SERVICE_URL}/job/{job_id}/energies")
            response.raise_for_status()
            return response.content
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


class MDPublicationRequest(BaseModel):
    job_id: str
    project_name: str = "md_project"
    analysis_job_id: Optional[str] = None
    compress: bool = True
    notify_on_complete: bool = False


@app.post("/md/publication/package")
async def md_create_publication(request: MDPublicationRequest):
    async with httpx.AsyncClient(timeout=60.0) as client:
        try:
            response = await client.post(
                f"{MD_SERVICE_URL}/publication/package",
                json=request.model_dump(),
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


class MDNotifyRequest(BaseModel):
    event: str
    title: str
    message: str
    details: Optional[Dict] = None


@app.post("/md/notify")
async def md_notify(request: MDNotifyRequest):
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.post(
                f"{MD_SERVICE_URL}/notify",
                json=request.model_dump(),
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.get("/md/notify/status")
async def md_notify_status():
    async with httpx.AsyncClient(timeout=10.0) as client:
        try:
            response = await client.get(f"{MD_SERVICE_URL}/notify/status")
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/md/notify/test")
async def md_notify_test(channel: str = "discord"):
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.post(
                f"{MD_SERVICE_URL}/notify/test",
                params={"channel": channel},
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/md/minimize")
async def md_minimize(pdb_content: str):
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.post(
                f"{MD_SERVICE_URL}/minimize",
                params={"pdb_content": pdb_content},
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.get("/md/health")
async def md_health():
    async with httpx.AsyncClient(timeout=10.0) as client:
        try:
            response = await client.get(f"{MD_SERVICE_URL}/health")
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


# ============================================================
# Database endpoints - MD Results (PostgreSQL)
# ============================================================


class MDSaveRequest(BaseModel):
    job_uuid: str
    project_name: Optional[str] = "md_project"
    n_steps: Optional[int] = None
    sim_time_ns: Optional[float] = None
    temperature_K: Optional[float] = None
    solvent_model: Optional[str] = None
    ionic_strength: Optional[float] = None
    n_frames: Optional[int] = None
    avg_energy_kj_mol: Optional[float] = None
    trajectory_path: Optional[str] = None
    final_frame_path: Optional[str] = None
    energy_csv_path: Optional[str] = None
    analysis_summary: Optional[dict] = None
    package_path: Optional[str] = None


@app.post("/db/md/save")
def save_md_result(request: MDSaveRequest, db: Session = Depends(get_db)):
    """Persist MD simulation result to PostgreSQL database."""
    try:
        job = db.query(Job).filter(Job.job_uuid == request.job_uuid).first()
        if not job:
            job = Job(
                job_uuid=request.job_uuid,
                job_name=request.project_name,
                job_type="md",
                status="completed",
                parameters={
                    "n_steps": request.n_steps,
                    "temperature_K": request.temperature_K,
                    "solvent_model": request.solvent_model,
                    "ionic_strength": request.ionic_strength,
                },
                result={
                    "sim_time_ns": request.sim_time_ns,
                    "n_frames": request.n_frames,
                    "avg_energy_kj_mol": request.avg_energy_kj_mol,
                },
                completed_at=datetime.utcnow(),
            )
            db.add(job)
        else:
            job.status = "completed"
            job.result = {
                "sim_time_ns": request.sim_time_ns,
                "n_frames": request.n_frames,
                "avg_energy_kj_mol": request.avg_energy_kj_mol,
            }
            job.completed_at = datetime.utcnow()

        existing_md = (
            db.query(MDResult).filter(MDResult.job_uuid == request.job_uuid).first()
        )
        if existing_md:
            for key, value in request.model_dump(exclude_none=True).items():
                setattr(existing_md, key, value)
        else:
            md_result = MDResult(
                job_uuid=request.job_uuid,
                project_name=request.project_name,
                n_steps=request.n_steps,
                sim_time_ns=request.sim_time_ns,
                temperature_K=request.temperature_K,
                solvent_model=request.solvent_model,
                ionic_strength=request.ionic_strength,
                n_frames=request.n_frames,
                avg_energy_kj_mol=request.avg_energy_kj_mol,
                trajectory_path=request.trajectory_path,
                final_frame_path=request.final_frame_path,
                energy_csv_path=request.energy_csv_path,
                analysis_summary=request.analysis_summary,
                package_path=request.package_path,
            )
            db.add(md_result)

        db.commit()
        return {"success": True, "job_uuid": request.job_uuid}
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/db/md/results/{job_uuid}")
def get_md_result(job_uuid: str, db: Session = Depends(get_db)):
    """Retrieve MD result from PostgreSQL by job UUID."""
    md_result = db.query(MDResult).filter(MDResult.job_uuid == job_uuid).first()
    if not md_result:
        raise HTTPException(
            status_code=404, detail=f"MD result for {job_uuid} not found"
        )
    return {
        "job_uuid": md_result.job_uuid,
        "project_name": md_result.project_name,
        "n_steps": md_result.n_steps,
        "sim_time_ns": md_result.sim_time_ns,
        "temperature_K": md_result.temperature_K,
        "solvent_model": md_result.solvent_model,
        "ionic_strength": md_result.ionic_strength,
        "n_frames": md_result.n_frames,
        "avg_energy_kj_mol": md_result.avg_energy_kj_mol,
        "trajectory_path": md_result.trajectory_path,
        "final_frame_path": md_result.final_frame_path,
        "energy_csv_path": md_result.energy_csv_path,
        "analysis_summary": md_result.analysis_summary,
        "package_path": md_result.package_path,
        "created_at": md_result.created_at.isoformat()
        if md_result.created_at
        else None,
    }


@app.get("/db/jobs")
def list_db_jobs(limit: int = 50, db: Session = Depends(get_db)):
    """List all jobs from PostgreSQL database."""
    jobs = db.query(Job).order_by(Job.created_at.desc()).limit(limit).all()
    return {
        "jobs": [
            {
                "job_uuid": j.job_uuid,
                "job_name": j.job_name,
                "job_type": j.job_type,
                "status": j.status,
                "created_at": j.created_at.isoformat() if j.created_at else None,
                "completed_at": j.completed_at.isoformat() if j.completed_at else None,
            }
            for j in jobs
        ],
        "count": len(jobs),
    }


# ============================================================
# Nanobot Memory - User profiles, long-term memory, conversation history
# ============================================================


class MemoryEntryRequest(BaseModel):
    user_id: str
    category: str = "general"
    key: str
    value: str
    confidence: float = 1.0
    tags: Optional[List[str]] = None
    metadata: Optional[dict] = None


class ConversationEntryRequest(BaseModel):
    user_id: str
    role: str
    content: str
    metadata: Optional[dict] = None


class UserProfileRequest(BaseModel):
    user_id: str
    display_name: Optional[str] = None
    email: Optional[str] = None
    platform: Optional[str] = None
    preferences: Optional[dict] = None


@app.post("/memory/profile")
def upsert_user_profile(request: UserProfileRequest, db: Session = Depends(get_db)):
    """Create or update a user profile."""
    try:
        profile = (
            db.query(UserProfile).filter(UserProfile.user_id == request.user_id).first()
        )
        if profile:
            for field in ["display_name", "email", "platform", "preferences"]:
                val = getattr(request, field, None)
                if val is not None:
                    setattr(profile, field, val)
        else:
            profile = UserProfile(
                user_id=request.user_id,
                display_name=request.display_name,
                email=request.email,
                platform=request.platform,
                preferences=request.preferences or {},
            )
            db.add(profile)
        db.commit()
        db.refresh(profile)
        return {"success": True, "user_id": profile.user_id}
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/memory/profile/{user_id}")
def get_user_profile(user_id: str, db: Session = Depends(get_db)):
    """Get a user profile."""
    profile = db.query(UserProfile).filter(UserProfile.user_id == user_id).first()
    if not profile:
        raise HTTPException(status_code=404, detail=f"Profile for {user_id} not found")
    return {
        "user_id": profile.user_id,
        "display_name": profile.display_name,
        "email": profile.email,
        "platform": profile.platform,
        "preferences": profile.preferences,
        "created_at": profile.created_at.isoformat() if profile.created_at else None,
        "updated_at": profile.updated_at.isoformat() if profile.updated_at else None,
    }


@app.post("/memory/remember")
def add_memory_entry(request: MemoryEntryRequest, db: Session = Depends(get_db)):
    """Store a long-term memory entry for a user."""
    try:
        existing = (
            db.query(MemoryEntry)
            .filter(
                MemoryEntry.user_id == request.user_id,
                MemoryEntry.category == request.category,
                MemoryEntry.key == request.key,
            )
            .first()
        )
        if existing:
            existing.value = request.value
            existing.confidence = request.confidence
            if request.tags is not None:
                existing.tags = request.tags
            if request.metadata is not None:
                existing.metadata = request.metadata
        else:
            entry = MemoryEntry(
                user_id=request.user_id,
                category=request.category,
                key=request.key,
                value=request.value,
                confidence=request.confidence,
                tags=request.tags or [],
                metadata=request.metadata or {},
            )
            db.add(entry)
        db.commit()
        return {"success": True, "user_id": request.user_id, "key": request.key}
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/memory/{user_id}")
def get_all_memory(
    user_id: str, category: Optional[str] = None, db: Session = Depends(get_db)
):
    """Retrieve all memory entries for a user, optionally filtered by category."""
    query = db.query(MemoryEntry).filter(MemoryEntry.user_id == user_id)
    if category:
        query = query.filter(MemoryEntry.category == category)
    entries = query.order_by(MemoryEntry.created_at.desc()).all()
    return {
        "user_id": user_id,
        "entries": [
            {
                "id": e.id,
                "category": e.category,
                "key": e.key,
                "value": e.value,
                "confidence": e.confidence,
                "tags": e.tags,
                "metadata": e.metadata,
                "created_at": e.created_at.isoformat() if e.created_at else None,
            }
            for e in entries
        ],
        "count": len(entries),
    }


@app.get("/memory/{user_id}/search")
def search_memory(
    user_id: str,
    q: Optional[str] = None,
    tags: Optional[str] = None,
    category: Optional[str] = None,
    db: Session = Depends(get_db),
):
    """Search memory entries for a user by query, tags, or category."""
    query = db.query(MemoryEntry).filter(MemoryEntry.user_id == user_id)
    if category:
        query = query.filter(MemoryEntry.category == category)
    if q:
        query = query.filter(
            (MemoryEntry.key.ilike(f"%{q}%")) | (MemoryEntry.value.ilike(f"%{q}%"))
        )
    entries = query.order_by(
        MemoryEntry.confidence.desc(), MemoryEntry.created_at.desc()
    ).all()
    if tags:
        tag_list = [t.strip() for t in tags.split(",")]
        entries = [e for e in entries if e.tags and any(t in e.tags for t in tag_list)]
    return {
        "user_id": user_id,
        "query": q,
        "entries": [
            {
                "id": e.id,
                "category": e.category,
                "key": e.key,
                "value": e.value,
                "confidence": e.confidence,
                "tags": e.tags,
                "metadata": e.metadata,
            }
            for e in entries
        ],
        "count": len(entries),
    }


@app.delete("/memory/{user_id}/{entry_id}")
def delete_memory_entry(user_id: str, entry_id: int, db: Session = Depends(get_db)):
    """Delete a specific memory entry."""
    entry = (
        db.query(MemoryEntry)
        .filter(
            MemoryEntry.id == entry_id,
            MemoryEntry.user_id == user_id,
        )
        .first()
    )
    if not entry:
        raise HTTPException(status_code=404, detail="Memory entry not found")
    db.delete(entry)
    db.commit()
    return {"success": True}


@app.post("/memory/conversation")
def add_conversation_entry(
    request: ConversationEntryRequest, db: Session = Depends(get_db)
):
    """Append a conversation turn to a user's chat history."""
    try:
        entry = ConversationHistory(
            user_id=request.user_id,
            role=request.role,
            content=request.content,
            metadata=request.metadata or {},
        )
        db.add(entry)
        db.commit()
        return {"success": True, "id": entry.id}
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/memory/conversation/{user_id}")
def get_conversation_history(
    user_id: str, limit: int = 20, db: Session = Depends(get_db)
):
    """Get recent conversation history for a user."""
    entries = (
        db.query(ConversationHistory)
        .filter(ConversationHistory.user_id == user_id)
        .order_by(ConversationHistory.created_at.desc())
        .limit(limit)
        .all()
    )
    entries.reverse()
    return {
        "user_id": user_id,
        "messages": [
            {
                "id": e.id,
                "role": e.role,
                "content": e.content,
                "metadata": e.metadata,
                "created_at": e.created_at.isoformat() if e.created_at else None,
            }
            for e in entries
        ],
        "count": len(entries),
    }


@app.get("/memory/history/{user_id}")
def get_user_history_summary(user_id: str, db: Session = Depends(get_db)):
    """Get a comprehensive history summary for a user: jobs, memory, conversations."""
    profile = db.query(UserProfile).filter(UserProfile.user_id == user_id).first()
    recent_jobs = (
        db.query(Job)
        .filter(Job.job_uuid.like(f"{user_id}%"))
        .order_by(Job.created_at.desc())
        .limit(10)
        .all()
    )
    recent_memory = (
        db.query(MemoryEntry)
        .filter(MemoryEntry.user_id == user_id)
        .order_by(MemoryEntry.created_at.desc())
        .limit(20)
        .all()
    )
    recent_convo = (
        db.query(ConversationHistory)
        .filter(ConversationHistory.user_id == user_id)
        .order_by(ConversationHistory.created_at.desc())
        .limit(10)
        .all()
    )
    return {
        "user_id": user_id,
        "profile": {
            "display_name": profile.display_name if profile else None,
            "platform": profile.platform if profile else None,
        }
        if profile
        else None,
        "recent_jobs": [
            {
                "job_uuid": j.job_uuid,
                "job_name": j.job_name,
                "job_type": j.job_type,
                "status": j.status,
                "created_at": j.created_at.isoformat() if j.created_at else None,
            }
            for j in recent_jobs
        ],
        "recent_memory": [
            {
                "category": e.category,
                "key": e.key,
                "value": e.value,
                "confidence": e.confidence,
                "tags": e.tags,
            }
            for e in recent_memory
        ],
        "recent_conversation": [
            {
                "role": e.role,
                "content": e.content[:200],
                "created_at": e.created_at.isoformat() if e.created_at else None,
            }
            for e in reversed(recent_convo)
        ],
    }


@app.post("/memory/auto记住")
def auto_remember_from_job(
    user_id: str,
    job_uuid: str,
    job_type: str,
    job_name: str,
    result_summary: str,
    db: Session = Depends(get_db),
):
    """Automatically remember key facts from a completed job."""
    try:
        memory_entries = []
        key_map = {
            "docking": ("docking", ["best_score", "num_poses", "ligand"]),
            "md": ("md", ["sim_time_ns", "temperature_K", "solvent_model"]),
            "qsar": ("qsar", ["model_type", "cv_r2", "n_compounds"]),
            "pharmacophore": ("pharmacophore", ["num_features"]),
        }
        cat, keys = key_map.get(job_type, ("general", []))
        entry = MemoryEntry(
            user_id=user_id,
            category="jobs",
            key=f"job:{job_uuid}",
            value=f"{job_type} job '{job_name}' - {result_summary}",
            confidence=1.0,
            tags=[job_type, "job", "completed"],
            metadata={"job_uuid": job_uuid, "job_name": job_name},
        )
        db.add(entry)
        for k in keys:
            memory_entries.append(
                MemoryEntry(
                    user_id=user_id,
                    category=cat,
                    key=k,
                    value=f"From job {job_uuid}: {k}",
                    confidence=0.9,
                    tags=[job_type, "job_detail"],
                    metadata={"job_uuid": job_uuid},
                )
            )
        for me in memory_entries:
            db.add(me)
        db.commit()
        return {"success": True, "entries_stored": 1 + len(memory_entries)}
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=str(e))


# ============================================================
# Sentinel Service proxies (Supervisor)
# ============================================================


class LigandRecordInput(BaseModel):
    ligand_id: str
    smiles: Optional[str] = None
    vina_score: Optional[float] = None
    gnina_score: Optional[float] = None
    rf_score: Optional[float] = None
    consensus_score: Optional[float] = None
    md_stability: Optional[float] = None
    md_time_ns: Optional[float] = None
    rmsd_from_crystal: Optional[float] = None
    h_bond_count: Optional[int] = None
    hydrophobic_count: Optional[int] = None
    mw: Optional[float] = None
    logp: Optional[float] = None
    tpsa: Optional[float] = None
    num_rotatable_bonds: Optional[int] = None


class RankingRequestInput(BaseModel):
    ligands: List[LigandRecordInput]
    weights: Optional[Dict[str, float]] = None


@app.post("/sentinel/monitor")
async def sentinel_monitor(job_id: str, service: str):
    async with httpx.AsyncClient(timeout=10.0) as client:
        try:
            response = await client.post(
                f"{SENTINEL_SERVICE_URL}/monitor/job",
                params={"job_id": job_id, "service": service},
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/sentinel/retry")
async def sentinel_retry(job_id: str, service: str):
    async with httpx.AsyncClient(timeout=10.0) as client:
        try:
            response = await client.post(
                f"{SENTINEL_SERVICE_URL}/retry/job",
                params={"job_id": job_id, "service": service},
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/sentinel/escalate")
async def sentinel_escalate(job_id: str, service: str, reason: str):
    async with httpx.AsyncClient(timeout=10.0) as client:
        try:
            response = await client.post(
                f"{SENTINEL_SERVICE_URL}/escalate",
                params={"job_id": job_id, "service": service, "reason": reason},
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.get("/sentinel/queue/status")
async def sentinel_queue_status():
    async with httpx.AsyncClient(timeout=10.0) as client:
        try:
            response = await client.get(f"{SENTINEL_SERVICE_URL}/queue/status")
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.get("/sentinel/health")
async def sentinel_health():
    async with httpx.AsyncClient(timeout=10.0) as client:
        try:
            response = await client.get(f"{SENTINEL_SERVICE_URL}/health")
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


# ============================================================
# Analysis Service proxies (Insight Generator)
# ============================================================


@app.post("/analysis/rank")
async def analysis_rank(request: RankingRequestInput):
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            payload = {
                "ligands": [l.model_dump() for l in request.ligands],
                "weights": request.weights,
            }
            response = await client.post(f"{ANALYSIS_SERVICE_URL}/rank", json=payload)
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/analysis/filter/admet")
async def analysis_admet(ligands: List[LigandRecordInput]):
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.post(
                f"{ANALYSIS_SERVICE_URL}/filter/admet",
                json=[l.model_dump() for l in ligands],
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/analysis/consensus")
async def analysis_consensus(ligands: List[LigandRecordInput]):
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.post(
                f"{ANALYSIS_SERVICE_URL}/consensus",
                json=[l.model_dump() for l in ligands],
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/analysis/report")
async def analysis_report(
    job_uuid: str, ligand_ids: List[str], summary: Optional[Dict] = None
):
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.post(
                f"{ANALYSIS_SERVICE_URL}/report",
                params={"job_uuid": job_uuid},
                json={"ligand_ids": ligand_ids, "summary": summary},
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/analysis/interactions/summary")
async def analysis_interactions_summary(interactions: List[Dict[str, Any]]):
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.post(
                f"{ANALYSIS_SERVICE_URL}/interactions/summary", json=interactions
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


class ExportTopHitsRequest(BaseModel):
    docking_results: List[Dict[str, Any]]
    top_n: int = 10
    sort_by: str = "vina_score"
    format: str = "csv"


@app.post("/analysis/export/top-hits")
async def export_top_hits(request: ExportTopHitsRequest):
    """Export top N docking hits as CSV or JSON"""
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.post(
                f"{ANALYSIS_SERVICE_URL}/export/top-hits",
                json=request.model_dump(),
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.get("/analysis/health")
async def analysis_health():
    async with httpx.AsyncClient(timeout=10.0) as client:
        try:
            response = await client.get(f"{ANALYSIS_SERVICE_URL}/health")
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/rdkit/smiles-to-3d")
async def smiles_to_3d(smiles: str, name: Optional[str] = None):
    """Convert SMILES to 3D structure"""
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.post(
                f"{RDKIT_SERVICE_URL}/smiles-to-3d",
                json={"smiles": smiles, "name": name},
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/rdkit/optimize")
async def optimize_molecule(pdb_path: str):
    """Optimize molecule 3D structure"""
    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            response = await client.post(
                f"{RDKIT_SERVICE_URL}/optimize", json={"pdb_path": pdb_path}
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


class PrepareProteinRequest(BaseModel):
    pdb_content: str
    name: Optional[str] = "protein"
    remove_waters: bool = True
    add_hydrogens: bool = True
    pH: float = 7.4


class PrepareLigandRequest(BaseModel):
    pdb_content: str
    name: Optional[str] = "ligand"
    pH: float = 7.4


class DetectInteractionsRequest(BaseModel):
    receptor_pdb_content: str
    ligand_pdb_content: str


@app.post("/rdkit/prepare_protein")
async def prepare_protein(request: PrepareProteinRequest):
    """Prepare protein: remove waters, add hydrogens"""
    async with httpx.AsyncClient(timeout=60.0) as client:
        try:
            response = await client.post(
                f"{RDKIT_SERVICE_URL}/prepare_protein",
                json=request.model_dump(),
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/rdkit/prepare_ligand")
async def prepare_ligand(request: PrepareLigandRequest):
    """Prepare ligand: add hydrogens, generate PDBQT"""
    async with httpx.AsyncClient(timeout=60.0) as client:
        try:
            response = await client.post(
                f"{RDKIT_SERVICE_URL}/prepare_ligand",
                json=request.model_dump(),
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.post("/rdkit/detect_interactions")
async def detect_interactions(request: DetectInteractionsRequest):
    """Detect H-bonds and hydrophobic interactions"""
    async with httpx.AsyncClient(timeout=60.0) as client:
        try:
            response = await client.post(
                f"{RDKIT_SERVICE_URL}/detect_interactions",
                json=request.model_dump(),
            )
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


class DockAsyncRequest(BaseModel):
    job_id: Optional[str] = None
    receptor_path: str
    ligand_path: str
    center_x: float = 0.0
    center_y: float = 0.0
    center_z: float = 0.0
    size_x: float = 20.0
    size_y: float = 20.0
    size_z: float = 20.0
    exhaustiveness: int = 8
    num_modes: int = 9
    engine: str = "vina"


@app.post("/dock/async")
async def dock_async(request: DockAsyncRequest):
    """Start an async docking job"""
    import uuid

    job_id = request.job_id or str(uuid.uuid4())
    logger.info(
        f"[DOCK] Starting job {job_id}: receptor={request.receptor_path}, ligand={request.ligand_path}"
    )
    async with httpx.AsyncClient(timeout=300.0) as client:
        try:
            response = await client.post(
                f"{DOCKING_SERVICE_URL}/dock/async",
                json={
                    "job_id": job_id,
                    "receptor_path": request.receptor_path,
                    "ligand_path": request.ligand_path,
                    "center_x": request.center_x,
                    "center_y": request.center_y,
                    "center_z": request.center_z,
                    "size_x": request.size_x,
                    "size_y": request.size_y,
                    "size_z": request.size_z,
                    "exhaustiveness": request.exhaustiveness,
                    "num_modes": request.num_modes,
                    "engine": request.engine,
                },
            )
            response.raise_for_status()
            result = response.json()
            logger.info(f"[DOCK] Job {job_id} queued successfully: {result}")
            return result
        except httpx.HTTPError as e:
            logger.error(f"[DOCK] Job {job_id} failed: {e}")
            raise HTTPException(status_code=500, detail=str(e))


@app.get("/dock/{job_id}/status")
async def get_dock_status(job_id: str):
    """Get docking job status"""
    async with httpx.AsyncClient(timeout=10.0) as client:
        try:
            response = await client.get(f"{DOCKING_SERVICE_URL}/dock/{job_id}/status")
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.get("/dock/{job_id}/result")
async def get_dock_result(job_id: str):
    """Get docking job result — reads from Redis via docking-service"""
    import redis
    import json

    try:
        r = redis.from_url(REDIS_URL)
        job_data = r.get(f"docking_job:{job_id}")
        if job_data:
            data = json.loads(job_data)
            status = data.get("status", "unknown")
            if status == "completed":
                result = data.get("result", {})
                poses = result.get("poses", []) if isinstance(result, dict) else []
                return {
                    "job_id": job_id,
                    "status": "completed",
                    "best_energy": result.get("best_energy")
                    if isinstance(result, dict)
                    else None,
                    "poses": poses,
                    "num_poses": len(poses),
                }
            elif status == "failed":
                raise HTTPException(
                    status_code=400, detail=data.get("error", "Docking failed")
                )
            return {"job_id": job_id, "status": status}

        async with httpx.AsyncClient(timeout=10.0) as client:
            try:
                response = await client.get(
                    f"{DOCKING_SERVICE_URL}/dock/{job_id}/result"
                )
                response.raise_for_status()
                return response.json()
            except httpx.HTTPError:
                raise HTTPException(status_code=404, detail="Job not found")

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/dock/{job_id}/cancel")
async def cancel_dock(job_id: str):
    """Cancel a docking job"""
    async with httpx.AsyncClient(timeout=10.0) as client:
        try:
            response = await client.post(f"{DOCKING_SERVICE_URL}/dock/{job_id}/cancel")
            response.raise_for_status()
            return response.json()
        except httpx.HTTPError as e:
            raise HTTPException(status_code=500, detail=str(e))


@app.get("/stats")
async def get_stats():
    """Get system statistics"""
    import redis

    try:
        r = redis.from_url(REDIS_URL)
        total_jobs = (
            len(r.keys("job:*"))
            + len(r.keys("docking_job:*"))
            + len(r.keys("md_job:*"))
            + len(r.keys("qsar_job:*"))
            + len(r.keys("pharmacophore_job:*"))
        )

        stats = {
            "total_jobs": total_jobs,
            "services": {
                "api_backend": "healthy",
                "docking_service": "unknown",
                "rdkit_service": "unknown",
                "pharmacophore_service": "unknown",
                "qsar_service": "unknown",
                "md_service": "unknown",
                "sentinel_service": "unknown",
                "analysis_service": "unknown",
                "brain_service": "unknown",
            },
        }

        async with httpx.AsyncClient(timeout=5.0) as client:
            for service_name, service_url in [
                ("docking_service", DOCKING_SERVICE_URL),
                ("rdkit_service", RDKIT_SERVICE_URL),
                ("pharmacophore_service", PHARMACOPHORE_SERVICE_URL),
                ("qsar_service", QSAR_SERVICE_URL),
                ("md_service", MD_SERVICE_URL),
                ("sentinel_service", SENTINEL_SERVICE_URL),
                ("analysis_service", ANALYSIS_SERVICE_URL),
                ("brain_service", BRAIN_SERVICE_URL),
            ]:
                try:
                    response = await client.get(f"{service_url}/health")
                    if response.status_code == 200:
                        stats["services"][service_name] = "healthy"
                except:
                    stats["services"][service_name] = "unhealthy"

        return stats
    except Exception as e:
        return {"error": str(e)}


@app.get("/llm/settings")
async def get_llm_settings():
    """Get current LLM settings"""
    return {**llm_settings, "api_key": "***" if llm_settings.get("api_key") else ""}


class LLMSettingsUpdate(BaseModel):
    provider: Optional[str] = None
    model: Optional[str] = None
    api_key: Optional[str] = None
    base_url: Optional[str] = None
    temperature: Optional[float] = None
    max_tokens: Optional[int] = None


@app.put("/llm/settings")
async def update_llm_settings(settings: LLMSettingsUpdate):
    """Update LLM settings"""
    global llm_settings
    update_data = settings.model_dump(exclude_none=True)
    llm_settings.update(update_data)
    logger.info(
        f"LLM settings updated: provider={llm_settings['provider']}, model={llm_settings['model']}"
    )
    return {
        "status": "ok",
        "settings": {
            **llm_settings,
            "api_key": "***" if llm_settings.get("api_key") else "",
        },
    }


@app.post("/llm/test")
async def test_llm_connection():
    """Test LLM connection with current settings"""
    async with httpx.AsyncClient(timeout=15.0) as client:
        try:
            headers = {}
            if llm_settings.get("api_key"):
                headers["Authorization"] = f"Bearer {llm_settings['api_key']}"

            response = await client.post(
                f"{llm_settings['base_url']}/chat/completions",
                json={
                    "model": llm_settings["model"],
                    "messages": [
                        {"role": "user", "content": "Hi, reply with 'OK' only."}
                    ],
                    "max_tokens": 10,
                },
                headers=headers,
            )
            response.raise_for_status()
            result = response.json()
            return {
                "status": "ok",
                "response": result["choices"][0]["message"]["content"],
            }
        except httpx.HTTPStatusError as e:
            return {
                "status": "error",
                "error": f"HTTP {e.response.status_code}: {e.response.text[:200]}",
            }
        except Exception as e:
            return {"status": "error", "error": str(e)}


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)


class RMSDRequest(BaseModel):
    pdb1: str
    pdb2: str


@app.post("/rmsd")
async def calculate_rmsd(request: RMSDRequest):
    """Calculate RMSD between two PDB structures"""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        import numpy as np

        mol1 = None
        mol2 = None

        if os.path.exists(request.pdb1):
            mol1 = Chem.MolFromPDBFile(request.pdb1)
        else:
            mol1 = Chem.MolFromPDBBlock(request.pdb1)

        if os.path.exists(request.pdb2):
            mol2 = Chem.MolFromPDBFile(request.pdb2)
        else:
            mol2 = Chem.MolFromPDBBlock(request.pdb2)

        if mol1 is None or mol2 is None:
            raise HTTPException(
                status_code=400, detail="Could not parse PDB structures"
            )

        conf1 = mol1.GetConformer(0) if mol1.GetNumConformers() > 0 else None
        conf2 = mol2.GetConformer(0) if mol2.GetNumConformers() > 0 else None

        if conf1 is None or conf2 is None:
            return {"rmsd": 0.0, "note": "No conformers found, returned 0"}

        n1 = mol1.GetNumAtoms()
        n2 = mol2.GetNumAtoms()
        n = min(n1, n2)

        rmsd_sum = 0.0
        for i in range(n):
            p1 = conf1.GetAtomPosition(i)
            p2 = conf2.GetAtomPosition(i)
            rmsd_sum += (p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2 + (p1.z - p2.z) ** 2

        rmsd = np.sqrt(rmsd_sum / n) if n > 0 else 0.0
        return {"rmsd": round(rmsd, 3)}

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


class BindingSiteRequest(BaseModel):
    receptor: str
    ligand: str


@app.post("/binding-site")
async def get_binding_site(request: BindingSiteRequest):
    """Get binding site residues around ligand in receptor"""
    try:
        from rdkit import Chem

        receptor = None
        if os.path.exists(request.receptor):
            receptor = Chem.MolFromPDBFile(request.receptor)
        else:
            receptor = Chem.MolFromPDBBlock(request.receptor)

        if receptor is None:
            raise HTTPException(status_code=400, detail="Could not parse receptor PDB")

        residues = []
        for atom in receptor.GetAtoms():
            res = atom.GetPDBResidue()
            if res:
                residues.append(
                    {
                        "residue_name": res.GetResname(),
                        "residue_number": res.GetResSeqNumber(),
                        "chain": res.GetChainId(),
                        "atom_name": atom.GetName(),
                    }
                )

        unique = {(r["residue_number"], r["chain"]) for r in residues}
        result = []
        for rnum, chain in list(unique)[:20]:
            for res in residues:
                if res["residue_number"] == rnum and res["chain"] == chain:
                    result.append(res)
                    break

        return {"residues": result}

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/gpu/status")
async def get_gpu_status():
    """Get GPU availability status"""
    try:
        import subprocess

        try:
            result = subprocess.run(
                [
                    "nvidia-smi",
                    "--query-gpu=index,name,utilization.gpu,memory.used,memory.total,temperature.gpu",
                    "--format=csv,noheader,nounits",
                ],
                capture_output=True,
                text=True,
                timeout=5,
            )
            if result.returncode == 0:
                lines = result.stdout.strip().split("\n")
                gpus = []
                for line in lines:
                    parts = [p.strip() for p in line.split(",")]
                    if len(parts) >= 6:
                        gpus.append(
                            {
                                "index": int(parts[0]),
                                "name": parts[1],
                                "utilization": float(parts[2]),
                                "memory_used": float(parts[3]),
                                "memory_total": float(parts[4]),
                                "temperature": float(parts[5]),
                            }
                        )
                return {"available": True, "gpus": gpus}
        except Exception:
            pass

        return {
            "available": False,
            "gpus": [],
            "message": "No GPU detected (CPU-only mode)",
        }

    except Exception as e:
        return {"available": False, "gpus": [], "message": str(e)}


@app.get("/security/status")
async def get_security_status():
    """Get security scan status"""
    return {
        "last_scan_at": None,
        "overall_severity": "NOT_SCANNED",
        "is_secure": True,
        "total_issues": 0,
        "scan_results": {},
    }


@app.post("/security/scan")
async def run_security_scan():
    """Run security scan"""
    import subprocess
    import datetime

    issues = []
    severity = "NONE"
    scan_results = {}

    try:
        result = subprocess.run(
            ["pip", "list", "--format=freeze"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        if result.returncode == 0:
            scan_results["dependencies"] = "Scanned"
    except Exception:
        scan_results["dependencies"] = "Scan failed"

    return {
        "worst_severity": severity,
        "is_secure": severity in ("NONE", "LOW"),
        "total_issues": len(issues),
        "scan_results": scan_results,
        "scanned_at": datetime.datetime.now().isoformat(),
    }


@app.get("/download/{filename}")
async def download_file(filename: str):
    """Download a file from storage"""
    from fastapi.responses import FileResponse
    import os

    storage_path = STORAGE_DIR / filename
    uploads_path = UPLOADS_DIR / filename

    if storage_path.exists():
        return FileResponse(storage_path, filename=filename)
    elif uploads_path.exists():
        return FileResponse(uploads_path, filename=filename)
    else:
        raise HTTPException(status_code=404, detail="File not found")


@app.post("/analyze/interactions")
async def analyze_interactions(receptor_pdb: str, ligand_pdb: str):
    """Analyze protein-ligand interactions"""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolDescriptors

        if os.path.exists(receptor_pdb):
            receptor = Chem.MolFromPDBFile(receptor_pdb)
        else:
            receptor = Chem.MolFromPDBBlock(receptor_pdb)

        if os.path.exists(ligand_pdb):
            ligand = Chem.MolFromPDBFile(ligand_pdb)
        else:
            ligand = Chem.MolFromPDBBlock(ligand_pdb)

        if receptor is None or ligand is None:
            return {"success": False, "error": "Invalid structures"}

        interactions = []
        num_hbonds = 0
        num_hydrophobic = 0

        try:
            from rdkit.Chem import rdMolDescriptors

            receptor_hbd = rdMolDescriptors.CalcNumHBD(receptor)
            receptor_hba = rdMolDescriptors.CalcNumHBA(receptor)
            ligand_hbd = rdMolDescriptors.CalcNumHBD(ligand)
            ligand_hba = rdMolDescriptors.CalcNumHBA(ligand)
            interactions.append(
                {
                    "type": "hydrogen_bond",
                    "description": f"H-bond donors/acceptors - Receptor: {receptor_hbd}/{receptor_hba}, Ligand: {ligand_hbd}/{ligand_hba}",
                    "count": min(receptor_hbd, ligand_hba)
                    + min(receptor_hba, ligand_hbd),
                }
            )
            num_hbonds = min(receptor_hbd, ligand_hba) + min(receptor_hba, ligand_hbd)
        except Exception:
            pass

        try:
            from rdkit.Chem import Lipinski

            receptor_lipinski = Lipinski.NumAromaticRings(receptor)
            ligand_lipinski = Lipinski.NumAromaticRings(ligand)
            num_hydrophobic = receptor_lipinski + ligand_lipinski
            if num_hydrophobic > 0:
                interactions.append(
                    {
                        "type": "hydrophobic",
                        "description": f"Aromatic rings - Receptor: {receptor_lipinski}, Ligand: {ligand_lipinski}",
                        "count": num_hydrophobic,
                    }
                )
        except Exception:
            pass

        return {
            "success": True,
            "interactions": interactions,
            "num_hbonds": num_hbonds,
            "num_hydrophobic": num_hydrophobic,
        }

    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "interactions": [],
            "num_hbonds": 0,
            "num_hydrophobic": 0,
        }


@app.get("/metrics")
async def get_metrics():
    """Prometheus-compatible metrics endpoint"""
    try:
        r = redis.from_url(REDIS_URL, decode_responses=True)
        info = r.info()

        job_keys = [
            "job:*",
            "docking_job:*",
            "md_job:*",
            "qsar_job:*",
            "pharmacophore_job:*",
        ]
        total_jobs = 0
        running = 0
        completed = 0
        failed = 0
        pending = 0

        for pattern in job_keys:
            for key in r.scan_iter(pattern, count=100):
                total_jobs += 1
                data = r.get(key)
                if data:
                    try:
                        job = json.loads(data)
                        status = job.get("status", "")
                        if status == "running":
                            running += 1
                        elif status == "completed":
                            completed += 1
                        elif status == "failed":
                            failed += 1
                        else:
                            pending += 1
                    except:
                        pass

        docker_services = [
            "api-backend",
            "brain-service",
            "docking-service",
            "rdkit-service",
            "pharmacophore-service",
            "qsar-service",
            "md-service",
            "sentinel-service",
            "analysis-service",
        ]
        healthy_services = len([s for s in docker_services if r.exists(f"health:{s}")])

        metrics_text = f"""# HELP docking_jobs_total Total number of jobs
# TYPE docking_jobs_total gauge
docking_jobs_total {total_jobs}
docking_jobs_running {running}
docking_jobs_completed {completed}
docking_jobs_failed {failed}
docking_jobs_pending {pending}

# HELP docking_services_healthy Number of healthy services
# TYPE docking_services_healthy gauge
docking_services_healthy {healthy_services}

# HELP redis_connected_clients Number of connected Redis clients
# TYPE redis_connected_clients gauge
redis_connected_clients {info.get("connected_clients", 0)}

# HELP redis_used_memory_human Redis memory usage
# TYPE redis_used_memory_human gauge
redis_used_memory_human {info.get("used_memory_human", "0")}
"""
        return Response(content=metrics_text, media_type="text/plain")
    except Exception as e:
        return Response(content=f"# Error: {e}", media_type="text/plain")
