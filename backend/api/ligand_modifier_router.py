"""
Ligand Modifier API Router
FastAPI endpoints for ligand modification pipeline.
"""
import uuid
import logging
from typing import List, Optional, Literal
from datetime import datetime

from fastapi import APIRouter, BackgroundTasks, HTTPException
from pydantic import BaseModel

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/ligand-modifier", tags=["ligand-modifier"])


class LigandModifierRequest(BaseModel):
    parent_smiles: str
    receptor_pdb: str
    mode: Literal["similarity_search", "prompt_based", "autonomous"]
    prompt: Optional[str] = None
    database: Literal["pubchem", "chembl"] = "pubchem"
    similarity_threshold: float = 0.85
    max_variants: int = 50
    docking_exhaustiveness: int = 8


class LigandModifierJob:
    def __init__(self, job_id: str):
        self.job_id = job_id
        self.status = "queued"
        self.results: List[dict] = []
        self.error: Optional[str] = None
        self.progress: float = 0.0
        self.created_at = datetime.now().isoformat()


JOBS: dict = {}


def _run_pipeline(req: LigandModifierRequest, job: LigandModifierJob):
    """Execute the full ligand modification pipeline."""
    try:
        from crew.tools.ligand_modification_tools import (
            SimilaritySearchRequest,
            ModificationPlan,
            fetch_similar_ligands,
            parse_modification_prompt,
            apply_rdkit_transformations,
            validate_and_filter,
        )

        # 1. Parse prompt or build plan
        job.status = "parsing"
        job.progress = 0.1

        if req.mode == "similarity_search":
            plan = ModificationPlan(mode="similarity_search", max_variants=req.max_variants)
        elif req.mode == "prompt_based" and req.prompt:
            plan = parse_modification_prompt(req.prompt, req.parent_smiles)
            plan.mode = "prompt_based"
            plan.max_variants = req.max_variants
        else:
            plan = ModificationPlan(
                mode="autonomous",
                strategy="optimize_binding",
                max_variants=req.max_variants,
                docking_params={"exhaustiveness": req.docking_exhaustiveness},
            )

        # 2. Fetch or generate variants
        job.status = "fetching" if req.mode == "similarity_search" else "transforming"
        job.progress = 0.3

        if req.mode == "similarity_search":
            similar = fetch_similar_ligands(
                SimilaritySearchRequest(
                    query_smiles=req.parent_smiles,
                    database=req.database,
                    similarity_threshold=req.similarity_threshold,
                    max_results=req.max_variants,
                )
            )
            variants = [
                {
                    "parent_smiles": req.parent_smiles,
                    "modified_smiles": s["smiles"],
                    "applied_transform": "database_fetch",
                    "source": "database",
                }
                for s in similar
            ]
        else:
            variants = apply_rdkit_transformations([req.parent_smiles], plan)

        # 3. Validate and filter
        job.status = "validating"
        job.progress = 0.6

        filtered = validate_and_filter(variants, plan.property_constraints)

        if not filtered:
            job.status = "completed"
            job.progress = 1.0
            job.results = []
            return

        # 4. Calculate properties for results (docking is optional)
        job.status = "docking"
        job.progress = 0.8

        results = []
        for v in filtered[: req.max_variants]:
            results.append(
                {
                    "parent_smiles": v.get("parent_smiles", req.parent_smiles),
                    "modified_smiles": v["modified_smiles"],
                    "applied_transform": v.get("applied_transform", "unknown"),
                    "properties": v.get("properties", {}),
                    "docking_score": None,
                    "delta_score": None,
                    "source": v.get("source", "transformation"),
                }
            )

        # Sort by MW as default ranking
        results.sort(key=lambda x: x["properties"].get("mw", 999))

        job.results = results
        job.status = "completed"
        job.progress = 1.0

    except Exception as e:
        job.status = "failed"
        job.error = str(e)
        logger.exception(f"Ligand modifier pipeline failed: {e}")


@router.post("/optimize")
async def start_optimization(req: LigandModifierRequest, bg_tasks: BackgroundTasks):
    """Start ligand modification job."""
    job_id = f"ligmod_{uuid.uuid4().hex[:8]}"
    job = LigandModifierJob(job_id)
    JOBS[job_id] = job

    bg_tasks.add_task(_run_pipeline, req, job)

    return {"job_id": job_id, "status": "queued", "message": "Optimization started"}


@router.get("/status/{job_id}")
def get_job_status(job_id: str):
    """Get job status and results."""
    if job_id not in JOBS:
        raise HTTPException(404, "Job not found")

    job = JOBS[job_id]
    return {
        "job_id": job.job_id,
        "status": job.status,
        "progress": job.progress,
        "results": job.results,
        "error": job.error,
    }


@router.delete("/cancel/{job_id}")
def cancel_job(job_id: str):
    """Cancel a running job."""
    if job_id not in JOBS:
        raise HTTPException(404, "Job not found")

    job = JOBS[job_id]
    if job.status in ["completed", "failed"]:
        return {"message": "Job already finished"}

    job.status = "cancelled"
    return {"message": "Job cancelled"}
