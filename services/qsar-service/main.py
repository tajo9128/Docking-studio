"""
QSAR Service - ML-based QSAR modeling
"""

import os
import logging
import asyncio
import uuid
from pathlib import Path
from typing import Optional, List, Dict, Any
from datetime import datetime
from contextlib import asynccontextmanager
import redis

from fastapi import FastAPI, HTTPException, BackgroundTasks, UploadFile, File, Form
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

from .descriptors import (
    descriptors_for_smiles_list,
    descriptors_upload,
    DESCRIPTOR_GROUPS,
)
from .qsar_engine import train_model, predict_single, predict_batch
from .model_store import list_models, get_model, delete_model

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("qsar-service")

REDIS_URL = os.getenv("REDIS_URL", "redis://redis:6379")
STORAGE_DIR = Path("/app/storage/qsar")
STORAGE_DIR.mkdir(parents=True, exist_ok=True)
TRAINING_TIMEOUT = int(os.getenv("TRAINING_TIMEOUT", "300"))

try:
    redis_client = redis.from_url(REDIS_URL, decode_responses=True)
    redis_client.ping()
    redis_available = True
except Exception:
    redis_client = None
    redis_available = False
    logger.warning(
        "Redis not available — training jobs will not persist across restarts"
    )


class DescriptorsRequest(BaseModel):
    smiles: List[str]
    groups: Optional[List[str]] = None


class TrainRequest(BaseModel):
    X: List[List[float]]
    y: List[float]
    feature_names: List[str]
    model_type: str
    model_name: str
    activity_column: str
    descriptor_groups: List[str]
    cv_folds: int = 5
    model_params: Optional[Dict[str, Any]] = None


class PredictRequest(BaseModel):
    model_id: str
    smiles: str


class BatchPredictRequest(BaseModel):
    model_id: str
    smiles_list: List[str]


def _set_job_status(
    job_id: str, status: str, result: Optional[Any] = None, error: Optional[str] = None
):
    if redis_client:
        import json

        data = {
            "status": status,
            "updated_at": datetime.now().isoformat(),
        }
        if result is not None:
            data["result"] = result
        if error is not None:
            data["error"] = error
        redis_client.setex(f"qsar_job:{job_id}", 3600, json.dumps(data))


@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("=" * 60)
    logger.info("QSAR Service starting up...")
    logger.info(f"Redis available: {redis_available}")
    logger.info("=" * 60)
    yield
    logger.info("QSAR Service shutting down...")


app = FastAPI(
    title="QSAR Service API",
    description="ML-based QSAR modeling with RDKit descriptors",
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


@app.get("/health")
def health():
    return {
        "status": "healthy",
        "service": "qsar-service",
        "redis": redis_available,
        "timestamp": datetime.now().isoformat(),
    }


@app.get("/")
def root():
    return {"service": "qsar-service", "version": "2.0.0"}


@app.get("/descriptor-groups")
def available_descriptor_groups():
    return {
        "groups": list(DESCRIPTOR_GROUPS.keys()),
        "descriptors": {k: DESCRIPTOR_GROUPS[k] for k in DESCRIPTOR_GROUPS},
    }


@app.post("/descriptors")
def calculate_descriptors_endpoint(request: DescriptorsRequest):
    """Calculate descriptors for a list of SMILES"""
    try:
        descriptors, valid_smiles, failed = descriptors_for_smiles_list(
            request.smiles,
            request.groups,
        )
        return {
            "success": True,
            "descriptors": descriptors,
            "valid_smiles": valid_smiles,
            "failed_smiles": failed,
            "n_valid": len(valid_smiles),
            "n_failed": len(failed),
        }
    except Exception as e:
        logger.error(f"Descriptor calculation error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/descriptors/upload")
def descriptors_upload_endpoint(
    file: UploadFile = File(...),
    smiles_col: str = Form(...),
    activity_col: str = Form(...),
    groups: Optional[str] = Form(None),
):
    """Upload CSV dataset and calculate descriptors"""
    try:
        content = file.file.read().decode("utf-8")
        group_list = groups.split(",") if groups else ["all"]

        result = descriptors_upload(
            csv_content=content,
            smiles_col=smiles_col,
            activity_col=activity_col,
            descriptor_groups=group_list,
        )
        return result
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"CSV upload error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


def _run_training(job_id: str, train_data: TrainRequest):
    try:
        _set_job_status(job_id, "running")
        result = train_model(
            X=train_data.X,
            y=train_data.y,
            feature_names=train_data.feature_names,
            model_type=train_data.model_type,
            model_name=train_data.model_name,
            activity_column=train_data.activity_column,
            descriptor_groups=train_data.descriptor_groups,
            cv_folds=train_data.cv_folds,
            model_params=train_data.model_params,
        )
        _set_job_status(job_id, "completed", result=result)
    except Exception as e:
        logger.error(f"Training job {job_id} failed: {e}")
        _set_job_status(job_id, "failed", error=str(e))


@app.post("/train")
def train(
    request: TrainRequest,
    background_tasks: BackgroundTasks,
):
    """Start async model training"""
    job_id = str(uuid.uuid4())[:8]
    _set_job_status(job_id, "pending")
    background_tasks.add_task(_run_training, job_id, request)
    return {"job_id": job_id, "status": "pending", "message": "Training started"}


@app.get("/train/{job_id}/status")
def training_status(job_id: str):
    """Poll training job status"""
    if not redis_client:
        raise HTTPException(status_code=503, detail="Redis not available")
    import json

    data = redis_client.get(f"qsar_job:{job_id}")
    if data is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")
    return json.loads(data)


@app.get("/train/{job_id}/results")
def training_results(job_id: str):
    """Fetch completed training results"""
    if not redis_client:
        raise HTTPException(status_code=503, detail="Redis not available")
    import json

    data = redis_client.get(f"qsar_job:{job_id}")
    if data is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")
    result = json.loads(data)
    if result["status"] != "completed":
        return {"status": result["status"], "error": result.get("error")}
    return {"status": "completed", "result": result.get("result")}


@app.post("/predict")
def predict(request: PredictRequest):
    """Predict activity for a single SMILES"""
    try:
        result = predict_single(request.model_id, request.smiles)
        return {"success": True, **result}
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Prediction error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/predict/batch")
def predict_batch_endpoint(request: BatchPredictRequest):
    """Predict activity for multiple SMILES"""
    try:
        result = predict_batch(request.model_id, request.smiles_list)
        return {"success": True, **result}
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Batch prediction error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/models")
def models_list():
    """List all saved models"""
    return {"models": list_models()}


@app.get("/models/{model_id}")
def model_info(model_id: str):
    """Get model metadata"""
    model = get_model(model_id)
    if model is None:
        raise HTTPException(status_code=404, detail=f"Model {model_id} not found")
    return model


@app.delete("/models/{model_id}")
def model_delete(model_id: str):
    """Delete a saved model"""
    deleted = delete_model(model_id)
    if not deleted:
        raise HTTPException(status_code=404, detail=f"Model {model_id} not found")
    return {"success": True, "model_id": model_id}


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8005)
