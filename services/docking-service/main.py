"""
Docking Service - AutoDock Vina execution API
"""

import os
import logging
from typing import Optional
from datetime import datetime
from contextlib import asynccontextmanager

import redis
import json

from fastapi import FastAPI, HTTPException, Form, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("docking-service")

STORAGE_DIR = os.environ.get("STORAGE_DIR", "/app/storage")
os.makedirs(STORAGE_DIR, exist_ok=True)


class DockingRequest(BaseModel):
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


class DockingJobRequest(BaseModel):
    job_id: str
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


@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("Docking Service starting up...")
    yield
    logger.info("Docking Service shutting down...")


app = FastAPI(
    title="Docking Service API",
    description="AutoDock Vina molecular docking",
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
        "service": "docking-service",
        "timestamp": datetime.now().isoformat(),
    }


@app.get("/")
def root():
    return {"service": "docking-service", "version": "2.0.0"}


@app.get("/vina/check")
def check_vina():
    """Check if Vina is available"""
    from docking_engine import check_vina as _check_vina

    available = _check_vina()
    return {"vina_available": available}


@app.get("/gnina/check")
def check_gnina():
    """Check if GNINA is available"""
    from docking_engine import check_gnina as _check_gnina

    available = _check_gnina()
    return {"gnina_available": available}


@app.post("/dock/gnina")
def run_gnina_docking(request: DockingRequest):
    """Run GNINA docking with CNN scoring"""
    try:
        from docking_engine import run_gnina_docking as _run_gnina_docking

        result = _run_gnina_docking(
            receptor_path=request.receptor_path,
            ligand_path=request.ligand_path,
            center_x=request.center_x,
            center_y=request.center_y,
            center_z=request.center_z,
            size_x=request.size_x,
            size_y=request.size_y,
            size_z=request.size_z,
            exhaustiveness=request.exhaustiveness,
            num_modes=request.num_modes,
            output_dir=STORAGE_DIR,
        )

        return result

    except Exception as e:
        logger.error(f"GNINA docking error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/dock")
def run_docking(request: DockingRequest):
    """Run molecular docking"""
    try:
        from docking_engine import run_docking as _run_docking

        result = _run_docking(
            receptor_path=request.receptor_path,
            ligand_path=request.ligand_path,
            engine=request.engine,
            center_x=request.center_x,
            center_y=request.center_y,
            center_z=request.center_z,
            size_x=request.size_x,
            size_y=request.size_y,
            size_z=request.size_z,
            exhaustiveness=request.exhaustiveness,
            num_modes=request.num_modes,
            output_dir=STORAGE_DIR,
        )

        return result

    except Exception as e:
        logger.error(f"Docking error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/dock/async")
def start_async_docking(request: DockingJobRequest, background_tasks: BackgroundTasks):
    """Start an asynchronous docking job"""
    try:
        from docking_engine import run_docking as _run_docking

        redis_url = os.environ.get("REDIS_URL", "redis://localhost:6379")
        r = redis.from_url(redis_url)

        job_data = {
            "job_id": request.job_id,
            "status": "queued",
            "receptor_path": request.receptor_path,
            "ligand_path": request.ligand_path,
            "engine": request.engine,
            "center": json.dumps([request.center_x, request.center_y, request.center_z]),
            "size": json.dumps([request.size_x, request.size_y, request.size_z]),
            "exhaustiveness": str(request.exhaustiveness),
            "num_modes": str(request.num_modes),
        }

        r.hset(f"docking_job:{request.job_id}", mapping=job_data)
        r.rpush("docking_queue", request.job_id)

        background_tasks.add_task(process_docking_job, request.job_id, job_data)

        return {"job_id": request.job_id, "status": "queued"}

    except Exception as e:
        logger.error(f"Failed to queue docking job: {e}")
        raise HTTPException(status_code=500, detail=str(e))


def process_docking_job(job_id: str, job_data: dict):
    """Background task to process docking job"""
    try:
        redis_url = os.environ.get("REDIS_URL", "redis://localhost:6379")
        r = redis.from_url(redis_url)

        def is_cancelled():
            job = r.hgetall(f"docking_job:{job_id}")
            return job.get("status") == "cancelled"

        if is_cancelled():
            logger.info(f"Job {job_id} was cancelled before starting")
            return

        r.hset(f"docking_job:{job_id}", mapping={**job_data, "status": "running"})

        from docking_engine import run_docking as _run_docking

        result = _run_docking(
            receptor_path=job_data["receptor_path"],
            ligand_path=job_data["ligand_path"],
            engine=job_data["engine"],
            center_x=job_data["center"][0],
            center_y=job_data["center"][1],
            center_z=job_data["center"][2],
            size_x=job_data["size"][0],
            size_y=job_data["size"][1],
            size_z=job_data["size"][2],
            exhaustiveness=job_data["exhaustiveness"],
            num_modes=job_data["num_modes"],
            output_dir=STORAGE_DIR,
        )

        if is_cancelled():
            logger.info(f"Job {job_id} was cancelled during execution")
            return

        if result.get("success"):
            r.hset(
                f"docking_job:{job_id}",
                mapping={**job_data, "status": "completed", "result": json.dumps(result)},
            )
        else:
            r.hset(
                f"docking_job:{job_id}",
                mapping={**job_data, "status": "failed", "error": result.get("error")},
            )

    except Exception as e:
        logger.error(f"Docking job {job_id} failed: {e}")
        r.hset(
            f"docking_job:{job_id}",
            mapping={**job_data, "status": "failed", "error": str(e)},
        )


@app.get("/dock/{job_id}/status")
def get_docking_status(job_id: str):
    """Get docking job status"""
    try:
        redis_url = os.environ.get("REDIS_URL", "redis://localhost:6379")
        r = redis.from_url(redis_url)

        data = r.hgetall(f"docking_job:{job_id}")
        if data:
            result = dict(data)
            if "center" in result:
                result["center"] = json.loads(result["center"])
            if "size" in result:
                result["size"] = json.loads(result["size"])
            if "exhaustiveness" in result:
                result["exhaustiveness"] = int(result["exhaustiveness"])
            if "num_modes" in result:
                result["num_modes"] = int(result["num_modes"])
            if "result" in result:
                result["result"] = json.loads(result["result"])
            return result
        return {"job_id": job_id, "status": "not_found"}

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/dock/{job_id}/result")
def get_docking_result(job_id: str):
    """Get docking job result"""
    try:
        redis_url = os.environ.get("REDIS_URL", "redis://localhost:6379")
        r = redis.from_url(redis_url)

        data = r.hgetall(f"docking_job:{job_id}")
        if data:
            job_data = dict(data)
            if job_data.get("status") == "completed":
                return json.loads(job_data.get("result", "{}"))
            elif job_data.get("status") == "failed":
                raise HTTPException(
                    status_code=400, detail=job_data.get("error", "Job failed")
                )
            return {"status": job_data.get("status")}
        raise HTTPException(status_code=404, detail="Job not found")

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/dock/{job_id}/cancel")
async def cancel_docking(job_id: str):
    """Cancel a running docking job"""
    redis_url = os.environ.get("REDIS_URL", "redis://localhost:6379")
    r = redis.from_url(redis_url)
    job_key = f"docking_job:{job_id}"

    if not r.exists(job_key):
        raise HTTPException(status_code=404, detail="Job not found")

    job_data = r.hgetall(job_key)
    current_status = job_data.get("status", "unknown")

    if current_status in ["completed", "failed", "cancelled"]:
        return {"job_id": job_id, "status": current_status}

    r.hset(job_key, mapping={"status": "cancelled"})
    r.lpush("docking_cancelled", job_id)

    return {"job_id": job_id, "status": "cancelled"}


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8002)
