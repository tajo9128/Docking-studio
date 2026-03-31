from fastapi import FastAPI, UploadFile, File, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from typing import Optional, List
from pathlib import Path
import asyncio
import subprocess
import uuid
import json
import os
from datetime import datetime

app = FastAPI(title="BioDockify", version="2.3.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

DATA_DIR = Path("data")
JOBS_DIR = DATA_DIR / "jobs"
JOBS_DIR.mkdir(parents=True, exist_ok=True)

class DockingJob(BaseModel):
    receptor_file: str
    ligand_smiles: Optional[str] = None
    ligand_file: Optional[str] = None
    exhaustiveness: int = 32
    n_poses: int = 10

class MDJob(BaseModel):
    pdb_file: str
    forcefield: str = "amber14-all.xml"
    duration_ns: float = 10.0
    temperature_k: float = 300.0

class AIRequest(BaseModel):
    message: str
    provider: str = "openai"
    model: Optional[str] = None

class AIConfig(BaseModel):
    provider: str
    api_key: Optional[str] = None
    base_url: Optional[str] = None
    model: str

ai_config_store = {}

@app.get("/")
async def root():
    return {"message": "BioDockify API", "version": "2.3.0"}

@app.get("/health")
async def health():
    return {"status": "healthy", "timestamp": datetime.now().isoformat()}

@app.post("/api/docking/jobs")
async def create_docking_job(job: DockingJob):
    job_id = str(uuid.uuid4())[:8]
    job_dir = JOBS_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)
    
    return {"job_id": job_id, "status": "created", "message": "Docking job submitted"}

@app.get("/api/docking/jobs/{job_id}")
async def get_docking_job(job_id: str):
    job_file = JOBS_DIR / job_id / "result.json"
    if job_file.exists():
        return json.loads(job_file.read_text())
    return {"job_id": job_id, "status": "pending"}

@app.post("/api/md/jobs")
async def create_md_job(job: MDJob):
    job_id = str(uuid.uuid4())[:8]
    return {"job_id": job_id, "status": "created"}

@app.get("/api/md/gpu-info")
async def get_gpu_info():
    return {
        "cuda_available": False,
        "opencl_available": False,
        "platform": "CPU",
        "message": "GPU detection will run inside container"
    }

@app.post("/api/ai/chat")
async def chat(request: AIRequest):
    return {
        "response": f"BioDockify AI using {request.provider}: {request.message}",
        "provider": request.provider
    }

@app.put("/api/ai/config")
async def update_ai_config(config: AIConfig):
    ai_config_store[config.provider] = config.dict()
    return {"status": "updated", "provider": config.provider}

@app.get("/api/ai/config/{provider}")
async def get_ai_config(provider: str):
    if provider in ai_config_store:
        return ai_config_store[provider]
    return {"provider": provider, "configured": False}

@app.get("/api/ai/providers")
async def list_providers():
    return {
        "providers": [
            {"id": "openai", "name": "OpenAI"},
            {"id": "claude", "name": "Claude"},
            {"id": "gemini", "name": "Gemini"},
            {"id": "mistral", "name": "Mistral"},
            {"id": "deepseek", "name": "DeepSeek"},
            {"id": "qwen", "name": "Qwen"},
            {"id": "siliconflow", "name": "SiliconFlow"},
            {"id": "openrouter", "name": "OpenRouter"},
            {"id": "ollama", "name": "Ollama"}
        ]
    }

@app.get("/api/stats")
async def get_stats():
    return {
        "total_jobs": 0,
        "completed_jobs": 0,
        "active_jobs": 0
    }

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
