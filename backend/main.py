from fastapi import FastAPI, HTTPException, UploadFile, File, Form
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from typing import Optional, List
import os
import uuid
import json
import logging
import time
import threading
from datetime import datetime
from analysis import (
    analyze_pose, 
    calculate_rmsd, 
    calculate_rmsd_files,
    calculate_advanced_interactions,
    get_binding_site_residues
)
from db import (
    init_db, 
    create_job, 
    update_job_status, 
    add_docking_result,
    add_interaction,
    get_job, 
    get_all_jobs,
    get_docking_results,
    get_interactions,
    delete_job
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)-8s | %(name)s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger("docking-studio")

init_db()
logger.info("Docking Studio backend initialized")

app = FastAPI(
    title="Docking Studio API",
    description="Backend API for molecular docking analysis",
    version="1.0.0"
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

STORAGE_DIR = os.path.join(os.path.dirname(__file__), "storage")
os.makedirs(STORAGE_DIR, exist_ok=True)

# Setup static files
STATIC_DIR = os.path.join(os.path.dirname(__file__), "static")
if os.path.exists(STATIC_DIR):
    app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")


@app.get("/")
async def root():
    """Serve the web interface"""
    index_path = os.path.join(STATIC_DIR, "index.html")
    if os.path.exists(index_path):
        return FileResponse(index_path)
    return {"message": "Docking Studio API", "version": "1.0.0"}


@app.on_event("startup")
async def startup_event():
    """Print startup information when server starts"""
    print("")
    print("=" * 60)
    print("  🧬 Docking Studio - Backend Started")
    print("=" * 60)
    print("")
    print("  📚 API Documentation (Swagger UI):")
    print("     ➤ http://localhost:8000/docs")
    print("")
    print("  📖 Alternative API Docs (ReDoc):")
    print("     ➤ http://localhost:8000/redoc")
    print("")
    print("  ✅ Health Check:")
    print("     ➤ http://localhost:8000/health")
    print("")
    print("  🔐 Security Status:")
    print("     ➤ http://localhost:8000/security/status")
    print("")
    print("  🤖 Ollama AI (if enabled):")
    print("     ➤ http://localhost:11434")
    print("")
    print("=" * 60)
    print("  🎯 Quick Start for Students:")
    print("     1. Open http://localhost:8000/docs in browser")
    print("     2. Read the API documentation")
    print("     3. Try the /dock endpoints")
    print("=" * 60)
    print("")





class PoseRequest(BaseModel):
    receptor: str
    ligand: str


class JobRequest(BaseModel):
    job_name: str
    receptor_path: str
    ligand_path: str
    engine: str = "vina"


class RMSDRequest(BaseModel):
    pdb1: str
    pdb2: str


@app.get("/health")
def health():
    logger.debug("Health check requested")
    return {"status": "healthy", "timestamp": datetime.now().isoformat()}


@app.post("/analyze")
def analyze(req: PoseRequest):
    """Analyze pose interactions"""
    logger.info("Pose analysis requested")
    try:
        result = analyze_pose(req.receptor, req.ligand)
        logger.info(f"Analysis complete: {result.get('hbond_count', 0)} hbonds, {result.get('hydrophobic_count', 0)} hydrophobic")
        return result
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/analyze/advanced")
def analyze_advanced(req: PoseRequest):
    """Advanced interaction analysis"""
    logger.info("Advanced analysis requested")
    try:
        result = calculate_advanced_interactions(req.receptor, req.ligand)
        logger.info(f"Advanced analysis complete")
        return result
    except Exception as e:
        logger.error(f"Advanced analysis failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/rmsd")
def rmsd(req: RMSDRequest):
    """Calculate RMSD between two PDB structures"""
    try:
        value = calculate_rmsd(req.pdb1, req.pdb2)
        return {"rmsd": value}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/rmsd/file")
def rmsd_file(pdb1_path: str = Form(...), pdb2_path: str = Form(...)):
    """Calculate RMSD between two PDB files"""
    try:
        value = calculate_rmsd_files(pdb1_path, pdb2_path)
        return {"rmsd": value}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/binding-site")
def binding_site(req: PoseRequest):
    """Get binding site residues"""
    try:
        residues = get_binding_site_residues(req.receptor, req.ligand)
        return {"residues": residues}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/jobs")
def create_new_job(job: JobRequest):
    """Create a new docking job"""
    logger.info(f"Creating job: {job.job_name} (engine: {job.engine})")
    job_uuid = str(uuid.uuid4())
    
    success = create_job(
        job_uuid=job_uuid,
        job_name=job.job_name,
        receptor_file=job.receptor_path,
        ligand_file=job.ligand_path,
        engine=job.engine
    )
    
    if success:
        logger.info(f"Job created successfully: {job_uuid}")
        return {"job_uuid": job_uuid, "status": "created"}
    else:
        logger.error(f"Failed to create job: {job.job_name}")
        raise HTTPException(status_code=500, detail="Failed to create job")


@app.get("/jobs")
def list_jobs(limit: int = 50):
    """List all jobs"""
    jobs = get_all_jobs(limit)
    return {"jobs": jobs}


@app.get("/jobs/{job_uuid}")
def get_job_info(job_uuid: str):
    """Get job details"""
    job = get_job(job_uuid)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    return job


@app.post("/jobs/{job_uuid}/status")
def update_job(job_uuid: str, status: str, binding_energy: Optional[float] = None, 
             confidence: Optional[float] = None):
    """Update job status"""
    success = update_job_status(job_uuid, status, binding_energy, confidence)
    if success:
        return {"status": "updated"}
    else:
        raise HTTPException(status_code=500, detail="Failed to update job")


@app.post("/jobs/{job_uuid}/results")
def add_result(job_uuid: str, pose_id: int, ligand_name: str, 
              vina_score: Optional[float] = None,
              gnina_score: Optional[float] = None,
              rf_score: Optional[float] = None,
              pdb_data: Optional[str] = None):
    """Add docking result for a pose"""
    success = add_docking_result(
        job_uuid, pose_id, ligand_name,
        vina_score, gnina_score, rf_score, pdb_data
    )
    if success:
        return {"status": "added"}
    else:
        raise HTTPException(status_code=500, detail="Failed to add result")


@app.get("/jobs/{job_uuid}/results")
def list_results(job_uuid: str):
    """Get all docking results for a job"""
    results = get_docking_results(job_uuid)
    return {"results": results}


@app.get("/jobs/{job_uuid}/interactions")
def list_interactions(job_uuid: str, pose_id: Optional[int] = None):
    """Get interactions for a job"""
    interactions = get_interactions(job_uuid, pose_id)
    return {"interactions": interactions}


@app.post("/jobs/{job_uuid}/interactions")
def add_job_interaction(job_uuid: str, pose_id: int, interaction_type: str,
                       atom_a: str, atom_b: str, distance: float):
    """Add interaction data"""
    success = add_interaction(job_uuid, pose_id, interaction_type, atom_a, atom_b, distance)
    if success:
        return {"status": "added"}
    else:
        raise HTTPException(status_code=500, detail="Failed to add interaction")


@app.delete("/jobs/{job_uuid}")
def remove_job(job_uuid: str):
    """Delete a job"""
    success = delete_job(job_uuid)
    if success:
        return {"status": "deleted"}
    else:
        raise HTTPException(status_code=500, detail="Failed to delete job")


@app.post("/upload")
async def upload_file(file: UploadFile = File(...)):
    """Upload a file"""
    file_path = os.path.join(STORAGE_DIR, file.filename)
    
    with open(file_path, "wb") as f:
        content = await file.read()
        f.write(content)
    
    return {"filename": file.filename, "path": file_path}


@app.get("/download/{filename}")
def download_file(filename: str):
    """Download a file"""
    file_path = os.path.join(STORAGE_DIR, filename)
    
    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="File not found")
    
    with open(file_path, "r") as f:
        content = f.read()
    
    return {"filename": filename, "content": content}


@app.get("/security/status")
def security_status():
    """Get current security status"""
    try:
        from security.monitor import SecurityMonitor
        monitor = SecurityMonitor()
        status = monitor.get_latest_status()
        
        if status:
            return status
        else:
            return {
                "last_scan_at": None,
                "overall_severity": "NOT_SCANNED",
                "is_secure": True,
                "total_issues": 0,
                "scan_results": {}
            }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/security/scan")
def run_security_scan():
    """Run full security scan"""
    try:
        from security.monitor import SecurityMonitor
        monitor = SecurityMonitor()
        summary = monitor.run_full_scan()
        return summary
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/security/reports")
def security_reports(limit: int = 10):
    """Get security reports"""
    try:
        from security.monitor import SecurityMonitor
        monitor = SecurityMonitor()
        reports = monitor.get_recent_reports(limit)
        return {"reports": reports}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/security/issues")
def security_issues(scan_type: Optional[str] = None):
    """Get detailed security issues"""
    try:
        from security.monitor import SecurityMonitor
        monitor = SecurityMonitor()
        issues = monitor.get_security_issues(scan_type)
        return {"issues": issues}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


class ChatRequest(BaseModel):
    message: str


@app.post("/chat")
def chat(req: ChatRequest):
    """Chat with AI assistant (Ollama or offline fallback)"""
    logger.info(f"Chat request received")
    try:
        from ai.llm_router import get_router
        router = get_router()
        result = router.chat(req.message)
        logger.info(f"Chat response: provider={result.get('provider')}, available={result.get('available')}")
        return result
    except Exception as e:
        logger.error(f"Chat failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/chat/status")
def chat_status():
    """Get chat provider status"""
    try:
        from ai.llm_router import get_router
        router = get_router()
        status = {
            "provider": router.provider,
            "ollama_available": router.detect_ollama(),
            "models": router.get_available_models() if router.provider == "ollama" else []
        }
        logger.debug(f"Chat status: {status}")
        return status
    except Exception as e:
        logger.error(f"Chat status failed: {e}")
        return {
            "provider": "offline",
            "ollama_available": False,
            "error": str(e)
        }


class DockingProgress:
    """Thread-safe docking progress tracker"""
    _jobs: dict = {}
    _lock = threading.Lock()
    
    @classmethod
    def start_job(cls, job_id: str, total: int = 100):
        with cls._lock:
            cls._jobs[job_id] = {"progress": 0, "total": total, "status": "running", "message": "Initializing..."}
    
    @classmethod
    def update_progress(cls, job_id: str, progress: int, message: str = ""):
        with cls._lock:
            if job_id in cls._jobs:
                cls._jobs[job_id]["progress"] = min(progress, 100)
                cls._jobs[job_id]["message"] = message
    
    @classmethod
    def set_status(cls, job_id: str, status: str, message: str = ""):
        with cls._lock:
            if job_id in cls._jobs:
                cls._jobs[job_id]["status"] = status
                cls._jobs[job_id]["message"] = message
    
    @classmethod
    def get_progress(cls, job_id: str):
        with cls._lock:
            return cls._jobs.get(job_id, {"progress": 0, "total": 100, "status": "unknown", "message": ""})
    
    @classmethod
    def clear_job(cls, job_id: str):
        with cls._lock:
            if job_id in cls._jobs:
                del cls._jobs[job_id]


@app.post("/dock/start")
def start_docking_job(job_id: str = Form(...), total_ligands: int = Form(10)):
    """Start a docking job and return job ID"""
    logger.info(f"Starting docking job: {job_id} with {total_ligands} ligands")
    DockingProgress.start_job(job_id, total_ligands)
    
    def simulate_docking():
        """Simulate docking progress for demonstration"""
        for i in range(0, 101, 5):
            progress_data = DockingProgress.get_progress(job_id)
            if progress_data["status"] == "cancelled":
                break
            
            DockingProgress.update_progress(job_id, i, f"Processing ligand {i*total_ligands//100}/{total_ligands}")
            time.sleep(0.3)
        
        DockingProgress.set_status(job_id, "completed", "Docking complete!")
        logger.info(f"Docking job completed: {job_id}")
    
    thread = threading.Thread(target=simulate_docking, daemon=True)
    thread.start()
    
    return {"job_id": job_id, "status": "started"}


@app.post("/dock/{job_id}/cancel")
def cancel_docking_job(job_id: str):
    """Cancel a running docking job"""
    DockingProgress.set_status(job_id, "cancelled", "Job cancelled by user")
    logger.info(f"Docking job cancelled: {job_id}")
    return {"job_id": job_id, "status": "cancelled"}


@app.get("/dock/{job_id}/stream")
def stream_docking_progress(job_id: str):
    """Stream docking progress using Server-Sent Events (SSE)"""
    def generate():
        last_progress = -1
        while True:
            progress_data = DockingProgress.get_progress(job_id)
            current_progress = progress_data["progress"]
            
            if current_progress != last_progress:
                last_progress = current_progress
                event_data = json.dumps({
                    "progress": current_progress,
                    "status": progress_data["status"],
                    "message": progress_data["message"]
                })
                yield f"data: {event_data}\n\n"
            
            if progress_data["status"] in ["completed", "cancelled", "failed"]:
                break
            
            time.sleep(0.2)
    
    return StreamingResponse(generate(), media_type="text/event-stream")


@app.get("/dock/{job_id}/status")
def get_docking_status(job_id: str):
    """Get current docking job status"""
    progress_data = DockingProgress.get_progress(job_id)
    return progress_data


@app.get("/gpu/status")
def get_gpu_status():
    """Get GPU status using nvidia-smi"""
    try:
        import subprocess
        result = subprocess.run(
            ['nvidia-smi', '--query-gpu=index,name,utilization.gpu,memory.used,memory.total,temperature.gpu',
             '--format=csv,noheader,nounits'],
            capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0:
            lines = result.stdout.strip().split('\n')
            gpus = []
            for line in lines:
                parts = [p.strip() for p in line.split(',')]
                if len(parts) >= 6:
                    gpus.append({
                        "index": int(parts[0]),
                        "name": parts[1],
                        "utilization": int(parts[2]),
                        "memory_used": int(parts[3]),
                        "memory_total": int(parts[4]),
                        "temperature": int(parts[5])
                    })
            return {"available": True, "gpus": gpus}
    except FileNotFoundError:
        pass
    except Exception as e:
        logger.warning(f"GPU detection error: {e}")
    
    return {"available": False, "gpus": [], "message": "No GPU detected"}


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
