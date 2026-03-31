# BioDockify Dockerfile
# Single container for molecular docking and MD simulations

FROM python:3.11-slim

LABEL maintainer="BioDockify"
LABEL description="Molecular Docking Studio with AI Assistant"

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PATH="/opt/conda/bin:$PATH"

WORKDIR /app

RUN apt-get update && apt-get install -y \
    wget \
    curl \
    gnupg2 \
    software-properties-common \
    && rm -rf /var/lib/apt/lists/*

RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/opt/conda/bin:$PATH"

RUN conda install -y -c conda-forge \
    rdkit \
    openmm \
    openbabel \
    ambertools \
    cudatoolkit=11.8 \
    -c defaults \
    && conda clean -a

RUN pip install --no-cache-dir \
    fastapi==0.109.2 \
    uvicorn[standard]==0.27.1 \
    python-multipart==0.0.9 \
    aiofiles==23.2.1 \
    pydantic==2.6.1 \
    pydantic-settings==2.1.0 \
    httpx==0.26.0 \
    pandas==2.2.1

RUN wget -q https://vina.scripps.edu/wp-content/uploads/sites/55/2022/12/autodock_vina_1_2_3_linux_x86_64.tar.gz && \
    tar -xzf autodock_vina_1_2_3_linux_x86_64.tar.gz -C /usr/local/bin && \
    rm autodock_vina_1_2_3_linux_x86_64.tar.gz

COPY backend/ /app/backend/
COPY frontend/dist/ /app/frontend/dist/

RUN mkdir -p /app/data/jobs

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/health || exit 1

CMD ["python", "-c", """
import uvicorn
from fastapi import FastAPI, UploadFile, File, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from pathlib import Path
import uuid
import json
import os
import shutil
from datetime import datetime

app = FastAPI(title='BioDockify', version='2.3.0')

app.add_middleware(
    CORSMiddleware,
    allow_origins=['*'],
    allow_credentials=True,
    allow_methods=['*'],
    allow_headers=['*'],
)

DATA_DIR = Path('/app/data')
JOBS_DIR = DATA_DIR / 'jobs'
JOBS_DIR.mkdir(parents=True, exist_ok=True)

frontend_dist = Path('/app/frontend/dist')

try:
    app.mount('/static', StaticFiles(directory=str(frontend_dist)), name='static')
except:
    pass

@app.get('/')
async def root():
    index_path = frontend_dist / 'index.html'
    if index_path.exists():
        return FileResponse(str(index_path))
    return {'message': 'BioDockify API', 'version': '2.3.0', 'frontend': 'not built'}

@app.get('/health')
async def health():
    return {'status': 'healthy', 'timestamp': datetime.utcnow().isoformat()}

@app.get('/api/stats')
async def get_stats():
    total = len(list(JOBS_DIR.glob('*.json')))
    completed = len([f for f in JOBS_DIR.glob('*.json') if 'completed' in f.read_text()])
    return {'total_jobs': total, 'completed_jobs': completed, 'active_jobs': total - completed}

@app.get('/api/md/gpu-info')
async def gpu_info():
    import subprocess
    try:
        result = subprocess.run(['nvidia-smi'], capture_output=True, timeout=5)
        if result.returncode == 0:
            return {'platform': 'CUDA', 'cuda_available': True, 'opencl_available': False, 'message': 'GPU detected'}
    except:
        pass
    try:
        result = subprocess.run(['which', 'clinfo'], capture_output=True)
        if result.returncode == 0:
            return {'platform': 'OpenCL', 'cuda_available': False, 'opencl_available': True, 'message': 'OpenCL available'}
    except:
        pass
    return {'platform': 'CPU', 'cuda_available': False, 'opencl_available': False, 'message': 'Running on CPU'}

@app.get('/api/ai/providers')
async def list_providers():
    return {'providers': [
        {'id': 'openai', 'name': 'OpenAI'},
        {'id': 'claude', 'name': 'Claude'},
        {'id': 'gemini', 'name': 'Gemini'},
        {'id': 'mistral', 'name': 'Mistral'},
        {'id': 'deepseek', 'name': 'DeepSeek'},
        {'id': 'qwen', 'name': 'Qwen'},
        {'id': 'siliconflow', 'name': 'SiliconFlow'},
        {'id': 'openrouter', 'name': 'OpenRouter'},
        {'id': 'ollama', 'name': 'Ollama'}
    ]}

@app.post('/api/docking/jobs')
async def create_docking_job():
    job_id = str(uuid.uuid4())[:8]
    job_info = {
        'job_id': job_id,
        'type': 'docking',
        'status': 'created',
        'created_at': datetime.utcnow().isoformat()
    }
    (JOBS_DIR / f'{job_id}.json').write_text(json.dumps(job_info, indent=2))
    return job_info

@app.get('/api/docking/jobs')
async def list_docking_jobs():
    jobs = []
    for f in JOBS_DIR.glob('*.json'):
        job = json.loads(f.read_text())
        if job.get('type') == 'docking':
            jobs.append(job)
    return {'jobs': sorted(jobs, key=lambda x: x.get('created_at', ''), reverse=True)}

@app.get('/api/docking/jobs/{job_id}')
async def get_docking_job(job_id: str):
    job_file = JOBS_DIR / f'{job_id}.json'
    if job_file.exists():
        return json.loads(job_file.read_text())
    return {'job_id': job_id, 'status': 'not_found'}

@app.post('/api/md/jobs')
async def create_md_job():
    job_id = str(uuid.uuid4())[:8]
    job_info = {
        'job_id': job_id,
        'type': 'md',
        'status': 'created',
        'created_at': datetime.utcnow().isoformat()
    }
    (JOBS_DIR / f'{job_id}.json').write_text(json.dumps(job_info, indent=2))
    return job_info

@app.post('/api/ai/chat')
async def chat(request: dict):
    provider = request.get('provider', 'demo')
    message = request.get('message', '')
    
    responses = {
        'openai': f'OpenAI GPT processing: {message[:50]}...',
        'claude': f'Claude analyzing: {message[:50]}...',
        'gemini': f'Gemini insights for: {message[:50]}...',
        'demo': f'BioDockify AI (Demo): You asked about {message[:30]}...'
    }
    
    return {'response': responses.get(provider, responses['demo']), 'provider': provider}

@app.put('/api/ai/config')
async def update_config(config: dict):
    provider = config.get('provider', 'unknown')
    return {'status': 'updated', 'provider': provider}

@app.get('/api/ai/config/{provider}')
async def get_config(provider: str):
    return {'provider': provider, 'configured': False}

uvicorn.run(app, host='0.0.0.0', port=8000)
"""]
