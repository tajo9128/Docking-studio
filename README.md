# BioDockify

Production-ready molecular docking software with Discovery Studio-inspired UI.

## Features

- **Molecular Docking**: AutoDock Vina integration with RDKit for ligand preparation
- **Molecular Dynamics**: OpenMM with automatic GPU detection (CUDA > OpenCL > CPU)
- **BioDockify AI**: Multi-provider AI assistant (Claude, Gemini, OpenRouter, Mistral, etc.)
- **Single Container Deployment**: Runs at localhost:8000

## Quick Start

```bash
# Build and run
docker build -t biodockify -f Dockerfile.single .
docker run -p 8000:8000 biodockify
```

Open http://localhost:8000 in your browser.

## Development

```bash
# Backend
cd backend
pip install -r requirements.txt
uvicorn main:app --reload

# Frontend
cd frontend
npm install
npm run dev
```

## Architecture

- **Backend**: FastAPI with async support
- **Frontend**: React with TypeScript
- **Docking**: AutoDock Vina via subprocess
- **MD**: OpenMM with platform auto-detection
- **AI**: Channel-based routing to multiple LLM providers
