# BioDockify Docking Studio AI

<p align="center">
  <img src="https://img.shields.io/badge/Version-2.3.5-blue.svg" alt="Version">
  <img src="https://img.shields.io/badge/Python-3.11-green.svg" alt="Python">
  <img src="https://img.shields.io/badge/License-MIT-purple.svg" alt="License">
</p>

> **AI-Powered Autonomous Drug Discovery Platform** — runs at `http://localhost:8000`

An intelligent molecular docking platform with Discovery Studio-inspired UI, AI-powered molecule optimization, and automated drug discovery workflows.

## Features

### 🧬 ChemDraw - Molecule Editor
- Draw and analyze molecules with real-time 2D/3D visualization
- SMILES input with structure validation
- 12 pre-loaded FDA-approved drugs

### 🤖 AI Optimization
- AI-powered molecular modification
- Bioisosteric replacement
- Halogen, OH, NH2 group addition
- Aromatic ring expansion
- Flexibility reduction

### 🔬 Drug-like Analysis
- Lipinski Rule of 5 compliance
- MW, LogP, HBD, HBA calculations
- TPSA (Topological Polar Surface Area)
- Rotatable bonds analysis

### 🧪 Molecular Docking (Smart Energy-Based Routing)
- **AutoDock Vina 1.2** — High-precision empirical scoring function
- **GNINA** — CNN-based deep learning scoring (CNN, CNN affinity, CNN pose)
- **RF (Random Forest) Score** — Machine learning based rescoring
- **Smart Routing:**
  - Energy ≤ -5.0 kcal/mol → **Vina only** → Returns log, docking, grid files
  - Energy > -5.0 kcal/mol → **GNINA + RF** → Returns all files (Vina + GNINA outputs)
- RDKit-only file preparation (no Meeko/OpenBabel)
- Flexible ligand & rigid receptor docking
- Grid box auto-detection
- Real-time job tracking & pose visualization
- **Downloadable output files:** Log, Docking PDBQT, Grid configuration

### 📊 Molecular Dynamics
- OpenMM simulation engine
- GPU acceleration (CUDA/OpenCL)
- Automatic platform detection

### 🖥️ 3D Viewer
- NGL Viewer integration
- Ball-and-stick visualization
- Interactive rotation & zoom

### 💬 AI Assistant
- Multi-provider support (OpenAI, Claude, Gemini, DeepSeek, etc.)
- Drug discovery insights
- Molecule analysis suggestions

## Quick Start

```bash
# Pull and run
docker pull tajo9128/biodockify-studio-ai:latest
docker run -p 8000:8000 tajo9128/biodockify-studio-ai:latest

# Or build locally
git clone https://github.com/tajo9128/BioDockify-Docking-Studio-AI.git
cd BioDockify-Docking-Studio-AI
docker build -f Dockerfile.single -t biodockify-studio-ai .
docker run -p 8000:8000 biodockify-studio-ai
```

Then open **http://localhost:8000** in your browser.

## AI Providers

Configure API keys in Settings:

| Provider | Models |
|----------|--------|
| OpenAI | GPT-4o, GPT-4o-mini |
| Claude | Claude-3-Opus, Claude-3-Sonnet |
| Gemini | Gemini Pro |
| Mistral | Mistral Large, Mistral Medium |
| DeepSeek | DeepSeek Chat |
| Qwen | Qwen Turbo, Qwen Plus |
| SiliconFlow | Qwen, Llama models |
| OpenRouter | Multiple providers |
| Ollama | Local models |

## Molecule Library

Pre-loaded with 12 FDA-approved drugs:

| Drug | SMILES | Use |
|------|--------|-----|
| Aspirin | CC(=O)Oc1ccccc1C(=O)O | Anti-inflammatory |
| Caffeine | Cn1cnc2c1c(=O)n(c(=O)n2C)C | Stimulant |
| Glucose | OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O | Energy metabolism |
| Ibuprofen | CC(C)Cc1ccc(cc1)C(C)C(=O)O | Analgesic |
| Morphine | CN1CCc2c(O)ccc(c2C1)C(O)=O | Analgesic |
| Metformin | CN(C)N=C(N)N | Diabetes |
| Warfarin | CC(=O)OC(Cc1c(O)c2ccccc2oc1=O)C(c1ccccc1)=O | Anticoagulant |
| Sildenafil | CCCC1=C2N(C(=O)N1CCC)CCCC2c3ccc(cc3)S(=O)(=O)N | PDE5 inhibitor |

## Drug-like Properties

| Property | Rule | Description |
|----------|------|-------------|
| MW | < 500 Da | Molecular weight |
| LogP | < 5 | Lipophilicity |
| HBD | ≤ 5 | Hydrogen bond donors |
| HBA | ≤ 10 | Hydrogen bond acceptors |
| TPSA | < 140 Å² | Topological polar surface area |
| Rotatable | ≤ 10 | Flexible bonds |

## Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    BioDockify Studio AI v2.3.5                           │
├─────────────────────────────────────────────────────────────────────────┤
│  Frontend (SPA - embedded HTML/JS)                                       │
│  ├── ChemDraw Panel (smiles-drawer + Ketcher)                          │
│  ├── 3D Viewer (NGL Viewer)                                            │
│  ├── Properties Panel (Lipinski Rule of 5)                             │
│  ├── AI Suggestions Panel (RDKit drug-likeness)                       │
│  ├── Molecular Optimization (mutation strategies)                       │
│  ├── Docking Panel (Vina config, receptor/ligand upload)               │
│  ├── MD Panel (OpenMM params, GPU info, trajectory view)               │
│  ├── Results Panel (scores, trajectories, analysis)                    │
│  └── Job Management (auto-refresh, queue, history)                     │
├─────────────────────────────────────────────────────────────────────────┤
│  Backend (FastAPI)                                                      │
│  ├── GET  /api/stats                       (job statistics)            │
│  ├── GET  /api/health                     (health check)               │
│  ├── GET  /api/md/gpu-info                (CUDA/CPU/OpenCL detect)    │
│  ├── GET  /api/ai/providers               (LLM providers list)         │
│  ├── POST /api/ai/test                     (test API key)             │
│  ├── POST /api/ai/chat                     (AI chat)                  │
│  ├── POST /api/docking/jobs               (create docking job)        │
│  ├── GET  /api/docking/jobs               (list all jobs)            │
│  ├── GET  /api/docking/jobs/{id}          (get job details)           │
│  ├── POST /api/docking/run                 (run AutoDock Vina)         │
│  ├── POST /api/docking/gnina              (run GNINA CNN docking)     │
│  ├── POST /api/docking/rescore             (RF/CNN rescoring)          │
│  ├── GET  /api/docking/results/{id}        (docking scores/poses)     │
│  ├── POST /api/chem/properties             (RDKit molecular calc)      │
│  ├── POST /api/chem/suggestions            (drug-likeness analysis)    │
│  ├── POST /api/chem/dock                   (prepare + create job)     │
│  ├── GET  /api/chem/3d/{id}               (PDB structure for viewer)   │
│  ├── POST /api/md/simulate                (run OpenMM MD)             │
│  ├── GET  /api/md/results/{id}            (MD trajectory data)       │
│  └── POST /api/chem/optimize              (AI molecular optimization)  │
├─────────────────────────────────────────────────────────────────────────┤
│  Container (python:3.11-slim)                                           │
│  ├── rdkit-pypi==2022.9.5     (properties, prep, analysis)            │
│  ├── autodock-vina            (empirical & vinardo scoring)          │
│  ├── gnina                    (CNN/RF deep learning docking)         │
│  ├── openmm                   (MD simulations, GPU-accelerated)     │
│  ├── fastapi + uvicorn         (web framework)                        │
│  ├── biopython                (protein structure handling)            │
│  ├── meeko                    (ligand preparation for Vina)           │
│  ├── Storage: /app/data/jobs/{id}/{pdb,pdbqt,xtc,dcd}               │
│  └── Port: 8000                                                      │
└─────────────────────────────────────────────────────────────────────────┘
```

## Docker Images

| Image | Description |
|-------|-------------|
| `tajo9128/biodockify-studio-ai:latest` | Latest release |
| `tajo9128/biodockify-studio-ai:v2.3.5` | Versioned release |
| `tajo9128/biodockify-studio-ai:full` | Full stack (Vina/GNINA/RF) |

## Development

```bash
# Local development (single container)
pip install -r requirements.txt
python app.py

# Frontend (if using separate frontend)
cd frontend
npm install
npm run dev
```

## License

MIT License - See LICENSE file for details.

---

<p align="center">
  <strong>Biodockify Studio AI</strong><br>
  AI-Powered Autonomous Drug Discovery Platform
</p>
