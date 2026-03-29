# 🧬 Docking Studio v2.0

**Local-first AI Drug Discovery Studio** — runs in your browser at `http://localhost:8000`.

Built with microservices architecture + Nanobot Brain AI agent. One-command startup via Docker Desktop.

---

## ⚡ Quick Start

```bash
# 1. Install Docker Desktop → https://www.docker.com/products/docker-desktop
# 2. Clone and start
git clone https://github.com/tajo9128/Docking-studio
cd Docking-studio
docker compose up

# 3. Open browser → http://localhost:8000
```

---

## 🧠 What is Nanobot Brain?

Nanobot is your AI drug discovery assistant. Unlike traditional software, Nanobot:

- **Understands natural language** — "Find best inhibitor for protein X"
- **Plans workflows** — Creates execution plans using available tools
- **Orchestrates pipelines** — Pharmacophore → Screen → Dock → Rank automatically
- **Learns from results** — Improves recommendations based on screening data

---

## 🧩 Plugin System (Skills/Tools)

Nanobot Brain has pluggable tools organized by category:

### 🧪 Cheminformatics
| Tool | Description |
|------|-------------|
| `smiles_to_3d` | Convert SMILES to 3D structure (RDKit) |
| `convert_format` | Convert between PDB, SDF, mol, pdbqt |
| `optimize_molecule` | MMFF force field optimization |

### 🔬 Docking
| Tool | Description |
|------|-------------|
| `dock_ligand` | AutoDock Vina molecular docking |
| `run_batch_docking` | Batch docking for compound libraries |

### 🎯 Pharmacophore
| Tool | Description |
|------|-------------|
| `generate_pharmacophore` | RDKit-based pharmacophore from receptor/ligand |
| `screen_library` | Virtual screening against pharmacophore model |

### 🌐 Integrations
| Tool | Description |
|------|-------------|
| `fetch_protein` | Fetch from PDB (e.g., "1ABC") |
| `fetch_compounds` | Fetch from PubChem by CID |
| `search_compounds` | Search PubChem by name/formula |
| `similarity_search` | Tanimoto similarity search |

### 📊 Analysis
| Tool | Description |
|------|-------------|
| `analyze_interactions` | H-bonds, hydrophobic, π-π stacking |
| `predict_binding` | ML-based binding affinity prediction |
| `compare_ligands` | Rank ligands by docking score |

---

## 🏗️ Architecture

```
Docker Desktop
│
├── gateway (nginx:3000)    → React UI + routing
├── api-backend (8000)      → Main API gateway
├── brain-service (8001)     → Nanobot AI planner
├── docking-service (8002)  → AutoDock Vina
├── rdkit-service (8003)    → Molecule processing
├── pharmacophore-service (8004) → Pharmacophore modeling
├── redis-worker            → Celery job processor
├── redis (6379)            → Job queue broker
├── postgres (5432)         → Metadata storage
└── ollama (11434)          → AI (optional)
```

**Brain = Planner Only** — Brain doesn't run docking directly. It creates plans and queues jobs via the API.

---

## 🚀 Autonomous Pipelines

Nanobot can run full drug discovery workflows:

```python
# Example: Virtual screening pipeline
1. fetch_protein("1ABC")           # Get target from PDB
2. generate_pharmacophore()         # Find key interaction features  
3. search_compounds("aspirin")     # Find candidate ligands
4. smiles_to_3d()                  # Convert to 3D
5. screen_library()                 # Filter by pharmacophore
6. dock_ligand()                   # Dock top candidates
7. analyze_interactions()           # Understand binding
8. rank_results()                  # Sort by score
```

---

## ✨ Features

| Feature | Description |
|---------|-------------|
| **3D Viewer** | NGL Viewer for protein-ligand complexes |
| **Vina Docking** | Physics-based AutoDock Vina |
| **GNINA CNN** | Deep learning scoring (optional GPU) |
| **Pharmacophore** | RDKit-based model generation & screening |
| **Nanobot AI** | Natural language drug discovery assistant |
| **Batch Processing** | Redis queue for large libraries |
| **Persistent Storage** | PostgreSQL + volumes |

---

## 🖥️ System Requirements

| | Minimum | Recommended |
|--|---------|-------------|
| RAM | 4 GB | 8+ GB |
| CPU | 2 cores | 4+ cores |
| GPU | Optional | NVIDIA GPU |
| Storage | 10 GB | 20 GB |

---

## 🛠️ Common Commands

```bash
# Start
docker compose up

# Start with GPU (NVIDIA)
docker compose up -d

# Stop
docker compose down

# View logs
docker compose logs -f

# Rebuild
docker compose up -d --build

# Clean start
docker compose down -v
docker compose up
```

---

## 🤖 Adding New Tools

Extend Nanobot by adding Python tools:

```python
# services/brain-service/tools/my_tool.py
from tools import BaseTool, ToolInput, ToolOutput

class MyTool(BaseTool):
    name = "my_tool"
    description = "Description of what it does"
    category = "custom"
    
    async def execute(self, input_data: ToolInput) -> ToolOutput:
        # Your logic here
        return ToolOutput(success=True, data={...})

# Register in tools/__init__.py
registry.register(MyTool())
```

---

## 📚 Documentation

- [QUICKSTART.md](QUICKSTART.md) — Step-by-step guide
- [docs/](docs/) — Additional documentation

---

## ⚠️ Disclaimer

For research and educational purposes only. Not validated for clinical use.

---

**Version:** v2.0.0 | Apache License 2.0
