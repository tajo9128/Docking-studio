# 🧬 Docking Studio

**A web application — runs in your browser at `http://localhost:8000`.**

No desktop app. No installation needed. Just open Docker, run the start script, and dock in your browser.

Built with AutoDock Vina + GNINA deep learning for accurate binding affinity prediction.

---

## ⚡ Quick Start

### Windows
```bat
1. Install Docker Desktop → https://www.docker.com/products/docker-desktop
2. Open start.bat
3. Open browser → http://localhost:8000
```

### Mac / Linux
```bash
1. Install Docker Desktop → https://www.docker.com/products/docker-desktop
2. Run: chmod +x start.sh && ./start.sh
3. Open browser → http://localhost:8000
```

---

## ✨ Features

| Feature | Description |
|---------|-------------|
| **3D Viewer** | Visualize protein-ligand complexes with 3Dmol.js |
| **Vina Docking** | Physics-based AutoDock Vina scoring |
| **GNINA CNN** | Deep learning CNN scoring for better accuracy |
| **RF-Score** | Random Forest consensus scoring |
| **Tri-Score Protocol** | Combines Vina + GNINA + RF-Score |
| **Pose Analysis** | RMSD, binding interactions, pharmacophores |
| **AI Assistant** | Optional Ollama-powered chat (auto-detected) |
| **Security Scanner** | Trivy, Bandit, Safety dependency checks |

---

## 🖥️ System Requirements

| | Minimum | Recommended |
|--|---------|-------------|
| RAM | 4 GB | 8+ GB |
| CPU | 2 cores | 4+ cores |
| GPU | Optional | NVIDIA GPU (for GNINA) |
| Storage | 5 GB | 10 GB |
| OS | Windows 10+, macOS 11+, Ubuntu 20.04+ | |

**Docker Desktop is required.** Get it free at [docker.com](https://www.docker.com/products/docker-desktop)

---

## 🚀 How It Works

```
┌─────────────────────────────────────────────────────────────┐
│  Your Browser (http://localhost:8000)                       │
│                                                             │
│  React SPA ── SSE/Rest ──→ FastAPI Backend (Docker)        │
│                                ↓                             │
│                          Vina + GNINA + RDKit               │
│                                ↓                             │
│                          SQLite Job Storage                 │
└─────────────────────────────────────────────────────────────┘

Optional: Ollama (port 11434) ── AI Chat Assistant
```

**Ollama is optional.** The docking studio works 100% offline without it.

---

## 📁 Project Structure

```
Docking-studio/
├── backend/            # FastAPI API (runs in Docker container)
│   ├── main.py         # API endpoints
│   ├── analysis.py     # RMSD, interactions, binding sites
│   └── db.py           # SQLite job storage
├── frontend/           # React SPA source (TypeScript)
├── docker-compose.yml  # Container orchestration
├── Dockerfile          # Multi-stage build (Node builds frontend, Python serves it)
├── start.sh            # Easy start script (Mac/Linux)
├── start.bat           # Easy start script (Windows)
└── QUICKSTART.md       # Step-by-step guide for students
```

**No desktop app.** Everything runs in the browser at `http://localhost:8000`.

---

## 🛠️ Common Commands

```bash
# Start (from project directory)
./start.sh          # Mac/Linux
start.bat           # Windows

# Stop
docker compose down

# View logs
docker compose logs -f backend

# Restart fresh
docker compose down -v
./start.sh

# Update to latest version
git pull origin main
docker compose up -d --build
```

---

## 🔧 Troubleshooting

### "Docker is not running"
Start Docker Desktop and wait for the whale icon to say "running".

### "Port 8000 is already in use"
```bash
docker compose down
# or
docker stop $(docker ps -q --filter "publish=8000")
```

### "Backend won't start"
```bash
docker compose logs backend
```
If it mentions memory, lower limits in `docker-compose.yml`.

### First build is slow
Normal. Docker downloads Vina, GNINA, RDKit (~3-5 GB). Subsequent builds use cache.

---

## 📚 Documentation

- [QUICKSTART.md](QUICKSTART.md) — Step-by-step for students
- [SECURITY.md](SECURITY.md) — Security scanning guide
- [docs/troubleshooting.md](docs/troubleshooting.md) — Common issues

---

## 🔬 Scientific Features

**Tri-Score Protocol:**
- **Vina Score** — Physics-based grid scoring
- **GNINA CNN** — Convolutional neural network (3D binding pose)
- **RF-Score** — Machine learning protein-ligand scoring

**Analysis Tools:**
- H-bond, hydrophobic, π-π stacking, salt bridge detection
- RMSD calculation for pose comparison
- Binding site residue identification
- Interaction diagram generation

---

## 🤖 Optional: Enable AI

Install [Ollama](https://ollama.com) for AI-powered docking interpretation.

```bash
# Install Ollama, then:
ollama pull llama3
```

Docking Studio auto-detects Ollama. Works fully offline without it.

---

## ⚠️ Disclaimer

For research and educational purposes only. Not validated for clinical or medical decision-making.

---

**Version:** v1.3.4 | Apache License 2.0
