# Docking Studio

Docking Studio is a professional desktop molecular docking platform built with:

* PyQt6 scientific desktop UI
* Dockerized backend (FastAPI)
* AutoDock Vina + GNINA integration
* Optional Ollama AI assistant (auto-detected)
* Agent Zero orchestration engine
* Security monitoring layer
* Multi-user workspace
* Plugin architecture

---

# 🚀 Features

* Real-time docking progress streaming
* Multi-pose 3D visualization
* 2D interaction diagram panel
* Pose clustering (RMSD-based)
* Binding pocket heatmap
* GNINA CNN heatmap overlay
* MM-GBSA energy estimation panel
* AI-assisted docking interpretation (optional)
* GPU utilization monitor
* Security status monitoring
* Scientific report export (PDF)
* Plugin extension system

---

# 🧱 System Architecture

```
PyQt6 Desktop
      ↓
Docker Backend (FastAPI)
      ↓
Docking Engine (Vina + GNINA)
      ↓
SQLite Job Storage
```

Optional:

```
Ollama (auto-detected)
```

Docking Studio works fully without AI.

---

# 📦 System Requirements

## Minimum

* 8 GB RAM
* 4 CPU cores
* Docker Desktop installed
* Python 3.10+

## Recommended

* 16–32 GB RAM
* NVIDIA GPU (for GNINA CNN acceleration)
* 8+ CPU cores

---

# 🐳 Installation (Backend)

### 1️⃣ Install Docker

Download and install Docker Desktop:
[https://www.docker.com/products/docker-desktop](https://www.docker.com/products/docker-desktop)

Verify:

```
docker --version
docker compose version
```

---

### 2️⃣ Clone Repository

```
git clone https://github.com/tajo9128/Docking-studio.git
cd Docking-studio
```

---

### 3️⃣ Configure Environment

Copy example environment:

```
cp .env.example .env
```

Edit `.env` if needed:

```
AI_MODE=auto
ALLOW_AI=true
OLLAMA_URL=http://host.docker.internal:11434
OLLAMA_MODEL=llama3
```

---

### 4️⃣ Start Backend

```
docker compose up -d --build
```

Verify backend:

```
http://localhost:8000/docs
```

---

# 🖥 Frontend Installation (PyQt6)

### Install Python Dependencies

```
pip install -r requirements.txt
```

### Run Desktop Application

```
python -m src.biodockify_main
```

---

# 🧠 Optional: Enable Ollama AI

Install Ollama:

[https://ollama.com](https://ollama.com)

Start Ollama normally.

Docking Studio will automatically detect it.

No manual configuration required.

If Ollama is not running → software works in offline deterministic mode.

---

# 🔐 Security Monitoring

Docking Studio includes:

* Trivy container scanning
* Bandit Python static analysis
* Safety dependency scanning
* Gitleaks secret detection

Check status:

```
http://localhost:8000/security/status
```

---

# 📊 Running a Docking Job

1. Load receptor (PDB / PDBQT)
2. Load ligand (SDF / PDBQT / SMILES)
3. Define grid or auto-detect from bound ligand
4. Click **Run Docking**
5. Monitor progress in real time
6. Analyze poses in viewer
7. Export report

---

# 📄 Exporting Scientific Report

Report includes:

* Docking scores
* GNINA CNN score
* MM-GBSA estimation
* Interaction table
* Heatmap visualization
* Pose image snapshot
* Security status
* AI interpretation (if enabled)

---

# 🔌 Plugin System

Add new plugins in:

```
plugins/
```

Each plugin must implement:

```
class DockingPlugin:
    name = "Plugin Name"
    def run(self, job_data):
        pass
```

Plugins load automatically at startup.

---

# 👥 Multi-User Workspace

Each user can:

* Create projects
* Store docking jobs
* Reopen previous results
* Maintain independent datasets

SQLite database persists between sessions.

---

# 🛠 Development Mode

Rebuild backend clean:

```
docker compose down -v
docker compose build --no-cache
docker compose up -d
```

---

# 📜 License

This project is licensed under the Apache License 2.0.

See LICENSE file for details.

---

# ⚠ Disclaimer

This software is provided for research and educational purposes.
No warranty is provided.
Not intended for clinical or medical decision making.

---

# 📧 Security Reporting

Report vulnerabilities privately through GitHub Security Advisories.

---

# 📌 Version

Current version: v1.2.3 (see `VERSION` file)

---

# 🏁 Quick Start

```bash
# Clone and start
git clone https://github.com/tajo9128/Docking-studio.git
cd Docking-studio
docker compose up -d

# Run frontend
pip install -r requirements.txt
python -m src.biodockify_main
```

---

# 📚 Documentation

- [Installation Guide](docs/installation.md)
- [User Guide](docs/user_guide.md)
- [Troubleshooting](docs/troubleshooting.md)
- [FAQ](docs/faq.md)

---

# 🙏 Acknowledgments

* AutoDock Vina
* GNINA
* RDKit
* ODDT
* PyQt6
* FastAPI
* Ollama
