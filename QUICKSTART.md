# 🚀 Quick Start Guide

## ⚠️ Before You Begin: Install Docker (Required)

Docking Studio runs inside Docker containers. **Docker Desktop must be installed first.**

### Download Docker Desktop

| OS | Download Link | Size |
|----|--------------|------|
| Windows | [docker.com/desktop/windows](https://www.docker.com/products/docker-desktop) | ~500 MB |
| Mac (Apple Silicon) | [docker.com/desktop/mac/apple-silicon](https://www.docker.com/products/docker-desktop) | ~500 MB |
| Mac (Intel) | [docker.com/desktop/mac/intel](https://www.docker.com/products/docker-desktop) | ~500 MB |
| Ubuntu | [docker.com/desktop/linux/ubuntu](https://www.docker.com/products/docker-desktop) | ~300 MB |

> **Windows users:** After installing, **restart your computer** — don't just log out/in.

### Start Docker Desktop

After installation, find and open **Docker Desktop** from your Applications/Start Menu.

Wait until you see the whale icon with green lights:

```
✅ Docker Desktop - Running
```

### Verify Docker is Ready

Open a **new** terminal and run:

```bash
docker --version
```

You should see something like `Docker version 27.x.x`. If you get an error, Docker Desktop is not running yet.

---

## Step 1: Clone or Download This Project

**If you have Git:**
```bash
git clone https://github.com/tajo9128/Docking-studio.git
cd Docking-studio
```

**If you don't have Git:**
→ Click the green **"Code"** button on GitHub
→ Click **"Download ZIP"**
→ Extract the ZIP file
→ Open terminal/command prompt in the extracted folder

---

## Step 2: Start Docking Studio

### Windows: Open `start.bat`

Double-click `start.bat` in the project folder.

Or open **Command Prompt** in the project folder and run:
```bat
start.bat
```

### Mac / Linux: Run `start.sh`

Open **Terminal** and run:
```bash
chmod +x start.sh
./start.sh
```

---

## Step 3: Wait and Open Browser

The script will show progress as it builds and starts the application.

**First time (takes 5-10 minutes):**
```
Building... (downloading Vina, GNINA, RDKit for the first time)
```

**After first run (takes ~10 seconds):**
```
Starting...
Waiting for backend...
✅ Docking Studio is ready!
```

Once ready, **open your browser** and go to:

### 🌐 [http://localhost:8000](http://localhost:8000)

---

## Quick Reference

| Task | Command |
|------|---------|
| Start (from project folder) | `start.bat` or `./start.sh` |
| Stop everything | `docker compose down` |
| View live logs | `docker compose logs -f backend` |
| Restart completely fresh | `docker compose down -v && start.bat` |
| Update to latest version | `git pull && docker compose up -d --build` |

---

## Port Reference

| URL | What it is |
|----|-----------|
| [http://localhost:8000](http://localhost:8000) | **Main Web UI** — your docking workspace |
| [http://localhost:8000/docs](http://localhost:8000/docs) | API Documentation — for developers |

---

## Troubleshooting

### "Docker is not running" / "command not found: docker"

→ **Docker Desktop is not started.**
→ Find "Docker Desktop" in your Applications/Start Menu and open it.
→ Wait 30 seconds for it to fully start.
→ Then run `start.bat` or `./start.sh` again.

---

### "Port 8000 is already in use"

Something else is using port 8000. Stop it with:
```bash
docker compose down
```
Then try starting again.

---

### "Backend failed to start" / containers keep restarting

Check what's wrong:
```bash
docker compose logs backend
```

If it mentions memory errors, open `docker-compose.yml` and lower these numbers:
```yaml
memory: 4G   # ← change to 2G if you have less RAM
```

---

### "First build failed" / "npm/node error"

You need **Node.js** installed for the frontend build to work inside Docker. This is usually pre-installed on Mac/Linux. On Windows, make sure you installed Docker Desktop (which includes Node.js in recent versions).

If you continue to have issues:
```bash
docker compose down -v
docker system prune -f
docker compose up -d --build
```

---

### The browser shows "This site can't be reached"

Wait 60 seconds. The backend takes time to start on first run.

If still not working:
```bash
docker compose logs backend --tail=20
```

---

## ⚠️ Common Student Problems

| Problem | Solution |
|---------|----------|
| "I don't have Docker" | Install Docker Desktop (free) from docker.com |
| "Docker installed but won't start" | Restart your computer. Enable WSL2 if prompted on Windows. |
| "I have limited RAM (4GB)" | Edit `docker-compose.yml` and change `memory: 4G` to `memory: 2G` |
| "I'm on a school computer" | Ask IT to enable Docker / Virtualization. |
| "I can't install software" | Use a cloud platform like Google Colab or GitHub Codespaces instead. |
| "Build keeps failing" | Run `docker system prune -f` to clear Docker cache, then try again. |

---

## What's Inside?

Once running, explore these pages:

| Page | What it does |
|------|-------------|
| **Dashboard** | System health, GPU status, quick actions |
| **New Docking** | Upload receptor + ligand → configure grid → run docking |
| **Job Queue** | See all running and completed jobs |
| **Results** | Binding scores, interactions, pose analysis |
| **3D Viewer** | Rotate, zoom, and screenshot molecular structures |
| **RMSD Analysis** | Compare binding poses |
| **Interactions** | H-bonds, hydrophobic, π-π stacking, salt bridges |
| **AI Assistant** | Ask questions about your results (needs Ollama) |
| **Security** | Container and dependency vulnerability scanner |
| **Settings** | Configure AI provider, Docker settings |

---

**Happy Docking! 🧬**
