#!/usr/bin/env python3
"""
BioDockify Studio AI - Localhost Launcher
Serves ALL backend API routes + React frontend from a single process.

All features included:
  Docking (Vina/GNINA/RF) • Composite Scoring • Flexibility • Constraints
  ChemDraw (cleanup/IUPAC/InChI/conformers/export) • ADMET • QSAR
  Pharmacophore (hypothesis/exclusion volumes) • Interaction Analysis
  Molecular Dynamics (equilibration/checkpoint/MM-GBSA)
  Benchmarking (PDBbind/RMSD/enrichment)
  AI Chat (multi-provider LLM) • CrewAI (validated tools/memory/meta-learning)
  Active Learning (Bayesian optimization) • NL-to-DAG Compiler
  Critique Agent • Knowledge Graph (Memgraph)
  Meta-Parameter Learning • Classroom/Assignment System

Usage:
  python app/launcher.py              # Auto-port, open browser
  python app/launcher.py --port 8000  # Specific port
  python app/launcher.py --no-browser # Don't open browser
  python app/launcher.py --debug      # Debug logging
"""
import sys
import os
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import socket
import threading
import time
import logging
import argparse
import webbrowser
from pathlib import Path

# ─── Path resolution ───────────────────────────────────────────────────────
if getattr(sys, 'frozen', False):
    BASE_DIR = Path(sys._MEIPASS)
    IS_FROZEN = True
else:
    BASE_DIR = Path(__file__).resolve().parent.parent
    IS_FROZEN = False

# ─── Logging ────────────────────────────────────────────────────────────────
log_level = logging.DEBUG if os.getenv("BIODOCKIFY_DEBUG") else logging.INFO
logging.basicConfig(
    level=log_level,
    format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
    datefmt="%H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger("biodockify")

# ─── Import backend with ALL routes ────────────────────────────────────────
backend_dir = BASE_DIR / "backend"
if str(backend_dir) not in sys.path:
    sys.path.insert(0, str(backend_dir))

from main import app  # noqa: E402 — FastAPI app with all 4295 lines of routes

# ─── Mount React frontend ──────────────────────────────────────────────────
from fastapi.staticfiles import StaticFiles  # noqa: E402
from fastapi.responses import FileResponse, HTMLResponse  # noqa: E402

FRONTEND_DIST = BASE_DIR / "frontend" / "dist"
FRONTEND_BUILD = BASE_DIR / "frontend" / "build"

frontend_dir = None
for d in [FRONTEND_DIST, FRONTEND_BUILD]:
    if d.exists() and (d / "index.html").exists():
        frontend_dir = d
        break

if frontend_dir:
    assets_dir = frontend_dir / "assets"
    if assets_dir.exists():
        app.mount("/assets", StaticFiles(directory=str(assets_dir)), name="assets")
    logger.info(f"✅ Frontend mounted from {frontend_dir}")
else:
    logger.warning("⚠️ Frontend not built. Run: cd frontend && npx vite build")

# ─── SPA catch-all (BEFORE any other catch-alls) ──────────────────────────
@app.get("/app/{path:path}")
@app.get("/studio/{path:path}")
async def spa_fallback(path: str):
    if frontend_dir:
        index = frontend_dir / "index.html"
        if index.exists():
            return FileResponse(str(index))
    return HTMLResponse(
        "<html><body style='background:#1a1a2e;color:#fff;padding:40px;font-family:sans-serif'>"
        "<h1>🧬 BioDockify Studio AI v4.0.0</h1>"
        "<p>API running. Frontend not built — run <code>cd frontend && npx vite build</code></p>"
        f"<p>Dir: {BASE_DIR}</p>"
        "<h3>API Endpoints:</h3><ul>"
        "<li><a href='/api/health'>/api/health</a></li>"
        "<li><a href='/docs'>/docs (Swagger)</a></li>"
        "<li><a href='/redoc'>/redoc</a></li>"
        "</ul></body></html>"
    )

# ─── Health check ──────────────────────────────────────────────────────────
@app.get("/api/health")
def health():
    return {
        "status": "ok",
        "version": "4.0.0",
        "mode": "desktop" if IS_FROZEN else "localhost",
        "frontend": "bundled" if frontend_dir else "not_built",
        "features": [
            "docking", "composite_scoring", "flexibility", "constraints",
            "chemdraw", "admet", "qsar", "pharmacophore", "interactions",
            "molecular_dynamics", "benchmarking", "ai_chat", "crewai",
            "active_learning", "nl_to_dag", "critique_agent", "knowledge_graph",
            "meta_param_learning",
        ],
    }

# ─── Hardware detection ────────────────────────────────────────────────────
def check_hardware() -> dict:
    try:
        import psutil
        cpu = psutil.cpu_count(logical=False) or 2
        ram = psutil.virtual_memory().total / (1024 ** 3)
    except ImportError:
        cpu, ram = 2, 4.0
    return {"cpu_cores": cpu, "ram_gb": round(ram, 1)}

# ─── Port finder ───────────────────────────────────────────────────────────
def find_free_port(start: int = 8000, max_attempts: int = 20) -> int:
    for port in range(start, start + max_attempts):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(("127.0.0.1", port)) != 0:
                return port
    raise RuntimeError(f"No free port in range {start}-{start + max_attempts}")

# ─── Browser opener ────────────────────────────────────────────────────────
def open_browser(port: int, delay: float = 2.5):
    time.sleep(delay)
    url = f"http://127.0.0.1:{port}"
    logger.info(f"🌐 Opening BioDockify in browser: {url}")
    try:
        webbrowser.open(url)
    except Exception as e:
        logger.warning(f"⚠️ Could not open browser: {e}")
        logger.info(f"👉 Please open {url} manually")

# ─── Main ──────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="BioDockify Studio AI Launcher")
    parser.add_argument("--port", type=int, default=None, help="Port (auto-detect if not set)")
    parser.add_argument("--no-browser", action="store_true", help="Don't auto-open browser")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    args = parser.parse_args()

    if args.debug:
        os.environ["BIODOCKIFY_DEBUG"] = "1"
        logging.getLogger().setLevel(logging.DEBUG)

    hw = check_hardware()
    logger.info(f"🚀 BioDockify Studio AI v4.0.0")
    logger.info(f"💻 Hardware: {hw['cpu_cores']} cores, {hw['ram_gb']}GB RAM")
    logger.info(f"📦 Mode: {'desktop' if IS_FROZEN else 'localhost'}")
    logger.info(f"🔧 Frontend: {'bundled' if frontend_dir else 'not_built'}")

    port = args.port or find_free_port(8000)
    os.environ["BIODOCKIFY_PORT"] = str(port)
    logger.info(f"🔗 Server: http://127.0.0.1:{port}")
    logger.info(f"📖 API Docs: http://127.0.0.1:{port}/docs")

    if not args.no_browser:
        threading.Thread(target=open_browser, args=(port,), daemon=True).start()

    logger.info("🎯 Press Ctrl+C to stop")
    try:
        uvicorn.run(app, host="127.0.0.1", port=port, log_level="info", access_log=False)
    except KeyboardInterrupt:
        logger.info("🛑 Shutting down")

if __name__ == "__main__":
    import uvicorn  # noqa: E402
    main()
