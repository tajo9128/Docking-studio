import uvicorn
from fastapi import FastAPI, UploadFile, File, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse
from pydantic import BaseModel
from pathlib import Path
import uuid
import json
from datetime import datetime

app = FastAPI(title='BioDockify', version='2.3.1')

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

HTML_CONTENT = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>BioDockify</title>
    <link rel="stylesheet" href="https://unpkg.com/smiles-drawer@2.2.1/dist/smiles-drawer.min.css">
    <script src="https://unpkg.com/smiles-drawer@2.2.1/dist/smiles-drawer.min.js"></script>
    <script src="https://unpkg.com/ketcher-core@2.6.2/dist/ketcher-core.min.js"></script>
    <script src="https://unpkg.com/ngl@2.1.0/dist/ngl.js"></script>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; background: #1a1a2e; color: #e8e8e8; min-height: 100vh; }
        .app { display: flex; flex-direction: column; min-height: 100vh; }
        .header { background: linear-gradient(135deg, #16213e 0%, #0f3460 100%); padding: 1rem 2rem; display: flex; align-items: center; justify-content: space-between; border-bottom: 1px solid #2a2a4a; }
        .header h1 { color: #00d9ff; font-size: 1.5rem; }
        .header-nav a { color: #a0a0a0; text-decoration: none; padding: 0.5rem 1rem; margin-left: 0.5rem; border-radius: 4px; }
        .header-nav a:hover, .header-nav a.active { color: #00d9ff; background: rgba(0,217,255,0.1); }
        .main { display: flex; flex: 1; }
        .editor-container { display: grid; grid-template-columns: 1fr 380px; gap: 1rem; height: calc(100vh - 200px); }
        .editor-main { background: #16213e; border-radius: 8px; border: 1px solid #2a2a4a; display: flex; flex-direction: column; overflow: hidden; }
        .editor-tabs { display: flex; background: #0f3460; border-bottom: 1px solid #2a2a4a; }
        .editor-tab { padding: 0.75rem 1.5rem; cursor: pointer; border: none; background: transparent; color: #a0a0a0; font-size: 0.875rem; }
        .editor-tab.active { color: #00d9ff; border-bottom: 2px solid #00d9ff; }
        .editor-toolbar { background: #0f3460; padding: 0.75rem; border-bottom: 1px solid #2a2a4a; display: flex; gap: 0.5rem; flex-wrap: wrap; }
        .editor-canvas { flex: 1; display: flex; align-items: center; justify-content: center; padding: 1rem; overflow: auto; }
        .viewer-3d { width: 100%; height: 100%; }
        .editor-sidebar { display: flex; flex-direction: column; gap: 0.75rem; overflow-y: auto; max-height: calc(100vh - 200px); }
        .smiles-input { width: 100%; padding: 0.75rem; background: #1a1a2e; border: 1px solid #2a2a4a; border-radius: 4px; color: #e8e8e8; font-family: monospace; resize: vertical; min-height: 60px; font-size: 0.875rem; }
        .preview-canvas { background: #1a1a2e; border-radius: 4px; padding: 0.5rem; min-height: 150px; display: flex; align-items: center; justify-content: center; }
        .preview-canvas canvas { max-width: 100%; }
        .prop-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 0.5rem; }
        .prop-item { background: #1a1a2e; padding: 0.5rem; border-radius: 4px; }
        .prop-label { font-size: 0.7rem; color: #a0a0a0; text-transform: uppercase; }
        .prop-value { font-size: 0.9rem; font-weight: 600; color: #00d9ff; }
        .suggestion-item { background: #1a1a2e; padding: 0.5rem 0.75rem; border-radius: 4px; margin-bottom: 0.5rem; font-size: 0.8rem; border-left: 3px solid #ffab00; }
        .suggestion-item.good { border-left-color: #00c853; }
        .suggestion-item.warning { border-left-color: #ffab00; }
        .suggestion-item.info { border-left-color: #00d9ff; }
        .tab-content { display: none; flex: 1; padding: 1rem; overflow: auto; }
        .tab-content.active { display: block; }
        .sidebar { width: 220px; background: #16213e; border-right: 1px solid #2a2a4a; padding: 1rem 0; }
        .sidebar a { color: #a0a0a0; text-decoration: none; padding: 0.75rem 1.5rem; display: flex; align-items: center; gap: 0.75rem; transition: all 0.2s; border-left: 3px solid transparent; }
        .sidebar a:hover { background: rgba(0,217,255,0.05); color: #e8e8e8; }
        .sidebar a.active { background: rgba(0,217,255,0.1); color: #00d9ff; border-left-color: #00d9ff; }
        .page { flex: 1; padding: 2rem; overflow-y: auto; }
        .page h2 { font-size: 1.75rem; margin-bottom: 0.5rem; }
        .page p { color: #a0a0a0; }
        .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1rem; margin: 2rem 0; }
        .stat-card { background: #16213e; border-radius: 8px; padding: 1.5rem; border: 1px solid #2a2a4a; }
        .stat-card h3 { font-size: 0.875rem; color: #a0a0a0; margin-bottom: 0.5rem; }
        .stat-card .value { font-size: 2rem; font-weight: 600; color: #00d9ff; }
        .card { background: #16213e; border-radius: 8px; padding: 1.5rem; border: 1px solid #2a2a4a; margin-bottom: 1rem; overflow: visible; }
        .card h3 { margin-bottom: 1rem; }
        .btn { padding: 0.75rem 1.5rem; border-radius: 6px; border: none; cursor: pointer; font-weight: 500; transition: all 0.2s; position: relative; z-index: 10; pointer-events: auto; }
        .btn-primary { background: #00d9ff; color: #1a1a2e; }
        .btn-primary:hover { background: #00b8d9; }
        .btn-secondary { background: #0f3460; color: #e8e8e8; border: 1px solid #2a2a4a; }
        .drop-zone { border: 2px dashed #2a2a4a; border-radius: 8px; padding: 3rem; text-align: center; cursor: pointer; transition: all 0.2s; }
        .drop-zone:hover { border-color: #00d9ff; background: rgba(0,217,255,0.05); }
        .form-group { margin-bottom: 1rem; }
        .form-group label { display: block; margin-bottom: 0.5rem; color: #a0a0a0; font-size: 0.875rem; }
        .form-group input, .form-group select { width: 100%; padding: 0.75rem; background: #1a1a2e; border: 1px solid #2a2a4a; border-radius: 4px; color: #e8e8e8; font-size: 1rem; }
        .grid-2 { display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; }
        .job-list { display: flex; flex-direction: column; gap: 0.5rem; }
        .job-item { background: #1a1a2e; border: 1px solid #2a2a4a; border-radius: 4px; padding: 1rem; display: flex; justify-content: space-between; cursor: pointer; }
        .job-item:hover { border-color: #00d9ff; }
        .job-item.active { border-color: #00d9ff; background: rgba(0,217,255,0.05); }
        .job-id { font-family: monospace; color: #00d9ff; }
        .job-status { font-size: 0.75rem; padding: 0.25rem 0.5rem; border-radius: 4px; }
        .status-dot { width: 8px; height: 8px; border-radius: 50%; display: inline-block; margin-right: 0.5rem; }
        .status-dot.green { background: #00c853; }
        .status-dot.yellow { background: #ffab00; }
        .status-bar { background: #16213e; padding: 0.5rem 2rem; border-top: 1px solid #2a2a4a; display: flex; justify-content: space-between; font-size: 0.75rem; color: #a0a0a0; }
        .ai-chat { display: flex; flex-direction: column; height: 400px; }
        .chat-messages { flex: 1; overflow-y: auto; padding: 1rem; background: #1a1a2e; border-radius: 4px; margin-bottom: 1rem; }
        .chat-message { margin-bottom: 1rem; padding: 0.75rem; border-radius: 4px; }
        .chat-message.user { background: #0f3460; margin-left: 2rem; }
        .chat-message.ai { background: #16213e; margin-right: 2rem; }
        .chat-input { display: flex; gap: 0.5rem; }
        .chat-input input { flex: 1; padding: 0.75rem; background: #1a1a2e; border: 1px solid #2a2a4a; border-radius: 4px; color: #e8e8e8; }
        .provider-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(250px, 1fr)); gap: 1rem; }
        .provider-card { background: #16213e; border: 1px solid #2a2a4a; border-radius: 8px; padding: 1.5rem; }
        .provider-card h4 { margin-bottom: 0.5rem; }
        @media (max-width: 768px) { .grid-2 { grid-template-columns: 1fr; } .sidebar { width: 60px; } .sidebar a span { display: none; } }
    </style>
</head>
<body>
    <div class="app">
        <header class="header">
            <h1>BioDockify</h1>
            <nav class="header-nav">
                <a href="#" onclick="showPage('dashboard')" class="active">Home</a>
                <a href="#" onclick="showPage('docking')">Docking</a>
                <a href="#" onclick="showPage('results')">Results</a>
                <a href="#" onclick="showPage('ai')">AI</a>
            </nav>
        </header>
        <div class="main">
            <aside class="sidebar">
                <a href="#" onclick="showPage('dashboard')" class="active"><span>🏠</span> Dashboard</a>
                <a href="#" onclick="showPage('docking')"><span>🧬</span> Docking</a>
                <a href="#" onclick="showPage('chemdraw')"><span>🖌️</span> ChemDraw</a>
                <a href="#" onclick="showPage('md')"><span>📊</span> MD Simulation</a>
                <a href="#" onclick="showPage('results')"><span>📈</span> Results</a>
                <a href="#" onclick="showPage('ai')"><span>🤖</span> AI Assistant</a>
                <a href="#" onclick="showPage('settings')"><span>⚙️</span> Settings</a>
            </aside>
            <main class="page" id="page-content"></main>
        </div>
        <footer class="status-bar">
            <span><span class="status-dot green"></span>BioDockify Backend</span>
            <span>v2.3.1</span>
        </footer>
    </div>
    <script>
        let currentPage = 'dashboard';
        let stats = { total_jobs: 0, completed_jobs: 0, active_jobs: 0 };
        let gpuInfo = { platform: 'Loading...' };
        let jobs = [];
        let messages = [];
        let selectedProvider = 'openai';
        
        const pages = {
            dashboard: () => `
                <h2>Dashboard</h2><p>BioDockify - Molecular Docking Studio</p>
                <div class="stats-grid">
                    <div class="stat-card"><h3>Total Jobs</h3><div class="value">${stats.total_jobs}</div></div>
                    <div class="stat-card"><h3>Completed</h3><div class="value">${stats.completed_jobs}</div></div>
                    <div class="stat-card"><h3>Active</h3><div class="value">${stats.active_jobs}</div></div>
                    <div class="stat-card"><h3>GPU Platform</h3><div class="value" style="font-size:1rem">${gpuInfo.platform}</div></div>
                </div>
                <div class="grid-2">
                    <div class="card"><h3>Quick Actions</h3>
                        <div style="display:flex;flex-direction:column;gap:0.5rem;margin-top:1rem">
                            <button class="btn btn-primary" onclick="showPage('docking')">New Docking Job</button>
                            <button class="btn btn-secondary" onclick="showPage('md')">Start MD Simulation</button>
                            <button class="btn btn-secondary" onclick="showPage('ai')">Ask BioDockify AI</button>
                        </div>
                    </div>
                    <div class="card"><h3>System Status</h3>
                        <div style="margin-top:1rem">
                            <div style="margin-bottom:0.5rem"><span class="status-dot green"></span>Backend API</div>
                            <div><span class="status-dot ${gpuInfo.platform !== 'CPU' ? 'green' : 'yellow'}"></span>GPU: ${gpuInfo.platform}</div>
                        </div>
                    </div>
                </div>`,
            docking: () => `
                <h2>Molecular Docking</h2><p>Run AutoDock Vina simulations</p>
                <div class="grid-2">
                    <div>
                        <div class="card"><h3>Receptor</h3>
                            <div class="drop-zone" onclick="document.getElementById('receptor-input').click()"><p>Drop PDB/PDBQT file</p></div>
                            <input type="file" id="receptor-input" accept=".pdb,.pdbqt" style="display:none">
                        </div>
                        <div class="card"><h3>Ligand</h3>
                            <div class="drop-zone" onclick="document.getElementById('ligand-input').click()"><p>Drop PDB/SDF file</p></div>
                            <input type="file" id="ligand-input" accept=".pdb,.pdbqt,.sdf" style="display:none">
                            <div style="text-align:center;color:#a0a0a0;margin:1rem 0">- OR -</div>
                            <div class="form-group"><label>SMILES String</label>
                                <input type="text" id="smiles-input" placeholder="CC(=O)Oc1ccccc1C(=O)O">
                            </div>
                        </div>
                        <div class="card"><h3>Parameters</h3>
                            <div class="form-group"><label>Exhaustiveness: <span id="exhaust-value">32</span></label>
                                <input type="range" min="1" max="64" value="32" oninput="document.getElementById('exhaust-value').textContent=this.value">
                            </div>
                            <div class="form-group"><label>Number of Poses: <span id="poses-value">10</span></label>
                                <input type="range" min="1" max="20" value="10" oninput="document.getElementById('poses-value').textContent=this.value">
                            </div>
                            <button id="start-docking-btn" class="btn btn-primary" style="width:100%;margin-top:1rem;position:relative;z-index:10">Start Docking</button>
                        </div>
                    </div>
                    <div class="card"><h3>Job Queue</h3>
                        <div class="job-list" id="job-list">
                            ${jobs.length === 0 ? '<p style="color:#a0a0a0">No jobs yet</p>' : ''}
                            ${jobs.map(j => '<div class="job-item"><span class="job-id">Job ' + j.job_id + '</span><span class="job-status">' + j.status + '</span></div>').join('')}
                        </div>
                    </div>
                </div>`,
            chemdraw: () => `
                <h2>ChemDraw - Structure Editor</h2><p>Draw molecules, analyze properties, and prepare for docking</p>
                <div class="editor-container">
                    <div class="editor-main">
                        <div class="editor-tabs">
                            <button class="editor-tab active" onclick="switchTab('draw')">2D Draw</button>
                            <button class="editor-tab" onclick="switchTab('viewer')">2D Viewer</button>
                            <button class="editor-tab" onclick="switchTab('viewer3d')">3D Viewer</button>
                        </div>
                        <div class="editor-toolbar">
                            <button class="btn btn-secondary" onclick="clearEditor()">Clear</button>
                            <button class="btn btn-secondary" onclick="copySmiles()">Copy SMILES</button>
                            <button class="btn btn-primary" onclick="dockMolecule()">🚀 Dock Molecule</button>
                            <button class="btn btn-secondary" onclick="analyzeMolecule()">🔬 Analyze</button>
                        </div>
                        <div id="tab-draw" class="tab-content active">
                            <div style="padding:1rem;text-align:center;color:#a0a0a0;">
                                <p style="font-size:2rem;">🖌️</p>
                                <p>SMILES-based structure editor</p>
                                <p style="font-size:0.875rem;margin-top:0.5rem">Enter SMILES in the sidebar to visualize</p>
                            </div>
                        </div>
                        <div id="tab-viewer" class="tab-content">
                            <div class="preview-canvas" id="preview-canvas" style="min-height:400px;width:100%;">
                                <canvas id="preview-smiles-canvas"></canvas>
                            </div>
                        </div>
                        <div id="tab-viewer3d" class="tab-content">
                            <div id="ngl-container" style="width:100%;height:400px;background:#1a1a2e;border-radius:4px;"></div>
                        </div>
                    </div>
                    <div class="editor-sidebar">
                        <div class="card">
                            <h3>SMILES Input</h3>
                            <textarea class="smiles-input" id="smiles-input" placeholder="Enter SMILES (e.g., CC(=O)Oc1ccccc1C(=O)O)">CC(=O)Oc1ccccc1C(=O)O</textarea>
                            <button class="btn btn-primary" style="width:100%;margin-top:0.5rem" onclick="analyzeMolecule()">🔬 Analyze</button>
                        </div>
                        <div class="card">
                            <h3>Properties (Lipinski)</h3>
                            <div class="prop-grid" id="prop-grid">
                                <div class="prop-item"><div class="prop-label">MW</div><div class="prop-value" id="prop-mw">-</div></div>
                                <div class="prop-item"><div class="prop-label">LogP</div><div class="prop-value" id="prop-logp">-</div></div>
                                <div class="prop-item"><div class="prop-label">HBD</div><div class="prop-value" id="prop-hbd">-</div></div>
                                <div class="prop-item"><div class="prop-label">HBA</div><div class="prop-value" id="prop-hba">-</div></div>
                                <div class="prop-item"><div class="prop-label">TPSA</div><div class="prop-value" id="prop-tpsa">-</div></div>
                                <div class="prop-item"><div class="prop-label">Rot. Bonds</div><div class="prop-value" id="prop-rot">-</div></div>
                                <div class="prop-item"><div class="prop-label">Formula</div><div class="prop-value" id="prop-formula" style="font-size:0.75rem">-</div></div>
                            </div>
                        </div>
                        <div class="card">
                            <h3>💡 AI Suggestions</h3>
                            <div id="ai-suggestions">
                                <p style="color:#a0a0a0;font-size:0.8rem">Click "Analyze" to get suggestions</p>
                            </div>
                        </div>
                        <div class="card">
                            <h3>Quick Examples</h3>
                            <div style="display:grid;grid-template-columns:1fr 1fr;gap:0.5rem;">
                                <button class="btn btn-secondary" style="font-size:0.7rem" onclick="loadExample('Aspirin','CC(=O)Oc1ccccc1C(=O)O')">Aspirin</button>
                                <button class="btn btn-secondary" style="font-size:0.7rem" onclick="loadExample('Caffeine','Cn1cnc2c1c(=O)n(c(=O)n2C)C')">Caffeine</button>
                                <button class="btn btn-secondary" style="font-size:0.7rem" onclick="loadExample('Glucose','OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O')">Glucose</button>
                                <button class="btn btn-secondary" style="font-size:0.7rem" onclick="loadExample('Benzene','c1ccccc1')">Benzene</button>
                                <button class="btn btn-secondary" style="font-size:0.7rem" onclick="loadExample('Ibuprofen','CC(C)Cc1ccc(cc1)C(C)C(=O)O')">Ibuprofen</button>
                                <button class="btn btn-secondary" style="font-size:0.7rem" onclick="loadExample('Morphine','CN1CCc2c(O)ccc(c2C1)C(O)=O')">Morphine</button>
                            </div>
                        </div>
                        <div class="card">
                            <h3>Docking Status</h3>
                            <div id="dock-status" style="font-size:0.8rem;color:#a0a0a0;">
                                <p>No active docking jobs</p>
                            </div>
                        </div>
                    </div>
                </div>`,
            md: () => `
                <h2>Molecular Dynamics</h2><p>Run OpenMM MD simulations with GPU acceleration</p>
                <div class="grid-2">
                    <div>
                        <div class="card"><h3>Structure File</h3>
                            <div class="drop-zone" onclick="document.getElementById('pdb-input').click()"><p>Drop PDB file</p></div>
                            <input type="file" id="pdb-input" accept=".pdb" style="display:none">
                        </div>
                        <div class="card"><h3>Parameters</h3>
                            <div class="form-group"><label>Force Field</label><select><option>AMBER14 All</option><option>AMBER14 TIP3P</option><option>CHARMM36</option></select></div>
                            <div class="form-group"><label>Duration (ns): <span id="duration-value">10</span></label>
                                <input type="range" min="1" max="100" value="10" oninput="document.getElementById('duration-value').textContent=this.value">
                            </div>
                            <div class="form-group"><label>Temperature (K): <span id="temp-value">300</span></label>
                                <input type="range" min="200" max="400" value="300" oninput="document.getElementById('temp-value').textContent=this.value">
                            </div>
                            <button class="btn btn-primary" style="width:100%">Start MD Simulation</button>
                        </div>
                    </div>
                    <div>
                        <div class="card"><h3>Active Simulations</h3><p style="color:#a0a0a0">No active simulations</p></div>
                        <div class="card"><h3>GPU Acceleration</h3>
                            <p style="color:#a0a0a0">OpenMM will detect: CUDA (NVIDIA) > OpenCL (AMD/Intel) > CPU</p>
                        </div>
                    </div>
                </div>`,
            results: () => `
                <h2>Results</h2><p>View docking and simulation results</p>
                <div class="grid-2">
                    <div class="card"><h3>All Jobs</h3>
                        <div class="job-list">
                            ${jobs.length === 0 ? '<p style="color:#a0a0a0">No jobs found</p>' : ''}
                            ${jobs.map(j => '<div class="job-item"><div><span class="job-id">' + j.job_id + '</span><div style="font-size:0.75rem;color:#a0a0a0">' + j.type + '</div></div><span class="job-status">' + j.status + '</span></div>').join('')}
                        </div>
                    </div>
                    <div class="card"><h3>Job Details</h3><p style="color:#a0a0a0">Select a job to view details</p></div>
                </div>`,
            ai: () => `
                <h2>BioDockify AI</h2><p>Ask about molecular docking and simulations</p>
                <div class="card">
                    <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:1rem">
                        <h3>Chat</h3>
                        <select value="openai" onchange="selectedProvider=this.value" style="padding:0.5rem;background:#1a1a2e;color:#e8e8e8;border:1px solid #2a2a4a;border-radius:4px">
                            <option value="openai">OpenAI</option><option value="claude">Claude</option><option value="gemini">Gemini</option>
                            <option value="mistral">Mistral</option><option value="deepseek">DeepSeek</option><option value="qwen">Qwen</option>
                            <option value="siliconflow">SiliconFlow</option><option value="openrouter">OpenRouter</option><option value="ollama">Ollama</option>
                        </select>
                    </div>
                    <div class="ai-chat">
                        <div class="chat-messages" id="chat-messages">
                            ${messages.length === 0 ? '<p style="color:#a0a0a0;text-align:center;margin-top:2rem">Ask me anything!</p>' : ''}
                            ${messages.map(m => '<div class="chat-message ' + m.role + '"><strong>' + (m.role === 'user' ? 'You' : 'BioDockify AI') + ':</strong> ' + m.content + '</div>').join('')}
                        </div>
                        <div class="chat-input">
                            <input type="text" id="chat-input" placeholder="Ask BioDockify AI..." onkeypress="if(event.key==='Enter')sendMessage()">
                            <button class="btn btn-primary" onclick="sendMessage()">Send</button>
                        </div>
                    </div>
                </div>`,
            settings: () => `
                <h2>Settings</h2><p>Configure AI providers</p>
                <div class="card"><h3>API Keys</h3>
                    <p style="color:#a0a0a0;margin:1rem 0">Configure your API keys. Keys are stored locally.</p>
                    <div class="provider-grid">
                        <div class="provider-card"><h4>DeepSeek</h4><p id="deepseek-status" style="color:#ffab00">Not configured</p>
                            <div class="form-group" style="margin-top:0.5rem"><input type="password" id="deepseek-key" placeholder="sk-..." style="width:100%;padding:0.5rem;background:#1a1a2e;border:1px solid #2a2a4a;border-radius:4px;color:#e8e8e8"></div>
                            <button class="btn btn-secondary" style="width:48%;margin-top:0.5rem" onclick="testProvider('deepseek')">Test</button>
                            <button class="btn btn-primary" style="width:48%;margin-top:0.5rem;float:right" onclick="saveProvider('deepseek')">Save</button>
                        </div>
                        <div class="provider-card"><h4>OpenAI</h4><p id="openai-status" style="color:#ffab00">Not configured</p>
                            <div class="form-group" style="margin-top:0.5rem"><input type="password" id="openai-key" placeholder="sk-..." style="width:100%;padding:0.5rem;background:#1a1a2e;border:1px solid #2a2a4a;border-radius:4px;color:#e8e8e8"></div>
                            <button class="btn btn-secondary" style="width:48%;margin-top:0.5rem" onclick="testProvider('openai')">Test</button>
                            <button class="btn btn-primary" style="width:48%;margin-top:0.5rem;float:right" onclick="saveProvider('openai')">Save</button>
                        </div>
                        <div class="provider-card"><h4>Claude</h4><p id="claude-status" style="color:#ffab00">Not configured</p>
                            <div class="form-group" style="margin-top:0.5rem"><input type="password" id="claude-key" placeholder="sk-ant-..." style="width:100%;padding:0.5rem;background:#1a1a2e;border:1px solid #2a2a4a;border-radius:4px;color:#e8e8e8"></div>
                            <button class="btn btn-secondary" style="width:48%;margin-top:0.5rem" onclick="testProvider('claude')">Test</button>
                            <button class="btn btn-primary" style="width:48%;margin-top:0.5rem;float:right" onclick="saveProvider('claude')">Save</button>
                        </div>
                    </div>
                </div>
                <div class="card"><h3>About</h3><p>BioDockify v2.3.1</p><p style="color:#a0a0a0;margin-top:0.5rem">Single-container molecular docking studio</p></div>`
        };
        
        function showPage(page) {
            currentPage = page;
            document.getElementById('page-content').innerHTML = pages[page]();
            document.querySelectorAll('.sidebar a, .header-nav a').forEach(a => a.classList.remove('active'));
            document.querySelectorAll('.sidebar a, .header-nav a').forEach(a => {
                if (a.textContent.toLowerCase().includes(page === 'dashboard' ? 'home' : page === 'md' ? 'simulation' : page)) a.classList.add('active');
            });
        }
        
        async function loadStats() { try { const res = await fetch('/api/stats'); stats = await res.json(); if (currentPage === 'dashboard') showPage('dashboard'); } catch(e) {} }
        async function loadGPUInfo() { try { const res = await fetch('/api/md/gpu-info'); gpuInfo = await res.json(); if (currentPage === 'dashboard') showPage('dashboard'); } catch(e) {} }
        async function loadJobs() { try { const res = await fetch('/api/docking/jobs'); const data = await res.json(); jobs = data.jobs || []; if (currentPage === 'docking' || currentPage === 'results') showPage(currentPage); } catch(e) {} }
        
        async function submitDockingJob() {
            console.log('Starting docking job...');
            try { 
                const res = await fetch('/api/docking/jobs', { method: 'POST' }); 
                const job = await res.json(); 
                console.log('Job created:', job);
                jobs.unshift(job); 
                showPage('docking'); 
            } catch(e) { 
                console.error('Error:', e); 
            }
        }
        
        document.addEventListener('DOMContentLoaded', function() {
            document.addEventListener('click', function(e) {
                if (e.target && (e.target.id === 'start-docking-btn' || e.target.textContent === 'Start Docking')) {
                    submitDockingJob();
                }
            });
        });
        
        async function sendMessage() {
            const input = document.getElementById('chat-input'); const text = input.value.trim();
            if (!text) return;
            messages.push({ role: 'user', content: text }); input.value = ''; showPage('ai');
            try {
                const res = await fetch('/api/ai/chat', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ message: text, provider: selectedProvider }) });
                const data = await res.json(); messages.push({ role: 'ai', content: data.response });
            } catch(e) { messages.push({ role: 'ai', content: 'Error: Could not connect' }); }
            showPage('ai');
        }
        
        function switchTab(tabName) {
            document.querySelectorAll('.editor-tab').forEach(t => t.classList.remove('active'));
            document.querySelectorAll('.tab-content').forEach(t => t.classList.remove('active'));
            document.querySelector(`.editor-tab[onclick="switchTab('${tabName}')"]`).classList.add('active');
            document.getElementById('tab-' + tabName).classList.add('active');
            if (tabName === 'viewer3d') init3DViewer();
        }
        
        function clearEditor() {
            document.getElementById('smiles-input').value = '';
            document.getElementById('prop-mw').textContent = '-';
            document.getElementById('prop-logp').textContent = '-';
            document.getElementById('prop-hbd').textContent = '-';
            document.getElementById('prop-hba').textContent = '-';
            document.getElementById('prop-tpsa').textContent = '-';
            document.getElementById('prop-rot').textContent = '-';
            document.getElementById('prop-formula').textContent = '-';
            document.getElementById('ai-suggestions').innerHTML = '<p style="color:#a0a0a0;font-size:0.8rem">Click "Analyze" to get suggestions</p>';
        }
        
        async function analyzeMolecule() {
            const smiles = document.getElementById('smiles-input').value.trim();
            if (!smiles) { alert('Please enter a SMILES string'); return; }
            
            try {
                const res = await fetch('/api/chem/properties', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ smiles: smiles })
                });
                const data = await res.json();
                
                if (data.error && data.fallback) {
                    alert('RDKit not available in this container. Using estimation only.');
                } else if (data.error) {
                    alert('Error: ' + data.error);
                    return;
                }
                
                document.getElementById('prop-mw').textContent = data.mw || '-';
                document.getElementById('prop-logp').textContent = data.logp || '-';
                document.getElementById('prop-hbd').textContent = data.hbd || '-';
                document.getElementById('prop-hba').textContent = data.hba || '-';
                document.getElementById('prop-tpsa').textContent = data.tpsa || '-';
                document.getElementById('prop-rot').textContent = data.rotatable_bonds || '-';
                document.getElementById('prop-formula').textContent = data.formula || '-';
                
                visualizeMolecule();
                
            } catch(e) { console.error('Analysis error:', e); }
            
            try {
                const res = await fetch('/api/chem/suggestions', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ smiles: smiles })
                });
                const data = await res.json();
                
                if (data.suggestions && data.suggestions.length > 0) {
                    const container = document.getElementById('ai-suggestions');
                    container.innerHTML = data.suggestions.map(s => {
                        let cls = 'info';
                        if (s.includes('good') || s.includes('Good') || s.includes('Ready')) cls = 'good';
                        if (s.includes('high') || s.includes('exceeds') || s.includes('low') || s.includes('poor')) cls = 'warning';
                        return `<div class="suggestion-item ${cls}">${s}</div>`;
                    }).join('');
                }
            } catch(e) { console.error('Suggestions error:', e); }
        }
        
        function visualizeMolecule() {
            const smiles = document.getElementById('smiles-input').value.trim();
            if (!smiles) return;
            
            try {
                const canvas = document.getElementById('preview-smiles-canvas');
                if (!canvas) return;
                canvas.width = 350;
                canvas.height = 280;
                const drawer = new SmilesDrawer.Drawer({ width: 350, height: 280, bondThickness: 2, bondLength: 30 });
                SmilesDrawer.parse(smiles, function(tree) {
                    drawer.draw(tree, canvas, 'light', false);
                });
            } catch(e) { console.error('Visualization error:', e); }
        }
        
        async function dockMolecule() {
            const smiles = document.getElementById('smiles-input').value.trim();
            if (!smiles) { alert('Please enter a SMILES string'); return; }
            
            const statusEl = document.getElementById('dock-status');
            statusEl.innerHTML = '<p style="color:#ffab00">Creating docking job...</p>';
            
            try {
                const res = await fetch('/api/chem/dock', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ smiles: smiles, receptor_id: 'default' })
                });
                const data = await res.json();
                
                if (data.error) {
                    statusEl.innerHTML = `<p style="color:#ff5252">Error: ${data.error}</p>`;
                } else {
                    statusEl.innerHTML = `<p style="color:#00c853">Job ${data.job_id} created!</p><p style="font-size:0.7rem;color:#a0a0a0">Go to Results page to view</p>`;
                }
            } catch(e) {
                statusEl.innerHTML = '<p style="color:#ff5252">Connection error</p>';
            }
        }
        
        let nglStage = null;
        
        function init3DViewer() {
            if (nglStage) return;
            
            try {
                const container = document.getElementById('ngl-container');
                if (!container) return;
                
                nglStage = NGL.createStage(container);
                
                const smiles = document.getElementById('smiles-input').value.trim();
                fetch('/api/chem/3d/test')
                    .then(r => r.json())
                    .then(data => {
                        const pdbBlob = new Blob([data.pdb], { type: 'chemical/x-pdb' });
                        nglStage.loadFile(pdbBlob, { ext: 'pdb', isOtf: false }).then(function(comp) {
                            comp.addRepresentation('ball-and-stick');
                            nglStage.autoView();
                        });
                    });
                
            } catch(e) { console.error('3D viewer error:', e); }
        }
        
        function copySmiles() {
            const smiles = document.getElementById('smiles-input').value.trim();
            if (smiles) { navigator.clipboard.writeText(smiles); }
        }
        
        function loadExample(name, smiles) {
            document.getElementById('smiles-input').value = smiles;
            analyzeMolecule();
        }
        
        async function testProvider(provider) {
            const keyInput = document.getElementById(provider + '-key');
            const statusEl = document.getElementById(provider + '-status');
            const apiKey = keyInput ? keyInput.value.trim() : '';
            
            if (!apiKey) {
                statusEl.textContent = 'Please enter API key';
                statusEl.style.color = '#ff5252';
                return;
            }
            
            statusEl.textContent = 'Testing...';
            statusEl.style.color = '#ffab00';
            
            try {
                const res = await fetch('/api/ai/test', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ provider: provider, api_key: apiKey })
                });
                const data = await res.json();
                
                if (data.success) {
                    statusEl.textContent = 'Connected!';
                    statusEl.style.color = '#00c853';
                } else {
                    statusEl.textContent = data.error || 'Failed';
                    statusEl.style.color = '#ff5252';
                }
            } catch(e) {
                statusEl.textContent = 'Connection failed';
                statusEl.style.color = '#ff5252';
            }
        }
        
        function saveProvider(provider) {
            const keyInput = document.getElementById(provider + '-key');
            const apiKey = keyInput ? keyInput.value.trim() : '';
            localStorage.setItem('api_key_' + provider, apiKey);
            const statusEl = document.getElementById(provider + '-status');
            if (apiKey) {
                statusEl.textContent = 'Saved';
                statusEl.style.color = '#00c853';
            }
        }
        
        loadStats(); loadGPUInfo(); loadJobs();
        setInterval(() => { loadStats(); loadGPUInfo(); loadJobs(); }, 30000);
    </script>
</body>
</html>'''

@app.get('/')
async def root():
    return HTMLResponse(content=HTML_CONTENT)

@app.get('/health')
async def health():
    return {'status': 'healthy', 'timestamp': datetime.utcnow().isoformat(), 'version': '2.3.1'}

@app.get('/api/stats')
async def get_stats():
    total = len(list(JOBS_DIR.glob('*.json')))
    completed = sum(1 for f in JOBS_DIR.glob('*.json') if 'completed' in f.read_text())
    return {'total_jobs': total, 'completed_jobs': completed, 'active_jobs': total - completed}

@app.get('/api/md/gpu-info')
async def gpu_info():
    import subprocess
    try:
        r = subprocess.run(['nvidia-smi'], capture_output=True, timeout=5)
        if r.returncode == 0: return {'platform': 'CUDA', 'cuda_available': True, 'opencl_available': False}
    except: pass
    return {'platform': 'CPU', 'cuda_available': False, 'opencl_available': False}

@app.get('/api/ai/providers')
async def providers():
    return {'providers': [
        {'id': 'openai', 'name': 'OpenAI'}, {'id': 'claude', 'name': 'Claude'}, {'id': 'gemini', 'name': 'Gemini'},
        {'id': 'mistral', 'name': 'Mistral'}, {'id': 'deepseek', 'name': 'DeepSeek'}, {'id': 'qwen', 'name': 'Qwen'},
        {'id': 'siliconflow', 'name': 'SiliconFlow'}, {'id': 'openrouter', 'name': 'OpenRouter'}, {'id': 'ollama', 'name': 'Ollama'}
    ]}

@app.post('/api/ai/test')
async def test_connection(req: dict):
    provider = req.get('provider', '')
    api_key = req.get('api_key', '')
    
    if not api_key:
        return {'success': False, 'error': 'No API key provided'}
    
    try:
        if provider == 'deepseek':
            import httpx
            async with httpx.AsyncClient() as client:
                resp = await client.post(
                    'https://api.deepseek.com/v1/chat/completions',
                    headers={'Authorization': f'Bearer {api_key}', 'Content-Type': 'application/json'},
                    json={'model': 'deepseek-chat', 'messages': [{'role': 'user', 'content': 'Hi'}]},
                    timeout=10.0
                )
                if resp.status_code == 200:
                    return {'success': True, 'message': 'Connection successful'}
                else:
                    return {'success': False, 'error': f'API error: {resp.status_code}'}
        else:
            return {'success': True, 'message': f'{provider} configured (basic validation)'}
    except Exception as e:
        return {'success': False, 'error': str(e)}

@app.post('/api/docking/jobs')
async def create_job():
    job_id = str(uuid.uuid4())[:8]
    job_info = {'job_id': job_id, 'type': 'docking', 'status': 'pending', 'created_at': datetime.utcnow().isoformat()}
    (JOBS_DIR / f'{job_id}.json').write_text(json.dumps(job_info))
    return job_info

@app.get('/api/docking/jobs')
async def list_jobs():
    all_jobs = [json.loads(f.read_text()) for f in JOBS_DIR.glob('*.json') if f.stat().st_size > 0]
    return {'jobs': sorted(all_jobs, key=lambda x: x.get('created_at', ''), reverse=True)}

@app.get('/api/docking/jobs/{job_id}')
async def get_job(job_id: str):
    f = JOBS_DIR / f'{job_id}.json'
    if f.exists(): return json.loads(f.read_text())
    return {'job_id': job_id, 'status': 'not_found'}

@app.post('/api/ai/chat')
async def chat(req: dict):
    return {'response': f'BioDockify AI ({req.get("provider", "demo")}): Configure API keys in Settings for full functionality.', 'provider': req.get('provider', 'demo')}

@app.post('/api/chem/properties')
async def get_properties(req: dict):
    smiles = req.get('smiles', '')
    if not smiles:
        return {'error': 'No SMILES provided'}
    
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Crippen, Lipinski
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}
        
        mol_h = Chem.AddHs(mol)
        atom_count = mol.GetNumAtoms()
        bond_count = mol.GetNumBonds()
        
        return {
            'smiles': smiles,
            'mw': round(Descriptors.MolWt(mol), 2),
            'logp': round(Crippen.MolLogP(mol), 2),
            'hbd': Lipinski.NumHDonors(mol),
            'hba': Lipinski.NumHAcceptors(mol),
            'tpsa': round(Descriptors.TPSA(mol), 2),
            'rotatable_bonds': Lipinski.NumRotatableBonds(mol),
            'aromatic_rings': Lipinski.NumAromaticRings(mol),
            'atoms': atom_count,
            'bonds': bond_count,
            'formula': Chem.rdMolDescriptors.CalcMolFormula(mol),
            'valid': True
        }
    except ImportError:
        return {'error': 'RDKit not available', 'fallback': True}
    except Exception as e:
        return {'error': str(e)}

@app.post('/api/chem/suggestions')
async def get_suggestions(req: dict):
    smiles = req.get('smiles', '')
    if not smiles:
        return {'suggestions': []}
    
    suggestions = []
    
    try:
        from rdkit import Chem
        from rdkit.Chem import Lipinski, Crippen
        
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = Lipinski.NumHDonors(mol)
            hba = Lipinski.NumHAcceptors(mol)
            rotatable = Lipinski.NumRotatableBonds(mol)
            
            if mw > 500:
                suggestions.append(f'Molecular weight ({mw:.1f}) is high. Consider reducing to <500 Da for better oral bioavailability.')
            if mw < 150:
                suggestions.append(f'Molecular weight ({mw:.1f}) is very low. May have rapid clearance.')
            
            if logp > 5:
                suggestions.append(f'LogP ({logp:.2f}) is high. May have poor solubility. Consider adding polar groups.')
            if logp < 0:
                suggestions.append(f'LogP ({logp:.2f}) is low. May have poor membrane penetration.')
            
            if hbd > 5:
                suggestions.append(f'H-bond donors ({hbd}) exceeds Lipinski rule of 5 limit (5).')
            if hba > 10:
                suggestions.append(f'H-bond acceptors ({hba}) exceeds Lipinski rule of 5 limit (10).')
            
            if rotatable > 10:
                suggestions.append(f'Rotatable bonds ({rotatable}) is high. May reduce binding specificity.')
            
            if 'c1ccccc1' in smiles or 'c1ccncc1' in smiles:
                suggestions.append('Contains aromatic ring - good for hydrophobic interactions with binding pocket.')
            
            if 'O' in smiles and 'N' not in smiles:
                suggestions.append('Consider adding nitrogen atoms to improve binding affinity.')
            
            if 'N' in smiles and 'O' not in smiles:
                suggestions.append('Consider adding oxygen atoms to improve solubility.')
            
            if mw < 400 and logp > 2 and logp < 4:
                suggestions.append('Good drug-like properties detected!')
            
            suggestions.append('Consider bioisosteric replacement for improved ADMET properties.')
            
    except ImportError:
        suggestions.append('RDKit not available - basic suggestions only.')
    except Exception as e:
        suggestions.append(f'Analysis complete. Error: {str(e)[:50]}')
    
    if len(suggestions) == 0:
        suggestions.append('Molecule appears within drug-like property ranges. Ready for docking!')
    
    return {'suggestions': suggestions, 'smiles': smiles}

@app.post('/api/chem/dock')
async def dock_molecule(req: dict):
    smiles = req.get('smiles', '')
    receptor_id = req.get('receptor_id', 'default')
    
    if not smiles:
        return {'error': 'No SMILES provided'}
    
    job_id = str(uuid.uuid4())[:8]
    job_dir = JOBS_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}
        
        mol_h = Chem.AddHs(mol)
        pdb_path = job_dir / 'ligand.pdb'
        Chem.MolToPDBFile(mol_h, str(pdb_path))
        
        job_info = {
            'job_id': job_id,
            'type': 'docking',
            'smiles': smiles,
            'receptor_id': receptor_id,
            'status': 'pending',
            'ligand_path': str(pdb_path),
            'created_at': datetime.utcnow().isoformat()
        }
        (job_dir / 'info.json').write_text(json.dumps(job_info, indent=2))
        (JOBS_DIR / f'{job_id}.json').write_text(json.dumps(job_info))
        
        return {'job_id': job_id, 'status': 'created', 'message': 'Docking job created'}
        
    except ImportError:
        job_info = {
            'job_id': job_id,
            'type': 'docking',
            'smiles': smiles,
            'status': 'created_no_rdkit',
            'created_at': datetime.utcnow().isoformat()
        }
        (JOBS_DIR / f'{job_id}.json').write_text(json.dumps(job_info))
        return {'job_id': job_id, 'status': 'created', 'message': 'Job created (RDKit not available for prep)'}
    except Exception as e:
        return {'error': str(e)}

@app.get('/api/chem/3d/{job_id}')
async def get_3d_structure(job_id: str):
    job_dir = JOBS_DIR / job_id
    pdb_path = job_dir / 'ligand.pdb'
    
    if pdb_path.exists():
        return {'pdb': pdb_path.read_text(), 'job_id': job_id}
    
    sample_pdb = """HETATM    1  C   UNL     1      -1.500   0.500   0.000  1.00  0.00           C
HETATM    2  C   UNL     1       0.500   0.500   0.000  1.00  0.00           C
HETATM    3  C   UNL     1       1.000   0.000   1.000  1.00  0.00           C
HETATM    4  C   UNL     1       0.000   0.000  -1.000  1.00  0.00           C
HETATM    5  O   UNL     1      -2.000   1.000  -1.000  1.00  0.00           O
HETATM    6  O   UNL     1      -2.000   0.000   1.000  1.00  0.00           O
CONECT    1    2    4    5
CONECT    2    1    3
CONECT    3    2
CONECT    4    1
CONECT    5    1
CONECT    6    1
END"""
    
    return {'pdb': sample_pdb, 'job_id': job_id, 'sample': True}

if __name__ == '__main__':
    uvicorn.run(app, host='0.0.0.0', port=8000)
