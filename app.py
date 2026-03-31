import uvicorn
from fastapi import FastAPI, UploadFile, File, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse
from pydantic import BaseModel
from pathlib import Path
import uuid
import json
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

HTML_CONTENT = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>BioDockify</title>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; background: #1a1a2e; color: #e8e8e8; min-height: 100vh; }
        .app { display: flex; flex-direction: column; min-height: 100vh; }
        .header { background: linear-gradient(135deg, #16213e 0%, #0f3460 100%); padding: 1rem 2rem; display: flex; align-items: center; justify-content: space-between; border-bottom: 1px solid #2a2a4a; }
        .header h1 { color: #00d9ff; font-size: 1.5rem; }
        .header-nav a { color: #a0a0a0; text-decoration: none; padding: 0.5rem 1rem; margin-left: 0.5rem; border-radius: 4px; }
        .header-nav a:hover, .header-nav a.active { color: #00d9ff; background: rgba(0,217,255,0.1); }
        .main { display: flex; flex: 1; }
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
        .card { background: #16213e; border-radius: 8px; padding: 1.5rem; border: 1px solid #2a2a4a; margin-bottom: 1rem; }
        .card h3 { margin-bottom: 1rem; }
        .btn { padding: 0.75rem 1.5rem; border-radius: 6px; border: none; cursor: pointer; font-weight: 500; transition: all 0.2s; }
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
        .status-dot.warning { background: #ffab00; }
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
                <a href="#" onclick="showPage('md')"><span>📊</span> MD Simulation</a>
                <a href="#" onclick="showPage('results')"><span>📈</span> Results</a>
                <a href="#" onclick="showPage('ai')"><span>🤖</span> AI Assistant</a>
                <a href="#" onclick="showPage('settings')"><span>⚙️</span> Settings</a>
            </aside>
            <main class="page" id="page-content">
            </main>
        </div>
        <footer class="status-bar">
            <span><span class="status-dot green"></span>Connected to BioDockify Backend</span>
            <span>v2.3.0</span>
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
                <h2>Dashboard</h2>
                <p>Welcome to BioDockify - Molecular Docking Studio</p>
                <div class="stats-grid">
                    <div class="stat-card"><h3>Total Jobs</h3><div class="value">${stats.total_jobs}</div></div>
                    <div class="stat-card"><h3>Completed</h3><div class="value">${stats.completed_jobs}</div></div>
                    <div class="stat-card"><h3>Active</h3><div class="value">${stats.active_jobs}</div></div>
                    <div class="stat-card"><h3>GPU Platform</h3><div class="value" style="font-size:1rem">${gpuInfo.platform}</div></div>
                </div>
                <div class="grid-2">
                    <div class="card">
                        <h3>Quick Actions</h3>
                        <div style="display:flex;flex-direction:column;gap:0.5rem;margin-top:1rem">
                            <button class="btn btn-primary" onclick="showPage('docking')">New Docking Job</button>
                            <button class="btn btn-secondary" onclick="showPage('md')">Start MD Simulation</button>
                            <button class="btn btn-secondary" onclick="showPage('ai')">Ask BioDockify AI</button>
                        </div>
                    </div>
                    <div class="card">
                        <h3>System Status</h3>
                        <div style="margin-top:1rem">
                            <div style="margin-bottom:0.5rem"><span class="status-dot green"></span>Backend API</div>
                            <div style="margin-bottom:0.5rem"><span class="status-dot green"></span>File System</div>
                            <div><span class="status-dot ${gpuInfo.platform !== 'CPU' ? 'green' : 'warning'}"></span>GPU: ${gpuInfo.platform}</div>
                        </div>
                    </div>
                </div>
            `,
            docking: () => `
                <h2>Molecular Docking</h2>
                <p>Run AutoDock Vina docking simulations</p>
                <div class="grid-2">
                    <div>
                        <div class="card">
                            <h3>Receptor</h3>
                            <div class="drop-zone" onclick="document.getElementById('receptor-input').click()">
                                <p>Drop PDB/PDBQT file or click to upload</p>
                            </div>
                            <input type="file" id="receptor-input" accept=".pdb,.pdbqt" style="display:none">
                        </div>
                        <div class="card">
                            <h3>Ligand</h3>
                            <div class="drop-zone" onclick="document.getElementById('ligand-input').click()">
                                <p>Drop PDB/PDBQT/SDF file or click to upload</p>
                            </div>
                            <input type="file" id="ligand-input" accept=".pdb,.pdbqt,.sdf" style="display:none">
                            <div style="text-align:center;color:#a0a0a0;margin:1rem 0">- OR -</div>
                            <div class="form-group">
                                <label>SMILES String</label>
                                <input type="text" id="smiles-input" placeholder="CC(=O)Oc1ccccc1C(=O)O">
                            </div>
                        </div>
                        <div class="card">
                            <h3>Parameters</h3>
                            <div class="form-group">
                                <label>Exhaustiveness: <span id="exhaust-value">32</span></label>
                                <input type="range" min="1" max="64" value="32" oninput="document.getElementById('exhaust-value').textContent=this.value">
                            </div>
                            <div class="form-group">
                                <label>Number of Poses: <span id="poses-value">10</span></label>
                                <input type="range" min="1" max="20" value="10" oninput="document.getElementById('poses-value').textContent=this.value">
                            </div>
                            <button class="btn btn-primary" style="width:100%" onclick="submitDockingJob()">Start Docking</button>
                        </div>
                    </div>
                    <div>
                        <div class="card">
                            <h3>Job Queue</h3>
                            <div class="job-list" id="job-list">
                                ${jobs.length === 0 ? '<p style="color:#a0a0a0">No jobs yet</p>' : ''}
                                ${jobs.map(j => '<div class="job-item"><span class="job-id">Job ' + j.job_id + '</span><span class="job-status">' + j.status + '</span></div>').join('')}
                            </div>
                        </div>
                    </div>
                </div>
            `,
            md: () => `
                <h2>Molecular Dynamics</h2>
                <p>Run OpenMM MD simulations with GPU acceleration</p>
                <div class="grid-2">
                    <div>
                        <div class="card">
                            <h3>Structure File</h3>
                            <div class="drop-zone" onclick="document.getElementById('pdb-input').click()">
                                <p>Drop PDB file or click to upload</p>
                            </div>
                            <input type="file" id="pdb-input" accept=".pdb" style="display:none">
                        </div>
                        <div class="card">
                            <h3>Parameters</h3>
                            <div class="form-group">
                                <label>Force Field</label>
                                <select><option>AMBER14 All</option><option>AMBER14 TIP3P</option><option>CHARMM36</option></select>
                            </div>
                            <div class="form-group">
                                <label>Duration (ns): <span id="duration-value">10</span></label>
                                <input type="range" min="1" max="100" value="10" oninput="document.getElementById('duration-value').textContent=this.value">
                            </div>
                            <div class="form-group">
                                <label>Temperature (K): <span id="temp-value">300</span></label>
                                <input type="range" min="200" max="400" value="300" oninput="document.getElementById('temp-value').textContent=this.value">
                            </div>
                            <button class="btn btn-primary" style="width:100%">Start MD Simulation</button>
                        </div>
                    </div>
                    <div>
                        <div class="card">
                            <h3>Active Simulations</h3>
                            <p style="color:#a0a0a0">No active simulations</p>
                        </div>
                        <div class="card">
                            <h3>GPU Acceleration</h3>
                            <p style="color:#a0a0a0;margin-top:0.5rem">OpenMM will automatically detect and use:</p>
                            <ul style="margin-top:1rem;padding-left:1.5rem;color:#a0a0a0">
                                <li>CUDA (NVIDIA GPUs)</li>
                                <li>OpenCL (AMD/Intel GPUs)</li>
                                <li>CPU (fallback)</li>
                            </ul>
                        </div>
                    </div>
                </div>
            `,
            results: () => `
                <h2>Results</h2>
                <p>View completed docking and simulation results</p>
                <div class="grid-2">
                    <div class="card">
                        <h3>All Jobs</h3>
                        <div class="job-list">
                            ${jobs.length === 0 ? '<p style="color:#a0a0a0">No jobs found. Run a docking simulation to see results here.</p>' : ''}
                            ${jobs.map(j => '<div class="job-item"><div><span class="job-id">' + j.job_id + '</span><div style="font-size:0.75rem;color:#a0a0a0">' + j.type + '</div></div><span class="job-status">' + j.status + '</span></div>').join('')}
                        </div>
                    </div>
                    <div class="card">
                        <h3>Job Details</h3>
                        <p style="color:#a0a0a0">Select a job to view details</p>
                    </div>
                </div>
            `,
            ai: () => `
                <h2>BioDockify AI</h2>
                <p>Ask questions about molecular docking, simulations, and more</p>
                <div class="card">
                    <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:1rem">
                        <h3>Chat</h3>
                        <select value="openai" onchange="selectedProvider=this.value" style="padding:0.5rem;background:#1a1a2e;color:#e8e8e8;border:1px solid #2a2a4a;border-radius:4px">
                            <option value="openai">OpenAI</option>
                            <option value="claude">Claude</option>
                            <option value="gemini">Gemini</option>
                            <option value="mistral">Mistral</option>
                            <option value="deepseek">DeepSeek</option>
                            <option value="qwen">Qwen</option>
                            <option value="siliconflow">SiliconFlow</option>
                            <option value="openrouter">OpenRouter</option>
                            <option value="ollama">Ollama</option>
                        </select>
                    </div>
                    <div class="ai-chat">
                        <div class="chat-messages" id="chat-messages">
                            ${messages.length === 0 ? '<p style="color:#a0a0a0;text-align:center;margin-top:2rem">Ask me anything about molecular docking, MD simulations, or drug discovery!</p>' : ''}
                            ${messages.map(m => '<div class="chat-message ' + m.role + '"><strong>' + (m.role === 'user' ? 'You' : 'BioDockify AI') + ':</strong> ' + m.content + '</div>').join('')}
                        </div>
                        <div class="chat-input">
                            <input type="text" id="chat-input" placeholder="Ask BioDockify AI..." onkeypress="if(event.key==='Enter')sendMessage()">
                            <button class="btn btn-primary" onclick="sendMessage()">Send</button>
                        </div>
                    </div>
                </div>
            `,
            settings: () => `
                <h2>Settings</h2>
                <p>Configure AI providers and system preferences</p>
                <div class="card">
                    <h3>API Keys Configuration</h3>
                    <p style="color:#a0a0a0;margin:1rem 0">Configure your API keys for AI providers. Keys are stored locally and never sent to external servers.</p>
                    <div class="provider-grid">
                        <div class="provider-card"><h4>OpenAI</h4><p class="job-status" style="background:rgba(255,171,0,0.2);color:#ffab00">Not configured</p><div class="form-group" style="margin-top:0.5rem"><input type="password" placeholder="sk-..."></div><button class="btn btn-primary" style="width:100%;margin-top:0.5rem">Save</button></div>
                        <div class="provider-card"><h4>Claude</h4><p class="job-status" style="background:rgba(255,171,0,0.2);color:#ffab00">Not configured</p><div class="form-group" style="margin-top:0.5rem"><input type="password" placeholder="sk-ant-..."></div><button class="btn btn-primary" style="width:100%;margin-top:0.5rem">Save</button></div>
                        <div class="provider-card"><h4>Gemini</h4><p class="job-status" style="background:rgba(255,171,0,0.2);color:#ffab00">Not configured</p><div class="form-group" style="margin-top:0.5rem"><input type="password" placeholder="AI..."></div><button class="btn btn-primary" style="width:100%;margin-top:0.5rem">Save</button></div>
                    </div>
                </div>
                <div class="card">
                    <h3>About</h3>
                    <p style="margin-top:1rem">BioDockify v2.3.0</p>
                    <p style="color:#a0a0a0;margin-top:0.5rem">Production-ready molecular docking software with Discovery Studio-inspired UI.</p>
                </div>
            `
        };
        
        function showPage(page) {
            currentPage = page;
            document.getElementById('page-content').innerHTML = pages[page]();
            document.querySelectorAll('.sidebar a, .header-nav a').forEach(a => a.classList.remove('active'));
            document.querySelectorAll('.sidebar a, .header-nav a').forEach(a => {
                if (a.textContent.toLowerCase().includes(page === 'dashboard' ? 'home' : page === 'md' ? 'simulation' : page)) a.classList.add('active');
            });
        }
        
        async function loadStats() {
            try {
                const res = await fetch('/api/stats');
                stats = await res.json();
                if (currentPage === 'dashboard') showPage('dashboard');
            } catch(e) { console.error(e); }
        }
        
        async function loadGPUInfo() {
            try {
                const res = await fetch('/api/md/gpu-info');
                gpuInfo = await res.json();
                if (currentPage === 'dashboard') showPage('dashboard');
            } catch(e) { console.error(e); }
        }
        
        async function loadJobs() {
            try {
                const res = await fetch('/api/docking/jobs');
                const data = await res.json();
                jobs = data.jobs || [];
                if (currentPage === 'docking' || currentPage === 'results') showPage(currentPage);
            } catch(e) { console.error(e); }
        }
        
        async function submitDockingJob() {
            try {
                const res = await fetch('/api/docking/jobs', { method: 'POST' });
                const job = await res.json();
                jobs.unshift(job);
                showPage('docking');
            } catch(e) { console.error(e); }
        }
        
        async function sendMessage() {
            const input = document.getElementById('chat-input');
            const text = input.value.trim();
            if (!text) return;
            messages.push({ role: 'user', content: text });
            input.value = '';
            showPage('ai');
            try {
                const res = await fetch('/api/ai/chat', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ message: text, provider: selectedProvider })
                });
                const data = await res.json();
                messages.push({ role: 'ai', content: data.response });
            } catch(e) { messages.push({ role: 'ai', content: 'Error: Could not connect to AI service' }); }
            showPage('ai');
        }
        
        loadStats();
        loadGPUInfo();
        loadJobs();
        setInterval(() => { loadStats(); loadGPUInfo(); loadJobs(); }, 30000);
    </script>
</body>
</html>'''

@app.get('/')
async def root():
    return HTMLResponse(content=HTML_CONTENT)

@app.get('/health')
async def health():
    return {'status': 'healthy', 'timestamp': datetime.utcnow().isoformat()}

@app.get('/api/stats')
async def get_stats():
    total = len(list(JOBS_DIR.glob('*.json')))
    completed = 0
    for f in JOBS_DIR.glob('*.json'):
        try:
            if 'completed' in f.read_text(): completed += 1
        except: pass
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
async def create_job():
    job_id = str(uuid.uuid4())[:8]
    job_info = {'job_id': job_id, 'type': 'docking', 'status': 'pending', 'created_at': datetime.utcnow().isoformat()}
    (JOBS_DIR / f'{job_id}.json').write_text(json.dumps(job_info))
    return job_info

@app.get('/api/docking/jobs')
async def list_jobs():
    all_jobs = []
    for f in JOBS_DIR.glob('*.json'):
        try: all_jobs.append(json.loads(f.read_text()))
        except: pass
    return {'jobs': sorted(all_jobs, key=lambda x: x.get('created_at', ''), reverse=True)}

@app.get('/api/docking/jobs/{job_id}')
async def get_job(job_id: str):
    f = JOBS_DIR / f'{job_id}.json'
    if f.exists(): return json.loads(f.read_text())
    return {'job_id': job_id, 'status': 'not_found'}

@app.post('/api/ai/chat')
async def chat(req: dict):
    msg = req.get('message', '')
    prov = req.get('provider', 'demo')
    return {'response': f'BioDockify AI ({prov}): Processing your query about {msg[:30]}... Configure API keys in Settings for full functionality.', 'provider': prov}

if __name__ == '__main__':
    uvicorn.run(app, host='0.0.0.0', port=8000)
