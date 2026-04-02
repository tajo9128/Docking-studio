# -*- mode: python ; coding: utf-8 -*-
"""
PyInstaller spec for BioDockify Studio AI — Windows build
Entry point: app/launcher.py (serves ALL backend routes + React frontend)
"""
import sys
from pathlib import Path

block_cipher = None

# Hidden imports: packages PyInstaller can't auto-detect
hidden_imports = [
    # RDKit
    'rdkit', 'rdkit.Chem', 'rdkit.Chem.AllChem', 'rdkit.Chem.rdMolDescriptors',
    'rdkit.Chem.rdMolTransforms', 'rdkit.Chem.MolStandardize',
    'rdkit.Chem.Descriptors', 'rdkit.Chem.Fingerprints',
    # OpenMM
    'openmm', 'openmm.app', 'openmm.unit', 'openmm.openmm',
    # CrewAI
    'crewai', 'crewai_tools', 'crewai.memory', 'crewai.task', 'crewai.agent',
    # Vector DB
    'chromadb', 'chromadb.api',
    # Scientific
    'numpy', 'numpy.core', 'numpy.linalg', 'scipy', 'scipy.spatial', 'scipy.stats',
    'pandas', 'mdtraj', 'mordred', 'xgboost',
    'sklearn', 'sklearn.ensemble', 'sklearn.gaussian_process', 'sklearn.metrics',
    # Web
    'uvicorn', 'uvicorn.main', 'uvicorn.config',
    'fastapi', 'fastapi.staticfiles', 'fastapi.middleware', 'fastapi.routing',
    'starlette', 'starlette.staticfiles', 'starlette.middleware',
    # Data
    'pydantic', 'pydantic_settings', 'pydantic.fields', 'pydantic.main',
    # HTTP
    'httpx', 'httpcore', 'requests', 'urllib3',
    # Utils
    'psutil', 'joblib', 'tqdm', 'python-dotenv', 'click', 'rich',
    # Crypto
    'secrets', 'hashlib', 'cryptography',
    # BioDockify modules
    'docking_engine', 'analysis', 'db', 'pharmacophore', 'classroom',
    'ai.llm_router', 'ai.config', 'ai.offline_engine',
    'crew.agents', 'crew.crews', 'crew.flows', 'crew.memory',
    'crew.meta_optimizer', 'crew.active_learning', 'crew.nl_compiler',
    'crew.critique_agent', 'crew.knowledge_graph',
    'crew.tools.docking_tools', 'crew.tools.chemistry_tools',
    'crew.tools.pharmacophore_tools', 'crew.tools.admet_tools',
    'crew.tools.analysis_tools', 'crew.tools.data_tools',
    'crew.tools.notification_tools', 'crew.tools.base',
]

# Data files to bundle
datas = [
    ('frontend/dist', 'frontend/dist'),
    ('backend/static', 'backend/static'),
]

a = Analysis(
    ['app/launcher.py'],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=hidden_imports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['tkinter', 'matplotlib.tests', 'pytest', 'torchvision', 'jupyter'],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='BioDockify',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='BioDockify',
)
