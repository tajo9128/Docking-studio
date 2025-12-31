"""
BioDockify Docking Studio - PyInstaller Build Specification for Windows
Production build for Windows 11 (x64)
Fixed for v1.0.36
"""

block_cipher = None

# ======================
# Clean Spec (v1.0.35 - Final Fix)
# ======================

from PyInstaller.utils.hooks import collect_all
from PyInstaller.building.build_main import Analysis, PYZ, EXE, COLLECT, MERGE

# Initialize variable to prevent NameError
pyqt_path = None

import os
import sys

import site

# CORRECTLY DEFINE tmp_ret NOW:
# Collect PyQt6
tmp_ret = collect_all('PyQt6')
datas_qt, binaries_qt, hiddenimports_qt = tmp_ret

# Manual Collection of Web Frameworks (Robust)
datas_web = []
hiddenimports_web = []

def collect_package_manual(package_name):
    """Manually find and collect package"""
    import importlib.util
    try:
        spec = importlib.util.find_spec(package_name)
        if spec and spec.submodule_search_locations:
            path = spec.submodule_search_locations[0]
            print(f"Manual collect: Found {package_name} at {path}")
            return [(path, package_name)]
        else:
            print(f"Manual collect: Could not find path for {package_name}")
            return []
    except ImportError:
        print(f"Manual collect: Failed to import {package_name}")
        return []

datas_web.extend(collect_package_manual('fastapi'))
datas_web.extend(collect_package_manual('uvicorn'))
datas_web.extend(collect_package_manual('starlette'))
datas_web.extend(collect_package_manual('email_validator'))
datas_web.extend(collect_package_manual('python_multipart'))

# Explicit hidden imports (Base)
base_hidden_imports = [
    'PyQt6.sip', 
    'fastapi',
    'uvicorn',
    'pydantic',
    'docker',
    'sqlalchemy',
    'starlette',
    'email_validator',
    'python_multipart',
    'src.database',
    'src.config',
    'src.utils.log_utils',
    'src.utils.docker_utils',
    'src', # Root package
    'ui',
    'ui.main_window',
    'ui.upload_widget',
    'ui.configuration_widget',
    'ui.progress_widget',
    'ui.results_widget',
    'ui.agent_zero_widget',
    # Alias with src prefix just in case
    'src.ui',
    'src.ui.main_window',
    'src.ui.upload_widget',
    'src.ui.configuration_widget',
    'src.ui.progress_widget',
    'src.ui.results_widget',
    'src.ui.agent_zero_widget'
]

# Merge all lists
hidden_imports = list(set(
    hiddenimports_qt + 
    hiddenimports_web +
    base_hidden_imports
))

binaries = binaries_qt

# ======================
# Data Files - FIXED v1.0.40
# ======================

# Find icon file
icon_paths = [
    'src/ui/styles/icon.ico',
    'biodockify_icon.ico',
    'icon.ico',
    os.path.abspath('icon.ico'),
]

icon_path = None
for path in icon_paths:
    if os.path.exists(path):
        icon_path = path
        print(f"Found icon: {path}")
        break

if not icon_path:
    print("WARNING: No icon file found, building without icon")
    icon_path = None

# Clean data files - NO duplicates
added_files = [
    ('src/templates/dock_vina.conf', 'templates'),
    ('LICENSE', '.'),
    ('src/ui/styles', 'ui/styles'),
] + datas_qt + datas_web

if icon_path:
    added_files.append((icon_path, '.'))

if pyqt_path:
    added_files.append((pyqt_path, 'PyQt6'))

# Analysis
a = Analysis(
    ['src/biodockify_main.py'],
    pathex=['.', 'src'],
    binaries=binaries,
    datas=added_files,
    hiddenimports=hidden_imports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        'tests',
        'tests.*',
        'setup',
        'setup.*',
        'conftest',
        'pytest',
        'pytest.*',
        'mypy',
        'black',
        'flake8'
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

# PYZ
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

# EXE
exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='BioDockify-Docking-Studio',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=icon_path  # ✅ Safe: None if not found
)

# COLLECT (creates the distribution)
coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='BioDockify-Docking-Studio',
)
