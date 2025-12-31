# -*- mode: python ; coding: utf-8 -*-
"""
BioDockify Docking Studio - PyInstaller Build Specification for Windows
Production build for Windows 11 (x64)
"""

block_cipher = None

# Custom data files to include
# Format: (Source Path, Destination Folder)
added_files = [
    ('src/templates/dock_vina.conf', 'templates'),
    ('LICENSE', '.'),
    ('src/ui/styles', 'ui/styles')
]

datas = []
binaries = []

from PyInstaller.utils.hooks import collect_all

# Hidden imports
hidden_imports = [
    'src.database',
    'src.config',
    'src.utils.log_utils'
]

# Collect all PyQt6 components (binaries, datas, hidden imports)
# This fixes "ModuleNotFoundError: No module named 'PyQt6'"
tmp_ret = collect_all('PyQt6')
datas += tmp_ret[0]
binaries += tmp_ret[1]
hidden_imports += tmp_ret[2]

a = Analysis(
    ['src/main.py'],
    pathex=[],
    datas=added_files + datas,
    hiddenimports=hidden_imports,
    binaries=binaries,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        'tests',
        'tests.*',
        'setup',
        'setup.*',
        'conftest'
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

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
    icon='src/ui/styles/icon.ico'
)
