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
    'src.utils.log_utils',
    'PyQt6',
    'PyQt6.QtCore',
    'PyQt6.QtGui',
    'PyQt6.QtWidgets',
    'PyQt6.QtWebEngineCore',
    'PyQt6.QtWebEngineWidgets',
    'PyQt6.QtNetwork',
    'PyQt6.QtPrintSupport',
    'PyQt6.sip',
]

# Note: PyInstaller 6.x should handle PyQt6 hooks automatically now that we upgraded requirements.
# We keep the explicit list just in case.

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
