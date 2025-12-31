# -*- mode: python ; coding: utf-8 -*-
"""
BioDockify Docking Studio - PyInstaller Build Specification for Windows
Production build for Windows 11 (x64)
"""

block_cipher = None

# ======================
# Force Copy Strategy (v1.0.11)
# ======================
import os
import PyQt6.QtCore
# Find the actual directory of PyQt6 on the build machine
qt_dir = os.path.dirname(PyQt6.QtCore.__file__) # e.g., site-packages/PyQt6/QtCore.pyd -> site-packages/PyQt6
qt_package_root = os.path.dirname(qt_dir) # site-packages/PyQt6 ? No, QtCore is inside PyQt6.
# If qt_dir is ".../site-packages/PyQt6", then we want to copy qt_dir to "PyQt6"
# Actually, QtCore.__file__ is usually .../PyQt6/QtCore.pyd (Windows) or .../PyQt6/QtCore.so
# So dirname is the PyQt6 folder.

print(f"DEBUG: Found PyQt6 at {qt_dir}")

# We manually add the entire folder as data
# Format: (Source, Dest)
# We want '.../PyQt6' -> 'PyQt6' in the dist folder
added_files = [
    ('src/templates/dock_vina.conf', 'templates'),
    ('LICENSE', '.'),
    ('src/ui/styles', 'ui/styles'),
    (qt_dir, 'PyQt6') 
]

# Standard hidden imports just to be safe
hidden_imports = [
    'PyQt6',
    'PyQt6.QtCore',
    'PyQt6.QtGui',
    'PyQt6.QtWidgets',
    'PyQt6.QtWebEngineWidgets',
    'PyQt6.QtNetwork',
    'PyQt6.QtPrintSupport',
    'fastapi',
    'uvicorn',
    'pydantic',
    'docker',
    'sqlalchemy',
    'src.database',
    'src.config',
    'src.utils.log_utils'
]

a = Analysis(
    ['src/main.py'],
    pathex=['.'],
    binaries=[],
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
