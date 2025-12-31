"""
BioDockify Docking Studio - PyInstaller Build Specification for Windows
Production build for Windows 11 (x64)
"""

block_cipher = None

# ======================
# Clean Spec (v1.0.12)
# ======================

# ======================
# Clean Spec (v1.0.17)
# ======================

# Initialize variable to prevent NameError
pyqt_path = None

import os
import site

# Attempt to find PyQt6 via direct import
try:
    import PyQt6
    pyqt_path = os.path.dirname(PyQt6.__file__)
    print(f"Spec DEBUG: Found PyQt6 via import at {pyqt_path}")
except ImportError:
    print("Spec WARNING: Could not import PyQt6. Falling back to site-packages scan.")
    try:
        site_dirs = site.getsitepackages()
        for d in site_dirs:
            candidate = os.path.join(d, 'PyQt6')
            if os.path.exists(candidate):
                pyqt_path = candidate
                print(f"Spec DEBUG: Found PyQt6 via scan at {pyqt_path}")
                break
    except Exception as e:
        print(f"Spec WARNING: Site scan failed: {e}")
except Exception as e:
    print(f"Spec WARNING: Unexpected error finding PyQt6: {e}")

# Explicit hidden imports (Option A)
hidden_imports = [
    'PyQt6',
    'PyQt6.QtCore',
    'PyQt6.QtGui',
    'PyQt6.QtWidgets',
    'PyQt6.QtWebEngineWidgets',
    'PyQt6.QtNetwork',
    'PyQt6.QtPrintSupport',
    'PyQt6.sip',  # often critical!
    'fastapi',
    'uvicorn',
    'pydantic',
    'docker',
    'sqlalchemy',
    'src.database',
    'src.config',
    'src.utils.log_utils'
]

# (Old import logic removed - handled at top of file)

# ======================
# Data Files
# ======================
added_files = [
    ('src/templates/dock_vina.conf', 'templates'),
    ('LICENSE', '.'),
    ('src/ui/styles', 'ui/styles')
]

if pyqt_path:
    added_files.append( (pyqt_path, 'PyQt6') )

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
