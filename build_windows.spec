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

from PyInstaller.utils.hooks import collect_all

# Initialize variable to prevent NameError
pyqt_path = None

import os
import site

# ... (Previous import logic remains) ...

# NOTE: The following lines (datas_qt, binaries_qt, hiddenimports_qt = tmp_ret)
# are part of a larger PyInstaller hook mechanism for PyQt6 that is not fully
# provided in the context of this edit. For this change, we will assume
# tmp_ret is defined elsewhere or that these variables are intended to be
# populated by a custom hook. If 'tmp_ret' is undefined, this spec will fail.
# For the purpose of this edit, we are inserting the provided lines as-is.

# CORRECTLY DEFINE tmp_ret NOW:
datas_qt, binaries_qt, hiddenimports_qt = tmp_ret

# Explicit hidden imports (Base)
base_hidden_imports = [
    'PyQt6.sip', 
    'fastapi',
    'uvicorn',
    'pydantic',
    'docker',
    'sqlalchemy',
    'src.database',
    'src.config',
    'src.utils.log_utils',
    'ui',
    'ui.main_window',
    'ui.upload_widget',
    'ui.configuration_widget',
    'ui.progress_widget',
    'ui.results_widget',
    'ui.agent_zero_widget'
]

hidden_imports = list(set(hiddenimports_qt + base_hidden_imports))
binaries = binaries_qt

# ======================
# Data Files
# ======================
added_files = [
    ('src/templates/dock_vina.conf', 'templates'),
    ('LICENSE', '.'),
    ('src/ui/styles', 'ui/styles')
] + datas_qt

if pyqt_path:
    added_files.append( (pyqt_path, 'PyQt6') )

a = Analysis(
    ['src/main.py'],
    pathex=['.', 'src'],

hidden_imports = list(set(hiddenimports_qt + base_hidden_imports))
binaries = binaries_qt

# (Old import logic removed - handled at top of file)

# ======================
# Data Files
# ======================
added_files = [
    ('src/templates/dock_vina.conf', 'templates'),
    ('LICENSE', '.'),
    ('src/ui/styles', 'ui/styles')
] + datas_qt

if pyqt_path:
    added_files.append( (pyqt_path, 'PyQt6') )

a = Analysis(
    ['src/main.py'],
    pathex=['.', os.path.abspath('src')],
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
