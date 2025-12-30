# -*- mode: python ; coding: utf-8 -*-
"""
BioDockify Docking Studio - PyInstaller Build Specification for macOS Apple Silicon
Production build for macOS Apple Silicon (ARM64)
"""

import sys
from PyInstaller.utils.hooks import collect_data_files, collect_submodules, collect_system_data_binaries
from PyInstaller.utils.hooks.qt import get_qt_binaries
from PyInstaller.building.build_main import Analysis

block_cipher = None
noencrypt = True

a = Analysis(['src/main.py'])

# Add hidden imports
a.datas += [('BioDockify/src/templates/dock_vina.conf', 'templates/dock_vina.conf', 'DATA')]
a.datas += [('BioDockify/LICENSE', 'LICENSE', 'DATA')]

# Add UI files
a.datas += [('BioDockify/src/ui/styles', 'ui/styles', 'DATA')]

# Add database initialization
hiddenimports = [
    'BioDockify.src.database',
    'BioDockify.src.config',
    'BioDockify.src.utils.log_utils'
]

# Collect PyQt6 binaries
pyqt6_binaries = get_qt_binaries('PyQt6')
a.binaries += pyqt6_binaries

# Exclude unused modules
excludes = [
    'BioDockify.tests',
    'BioDockify.tests.*',
    'BioDockify.*.tests.*',
    'BioDockify.*.test_*',
    'BioDockify.setup',
    'BioDockify.setup.*',
    'BioDockify.*.setup_*',
    'BioDockify.build_*',
    'BioDockify.*.build_*',
    'BioDockify.conftest',
    'BioDockify.*.conftest'
]

binaries_excludes = [
    'BioDockify/src/ui/styles/colors.py',
]

# Build configuration
app = BUNDLE(
    'BioDockify',
    name='BioDockify',
    icon='src/ui/styles/icon.icns',
    bundle_identifier='com.biodockify.studio',
    info_plist={
        'NSHighResolutionCapable': 'True',
        'LSBackgroundOnly': 'False'
    }
)
