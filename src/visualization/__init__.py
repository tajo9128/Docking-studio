"""
BioDockify Visualization Module
Professional molecular visualization platform for Docking Studio.

Modules:
- advanced_viewer: Enhanced 3D molecular viewer with 3Dmol.js
- interaction_detector: Detection and visualization of molecular interactions
- surface_engine: Surface rendering (VDW, SAS, Electrostatic)
- pharmacophore_extractor: Pharmacophore feature extraction
- electrostatics_overlay: Electrostatic potential visualization
- ai_heatmap: AI confidence score overlay
- pose_comparator: Split-view pose comparison with RMSD
- export_manager: Publication-ready export (PNG, SVG, 4K)
- scene_animator: Pose animation and rotation
- download_manager: Download and export docking results
"""

# Import core modules (no GUI dependency)
from .interaction_detector import InteractionDetector, InteractionType
from .surface_engine import SurfaceEngine, SurfaceType
from .pharmacophore_extractor import PharmacophoreExtractor, PharmacophoreFeature
from .electrostatics_overlay import ElectrostaticsOverlay
from .ai_heatmap import AIHeatmapGenerator
from .pose_comparator import PoseComparator, RMSDCalculator
from .export_manager import ExportManager, ExportFormat
from .scene_animator import SceneAnimator, AnimationMode
from .download_manager import DownloadManager, FileExporter, ResultFile, DockingResult

# Try importing GUI modules (require PyQt6 WebEngine)
try:
    from .advanced_viewer import AdvancedMolecularViewer
    GUI_AVAILABLE = True
except ImportError:
    GUI_AVAILABLE = False
    AdvancedMolecularViewer = None

__all__ = [
    'AdvancedMolecularViewer',
    'GUI_AVAILABLE',
    'InteractionDetector',
    'InteractionType',
    'SurfaceEngine',
    'SurfaceType',
    'PharmacophoreExtractor',
    'PharmacophoreFeature',
    'ElectrostaticsOverlay',
    'AIHeatmapGenerator',
    'PoseComparator',
    'RMSDCalculator',
    'ExportManager',
    'ExportFormat',
    'SceneAnimator',
    'AnimationMode',
    'DownloadManager',
    'FileExporter',
    'ResultFile',
    'DockingResult',
]
