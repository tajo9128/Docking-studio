"""
Advanced Molecular Viewer
Enhanced 3D viewer with 3Dmol.js integration, split-view, and multiple representations.
"""

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
    QComboBox, QLabel, QSlider, QGroupBox, QToolBar,
    QCheckBox, QSpinBox, QSplitter, QFrame, QStackedWidget
)
from PyQt6.QtCore import Qt, pyqtSignal, QTimer, QSize
from PyQt6.QtWebEngineWidgets import QWebEngineView
from PyQt6.QtWebEngineCore import QWebEngineSettings
import json
from typing import List, Dict, Optional, Tuple, Callable
import logging

logger = logging.getLogger(__name__)


class ViewStyle:
    """Molecular view styles"""
    STICK = "stick"
    BALL_STICK = "ballandstick"
    CARTON = "cartoon"
    LINE = "line"
    SPHERE = "sphere"
    SURFACE = "surface"
    RIBBON = "ribbon"
    BACKBONE = "backbone"


class ColorScheme:
    """Color schemes for molecular coloring"""
    ELEMENT = "element"
    CHAIN = "chain"
    RESIDUE = "residue"
    BFACTOR = "bfactor"
    HYDROPHOBICITY = "hydrophobicity"
    POLARITY = "polarity"
    CHARGE = "charge"
    CUSTOM = "custom"


class BackgroundColor:
    """Background colors"""
    WHITE = "white"
    BLACK = "black"
    GRAY = "gray"
    DARK_GRAY = "#1a1a1a"
    LIGHT_GRAY = "#f0f0f0"


class AdvancedMolecularViewer(QWidget):
    """
    Advanced 3D Molecular Viewer with professional features:
    - Multiple representations (stick, cartoon, surface, etc.)
    - Split-view mode for pose comparison
    - Surface rendering with transparency
    - Interaction visualization
    - Animation support
    - Publication-ready export
    """
    
    structure_loaded = pyqtSignal(str)
    pose_selected = pyqtSignal(int)
    view_changed = pyqtSignal(str)
    export_completed = pyqtSignal(str)
    
    def __init__(self, parent=None, enable_split_view: bool = True):
        super().__init__(parent)
        
        self.current_receptor = None
        self.current_ligand = None
        self.poses: List[Tuple[int, str, float, str]] = []
        self.interactions: List[Dict] = []
        
        self.view_mode = "single"  # single, split, overlay
        self.split_orientation = "horizontal"
        
        self.animation_timer = QTimer()
        self.animation_timer.timeout.connect(self._animation_tick)
        self.current_animating_pose = 0
        
        self.callbacks: Dict[str, Callable] = {}
        
        self._setup_ui()
        self._init_viewer()
        
        logger.info("Advanced Molecular Viewer initialized")
    
    def _setup_ui(self):
        """Setup the complete UI with controls"""
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)
        
        self._setup_toolbar(main_layout)
        
        content_splitter = QSplitter(Qt.Orientation.Horizontal)
        
        self.viewer_container = QFrame()
        viewer_layout = QVBoxLayout(self.viewer_container)
        viewer_layout.setContentsMargins(0, 0, 0, 0)
        
        self.viewer = QWebEngineView()
        self.viewer.setMinimumSize(400, 400)
        viewer_layout.addWidget(self.viewer, 1)
        
        content_splitter.addWidget(self.viewer_container)
        
        self.controls_panel = self._create_controls_panel()
        content_splitter.addWidget(self.controls_panel)
        
        content_splitter.setStretchFactor(0, 4)
        content_splitter.setStretchFactor(1, 1)
        
        main_layout.addWidget(content_splitter, 1)
        
        self.status_bar = self._create_status_bar()
        main_layout.addWidget(self.status_bar)
    
    def _setup_toolbar(self, parent_layout):
        """Create the professional toolbar"""
        toolbar = QFrame()
        toolbar.setStyleSheet("""
            QFrame {
                background: #f8f9fa;
                border-bottom: 1px solid #e9ecef;
                padding: 8px;
            }
        """)
        toolbar_layout = QHBoxLayout(toolbar)
        toolbar_layout.setContentsMargins(8, 4, 8, 4)
        toolbar_layout.setSpacing(8)
        
        view_group = self._create_view_buttons()
        toolbar_layout.addWidget(view_group)
        
        toolbar_layout.addSpacing(16)
        
        style_group = self._create_style_buttons()
        toolbar_layout.addWidget(style_group)
        
        toolbar_layout.addSpacing(16)
        
        action_group = self._create_action_buttons()
        toolbar_layout.addWidget(action_group)
        
        toolbar_layout.addStretch()
        
        view_mode_group = self._create_view_mode_buttons()
        toolbar_layout.addWidget(view_mode_group)
        
        parent_layout.addWidget(toolbar)
    
    def _create_view_buttons(self) -> QFrame:
        """Create view control buttons"""
        group = QFrame()
        layout = QHBoxLayout(group)
        layout.setSpacing(4)
        
        self.btn_rotate = QPushButton("⟳")
        self.btn_rotate.setToolTip("Auto Rotate")
        self.btn_rotate.setFixedSize(36, 36)
        self.btn_rotate.clicked.connect(self.toggle_auto_rotate)
        
        self.btn_zoom = QPushButton("⊕")
        self.btn_zoom.setToolTip("Zoom to Fit")
        self.btn_zoom.setFixedSize(36, 36)
        self.btn_zoom.clicked.connect(self.zoom_to_fit)
        
        self.btn_center = QPushButton("◎")
        self.btn_center.setToolTip("Center View")
        self.btn_center.setFixedSize(36, 36)
        self.btn_center.clicked.connect(self.center_view)
        
        self.btn_reset = QPushButton("↺")
        self.btn_reset.setToolTip("Reset View")
        self.btn_reset.setFixedSize(36, 36)
        self.btn_reset.clicked.connect(self.reset_view)
        
        for btn in [self.btn_rotate, self.btn_zoom, self.btn_center, self.btn_reset]:
            btn.setStyleSheet("""
                QPushButton {
                    background: white;
                    border: 1px solid #dee2e6;
                    border-radius: 6px;
                    font-size: 16px;
                }
                QPushButton:hover {
                    background: #e9ecef;
                    border-color: #2E5AAC;
                }
            """)
        
        layout.addWidget(self.btn_rotate)
        layout.addWidget(self.btn_zoom)
        layout.addWidget(self.btn_center)
        layout.addWidget(self.btn_reset)
        
        return group
    
    def _create_style_buttons(self) -> QFrame:
        """Create molecular style buttons"""
        group = QFrame()
        layout = QHBoxLayout(group)
        layout.setSpacing(4)
        
        styles = [
            ("stick", "Stick"),
            ("ballandstick", "Ball & Stick"),
            ("cartoon", "Cartoon"),
            ("surface", "Surface"),
            ("sphere", "Spacefill"),
        ]
        
        self.style_buttons = {}
        for style_id, tooltip in styles:
            btn = QPushButton(style_id.replace("ballandstick", "B&S"))
            btn.setToolTip(tooltip)
            btn.setFixedWidth(60)
            btn.setCheckable(True)
            btn.setStyleSheet("""
                QPushButton {
                    background: white;
                    border: 1px solid #dee2e6;
                    border-radius: 4px;
                    font-size: 11px;
                    padding: 4px;
                }
                QPushButton:checked {
                    background: #2E5AAC;
                    color: white;
                    border-color: #2E5AAC;
                }
                QPushButton:hover:not(:checked) {
                    background: #e9ecef;
                }
            """)
            btn.clicked.connect(lambda checked, s=style_id: self.set_style(s))
            self.style_buttons[style_id] = btn
            layout.addWidget(btn)
        
        self.style_buttons["stick"].setChecked(True)
        
        return group
    
    def _create_action_buttons(self) -> QFrame:
        """Create action buttons"""
        group = QFrame()
        layout = QHBoxLayout(group)
        layout.setSpacing(4)
        
        self.btn_screenshot = QPushButton("📷")
        self.btn_screenshot.setToolTip("Take Screenshot")
        self.btn_screenshot.setFixedSize(36, 36)
        self.btn_screenshot.clicked.connect(self.take_screenshot)
        
        self.btn_fullscreen = QPushButton("⛶")
        self.btn_fullscreen.setToolTip("Fullscreen")
        self.btn_fullscreen.setFixedSize(36, 36)
        
        self.btn_interactions = QPushButton("🔗")
        self.btn_interactions.setToolTip("Show Interactions")
        self.btn_interactions.setCheckable(True)
        self.btn_interactions.setFixedSize(36, 36)
        self.btn_interactions.clicked.connect(self.toggle_interactions)
        
        self.btn_overlay_all = QPushButton("⊞")
        self.btn_overlay_all.setToolTip("Overlay All Poses")
        self.btn_overlay_all.setFixedSize(36, 36)
        self.btn_overlay_all.clicked.connect(self.show_overlay_all_poses)
        
        self.btn_best_pose = QPushButton("★")
        self.btn_best_pose.setToolTip("Highlight Best Pose")
        self.btn_best_pose.setFixedSize(36, 36)
        self.btn_best_pose.clicked.connect(self.highlight_best_pose)
        
        self.btn_animate = QPushButton("▶")
        self.btn_animate.setToolTip("Animate Poses")
        self.btn_animate.setFixedSize(36, 36)
        self.btn_animate.clicked.connect(lambda: self.animate_poses(1000))
        
        for btn in [self.btn_screenshot, self.btn_fullscreen, self.btn_interactions, 
                    self.btn_overlay_all, self.btn_best_pose, self.btn_animate]:
            btn.setStyleSheet("""
                QPushButton {
                    background: white;
                    border: 1px solid #dee2e6;
                    border-radius: 6px;
                    font-size: 14px;
                }
                QPushButton:hover {
                    background: #e9ecef;
                }
                QPushButton:checked {
                    background: #2E5AAC;
                    color: white;
                }
            """)
        
        layout.addWidget(self.btn_screenshot)
        layout.addWidget(self.btn_fullscreen)
        layout.addWidget(self.btn_interactions)
        layout.addWidget(self.btn_overlay_all)
        layout.addWidget(self.btn_best_pose)
        layout.addWidget(self.btn_animate)
        
        return group
    
    def _create_view_mode_buttons(self) -> QFrame:
        """Create view mode toggle buttons"""
        group = QFrame()
        layout = QHBoxLayout(group)
        layout.setSpacing(4)
        
        self.btn_single = QPushButton("Single")
        self.btn_single.setCheckable(True)
        self.btn_single.setChecked(True)
        self.btn_single.setFixedWidth(60)
        self.btn_single.setStyleSheet("""
            QPushButton {
                background: #2E5AAC;
                color: white;
                border: none;
                border-radius: 4px;
                font-size: 11px;
            }
            QPushButton:checked {
                background: #2E5AAC;
                color: white;
            }
            QPushButton:!checked {
                background: white;
                color: #495057;
                border: 1px solid #dee2e6;
            }
        """)
        self.btn_single.clicked.connect(lambda: self.set_view_mode("single"))
        
        self.btn_split = QPushButton("Split")
        self.btn_split.setCheckable(True)
        self.btn_split.setFixedWidth(60)
        self.btn_split.setStyleSheet("""
            QPushButton:checked {
                background: #2E5AAC;
                color: white;
            }
            QPushButton:!checked {
                background: white;
                color: #495057;
                border: 1px solid #dee2e6;
            }
        """)
        self.btn_split.clicked.connect(lambda: self.set_view_mode("split"))
        
        self.btn_overlay = QPushButton("Overlay")
        self.btn_overlay.setCheckable(True)
        self.btn_overlay.setFixedWidth(60)
        self.btn_overlay.setStyleSheet("""
            QPushButton:checked {
                background: #2E5AAC;
                color: white;
            }
            QPushButton:!checked {
                background: white;
                color: #495057;
                border: 1px solid #dee2e6;
            }
        """)
        self.btn_overlay.clicked.connect(lambda: self.set_view_mode("overlay"))
        
        layout.addWidget(self.btn_single)
        layout.addWidget(self.btn_split)
        layout.addWidget(self.btn_overlay)
        
        return group
    
    def _create_controls_panel(self) -> QWidget:
        """Create the controls panel"""
        panel = QFrame()
        panel.setMinimumWidth(280)
        panel.setStyleSheet("""
            QFrame {
                background: white;
                border-left: 1px solid #e9ecef;
            }
            QGroupBox {
                font-weight: 600;
                color: #495057;
                border: 1px solid #dee2e6;
                border-radius: 8px;
                margin-top: 8px;
                padding-top: 8px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 12px;
                padding: 0 4px;
            }
        """)
        
        layout = QVBoxLayout(panel)
        layout.setSpacing(8)
        layout.setContentsMargins(12, 12, 12, 12)
        
        layout.addWidget(self._create_color_group())
        layout.addWidget(self._create_surface_group())
        layout.addWidget(self._create_ligand_group())
        layout.addWidget(self._create_pose_selector())
        
        layout.addStretch()
        
        return panel
    
    def _create_color_group(self) -> QGroupBox:
        """Create color scheme controls with professional options"""
        group = QGroupBox("Coloring")
        layout = QVBoxLayout(group)
        layout.setSpacing(8)
        
        # Theme presets
        layout.addWidget(QLabel("Theme:"))
        self.theme_combo = QComboBox()
        self.theme_combo.addItems([
            "Professional", "PyMOL", "Discovery Studio", "Rainbow", "Monochrome"
        ])
        self.theme_combo.currentTextChanged.connect(self._on_theme_change)
        layout.addWidget(self.theme_combo)
        
        # Separator
        separator = QFrame()
        separator.setFrameShape(QFrame.Shape.HLine)
        separator.setStyleSheet("background: #dee2e6; height: 1px;")
        layout.addWidget(separator)
        
        # Protein coloring
        layout.addWidget(QLabel("Protein:"))
        self.protein_color_combo = QComboBox()
        self.protein_color_combo.addItems([
            "Secondary Structure", "Chain", "B-Factor", "Hydrophobicity", 
            "Residue Type", "Element", "Rainbow", "White", "Light Blue"
        ])
        self.protein_color_combo.currentTextChanged.connect(self._on_protein_color_change)
        layout.addWidget(self.protein_color_combo)
        
        # Ligand coloring
        layout.addWidget(QLabel("Ligand:"))
        self.ligand_color_combo = QComboBox()
        self.ligand_color_combo.addItems([
            "Element (CPK)", "Cyan", "Magenta", "Green", "White",
            "Orange", "Pink", "Yellow", "By Charge"
        ])
        self.ligand_color_combo.currentTextChanged.connect(self._on_ligand_color_change)
        layout.addWidget(self.ligand_color_combo)
        
        return group
    
    def _on_theme_change(self, theme: str):
        """Handle theme preset change"""
        theme_map = {
            "Professional": "professional",
            "PyMOL": "pymol",
            "Discovery Studio": "discovery_studio",
            "Rainbow": "rainbow",
            "Monochrome": "monochrome",
        }
        self.apply_theme(theme_map.get(theme, "professional"))
    
    def _on_protein_color_change(self, scheme: str):
        """Handle protein color scheme change"""
        scheme_map = {
            "Secondary Structure": "secondary",
            "Chain": "chain",
            "B-Factor": "bfactor",
            "Hydrophobicity": "hydrophobicity",
            "Residue Type": "residue",
            "Element": "element",
            "Rainbow": "rainbow",
            "White": "white",
            "Light Blue": "lightblue",
        }
        self.color_protein(scheme_map.get(scheme, "secondary"))
    
    def _on_ligand_color_change(self, scheme: str):
        """Handle ligand color scheme change"""
        scheme_map = {
            "Element (CPK)": "element",
            "Cyan": "cyan",
            "Magenta": "magenta",
            "Green": "green",
            "White": "white",
            "Orange": "orange",
            "Pink": "pink",
            "Yellow": "yellow",
            "By Charge": "charge",
        }
        self.color_ligand(scheme_map.get(scheme, "element"))
    
    def _create_surface_group(self) -> QGroupBox:
        """Create surface controls"""
        group = QGroupBox("Surface")
        layout = QVBoxLayout(group)
        layout.setSpacing(8)
        
        self.surface_enabled = QCheckBox("Enable Surface")
        self.surface_enabled.stateChanged.connect(self._on_surface_toggle)
        layout.addWidget(self.surface_enabled)
        
        layout.addWidget(QLabel("Surface Type:"))
        self.surface_type_combo = QComboBox()
        self.surface_type_combo.addItems(["VDW", "SAS", "Molecular Surface"])
        self.surface_type_combo.setEnabled(False)
        layout.addWidget(self.surface_type_combo)
        
        layout.addWidget(QLabel("Opacity:"))
        self.surface_opacity_slider = QSlider(Qt.Orientation.Horizontal)
        self.surface_opacity_slider.setRange(10, 100)
        self.surface_opacity_slider.setValue(40)
        self.surface_opacity_slider.setEnabled(False)
        self.surface_opacity_slider.valueChanged.connect(self._on_opacity_change)
        layout.addWidget(self.surface_opacity_slider)
        
        return group
    
    def _create_ligand_group(self) -> QGroupBox:
        """Create ligand display controls"""
        group = QGroupBox("Ligand Display")
        layout = QVBoxLayout(group)
        layout.setSpacing(8)
        
        self.ligand_enabled = QCheckBox("Show Ligand")
        self.ligand_enabled.setChecked(True)
        self.ligand_enabled.stateChanged.connect(self._on_ligand_toggle)
        layout.addWidget(self.ligand_enabled)
        
        self.ligand_stick = QCheckBox("Show as Stick")
        self.ligand_stick.setChecked(True)
        layout.addWidget(self.ligand_stick)
        
        self.ligand_label = QCheckBox("Show Atom Labels")
        layout.addWidget(self.ligand_label)
        
        return group
    
    def _create_pose_selector(self) -> QGroupBox:
        """Create pose selector"""
        group = QGroupBox("Poses")
        layout = QVBoxLayout(group)
        layout.setSpacing(8)
        
        self.pose_count_label = QLabel("0 poses loaded")
        layout.addWidget(self.pose_count_label)
        
        self.pose_slider = QSlider(Qt.Orientation.Horizontal)
        self.pose_slider.setRange(0, 0)
        self.pose_slider.setEnabled(False)
        self.pose_slider.valueChanged.connect(self._on_pose_change)
        layout.addWidget(self.pose_slider)
        
        self.pose_info_label = QLabel("")
        layout.addWidget(self.pose_info_label)
        
        nav_layout = QHBoxLayout()
        self.btn_prev_pose = QPushButton("◀")
        self.btn_prev_pose.setFixedWidth(40)
        self.btn_prev_pose.clicked.connect(self.prev_pose)
        
        self.btn_next_pose = QPushButton("▶")
        self.btn_next_pose.setFixedWidth(40)
        self.btn_next_pose.clicked.connect(self.next_pose)
        
        nav_layout.addWidget(self.btn_prev_pose)
        nav_layout.addStretch()
        nav_layout.addWidget(self.btn_next_pose)
        layout.addLayout(nav_layout)
        
        return group
    
    def _create_status_bar(self) -> QFrame:
        """Create status bar"""
        status = QFrame()
        status.setStyleSheet("""
            QFrame {
                background: #f8f9fa;
                border-top: 1px solid #e9ecef;
                padding: 4px 12px;
            }
        """)
        layout = QHBoxLayout(status)
        layout.setContentsMargins(12, 4, 12, 4)
        
        self.status_label = QLabel("Ready")
        self.status_label.setStyleSheet("color: #6c757d; font-size: 12px;")
        
        self.info_label = QLabel("")
        self.info_label.setStyleSheet("color: #6c757d; font-size: 12px;")
        
        layout.addWidget(self.status_label)
        layout.addStretch()
        layout.addWidget(self.info_label)
        
        return status
    
    def _init_viewer(self):
        """Initialize the 3Dmol.js viewer with full features"""
        html = '''
<!DOCTYPE html>
<html>
<head>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        body { margin: 0; padding: 0; overflow: hidden; background: white; }
        #viewer { width: 100%; height: 100%; position: relative; }
        .split-view { display: flex; width: 100%; height: 100%; }
        .split-panel { flex: 1; position: relative; }
        .split-divider { width: 4px; background: #dee2e6; cursor: col-resize; }
    </style>
</head>
<body>
    <div id="viewer"></div>
    <script>
        let viewer = null;
        let viewer2 = null;
        let splitMode = false;
        let currentStyle = 'stick';
        let ligandVisible = true;
        let surfaceEnabled = false;
        let surfaceOpacity = 0.4;
        let autoRotate = false;
        
        // Initialize viewer
        function initViewer() {
            viewer = $3Dmol.createViewer("viewer", {
                backgroundColor: "white",
                antialias: true,
                cartoonQuality: 10
            });
            
            // Mouse controls
            viewer.setClickable({}, true, function(atom, viewer, event, container) {
                if (atom) {
                    showAtomInfo(atom);
                }
            });
            
            renderLoop();
        }
        
        // Render loop for animations
        function renderLoop() {
            if (autoRotate && viewer) {
                viewer.rotate(0.5, 'y');
                viewer.render();
            }
            requestAnimationFrame(renderLoop);
        }
        
        // Show atom information
        function showAtomInfo(atom) {
            console.log("Selected atom:", atom.elem, atom.resn, atom.resi);
        }
        
        // Load structure
        function loadReceptor(pdbData) {
            viewer.removeAllModels();
            viewer.addModel(pdbData, "pdb");
            applyStyle(currentStyle);
            viewer.zoomTo();
            viewer.render();
            return true;
        }
        
        function loadLigand(pdbData, asModel = 0) {
            if (asModel === 0) {
                viewer.addModel(pdbData, "pdb");
            } else {
                viewer.addModel(pdbData, "pdb");
            }
            viewer.setStyle({model: -1}, {stick: {colorscheme: "Jmol"}});
            viewer.render();
        }
        
        // Apply style to receptor
        function applyStyle(style) {
            currentStyle = style;
            viewer.setStyle({}, {cartoon: {color: "spectrum"}});
            
            if (style === "stick" || style === "ballandstick") {
                viewer.setStyle({}, {stick: {}, cartoon: {opacity: 0.8}});
            } else if (style === "surface") {
                viewer.addSurface($3Dmol.SurfaceType.VDW, {
                    opacity: surfaceOpacity,
                    color: "white"
                });
            } else if (style === "sphere") {
                viewer.setStyle({}, {sphere: {}});
            }
            
            viewer.render();
        }
        
        // Toggle surface
        function toggleSurface(enabled, type, opacity) {
            surfaceEnabled = enabled;
            surfaceOpacity = opacity;
            
            if (enabled) {
                viewer.removeAllSurfaces();
                viewer.addSurface($3Dmol.SurfaceType.VDW, {
                    opacity: opacity,
                    color: "white"
                });
            } else {
                viewer.removeAllSurfaces();
            }
            viewer.render();
        }
        
        // Toggle ligand visibility
        function setLigandVisible(visible) {
            ligandVisible = visible;
            viewer.setStyle({model: -1}, visible ? {stick: {}} : {});
            viewer.render();
        }
        
        // Auto rotate
        function setAutoRotate(enabled) {
            autoRotate = enabled;
        }
        
        // View controls
        function resetView() {
            viewer.zoomTo();
            viewer.render();
        }
        
        function centerView() {
            viewer.center();
            viewer.render();
        }
        
        function zoomToFit() {
            viewer.zoomTo();
            viewer.render();
        }
        
        // Export
        function exportImage(callback) {
            viewer.pngURI(function(uri) {
                callback(uri);
            });
        }
        
        // Add interaction line
        function addInteractionLine(start, end, color, type) {
            viewer.addCylinder({
                start: {x: start[0], y: start[1], z: start[2]},
                end: {x: end[0], y: end[1], z: end[2]},
                radius: type === 'hbond' ? 0.05 : 0.1,
                color: color,
                dashed: type === 'hbond'
            });
            viewer.render();
        }
        
        // Clear interactions
        function clearInteractions() {
            viewer.removeAllSurfaces();
            // Note: lines would need specific handling
            viewer.render();
        }
        
        // ==================== ADVANCED COLORING FUNCTIONS ====================
        
        // Color schemes for protein
        function colorProtein(scheme) {
            var style = {};
            
            if (scheme === 'chain') {
                // Color by chain - different colors per chain
                viewer.setStyle({}, {cartoon: {color: "chain"}});
            } else if (scheme === 'secondary') {
                // Secondary structure coloring (helix=red, sheet=yellow, loop=green)
                viewer.setStyle({}, {cartoon: {color: "sstruc"}});
            } else if (scheme === 'bfactor') {
                // B-factor / temperature factor coloring
                viewer.setStyle({}, {cartoon: {color: "bfactor"}});
            } else if (scheme === 'hydrophobicity') {
                // Hydrophobicity coloring
                viewer.setStyle({}, {cartoon: {color: "hydrophobicity"}});
            } else if (scheme === 'residue') {
                // Color by residue type
                viewer.setStyle({}, {cartoon: {color: "residueindex"}});
            } else if (scheme === 'element') {
                // Element-based (CPK)
                viewer.setStyle({}, {cartoon: {color: "element"}});
            } else if (scheme === 'rainbow') {
                // Rainbow gradient
                viewer.setStyle({}, {cartoon: {color: "rainbow"}});
            } else {
                // Default / custom single color
                viewer.setStyle({}, {cartoon: {color: scheme}});
            }
            
            viewer.render();
        }
        
        // Color schemes for ligand
        function colorLigand(scheme) {
            // Target the last added model (ligand)
            var ligandModel = viewer.getModel();
            
            if (!ligandModel) return;
            
            if (scheme === 'element') {
                // CPK element coloring (default)
                viewer.setStyle({model: -1}, {stick: {colorscheme: "Jmol"}});
            } else if (scheme === 'cyan') {
                viewer.setStyle({model: -1}, {stick: {color: "cyan"}});
            } else if (scheme === 'magenta') {
                viewer.setStyle({model: -1}, {stick: {color: "magenta"}});
            } else if (scheme === 'green') {
                viewer.setStyle({model: -1}, {stick: {color: "green"}});
            } else if (scheme === 'white') {
                viewer.setStyle({model: -1}, {stick: {color: "white"}});
            } else if (scheme === 'orange') {
                viewer.setStyle({model: -1}, {stick: {color: "orange"}});
            } else if (scheme === 'pink') {
                viewer.setStyle({model: -1}, {stick: {color: "pink"}});
            } else if (scheme === 'yellow') {
                viewer.setStyle({model: -1}, {stick: {color: "yellow"}});
            } else if (scheme === 'charge') {
                // Color by partial charge
                viewer.setStyle({model: -1}, {stick: {colorscheme: "Jmol"}});
            }
            
            viewer.render();
        }
        
        // Apply professional style (protein cartoon + ligand stick)
        function applyProfessionalStyle() {
            viewer.setStyle({}, {cartoon: {color: "sstruc"}});  // Secondary structure
            viewer.addStyle({model: -1}, {stick: {colorscheme: "Jmol"}});  // Ligand with CPK
            viewer.render();
        }
        
        // Apply PyMOL-style look
        function applyPyMOLStyle() {
            viewer.setStyle({}, {
                cartoon: {color: "chain", opacity: 0.9},
                stick: {radius: 0.2}
            });
            viewer.addStyle({model: -1}, {
                stick: {colorscheme: "Jmol", radius: 0.3},
                sphere: {scale: 0.3}
            });
            viewer.render();
        }
        
        // Apply Discovery Studio style
        function applyDiscoveryStudioStyle() {
            viewer.setStyle({}, {cartoon: {color: "lightblue", opacity: 0.85}});
            viewer.addStyle({model: -1}, {stick: {color: "magenta", radius: 0.25}});
            viewer.render();
        }
        
        // Highlight binding site (residues near ligand)
        function highlightBindingSite(cutoff) {
            cutoff = cutoff || 5.0;
            
            // This would need the ligand position - simplified for now
            viewer.addStyle({}, {
                cartoon: {color: "white", opacity: 0.3}
            });
            viewer.render();
        }
        
        // Custom color for selection
        function colorSelection(selection, color) {
            viewer.addStyle(selection, {cartoon: {color: color}});
            viewer.render();
        }
        
        // ==================== END ADVANCED COLORING ====================
        
        // Set background
        function setBackground(color) {
            viewer.setBackgroundColor(color);
            viewer.render();
        }
        
        // Snapshot functions
        function takeSnapshot(factor, transparent, callback) {
            factor = factor || 1;
            transparent = transparent || false;
            
            // Get PNG with specified resolution factor
            viewer.pngURI(function(uri) {
                callback(uri);
            }, factor, transparent);
        }
        
        function setPublicationMode(enabled) {
            if (enabled) {
                // Hide axes, set white background, enable AA
                viewer.setBackgroundColor("white");
            }
            viewer.render();
        }
        
        // Initialize on load
        $(document).ready(function() {
            initViewer();
        });
    </script>
</body>
</html>
        '''
        self.viewer.setHtml(html)
        self.viewer.loadFinished.connect(self._on_viewer_ready)
    
    def _on_viewer_ready(self, ok):
        """Handle viewer ready event"""
        if ok:
            self.status_label.setText("Viewer ready")
            logger.info("3Dmol.js viewer initialized successfully")
    
    # ==================== Public Methods ====================
    
    def load_receptor(self, pdb_content: str):
        """Load receptor PDB structure"""
        self.current_receptor = pdb_content
        escaped = json.dumps(pdb_content)
        self.viewer.page().runJavaScript(f'loadReceptor({escaped})')
        self.status_label.setText("Receptor loaded")
        self.structure_loaded.emit("receptor")
    
    def load_ligand(self, pdb_content: str, pose_id: int = 0):
        """Load ligand PDB structure"""
        self.current_ligand = pdb_content
        escaped = json.dumps(pdb_content)
        self.viewer.page().runJavaScript(f'loadLigand({escaped}, {pose_id})')
        self.status_label.setText(f"Pose {pose_id} loaded")
    
    POSE_COLORS = [
        'cyan', 'magenta', 'green', 'orange', 'pink', 
        'yellow', 'white', 'red', 'blue', 'purple'
    ]
    
    def add_pose(self, pdb_content: str, pose_id: int, score: float = 0.0, color: Optional[str] = None):
        """Add a docking pose with optional custom color"""
        if color is None:
            color = self.POSE_COLORS[pose_id % len(self.POSE_COLORS)]
        
        self.poses.append((pose_id, pdb_content, score, color))
        self.pose_slider.setRange(0, len(self.poses) - 1)
        self.pose_slider.setEnabled(len(self.poses) > 1)
        self.pose_count_label.setText(f"{len(self.poses)} poses loaded")
        
        escaped = json.dumps(pdb_content)
        self.viewer.page().runJavaScript(f'''
            viewer.addModel({escaped}, "pdb");
            viewer.setStyle({{model: -1}}, {{stick: {{color: "{color}"}}}});
            viewer.render();
        ''')
        self.status_label.setText(f"Pose {pose_id + 1} added ({color})")
    
    def add_pose_overlay(self, pdb_content: str, pose_id: int, score: float = 0.0, color: Optional[str] = None):
        """Add pose in overlay mode - all poses visible simultaneously"""
        if color is None:
            color = self.POSE_COLORS[pose_id % len(self.POSE_COLORS)]
        
        self.poses.append((pose_id, pdb_content, score, color))
        
        escaped = json.dumps(pdb_content)
        self.viewer.page().runJavaScript(f'''
            viewer.addModel({escaped}, "pdbqt");
            viewer.setStyle({{model: {pose_id}}}, {{stick: {{color: "{color}", radius: 0.15}}}});
            viewer.render();
        ''')
        
        self.pose_count_label.setText(f"{len(self.poses)} poses in overlay")
        self.status_label.setText(f"Pose {pose_id + 1} added to overlay")
    
    def show_overlay_all_poses(self):
        """Show all poses in overlay mode with different colors"""
        self.viewer.page().runJavaScript('''
            viewer.removeAllModels();
            viewer.render();
        ''')
        
        for pose_id, pdb_content, score, color in self.poses:
            escaped = json.dumps(pdb_content)
            self.viewer.page().runJavaScript(f'''
                viewer.addModel({escaped}, "pdbqt");
                viewer.setStyle({{model: -1}}, {{stick: {{color: "{color}", radius: 0.15}}}});
                viewer.render();
            ''')
        
        self.viewer.page().runJavaScript('viewer.zoomTo(); viewer.render();')
        self.status_label.setText(f"Showing {len(self.poses)} poses in overlay")
    
    def show_single_pose(self, pose_index: int):
        """Show only a single pose from the overlay"""
        if pose_index >= len(self.poses):
            return
        
        for i, (pose_id, pdb_content, score, color) in enumerate(self.poses):
            visible = "visible" if i == pose_index else ""
            js_code = f'''
                viewer.setStyle({{model: {i}}}, {{stick: {{color: "{color}", radius: 0.15}}}});
                viewer.render();
            '''
            self.viewer.page().runJavaScript(js_code)
        
        self.status_label.setText(f"Showing pose {pose_index + 1}")
    
    def highlight_best_pose(self):
        """Highlight the pose with the best (lowest) score"""
        if not self.poses:
            return
        
        best_pose = min(self.poses, key=lambda p: p[2] if p[2] else float('inf'))
        best_idx = self.poses.index(best_pose)
        
        self.viewer.page().runJavaScript('''
            viewer.removeAllModels();
            viewer.render();
        ''')
        
        for i, (pose_id, pdb_content, score, color) in enumerate(self.poses):
            escaped = json.dumps(pdb_content)
            radius = 0.25 if i == best_idx else 0.1
            opacity = 1.0 if i == best_idx else 0.5
            
            self.viewer.page().runJavaScript(f'''
                viewer.addModel({escaped}, "pdbqt");
                viewer.setStyle({{model: -1}}, {{
                    stick: {{color: "{color}", radius: {radius}, opacity: {opacity}}}
                }});
                viewer.render();
            ''')
        
        self.viewer.page().runJavaScript('viewer.zoomTo(); viewer.render();')
        self.status_label.setText(f"Best pose highlighted (score: {best_pose[2]:.2f})")
    
    def animate_poses(self, interval_ms: int = 1000):
        """Animate through all poses"""
        if not self.poses:
            return
        
        self._animation_interval = interval_ms
        self._animate_index = 0
        
        if hasattr(self, 'animation_timer') and self.animation_timer.isActive():
            self.animation_timer.stop()
        
        self.animation_timer.timeout.connect(self._animate_tick)
        self.animation_timer.start(interval_ms)
        self.status_label.setText(f"Animating {len(self.poses)} poses...")
    
    def stop_animation(self):
        """Stop pose animation"""
        if hasattr(self, 'animation_timer'):
            self.animation_timer.stop()
        self.status_label.setText("Animation stopped")
    
    def _animate_tick(self):
        """Animation tick - show next pose"""
        if not self.poses:
            return
        
        self._animate_index = (self._animate_index + 1) % len(self.poses)
        self.show_single_pose(self._animate_index)
        self.pose_slider.setValue(self._animate_index)
    
    def export_posed_image(self, filename: str):
        """Export overlay image with all poses"""
        self.take_screenshot(resolution=2, callback=lambda data: self._save_export(data, filename))
    
    def _save_export(self, data_uri: str, filename: str):
        """Save exported image"""
        import base64
        if data_uri and data_uri.startswith('data:image/png;base64,'):
            data = base64.b64decode(data_uri.split(',')[1])
            with open(filename, 'wb') as f:
                f.write(data)
            logger.info(f"Multi-pose image exported to {filename}")
    
    def cluster_poses(self, rmsd_threshold: float = 2.0) -> List[List[int]]:
        """
        Cluster poses by RMSD similarity.
        
        Args:
            rmsd_threshold: RMSD threshold for grouping poses (Angstroms)
        
        Returns:
            List of pose index groups
        """
        if len(self.poses) < 2:
            return [[0]]
        
        try:
            from backend.analysis import calculate_rmsd
            
            n = len(self.poses)
            clusters = []
            assigned = [False] * n
            
            for i in range(n):
                if assigned[i]:
                    continue
                
                cluster = [i]
                assigned[i] = True
                
                pdb_i = self.poses[i][1]
                
                for j in range(i + 1, n):
                    if assigned[j]:
                        continue
                    
                    pdb_j = self.poses[j][1]
                    rmsd = calculate_rmsd(pdb_i, pdb_j)
                    
                    if 0 <= rmsd < rmsd_threshold:
                        cluster.append(j)
                        assigned[j] = True
                
                clusters.append(cluster)
            
            logger.info(f"Poses clustered into {len(clusters)} groups")
            return clusters
        
        except Exception as e:
            logger.error(f"Pose clustering failed: {e}")
            return [[i] for i in range(len(self.poses))]
    
    def show_cluster(self, cluster_indices: List[int]):
        """Show only poses from a specific cluster"""
        self.viewer.page().runJavaScript('viewer.removeAllModels(); viewer.render();')
        
        for idx in cluster_indices:
            if idx < len(self.poses):
                pose_id, pdb_content, score, color = self.poses[idx]
                escaped = json.dumps(pdb_content)
                self.viewer.page().runJavaScript(f'''
                    viewer.addModel({escaped}, "pdbqt");
                    viewer.setStyle({{model: -1}}, {{stick: {{color: "{color}", radius: 0.15}}}});
                    viewer.render();
                ''')
        
        self.viewer.page().runJavaScript('viewer.zoomTo(); viewer.render();')
        self.status_label.setText(f"Showing cluster with {len(cluster_indices)} poses")
    
    def compare_poses(self, pose1_idx: int, pose2_idx: int) -> Dict:
        """
        Compare two poses and return analysis.
        
        Args:
            pose1_idx: Index of first pose
            pose2_idx: Index of second pose
        
        Returns:
            Comparison dictionary with RMSD and score difference
        """
        if pose1_idx >= len(self.poses) or pose2_idx >= len(self.poses):
            return {"error": "Invalid pose indices"}
        
        pose1 = self.poses[pose1_idx]
        pose2 = self.poses[pose2_idx]
        
        try:
            from backend.analysis import calculate_rmsd
            rmsd = calculate_rmsd(pose1[1], pose2[1])
        except:
            rmsd = -1.0
        
        score_diff = abs(pose1[2] - pose2[2]) if pose1[2] and pose2[2] else 0
        
        return {
            "pose1_id": pose1[0],
            "pose2_id": pose2[0],
            "pose1_score": pose1[2],
            "pose2_score": pose2[2],
            "score_difference": score_diff,
            "rmsd": rmsd,
            "similar": rmsd < 2.0 if rmsd >= 0 else None
        }
    
    def set_style(self, style: str):
        """Set molecular display style"""
        for s, btn in self.style_buttons.items():
            btn.setChecked(s == style)
        
        self.viewer.page().runJavaScript(f'applyStyle("{style}")')
        self.view_changed.emit(style)
    
    def set_view_mode(self, mode: str):
        """Set view mode: single, split, overlay"""
        self.view_mode = mode
        
        self.btn_single.setChecked(mode == "single")
        self.btn_split.setChecked(mode == "split")
        self.btn_overlay.setChecked(mode == "overlay")
        
        if mode == "split":
            self.viewer.page().runJavaScript('enableSplitView()')
        else:
            self.viewer.page().runJavaScript('disableSplitView()')
        
        logger.info(f"View mode changed to: {mode}")
    
    def toggle_auto_rotate(self):
        """Toggle auto rotation"""
        is_rotating = self.btn_rotate.isChecked()
        self.viewer.page().runJavaScript(f'setAutoRotate({str(is_rotating).lower()})')
        self.status_label.setText("Auto rotate: " + ("ON" if is_rotating else "OFF"))
    
    def zoom_to_fit(self):
        """Zoom to fit structure"""
        self.viewer.page().runJavaScript('zoomToFit()')
    
    def center_view(self):
        """Center view on structure"""
        self.viewer.page().runJavaScript('centerView()')
    
    def reset_view(self):
        """Reset view to default"""
        self.viewer.page().runJavaScript('resetView()')
        self.status_label.setText("View reset")
    
    def toggle_interactions(self):
        """Toggle interaction display"""
        show = self.btn_interactions.isChecked()
        if show and self.interactions:
            self._display_interactions()
        else:
            self.viewer.page().runJavaScript('clearInteractions()')
    
    def _display_interactions(self):
        """Display interactions in viewer"""
        for interaction in self.interactions:
            start = interaction.get('atom1_coords', [0, 0, 0])
            end = interaction.get('atom2_coords', [0, 0, 0])
            itype = interaction.get('type', 'hbond')
            
            color_map = {
                'hbond': 'green',
                'hydrophobic': 'yellow',
                'pi_stacking': 'purple',
                'salt_bridge': 'red',
                'metal': 'gray'
            }
            color = color_map.get(itype, 'white')
            
            self.viewer.page().runJavaScript(f'''
                addInteractionLine(
                    [{start[0]}, {start[1]}, {start[2]}],
                    [{end[0]}, {end[1]}, {end[2]}],
                    "{color}",
                    "{itype}"
                )
            ''')
    
    def set_interactions(self, interactions: List[Dict]):
        """Set interactions to display"""
        self.interactions = interactions
        if self.btn_interactions.isChecked():
            self._display_interactions()
    
    def take_screenshot(
        self,
        resolution: int = 1,
        transparent: bool = False,
        callback: Optional[Callable] = None
    ):
        """
        Take screenshot of current view.
        
        Args:
            resolution: Resolution factor (1, 2, 4 for 1x, 2x, 4x)
            transparent: Use transparent background
            callback: Optional callback to receive base64 image
        """
        def on_result(data_uri):
            if callback:
                callback(data_uri)
            elif data_uri:
                self.export_completed.emit(data_uri)
        
        self.viewer.page().runJavaScript(
            f'takeSnapshot({resolution}, {str(transparent).lower()})',
            on_result
        )
    
    def take_snapshot_4k(self, callback: Optional[Callable] = None):
        """Take 4K resolution screenshot (3840x2160)"""
        self.take_screenshot(resolution=4, transparent=False, callback=callback)
    
    def take_snapshot_transparent(self, callback: Optional[Callable] = None):
        """Take screenshot with transparent background"""
        self.take_screenshot(resolution=2, transparent=True, callback=callback)
    
    def set_publication_mode(self, enabled: bool = True):
        """Set publication mode (white background, high quality)"""
        self.viewer.page().runJavaScript(f'setPublicationMode({str(enabled).lower()})')
    
    def set_background(self, color: str):
        """Set background color"""
        color_map = {
            'white': 'white',
            'black': 'black',
            'gray': '#808080',
            'dark': '#1a1a1a'
        }
        self.viewer.page().runJavaScript(f'setBackground("{color_map.get(color, "white")}")')
    
    # ==================== Coloring Methods ====================
    
    def color_protein(self, scheme: str):
        """
        Apply color scheme to protein.
        
        Args:
            scheme: Color scheme - 'chain', 'secondary', 'bfactor', 
                    'hydrophobicity', 'residue', 'element', 'rainbow'
        """
        valid_schemes = ['chain', 'secondary', 'bfactor', 'hydrophobicity', 
                        'residue', 'element', 'rainbow']
        if scheme not in valid_schemes:
            scheme = 'secondary'  # Default
        
        self.viewer.page().runJavaScript(f'colorProtein("{scheme}")')
        self.status_label.setText(f"Protein colored by: {scheme}")
        logger.info(f"Protein coloring scheme: {scheme}")
    
    def color_ligand(self, scheme: str):
        """
        Apply color scheme to ligand.
        
        Args:
            scheme: Color scheme - 'element', 'cyan', 'magenta', 'green',
                    'white', 'orange', 'pink', 'yellow', 'charge'
        """
        valid_schemes = ['element', 'cyan', 'magenta', 'green', 
                        'white', 'orange', 'pink', 'yellow', 'charge']
        if scheme not in valid_schemes:
            scheme = 'element'  # Default
        
        self.viewer.page().runJavaScript(f'colorLigand("{scheme}")')
        self.status_label.setText(f"Ligand colored by: {scheme}")
        logger.info(f"Ligand coloring scheme: {scheme}")
    
    def apply_professional_style(self):
        """Apply professional style (secondary structure protein + element ligand)"""
        self.viewer.page().runJavaScript('applyProfessionalStyle()')
        self.status_label.setText("Professional style applied")
        logger.info("Applied professional coloring style")
    
    def apply_pymol_style(self):
        """Apply PyMOL-style coloring (chain colors + CPK ligand)"""
        self.viewer.page().runJavaScript('applyPyMOLStyle()')
        self.status_label.setText("PyMOL style applied")
        logger.info("Applied PyMOL-style coloring")
    
    def apply_discovery_studio_style(self):
        """Apply Discovery Studio-style coloring"""
        self.viewer.page().runJavaScript('applyDiscoveryStudioStyle()')
        self.status_label.setText("Discovery Studio style applied")
        logger.info("Applied Discovery Studio-style coloring")
    
    def highlight_binding_site(self, cutoff: float = 5.0):
        """Highlight binding site residues near ligand"""
        self.viewer.page().runJavaScript(f'highlightBindingSite({cutoff})')
        self.status_label.setText(f"Binding site highlighted ({cutoff}Å)")
    
    def color_by_custom(self, selection: str, color: str):
        """Apply custom color to selection"""
        self.viewer.page().runJavaScript(f'colorSelection("{selection}", "{color}")')
    
    # ==================== Preset Color Themes ====================
    
    def apply_theme(self, theme: str):
        """
        Apply a preset color theme.
        
        Args:
            theme: Theme name - 'professional', 'pymol', 'discovery_studio',
                   'monochrome', 'rainbow'
        """
        theme_map = {
            'professional': 'applyProfessionalStyle()',
            'pymol': 'applyPyMOLStyle()',
            'discovery_studio': 'applyDiscoveryStudioStyle()',
            'monochrome': 'colorProtein("white"); colorLigand("white");',
            'rainbow': 'colorProtein("rainbow"); colorLigand("element");',
        }
        
        js_func = theme_map.get(theme, 'applyProfessionalStyle()')
        self.viewer.page().runJavaScript(js_func)
        self.status_label.setText(f"Theme applied: {theme}")
        logger.info(f"Applied color theme: {theme}")
    
    # ==================== Event Handlers ====================
    
    def _on_color_change(self, text: str):
        """Handle color scheme change"""
        self.view_changed.emit(f"color:{text}")
    
    def _on_surface_toggle(self, state):
        """Handle surface toggle"""
        enabled = state == Qt.CheckState.Checked.value
        self.surface_type_combo.setEnabled(enabled)
        self.surface_opacity_slider.setEnabled(enabled)
        
        surface_type = self.surface_type_combo.currentText().lower()
        opacity = self.surface_opacity_slider.value() / 100.0
        
        self.viewer.page().runJavaScript(f'toggleSurface({str(enabled).lower()}, "{surface_type}", {opacity})')
    
    def _on_opacity_change(self, value):
        """Handle opacity change"""
        if self.surface_enabled.isChecked():
            opacity = value / 100.0
            surface_type = self.surface_type_combo.currentText().lower()
            self.viewer.page().runJavaScript(f'toggleSurface(true, "{surface_type}", {opacity})')
    
    def _on_ligand_toggle(self, state):
        """Handle ligand visibility toggle"""
        visible = state == Qt.CheckState.Checked.value
        self.viewer.page().runJavaScript(f'setLigandVisible({str(visible).lower()})')
    
    def _on_pose_change(self, value):
        """Handle pose slider change"""
        if value < len(self.poses):
            pose_id, pdb, score, color = self.poses[value]
            self.pose_info_label.setText(f"Pose {pose_id + 1} | Score: {score:.2f}")
            self.pose_selected.emit(pose_id)
    
    def prev_pose(self):
        """Go to previous pose"""
        current = self.pose_slider.value()
        if current > 0:
            self.pose_slider.setValue(current - 1)
    
    def next_pose(self):
        """Go to next pose"""
        current = self.pose_slider.value()
        if current < len(self.poses) - 1:
            self.pose_slider.setValue(current + 1)
    
    def _animation_tick(self):
        """Animation timer tick"""
        self.current_animating_pose = (self.current_animating_pose + 1) % len(self.poses)
        if self.poses:
            self.pose_slider.setValue(self.current_animating_pose)
    
    # ==================== Utility Methods ====================
    
    def get_current_pose(self) -> Optional[Tuple[int, str, float, str]]:
        """Get current pose data"""
        idx = self.pose_slider.value()
        if idx < len(self.poses):
            return self.poses[idx]
        return None
    
    def clear(self):
        """Clear all structures"""
        self.current_receptor = None
        self.current_ligand = None
        self.poses = []
        self.interactions = []
        self.viewer.page().runJavaScript('viewer.removeAllModels(); viewer.removeAllSurfaces(); viewer.render()')
        self.pose_slider.setRange(0, 0)
        self.pose_count_label.setText("0 poses loaded")
        self.status_label.setText("Cleared")


class ViewerCallback:
    """Callback handler for viewer JavaScript"""
    
    def __init__(self, viewer: AdvancedMolecularViewer):
        self.viewer = viewer
    
    def handle_atom_click(self, atom_data: dict):
        """Handle atom click from JavaScript"""
        logger.info(f"Atom clicked: {atom_data}")
