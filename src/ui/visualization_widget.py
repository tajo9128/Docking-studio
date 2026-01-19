"""
BioDockify Visualization Widget
3D Molecular Viewer using PyQt6 and py3Dmol/WebEngine
"""

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, 
    QComboBox, QLabel, QSlider, QGroupBox, QToolBar
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWebEngineWidgets import QWebEngineView
import json


class MolecularViewer3D(QWidget):
    """
    3D Molecular Viewer Widget with interactive controls.
    Uses WebEngine with 3Dmol.js for rendering.
    """
    
    structure_loaded = pyqtSignal(str)
    pose_selected = pyqtSignal(int)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.current_structure = None
        self.poses = []
        self.setup_ui()
        
    def setup_ui(self):
        """Initialize the UI components."""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Toolbar
        toolbar = self.create_toolbar()
        layout.addWidget(toolbar)
        
        # 3D Viewer (WebEngine)
        self.viewer = QWebEngineView()
        self.viewer.setMinimumSize(400, 400)
        layout.addWidget(self.viewer, 1)
        
        # Controls panel
        controls = self.create_controls()
        layout.addWidget(controls)
        
        # Initialize 3Dmol viewer
        self.init_viewer()
        
    def create_toolbar(self):
        """Create the toolbar with view controls."""
        toolbar = QToolBar()
        toolbar.setMovable(False)
        
        # View mode buttons
        self.btn_ball_stick = QPushButton("Ball & Stick")
        self.btn_ball_stick.clicked.connect(lambda: self.set_style("stick"))
        toolbar.addWidget(self.btn_ball_stick)
        
        self.btn_cartoon = QPushButton("Cartoon")
        self.btn_cartoon.clicked.connect(lambda: self.set_style("cartoon"))
        toolbar.addWidget(self.btn_cartoon)
        
        self.btn_surface = QPushButton("Surface")
        self.btn_surface.clicked.connect(lambda: self.set_style("surface"))
        toolbar.addWidget(self.btn_surface)
        
        self.btn_spacefill = QPushButton("Spacefill")
        self.btn_spacefill.clicked.connect(lambda: self.set_style("sphere"))
        toolbar.addWidget(self.btn_spacefill)
        
        toolbar.addSeparator()
        
        # Reset view
        self.btn_reset = QPushButton("Reset View")
        self.btn_reset.clicked.connect(self.reset_view)
        toolbar.addWidget(self.btn_reset)
        
        # Center
        self.btn_center = QPushButton("Center")
        self.btn_center.clicked.connect(self.center_view)
        toolbar.addWidget(self.btn_center)
        
        return toolbar
    
    def create_controls(self):
        """Create the control panel."""
        group = QGroupBox("Display Options")
        layout = QHBoxLayout(group)
        
        # Color scheme
        layout.addWidget(QLabel("Color:"))
        self.color_combo = QComboBox()
        self.color_combo.addItems([
            "Element", "Chain", "Residue", "B-Factor", "Hydrophobicity"
        ])
        self.color_combo.currentTextChanged.connect(self.change_color_scheme)
        layout.addWidget(self.color_combo)
        
        # Background
        layout.addWidget(QLabel("Background:"))
        self.bg_combo = QComboBox()
        self.bg_combo.addItems(["White", "Black", "Gray"])
        self.bg_combo.currentTextChanged.connect(self.change_background)
        layout.addWidget(self.bg_combo)
        
        # Opacity slider
        layout.addWidget(QLabel("Opacity:"))
        self.opacity_slider = QSlider(Qt.Orientation.Horizontal)
        self.opacity_slider.setRange(0, 100)
        self.opacity_slider.setValue(100)
        self.opacity_slider.valueChanged.connect(self.change_opacity)
        layout.addWidget(self.opacity_slider)
        
        layout.addStretch()
        
        return group
    
    def init_viewer(self):
        """Initialize the 3Dmol.js viewer."""
        html = '''
        <!DOCTYPE html>
        <html>
        <head>
            <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
            <style>
                body { margin: 0; padding: 0; overflow: hidden; }
                #viewer { width: 100%; height: 100vh; }
            </style>
        </head>
        <body>
            <div id="viewer"></div>
            <script>
                var viewer = $3Dmol.createViewer("viewer", {
                    backgroundColor: "white"
                });
                
                function loadStructure(pdb, format) {
                    viewer.removeAllModels();
                    viewer.addModel(pdb, format);
                    viewer.setStyle({}, {stick: {}});
                    viewer.zoomTo();
                    viewer.render();
                }
                
                function setStyle(style) {
                    var styleObj = {};
                    styleObj[style] = {};
                    viewer.setStyle({}, styleObj);
                    viewer.render();
                }
                
                function setBackground(color) {
                    viewer.setBackgroundColor(color);
                    viewer.render();
                }
                
                function resetView() {
                    viewer.zoomTo();
                    viewer.render();
                }
                
                function addSurface(opacity) {
                    viewer.addSurface($3Dmol.SurfaceType.VDW, {
                        opacity: opacity,
                        color: "white"
                    });
                    viewer.render();
                }
            </script>
        </body>
        </html>
        '''
        self.viewer.setHtml(html)
    
    def load_pdb(self, pdb_content: str):
        """Load a PDB structure into the viewer."""
        self.current_structure = pdb_content
        escaped = json.dumps(pdb_content)
        self.viewer.page().runJavaScript(f'loadStructure({escaped}, "pdb")')
        self.structure_loaded.emit("pdb")
    
    def load_mol2(self, mol2_content: str):
        """Load a MOL2 structure into the viewer."""
        self.current_structure = mol2_content
        escaped = json.dumps(mol2_content)
        self.viewer.page().runJavaScript(f'loadStructure({escaped}, "mol2")')
        self.structure_loaded.emit("mol2")
    
    def load_sdf(self, sdf_content: str):
        """Load an SDF structure into the viewer."""
        self.current_structure = sdf_content
        escaped = json.dumps(sdf_content)
        self.viewer.page().runJavaScript(f'loadStructure({escaped}, "sdf")')
        self.structure_loaded.emit("sdf")
    
    def set_style(self, style: str):
        """Set the display style."""
        self.viewer.page().runJavaScript(f'setStyle("{style}")')
    
    def change_color_scheme(self, scheme: str):
        """Change the color scheme."""
        scheme_map = {
            "Element": "default",
            "Chain": "chain",
            "Residue": "residue",
            "B-Factor": "b",
            "Hydrophobicity": "hydrophobicity"
        }
        color = scheme_map.get(scheme, "default")
        self.viewer.page().runJavaScript(
            f'viewer.setStyle({{}}, {{stick: {{colorscheme: "{color}"}}}});viewer.render()'
        )
    
    def change_background(self, color: str):
        """Change the background color."""
        color_map = {"White": "white", "Black": "black", "Gray": "#808080"}
        self.viewer.page().runJavaScript(f'setBackground("{color_map.get(color, "white")}")')
    
    def change_opacity(self, value: int):
        """Change surface opacity."""
        opacity = value / 100.0
        self.viewer.page().runJavaScript(f'addSurface({opacity})')
    
    def reset_view(self):
        """Reset the view to default."""
        self.viewer.page().runJavaScript('resetView()')
    
    def center_view(self):
        """Center the view on the structure."""
        self.viewer.page().runJavaScript('viewer.zoomTo();viewer.render()')
    
    def add_pose(self, pose_pdb: str, pose_id: int):
        """Add a docking pose to the viewer."""
        self.poses.append((pose_id, pose_pdb))
        escaped = json.dumps(pose_pdb)
        self.viewer.page().runJavaScript(f'''
            viewer.addModel({escaped}, "pdb");
            viewer.setStyle({{model: -1}}, {{stick: {{color: "cyan"}}}});
            viewer.render();
        ''')
    
    def highlight_interactions(self, interactions: list):
        """Highlight molecular interactions."""
        for interaction in interactions:
            if interaction['type'] == 'hbond':
                color = 'blue'
            elif interaction['type'] == 'hydrophobic':
                color = 'yellow'
            elif interaction['type'] == 'pi_stacking':
                color = 'green'
            else:
                color = 'gray'
            
            # Add cylinder to represent interaction
            start = interaction['atom1_coords']
            end = interaction['atom2_coords']
            self.viewer.page().runJavaScript(f'''
                viewer.addCylinder({{
                    start: {{x: {start[0]}, y: {start[1]}, z: {start[2]}}},
                    end: {{x: {end[0]}, y: {end[1]}, z: {end[2]}}},
                    radius: 0.1,
                    color: "{color}",
                    dashed: true
                }});
                viewer.render();
            ''')
    
    def export_image(self, filename: str):
        """Export the current view as an image."""
        self.viewer.page().runJavaScript('''
            var png = viewer.pngURI();
            // Return PNG data
        ''', lambda png: self._save_image(png, filename))
    
    def _save_image(self, png_data: str, filename: str):
        """Save PNG data to file."""
        import base64
        if png_data and png_data.startswith('data:image/png;base64,'):
            data = base64.b64decode(png_data.split(',')[1])
            with open(filename, 'wb') as f:
                f.write(data)
