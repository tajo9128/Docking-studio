"""
BioDockify Interaction Diagram Widget
2D Ligand Interaction Diagram using PyQt6 graphics
"""

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, 
    QLabel, QScrollArea, QFrame, QGraphicsView, 
    QGraphicsScene, QGraphicsEllipseItem, QGraphicsLineItem,
    QGraphicsTextItem, QToolBar, QFileDialog
)
from PyQt6.QtCore import Qt, QPointF, QRectF
from PyQt6.QtGui import (
    QPen, QBrush, QColor, QFont, QPainter, QPixmap
)
import math
from typing import List, Dict, Tuple


class InteractionDiagram(QWidget):
    """
    2D Ligand Interaction Diagram Widget.
    Shows protein-ligand interactions in a schematic 2D view.
    """
    
    # Interaction type colors
    COLORS = {
        'hbond': QColor(0, 100, 255),       # Blue
        'hydrophobic': QColor(255, 200, 0),  # Yellow
        'pi_stacking': QColor(0, 200, 100),  # Green
        'salt_bridge': QColor(255, 0, 100),  # Magenta
        'metal': QColor(150, 150, 150),      # Gray
        'halogen': QColor(0, 200, 200),      # Cyan
    }
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.interactions = []
        self.residues = []
        self.ligand_atoms = []
        self.setup_ui()
        
    def setup_ui(self):
        """Initialize the UI components."""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Toolbar
        toolbar = self.create_toolbar()
        layout.addWidget(toolbar)
        
        # Graphics view
        self.scene = QGraphicsScene()
        self.view = QGraphicsView(self.scene)
        self.view.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.view.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)
        self.view.setMinimumSize(400, 400)
        layout.addWidget(self.view)
        
        # Legend
        legend = self.create_legend()
        layout.addWidget(legend)
        
    def create_toolbar(self):
        """Create the toolbar."""
        toolbar = QToolBar()
        toolbar.setMovable(False)
        
        self.btn_zoom_in = QPushButton("Zoom In")
        self.btn_zoom_in.clicked.connect(self.zoom_in)
        toolbar.addWidget(self.btn_zoom_in)
        
        self.btn_zoom_out = QPushButton("Zoom Out")
        self.btn_zoom_out.clicked.connect(self.zoom_out)
        toolbar.addWidget(self.btn_zoom_out)
        
        self.btn_fit = QPushButton("Fit View")
        self.btn_fit.clicked.connect(self.fit_view)
        toolbar.addWidget(self.btn_fit)
        
        toolbar.addSeparator()
        
        self.btn_export = QPushButton("Export PNG")
        self.btn_export.clicked.connect(self.export_image)
        toolbar.addWidget(self.btn_export)
        
        return toolbar
    
    def create_legend(self):
        """Create the legend panel."""
        frame = QFrame()
        frame.setFrameShape(QFrame.Shape.StyledPanel)
        layout = QHBoxLayout(frame)
        
        for interaction_type, color in self.COLORS.items():
            # Color box
            color_label = QLabel()
            color_label.setFixedSize(16, 16)
            color_label.setStyleSheet(f"background-color: {color.name()}; border: 1px solid black;")
            layout.addWidget(color_label)
            
            # Label
            text_label = QLabel(interaction_type.replace('_', ' ').title())
            layout.addWidget(text_label)
            
            layout.addSpacing(10)
        
        layout.addStretch()
        return frame
    
    def set_data(self, ligand_atoms: List[Dict], residues: List[Dict], interactions: List[Dict]):
        """
        Set the data for the diagram.
        
        Args:
            ligand_atoms: List of ligand atoms with name, coords
            residues: List of protein residues with name, chain, number
            interactions: List of interactions with type, atoms
        """
        self.ligand_atoms = ligand_atoms
        self.residues = residues
        self.interactions = interactions
        self.draw_diagram()
    
    def draw_diagram(self):
        """Draw the interaction diagram."""
        self.scene.clear()
        
        if not self.ligand_atoms:
            return
        
        # Center position
        center_x, center_y = 300, 300
        ligand_radius = 80
        residue_radius = 200
        
        # Draw ligand in center
        ligand_ellipse = self.scene.addEllipse(
            center_x - ligand_radius/2, 
            center_y - ligand_radius/2,
            ligand_radius, ligand_radius,
            QPen(QColor(100, 100, 100), 2),
            QBrush(QColor(200, 200, 255))
        )
        
        # Label for ligand
        ligand_label = self.scene.addText("LIGAND", QFont("Arial", 10, QFont.Weight.Bold))
        ligand_label.setPos(center_x - 25, center_y - 10)
        
        # Draw residues around the ligand
        num_residues = len(self.residues)
        if num_residues == 0:
            return
            
        angle_step = 2 * math.pi / num_residues
        residue_positions = {}
        
        for i, residue in enumerate(self.residues):
            angle = i * angle_step - math.pi / 2  # Start from top
            x = center_x + residue_radius * math.cos(angle)
            y = center_y + residue_radius * math.sin(angle)
            
            residue_positions[residue.get('id', i)] = (x, y)
            
            # Determine residue color based on properties
            res_name = residue.get('name', 'UNK')
            color = self.get_residue_color(res_name)
            
            # Draw residue circle
            self.scene.addEllipse(
                x - 25, y - 25, 50, 50,
                QPen(QColor(0, 0, 0), 1),
                QBrush(color)
            )
            
            # Residue label
            label = f"{res_name}{residue.get('number', '')}"
            text = self.scene.addText(label, QFont("Arial", 8))
            text.setPos(x - 20, y - 8)
        
        # Draw interactions
        for interaction in self.interactions:
            int_type = interaction.get('type', 'hbond')
            res_id = interaction.get('residue_id')
            
            if res_id in residue_positions:
                res_x, res_y = residue_positions[res_id]
                
                # Draw line from ligand to residue
                pen = QPen(self.COLORS.get(int_type, QColor(128, 128, 128)), 2)
                
                if int_type == 'hbond':
                    pen.setStyle(Qt.PenStyle.DashLine)
                elif int_type == 'hydrophobic':
                    pen.setStyle(Qt.PenStyle.DotLine)
                else:
                    pen.setStyle(Qt.PenStyle.SolidLine)
                
                self.scene.addLine(
                    center_x, center_y,
                    res_x, res_y,
                    pen
                )
                
                # Add distance label
                distance = interaction.get('distance', 0)
                if distance > 0:
                    mid_x = (center_x + res_x) / 2
                    mid_y = (center_y + res_y) / 2
                    dist_text = self.scene.addText(f"{distance:.1f}Å", QFont("Arial", 7))
                    dist_text.setPos(mid_x, mid_y)
    
    def get_residue_color(self, res_name: str) -> QColor:
        """Get color based on residue type."""
        # Hydrophobic
        if res_name in ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO']:
            return QColor(255, 255, 150)  # Yellow
        # Polar
        elif res_name in ['SER', 'THR', 'CYS', 'TYR', 'ASN', 'GLN']:
            return QColor(150, 255, 150)  # Green
        # Charged positive
        elif res_name in ['LYS', 'ARG', 'HIS']:
            return QColor(150, 150, 255)  # Blue
        # Charged negative
        elif res_name in ['ASP', 'GLU']:
            return QColor(255, 150, 150)  # Red
        # Glycine
        elif res_name == 'GLY':
            return QColor(255, 200, 200)  # Pink
        else:
            return QColor(200, 200, 200)  # Gray
    
    def zoom_in(self):
        """Zoom in the view."""
        self.view.scale(1.2, 1.2)
    
    def zoom_out(self):
        """Zoom out the view."""
        self.view.scale(0.8, 0.8)
    
    def fit_view(self):
        """Fit the entire diagram in view."""
        self.view.fitInView(self.scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio)
    
    def export_image(self):
        """Export the diagram as PNG."""
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export Diagram", "", "PNG Images (*.png)"
        )
        if filename:
            self.save_to_file(filename)
    
    def save_to_file(self, filename: str):
        """Save the diagram to a file."""
        rect = self.scene.sceneRect()
        pixmap = QPixmap(int(rect.width()), int(rect.height()))
        pixmap.fill(Qt.GlobalColor.white)
        
        painter = QPainter(pixmap)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.scene.render(painter)
        painter.end()
        
        pixmap.save(filename)
    
    def clear(self):
        """Clear the diagram."""
        self.scene.clear()
        self.interactions = []
        self.residues = []
        self.ligand_atoms = []
