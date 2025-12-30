"""
BioDockify Docking Studio - Upload Widget
Handles file upload (receptor and ligand)
"""

from PyQt6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QFileDialog, QDragEnterEvent, QDropEvent
from PyQt6.QtCore import QMimeData, pyqtSignal, Qt
from PyQt6.QtGui import QPalette, QFrame
import logging
from pathlib import Path
import mimetypes

logger = logging.getLogger(__name__)

class UploadWidget(QWidget):
    """File upload widget for receptor and ligand"""
    
    # Signals
    receptor_uploaded = pyqtSignal(str, str)  # filepath, filename
    ligand_uploaded = pyqtSignal(str, str)  # filepath, filename
    upload_failed = pyqtSignal(str, str)  # upload_type, error_message
    
    def __init__(self):
        """Initialize upload widget"""
        super().__init__()
        self.receptor_path = None
        self.ligand_path = None
        self._setup_ui()
    
    def _setup_ui(self) -> None:
        """Setup UI components"""
        layout = QVBoxLayout(self)
        
        # Receptor section
        receptor_frame = QFrame()
        receptor_frame.setStyleSheet("border: 2px solid #E0E0E0; border-radius: 8px; padding: 10px;")
        
        receptor_label = QLabel("Receptor (Protein):")
        receptor_label.setStyleSheet("font-size: 14px; font-weight: bold; color: #2196F3;")
        
        self.receptor_file_label = QLabel("Drop file here or click to upload")
        self.receptor_file_label.setStyleSheet("padding: 10px; color: #666; font-style: italic;")
        
        self.receptor_browse_button = QPushButton("Browse")
        self.receptor_browse_button.setStyleSheet("background-color: #4CAF50; color: white; border: none; padding: 8px; border-radius: 5px;")
        self.receptor_browse_button.clicked.connect(self.browse_receptor)
        
        receptor_layout = QVBoxLayout(receptor_frame)
        receptor_layout.addWidget(receptor_label)
        receptor_layout.addWidget(self.receptor_file_label)
        receptor_layout.addWidget(self.receptor_browse_button)
        receptor_layout.addStretch()
        
        # Ligand section
        ligand_frame = QFrame()
        ligand_frame.setStyleSheet("border: 2px solid #E0E0E0; border-radius: 8px; padding: 10px;")
        
        ligand_label = QLabel("Ligand (Small Molecule):")
        ligand_label.setStyleSheet("font-size: 14px; font-weight: bold; color: #2196F3;")
        
        self.ligand_file_label = QLabel("Drop file here or click to upload")
        self.ligand_file_label.setStyleSheet("padding: 10px; color: #666; font-style: italic;")
        
        self.ligand_browse_button = QPushButton("Browse")
        self.ligand_browse_button.setStyleSheet("background-color: #4CAF50; color: white; border: none; padding: 8px; border-radius: 5px;")
        self.ligand_browse_button.clicked.connect(self.browse_ligand)
        
        ligand_layout = QVBoxLayout(ligand_frame)
        ligand_layout.addWidget(ligand_label)
        ligand_layout.addWidget(self.ligand_file_label)
        ligand_layout.addWidget(self.ligand_browse_button)
        ligand_layout.addStretch()
        
        # Add frames to main layout
        layout.addWidget(receptor_frame)
        layout.addWidget(ligand_frame)
        
        layout.addStretch()
        
        self.setLayout(layout)
        
        # Enable drag and drop
        self.setAcceptDrops(True)
    
    def browse_receptor(self) -> None:
        """Browse for receptor file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Receptor File",
            "",
            "PDB Files (*.pdb); PDBQT Files (*.pdbqt); All Files (*)"
        )
        
        if file_path:
            self._process_receptor_upload(file_path)
    
    def browse_ligand(self) -> None:
        """Browse for ligand file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Ligand File",
            "",
            "SDF Files (*.sdf); MOL2 Files (*.mol2); PDB Files (*.pdb); All Files (*)"
        )
        
        if file_path:
            self._process_ligand_upload(file_path)
    
    def _process_receptor_upload(self, file_path: str) -> None:
        """Process receptor file upload"""
        logger.info(f"Processing receptor upload: {file_path}")
        
        try:
            # Validate file exists
            if not Path(file_path).exists():
                error_message = "File not found"
                self.upload_failed.emit("receptor", error_message)
                self.receptor_file_label.setText(f"Error: {error_message}")
                self.receptor_file_label.setStyleSheet("padding: 10px; color: #FF4444; font-weight: bold;")
                return
            
            # Validate file size
            file_size = Path(file_path).stat().st_size
            max_size = 50 * 1024 * 1024  # 50 MB
            if file_size > max_size:
                error_message = f"File too large ({file_size / (1024*1024):.1f} MB). Maximum size is {max_size / (1024*1024):.1f} MB."
                self.upload_failed.emit("receptor", error_message)
                self.receptor_file_label.setText(f"Error: {error_message}")
                self.receptor_file_label.setStyleSheet("padding: 10px; color: #FF4444; font-weight: bold;")
                return
            
            # Get file type and validate
            file_type = self._get_file_type(file_path)
            if file_type not in ["pdb", "pdbqt"]:
                error_message = f"Invalid file type: {file_type}. Supported types: PDB, PDBQT."
                self.upload_failed.emit("receptor", error_message)
                self.receptor_file_label.setText(f"Error: {error_message}")
                self.receptor_file_label.setStyleSheet("padding: 10px; color: #FF4444; font-weight: bold;")
                return
            
            # Store file path
            self.receptor_path = file_path
            self.receptor_file_label.setText(Path(file_path).name)
            self.receptor_file_label.setStyleSheet("padding: 10px; color: #666;")
            
            # Emit signal
            self.receptor_uploaded.emit(file_path, Path(file_path).name)
            logger.info(f"Receptor uploaded successfully: {file_path}")
        
        except Exception as e:
            error_message = f"Failed to process receptor file: {str(e)}"
            self.upload_failed.emit("receptor", error_message)
            self.receptor_file_label.setText(f"Error: {error_message}")
            self.receptor_file_label.setStyleSheet("padding: 10px; color: #FF4444; font-weight: bold;")
            logger.error(error_message)
    
    def _process_ligand_upload(self, file_path: str) -> None:
        """Process ligand file upload"""
        logger.info(f"Processing ligand upload: {file_path}")
        
        try:
            # Validate file exists
            if not Path(file_path).exists():
                error_message = "File not found"
                self.upload_failed.emit("ligand", error_message)
                self.ligand_file_label.setText(f"Error: {error_message}")
                self.ligand_file_label.setStyleSheet("padding: 10px; color: #FF4444; font-weight: bold;")
                return
            
            # Validate file size
            file_size = Path(file_path).stat().st_size
            max_size = 50 * 1024 * 1024  # 50 MB
            if file_size > max_size:
                error_message = f"File too large ({file_size / (1024*1024):.1f} MB). Maximum size is {max_size / (1024*1024):.1f} MB."
                self.upload_failed.emit("ligand", error_message)
                self.ligand_file_label.setText(f"Error: {error_message}")
                self.ligand_file_label.setStyleSheet("padding: 10px; color: #FF4444; font-weight: bold;")
                return
            
            # Get file type and validate
            file_type = self._get_file_type(file_path)
            if file_type not in ["sdf", "mol2", "pdb"]:
                error_message = f"Invalid file type: {file_type}. Supported types: SDF, MOL2, PDB."
                self.upload_failed.emit("ligand", error_message)
                self.ligand_file_label.setText(f"Error: {error_message}")
                self.ligand_file_label.setStyleSheet("padding: 10px; color: #FF4444; font-weight: bold;")
                return
            
            # Store file path
            self.ligand_path = file_path
            self.ligand_file_label.setText(Path(file_path).name)
            self.ligand_file_label.setStyleSheet("padding: 10px; color: #666;")
            
            # Emit signal
            self.ligand_uploaded.emit(file_path, Path(file_path).name)
            logger.info(f"Ligand uploaded successfully: {file_path}")
        
        except Exception as e:
            error_message = f"Failed to process ligand file: {str(e)}"
            self.upload_failed.emit("ligand", error_message)
            self.ligand_file_label.setText(f"Error: {error_message}")
            self.ligand_file_label.setStyleSheet("padding: 10px; color: #FF4444; font-weight: bold;")
            logger.error(error_message)
    
    def _get_file_type(self, file_path: str) -> str:
        """Get file type from extension"""
        path_obj = Path(file_path)
        return path_obj.suffix[1:].lower()
    
    def dragEnterEvent(self, event: QDragEnterEvent) -> None:
        """Handle drag enter event"""
        event.acceptProposedAction()
        logger.debug("Drag enter event")
    
    def dropEvent(self, event: QDropEvent) -> None:
        """Handle drop event"""
        logger.debug("Drop event detected")
        
        # Process dropped files
        for url in event.mimeData().urls():
            file_path = url.toLocalFile()
            file_name = url.fileName()
            file_type = self._get_file_type(file_path)
            
            # Validate file
            if Path(file_path).exists():
                if file_type in ["pdb", "pdbqt"]:
                    self._process_receptor_upload(file_path)
                elif file_type in ["sdf", "mol2", "pdb"]:
                    self._process_ligand_upload(file_path)
            else:
                error_message = f"Unsupported file type: {file_type}. Expected PDB/PDBQT for receptor or SDF/MOL2/PDB for ligand."
                self.upload_failed.emit("unknown", error_message)
    
    def get_receptor_path(self) -> Optional[str]:
        """Get receptor file path"""
        return self.receptor_path
    
    def get_ligand_path(self) -> Optional[str]:
        """Get ligand file path"""
        return self.ligand_path
    
    def reset(self) -> None:
        """Reset upload widget"""
        self.receptor_path = None
        self.ligand_path = None
        self.receptor_file_label.setText("Drop file here or click to upload")
        self.receptor_file_label.setStyleSheet("padding: 10px; color: #666; font-style: italic;")
        self.ligand_file_label.setText("Drop file here or click to upload")
        self.ligand_file_label.setStyleSheet("padding: 10px; color: #666; font-style: italic;")
        
        logger.info("Upload widget reset")
