"""
BioDockify Docking Studio - Configuration Widget
Handles docking parameters configuration
"""

from PyQt6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QLabel, QSpinBox, QDoubleSpinBox, QComboBox, QCheckBox, QGroupBox
from PyQt6.QtCore import pyqtSignal
import logging

logger = logging.getLogger(__name__)

class ConfigurationWidget(QWidget):
    """Configuration widget for docking parameters"""
    
    # Signals
    parameters_changed = pyqtSignal(dict)  # Dictionary of parameters
    
    def __init__(self):
        """Initialize configuration widget"""
        super().__init__()
        self.parameters = {
            "center_x": 0.0,
            "center_y": 0.0,
            "center_z": 0.0,
            "size_x": 20.0,
            "size_y": 20.0,
            "size_z": 20.0,
            "exhaustiveness": 8,
            "num_modes": 9
        }
        self._setup_ui()
    
    def _setup_ui(self) -> None:
        """Setup UI components"""
        layout = QVBoxLayout(self)
        
        # Docking parameters group
        docking_params_group = QGroupBox("Docking Parameters")
        docking_params_group.setStyleSheet("border: 2px solid #E0E0E0; border-radius: 8px; padding: 10px;")
        
        # Search box center coordinates
        center_layout = QHBoxLayout()
        
        center_label = QLabel("Search Box Center (Å):")
        center_label.setStyleSheet("font-size: 12px; font-weight: bold; color: #2196F3;")
        center_layout.addWidget(center_label)
        
        self.center_x_input = QDoubleSpinBox()
        self.center_x_input.setRange(-100.0, 100.0)
        self.center_x_input.setDecimals(1)
        self.center_x_input.setValue(self.parameters["center_x"])
        self.center_x_input.setSingleStep(1.0)
        self.center_x_input.valueChanged.connect(self._on_parameter_changed)
        center_layout.addWidget(self.center_x_input)
        
        self.center_y_input = QDoubleSpinBox()
        self.center_y_input.setRange(-100.0, 100.0)
        self.center_y_input.setDecimals(1)
        self.center_y_input.setValue(self.parameters["center_y"])
        self.center_y_input.setSingleStep(1.0)
        self.center_y_input.valueChanged.connect(self._on_parameter_changed)
        center_layout.addWidget(self.center_y_input)
        
        self.center_z_input = QDoubleSpinBox()
        self.center_z_input.setRange(-100.0, 100.0)
        self.center_z_input.setDecimals(1)
        self.center_z_input.setValue(self.parameters["center_z"])
        self.center_z_input.setSingleStep(1.0)
        self.center_z_input.valueChanged.connect(self._on_parameter_changed)
        center_layout.addWidget(self.center_z_input)
        
        center_layout.addStretch()
        
        # Search box size
        size_layout = QHBoxLayout()
        
        size_label = QLabel("Search Box Size (Å):")
        size_label.setStyleSheet("font-size: 12px; font-weight: bold; color: #2196F3;")
        size_layout.addWidget(size_label)
        
        self.size_x_input = QDoubleSpinBox()
        self.size_x_input.setRange(10.0, 100.0)
        self.size_x_input.setDecimals(0)
        self.size_x_input.setValue(self.parameters["size_x"])
        self.size_x_input.setSingleStep(1.0)
        self.size_x_input.valueChanged.connect(self._on_parameter_changed)
        size_layout.addWidget(self.size_x_input)
        
        self.size_y_input = QDoubleSpinBox()
        self.size_y_input.setRange(10.0, 100.0)
        self.size_y_input.setDecimals(0)
        self.size_y_input.setValue(self.parameters["size_y"])
        self.size_y_input.setSingleStep(1.0)
        self.size_y_input.valueChanged.connect(self._on_parameter_changed)
        size_layout.addWidget(self.size_y_input)
        
        self.size_z_input = QDoubleSpinBox()
        self.size_z_input.setRange(10.0, 100.0)
        self.size_z_input.setDecimals(0)
        self.size_z_input.setValue(self.parameters["size_z"])
        self.size_z_input.setSingleStep(1.0)
        self.size_z_input.valueChanged.connect(self._on_parameter_changed)
        size_layout.addWidget(self.size_z_input)
        
        size_layout.addStretch()
        
        # Add to docking parameters group
        params_layout = QVBoxLayout(docking_params_group)
        params_layout.addLayout(center_layout)
        params_layout.addLayout(size_layout)
        
        # Sampling parameters group
        sampling_params_group = QGroupBox("Sampling Parameters")
        sampling_params_group.setStyleSheet("border: 2px solid #E0E0E0; border-radius: 8px; padding: 10px;")
        
        # Exhaustiveness
        exhaustiveness_layout = QHBoxLayout()
        
        exhaustiveness_label = QLabel("Exhaustiveness (1-32):")
        exhaustiveness_label.setStyleSheet("font-size: 12px; font-weight: bold; color: #2196F3;")
        exhaustiveness_layout.addWidget(exhaustiveness_label)
        
        self.exhaustiveness_input = QSpinBox()
        self.exhaustiveness_input.setRange(1, 32)
        self.exhaustiveness_input.setValue(self.parameters["exhaustiveness"])
        self.exhaustiveness_input.setSingleStep(1)
        self.exhaustiveness_input.valueChanged.connect(self._on_parameter_changed)
        exhaustiveness_layout.addWidget(self.exhaustiveness_input)
        
        exhaustiveness_layout.addStretch()
        
        # Number of modes
        num_modes_layout = QHBoxLayout()
        
        num_modes_label = QLabel("Number of Modes (1-20):")
        num_modes_label.setStyleSheet("font-size: 12px; font-weight: bold; color: #2196F3;")
        num_modes_layout.addWidget(num_modes_label)
        
        self.num_modes_input = QSpinBox()
        self.num_modes_input.setRange(1, 20)
        self.num_modes_input.setValue(self.parameters["num_modes"])
        self.num_modes_input.setSingleStep(1)
        self.num_modes_input.valueChanged.connect(self._on_parameter_changed)
        num_modes_layout.addWidget(self.num_modes_input)
        
        num_modes_layout.addStretch()
        
        # Add to sampling parameters group
        sampling_layout = QVBoxLayout(sampling_params_group)
        sampling_layout.addLayout(exhaustiveness_layout)
        sampling_layout.addLayout(num_modes_layout)
        
        # Advanced parameters
        advanced_group = QGroupBox("Advanced Parameters")
        advanced_group.setStyleSheet("border: 2px solid #E0E0E0; border-radius: 8px; padding: 10px;")
        
        # Flexible receptor
        flexible_receptor = QCheckBox("Enable Flexible Receptor Side Chains")
        flexible_receptor.setStyleSheet("color: #2196F3; padding: 5px;")
        
        # Auto parameters
        auto_params_layout = QVBoxLayout(advanced_group)
        auto_params_layout.addWidget(flexible_receptor)
        
        # Add all groups to main layout
        layout.addWidget(docking_params_group)
        layout.addWidget(sampling_params_group)
        layout.addWidget(advanced_group)
        
        layout.addStretch()
        
        self.setLayout(layout)
    
    def _on_parameter_changed(self, value) -> None:
        """Handle parameter change"""
        sender = self.sender()
        param_name = sender.objectName()
        
        if param_name == "center_x_input":
            self.parameters["center_x"] = value
        elif param_name == "center_y_input":
            self.parameters["center_y"] = value
        elif param_name == "center_z_input":
            self.parameters["center_z"] = value
        elif param_name == "size_x_input":
            self.parameters["size_x"] = value
        elif param_name == "size_y_input":
            self.parameters["size_y"] = value
        elif param_name == "size_z_input":
            self.parameters["size_z"] = value
        elif param_name == "exhaustiveness_input":
            self.parameters["exhaustiveness"] = value
        elif param_name == "num_modes_input":
            self.parameters["num_modes"] = value
        
        logger.debug(f"Parameter changed: {param_name} = {value}")
        
        # Emit signal
        self.parameters_changed.emit(self.parameters)
    
    def get_parameters(self) -> dict:
        """Get current parameters"""
        return self.parameters.copy()
    
    def set_parameters(self, params: dict) -> None:
        """Set parameters"""
        self.parameters.update(params)
        
        self.center_x_input.setValue(params.get("center_x", self.parameters["center_x"]))
        self.center_y_input.setValue(params.get("center_y", self.parameters["center_y"]))
        self.center_z_input.setValue(params.get("center_z", self.parameters["center_z"]))
        self.size_x_input.setValue(params.get("size_x", self.parameters["size_x"]))
        self.size_y_input.setValue(params.get("size_y", self.parameters["size_y"]))
        self.size_z_input.setValue(params.get("size_z", self.parameters["size_z"]))
        self.exhaustiveness_input.setValue(params.get("exhaustiveness", self.parameters["exhaustiveness"]))
        self.num_modes_input.setValue(params.get("num_modes", self.parameters["num_modes"]))
        
        logger.info(f"Parameters set: {params}")
    
    def reset(self) -> None:
        """Reset to default parameters"""
        self.parameters = {
            "center_x": 0.0,
            "center_y": 0.0,
            "center_z": 0.0,
            "size_x": 20.0,
            "size_y": 20.0,
            "size_z": 20.0,
            "exhaustiveness": 8,
            "num_modes": 9
        }
        self._update_ui()
        
        logger.info("Configuration widget reset")
    
    def _update_ui(self) -> None:
        """Update UI to match parameters"""
        self.center_x_input.blockSignals(True)
        self.center_y_input.blockSignals(True)
        self.center_z_input.blockSignals(True)
        self.size_x_input.blockSignals(True)
        self.size_y_input.blockSignals(True)
        self.size_z_input.blockSignals(True)
        self.exhaustiveness_input.blockSignals(True)
        self.num_modes_input.blockSignals(True)
        
        self.center_x_input.setValue(self.parameters["center_x"])
        self.center_y_input.setValue(self.parameters["center_y"])
        self.center_z_input.setValue(self.parameters["center_z"])
        self.size_x_input.setValue(self.parameters["size_x"])
        self.size_y_input.setValue(self.parameters["size_y"])
        self.size_z_input.setValue(self.parameters["size_z"])
        self.exhaustiveness_input.setValue(self.parameters["exhaustiveness"])
        self.num_modes_input.setValue(self.parameters["num_modes"])
        
        self.center_x_input.blockSignals(False)
        self.center_y_input.blockSignals(False)
        self.center_z_input.blockSignals(False)
        self.size_x_input.blockSignals(False)
        self.size_y_input.blockSignals(False)
        self.size_z_input.blockSignals(False)
        self.exhaustiveness_input.blockSignals(False)
        self.num_modes_input.blockSignals(False)
        
        logger.debug("UI updated with current parameters")
