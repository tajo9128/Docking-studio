"""
Base Docking Engine Interface
Abstract base class for all docking engines.
"""

from abc import ABC, abstractmethod
from typing import Dict, Any, List, Optional
from dataclasses import dataclass
from pathlib import Path


@dataclass
class DockingResult:
    """Container for docking results"""
    status: str
    binding_energy: Optional[float] = None
    num_modes: int = 0
    poses: List[Dict[str, Any]] = None
    gnina_cnn_score: Optional[float] = None
    gnina_cnn_affinity: Optional[float] = None
    logs: str = ""
    job_id: str = ""
    
    def __post_init__(self):
        if self.poses is None:
            self.poses = []


@dataclass
class DockingConfig:
    """Configuration for docking run"""
    receptor_file: str
    ligand_file: str
    center_x: float
    center_y: float
    center_z: float
    size_x: float = 22.0
    size_y: float = 22.0
    size_z: float = 22.0
    exhaustiveness: int = 8
    num_modes: int = 9
    energy_range: float = 3.0
    cpu: int = 0
    job_id: str = ""
    output_dir: str = "data"
    
    @property
    def grid_box(self) -> Dict[str, float]:
        return {
            "center_x": self.center_x,
            "center_y": self.center_y,
            "center_z": self.center_z,
            "size_x": self.size_x,
            "size_y": self.size_y,
            "size_z": self.size_z
        }


class DockingEngine(ABC):
    """Abstract base class for docking engines"""
    
    @property
    @abstractmethod
    def name(self) -> str:
        """Engine name"""
        pass
    
    @property
    @abstractmethod
    def requires_gpu(self) -> bool:
        """Whether this engine requires GPU"""
        pass
    
    @property
    @abstractmethod
    def supports_cnn_scoring(self) -> bool:
        """Whether this engine supports CNN scoring (GNINA)"""
        pass
    
    @abstractmethod
    def run(self, config: DockingConfig) -> DockingResult:
        """
        Run docking with given configuration.
        
        Args:
            config: Docking configuration
            
        Returns:
            DockingResult: Results from docking run
        """
        pass
    
    @abstractmethod
    def parse_output(self, output_file: str) -> Dict[str, Any]:
        """
        Parse docking output file.
        
        Args:
            output_file: Path to output file
            
        Returns:
            dict: Parsed results
        """
        pass
    
    def validate_config(self, config: DockingConfig) -> List[str]:
        """
        Validate configuration and return list of errors.
        
        Args:
            config: Configuration to validate
            
        Returns:
            list: List of error messages (empty if valid)
        """
        errors = []
        
        if not Path(config.receptor_file).exists():
            errors.append(f"Receptor file not found: {config.receptor_file}")
        
        if not Path(config.ligand_file).exists():
            errors.append(f"Ligand file not found: {config.ligand_file}")
        
        if config.size_x <= 0 or config.size_y <= 0 or config.size_z <= 0:
            errors.append("Grid size must be positive")
        
        return errors
