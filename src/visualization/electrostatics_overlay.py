"""
Electrostatics Overlay
Electrostatic potential visualization on molecular surfaces.
"""

from enum import Enum
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
import numpy as np
import logging

logger = logging.getLogger(__name__)


class ElectrostaticMethod(Enum):
    """Methods for electrostatic calculation"""
    APBS = "apbs"              # Adaptive Poisson-Boltzmann Solver
    DELPHI = "delphi"          # DelPhi
    GAUSSIAN = "gaussian"      # Gaussian
    APPROXIMATE = "approximate" # Quick approximation


@dataclass
class ElectrostaticMap:
    """Electrostatic potential map"""
    potentials: np.ndarray      # Potential values at each vertex
    vertices: np.ndarray         # Vertex positions
    colors: np.ndarray          # RGB colors (red-white-blue)
    min_potential: float        # Minimum potential (kT/e)
    max_potential: float        # Maximum potential (kT/e)


class ElectrostaticsOverlay:
    """
    Generates and displays electrostatic potential maps on molecular surfaces.
    
    Uses simplified charge model for quick visualization.
    For production use, integrates with APBS or DelPhi via Docker.
    """
    
    # Color scheme for potentials
    NEGATIVE_COLOR = [1.0, 0.2, 0.2]   # Red (negative)
    NEUTRAL_COLOR = [1.0, 1.0, 1.0]    # White (neutral)
    POSITIVE_COLOR = [0.2, 0.2, 1.0]   # Blue (positive)
    
    # Charge values for common atoms (partial charges)
    ATOM_CHARGES = {
        'N': 0.5, 'O': -0.5, 'S': 0.0, 'C': 0.0,
        'H': 0.2, 'P': 0.5, 'FE': 1.0, 'ZN': 0.5,
        'MG': 0.6, 'CA': 0.7, 'MN': 0.8, 'CO': 0.7,
    }
    
    def __init__(self, method: ElectrostaticMethod = ElectrostaticMethod.APPROXIMATE):
        """
        Initialize electrostatics calculator.
        
        Args:
            method: Calculation method to use
        """
        self.method = method
        self.current_map: Optional[ElectrostaticMap] = None
        self.scale_factor = 1.0
        
        logger.info(f"ElectrostaticsOverlay initialized (method: {method.value})")
    
    def calculate_potential(
        self,
        vertices: np.ndarray,
        atoms: List[Dict],
        charge_model: str = " Kollman"
    ) -> ElectrostaticMap:
        """
        Calculate electrostatic potential at each surface vertex.
        
        Args:
            vertices: Surface vertex positions (N x 3)
            atoms: List of atom dictionaries with x, y, z, elem, charge
            charge_model: Charge model to use
            
        Returns:
            ElectrostaticMap with calculated potentials
        """
        n_verts = len(vertices)
        potentials = np.zeros(n_verts)
        
        # Extract atom data
        atom_coords = np.array([(a['x'], a['y'], a['z']) for a in atoms])
        atom_charges = np.array([
            a.get('charge', self.ATOM_CHARGES.get(a.get('elem', 'C').upper(), 0.0))
            for a in atoms
        ])
        
        # Calculate potential at each vertex using Coulomb's law
        # V = sum(qi / ri) for all atoms
        logger.info(f"Calculating electrostatic potential for {n_verts} vertices, {len(atoms)} atoms")
        
        for i, vert in enumerate(vertices):
            # Distance to each atom
            distances = np.linalg.norm(atom_coords - vert, axis=1)
            
            # Avoid division by zero
            distances[distances < 0.1] = 0.1
            
            # Calculate potential (simplified Coulomb)
            potential = np.sum(atom_charges / distances)
            potentials[i] = potential
        
        # Normalize potentials
        max_pot = np.max(np.abs(potentials))
        if max_pot > 0:
            potentials = potentials / max_pot
        
        # Generate colors
        colors = self._potentials_to_colors(potentials)
        
        electrostatic_map = ElectrostaticMap(
            potentials=potentials,
            vertices=vertices,
            colors=colors,
            min_potential=np.min(potentials),
            max_potential=np.max(potentials)
        )
        
        self.current_map = electrostatic_map
        
        logger.info(f"Electrostatic potential calculated: range [{np.min(potentials):.2f}, {np.max(potentials):.2f}]")
        
        return electrostatic_map
    
    def _potentials_to_colors(self, potentials: np.ndarray) -> np.ndarray:
        """Convert potential values to RGB colors"""
        n = len(potentials)
        colors = np.zeros((n, 3))
        
        for i, pot in enumerate(potentials):
            if pot < 0:
                # Negative: interpolate white to red
                t = min(1.0, -pot)
                colors[i] = [
                    NEUTRAL_COLOR[0] * (1-t) + NEGATIVE_COLOR[0] * t,
                    NEUTRAL_COLOR[1] * (1-t) + NEGATIVE_COLOR[1] * t,
                    NEUTRAL_COLOR[2] * (1-t) + NEGATIVE_COLOR[2] * t,
                ]
            else:
                # Positive: interpolate white to blue
                t = min(1.0, pot)
                colors[i] = [
                    NEUTRAL_COLOR[0] * (1-t) + POSITIVE_COLOR[0] * t,
                    NEUTRAL_COLOR[1] * (1-t) + POSITIVE_COLOR[1] * t,
                    NEUTRAL_COLOR[2] * (1-t) + POSITIVE_COLOR[2] * t,
                ]
        
        return colors
    
    def get_3dmol_colors(self) -> List[str]:
        """Get colors as hex strings for 3Dmol.js"""
        if self.current_map is None:
            return []
        
        hex_colors = []
        for color in self.current_map.colors:
            r, g, b = int(color[0] * 255), int(color[1] * 255), int(color[2] * 255)
            hex_colors.append(f'#{r:02x}{g:02x}{b:02x}')
        
        return hex_colors
    
    def get_color_scale(self) -> Dict:
        """Get color scale information"""
        return {
            'negative': '#FF3333',
            'neutral': '#FFFFFF',
            'positive': '#3333FF',
            'min': float(np.min(self.current_map.potentials)) if self.current_map else 0,
            'max': float(np.max(self.current_map.potentials)) if self.current_map else 0,
        }
    
    def apply_to_surface(
        self,
        surface_vertices: np.ndarray,
        atoms: List[Dict]
    ) -> ElectrostaticMap:
        """Apply electrostatic calculation to surface"""
        return self.calculate_potential(surface_vertices, atoms)
    
    @staticmethod
    def create_gradient_legend() -> Dict:
        """Create gradient legend data for UI"""
        return {
            'type': 'gradient',
            'colors': [
                {'position': 0.0, 'color': '#FF3333', 'label': 'Negative (-)'},
                {'position': 0.5, 'color': '#FFFFFF', 'label': 'Neutral'},
                {'position': 1.0, 'color': '#3333FF', 'label': 'Positive (+)'},
            ],
            'unit': 'kT/e',
        }
    
    def set_scale_factor(self, factor: float):
        """Set scale factor for potential visualization"""
        self.scale_factor = max(0.1, min(10.0, factor))


class APBSElectrostatics:
    """APBS-based electrostatics calculation (requires Docker)"""
    
    def __init__(self, docker_manager=None):
        self.docker_manager = docker_manager
        self.apbs_available = False
    
    async def calculate_with_apbs(
        self,
        pdb_file: str,
        output_dir: str
    ) -> Dict:
        """Calculate electrostatics using APBS in Docker"""
        logger.info("APBS electrostatics calculation requested")
        
        # This would require proper APBS setup in Docker
        # Placeholder for now
        return {
            'status': 'not_implemented',
            'message': 'APBS integration requires Docker container setup'
        }
    
    def parse_apbs_output(self, output_file: str) -> np.ndarray:
        """Parse APBS output file"""
        # Parse potential values from APBS DX file
        logger.info(f"Parsing APBS output from {output_file}")
        return np.array([])


class ElectrostaticColoring:
    """Advanced electrostatic coloring modes"""
    
    @staticmethod
    def by_charge(atoms: List[Dict]) -> Dict[int, Tuple[float, float, float]]:
        """Color by partial charge"""
        colors = {}
        
        for i, atom in enumerate(atoms):
            charge = atom.get('charge', 0.0)
            
            if charge > 0:
                t = min(1.0, charge)
                colors[i] = (
                    1.0 - t,           # Red decreases
                    1.0 - t,           # Green decreases  
                    1.0 - 0.5 * t      # Blue moderate
                )
            elif charge < 0:
                t = min(1.0, -charge)
                colors[i] = (
                    1.0 - 0.5 * t,    # Red moderate
                    1.0 - t,          # Green decreases
                    1.0 - t           # Blue decreases
                )
            else:
                colors[i] = (0.8, 0.8, 0.8)  # Gray for neutral
        
        return colors
    
    @staticmethod
    def by_pI(protein_sequence: str) -> Dict[int, Tuple[float, float, float]]:
        """Color by predicted pI (simplified)"""
        colors = {}
        
        # Simple pKa values
        pKa = {
            'D': 3.9, 'E': 4.3, 'C': 8.3,
            'Y': 10.1, 'H': 6.0, 'K': 10.5, 'R': 12.5
        }
        
        for i, aa in enumerate(protein_sequence):
            if aa in pKa:
                pka = pKa[aa]
                # Color based on charge at pH 7
                if pka < 7:  # Acidic - negative
                    t = (7 - pka) / 7
                    colors[i] = (1.0, 1.0 - t, 1.0 - t)
                else:  # Basic - positive
                    t = (pka - 7) / 6
                    colors[i] = (1.0 - t, 1.0 - t, 1.0)
            else:
                colors[i] = (0.8, 0.8, 0.8)
        
        return colors
