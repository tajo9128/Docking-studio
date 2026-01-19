"""Molecular Engine - Handles Atoms and Bonds"""

from typing import List, Dict, Optional
import math

from ..spatial_hash import SpatialHashGrid
from ..logging_config import get_logger

logger = get_logger(__name__)

class MolecularEngine:
    """Molecular engine for handling atoms and bonds"""
    
    def __init__(self):
        self.atoms = []
        self.bonds = []
        self.spatial_grid: Optional[SpatialHashGrid] = None
        
        self.COVALENT_RADII = {
            'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66,
            'F': 0.57, 'P': 1.07, 'S': 1.05, 'Cl': 1.02,
            'Br': 1.20, 'I': 1.39, 'Fe': 1.32, 'Mg': 1.30,
            'Ca': 1.67, 'Mn': 1.39, 'Zn': 1.31,
        }
    
    def initialize(self, atoms: List[dict], bonds: List[dict] = None) -> None:
        """Initialize molecular engine with atoms and bonds"""
        logger.info(f"Initializing molecular engine with {len(atoms)} atoms")
        
        self.atoms = atoms
        self.bonds = bonds if bonds else []
        self.spatial_grid = SpatialHashGrid(atoms)
        
        if not self.bonds:
            logger.info("Detecting bonds...")
            self.bonds = self._detect_bonds_optimized()
    
    def _detect_bonds_optimized(self) -> List[dict]:
        """Detect bonds using spatial hashing (O(n) complexity)"""
        bonds = []
        seen_atom_pairs = set()
        
        for i in range(len(self.atoms)):
            atom1 = self.atoms[i]
            neighbors = self.spatial_grid.get_neighbors(i, self.atoms)
            
            for j in neighbors:
                if j <= i:
                    continue
                
                atom2 = self.atoms[j]
                pair_key = tuple(sorted([i, j]))
                
                if pair_key in seen_atom_pairs:
                    continue
                
                seen_atom_pairs.add(pair_key)
                distance = self._calculate_distance(atom1, atom2)
                
                r1 = self.COVALENT_RADII.get(atom1.get('element', 'C'), 0.76)
                r2 = self.COVALENT_RADII.get(atom2.get('element', 'C'), 0.76)
                covalent_distance = r1 + r2
                
                if distance > 0.5 and distance <= covalent_distance + 0.2:
                    ratio = distance / covalent_distance
                    bond_type = "single"
                    bond_order = 1
                    
                    if ratio <= 0.9:
                        bond_type = "triple"
                        bond_order = 3
                    elif ratio <= 0.95:
                        bond_type = "double"
                        bond_order = 2
                    elif ratio >= 1.0 and ratio <= 1.1:
                        bond_type = "aromatic"
                        bond_order = 1.5
                    
                    bonds.append({
                        'atom1_index': i,
                        'atom2_index': j,
                        'type': bond_type,
                        'order': bond_order,
                        'distance': distance,
                    })
        
        logger.info(f"Detected {len(bonds)} bonds")
        return bonds
    
    def _calculate_distance(self, atom1: dict, atom2: dict) -> float:
        """Calculate Euclidean distance between two atoms"""
        dx = atom1['x'] - atom2['x']
        dy = atom1['y'] - atom2['y']
        dz = atom1['z'] - atom2['z']
        return math.sqrt(dx**2 + dy**2 + dz**2)
    
    def get_bounding_box(self) -> Dict[str, float]:
        """Calculate bounding box of molecule"""
        if not self.atoms:
            return {'min_x': 0.0, 'max_x': 0.0, 'min_y': 0.0, 'max_y': 0.0, 'min_z': 0.0, 'max_z': 0.0}
        
        xs = [atom['x'] for atom in self.atoms]
        ys = [atom['y'] for atom in self.atoms]
        zs = [atom['z'] for atom in self.atoms]
        
        return {
            'min_x': min(xs), 'max_x': max(xs),
            'min_y': min(ys), 'max_y': max(ys),
            'min_z': min(zs), 'max_z': max(zs),
        }
    
    def get_center_of_mass(self) -> Dict[str, float]:
        """Calculate center of mass (geometric center)"""
        if not self.atoms:
            return {'x': 0.0, 'y': 0.0, 'z': 0.0}
        
        xs = [atom['x'] for atom in self.atoms]
        ys = [atom['y'] for atom in self.atoms]
        zs = [atom['z'] for atom in self.atoms]
        
        return {'x': sum(xs) / len(xs), 'y': sum(ys) / len(ys), 'z': sum(zs) / len(zs)}
