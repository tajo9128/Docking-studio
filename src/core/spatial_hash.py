"""Spatial Hash Grid for O(n) Complexity"""

import numpy as np
from typing import List, Set
from ..config import settings
from ..logging_config import get_logger

logger = get_logger(__name__)

class SpatialHashGrid:
    """Spatial hash grid for O(n) neighbor search complexity"""
    
    def __init__(self, atoms: List[dict]):
        """Initialize spatial hash grid"""
        self.cell_size = settings.SPATIAL_GRID_CELL_SIZE
        self.grid = {}
        self.seen_atom_pairs = set()
        
        # Calculate optimal cell size
        if atoms:
            coords = np.array([[atom['x'], atom['y'], atom['z']] for atom in atoms])
            dimensions = coords.max(axis=0) - coords.min(axis=0)
            max_dim = np.max(dimensions)
            
            cells_per_dim = max(1, int(np.ceil(len(atoms) / 1000)))
            cell_size = max_dim / cells_per_dim
            self.cell_size = max(cell_size, 5.0)
        
        self._build_grid(atoms)
    
    def _build_grid(self, atoms: List[dict]) -> None:
        """Build spatial hash grid"""
        for i, atom in enumerate(atoms):
            x, y, z = atom['x'], atom['y'], atom['z']
            cell_key = self._get_cell_key(x, y, z)
            
            if cell_key not in self.grid:
                self.grid[cell_key] = []
            self.grid[cell_key].append(i)
    
    def _get_cell_key(self, x: float, y: float, z: float) -> str:
        """Get cell key for coordinates"""
        cell_x = int(np.floor(x / self.cell_size))
        cell_y = int(np.floor(y / self.cell_size))
        cell_z = int(np.floor(z / self.cell_size))
        return f"{cell_x},{cell_y},{cell_z}"
    
    def get_neighbors(self, atom_index: int, atoms: List[dict]) -> List[int]:
        """Get neighbor atom indices (in same or adjacent cells)"""
        atom = atoms[atom_index]
        x, y, z = atom['x'], atom['y'], atom['z']
        neighbors = []
        
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                for dz in [-1, 0, 1]:
                    cell_key = self._get_cell_key(
                        x + dx * self.cell_size,
                        y + dy * self.cell_size,
                        z + dz * self.cell_size
                    )
                    cell_atoms = self.grid.get(cell_key, [])
                    neighbors.extend([a for a in cell_atoms if a != atom_index])
        
        return list(set(neighbors))
