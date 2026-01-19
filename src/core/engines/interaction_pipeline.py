"""Interaction Pipeline - Handles Scientific Analysis"""

from typing import List, Dict, Optional
import math

from ..spatial_hash import SpatialHashGrid
from ..logging_config import get_logger

logger = get_logger(__name__)

class AnalysisThresholds:
    """Analysis thresholds (literature-based)"""
    
    def __init__(self):
        self.HYDROGEN_BOND = {
            'min': 1.5, 'max': 2.5, 'angle_min': 120,
            'donors': ['N', 'O', 'S', 'F', 'Cl', 'Br', 'I'],
            'acceptors': ['N', 'O', 'S', 'F', 'Cl', 'Br', 'I'],
        }
        
        self.SALT_BRIDGE = {
            'distance_max': 4.0,
            'positive_residues': ['LYS', 'ARG', 'HIS', 'LYS+', 'ARG+', 'HIS+'],
            'negative_residues': ['ASP', 'GLU', 'ASP-', 'GLU-'],
        }
        
        self.VDW = {'min': 0.7, 'max': 1.1}
        
        self.VDW_RADII = {
            'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
            'F': 1.47, 'P': 1.80, 'S': 1.80, 'Cl': 1.75,
            'Br': 1.85, 'I': 1.98, 'Fe': 2.00, 'Mg': 1.73,
            'Ca': 2.31, 'Mn': 2.00, 'Zn': 1.39,
        }
    
    def dict(self) -> Dict:
        return {
            'hydrogen_bond': self.HYDROGEN_BOND,
            'salt_bridge': self.SALT_BRIDGE,
            'vdw': self.VDW,
        }

class InteractionPipeline:
    """Interaction pipeline for scientific analysis"""
    
    def __init__(self):
        self.thresholds = AnalysisThresholds()
    
    def analyze(self, atoms: List[dict], bonds: List[dict]) -> Dict[str, List[dict]]:
        """Analyze molecular interactions"""
        logger.info(f"Analyzing {len(atoms)} atoms")
        
        interactions = {'hydrogen_bonds': [], 'vdw_contacts': [], 'salt_bridges': []}
        seen_atom_pairs = set()
        grid = SpatialHashGrid(atoms)
        
        for i in range(len(atoms)):
            atom1 = atoms[i]
            neighbors = grid.get_neighbors(i, atoms)
            
            for j in neighbors:
                if j <= i:
                    continue
                
                atom2 = atoms[j]
                pair_key = tuple(sorted([i, j]))
                
                if pair_key in seen_atom_pairs:
                    continue
                
                seen_atom_pairs.add(pair_key)
                distance = self._calculate_distance(atom1, atom2)
                
                if self._is_hydrogen_bond(atom1, atom2, distance):
                    interactions['hydrogen_bonds'].append({
                        'atom1_index': i, 'atom2_index': j, 'distance': distance,
                        'angle': self._calculate_angle(atom1, atom2, atoms, i, j),
                        'atom1_residue': atom1.get('res_name', ''),
                        'atom1_residue_seq': atom1.get('res_seq', 0),
                        'atom2_residue': atom2.get('res_name', ''),
                        'atom2_residue_seq': atom2.get('res_seq', 0),
                        'confidence': 1.0,
                    })
                
                if self._is_salt_bridge(atom1, atom2, distance):
                    interactions['salt_bridges'].append({
                        'atom1_index': i, 'atom2_index': j, 'distance': distance,
                        'atom1_residue': atom1.get('res_name', ''),
                        'atom1_residue_seq': atom1.get('res_seq', 0),
                        'atom2_residue': atom2.get('res_name', ''),
                        'atom2_residue_seq': atom2.get('res_seq', 0),
                        'confidence': 1.0,
                    })
                
                if self._is_vdw_contact(atom1, atom2, distance):
                    interactions['vdw_contacts'].append({
                        'atom1_index': i, 'atom2_index': j, 'distance': distance,
                        'atom1_residue': atom1.get('res_name', ''),
                        'atom1_residue_seq': atom1.get('res_seq', 0),
                        'atom2_residue': atom2.get('res_name', ''),
                        'atom2_residue_seq': atom2.get('res_seq', 0),
                        'confidence': 1.0,
                    })
        
        logger.info(f"Found {len(interactions['hydrogen_bonds'])} H-bonds")
        return interactions
    
    def _calculate_distance(self, atom1: dict, atom2: dict) -> float:
        dx = atom1['x'] - atom2['x']
        dy = atom1['y'] - atom2['y']
        dz = atom1['z'] - atom2['z']
        return math.sqrt(dx**2 + dy**2 + dz**2)
    
    def _is_hydrogen_bond(self, atom1: dict, atom2: dict, distance: float) -> bool:
        if distance < self.thresholds.HYDROGEN_BOND['min'] or distance > self.thresholds.HYDROGEN_BOND['max']:
            return False
        
        if atom1.get('element', '') != 'H' and atom2.get('element', '') != 'H':
            return False
        
        return True
    
    def _is_salt_bridge(self, atom1: dict, atom2: dict, distance: float) -> bool:
        if distance > self.thresholds.SALT_BRIDGE['distance_max']:
            return False
        
        res1 = atom1.get('res_name', '')
        res2 = atom2.get('res_name', '')
        
        is_pos1 = any(res1.startswith(r) for r in self.thresholds.SALT_BRIDGE['positive_residues'])
        is_neg1 = any(res1.startswith(r) for r in self.thresholds.SALT_BRIDGE['negative_residues'])
        is_pos2 = any(res2.startswith(r) for r in self.thresholds.SALT_BRIDGE['positive_residues'])
        is_neg2 = any(res2.startswith(r) for r in self.thresholds.SALT_BRIDGE['negative_residues'])
        
        return (is_pos1 and is_neg2) or (is_pos2 and is_neg1)
    
    def _is_vdw_contact(self, atom1: dict, atom2: dict, distance: float) -> bool:
        r1 = self.thresholds.VDW_RADII.get(atom1.get('element', 'C'), 1.70)
        r2 = self.thresholds.VDW_RADII.get(atom2.get('element', 'C'), 1.70)
        vdw_sum = r1 + r2
        
        min_dist = self.thresholds.VDW['min'] * vdw_sum
        max_dist = self.thresholds.VDW['max'] * vdw_sum
        
        return min_dist <= distance <= max_dist
    
    def _calculate_angle(self, atom1: dict, atom2: dict, atoms: List[dict], i: int, j: int) -> Optional[float]:
        return None  # Simplified for now
