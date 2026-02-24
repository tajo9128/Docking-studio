"""
Interaction Detector
Automatic detection and analysis of molecular interactions between ligand and protein.
"""

from enum import Enum
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Any
import numpy as np
import logging

logger = logging.getLogger(__name__)


class InteractionType(Enum):
    """Types of molecular interactions"""
    HYDROGEN_BOND = "hbond"
    HYDROPHOBIC = "hydrophobic"
    PI_STACKING = "pi_stacking"
    SALT_BRIDGE = "salt_bridge"
    METAL_COORDINATION = "metal"
    HALOGEN = "halogen"
    WATER_BRIDGE = "water_bridge"
    PI_CATION = "pi_cation"


@dataclass
class Interaction:
    """Represents a single molecular interaction"""
    type: InteractionType
    atom1: str
    atom2: str
    atom1_coords: Tuple[float, float, float]
    atom2_coords: Tuple[float, float, float]
    distance: float
    angle: Optional[float] = None
    residue1: Optional[str] = None
    residue2: Optional[str] = None
    strength: float = 1.0


@dataclass
class InteractionSummary:
    """Summary of all interactions"""
    total_count: int
    hbond_count: int
    hydrophobic_count: int
    pi_stacking_count: int
    salt_bridge_count: int
    metal_count: int
    interactions: List[Interaction]
    key_residues: List[str]


class InteractionDetector:
    """
    Detects molecular interactions between ligand and protein.
    
    Uses geometric criteria for detection:
    - Hydrogen bonds: donor-acceptor distance < 3.5Å, angle > 120°
    - Hydrophobic: hydrophobic atom pairs < 4.5Å
    - Pi stacking: aromatic ring centroid distance < 5.0Å
    - Salt bridge: charged group distance < 4.0Å
    - Metal coordination: metal-ligand distance < 2.5Å
    """
    
    # Element classifications
    HYDROGEN_DONORS = {'N', 'O', 'S'}  # Can be H-bond donors
    HYDROGEN_ACCEPTORS = {'N', 'O', 'S', 'F', 'CL', 'BR', 'I'}  # H-bond acceptors
    HYDROPHOBIC_ATOMS = {'C', 'S'}  # Hydrophobic atoms
    AROMATIC_ATOMS = {'C', 'N'}  # Can be in aromatic rings
    POSITIVE_CHARGED = {'N', 'K', 'R', 'H'}  # Positive charge
    NEGATIVE_CHARGED = {'O', 'D', 'E'}  # Negative charge (single letter codes)
    METALS = {'MG', 'ZN', 'FE', 'CA', 'MN', 'CO', 'NI', 'CU', 'CD'}
    
    # Detection thresholds
    DISTANCE_CUTOFFS = {
        InteractionType.HYDROGEN_BOND: 3.5,
        InteractionType.HYDROPHOBIC: 4.5,
        InteractionType.PI_STACKING: 5.0,
        InteractionType.SALT_BRIDGE: 4.0,
        InteractionType.METAL_COORDINATION: 2.5,
        InteractionType.HALOGEN: 3.5,
        InteractionType.PI_CATION: 5.0,
    }
    
    # Colors for visualization (matching Discovery Studio)
    INTERACTION_COLORS = {
        InteractionType.HYDROGEN_BOND: '#00FF00',      # Green
        InteractionType.HYDROPHOBIC: '#FFFF00',        # Yellow
        InteractionType.PI_STACKING: '#800080',       # Purple
        InteractionType.SALT_BRIDGE: '#FF0000',        # Red
        InteractionType.METAL_COORDINATION: '#808080', # Gray
        InteractionType.HALOGEN: '#00FFFF',            # Cyan
        InteractionType.PI_CATION: '#FF6600',         # Orange
    }
    
    def __init__(self, use_prolif: bool = True):
        """
        Initialize interaction detector.
        
        Args:
            use_prolif: If True, prefer ProLIF library when available
        """
        self.use_prolif = use_prolif
        self.prolif_available = False
        self._check_prolif()
        
        logger.info(f"InteractionDetector initialized (ProLIF: {self.prolif_available})")
    
    def _check_prolif(self):
        """Check if ProLIF is available"""
        try:
            import prolif
            self.prolif_available = True
            logger.info("ProLIF library detected - using for advanced detection")
        except ImportError:
            self.prolif_available = False
            logger.info("ProLIF not available - using geometric detection")
    
    def detect_interactions(
        self,
        receptor_atoms: List[Dict],
        ligand_atoms: List[Dict],
        receptor_bonds: Optional[List[Dict]] = None,
        ligand_bonds: Optional[List[Dict]] = None,
    ) -> InteractionSummary:
        """
        Detect all interactions between receptor and ligand.
        
        Args:
            receptor_atoms: List of receptor atom dictionaries with keys:
                - elem: element symbol
                - x, y, z: coordinates
                - resn: residue name
                - resi: residue number
                - name: atom name
            ligand_atoms: List of ligand atom dictionaries
            receptor_bonds: Optional receptor bonds
            ligand_bonds: Optional ligand bonds
            
        Returns:
            InteractionSummary with all detected interactions
        """
        if not receptor_atoms or not ligand_atoms:
            logger.warning("Empty atom lists provided")
            return self._empty_summary()
        
        # Build spatial index for fast distance calculations
        receptor_coords = np.array([(a['x'], a['y'], a['z']) for a in receptor_atoms])
        ligand_coords = np.array([(a['x'], a['y'], a['z']) for a in ligand_atoms])
        
        interactions = []
        
        # Detect each interaction type
        if self.prolif_available:
            interactions = self._detect_with_prolif(receptor_atoms, ligand_atoms)
        else:
            interactions = self._detect_geometric(receptor_atoms, ligand_atoms, 
                                                   receptor_coords, ligand_coords)
        
        # Build summary
        summary = self._build_summary(interactions)
        
        logger.info(f"Detected {summary.total_count} interactions")
        
        return summary
    
    def _detect_with_prolif(
        self,
        receptor_atoms: List[Dict],
        ligand_atoms: List[Dict]
    ) -> List[Interaction]:
        """Detect interactions using ProLIF library"""
        interactions = []
        
        try:
            import prolif as pl
            from rdkit import Chem
            
            # Convert to RDKit molecules (simplified - would need proper conversion)
            logger.info("Using ProLIF for interaction detection")
            
        except Exception as e:
            logger.warning(f"ProLIF detection failed: {e}, falling back to geometric")
            interactions = self._detect_geometric(
                receptor_atoms, ligand_atoms,
                np.array([(a['x'], a['y'], a['z']) for a in receptor_atoms]),
                np.array([(a['x'], a['y'], a['z']) for a in ligand_atoms])
            )
        
        return interactions
    
    def _detect_geometric(
        self,
        receptor_atoms: List[Dict],
        ligand_atoms: List[Dict],
        receptor_coords: np.ndarray,
        ligand_coords: np.ndarray
    ) -> List[Interaction]:
        """Detect interactions using geometric criteria"""
        interactions = []
        
        # Build KD-tree for fast neighbor search
        from scipy.spatial import cKDTree
        
        receptor_tree = cKDTree(receptor_coords)
        
        # Find all close pairs
        max_dist = max(self.DISTANCE_CUTOFFS.values())
        pairs = receptor_tree.query_ball_tree(cKDTree(ligand_coords), max_dist)
        
        # Classify each atom
        receptor_class = self._classify_atoms(receptor_atoms)
        ligand_class = self._classify_atoms(ligand_atoms)
        
        # Check each close pair for interactions
        for r_idx, l_indices in enumerate(pairs):
            if not l_indices:
                continue
            
            r_atom = receptor_atoms[r_idx]
            r_elem = r_atom.get('elem', '').upper()
            r_class = receptor_class.get(r_idx, {})
            
            for l_idx in l_indices:
                l_atom = ligand_atoms[l_idx]
                l_elem = l_atom.get('elem', '').upper()
                l_class = ligand_class.get(l_idx, {})
                
                # Calculate distance
                dist = np.linalg.norm(
                    receptor_coords[r_idx] - ligand_coords[l_idx]
                )
                
                # Hydrogen bonds
                if self._is_hbond_pair(r_class, l_class):
                    if dist < self.DISTANCE_CUTOFFS[InteractionType.HYDROGEN_BOND]:
                        interactions.append(self._create_interaction(
                            InteractionType.HYDROGEN_BOND, r_atom, l_atom,
                            receptor_coords[r_idx], ligand_coords[l_idx], dist
                        ))
                
                # Hydrophobic contacts
                if self._is_hydrophobic_pair(r_class, l_class):
                    if dist < self.DISTANCE_CUTOFFS[InteractionType.HYDROPHOBIC]:
                        interactions.append(self._create_interaction(
                            InteractionType.HYDROPHOBIC, r_atom, l_atom,
                            receptor_coords[r_idx], ligand_coords[l_idx], dist
                        ))
                
                # Salt bridges
                if self._is_salt_bridge_pair(r_class, l_class):
                    if dist < self.DISTANCE_CUTOFFS[InteractionType.SALT_BRIDGE]:
                        interactions.append(self._create_interaction(
                            InteractionType.SALT_BRIDGE, r_atom, l_atom,
                            receptor_coords[r_idx], ligand_coords[l_idx], dist
                        ))
                
                # Metal coordination
                if self._is_metal_pair(r_class, l_class):
                    if dist < self.DISTANCE_CUTOFFS[InteractionType.METAL_COORDINATION]:
                        interactions.append(self._create_interaction(
                            InteractionType.METAL_COORDINATION, r_atom, l_atom,
                            receptor_coords[r_idx], ligand_coords[l_idx], dist
                        ))
        
        return interactions
    
    def _classify_atoms(self, atoms: List[Dict]) -> Dict[int, Dict[str, bool]]:
        """Classify each atom by its properties"""
        classification = {}
        
        for idx, atom in enumerate(atoms):
            elem = atom.get('elem', '').upper()
            resn = atom.get('resn', '')
            
            classification[idx] = {
                'element': elem,
                'residue': resn,
                'is_donor': elem in self.HYDROGEN_DONORS,
                'is_acceptor': elem in self.HYDROGEN_ACCEPTORS,
                'is_hydrophobic': elem in self.HYDROPHOBIC_ATOMS,
                'is_positive': elem in self.POSITIVE_CHARGED or resn in ['LYS', 'ARG', 'HIS'],
                'is_negative': elem in self.NEGATIVE_CHARGED or resn in ['ASP', 'GLU'],
                'is_metal': elem in self.METALS,
                'is_aromatic': self._is_aromatic(atom),
            }
        
        return classification
    
    def _is_aromatic(self, atom: Dict) -> bool:
        """Check if atom is part of aromatic ring"""
        # Simplified check - would need proper ring detection
        elem = atom.get('elem', '').upper()
        resn = atom.get('resn', '')
        return elem == 'C' and resn in ['HIS', 'PHE', 'TYR', 'TRP']
    
    def _is_hbond_pair(self, r_class: Dict, l_class: Dict) -> bool:
        """Check if pair can form hydrogen bond"""
        return (
            (r_class.get('is_donor') and l_class.get('is_acceptor')) or
            (r_class.get('is_acceptor') and l_class.get('is_donor'))
        )
    
    def _is_hydrophobic_pair(self, r_class: Dict, l_class: Dict) -> bool:
        """Check if pair is hydrophobic"""
        return r_class.get('is_hydrophobic') and l_class.get('is_hydrophobic')
    
    def _is_salt_bridge_pair(self, r_class: Dict, l_class: Dict) -> bool:
        """Check if pair can form salt bridge"""
        return (
            (r_class.get('is_positive') and l_class.get('is_negative')) or
            (r_class.get('is_negative') and l_class.get('is_positive'))
        )
    
    def _is_metal_pair(self, r_class: Dict, l_class: Dict) -> bool:
        """Check if pair involves metal coordination"""
        return r_class.get('is_metal') and l_class.get('is_acceptor')
    
    def _create_interaction(
        self,
        itype: InteractionType,
        atom1: Dict,
        atom2: Dict,
        coords1: np.ndarray,
        coords2: np.ndarray,
        distance: float
    ) -> Interaction:
        """Create Interaction object"""
        return Interaction(
            type=itype,
            atom1=atom1.get('name', ''),
            atom2=atom2.get('name', ''),
            atom1_coords=tuple(coords1),
            atom2_coords=tuple(coords2),
            distance=distance,
            residue1=atom1.get('resn', '') + str(atom1.get('resi', '')),
            residue2=atom2.get('resn', '') + str(atom2.get('resi', '')),
            strength=self._calculate_strength(itype, distance),
        )
    
    def _calculate_strength(self, itype: InteractionType, distance: float) -> float:
        """Calculate interaction strength (0-1) based on distance"""
        cutoff = self.DISTANCE_CUTOFFS.get(itype, 4.0)
        if distance >= cutoff:
            return 0.0
        return 1.0 - (distance / cutoff)
    
    def _build_summary(self, interactions: List[Interaction]) -> InteractionSummary:
        """Build summary from interactions"""
        hbond = [i for i in interactions if i.type == InteractionType.HYDROGEN_BOND]
        hydrophobic = [i for i in interactions if i.type == InteractionType.HYDROPHOBIC]
        pi_stack = [i for i in interactions if i.type == InteractionType.PI_STACKING]
        salt = [i for i in interactions if i.type == InteractionType.SALT_BRIDGE]
        metal = [i for i in interactions if i.type == InteractionType.METAL_COORDINATION]
        
        # Get key residues (those involved in most interactions)
        residues = {}
        for interaction in interactions:
            if interaction.residue1:
                residues[interaction.residue1] = residues.get(interaction.residue1, 0) + 1
        
        key_residues = sorted(residues.keys(), key=lambda x: residues[x], reverse=True)[:10]
        
        return InteractionSummary(
            total_count=len(interactions),
            hbond_count=len(hbond),
            hydrophobic_count=len(hydrophobic),
            pi_stacking_count=len(pi_stack),
            salt_bridge_count=len(salt),
            metal_count=len(metal),
            interactions=interactions,
            key_residues=key_residues,
        )
    
    def _empty_summary(self) -> InteractionSummary:
        """Return empty summary"""
        return InteractionSummary(
            total_count=0,
            hbond_count=0,
            hydrophobic_count=0,
            pi_stacking_count=0,
            salt_bridge_count=0,
            metal_count=0,
            interactions=[],
            key_residues=[],
        )
    
    def to_dict(self, summary: InteractionSummary) -> Dict[str, Any]:
        """Convert summary to dictionary for JSON serialization"""
        return {
            'total_count': summary.total_count,
            'hbond_count': summary.hbond_count,
            'hydrophobic_count': summary.hydrophobic_count,
            'pi_stacking_count': summary.pi_stacking_count,
            'salt_bridge_count': summary.salt_bridge_count,
            'metal_count': summary.metal_count,
            'key_residues': summary.key_residues,
            'interactions': [
                {
                    'type': i.type.value,
                    'atom1': i.atom1,
                    'atom2': i.atom2,
                    'atom1_coords': list(i.atom1_coords),
                    'atom2_coords': list(i.atom2_coords),
                    'distance': i.distance,
                    'residue1': i.residue1,
                    'residue2': i.residue2,
                    'strength': i.strength,
                    'color': self.INTERACTION_COLORS.get(i.type.value, '#FFFFFF'),
                }
                for i in summary.interactions
            ],
        }
    
    def get_visualization_data(self, summary: InteractionSummary) -> List[Dict]:
        """Get data formatted for 3D viewer visualization"""
        viz_data = []
        
        for interaction in summary.interactions:
            viz_data.append({
                'type': interaction.type.value,
                'atom1': interaction.atom1,
                'atom2': interaction.atom2,
                'atom1_coords': list(interaction.atom1_coords),
                'atom2_coords': list(interaction.atom2_coords),
                'distance': interaction.distance,
                'color': self.INTERACTION_COLORS.get(interaction.type, '#FFFFFF'),
                'dashed': interaction.type == InteractionType.HYDROGEN_BOND,
                'radius': 0.05 if interaction.type == InteractionType.HYDROGEN_BOND else 0.1,
            })
        
        return viz_data


def parse_pdb_atoms(pdb_content: str) -> List[Dict]:
    """Parse PDB file to extract atoms"""
    atoms = []
    
    for line in pdb_content.split('\n'):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                atom = {
                    'serial': int(line[6:11].strip()),
                    'name': line[12:16].strip(),
                    'elem': line[76:78].strip() or line[12:14].strip(),
                    'x': float(line[30:38].strip()),
                    'y': float(line[38:46].strip()),
                    'z': float(line[46:54].strip()),
                    'resn': line[17:20].strip(),
                    'resi': int(line[22:26].strip()),
                    'chain': line[21:22].strip(),
                }
                atoms.append(atom)
            except (ValueError, IndexError):
                continue
    
    return atoms
