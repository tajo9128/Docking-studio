"""
Pharmacophore Extractor
Automatic extraction and visualization of pharmacophore features.
"""

from enum import Enum
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple, Set
import numpy as np
import logging

logger = logging.getLogger(__name__)


class PharmacophoreFeature(Enum):
    """Types of pharmacophore features"""
    HYDROGEN_BOND_DONOR = "HBD"
    HYDROGEN_BOND_ACCEPTOR = "HBA"
    AROMATIC_RING = "AR"
    HYDROPHOBIC = "HYDRO"
    POSITIVE_CHARGE = "POS"
    NEGATIVE_CHARGE = "NEG"
    IONIZABLE = "ION"
    METAL = "METAL"


@dataclass
class PharmacophorePoint:
    """A single pharmacophore feature point"""
    feature_type: PharmacophoreFeature
    center: Tuple[float, float, float]
    atoms: List[str]
    residues: List[str]
    strength: float  # 0-1 scale
    enabled: bool = True


@dataclass
class PharmacophoreModel:
    """Complete pharmacophore model"""
    features: List[PharmacophorePoint]
    ligand_name: str
    created_from: str  # 'auto' or 'manual'


class PharmacophoreExtractor:
    """
    Extracts pharmacophore features from molecular structures.
    
    Features detected:
    - Hydrogen bond donors (HBD)
    - Hydrogen bond acceptors (HBA)
    - Aromatic rings (AR)
    - Hydrophobic regions
    - Charged groups (positive/negative)
    - Metal binding sites
    """
    
    # Feature colors for visualization (matching common standards)
    FEATURE_COLORS = {
        PharmacophoreFeature.HYDROGEN_BOND_DONOR: "#00FF00",    # Green
        PharmacophoreFeature.HYDROGEN_BOND_ACCEPTOR: "#0000FF", # Blue
        PharmacophoreFeature.AROMATIC_RING: "#800080",          # Purple
        PharmacophoreFeature.HYDROPHOBIC: "#FFFF00",           # Yellow
        PharmacophoreFeature.POSITIVE_CHARGE: "#FF0000",       # Red
        PharmacophoreFeature.NEGATIVE_CHARGE: "#FF6600",        # Orange
        PharmacophoreFeature.IONIZABLE: "#00FFFF",             # Cyan
        PharmacophoreFeature.METAL: "#808080",                 # Gray
    }
    
    # Feature shapes for visualization
    FEATURE_SHAPES = {
        PharmacophoreFeature.HYDROGEN_BOND_DONOR: "sphere",
        PharmacophoreFeature.HYDROGEN_BOND_ACCEPTOR: "sphere",
        PharmacophoreFeature.AROMATIC_RING: "ring",
        PharmacophoreFeature.HYDROPHOBIC: "sphere",
        PharmacophoreFeature.POSITIVE_CHARGE: "sphere",
        PharmacophoreFeature.NEGATIVE_CHARGE: "sphere",
        PharmacophoreFeature.IONIZABLE: "sphere",
        PharmacophoreFeature.METAL: "sphere",
    }
    
    # Feature radii
    FEATURE_RADII = {
        PharmacophoreFeature.HYDROGEN_BOND_DONOR: 0.5,
        PharmacophoreFeature.HYDROGEN_BOND_ACCEPTOR: 0.5,
        PharmacophoreFeature.AROMATIC_RING: 0.8,
        PharmacophoreFeature.HYDROPHOBIC: 0.6,
        PharmacophoreFeature.POSITIVE_CHARGE: 0.5,
        PharmacophoreFeature.NEGATIVE_CHARGE: 0.5,
        PharmacophoreFeature.IONIZABLE: 0.5,
        PharmacophoreFeature.METAL: 0.6,
    }
    
    def __init__(self):
        """Initialize pharmacophore extractor"""
        self.current_model: Optional[PharmacophoreModel] = None
        logger.info("PharmacophoreExtractor initialized")
    
    def extract_from_ligand(
        self,
        atoms: List[Dict],
        bonds: Optional[List[Dict]] = None,
        name: str = "ligand"
    ) -> PharmacophoreModel:
        """
        Extract pharmacophore features from ligand.
        
        Args:
            atoms: List of atom dictionaries
            bonds: Optional bond list
            name: Name for the pharmacophore model
            
        Returns:
            PharmacophoreModel with extracted features
        """
        features = []
        
        # Classify atoms
        atom_classes = self._classify_atoms(atoms)
        
        # Extract each feature type
        features.extend(self._find_hbd(atoms, atom_classes))
        features.extend(self._find_hba(atoms, atom_classes))
        features.extend(self._find_aromatic(atoms, atom_classes))
        features.extend(self._find_hydrophobic(atoms, atom_classes))
        features.extend(self._find_charged(atoms, atom_classes))
        features.extend(self._find_metal(atoms, atom_classes))
        
        model = PharmacophoreModel(
            features=features,
            ligand_name=name,
            created_from="auto"
        )
        
        self.current_model = model
        
        logger.info(f"Extracted {len(features)} pharmacophore features from {name}")
        
        return model
    
    def _classify_atoms(self, atoms: List[Dict]) -> Dict[int, Dict[str, bool]]:
        """Classify atoms by their properties"""
        classification = {}
        
        hbd_elements = {'N', 'O', 'S'}  # Can be H-bond donors
        hba_elements = {'N', 'O', 'S', 'F', 'CL', 'BR', 'I'}  # H-bond acceptors
        hydrophobic = {'C', 'S'}
        positive = {'N', 'K', 'R', 'H', 'NA', 'CA'}  # Simplified
        negative = {'O', 'D', 'E', 'CL', 'BR', 'F', 'I', 'S'}
        metals = {'MG', 'ZN', 'FE', 'CA', 'MN', 'CO', 'NI', 'CU', 'CD'}
        
        aromatic_residues = {'HIS', 'PHE', 'TYR', 'TRP'}
        
        for idx, atom in enumerate(atoms):
            elem = atom.get('elem', '').upper()
            resn = atom.get('resn', '')
            
            classification[idx] = {
                'is_hbd': elem in hbd_elements and atom.get('hydrogens', 0) > 0,
                'is_hba': elem in hba_elements,
                'is_hydrophobic': elem in hydrophobic,
                'is_positive': elem in positive or resn in ['LYS', 'ARG', 'HIS'],
                'is_negative': elem in negative or resn in ['ASP', 'GLU'],
                'is_metal': elem in metals,
                'is_aromatic': resn in aromatic_residues,
                'element': elem,
                'residue': resn,
                'name': atom.get('name', ''),
            }
        
        return classification
    
    def _find_hbd(
        self,
        atoms: List[Dict],
        classes: Dict[int, Dict]
    ) -> List[PharmacophorePoint]:
        """Find hydrogen bond donors"""
        hbd_atoms = [
            i for i, c in classes.items()
            if c.get('is_hbd', False) and c.get('hydrogens', 0) > 0
        ]
        
        features = []
        for idx in hbd_atoms:
            atom = atoms[idx]
            # Find attached hydrogens
            hydrogens = [
                a for a in atoms
                if a.get('name', '').startswith('H') and
                self._is_connected(atom, a, atoms)
            ]
            
            # Use hydrogen positions as feature points
            for h in hydrogens:
                features.append(PharmacophorePoint(
                    feature_type=PharmacophoreFeature.HYDROGEN_BOND_DONOR,
                    center=(h['x'], h['y'], h['z']),
                    atoms=[atom.get('name', ''), h.get('name', '')],
                    residues=[atom.get('resn', '') + str(atom.get('resi', ''))],
                    strength=0.8,
                ))
            
            # If no hydrogens, use heavy atom
            if not hydrogens:
                features.append(PharmacophorePoint(
                    feature_type=PharmacophoreFeature.HYDROGEN_BOND_DONOR,
                    center=(atom['x'], atom['y'], atom['z']),
                    atoms=[atom.get('name', '')],
                    residues=[atom.get('resn', '') + str(atom.get('resi', ''))],
                    strength=0.6,
                ))
        
        return features
    
    def _find_hba(
        self,
        atoms: List[Dict],
        classes: Dict[int, Dict]
    ) -> List[PharmacophorePoint]:
        """Find hydrogen bond acceptors"""
        hba_atoms = [
            i for i, c in classes.items()
            if c.get('is_hba', False)
        ]
        
        features = []
        for idx in hba_atoms:
            atom = atoms[idx]
            features.append(PharmacophorePoint(
                feature_type=PharmacophoreFeature.HYDROGEN_BOND_ACCEPTOR,
                center=(atom['x'], atom['y'], atom['z']),
                atoms=[atom.get('name', '')],
                residues=[atom.get('resn', '') + str(atom.get('resi', ''))],
                strength=0.8,
            ))
        
        return features
    
    def _find_aromatic(
        self,
        atoms: List[Dict],
        classes: Dict[int, Dict]
    ) -> List[PharmacophorePoint]:
        """Find aromatic rings"""
        aromatic_atoms = [
            i for i, c in classes.items()
            if c.get('is_aromatic', False)
        ]
        
        if not aromatic_atoms:
            return []
        
        # Group by residue
        by_residue: Dict[str, List] = {}
        for idx in aromatic_atoms:
            atom = atoms[idx]
            res = atom.get('resn', '') + str(atom.get('resi', ''))
            if res not in by_residue:
                by_residue[res] = []
            by_residue[res].append(idx)
        
        features = []
        for res, indices in by_residue.items():
            # Calculate centroid
            coords = np.array([
                (atoms[i]['x'], atoms[i]['y'], atoms[i]['z'])
                for i in indices
            ])
            centroid = coords.mean(axis=0)
            
            features.append(PharmacophorePoint(
                feature_type=PharmacophoreFeature.AROMATIC_RING,
                center=tuple(centroid),
                atoms=[atoms[i].get('name', '') for i in indices],
                residues=[res],
                strength=0.9,
            ))
        
        return features
    
    def _find_hydrophobic(
        self,
        atoms: List[Dict],
        classes: Dict[int, Dict]
    ) -> List[PharmacophorePoint]:
        """Find hydrophobic regions"""
        hydro_atoms = [
            i for i, c in classes.items()
            if c.get('is_hydrophobic', False)
        ]
        
        # Group nearby hydrophobic atoms
        clusters = self._cluster_atoms(atoms, hydro_atoms, cutoff=2.0)
        
        features = []
        for cluster in clusters:
            coords = np.array([
                (atoms[i]['x'], atoms[i]['y'], atoms[i]['z'])
                for i in cluster
            ])
            centroid = coords.mean(axis=0)
            
            features.append(PharmacophorePoint(
                feature_type=PharmacophoreFeature.HYDROPHOBIC,
                center=tuple(centroid),
                atoms=[atoms[i].get('name', '') for i in cluster],
                residues=list(set(
                    atoms[i].get('resn', '') + str(atoms[i].get('resi', ''))
                    for i in cluster
                )),
                strength=min(0.5 + 0.1 * len(cluster), 1.0),
            ))
        
        return features
    
    def _find_charged(
        self,
        atoms: List[Dict],
        classes: Dict[int, Dict]
    ) -> List[PharmacophorePoint]:
        """Find charged groups"""
        features = []
        
        # Positive charges
        pos_atoms = [
            i for i, c in classes.items()
            if c.get('is_positive', False)
        ]
        
        for idx in pos_atoms:
            atom = atoms[idx]
            features.append(PharmacophorePoint(
                feature_type=PharmacophoreFeature.POSITIVE_CHARGE,
                center=(atom['x'], atom['y'], atom['z']),
                atoms=[atom.get('name', '')],
                residues=[atom.get('resn', '') + str(atom.get('resi', ''))],
                strength=0.9,
            ))
        
        # Negative charges
        neg_atoms = [
            i for i, c in classes.items()
            if c.get('is_negative', False)
        ]
        
        for idx in neg_atoms:
            atom = atoms[idx]
            features.append(PharmacophorePoint(
                feature_type=PharmacophoreFeature.NEGATIVE_CHARGE,
                center=(atom['x'], atom['y'], atom['z']),
                atoms=[atom.get('name', '')],
                residues=[atom.get('resn', '') + str(atom.get('resi', ''))],
                strength=0.9,
            ))
        
        return features
    
    def _find_metal(
        self,
        atoms: List[Dict],
        classes: Dict[int, Dict]
    ) -> List[PharmacophorePoint]:
        """Find metal binding sites"""
        metal_atoms = [
            i for i, c in classes.items()
            if c.get('is_metal', False)
        ]
        
        features = []
        for idx in metal_atoms:
            atom = atoms[idx]
            features.append(PharmacophorePoint(
                feature_type=PharmacophoreFeature.METAL,
                center=(atom['x'], atom['y'], atom['z']),
                atoms=[atom.get('name', '')],
                residues=[atom.get('resn', '') + str(atom.get('resi', ''))],
                strength=1.0,
            ))
        
        return features
    
    def _is_connected(self, atom1: Dict, atom2: Dict, all_atoms: List[Dict]) -> bool:
        """Check if two atoms are connected (simplified)"""
        dist = np.sqrt(
            (atom1['x'] - atom2['x'])**2 +
            (atom1['y'] - atom2['y'])**2 +
            (atom1['z'] - atom2['z'])**2
        )
        return dist < 1.5
    
    def _cluster_atoms(
        self,
        atoms: List[Dict],
        indices: List[int],
        cutoff: float = 2.0
    ) -> List[List[int]]:
        """Cluster atoms by distance"""
        if not indices:
            return []
        
        coords = np.array([
            (atoms[i]['x'], atoms[i]['y'], atoms[i]['z'])
            for i in indices
        ])
        
        # Simple clustering
        clusters = []
        used = set()
        
        for i, idx in enumerate(indices):
            if idx in used:
                continue
            
            cluster = [idx]
            used.add(idx)
            
            for j, other_idx in enumerate(indices[i+1:], i+1):
                if other_idx in used:
                    continue
                
                dist = np.linalg.norm(coords[i] - coords[j])
                if dist < cutoff:
                    cluster.append(other_idx)
                    used.add(other_idx)
            
            clusters.append(cluster)
        
        return clusters
    
    def get_visualization_data(self, model: Optional[PharmacophoreModel] = None) -> List[Dict]:
        """Get data for 3D viewer visualization"""
        if model is None:
            model = self.current_model
        
        if model is None:
            return []
        
        viz_data = []
        
        for feature in model.features:
            viz_data.append({
                'type': feature.feature_type.value,
                'center': list(feature.center),
                'color': self.FEATURE_COLORS.get(feature.feature_type, '#808080'),
                'shape': self.FEATURE_SHAPES.get(feature.feature_type, 'sphere'),
                'radius': self.FEATURE_RADII.get(feature.feature_type, 0.5),
                'strength': feature.strength,
                'atoms': feature.atoms,
                'residues': feature.residues,
            })
        
        return viz_data
    
    def get_2d_coordinates(self, model: Optional[PharmacophoreModel] = None) -> Dict:
        """Get 2D coordinates for RDKit drawing"""
        if model is None:
            model = self.current_model
        
        if model is None:
            return {}
        
        # Simplified 2D mapping
        coords = {}
        for i, feature in enumerate(model.features):
            angle = 2 * np.pi * i / len(model.features)
            r = 2.0
            coords[feature.feature_type.value] = (r * np.cos(angle), r * np.sin(angle))
        
        return coords
    
    def filter_by_type(
        self,
        feature_types: Set[PharmacophoreFeature],
        model: Optional[PharmacophoreModel] = None
    ) -> List[PharmacophorePoint]:
        """Filter features by type"""
        if model is None:
            model = self.current_model
        
        if model is None:
            return []
        
        return [f for f in model.features if f.feature_type in feature_types]
    
    def get_feature_summary(self, model: Optional[PharmacophoreModel] = None) -> Dict:
        """Get summary of features"""
        if model is None:
            model = self.current_model
        
        if model is None:
            return {}
        
        summary = {}
        for feature in model.features:
            ftype = feature.feature_type.value
            if ftype not in summary:
                summary[ftype] = 0
            summary[ftype] += 1
        
        return summary
