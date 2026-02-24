"""
Feature Extraction Module
Extracts interaction fingerprints using ODDT for ML-based scoring.
"""

import logging
from typing import Dict, Any, List, Optional
import numpy as np

logger = logging.getLogger(__name__)

oddt_available = False

try:
    import oddt
    from oddt import toolkit
    try:
        from oddt.fingerprints import InteractionFingerprint
    except ImportError:
        from oddt.scoring import InteractionFingerprint
    oddt_available = True
except ImportError:
    oddt_available = False
    logger.warning("ODDT not available - some features disabled")


class FeatureExtractor:
    """
    Extracts molecular features for ML-based scoring.
    Uses ODDT interaction fingerprints when available.
    """
    
    def __init__(self):
        self.oddt_available = oddt_available
    
    def extract_interaction_fingerprint(
        self,
        receptor_path: str,
        ligand_path: str
    ) -> Optional[np.ndarray]:
        """
        Extract 47-dimensional interaction fingerprint.
        
        Args:
            receptor_path: Path to receptor file (PDBQT)
            ligand_path: Path to ligand file (PDBQT)
            
        Returns:
            np.ndarray: 47-dim feature vector or None on failure
        """
        if not self.oddt_available:
            logger.warning("ODDT not available - returning placeholder features")
            return self._placeholder_features()
        
        try:
            rec = next(oddt.toolkit.readfile('pdbqt', receptor_path))
            lig = next(oddt.toolkit.readfile('pdbqt', ligand_path))
            rec.protein = True
            
            ifp = InteractionFingerprint(lig, rec)
            return np.array(ifp).flatten()
            
        except Exception as e:
            logger.error(f"Feature extraction failed: {e}")
            return self._placeholder_features()
    
    def _placeholder_features(self) -> np.ndarray:
        """Return placeholder features when ODDT unavailable"""
        return np.zeros(47)
    
    def build_feature_vector(
        self,
        vina_affinity: float,
        gnina_cnn_score: Optional[float] = None,
        gnina_cnn_affinity: Optional[float] = None,
        receptor_path: Optional[str] = None,
        ligand_path: Optional[str] = None
    ) -> np.ndarray:
        """
        Build complete feature vector for ML scoring.
        
        Features:
        - 5 Vina components (if available)
        - 2 GNINA scores (CNNscore, CNNaffinity)
        - 47 ODDT interaction fingerprint
        
        Total: 54 features
        
        Args:
            vina_affinity: Vina binding affinity (kcal/mol)
            gnina_cnn_score: GNINA CNN score (0-1)
            gnina_cnn_affinity: GNINA predicted affinity (pKd)
            receptor_path: Path to receptor for IFP
            ligand_path: Path to ligand for IFP
            
        Returns:
            np.ndarray: 54-dim feature vector
        """
        features = []
        
        features.append(vina_affinity)
        features.extend([0] * 4)
        
        if gnina_cnn_score is not None:
            features.append(gnina_cnn_score)
        else:
            features.append(0)
        
        if gnina_cnn_affinity is not None:
            features.append(gnina_cnn_affinity)
        else:
            features.append(0)
        
        if receptor_path and ligand_path:
            ifp = self.extract_interaction_fingerprint(receptor_path, ligand_path)
            if ifp is not None:
                features.extend(ifp.tolist())
            else:
                features.extend(self._placeholder_features().tolist())
        else:
            features.extend(self._placeholder_features().tolist())
        
        return np.array(features)
    
    def extract_batch_features(
        self,
        receptor_path: str,
        ligand_paths: List[str]
    ) -> List[np.ndarray]:
        """
        Extract features for multiple ligands.
        
        Args:
            receptor_path: Path to receptor
            ligand_paths: List of ligand file paths
            
        Returns:
            list: List of feature vectors
        """
        features_list = []
        
        for ligand_path in ligand_paths:
            features = self.extract_interaction_fingerprint(
                receptor_path,
                ligand_path
            )
            features_list.append(features)
        
        return features_list


def smiles_to_fingerprint(smiles: str, n_bits: int = 2048) -> Optional[np.ndarray]:
    """
    Convert SMILES to Morgan fingerprint (ECFP4).
    
    Args:
        smiles: SMILES string
        n_bits: Number of bits in fingerprint
        
    Returns:
        np.ndarray: Fingerprint or None on failure
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
        return np.array(fp)
        
    except Exception as e:
        logger.error(f"SMILES to fingerprint failed: {e}")
        return None
