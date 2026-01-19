"""
BioDockify ADMET Predictor
Absorption, Distribution, Metabolism, Excretion, Toxicity Predictions
"""

from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
import logging
import math

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, AllChem, Crippen
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

logger = logging.getLogger(__name__)


@dataclass
class ADMETResult:
    """Complete ADMET prediction result."""
    absorption: Dict[str, Any] = field(default_factory=dict)
    distribution: Dict[str, Any] = field(default_factory=dict)
    metabolism: Dict[str, Any] = field(default_factory=dict)
    excretion: Dict[str, Any] = field(default_factory=dict)
    toxicity: Dict[str, Any] = field(default_factory=dict)
    overall_score: float = 0.0
    warnings: List[str] = field(default_factory=list)
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            'absorption': self.absorption,
            'distribution': self.distribution,
            'metabolism': self.metabolism,
            'excretion': self.excretion,
            'toxicity': self.toxicity,
            'overall_score': self.overall_score,
            'warnings': self.warnings
        }


class ADMETPredictor:
    """
    ADMET Property Predictor.
    
    Predicts pharmacokinetic and toxicity properties of molecules
    based on molecular descriptors and empirical rules.
    """
    
    def __init__(self):
        """Initialize the ADMET predictor."""
        if not RDKIT_AVAILABLE:
            logger.warning("RDKit not available. ADMET predictions will be limited.")
    
    def predict(self, molecule) -> ADMETResult:
        """
        Predict ADMET properties for a molecule.
        
        Args:
            molecule: RDKit Mol object or SMILES string
            
        Returns:
            ADMETResult with all predictions
        """
        if not RDKIT_AVAILABLE:
            return ADMETResult(warnings=["RDKit not available"])
        
        # Convert SMILES to Mol if needed
        if isinstance(molecule, str):
            mol = Chem.MolFromSmiles(molecule)
            if mol is None:
                return ADMETResult(warnings=["Invalid SMILES string"])
        else:
            mol = molecule
        
        # Calculate descriptors
        descriptors = self._calculate_descriptors(mol)
        
        # Make predictions
        absorption = self._predict_absorption(mol, descriptors)
        distribution = self._predict_distribution(mol, descriptors)
        metabolism = self._predict_metabolism(mol, descriptors)
        excretion = self._predict_excretion(mol, descriptors)
        toxicity = self._predict_toxicity(mol, descriptors)
        
        # Calculate overall score
        overall = self._calculate_overall_score(
            absorption, distribution, metabolism, excretion, toxicity
        )
        
        # Collect warnings
        warnings = []
        if absorption.get('caco2_permeability', '') == 'low':
            warnings.append("Low intestinal absorption predicted")
        if distribution.get('bbb_penetration') == 'no':
            warnings.append("Poor blood-brain barrier penetration")
        if toxicity.get('ames_mutagenicity') == 'positive':
            warnings.append("Potential Ames mutagenicity")
        if toxicity.get('herg_inhibition') == 'high':
            warnings.append("High hERG inhibition risk - cardiotoxicity concern")
        
        return ADMETResult(
            absorption=absorption,
            distribution=distribution,
            metabolism=metabolism,
            excretion=excretion,
            toxicity=toxicity,
            overall_score=overall,
            warnings=warnings
        )
    
    def _calculate_descriptors(self, mol) -> Dict[str, float]:
        """Calculate molecular descriptors."""
        return {
            'mw': Descriptors.MolWt(mol),
            'logp': Crippen.MolLogP(mol),
            'hba': Descriptors.NumHAcceptors(mol),
            'hbd': Descriptors.NumHDonors(mol),
            'tpsa': Descriptors.TPSA(mol),
            'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
            'rings': Descriptors.RingCount(mol),
            'aromatic_rings': Descriptors.NumAromaticRings(mol),
            'heavy_atoms': Descriptors.HeavyAtomCount(mol),
            'fraction_sp3': Descriptors.FractionCSP3(mol),
            'molar_refractivity': Crippen.MolMR(mol),
        }
    
    def _predict_absorption(self, mol, desc: Dict) -> Dict[str, Any]:
        """Predict absorption properties."""
        # Caco-2 permeability estimation
        # Based on TPSA and LogP
        tpsa = desc['tpsa']
        logp = desc['logp']
        
        if tpsa < 60 and 1 < logp < 4:
            caco2 = 'high'
            caco2_score = 0.9
        elif tpsa < 100 and 0 < logp < 5:
            caco2 = 'medium'
            caco2_score = 0.6
        else:
            caco2 = 'low'
            caco2_score = 0.3
        
        # Human Intestinal Absorption
        # Based on Lipinski and TPSA
        mw = desc['mw']
        hba = desc['hba']
        hbd = desc['hbd']
        
        if mw < 500 and logp < 5 and hba < 10 and hbd < 5 and tpsa < 140:
            hia = 'high (>80%)'
            hia_score = 0.85
        elif tpsa < 200:
            hia = 'medium (50-80%)'
            hia_score = 0.65
        else:
            hia = 'low (<50%)'
            hia_score = 0.3
        
        # Oral bioavailability
        # Rule of 5 violations
        violations = 0
        if mw > 500: violations += 1
        if logp > 5: violations += 1
        if hba > 10: violations += 1
        if hbd > 5: violations += 1
        
        if violations == 0:
            bioavailability = 'good'
            bio_score = 0.8
        elif violations == 1:
            bioavailability = 'moderate'
            bio_score = 0.5
        else:
            bioavailability = 'poor'
            bio_score = 0.2
        
        return {
            'caco2_permeability': caco2,
            'caco2_score': caco2_score,
            'human_intestinal_absorption': hia,
            'hia_score': hia_score,
            'oral_bioavailability': bioavailability,
            'bioavailability_score': bio_score,
            'lipinski_violations': violations,
        }
    
    def _predict_distribution(self, mol, desc: Dict) -> Dict[str, Any]:
        """Predict distribution properties."""
        tpsa = desc['tpsa']
        logp = desc['logp']
        mw = desc['mw']
        
        # Blood-Brain Barrier penetration
        # Based on TPSA and MW
        if tpsa < 70 and mw < 450 and logp > 1:
            bbb = 'yes'
            bbb_score = 0.8
        elif tpsa < 90 and mw < 500:
            bbb = 'possible'
            bbb_score = 0.5
        else:
            bbb = 'no'
            bbb_score = 0.2
        
        # Plasma Protein Binding
        if logp > 4:
            ppb = 'high (>90%)'
            ppb_score = 0.3  # High binding is often undesirable
        elif logp > 2:
            ppb = 'medium (50-90%)'
            ppb_score = 0.6
        else:
            ppb = 'low (<50%)'
            ppb_score = 0.8
        
        # Volume of distribution estimation
        if logp > 3:
            vd = 'high'
        elif logp > 0:
            vd = 'medium'
        else:
            vd = 'low'
        
        return {
            'bbb_penetration': bbb,
            'bbb_score': bbb_score,
            'plasma_protein_binding': ppb,
            'ppb_score': ppb_score,
            'volume_of_distribution': vd,
        }
    
    def _predict_metabolism(self, mol, desc: Dict) -> Dict[str, Any]:
        """Predict metabolism properties."""
        mw = desc['mw']
        logp = desc['logp']
        
        # CYP inhibition likelihood (simplified)
        aromatic_rings = desc['aromatic_rings']
        
        cyp_risk = {}
        for cyp in ['CYP1A2', 'CYP2C9', 'CYP2C19', 'CYP2D6', 'CYP3A4']:
            # Simplified heuristic
            if aromatic_rings > 2 and logp > 3:
                cyp_risk[cyp] = 'high'
            elif aromatic_rings > 1:
                cyp_risk[cyp] = 'medium'
            else:
                cyp_risk[cyp] = 'low'
        
        # Metabolic stability
        if desc['fraction_sp3'] > 0.4 and desc['rotatable_bonds'] < 5:
            stability = 'high'
            stab_score = 0.8
        elif desc['rotatable_bonds'] < 10:
            stability = 'medium'
            stab_score = 0.5
        else:
            stability = 'low'
            stab_score = 0.3
        
        return {
            'cyp_inhibition': cyp_risk,
            'metabolic_stability': stability,
            'stability_score': stab_score,
        }
    
    def _predict_excretion(self, mol, desc: Dict) -> Dict[str, Any]:
        """Predict excretion properties."""
        mw = desc['mw']
        logp = desc['logp']
        tpsa = desc['tpsa']
        
        # Renal clearance
        if mw < 350 and logp < 2 and tpsa > 50:
            renal = 'high'
            renal_score = 0.8
        elif mw < 500:
            renal = 'medium'
            renal_score = 0.5
        else:
            renal = 'low'
            renal_score = 0.3
        
        # Half-life estimation (very rough)
        if logp > 4 and mw > 500:
            half_life = 'long (>12h)'
        elif logp > 2:
            half_life = 'medium (4-12h)'
        else:
            half_life = 'short (<4h)'
        
        return {
            'renal_clearance': renal,
            'renal_score': renal_score,
            'half_life_estimate': half_life,
        }
    
    def _predict_toxicity(self, mol, desc: Dict) -> Dict[str, Any]:
        """Predict toxicity properties."""
        logp = desc['logp']
        mw = desc['mw']
        tpsa = desc['tpsa']
        aromatic = desc['aromatic_rings']
        
        # hERG inhibition (cardiotoxicity)
        if logp > 4 and aromatic > 1:
            herg = 'high'
            herg_score = 0.2
        elif logp > 3:
            herg = 'medium'
            herg_score = 0.5
        else:
            herg = 'low'
            herg_score = 0.8
        
        # Ames mutagenicity (simplified)
        # Check for nitro groups, anilines, etc.
        smiles = Chem.MolToSmiles(mol)
        if '[N+](=O)[O-]' in smiles or 'c1ccc(N)cc1' in smiles:
            ames = 'positive'
            ames_score = 0.2
        else:
            ames = 'negative'
            ames_score = 0.8
        
        # Hepatotoxicity
        if mw > 600 and logp > 4:
            hepato = 'high'
            hepato_score = 0.3
        elif mw > 500:
            hepato = 'medium'
            hepato_score = 0.5
        else:
            hepato = 'low'
            hepato_score = 0.8
        
        # Acute oral toxicity (LD50 class)
        if logp > 5:
            acute = 'III (slightly toxic)'
        elif logp > 3:
            acute = 'IV (practically non-toxic)'
        else:
            acute = 'V (non-toxic)'
        
        return {
            'herg_inhibition': herg,
            'herg_score': herg_score,
            'ames_mutagenicity': ames,
            'ames_score': ames_score,
            'hepatotoxicity': hepato,
            'hepato_score': hepato_score,
            'acute_toxicity_class': acute,
        }
    
    def _calculate_overall_score(self, absorption, distribution, metabolism, 
                                  excretion, toxicity) -> float:
        """Calculate overall ADMET score (0-1)."""
        scores = []
        
        scores.append(absorption.get('caco2_score', 0.5))
        scores.append(absorption.get('hia_score', 0.5))
        scores.append(absorption.get('bioavailability_score', 0.5))
        scores.append(distribution.get('bbb_score', 0.5))
        scores.append(distribution.get('ppb_score', 0.5))
        scores.append(metabolism.get('stability_score', 0.5))
        scores.append(excretion.get('renal_score', 0.5))
        scores.append(toxicity.get('herg_score', 0.5))
        scores.append(toxicity.get('ames_score', 0.5))
        scores.append(toxicity.get('hepato_score', 0.5))
        
        return sum(scores) / len(scores)
    
    def predict_smiles(self, smiles: str) -> ADMETResult:
        """Convenience method to predict from SMILES."""
        return self.predict(smiles)
    
    def predict_batch(self, molecules: List) -> List[ADMETResult]:
        """Predict ADMET for multiple molecules."""
        return [self.predict(mol) for mol in molecules]
