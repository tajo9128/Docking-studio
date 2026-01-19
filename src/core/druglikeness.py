"""
BioDockify Drug-likeness Calculator
Lipinski, Veber, Lead-likeness, Fragment-likeness rules
"""

from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
import logging

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, Crippen
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

logger = logging.getLogger(__name__)


@dataclass
class DrugLikenessResult:
    """Drug-likeness analysis result."""
    smiles: str
    descriptors: Dict[str, float] = field(default_factory=dict)
    lipinski: Dict[str, Any] = field(default_factory=dict)
    veber: Dict[str, Any] = field(default_factory=dict)
    ghose: Dict[str, Any] = field(default_factory=dict)
    lead_likeness: Dict[str, Any] = field(default_factory=dict)
    fragment_likeness: Dict[str, Any] = field(default_factory=dict)
    overall_druglikeness: str = 'unknown'
    score: float = 0.0
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            'smiles': self.smiles,
            'descriptors': self.descriptors,
            'lipinski': self.lipinski,
            'veber': self.veber,
            'ghose': self.ghose,
            'lead_likeness': self.lead_likeness,
            'fragment_likeness': self.fragment_likeness,
            'overall_druglikeness': self.overall_druglikeness,
            'score': self.score
        }


class DrugLikenessCalculator:
    """
    Drug-likeness Calculator.
    
    Evaluates molecules against multiple drug-likeness rules:
    - Lipinski's Rule of 5
    - Veber's Rules
    - Ghose Filter
    - Lead-likeness
    - Fragment-likeness (Rule of 3)
    """
    
    def __init__(self):
        """Initialize the calculator."""
        if not RDKIT_AVAILABLE:
            logger.warning("RDKit not available. Drug-likeness calculations limited.")
    
    def calculate(self, molecule) -> DrugLikenessResult:
        """
        Calculate drug-likeness for a molecule.
        
        Args:
            molecule: RDKit Mol object or SMILES string
            
        Returns:
            DrugLikenessResult with all analyses
        """
        if not RDKIT_AVAILABLE:
            return DrugLikenessResult(
                smiles=str(molecule) if isinstance(molecule, str) else "",
                overall_druglikeness='unknown'
            )
        
        # Convert SMILES to Mol if needed
        if isinstance(molecule, str):
            smiles = molecule
            mol = Chem.MolFromSmiles(molecule)
            if mol is None:
                return DrugLikenessResult(smiles=smiles, overall_druglikeness='invalid')
        else:
            mol = molecule
            smiles = Chem.MolToSmiles(mol)
        
        # Calculate descriptors
        descriptors = self._calculate_descriptors(mol)
        
        # Apply different filters
        lipinski = self._check_lipinski(descriptors)
        veber = self._check_veber(descriptors)
        ghose = self._check_ghose(descriptors)
        lead = self._check_lead_likeness(descriptors)
        fragment = self._check_fragment_likeness(descriptors)
        
        # Determine overall drug-likeness
        overall, score = self._evaluate_overall(lipinski, veber, ghose)
        
        return DrugLikenessResult(
            smiles=smiles,
            descriptors=descriptors,
            lipinski=lipinski,
            veber=veber,
            ghose=ghose,
            lead_likeness=lead,
            fragment_likeness=fragment,
            overall_druglikeness=overall,
            score=score
        )
    
    def _calculate_descriptors(self, mol) -> Dict[str, float]:
        """Calculate all relevant molecular descriptors."""
        return {
            'molecular_weight': Descriptors.MolWt(mol),
            'logp': Crippen.MolLogP(mol),
            'hba': Descriptors.NumHAcceptors(mol),
            'hbd': Descriptors.NumHDonors(mol),
            'tpsa': Descriptors.TPSA(mol),
            'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
            'rings': Descriptors.RingCount(mol),
            'aromatic_rings': Descriptors.NumAromaticRings(mol),
            'heavy_atoms': Descriptors.HeavyAtomCount(mol),
            'molar_refractivity': Crippen.MolMR(mol),
            'fraction_sp3': Descriptors.FractionCSP3(mol),
            'num_atoms': mol.GetNumAtoms(),
        }
    
    def _check_lipinski(self, desc: Dict) -> Dict[str, Any]:
        """
        Check Lipinski's Rule of 5.
        
        Rules:
        - MW <= 500
        - LogP <= 5
        - HBA <= 10
        - HBD <= 5
        """
        mw = desc['molecular_weight']
        logp = desc['logp']
        hba = desc['hba']
        hbd = desc['hbd']
        
        violations = []
        if mw > 500:
            violations.append(f"MW ({mw:.1f}) > 500")
        if logp > 5:
            violations.append(f"LogP ({logp:.2f}) > 5")
        if hba > 10:
            violations.append(f"HBA ({hba}) > 10")
        if hbd > 5:
            violations.append(f"HBD ({hbd}) > 5")
        
        return {
            'passed': len(violations) <= 1,  # Allow 1 violation
            'violations': len(violations),
            'details': violations,
            'mw_ok': mw <= 500,
            'logp_ok': logp <= 5,
            'hba_ok': hba <= 10,
            'hbd_ok': hbd <= 5,
        }
    
    def _check_veber(self, desc: Dict) -> Dict[str, Any]:
        """
        Check Veber's Rules for oral bioavailability.
        
        Rules:
        - Rotatable bonds <= 10
        - TPSA <= 140 Å²
        """
        rotbonds = desc['rotatable_bonds']
        tpsa = desc['tpsa']
        
        violations = []
        if rotbonds > 10:
            violations.append(f"Rotatable bonds ({rotbonds}) > 10")
        if tpsa > 140:
            violations.append(f"TPSA ({tpsa:.1f}) > 140")
        
        return {
            'passed': len(violations) == 0,
            'violations': len(violations),
            'details': violations,
            'rotbonds_ok': rotbonds <= 10,
            'tpsa_ok': tpsa <= 140,
        }
    
    def _check_ghose(self, desc: Dict) -> Dict[str, Any]:
        """
        Check Ghose Filter.
        
        Rules:
        - 160 <= MW <= 480
        - -0.4 <= LogP <= 5.6
        - 40 <= MR <= 130
        - 20 <= Atoms <= 70
        """
        mw = desc['molecular_weight']
        logp = desc['logp']
        mr = desc['molar_refractivity']
        atoms = desc['num_atoms']
        
        violations = []
        if not (160 <= mw <= 480):
            violations.append(f"MW ({mw:.1f}) not in [160, 480]")
        if not (-0.4 <= logp <= 5.6):
            violations.append(f"LogP ({logp:.2f}) not in [-0.4, 5.6]")
        if not (40 <= mr <= 130):
            violations.append(f"MR ({mr:.1f}) not in [40, 130]")
        if not (20 <= atoms <= 70):
            violations.append(f"Atoms ({atoms}) not in [20, 70]")
        
        return {
            'passed': len(violations) == 0,
            'violations': len(violations),
            'details': violations,
        }
    
    def _check_lead_likeness(self, desc: Dict) -> Dict[str, Any]:
        """
        Check Lead-likeness.
        
        Rules:
        - 250 <= MW <= 350
        - LogP <= 3.5
        - Rotatable bonds <= 7
        """
        mw = desc['molecular_weight']
        logp = desc['logp']
        rotbonds = desc['rotatable_bonds']
        
        violations = []
        if not (250 <= mw <= 350):
            violations.append(f"MW ({mw:.1f}) not in [250, 350]")
        if logp > 3.5:
            violations.append(f"LogP ({logp:.2f}) > 3.5")
        if rotbonds > 7:
            violations.append(f"Rotatable bonds ({rotbonds}) > 7")
        
        return {
            'passed': len(violations) == 0,
            'violations': len(violations),
            'details': violations,
            'suitable_for_optimization': len(violations) <= 1,
        }
    
    def _check_fragment_likeness(self, desc: Dict) -> Dict[str, Any]:
        """
        Check Fragment-likeness (Rule of 3).
        
        Rules:
        - MW <= 300
        - LogP <= 3
        - HBA <= 3
        - HBD <= 3
        - Rotatable bonds <= 3
        """
        mw = desc['molecular_weight']
        logp = desc['logp']
        hba = desc['hba']
        hbd = desc['hbd']
        rotbonds = desc['rotatable_bonds']
        
        violations = []
        if mw > 300:
            violations.append(f"MW ({mw:.1f}) > 300")
        if logp > 3:
            violations.append(f"LogP ({logp:.2f}) > 3")
        if hba > 3:
            violations.append(f"HBA ({hba}) > 3")
        if hbd > 3:
            violations.append(f"HBD ({hbd}) > 3")
        if rotbonds > 3:
            violations.append(f"Rotatable bonds ({rotbonds}) > 3")
        
        return {
            'passed': len(violations) == 0,
            'violations': len(violations),
            'details': violations,
            'is_fragment': len(violations) == 0,
        }
    
    def _evaluate_overall(self, lipinski: Dict, veber: Dict, 
                          ghose: Dict) -> tuple:
        """Evaluate overall drug-likeness."""
        score = 0.0
        
        # Lipinski weight: 40%
        if lipinski['passed']:
            score += 0.4
        elif lipinski['violations'] == 2:
            score += 0.2
        
        # Veber weight: 30%
        if veber['passed']:
            score += 0.3
        elif veber['violations'] == 1:
            score += 0.15
        
        # Ghose weight: 30%
        if ghose['passed']:
            score += 0.3
        elif ghose['violations'] <= 2:
            score += 0.15
        
        # Determine category
        if score >= 0.8:
            category = 'excellent'
        elif score >= 0.6:
            category = 'good'
        elif score >= 0.4:
            category = 'moderate'
        elif score >= 0.2:
            category = 'poor'
        else:
            category = 'very_poor'
        
        return category, score
    
    def calculate_smiles(self, smiles: str) -> DrugLikenessResult:
        """Convenience method for SMILES input."""
        return self.calculate(smiles)
    
    def calculate_batch(self, molecules: List) -> List[DrugLikenessResult]:
        """Calculate drug-likeness for multiple molecules."""
        return [self.calculate(mol) for mol in molecules]
