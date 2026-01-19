"""
BioDockify PAINS Filter
Pan-Assay Interference Compounds (PAINS) detection using SMARTS patterns
"""

from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
import logging

try:
    from rdkit import Chem
    from rdkit.Chem import FilterCatalog
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

logger = logging.getLogger(__name__)


@dataclass
class PAINSResult:
    """Result of PAINS filter analysis."""
    passed: bool
    alerts: List[str]
    matched_patterns: List[str]
    severity: str  # 'none', 'low', 'medium', 'high'


class PAINSFilter:
    """
    PAINS (Pan-Assay Interference Compounds) Filter.
    
    Detects compounds that frequently show false positive activity 
    in high-throughput screens due to assay interference.
    """
    
    # PAINS SMARTS patterns (subset of common patterns)
    PAINS_SMARTS = {
        # Quinones
        'quinone_A': '[#6]1([#6]=[#6][#6](=[O])[#6]=[#6]1=[O])',
        
        # Catechols
        'catechol_A': 'c1c([OH])c([OH])ccc1',
        
        # Rhodanines
        'rhodanine_A': '[#6]1(=[#16])[#7][#6](=[O])[#16][#6]1',
        
        # Michael acceptors
        'michael_acceptor_1': '[#6]=[#6][#6](=[O])[!#7]',
        
        # Aldehydes
        'aldehyde': '[CH1](=O)',
        
        # Thiols
        'thiol': '[SH]',
        
        # Alpha-halo carbonyl
        'alpha_halo_carbonyl': '[F,Cl,Br,I][CH2][C](=[O])',
        
        # Aliphatic nitro
        'aliphatic_nitro': '[CH2][N+](=[O])[O-]',
        
        # Azo compounds
        'azo': '[N]=[N]',
        
        # Peroxides
        'peroxide': '[O][O]',
        
        # Acyl hydrazines
        'acyl_hydrazine': '[#6](=[O])[#7][#7]',
        
        # Sulfonyl halides
        'sulfonyl_halide': '[S](=[O])(=[O])[F,Cl,Br,I]',
        
        # Isocyanates
        'isocyanate': '[N]=[C]=[O]',
        
        # Thioureas
        'thiourea': '[#7][#6](=[S])[#7]',
        
        # Hydroxamic acids
        'hydroxamic_acid': '[#6](=[O])[#7][OH]',
        
        # Enones
        'enone': '[#6]=[#6][#6](=[O])[#6,#1]',
        
        # Activated esters
        'activated_ester': '[#6](=[O])[O][c]1[n,o,s]',
    }
    
    def __init__(self):
        """Initialize the PAINS filter."""
        self.patterns = {}
        self._compile_patterns()
    
    def _compile_patterns(self):
        """Compile SMARTS patterns."""
        if not RDKIT_AVAILABLE:
            logger.warning("RDKit not available. PAINS filtering disabled.")
            return
        
        for name, smarts in self.PAINS_SMARTS.items():
            try:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern:
                    self.patterns[name] = pattern
            except Exception as e:
                logger.warning(f"Failed to compile SMARTS pattern {name}: {e}")
    
    def filter(self, molecule) -> PAINSResult:
        """
        Filter a molecule for PAINS patterns.
        
        Args:
            molecule: RDKit Mol object or SMILES string
            
        Returns:
            PAINSResult with pass/fail status and matched patterns
        """
        if not RDKIT_AVAILABLE:
            return PAINSResult(
                passed=True,
                alerts=["RDKit not available - PAINS check skipped"],
                matched_patterns=[],
                severity='none'
            )
        
        # Convert SMILES to Mol if needed
        if isinstance(molecule, str):
            mol = Chem.MolFromSmiles(molecule)
            if mol is None:
                return PAINSResult(
                    passed=False,
                    alerts=["Invalid SMILES string"],
                    matched_patterns=[],
                    severity='high'
                )
        else:
            mol = molecule
        
        # Check for PAINS patterns
        matched = []
        alerts = []
        
        for name, pattern in self.patterns.items():
            if mol.HasSubstructMatch(pattern):
                matched.append(name)
                alerts.append(f"PAINS alert: {name.replace('_', ' ').title()}")
        
        # Determine severity
        if len(matched) == 0:
            severity = 'none'
            passed = True
        elif len(matched) <= 1:
            severity = 'low'
            passed = True  # Borderline - warn but pass
        elif len(matched) <= 3:
            severity = 'medium'
            passed = False
        else:
            severity = 'high'
            passed = False
        
        return PAINSResult(
            passed=passed,
            alerts=alerts,
            matched_patterns=matched,
            severity=severity
        )
    
    def filter_smiles(self, smiles: str) -> PAINSResult:
        """
        Filter a SMILES string for PAINS patterns.
        
        Args:
            smiles: SMILES string
            
        Returns:
            PAINSResult
        """
        return self.filter(smiles)
    
    def filter_batch(self, molecules: List) -> List[PAINSResult]:
        """
        Filter multiple molecules for PAINS patterns.
        
        Args:
            molecules: List of RDKit Mol objects or SMILES strings
            
        Returns:
            List of PAINSResult objects
        """
        return [self.filter(mol) for mol in molecules]
    
    def get_pattern_description(self, pattern_name: str) -> str:
        """Get a description for a PAINS pattern."""
        descriptions = {
            'quinone_A': 'Quinone - redox cycling, protein reactivity',
            'catechol_A': 'Catechol - oxidation, metal chelation',
            'rhodanine_A': 'Rhodanine - promiscuous binding',
            'michael_acceptor_1': 'Michael acceptor - covalent protein modification',
            'aldehyde': 'Aldehyde - Schiff base formation',
            'thiol': 'Thiol - disulfide formation, redox',
            'alpha_halo_carbonyl': 'Alpha-halo carbonyl - alkylation',
            'aliphatic_nitro': 'Aliphatic nitro - redox, toxicity',
            'azo': 'Azo compound - non-specific binding',
            'peroxide': 'Peroxide - oxidative damage',
            'acyl_hydrazine': 'Acyl hydrazine - reactive',
            'sulfonyl_halide': 'Sulfonyl halide - electrophilic',
            'isocyanate': 'Isocyanate - protein modification',
            'thiourea': 'Thiourea - non-specific binding',
            'hydroxamic_acid': 'Hydroxamic acid - metal chelation',
            'enone': 'Enone - Michael acceptor',
            'activated_ester': 'Activated ester - acylation',
        }
        return descriptions.get(pattern_name, 'Unknown pattern')
