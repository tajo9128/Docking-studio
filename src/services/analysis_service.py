"""Analysis Service - Orchestrates Molecular Analysis"""

from typing import Optional, Dict, List, Any
from datetime import datetime
import logging

logger = logging.getLogger(__name__)


class AnalysisService:
    """Analysis service for molecular interactions"""
    
    def __init__(self):
        self.engines_available = False
        try:
            from ..core.engines.molecular_engine import MolecularEngine
            from ..core.engines.interaction_pipeline import InteractionPipeline
            self.molecular_engine = MolecularEngine()
            self.interaction_pipeline = InteractionPipeline()
            self.engines_available = True
        except ImportError as e:
            logger.warning(f"Molecular engines not available: {e}")
            self.molecular_engine = None
            self.interaction_pipeline = None
    
    def analyze_interactions(self, atoms: List[Dict], bonds: Optional[List[Dict]] = None) -> Dict[str, Any]:
        """
        Analyze molecular interactions (H-bonds, VdW, Salt Bridges)
        
        Args:
            atoms: List of atom dictionaries with x, y, z coordinates
            bonds: Optional list of bonds
            
        Returns:
            Dictionary with interaction analysis results
        """
        logger.info(f"Analyzing interactions for {len(atoms)} atoms")
        
        start_time = datetime.now()
        
        if not atoms or len(atoms) < 2:
            return {
                'hydrogen_bonds': [],
                'vdw_contacts': [],
                'salt_bridges': [],
                'total_interactions': 0,
                'processing_time_ms': 0,
                'atom_count': len(atoms) if atoms else 0,
            }
        
        if not self.engines_available:
            logger.warning("No molecular engines available, returning empty results")
            return {
                'hydrogen_bonds': [],
                'vdw_contacts': [],
                'salt_bridges': [],
                'total_interactions': 0,
                'processing_time_ms': 0,
                'atom_count': len(atoms),
                'warning': 'Molecular engines not available'
            }
        
        try:
            # Initialize molecular engine
            self.molecular_engine.initialize(atoms, bonds or [])
            
            # Run interaction pipeline
            results = self.interaction_pipeline.analyze(atoms, bonds or [])
            
            processing_time = (datetime.now() - start_time).total_seconds() * 1000
            
            return {
                'hydrogen_bonds': results.get('hydrogen_bonds', []),
                'vdw_contacts': results.get('vdw_contacts', []),
                'salt_bridges': results.get('salt_bridges', []),
                'total_interactions': sum(len(v) for v in results.values()),
                'processing_time_ms': processing_time,
                'atom_count': len(atoms),
                'bond_count': len(bonds) if bonds else len(self.molecular_engine.bonds),
            }
            
        except Exception as e:
            logger.error(f"Analysis failed: {e}", exc_info=True)
            return {
                'hydrogen_bonds': [],
                'vdw_contacts': [],
                'salt_bridges': [],
                'total_interactions': 0,
                'processing_time_ms': 0,
                'atom_count': len(atoms),
                'error': str(e)
            }
    
    def get_summary(self, results: Dict[str, Any]) -> str:
        """Get a human-readable summary of the analysis results"""
        hbonds = len(results.get('hydrogen_bonds', []))
        vdw = len(results.get('vdw_contacts', []))
        salt = len(results.get('salt_bridges', []))
        
        return f"Found {hbonds} hydrogen bonds, {vdw} VdW contacts, and {salt} salt bridges"
