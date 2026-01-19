"""Analysis Service - Orchestrates Molecular Analysis"""

from typing import Optional
from datetime import datetime

from ..database import Structure, Interaction, get_db
from ..schemas import AnalysisResponse, AnalysisMetadata, HydrogenBond, VDWContact, SaltBridge
from ..logging_config import get_logger
from ..core.engines.molecular_engine import MolecularEngine
from ..core.engines.interaction_pipeline import InteractionPipeline
from ..core.exceptions import AnalysisException
from ..core.utils import PerformanceTimer, get_current_time_ms

logger = get_logger(__name__)

class AnalysisService:
    """Analysis service for molecular interactions"""
    
    def __init__(self):
        self.molecular_engine = MolecularEngine()
        self.interaction_pipeline = InteractionPipeline()
    
    async def analyze_interactions(self, structure_id: str, options: Optional[dict] = None) -> AnalysisResponse:
        """Analyze molecular interactions (H-bonds, VdW, Salt Bridges)"""
        logger.info(f"Analyzing interactions: {structure_id}")
        
        start_time = get_current_time_ms()
        
        async with get_db() as db:
            structure = await db.get(Structure, structure_id)
            
            if not structure or not structure.parsed_data:
                raise AnalysisException(
                    message="Structure not found or not parsed",
                    code="STRUCTURE_NOT_FOUND"
                )
            
            atoms_data = structure.parsed_data.get('atoms', [])
            bonds_data = structure.parsed_data.get('bonds', [])
            
            if not atoms_data:
                raise AnalysisException(message="No atoms found in structure", code="NO_ATOMS")
            
            if len(atoms_data) < 2:
                return AnalysisResponse(
                    structure_id=structure_id,
                    hydrogen_bonds=[],
                    vdw_contacts=[],
                    salt_bridges=[],
                    total_interactions=0,
                    metadata=AnalysisMetadata(
                        processing_time_ms=0,
                        atom_count=len(atoms_data),
                        bond_count=0,
                        algorithm="skipped",
                        thresholds={},
                    ),
                    stage="analyzed",
                    timestamp=datetime.now().isoformat(),
                )
            
            try:
                with PerformanceTimer("Interaction Analysis"):
                    self.molecular_engine.initialize(atoms_data, bonds_data)
                    interaction_results = self.interaction_pipeline.analyze(atoms_data, bonds_data)
                    
                    hydrogen_bonds = []
                    vdw_contacts = []
                    salt_bridges = []
                    
                    for hb in interaction_results.get('hydrogen_bonds', []):
                        hydrogen_bonds.append(HydrogenBond(
                            atom1_index=hb['atom1_index'],
                            atom2_index=hb['atom2_index'],
                            distance=hb['distance'],
                            angle=hb.get('angle'),
                            atom1_residue=hb['atom1_residue'],
                            atom1_residue_seq=hb['atom1_residue_seq'],
                            atom2_residue=hb['atom2_residue'],
                            atom2_residue_seq=hb['atom2_residue_seq'],
                            confidence=hb['confidence'],
                            is_predicted=False,
                        ))
                    
                    for vdw in interaction_results.get('vdw_contacts', []):
                        vdw_contacts.append(VDWContact(
                            atom1_index=vdw['atom1_index'],
                            atom2_index=vdw['atom2_index'],
                            distance=vdw['distance'],
                            atom1_residue=vdw['atom1_residue'],
                            atom1_residue_seq=vdw['atom1_residue_seq'],
                            atom2_residue=vdw['atom2_residue'],
                            atom2_residue_seq=vdw['atom2_residue_seq'],
                            confidence=vdw['confidence'],
                            is_predicted=False,
                        ))
                    
                    for sb in interaction_results.get('salt_bridges', []):
                        salt_bridges.append(SaltBridge(
                            atom1_index=sb['atom1_index'],
                            atom2_index=sb['atom2_index'],
                            distance=sb['distance'],
                            atom1_residue=sb['atom1_residue'],
                            atom1_residue_seq=sb['atom1_residue_seq'],
                            atom2_residue=sb['atom2_residue'],
                            atom2_residue_seq=sb['atom2_residue_seq'],
                            confidence=sb['confidence'],
                            is_predicted=False,
                        ))
                    
                    # Save to database
                    for hb in hydrogen_bonds:
                        db.add(Interaction(
                            structure_id=structure_id,
                            interaction_type="hydrogen_bond",
                            atom1_index=hb.atom1_index,
                            atom2_index=hb.atom2_index,
                            distance=hb.distance,
                            atom1_residue=hb.atom1_residue,
                            atom1_residue_seq=hb.atom1_residue_seq,
                            atom2_residue=hb.atom2_residue,
                            atom2_residue_seq=hb.atom2_residue_seq,
                            confidence=hb.confidence,
                            metadata={'angle': hb.angle},
                        ))
                    
                    total_interactions = len(hydrogen_bonds) + len(vdw_contacts) + len(salt_bridges)
                    
                    structure.analysis_data = {
                        'hydrogen_bonds': len(hydrogen_bonds),
                        'vdw_contacts': len(vdw_contacts),
                        'salt_bridges': len(salt_bridges),
                        'total_interactions': total_interactions,
                    }
                    
                    await db.commit()
                    
                    processing_time = get_current_time_ms() - start_time
                    
                    logger.info(f"Analysis complete: {structure_id}")
                    
                    return AnalysisResponse(
                        structure_id=structure_id,
                        hydrogen_bonds=hydrogen_bonds,
                        vdw_contacts=vdw_contacts,
                        salt_bridges=salt_bridges,
                        total_interactions=total_interactions,
                        metadata=AnalysisMetadata(
                            processing_time_ms=processing_time,
                            atom_count=len(atoms_data),
                            bond_count=len(bonds_data) if bonds_data else 0,
                            algorithm="O(n) spatial hash grid",
                            thresholds=self.interaction_pipeline.thresholds.dict(),
                        ),
                        stage="analyzed",
                        timestamp=datetime.now().isoformat(),
                    )
                    
            except Exception as e:
                logger.error(f"Analysis failed: {structure_id}", exc_info=True)
                raise AnalysisException(message=f"Failed to analyze: {str(e)}", code="ANALYSIS_ERROR")
