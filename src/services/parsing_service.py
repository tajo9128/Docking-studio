"""Parsing Service - Orchestrates Structure Parsing"""

from typing import Optional, List

from ..core.parsers.pdb_parser import PDBParser
from ..core.parsers.sdf_parser import SDFParser
from ..core.parsers.mol2_parser import MOL2Parser
from ..database import Structure, get_db
from ..schemas import StructureParseResponse, AtomModel, BondModel, StructureMetadata
from ..logging_config import get_logger
from ..core.validators import StructureValidator
from ..core.exceptions import ParseException
from ..core.utils import PerformanceTimer

logger = get_logger(__name__)

class ParsingService:
    """Parsing service for structure files"""
    
    PARSERS = {
        'pdb': PDBParser,
        'pdbqt': PDBParser,
        'sdf': SDFParser,
        'mol2': MOL2Parser,
        'mol': SDFParser,
        'sd': SDFParser,
        'mcif': PDBParser,
        'mmcif': PDBParser,
    }
    
    def __init__(self):
        self.parsers = {}
        for ext, parser_class in self.PARSERS.items():
            self.parsers[ext] = parser_class(strict_mode=True)
    
    async def parse_structure(self, structure_id: str, content: str, filename: str) -> StructureParseResponse:
        """Parse structure file and save to database"""
        logger.info(f"Parsing structure: {filename}")
        
        file_ext = filename.split('.')[-1].lower() if '.' in filename else 'pdb'
        parser = self.parsers.get(file_ext)
        
        if not parser:
            raise ParseException(message=f"Unsupported file type: {file_ext}", code="UNSUPPORTED_FILE_TYPE")
        
        try:
            with PerformanceTimer("Parsing"):
                parse_result = await parser.parse(content)
                
                atoms = []
                bonds = []
                
                for i, atom in enumerate(parse_result.atoms):
                    atoms.append(AtomModel(
                        index=i,
                        serial=atom.get('serial', i),
                        name=atom.get('name', ''),
                        alt_loc=atom.get('alt_loc', ''),
                        res_name=atom.get('res_name', ''),
                        chain_id=atom.get('chain_id', ''),
                        res_seq=atom.get('res_seq', 0),
                        i_code=atom.get('i_code', ''),
                        x=atom['x'],
                        y=atom['y'],
                        z=atom['z'],
                        occupancy=atom.get('occupancy', 1.0),
                        temp_factor=atom.get('temp_factor', 0.0),
                        element=atom.get('element', 'C'),
                        charge=atom.get('charge', 0.0),
                    ))
                
                if parse_result.bonds:
                    for bond in parse_result.bonds:
                        bonds.append(BondModel(
                            atom1_index=bond['atom1_index'],
                            atom2_index=bond['atom2_index'],
                            type=bond['type'],
                            order=bond.get('order', 1),
                            distance=bond['distance'],
                        ))
                
                chains = set()
                for atom in atoms:
                    if atom.chain_id:
                        chains.add(atom.chain_id)
                
                metadata = StructureMetadata(
                    file_name=filename,
                    file_size=len(content),
                    atom_count=len(atoms),
                    bond_count=len(bonds),
                    chain_count=len(chains),
                    model_count=len(parse_result.models) if parse_result.models else 1,
                    title=parse_result.metadata.get('title'),
                    experimental_technique=parse_result.metadata.get('experimental_technique'),
                    resolution=parse_result.metadata.get('resolution'),
                    warnings=parse_result.warnings,
                )
                
                if not StructureValidator.validate_atom_count(metadata.atom_count):
                    raise ParseException(
                        message=f"Structure contains {metadata.atom_count} atoms, exceeds maximum",
                        code="ATOM_COUNT_EXCEEDED"
                    )
                
                async with get_db() as db:
                    structure = await db.get(Structure, structure_id)
                    if not structure:
                        raise ParseException(message="Structure not found", code="STRUCTURE_NOT_FOUND")
                    
                    structure.parsed_data = {
                        'atoms': [atom.dict() for atom in atoms],
                        'bonds': [bond.dict() for bond in bonds],
                        'metadata': metadata.dict(),
                    }
                    structure.atom_count = metadata.atom_count
                    structure.bond_count = metadata.bond_count
                    
                    await db.commit()
                
                logger.info(f"Structure parsed successfully: {filename}")
                
                return StructureParseResponse(
                    structure_id=structure_id,
                    metadata=metadata,
                    atoms=atoms,
                    bonds=bonds if bonds else None,
                    stage="parsed",
                    timestamp="",
                )
                
        except Exception as e:
            logger.error(f"Failed to parse structure: {filename}", exc_info=True)
            raise ParseException(message=f"Failed to parse structure: {str(e)}", code="PARSE_ERROR")
