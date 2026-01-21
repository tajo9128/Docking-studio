"""Parsing Service - Orchestrates Structure Parsing"""

from typing import Optional, List, Dict, Any
import logging

logger = logging.getLogger(__name__)


class ParsingService:
    """Parsing service for structure files"""
    
    SUPPORTED_FORMATS = ['pdb', 'pdbqt', 'sdf', 'mol2', 'mol', 'sd', 'mmcif']
    
    def __init__(self):
        self.parsers = {}
        self._initialize_parsers()
    
    def _initialize_parsers(self):
        """Initialize available parsers"""
        try:
            from ..core.parsers.pdb_parser import PDBParser
            from ..core.parsers.sdf_parser import SDFParser
            from ..core.parsers.mol2_parser import MOL2Parser
            
            self.parsers = {
                'pdb': PDBParser(strict_mode=True),
                'pdbqt': PDBParser(strict_mode=True),
                'sdf': SDFParser(strict_mode=True),
                'mol2': MOL2Parser(strict_mode=True),
                'mol': SDFParser(strict_mode=True),
                'sd': SDFParser(strict_mode=True),
                'mmcif': PDBParser(strict_mode=True),
            }
        except ImportError as e:
            logger.warning(f"Some parsers not available: {e}")
    
    async def parse_structure(self, content: str, filename: str) -> Dict[str, Any]:
        """
        Parse structure file content
        
        Args:
            content: File content as string
            filename: Original filename
            
        Returns:
            Dictionary with parsed atoms, bonds, and metadata
        """
        logger.info(f"Parsing structure: {filename}")
        
        file_ext = filename.split('.')[-1].lower() if '.' in filename else 'pdb'
        parser = self.parsers.get(file_ext)
        
        if not parser:
            return {
                'error': f"Unsupported file type: {file_ext}",
                'atoms': [],
                'bonds': [],
                'metadata': {'filename': filename}
            }
        
        try:
            parse_result = await parser.parse(content)
            
            atoms = []
            for i, atom in enumerate(parse_result.atoms):
                atoms.append({
                    'index': i,
                    'serial': atom.get('serial', i),
                    'name': atom.get('name', ''),
                    'res_name': atom.get('res_name', ''),
                    'chain_id': atom.get('chain_id', ''),
                    'res_seq': atom.get('res_seq', 0),
                    'x': atom.get('x', 0.0),
                    'y': atom.get('y', 0.0),
                    'z': atom.get('z', 0.0),
                    'element': atom.get('element', 'C'),
                })
            
            bonds = []
            if hasattr(parse_result, 'bonds') and parse_result.bonds:
                for bond in parse_result.bonds:
                    bonds.append({
                        'atom1_index': bond.get('atom1_index'),
                        'atom2_index': bond.get('atom2_index'),
                        'type': bond.get('type', 'single'),
                        'order': bond.get('order', 1),
                    })
            
            metadata = {
                'filename': filename,
                'file_size': len(content),
                'atom_count': len(atoms),
                'bond_count': len(bonds),
            }
            
            if hasattr(parse_result, 'metadata'):
                metadata.update(parse_result.metadata)
            
            logger.info(f"Parsed {len(atoms)} atoms from {filename}")
            
            return {
                'atoms': atoms,
                'bonds': bonds,
                'metadata': metadata,
            }
            
        except Exception as e:
            logger.error(f"Failed to parse structure: {filename}", exc_info=True)
            return {
                'error': str(e),
                'atoms': [],
                'bonds': [],
                'metadata': {'filename': filename}
            }
    
    def get_supported_formats(self) -> List[str]:
        """Get list of supported file formats"""
        return self.SUPPORTED_FORMATS
