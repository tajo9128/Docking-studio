"""
Format detection and conversion utilities for protein/ligand inputs.
Supports: PDB, PDBQT, MOL, MOL2, SDF, SMILES, InChI
"""
import re
import logging
from enum import Enum
from pathlib import Path
from typing import Union, Optional
from rdkit import Chem

logger = logging.getLogger(__name__)


class MoleculeType(Enum):
    PROTEIN = "protein"
    LIGAND = "ligand"
    UNKNOWN = "unknown"


class InputFormat(Enum):
    SMILES = "smiles"
    INCHI = "inchi"
    PDB = "pdb"
    PDBQT = "pdbqt"
    MOL = "mol"
    MOL2 = "mol2"
    SDF = "sdf"
    PDB_ID = "pdb_id"
    UNKNOWN = "unknown"


def detect_format(content: Union[str, Path],
                  filename: Optional[str] = None,
                  mol_type: MoleculeType = MoleculeType.UNKNOWN) -> InputFormat:
    """Detect input format from content string or file path."""
    if filename:
        ext = Path(filename).suffix.lower().lstrip('.')
        ext_map = {
            'smi': InputFormat.SMILES, 'smiles': InputFormat.SMILES,
            'inchi': InputFormat.INCHI,
            'pdb': InputFormat.PDB, 'pdbqt': InputFormat.PDBQT,
            'mol': InputFormat.MOL, 'sdf': InputFormat.SDF,
            'mol2': InputFormat.MOL2,
        }
        if ext in ext_map:
            return ext_map[ext]

    if isinstance(content, Path):
        try:
            with open(content, 'r') as f:
                content = f.read(2000)
        except Exception as e:
            logger.error(f"Failed to read file {content}: {e}")
            return InputFormat.UNKNOWN

    content_stripped = content.strip()

    if mol_type == MoleculeType.LIGAND:
        if len(content_stripped) < 500 and '\n' not in content_stripped:
            if re.match(r'^[A-Za-z0-9@+\-\[\]\(\)\\%=#]+$', content_stripped):
                if any(p in content_stripped for p in ['=', '#', '(', ')', '@', '%']):
                    return InputFormat.SMILES
                if re.match(r'^[A-Za-z]+$', content_stripped) and len(content_stripped) <= 50:
                    return InputFormat.SMILES

    if content_stripped.startswith('InChI='):
        return InputFormat.INCHI

    if any(line.startswith(('ATOM  ', 'HETATM', 'TER', 'END')) for line in content_stripped.split('\n')[:20]):
        first_atom = next((l for l in content_stripped.split('\n') if l.startswith(('ATOM', 'HETATM'))), None)
        if first_atom and len(first_atom) >= 78:
            atom_type_field = first_atom[76:78].strip()
            valid_ad4_types = {'C', 'A', 'N', 'NA', 'OA', 'S', 'SA', 'P', 'F', 'CL', 'BR', 'I', 'HD'}
            if atom_type_field in valid_ad4_types:
                return InputFormat.PDBQT
        return InputFormat.PDB

    if 'V2000' in content_stripped or 'V3000' in content_stripped:
        if '$$$$' in content_stripped:
            return InputFormat.SDF
        return InputFormat.MOL

    if '@<TRIPOS>' in content_stripped:
        return InputFormat.MOL2

    if mol_type == MoleculeType.PROTEIN:
        if re.match(r'^[0-9][A-Za-z0-9]{3}$', content_stripped.upper()):
            return InputFormat.PDB_ID

    return InputFormat.UNKNOWN


def is_valid_smiles(smiles: str) -> bool:
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except Exception:
        return False


def is_valid_inchi(inchi: str) -> bool:
    try:
        mol = Chem.MolFromInchi(inchi)
        return mol is not None
    except Exception:
        return False
