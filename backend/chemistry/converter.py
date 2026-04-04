"""
Universal converter: Any input format → RDKit Mol → PDBQT
Uses: RDKit, Open Babel, Meeko (in order of preference)
"""
import tempfile
import os
import logging
from pathlib import Path
from typing import Union, Optional
from rdkit import Chem
from rdkit.Chem import AllChem

from chemistry.formats import InputFormat, MoleculeType, detect_format, is_valid_smiles, is_valid_inchi

logger = logging.getLogger(__name__)


class ConversionError(Exception):
    """Raised when format conversion fails."""
    pass


class UniversalConverter:
    """Convert any supported input format to RDKit Mol, then to PDBQT."""

    def __init__(self, pH: float = 7.4, add_hydrogens: bool = True):
        self.pH = pH
        self.add_hydrogens = add_hydrogens

    def to_mol(self,
               content: Union[str, Path],
               fmt: Optional[Union[str, InputFormat]] = None,
               mol_type: MoleculeType = MoleculeType.UNKNOWN,
               filename: Optional[str] = None) -> Chem.Mol:
        """Convert any supported format to RDKit Mol object."""
        if fmt is None:
            fmt = detect_format(content, filename, mol_type)
        elif isinstance(fmt, str):
            fmt = InputFormat(fmt.lower())

        logger.debug(f"Converting {fmt.value} → RDKit Mol (type: {mol_type.value})")

        try:
            if fmt == InputFormat.SMILES:
                return self._from_smiles(str(content).strip())
            elif fmt == InputFormat.INCHI:
                return self._from_inchi(str(content).strip())
            elif fmt == InputFormat.PDB:
                return self._from_pdb(content, mol_type)
            elif fmt == InputFormat.PDBQT:
                return self._from_pdbqt(content)
            elif fmt == InputFormat.MOL:
                return self._from_mol(content)
            elif fmt == InputFormat.MOL2:
                return self._from_mol2(content)
            elif fmt == InputFormat.SDF:
                return self._from_sdf(content)
            elif fmt == InputFormat.PDB_ID and mol_type == MoleculeType.PROTEIN:
                return self._from_pdb_id(str(content).strip().upper())
            else:
                raise ConversionError(f"Unsupported format: {fmt.value}")
        except Exception as e:
            raise ConversionError(f"Failed to convert {fmt.value} to Mol: {e}")

    def _from_smiles(self, smiles: str) -> Chem.Mol:
        if not is_valid_smiles(smiles):
            raise ConversionError(f"Invalid SMILES: {smiles[:50]}...")
        mol = Chem.MolFromSmiles(smiles)
        if self.add_hydrogens:
            mol = Chem.AddHs(mol)
        return mol

    def _from_inchi(self, inchi: str) -> Chem.Mol:
        if not is_valid_inchi(inchi):
            raise ConversionError(f"Invalid InChI: {inchi[:50]}...")
        mol = Chem.MolFromInchi(inchi)
        if mol is None:
            raise ConversionError("Failed to parse InChI")
        if self.add_hydrogens:
            mol = Chem.AddHs(mol)
        return mol

    def _from_pdb(self, content: Union[str, Path], mol_type: MoleculeType) -> Chem.Mol:
        if isinstance(content, Path):
            with open(content, 'r') as f:
                content = f.read()
        mol = Chem.MolFromPDBBlock(content, removeHs=False)
        if mol is None:
            raise ConversionError("Failed to parse PDB with RDKit")
        if mol_type == MoleculeType.LIGAND:
            ligand_block = self._extract_ligand_from_pdb(content)
            if ligand_block:
                mol = Chem.MolFromPDBBlock(ligand_block, removeHs=False)
                if mol is None:
                    raise ConversionError("Failed to parse ligand from PDB")
        if self.add_hydrogens:
            mol = Chem.AddHs(mol)
        return mol

    def _extract_ligand_from_pdb(self, pdb_content: str) -> Optional[str]:
        lines = pdb_content.split('\n')
        ligand_lines = [l for l in lines if l.startswith('HETATM')]
        if not ligand_lines:
            return None
        return '\n'.join(ligand_lines) + '\nTER\nEND\n'

    def _from_pdbqt(self, content: Union[str, Path]) -> Chem.Mol:
        if isinstance(content, Path):
            with open(content, 'r') as f:
                content = f.read()
        pdb_lines = []
        for line in content.split('\n'):
            if line.startswith(('ATOM', 'HETATM')) and len(line) >= 76:
                pdb_lines.append(line[:76] + '  \n')
            else:
                pdb_lines.append(line)
        pdb_content = '\n'.join(pdb_lines)
        return self._from_pdb(pdb_content, MoleculeType.UNKNOWN)

    def _from_mol(self, content: Union[str, Path]) -> Chem.Mol:
        if isinstance(content, Path):
            with open(content, 'r') as f:
                content = f.read()
        mol = Chem.MolFromMolBlock(content, removeHs=False)
        if mol is None:
            raise ConversionError("Failed to parse MOL file")
        if self.add_hydrogens:
            mol = Chem.AddHs(mol)
        return mol

    def _from_mol2(self, content: Union[str, Path]) -> Chem.Mol:
        if isinstance(content, Path):
            with open(content, 'r') as f:
                content = f.read()
        mol = Chem.MolFromMol2Block(content, removeHs=False)
        if mol is not None:
            if self.add_hydrogens:
                mol = Chem.AddHs(mol)
            return mol
        return self._from_file_via_obabel(content, 'mol2')

    def _from_sdf(self, content: Union[str, Path]) -> Chem.Mol:
        if isinstance(content, Path):
            with open(content, 'r') as f:
                content = f.read()
        supplier = Chem.SDMolSupplier()
        supplier.SetData(content.encode(), len(content), True, False)
        mol = None
        for m in supplier:
            if m is not None:
                mol = m
                break
        if mol is None:
            raise ConversionError("Failed to parse SDF file")
        if self.add_hydrogens:
            mol = Chem.AddHs(mol)
        return mol

    def _from_pdb_id(self, pdb_id: str) -> Chem.Mol:
        import requests
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(url, timeout=30)
        if response.status_code != 200:
            raise ConversionError(f"Failed to fetch PDB {pdb_id}: HTTP {response.status_code}")
        return self._from_pdb(response.text, MoleculeType.PROTEIN)

    def _from_file_via_obabel(self, content: str, fmt: str) -> Chem.Mol:
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix=f'.{fmt}', delete=False) as f:
            f.write(content)
            temp_path = f.name
        try:
            try:
                from openbabel import pybel
                obmol = next(pybel.readfile(fmt, temp_path))
                smiles = obmol.write('smiles').strip()
                return self._from_smiles(smiles)
            except ImportError:
                import subprocess
                result = subprocess.run(
                    ['obabel', temp_path, '-o', 'smi', '-O', '/dev/stdout'],
                    capture_output=True, text=True, timeout=30
                )
                if result.returncode == 0:
                    smiles = result.stdout.strip().split()[0]
                    return self._from_smiles(smiles)
                raise ConversionError(f"Open Babel failed: {result.stderr}")
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)
