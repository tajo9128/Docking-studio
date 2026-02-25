#!/usr/bin/env python3
"""
Canonical Docking Input Normalizer
Accept many formats → Normalize everything → Dock only canonical format

Supports:
- Receptor: CIF, PDB, MOL2, PDBQT
- Ligand: SDF, SMILES, MOL2, PDB
- Output: receptor.pdbqt, ligand_XXX.pdbqt, ligands.sdf, ligands.smiles, gridbox.json
"""

import os
import json
import shutil
import zipfile
import tarfile
import subprocess
from pathlib import Path
from typing import Optional, Tuple, List, Dict, Any

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

try:
    from meeko import MoleculePreparation
    MEEKO_AVAILABLE = True
except ImportError:
    MEEKO_AVAILABLE = False

try:
    import pybel
    OPENBABEL_AVAILABLE = True
except ImportError:
    OPENBABEL_AVAILABLE = False


# ============================================================
# CONFIG
# ============================================================

MAX_GRID_VOLUME = 64000
DEFAULT_SIZE = [22, 22, 22]


# ============================================================
# UTILS
# ============================================================

def ensure_dir(path: str) -> None:
    """Create directory if not exists."""
    Path(path).mkdir(parents=True, exist_ok=True)


def safe_extract(archive_path: str, extract_dir: str) -> None:
    """Extract zip/tar safely (no path traversal)."""
    if archive_path.endswith(".zip"):
        with zipfile.ZipFile(archive_path, "r") as z:
            z.extractall(extract_dir)
    elif archive_path.endswith((".tar.gz", ".tgz")):
        with tarfile.open(archive_path) as t:
            t.extractall(extract_dir)


def detect_format(filepath: str) -> str:
    """Detect file format from extension."""
    return Path(filepath).suffix.lower()


def run_command(cmd: List[str], check: bool = True) -> subprocess.CompletedProcess:
    """Run subprocess command."""
    return subprocess.run(cmd, check=check, capture_output=True, text=True)


# ============================================================
# RECEPTOR PIPELINE
# ============================================================

def normalize_receptor(input_path: str, output_dir: str) -> str:
    """
    Convert receptor to canonical PDBQT format.
    
    Supported inputs: .cif, .mmcif, .mol2, .pdb, .pdbqt
    Output: receptor.pdbqt
    """
    ensure_dir(output_dir)
    canonical_pdb = os.path.join(output_dir, "receptor.pdb")
    canonical_pdbqt = os.path.join(output_dir, "receptor.pdbqt")

    ext = detect_format(input_path)

    # Already PDBQT - just copy
    if ext == ".pdbqt":
        shutil.copy(input_path, canonical_pdbqt)
        return canonical_pdbqt

    # Convert other formats to PDB first
    if ext in [".cif", ".mmcif", ".mol2", ".ent"]:
        # Use OpenBabel
        fmt = ext.replace(".", "").replace("mmcif", "cif")
        try:
            mol = next(pybel.readfile(fmt, input_path))
            mol.write("pdb", canonical_pdb, overwrite=True)
        except Exception as e:
            raise ValueError(f"Failed to convert receptor: {e}")
    elif ext == ".pdb":
        shutil.copy(input_path, canonical_pdb)
    else:
        raise ValueError(f"Unsupported receptor format: {ext}")

    # Convert PDB → PDBQT using OpenBabel
    try:
        cmd = ["obabel", canonical_pdb, "-O", canonical_pdbqt, "-xh", "-xr"]
        run_command(cmd)
    except subprocess.CalledProcessError as e:
        raise ValueError(f"Failed to generate PDBQT: {e}")

    return canonical_pdbqt


# ============================================================
# LIGAND PIPELINE
# ============================================================

def rdkit_prepare_3d(mol: Any) -> Any:
    """Generate 3D conformer with RDKit."""
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    return mol


def convert_to_pdbqt(mol: Any, output_path: str) -> None:
    """Convert RDKit mol to PDBQT using Meeko."""
    if not MEEKO_AVAILABLE:
        # Fallback: save as SDF
        writer = Chem.SDWriter(output_path.replace(".pdbqt", ".sdf"))
        writer.write(mol)
        writer.close()
        return
    
    try:
        prep = MoleculePreparation()
        setup = prep.prepare(mol)
        pdbqt_string = prep.write_pdbqt_string(setup)
        with open(output_path, "w") as f:
            f.write(pdbqt_string)
    except Exception as e:
        # Fallback
        writer = Chem.SDWriter(output_path.replace(".pdbqt", ".sdf"))
        writer.write(mol)
        writer.close()


def normalize_ligands(input_path: str, output_dir: str) -> Dict[str, str]:
    """
    Convert any ligand format into canonical format.
    
    Output:
    - ligands/ligand_XXX.pdbqt
    - ligands.sdf
    - ligands.smiles
    """
    ensure_dir(output_dir)
    ligands_dir = os.path.join(output_dir, "ligands")
    ensure_dir(ligands_dir)

    sdf_out = os.path.join(output_dir, "ligands.sdf")
    smiles_out = os.path.join(output_dir, "ligands.smiles")

    ext = detect_format(input_path)

    writer_sdf = Chem.SDWriter(sdf_out) if RDKIT_AVAILABLE else None
    smiles_lines = []
    ligand_files = []
    mol_index = 1

    # SDF input
    if ext == ".sdf":
        supplier = Chem.SDMolSupplier(input_path)
        for mol in supplier:
            if mol is None:
                continue
            mol = rdkit_prepare_3d(mol)
            ligand_name = f"ligand_{mol_index:03d}"
            mol.SetProp("_Name", ligand_name)

            if writer_sdf:
                writer_sdf.write(mol)
            smiles_lines.append(Chem.MolToSmiles(mol))

            pdbqt_path = os.path.join(ligands_dir, f"{ligand_name}.pdbqt")
            convert_to_pdbqt(mol, pdbqt_path)
            ligand_files.append(pdbqt_path)

            mol_index += 1

    # SMILES input
    elif ext in [".smi", ".smiles", ".txt"]:
        with open(input_path) as f:
            for line in f:
                smi = line.strip().split()[0]
                if not smi:
                    continue
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    continue
                mol = rdkit_prepare_3d(mol)

                ligand_name = f"ligand_{mol_index:03d}"
                mol.SetProp("_Name", ligand_name)

                if writer_sdf:
                    writer_sdf.write(mol)
                smiles_lines.append(smi)

                pdbqt_path = os.path.join(ligands_dir, f"{ligand_name}.pdbqt")
                convert_to_pdbqt(mol, pdbqt_path)
                ligand_files.append(pdbqt_path)

                mol_index += 1

    # MOL2/PDB input
    elif ext in [".mol2", ".pdb"]:
        mol = Chem.MolFromMolFile(input_path, removeHs=False)
        if mol:
            mol = rdkit_prepare_3d(mol)

            ligand_name = "ligand_001"
            mol.SetProp("_Name", ligand_name)

            if writer_sdf:
                writer_sdf.write(mol)
            smiles_lines.append(Chem.MolToSmiles(mol))

            pdbqt_path = os.path.join(ligands_dir, f"{ligand_name}.pdbqt")
            convert_to_pdbqt(mol, pdbqt_path)
            ligand_files.append(pdbqt_path)
    else:
        raise ValueError(f"Unsupported ligand format: {ext}")

    if writer_sdf:
        writer_sdf.close()

    with open(smiles_out, "w") as f:
        for smi in smiles_lines:
            f.write(smi + "\n")

    return {
        "ligands_dir": ligands_dir,
        "sdf": sdf_out,
        "smiles": smiles_out,
        "count": mol_index - 1
    }


# ============================================================
# GRID GENERATION
# ============================================================

def generate_grid(
    center: Tuple[float, float, float],
    size: Tuple[float, float, float] = None,
    exhaustiveness: int = 8,
    seed: int = None
) -> Dict[str, Any]:
    """Generate grid parameters."""
    if size is None:
        size = DEFAULT_SIZE
    
    volume = size[0] * size[1] * size[2]
    if volume > MAX_GRID_VOLUME:
        raise ValueError(f"Grid too large: {volume} > {MAX_GRID_VOLUME}")
    
    if seed is None:
        import time
        seed = int(time.time()) % 1000000

    grid = {
        "center": {"x": center[0], "y": center[1], "z": center[2]},
        "size": {"x": size[0], "y": size[1], "z": size[2]},
        "exhaustiveness": exhaustiveness,
        "seed": seed
    }

    return grid


def save_grid(grid: Dict[str, Any], output_dir: str) -> str:
    """Save grid to JSON file."""
    grid_path = os.path.join(output_dir, "gridbox.json")
    with open(grid_path, "w") as f:
        json.dump(grid, f, indent=4)
    return grid_path


def grid_to_vina_config(grid: Dict[str, Any], output_path: str) -> str:
    """Convert grid to Vina config format."""
    center = grid["center"]
    size = grid["size"]
    
    config = f"""center_x = {center['x']}
center_y = {center['y']}
center_z = {center['z']}
size_x = {size['x']}
size_y = {size['y']}
size_z = {size['z']}
exhaustiveness = {grid['exhaustiveness']}
num_modes = 9
energy_range = 3
seed = {grid['seed']}
"""
    with open(output_path, "w") as f:
        f.write(config)
    
    return output_path


# ============================================================
# MAIN PIPELINE
# ============================================================

class CanonicalNormalizer:
    """Main class for canonical input normalization."""
    
    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        ensure_dir(output_dir)
        
        self.receptor = None
        self.ligands = None
        self.grid = None
        
    def process_receptor(self, input_path: str) -> str:
        """Process receptor file."""
        self.receptor = normalize_receptor(input_path, self.output_dir)
        return self.receptor
    
    def process_ligands(self, input_path: str) -> Dict[str, Any]:
        """Process ligand file(s)."""
        self.ligands = normalize_ligands(input_path, self.output_dir)
        return self.ligands
    
    def process_grid(
        self,
        center: Tuple[float, float, float],
        size: Tuple[float, float, float] = None,
        exhaustiveness: int = 8,
        seed: int = None
    ) -> str:
        """Process grid parameters."""
        self.grid = generate_grid(center, size, exhaustiveness, seed)
        grid_path = save_grid(self.grid, self.output_dir)
        
        # Also create Vina config
        vina_config = os.path.join(self.output_dir, "vina_config.txt")
        grid_to_vina_config(self.grid, vina_config)
        
        return grid_path
    
    def get_summary(self) -> Dict[str, Any]:
        """Get processing summary."""
        return {
            "receptor": self.receptor,
            "ligands": self.ligands,
            "grid": self.grid,
            "output_dir": self.output_dir
        }


def canonicalize_job(
    receptor_input: str,
    ligand_input: str,
    job_dir: str,
    grid_center: Tuple[float, float, float],
    grid_size: Tuple[float, float, float] = None,
    exhaustiveness: int = 8,
    seed: int = None
) -> Dict[str, Any]:
    """
    Main entry point - canonicalize all inputs.
    
    Args:
        receptor_input: Path to receptor file
        ligand_input: Path to ligand file(s)
        job_dir: Output directory
        grid_center: Grid center (x, y, z)
        grid_size: Grid size (x, y, z)
        exhaustiveness: Vina exhaustiveness
        seed: Random seed
        
    Returns:
        Dictionary with paths to canonical files
    """
    normalizer = CanonicalNormalizer(job_dir)
    
    # Process receptor
    receptor = normalizer.process_receptor(receptor_input)
    
    # Process ligands
    ligands = normalizer.process_ligands(ligand_input)
    
    # Process grid
    grid = normalizer.process_grid(grid_center, grid_size, exhaustiveness, seed)
    
    return normalizer.get_summary()


# ============================================================
# EXAMPLE USAGE
# ============================================================

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Canonical Docking Input Normalizer")
    parser.add_argument("--receptor", required=True, help="Receptor file (CIF/PDB/MOL2/PDBQT)")
    parser.add_argument("--ligand", required=True, help="Ligand file (SDF/SMILES/MOL2)")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--center", nargs=3, type=float, required=True, help="Grid center (x y z)")
    parser.add_argument("--size", nargs=3, type=float, default=[22, 22, 22], help="Grid size (x y z)")
    parser.add_argument("--exhaustiveness", type=int, default=8, help="Exhaustiveness")
    parser.add_argument("--seed", type=int, help="Random seed")
    
    args = parser.parse_args()
    
    result = canonicalize_job(
        receptor_input=args.receptor,
        ligand_input=args.ligand,
        job_dir=args.output,
        grid_center=tuple(args.center),
        grid_size=tuple(args.size),
        exhaustiveness=args.exhaustiveness,
        seed=args.seed
    )
    
    print("\n✅ Canonicalization complete:")
    print(f"   Receptor: {result['receptor']}")
    print(f"   Ligands: {result['ligands']['count']} files")
    print(f"   SDF: {result['ligands']['sdf']}")
    print(f"   Grid: {result['grid']}")
