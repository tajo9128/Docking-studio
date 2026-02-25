import numpy as np
from typing import List, Dict, Any, Optional, Tuple
import logging

logger = logging.getLogger(__name__)

try:
    from Bio.PDB import PDBParser, Selection, Atom
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    logger.warning("BioPython not available - some features disabled")


def get_atoms_from_pdb(pdb_data: str) -> np.ndarray:
    """Extract atom coordinates from PDB data string"""
    if not BIOPYTHON_AVAILABLE:
        return np.array([])
    
    import tempfile
    import os
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write(pdb_data)
        temp_path = f.name
    
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("mol", temp_path)
        
        atoms = []
        for atom in structure.get_atoms():
            atoms.append(atom.get_coord())
        
        return np.array(atoms)
    finally:
        os.unlink(temp_path)


def get_atoms_from_file(pdb_path: str) -> np.ndarray:
    """Extract atom coordinates from PDB file"""
    if not BIOPYTHON_AVAILABLE:
        return np.array([])
    
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("mol", pdb_path)
        
        atoms = []
        for atom in structure.get_atoms():
            atoms.append(atom.get_coord())
        
        return np.array(atoms)
    except Exception as e:
        logger.error(f"Error parsing PDB file: {e}")
        return np.array([])


def calculate_rmsd(pdb1_data: str, pdb2_data: str) -> float:
    """
    Calculate RMSD between two PDB structures.
    
    Args:
        pdb1_data: First PDB structure (string)
        pdb2_data: Second PDB structure (string)
    
    Returns:
        RMSD value in Angstroms
    """
    atoms1 = get_atoms_from_pdb(pdb1_data)
    atoms2 = get_atoms_from_pdb(pdb2_data)
    
    if len(atoms1) == 0 or len(atoms2) == 0:
        logger.warning("Could not extract atoms from PDB data")
        return -1.0
    
    min_len = min(len(atoms1), len(atoms2))
    atoms1 = atoms1[:min_len]
    atoms2 = atoms2[:min_len]
    
    rmsd = float(np.sqrt(((atoms1 - atoms2) ** 2).sum() / min_len))
    return rmsd


def calculate_rmsd_files(pdb1_path: str, pdb2_path: str) -> float:
    """Calculate RMSD between two PDB files"""
    atoms1 = get_atoms_from_file(pdb1_path)
    atoms2 = get_atoms_from_file(pdb2_path)
    
    if len(atoms1) == 0 or len(atoms2) == 0:
        return -1.0
    
    min_len = min(len(atoms1), len(atoms2))
    atoms1 = atoms1[:min_len]
    atoms2 = atoms2[:min_len]
    
    return float(np.sqrt(((atoms1 - atoms2) ** 2).sum() / min_len))


def detect_hbonds(receptor_atoms: np.ndarray, ligand_atoms: np.ndarray, 
                cutoff: float = 3.5) -> List[Tuple[int, int, float]]:
    """
    Detect hydrogen bonds between receptor and ligand.
    
    Args:
        receptor_atoms: Receptor atom coordinates
        ligand_atoms: Ligand atom coordinates  
        cutoff: Distance cutoff in Angstroms
    
    Returns:
        List of (receptor_idx, ligand_idx, distance) tuples
    """
    hbonds = []
    
    for i, ra in enumerate(receptor_atoms):
        for j, la in enumerate(ligand_atoms):
            dist = np.linalg.norm(ra - la)
            if dist < cutoff:
                hbonds.append((i, j, float(dist)))
    
    return hbonds


def detect_hydrophobic_contacts(receptor_atoms: np.ndarray, ligand_atoms: np.ndarray,
                                cutoff: float = 4.5) -> List[Tuple[int, int, float]]:
    """Detect hydrophobic contacts"""
    contacts = []
    
    for i, ra in enumerate(receptor_atoms):
        for j, la in enumerate(ligand_atoms):
            dist = np.linalg.norm(ra - la)
            if dist < cutoff:
                contacts.append((i, j, float(dist)))
    
    return contacts


def analyze_pose(receptor_pdb: str, ligand_pdb: str, 
                receptor_atoms: Optional[np.ndarray] = None,
                ligand_atoms: Optional[np.ndarray] = None) -> Dict[str, Any]:
    """
    Comprehensive pose analysis.
    
    Args:
        receptor_pdb: Receptor PDB data or path
        ligand_pdb: Ligand PDB data or path
        receptor_atoms: Pre-computed receptor atoms (optional)
        ligand_atoms: Pre-computed ligand atoms (optional)
    
    Returns:
        Dictionary with analysis results
    """
    if receptor_atoms is None:
        if '\n' in receptor_pdb:
            receptor_atoms = get_atoms_from_pdb(receptor_pdb)
        else:
            receptor_atoms = get_atoms_from_file(receptor_pdb)
    
    if ligand_atoms is None:
        if '\n' in ligand_pdb:
            ligand_atoms = get_atoms_from_pdb(ligand_pdb)
        else:
            ligand_atoms = get_atoms_from_file(ligand_pdb)
    
    if len(receptor_atoms) == 0 or len(ligand_atoms) == 0:
        return {
            "hbonds": [],
            "hydrophobic_contacts": [],
            "total_contacts": 0,
            "binding_score": 0.0
        }
    
    hbonds = detect_hbonds(receptor_atoms, ligand_atoms)
    hydrophobic = detect_hydrophobic_contacts(receptor_atoms, ligand_atoms)
    
    avg_hbond_dist = np.mean([d for _, _, d in hbonds]) if hbonds else 0.0
    
    return {
        "hbonds": hbonds,
        "hbond_count": len(hbonds),
        "hydrophobic_contacts": hydrophobic,
        "hydrophobic_count": len(hydrophobic),
        "total_contacts": len(hbonds) + len(hydrophobic),
        "avg_hbond_distance": avg_hbond_dist,
        "binding_score": len(hbonds) * 0.5 + len(hydrophobic) * 0.2
    }


def calculate_advanced_interactions(receptor_pdb: str, ligand_pdb: str) -> Dict[str, Any]:
    """
    Advanced interaction analysis including detailed geometry.
    
    Returns:
        Detailed interaction analysis
    """
    if not BIOPYTHON_AVAILABLE:
        return {
            "error": "BioPython not available",
            "hbonds": [],
            "hydrophobic": [],
            "pi_stacking": [],
            "salt_bridges": []
        }
    
    import tempfile
    import os
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write(receptor_pdb)
        rec_path = f.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write(ligand_pdb)
        lig_path = f.name
    
    try:
        parser = PDBParser(QUIET=True)
        rec = parser.get_structure("rec", rec_path)
        lig = parser.get_structure("lig", lig_path)
        
        rec_atoms = np.array([a.get_coord() for a in rec.get_atoms()])
        lig_atoms = np.array([a.get_coord() for a in lig.get_atoms()])
        
        hbonds = detect_hbonds(rec_atoms, lig_atoms, 3.5)
        hydrophobic = detect_hydrophobic_contacts(rec_atoms, lig_atoms, 4.5)
        
        return {
            "hbonds": [
                {"receptor_idx": r, "ligand_idx": l, "distance": d}
                for r, l, d in hbonds
            ],
            "hydrophobic": [
                {"receptor_idx": r, "ligand_idx": l, "distance": d}
                for r, l, d in hydrophobic
            ],
            "pi_stacking": [],
            "salt_bridges": [],
            "total_hbonds": len(hbonds),
            "total_hydrophobic": len(hydrophobic)
        }
    finally:
        os.unlink(rec_path)
        os.unlink(lig_path)


def get_binding_site_residues(receptor_pdb: str, ligand_pdb: str, 
                              cutoff: float = 5.0) -> List[Dict[str, Any]]:
    """Get binding site residues within cutoff of ligand"""
    if not BIOPYTHON_AVAILABLE:
        return []
    
    import tempfile
    import os
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write(receptor_pdb)
        rec_path = f.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write(ligand_pdb)
        lig_path = f.name
    
    try:
        parser = PDBParser(QUIET=True)
        rec = parser.get_structure("rec", rec_path)
        lig = parser.get_structure("lig", lig_path)
        
        lig_atoms = np.array([a.get_coord() for a in lig.get_atoms()])
        
        binding_residues = []
        
        for residue in rec.get_residues():
            for atom in residue:
                atom_coord = atom.get_coord()
                
                for la in lig_atoms:
                    dist = np.linalg.norm(atom_coord - la)
                    if dist < cutoff:
                        binding_residues.append({
                            "residue_name": residue.resname,
                            "residue_number": residue.id[1],
                            "chain": residue.parent().id,
                            "atom_name": atom.name
                        })
                        break
        
        unique_residues = []
        seen = set()
        for r in binding_residues:
            key = (r['residue_name'], r['residue_number'], r['chain'])
            if key not in seen:
                seen.add(key)
                unique_residues.append(r)
        
        return unique_residues
    finally:
        os.unlink(rec_path)
        os.unlink(lig_path)
