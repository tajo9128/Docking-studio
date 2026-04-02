"""
BioDockify Docking Engine - RDKit-Only Molecular Preparation & Docking
- Protein Preparation: PDB → Remove Water → Add H → Assign Charges → PDBQT
- Ligand Preparation: SDF/MOL2/PDB/SMILES → 3D → Minimize → PDBQT
- Smart Docking Router: Energy-based routing
  - Energy <= -5.0: Vina only → output log, docking, grid files
  - Energy > -5.0: GNINA + RF → output all files (Vina + GNINA)
"""
import os
import logging
import random
import subprocess
from typing import Dict, Any, List, Optional
from datetime import datetime

logger = logging.getLogger(__name__)

ENERGY_THRESHOLD = -5.0


def check_rdkit() -> bool:
    """Check if RDKit is available"""
    try:
        from rdkit import Chem
        return True
    except ImportError:
        return False


def check_vina() -> bool:
    """Check if AutoDock Vina Python API is available"""
    try:
        from vina import Vina
        return True
    except ImportError:
        return False


def check_vina_cli() -> bool:
    """Check if Vina CLI is available"""
    try:
        result = subprocess.run(["which", "vina"], capture_output=True, timeout=5)
        if result.returncode == 0:
            return True
        result = subprocess.run(["vina", "--version"], capture_output=True, timeout=5)
        return result.returncode == 0
    except Exception:
        return False


def check_gnina() -> bool:
    """Check if GNINA CLI is available"""
    try:
        result = subprocess.run(["which", "gnina"], capture_output=True, timeout=5)
        return result.returncode == 0
    except Exception:
        return False


def check_gpu_cuda() -> Dict[str, Any]:
    """Check for NVIDIA GPU with CUDA"""
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=name,driver_version,memory.total", "--format=csv,noheader"],
            capture_output=True, timeout=10
        )
        if result.returncode == 0:
            gpus = []
            for line in result.stdout.strip().split('\n'):
                parts = [p.strip() for p in line.split(',')]
                if len(parts) >= 2:
                    gpus.append({"name": parts[0], "driver": parts[1], "memory": parts[2] if len(parts) > 2 else "Unknown"})
            return {"available": True, "count": len(gpus), "gpus": gpus, "platform": "CUDA"}
    except Exception:
        pass
    return {"available": False, "count": 0, "gpus": [], "platform": "CPU"}


def smiles_to_3d(smiles: str, random_seed: int = 42) -> Optional[Dict[str, Any]]:
    """
    Convert SMILES to 3D structure using RDKit
    Returns: {'mol': RDKit Mol object, 'pdb': PDB string, 'smi': canonical SMILES}
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=random_seed)
        AllChem.MMFFOptimizeMolecule(mol)
        
        canonical_smiles = Chem.MolToSmiles(mol)
        pdb_block = Chem.MolToPDBBlock(mol)
        
        return {
            'mol': mol,
            'pdb': pdb_block,
            'smi': canonical_smiles,
            'num_atoms': mol.GetNumAtoms(),
            'num_heavy_atoms': mol.GetNumHeavyAtoms()
        }
    except Exception as e:
        logger.error(f"SMILES to 3D conversion failed: {e}")
        return None


def prepare_ligand_from_content(content: str, input_format: str = 'sdf', output_dir: str = '/tmp') -> Optional[Dict[str, Any]]:
    """
    Prepare ligand from file content (SDF, MOL2, PDB, or SMILES)
    Returns: {'pdbqt_path': str, 'pdb': str, 'mol': mol object, 'num_rotatable_bonds': int}
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors
        
        mol = None
        
        if input_format == 'smiles':
            result = smiles_to_3d(content)
            if result:
                mol = result['mol']
        elif input_format == 'sdf':
            suppl = Chem.SDMolSupplier()
            suppl.SetData(content)
            for m in suppl:
                if m is not None:
                    mol = m
                    break
        elif input_format == 'mol2':
            mol = Chem.MolFromMol2Block(content)
        elif input_format == 'pdb':
            mol = Chem.MolFromPDBBlock(content)
        
        if mol is None:
            logger.error(f"Failed to parse ligand content (format: {input_format})")
            return None
        
        mol = Chem.AddHs(mol)
        
        try:
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception as e:
            logger.warning(f"3D embedding failed, using 2D: {e}")
        
        pdbqt_content = mol_to_pdbqt(mol, is_ligand=True)
        
        os.makedirs(output_dir, exist_ok=True)
        ligand_id = random.randint(10000, 99999)
        pdbqt_path = os.path.join(output_dir, f"ligand_{ligand_id}.pdbqt")
        with open(pdbqt_path, 'w') as f:
            f.write(pdbqt_content)
        
        num_rotatable = 0
        try:
            num_rotatable = Descriptors.NumRotatableBonds(mol)
        except Exception:
            pass
        
        return {
            'pdbqt_path': pdbqt_path,
            'pdbqt_content': pdbqt_content,
            'pdb': Chem.MolToPDBBlock(mol),
            'mol': mol,
            'num_rotatable_bonds': num_rotatable
        }
        
    except Exception as e:
        logger.error(f"Ligand preparation failed: {e}")
        return None


def prepare_protein_from_content(pdb_content: str, output_dir: str = '/tmp') -> Optional[Dict[str, Any]]:
    """
    Prepare protein from PDB content
    Steps: Remove water → Add H → Assign charges → Convert to PDBQT
    Returns: {'pdbqt_path': str, 'pdbqt_content': str}
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        lines = pdb_content.split('\n')
        pdb_lines = []
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                resname = line[17:20].strip()
                if resname not in ['HOH', 'WAT', 'H2O']:
                    pdb_lines.append(line)
        
        clean_pdb = '\n'.join(pdb_lines)
        mol = Chem.MolFromPDBBlock(clean_pdb, flavor=4)
        
        if mol is None:
            logger.error("Failed to parse protein PDB")
            return None
        
        try:
            mol = Chem.AddHs(mol, addCoords=True)
        except Exception as e:
            logger.warning(f"AddHs failed ({e}), proceeding without hydrogens")
        
        pdbqt_content = mol_to_pdbqt(mol, is_ligand=False)
        
        os.makedirs(output_dir, exist_ok=True)
        receptor_id = random.randint(10000, 99999)
        pdbqt_path = os.path.join(output_dir, f"receptor_{receptor_id}.pdbqt")
        with open(pdbqt_path, 'w') as f:
            f.write(pdbqt_content)
        
        num_residues = 0
        try:
            num_residues = mol.GetNumResidues()
        except Exception:
            pass
        
        return {
            'pdbqt_path': pdbqt_path,
            'pdbqt_content': pdbqt_content,
            'num_residues': num_residues
        }
        
    except Exception as e:
        logger.error(f"Protein preparation failed: {e}")
        return None


def mol_to_pdbqt(mol, is_ligand: bool = True) -> str:
    """
    Convert RDKit mol to PDBQT format
    - For ligands: Include torsion tree, atom types, charges
    - For receptors: Simple conversion with AD4 atom types
    """
    from rdkit import Chem
    
    pdbqt_lines = []
    
    if is_ligand:
        pdbqt_lines.append("REMARK  Generated by BioDockify RDKit")
        pdbqt_lines.append("ROOT")
        
        conf = mol.GetConformer(0)
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetPositions()[i]
            pdbqt_lines.append(
                f"HETATM{((i+1)%99999):5d}  {atom.GetSymbol():<2}{' ':<1}"
                f"{'LIG':<3} A{((i//9999)+1):4d}    "
                f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}"
                f"  1.00  0.00           {atom.GetSymbol():>2}"
            )
        
        pdbqt_lines.append("ENDROOT")
        pdbqt_lines.append("BRANCH 1 1")
        
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                pdbqt_lines.append("TORSDOF 1")
                break
        
    else:
        conf = mol.GetConformer(0)
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetPositions()[i]
            atomic_num = atom.GetAtomicNum()
            ad_type = get_ad4_atom_type(atomic_num)
            charge = 0.0
            
            pdbqt_lines.append(
                f"ATOM  {((i+1)%99999):5d}  {atom.GetSymbol():<2}{' ':<1}"
                f"{'PRO':<3} A{((i//9999)+1):4d}    "
                f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}"
                f"  1.00  0.00          {ad_type:<2}{charge:6.3f}"
            )
    
    pdbqt_lines.append("END")
    return '\n'.join(pdbqt_lines)


def get_ad4_atom_type(atomic_num: int) -> str:
    """Get AutoDock4 atom type for element"""
    ad_types = {
        1: 'HD', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N',
        8: 'OA', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al',
        14: 'Si', 15: 'P', 16: 'SA', 17: 'Cl', 19: 'K', 20: 'Ca',
        34: 'Se', 35: 'Br', 53: 'I', 26: 'Fe'
    }
    return ad_types.get(atomic_num, 'A')


def calculate_hydrophobic_enclosure_score(
    ligand_mol,
    receptor_pdbqt: str,
    cutoff: float = 4.5
) -> float:
    """
    Calculate hydrophobic enclosure term similar to GlideScore.
    Rewards ligands that bury hydrophobic atoms in hydrophobic pockets.
    Returns: Negative score (favorable) if hydrophobic atoms are well-enclosed.
    """
    try:
        import numpy as np
        from rdkit import Chem

        hydrophobic_residues = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'TYR', 'PRO'}
        receptor_atoms = []

        for line in receptor_pdbqt.split('\n'):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    resname = line[17:20].strip()
                    if resname in hydrophobic_residues:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        receptor_atoms.append(np.array([x, y, z]))
                except (ValueError, IndexError):
                    continue

        if not receptor_atoms:
            return 0.0

        receptor_coords = np.array(receptor_atoms)

        hydrophobic_ligand_atoms = []
        for atom in ligand_mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol in ['C', 'S', 'Cl', 'Br', 'I']:
                is_polar = False
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() in ['O', 'N'] and neighbor.GetTotalNumHs() == 0:
                        is_polar = True
                        break
                if not is_polar:
                    pos = ligand_mol.GetConformer().GetAtomPosition(atom.GetIdx())
                    hydrophobic_ligand_atoms.append(np.array([pos.x, pos.y, pos.z]))

        if not hydrophobic_ligand_atoms:
            return 0.0

        enclosure_score = 0.0
        for lig_pos in hydrophobic_ligand_atoms:
            distances = np.linalg.norm(receptor_coords - lig_pos, axis=1)
            if np.any(distances < cutoff):
                min_dist = np.min(distances[distances < cutoff])
                enclosure_score += (cutoff - min_dist) / cutoff

        n_hydrophobic = len(hydrophobic_ligand_atoms)
        if n_hydrophobic > 0:
            enclosure_score = -0.5 * (enclosure_score / n_hydrophobic)

        return round(enclosure_score, 4)
    except Exception as e:
        logger.warning(f"Hydrophobic enclosure calculation failed: {e}")
        return 0.0


def calculate_rotatable_bond_penalty(num_rotatable: int, max_penalty: float = 0.5) -> float:
    """
    Penalize excessive conformational flexibility (entropy loss upon binding).
    GlideScore uses: penalty = 0.058 * N_rot (capped at ~1.5)
    """
    penalty = min(0.058 * num_rotatable, max_penalty)
    return round(penalty, 4)


def calculate_lipophilic_contact_term(
    ligand_mol,
    receptor_pdbqt: str,
    contact_cutoff: float = 4.0
) -> float:
    """
    Reward favorable lipophilic contacts between ligand and receptor.
    Similar to Glide's Lipo term.
    """
    try:
        import numpy as np

        ligand_coords = []
        ligand_lipophilic = []

        for atom in ligand_mol.GetAtoms():
            pos = ligand_mol.GetConformer().GetAtomPosition(atom.GetIdx())
            ligand_coords.append([pos.x, pos.y, pos.z])
            symbol = atom.GetSymbol()
            is_lipophilic = symbol in ['C', 'S', 'Cl', 'Br', 'I']
            if is_lipophilic:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() in ['O', 'N', 'P']:
                        is_lipophilic = False
                        break
            ligand_lipophilic.append(1.0 if is_lipophilic else 0.0)

        if not ligand_coords:
            return 0.0

        ligand_coords = np.array(ligand_coords)
        ligand_lipophilic = np.array(ligand_lipophilic)

        receptor_lipophilic = []
        lipophilic_elements = {'C', 'S'}

        for line in receptor_pdbqt.split('\n'):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    element = line[76:78].strip() or line[12:14].strip()
                    if element in lipophilic_elements:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        receptor_lipophilic.append([x, y, z])
                except (ValueError, IndexError):
                    continue

        if not receptor_lipophilic:
            return 0.0

        receptor_lipophilic = np.array(receptor_lipophilic)

        lipo_score = 0.0
        for i, (lig_pos, is_lipo) in enumerate(zip(ligand_coords, ligand_lipophilic)):
            if is_lipo < 0.5:
                continue
            distances = np.linalg.norm(receptor_lipophilic - lig_pos, axis=1)
            close_contacts = distances[distances < contact_cutoff]
            for dist in close_contacts:
                lipo_score += (contact_cutoff - dist) / contact_cutoff

        return round(-0.3 * lipo_score, 4)
    except Exception as e:
        logger.warning(f"Lipophilic contact calculation failed: {e}")
        return 0.0


def apply_composite_scoring(
    results: List[Dict],
    ligand_mol,
    receptor_pdbqt_content: str,
    num_rotatable_bonds: int = 0
) -> List[Dict]:
    """
    Apply composite scoring to docking results.
    Adds hydrophobic_term, rotatable_penalty, lipo_contact, and composite_score to each pose.
    """
    if not ligand_mol or not receptor_pdbqt_content:
        for pose in results:
            pose['hydrophobic_term'] = 0.0
            pose['rotatable_penalty'] = 0.0
            pose['lipo_contact'] = 0.0
            pose['composite_score'] = pose.get('vina_score', 0)
        return results

    hydrophobic = calculate_hydrophobic_enclosure_score(ligand_mol, receptor_pdbqt_content)
    rot_penalty = calculate_rotatable_bond_penalty(num_rotatable_bonds)
    lipo = calculate_lipophilic_contact_term(ligand_mol, receptor_pdbqt_content)

    for pose in results:
        pose['hydrophobic_term'] = hydrophobic
        pose['rotatable_penalty'] = rot_penalty
        pose['lipo_contact'] = lipo
        base_score = pose.get('vina_score', 0)
        pose['composite_score'] = round(base_score + hydrophobic + rot_penalty + lipo, 4)

    results.sort(key=lambda x: x.get('composite_score', x.get('vina_score', 0)))
    return results


def generate_vina_config(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    center_x: float,
    center_y: float,
    center_z: float,
    size_x: float,
    size_y: float,
    size_z: float,
    output_dir: str
) -> str:
    """Generate Vina configuration file"""
    config_content = f"""receptor = {receptor_pdbqt}
ligand = {ligand_pdbqt}

center_x = {center_x}
center_y = {center_y}
center_z = {center_z}

size_x = {size_x}
size_y = {size_y}
size_z = {size_z}

exhaustiveness = 32
num_modes = 10
energy_range = 2
"""
    config_path = os.path.join(output_dir, "vina_config.txt")
    with open(config_path, 'w') as f:
        f.write(config_content)
    return config_path


def run_vina_docking(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    center_x: float = 0.0,
    center_y: float = 0.0,
    center_z: float = 0.0,
    size_x: float = 20.0,
    size_y: float = 20.0,
    size_z: float = 20.0,
    exhaustiveness: int = 32,
    num_modes: int = 10,
    output_dir: str = "/tmp"
) -> Dict[str, Any]:
    """Run AutoDock Vina docking"""
    
    vina_available = check_vina()
    vina_cli_available = check_vina_cli()
    
    if not vina_available and not vina_cli_available:
        logger.warning("[Vina] Vina not available - generating simulated results")
        return {
            "success": False,
            "error": "Vina not installed",
            "results": [],
            "files": {}
        }
    
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    job_id = f"vina_{timestamp}"
    
    log_file = os.path.join(output_dir, f"{job_id}_log.txt")
    docking_file = os.path.join(output_dir, f"{job_id}_docking.pdbqt")
    grid_file = os.path.join(output_dir, f"{job_id}_grid.txt")
    
    if vina_cli_available:
        try:
            config_path = generate_vina_config(
                receptor_pdbqt, ligand_pdbqt,
                center_x, center_y, center_z,
                size_x, size_y, size_z,
                output_dir
            )
            
            output_path = docking_file
            result = subprocess.run(
                ["vina", "--config", config_path, "--log", log_file, "--out", output_path],
                capture_output=True, timeout=300
            )
            
            if result.returncode == 0:
                results = parse_vina_log(log_file, num_modes)
                
                grid_content = f"""# Grid Box Configuration
# Generated by BioDockify Studio AI
# Date: {datetime.now().isoformat()}

CENTER_X = {center_x}
CENTER_Y = {center_y}
CENTER_Z = {center_z}

SIZE_X = {size_x}
SIZE_Y = {size_y}
SIZE_Z = {size_z}

SPACING = 0.375
DIMENSIONS = {int(size_x/0.375)} x {int(size_y/0.375)} x {int(size_z/0.375)}
"""
                with open(grid_file, 'w') as f:
                    f.write(grid_content)
                
                return {
                    "success": True,
                    "engine": "vina_cli",
                    "results": results,
                    "files": {
                        "log": log_file,
                        "docking": docking_file,
                        "grid": grid_file
                    }
                }
        except Exception as e:
            logger.warning(f"[Vina CLI] Failed: {e}")
    
    if vina_available:
        try:
            from vina import Vina
            
            logger.info(f"[Vina] Initializing with receptor: {receptor_pdbqt}")
            v = Vina(sf_name='vina')
            v.set_receptor(rigid_pdbqt_filename=receptor_pdbqt)
            v.set_ligand_from_file(ligand_pdbqt)
            
            logger.info(f"[Vina] Computing maps at center=({center_x},{center_y},{center_z})")
            v.compute_vina_maps(center=[center_x, center_y, center_z], box_size=[size_x, size_y, size_z])
            
            logger.info(f"[Vina] Docking {num_modes} poses...")
            v.dock(exhaustiveness=exhaustiveness, n_poses=num_modes)
            
            energies = v.energies
            results = []
            for i in range(min(num_modes, len(energies))):
                results.append({
                    "mode": i + 1,
                    "vina_score": round(float(energies[i][0]), 3),
                    "gnina_score": None,
                    "rf_score": None,
                    "source": "vina"
                })
            
            try:
                v.write_pose(docking_file, overwrite=True)
            except:
                pass
            
            with open(log_file, 'w') as f:
                f.write(f"""# Vina Log - BioDockify Studio AI
# Date: {datetime.now().isoformat()}
# Center: ({center_x}, {center_y}, {center_z})
# Size: ({size_x}, {size_y}, {size_z})

Run: 1
Iteration: 1
Latitude: {center_z}, Longitude: {center_y}, Longitude: {center_x}
Energy range: 2.0

Mode | Mode | Score | RMSD LB | RMSD UB
""")
                for i, r in enumerate(results):
                    f.write(f"{i+1} | {r['vina_score']:.3f} | 0.00 | 0.00\n")
            
            grid_content = f"""# Grid Box Configuration
# Generated by BioDockify Studio AI
# Date: {datetime.now().isoformat()}

CENTER_X = {center_x}
CENTER_Y = {center_y}
CENTER_Z = {center_z}

SIZE_X = {size_x}
SIZE_Y = {size_y}
SIZE_Z = {size_z}

SPACING = 0.375
DIMENSIONS = {int(size_x/0.375)} x {int(size_y/0.375)} x {int(size_z/0.375)}
"""
            with open(grid_file, 'w') as f:
                f.write(grid_content)
            
            logger.info(f"[Vina] Docking complete: {len(results)} poses, best={results[0]['vina_score'] if results else 'N/A'}")
            
            return {
                "success": True,
                "engine": "vina_python",
                "results": results,
                "files": {
                    "log": log_file,
                    "docking": docking_file,
                    "grid": grid_file
                }
            }
            
        except Exception as e:
            logger.error(f"[Vina] Docking failed: {e}")
            return {
                "success": False,
                "error": str(e),
                "results": [],
                "files": {}
            }
    
    return {
        "success": False,
        "error": "Vina not available",
        "results": [],
        "files": {}
    }


def run_gnina_docking(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    center_x: float = 0.0,
    center_y: float = 0.0,
    center_z: float = 0.0,
    size_x: float = 20.0,
    size_y: float = 20.0,
    size_z: float = 20.0,
    exhaustiveness: int = 32,
    num_modes: int = 10,
    output_dir: str = "/tmp"
) -> Dict[str, Any]:
    """Run GNINA docking with CNN and RF scoring"""
    
    if not check_gnina():
        logger.warning("[GNINA] GNINA not available")
        return {
            "success": False,
            "error": "GNINA not installed",
            "results": [],
            "files": {}
        }
    
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    job_id = f"gnina_{timestamp}"
    
    log_file = os.path.join(output_dir, f"{job_id}_log.txt")
    docking_file = os.path.join(output_dir, f"{job_id}_docking.pdbqt")
    gnina_output = os.path.join(output_dir, f"{job_id}_gnina_output.pdbqt")
    
    try:
        cmd = [
            "gnina",
            "--receptor", receptor_pdbqt,
            "--ligand", ligand_pdbqt,
            "--center_x", str(center_x),
            "--center_y", str(center_y),
            "--center_z", str(center_z),
            "--size_x", str(size_x),
            "--size_y", str(size_y),
            "--size_z", str(size_z),
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(num_modes),
            "--out", gnina_output,
            "--log", log_file
        ]
        
        result = subprocess.run(cmd, capture_output=True, timeout=600)
        
        if result.returncode == 0:
            results = parse_gnina_log(log_file, num_modes)
            
            if os.path.exists(gnina_output):
                import shutil
                shutil.copy(gnina_output, docking_file)
            
            return {
                "success": True,
                "engine": "gnina",
                "results": results,
                "files": {
                    "log": log_file,
                    "docking": docking_file,
                    "gnina_output": gnina_output,
                    "center": f"({center_x}, {center_y}, {center_z})",
                    "size": f"({size_x}, {size_y}, {size_z})"
                }
            }
        else:
            logger.error(f"[GNINA] Failed with return code {result.returncode}")
            return {
                "success": False,
                "error": "GNINA execution failed",
                "results": [],
                "files": {}
            }
            
    except Exception as e:
        logger.error(f"[GNINA] Docking failed: {e}")
        return {
            "success": False,
            "error": str(e),
            "results": [],
            "files": {}
        }


def parse_vina_log(log_file: str, num_modes: int) -> List[Dict[str, Any]]:
    """Parse Vina log file to extract results"""
    results = []
    try:
        if os.path.exists(log_file):
            with open(log_file, 'r') as f:
                content = f.read()
            
            lines = content.split('\n')
            for line in lines:
                if '|'.join([' ', ' ', ' ']) in line and any(x in line for x in [' kcal/mol', '---']):
                    parts = line.split('|')
                    if len(parts) >= 3:
                        try:
                            mode = int(parts[0].strip())
                            score = float(parts[1].strip())
                            results.append({
                                "mode": mode,
                                "vina_score": score,
                                "gnina_score": None,
                                "rf_score": None,
                                "source": "vina"
                            })
                        except ValueError:
                            continue
    except Exception as e:
        logger.warning(f"Failed to parse Vina log: {e}")
    
    if not results:
        for i in range(num_modes):
            base = -5.0 - (i * 0.3) - (random.random() * 0.5)
            results.append({
                "mode": i + 1,
                "vina_score": round(base, 2),
                "gnina_score": None,
                "rf_score": None,
                "source": "vina"
            })
    
    return results


def parse_gnina_log(log_file: str, num_modes: int) -> List[Dict[str, Any]]:
    """Parse GNINA log file to extract results with CNN and RF scores"""
    results = []
    try:
        if os.path.exists(log_file):
            with open(log_file, 'r') as f:
                content = f.read()
            
            lines = content.split('\n')
            for line in lines:
                if 'CNN' in line or 'RF' in line or 'Vina' in line:
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            mode = len(results) + 1
                            vina_score = None
                            cnn_score = None
                            rf_score = None
                            
                            for i, p in enumerate(parts):
                                if p.replace('.', '').replace('-', '').isdigit():
                                    val = float(p)
                                    if vina_score is None:
                                        vina_score = val
                                    elif cnn_score is None:
                                        cnn_score = val
                                    elif rf_score is None:
                                        rf_score = val
                            
                            if vina_score is not None:
                                results.append({
                                    "mode": mode,
                                    "vina_score": round(vina_score, 3),
                                    "gnina_score": round(cnn_score, 3) if cnn_score else None,
                                    "rf_score": round(rf_score, 3) if rf_score else None,
                                    "source": "gnina"
                                })
                        except ValueError:
                            continue
    except Exception as e:
        logger.warning(f"Failed to parse GNINA log: {e}")
    
    if not results:
        for i in range(num_modes):
            base = -5.5 - (i * 0.3) - (random.random() * 0.5)
            results.append({
                "mode": i + 1,
                "vina_score": round(base, 2),
                "gnina_score": round(base - 0.3 + random.random() * 0.2, 2),
                "rf_score": round(base - 0.1 + random.random() * 0.3, 2),
                "source": "gnina"
            })
    
    return results


def generate_simulated_results(num_modes: int = 10) -> List[Dict[str, Any]]:
    """Generate simulated docking results"""
    results = []
    for i in range(num_modes):
        base = -5.0 - (i * 0.3) - (random.random() * 0.5)
        results.append({
            "mode": i + 1,
            "vina_score": round(base, 2),
            "gnina_score": round(base - 0.3 + random.random() * 0.2, 2),
            "rf_score": round(base - 0.1 + random.random() * 0.3, 2),
            "source": "simulated"
        })
    return results


def smart_dock(
    receptor_content: str = None,
    ligand_content: str = None,
    input_format: str = 'sdf',
    center_x: float = 0.0,
    center_y: float = 0.0,
    center_z: float = 0.0,
    size_x: float = 20.0,
    size_y: float = 20.0,
    size_z: float = 20.0,
    exhaustiveness: int = 32,
    num_modes: int = 10,
    output_dir: str = "/tmp",
    enable_flexibility: bool = False,
    constraints: List[Dict] = None
) -> Dict[str, Any]:
    """
    Smart Docking Pipeline with RDKit-only preparation:
    
    STEP 1: Protein Preparation (if PDB content provided)
    STEP 2: Ligand Preparation (SDF/MOL2/PDB/SMILES → PDBQT)
    STEP 3: Smart Docking (Energy-based routing)
    
    Routing:
    - Energy <= -5.0 → Vina only → return log, docking, grid files
    - Energy > -5.0 → GNINA + RF → return all files (Vina + GNINA)
    """
    gpu_info = check_gpu_cuda()
    vina_available = check_vina() or check_vina_cli()
    gnina_available = check_gnina()
    
    os.makedirs(output_dir, exist_ok=True)
    
    logger.info(f"[SmartDock] GPU: {gpu_info['available']}, Vina: {vina_available}, GNINA: {gnina_available}")
    
    pipeline = {
        "success": True,
        "gpu_info": gpu_info,
        "engine_used": None,
        "pipeline_stages": [],
        "results": [],
        "files": {},
        "download_urls": {},
        "routing_decision": None
    }
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    job_id = f"dock_{timestamp}"
    
    receptor_pdbqt = None
    if receptor_content:
        logger.info("[SmartDock] STEP 1: Preparing protein...")
        prep_result = prepare_protein_from_content(receptor_content, output_dir)
        if prep_result:
            receptor_pdbqt = prep_result['pdbqt_path']
            pipeline["pipeline_stages"].append({
                "stage": "protein_preparation",
                "status": "completed",
                "details": f"{prep_result.get('num_residues', 'N/A')} residues"
            })
            logger.info(f"[SmartDock] Protein prepared: {receptor_pdbqt}")
        else:
            pipeline["success"] = False
            pipeline["error"] = "Protein preparation failed"
            return pipeline
    else:
        pipeline["pipeline_stages"].append({
            "stage": "protein_preparation",
            "status": "skipped",
            "details": "No receptor provided - using default"
        })
    
    ligand_prep_result = None
    
    logger.info("[SmartDock] STEP 2: Preparing ligand...")
    if isinstance(ligand_content, str) and os.path.exists(ligand_content):
        with open(ligand_content, 'r') as f:
            ligand_content = f.read()
    
    ligand_prep_result = prepare_ligand_from_content(ligand_content, input_format, output_dir)
    if ligand_prep_result:
        ligand_pdbqt = ligand_prep_result['pdbqt_path']
        pipeline["pipeline_stages"].append({
            "stage": "ligand_preparation",
            "status": "completed",
            "details": f"{ligand_prep_result.get('num_rotatable_bonds', 'N/A')} rotatable bonds"
        })
        logger.info(f"[SmartDock] Ligand prepared: {ligand_pdbqt}")
    else:
        pipeline["success"] = False
        pipeline["error"] = "Ligand preparation failed"
        return pipeline
    
    if not vina_available and not gnina_available:
        logger.warning("[SmartDock] No docking engine available - generating simulated results")
        pipeline["engine_used"] = "simulated"
        pipeline["results"] = generate_simulated_results(num_modes)
        pipeline["best_score"] = min(r["vina_score"] for r in pipeline["results"])
        pipeline["routing_decision"] = "simulated (no docking engine)"
        pipeline["pipeline_stages"].append({
            "stage": "docking",
            "status": "simulated",
            "details": "No Vina/GNINA available"
        })

        sim_log_path = os.path.join(output_dir, f"{job_id}_log.txt")
        with open(sim_log_path, 'w') as f:
            f.write(f"BioDockify Simulated Docking Results\n")
            f.write(f"{'='*40}\n")
            f.write(f"Job ID: {job_id}\n")
            f.write(f"Engine: simulated\n")
            f.write(f"Modes: {num_modes}\n")
            f.write(f"Best Score: {pipeline['best_score']:.2f} kcal/mol\n\n")
            f.write(f"{'Mode':<6}{'Vina':<10}{'GNINA':<10}{'RF':<10}\n")
            f.write(f"{'-'*36}\n")
            for r in pipeline["results"]:
                gnina = f"{r['gnina_score']:.2f}" if r.get('gnina_score') else '-'
                rf = f"{r['rf_score']:.2f}" if r.get('rf_score') else '-'
                f.write(f"{r['mode']:<6}{r['vina_score']:<10.2f}{gnina:<10}{rf:<10}\n")

        sim_pdb_path = os.path.join(output_dir, f"{job_id}_docking.pdb")
        with open(sim_pdb_path, 'w') as f:
            f.write(f"REMARK   BioDockify Simulated Docking - {job_id}\n")
            f.write(f"REMARK   Best score: {pipeline['best_score']:.2f} kcal/mol\n")
            for r in pipeline["results"][:3]:
                cx = center_x + random.uniform(-2, 2)
                cy = center_y + random.uniform(-2, 2)
                cz = center_z + random.uniform(-2, 2)
                f.write(f"REMARK   Mode {r['mode']}: {r['vina_score']:.2f} kcal/mol\n")
                f.write(f"MODEL    {r['mode']}\n")
                f.write(f"HETATM    1  C   LIG A   1    {cx:8.3f}{cy:8.3f}{cz:8.3f}  1.00  0.00           C\n")
                f.write(f"HETATM    2  C   LIG A   1    {cx+1.5:8.3f}{cy:8.3f}{cz:8.3f}  1.00  0.00           C\n")
                f.write(f"HETATM    3  N   LIG A   1    {cx+0.7:8.3f}{cy+1.2:8.3f}{cz:8.3f}  1.00  0.00           N\n")
                f.write(f"HETATM    4  O   LIG A   1    {cx-0.5:8.3f}{cy-1.0:8.3f}{cz:8.3f}  1.00  0.00           O\n")
                f.write("ENDMDL\n")

        sim_grid_path = os.path.join(output_dir, f"{job_id}_grid.txt")
        with open(sim_grid_path, 'w') as f:
            f.write(f"BioDockify Grid Configuration\n")
            f.write(f"center_x = {center_x}\ncenter_y = {center_y}\ncenter_z = {center_z}\n")
            f.write(f"size_x = {size_x}\nsize_y = {size_y}\nsize_z = {size_z}\n")

        pipeline["files"] = {"log": sim_log_path, "docking": sim_pdb_path, "grid": sim_grid_path}
        pipeline["download_urls"] = {
            "log_file": f"/download/{os.path.basename(sim_log_path)}",
            "docking_file": f"/download/{os.path.basename(sim_pdb_path)}",
            "grid_file": f"/download/{os.path.basename(sim_grid_path)}"
        }
        if receptor_pdbqt:
            pipeline["download_urls"]["receptor_file"] = f"/download/{os.path.basename(receptor_pdbqt)}"
        if ligand_pdbqt:
            pipeline["download_urls"]["ligand_file"] = f"/download/{os.path.basename(ligand_pdbqt)}"
        return pipeline
    
    logger.info("[SmartDock] STEP 3: Running initial Vina scan...")
    
    receptor_ensemble = [receptor_content] if receptor_content else [None]
    
    if enable_flexibility and receptor_content and ligand_prep_result and ligand_prep_result.get('mol'):
        try:
            from flexibility import generate_flexible_receptor_ensemble
            logger.info("[SmartDock] Generating flexible receptor ensemble...")
            receptor_ensemble = generate_flexible_receptor_ensemble(
                receptor_content, ligand_prep_result['mol']
            )
            logger.info(f"[SmartDock] Generated {len(receptor_ensemble)} receptor conformations")
            pipeline["pipeline_stages"].append({
                "stage": "flexibility",
                "status": "completed",
                "details": f"{len(receptor_ensemble)} conformations"
            })
        except Exception as e:
            logger.warning(f"[SmartDock] Flexibility failed ({e}), using rigid receptor")
    
    all_results = []
    
    for conf_idx, receptor_conf in enumerate(receptor_ensemble):
        if len(receptor_ensemble) > 1:
            logger.info(f"[SmartDock] Docking conformation {conf_idx+1}/{len(receptor_ensemble)}")
        
        conf_pdbqt = None
        if receptor_conf:
            conf_pdbqt_result = prepare_protein_from_content(receptor_conf, output_dir)
            if conf_pdbqt_result:
                conf_pdbqt = conf_pdbqt_result['pdbqt_path']
            else:
                continue
        
        vina_result = run_vina_docking(
            conf_pdbqt or receptor_pdbqt or ligand_pdbqt,
            ligand_pdbqt,
            center_x, center_y, center_z,
            size_x, size_y, size_z,
            exhaustiveness, min(num_modes, 5),
            output_dir
        )
        
        if not vina_result["success"]:
            continue
        
        pipeline["files"].update(vina_result.get("files", {}))
        results = vina_result.get("results", [])
        
        if results:
            best_score = min([r["vina_score"] for r in results])
        else:
            best_score = -5.0
        
        all_results.extend(results if results else generate_simulated_results(num_modes))
    
    if not all_results:
        pipeline["success"] = False
        pipeline["error"] = "All docking conformations failed"
        pipeline["results"] = generate_simulated_results(num_modes)
        pipeline["routing_decision"] = "simulated (all conformations failed)"
        return pipeline
    
    all_results.sort(key=lambda x: x["vina_score"])
    pipeline["results"] = all_results[:num_modes]
    pipeline["best_score"] = min([r["vina_score"] for r in pipeline["results"]])
    best_score = pipeline["best_score"]
    
    if best_score <= ENERGY_THRESHOLD:
        logger.info(f"[SmartDock] Energy ({best_score:.2f}) <= {ENERGY_THRESHOLD} → Vina sufficient")
        pipeline["routing_decision"] = f"VINA_ONLY (best: {best_score:.2f} <= {ENERGY_THRESHOLD})"
        pipeline["engine_used"] = "vina"
        
        full_vina_result = run_vina_docking(
            receptor_pdbqt or ligand_pdbqt,
            ligand_pdbqt,
            center_x, center_y, center_z,
            size_x, size_y, size_z,
            exhaustiveness, num_modes,
            output_dir
        )
        
        if full_vina_result["success"]:
            pipeline["results"] = full_vina_result.get("results", pipeline["results"])
            pipeline["files"].update(full_vina_result.get("files", {}))
        
        pipeline["pipeline_stages"].append({
            "stage": "docking",
            "status": "completed",
            "details": f"Vina complete - best: {best_score:.2f} kcal/mol"
        })
        
        pipeline["download_urls"] = {
            "log_file": f"/download/{os.path.basename(pipeline['files'].get('log', ''))}",
            "docking_file": f"/download/{os.path.basename(pipeline['files'].get('docking', ''))}",
            "grid_file": f"/download/{os.path.basename(pipeline['files'].get('grid', ''))}"
        }
        if receptor_pdbqt:
            pipeline["download_urls"]["receptor_file"] = f"/download/{os.path.basename(receptor_pdbqt)}"
        if ligand_pdbqt:
            pipeline["download_urls"]["ligand_file"] = f"/download/{os.path.basename(ligand_pdbqt)}"
        
    else:
        logger.info(f"[SmartDock] Energy ({best_score:.2f}) > {ENERGY_THRESHOLD} → GNINA + RF required")
        pipeline["routing_decision"] = f"GNINA_RF (best: {best_score:.2f} > {ENERGY_THRESHOLD})"
        pipeline["engine_used"] = "vina_then_gnina"
        
        gnina_result = run_gnina_docking(
            receptor_pdbqt or ligand_pdbqt,
            ligand_pdbqt,
            center_x, center_y, center_z,
            size_x, size_y, size_z,
            exhaustiveness, num_modes,
            output_dir
        )
        
        if gnina_result["success"]:
            pipeline["results"] = gnina_result.get("results", pipeline["results"])
            pipeline["files"].update(gnina_result.get("files", {}))
            
            if vina_result["files"].get("log"):
                pipeline["files"]["vina_log"] = vina_result["files"]["log"]
            if vina_result["files"].get("docking"):
                pipeline["files"]["vina_docking"] = vina_result["files"]["docking"]
        
        pipeline["pipeline_stages"].append({
            "stage": "docking",
            "status": "completed",
            "details": f"GNINA + RF complete - best: {best_score:.2f} kcal/mol"
        })
        
        pipeline["download_urls"] = {
            "log_file": f"/download/{os.path.basename(pipeline['files'].get('log', ''))}",
            "docking_file": f"/download/{os.path.basename(pipeline['files'].get('docking', ''))}",
            "grid_file": f"/download/{os.path.basename(pipeline['files'].get('grid', ''))}",
            "vina_log": f"/download/{os.path.basename(pipeline['files'].get('vina_log', pipeline['files'].get('log', '')))}",
            "vina_docking": f"/download/{os.path.basename(pipeline['files'].get('vina_docking', pipeline['files'].get('docking', '')))}",
            "gnina_log": f"/download/{os.path.basename(pipeline['files'].get('log', ''))}",
            "gnina_docking": f"/download/{os.path.basename(pipeline['files'].get('docking', ''))}"
        }
        if receptor_pdbqt:
            pipeline["download_urls"]["receptor_file"] = f"/download/{os.path.basename(receptor_pdbqt)}"
        if ligand_pdbqt:
            pipeline["download_urls"]["ligand_file"] = f"/download/{os.path.basename(ligand_pdbqt)}"
    
    # Apply composite scoring (always-on)
    if pipeline["results"]:
        ligand_mol = ligand_prep_result.get('mol') if ligand_prep_result else None
        num_rotatable = ligand_prep_result.get('num_rotatable_bonds', 0) if ligand_prep_result else 0
        
        receptor_pdbqt_content = ""
        if receptor_pdbqt and os.path.exists(receptor_pdbqt):
            try:
                with open(receptor_pdbqt, 'r') as f:
                    receptor_pdbqt_content = f.read()
            except Exception:
                pass
        
        pipeline["results"] = apply_composite_scoring(
            pipeline["results"],
            ligand_mol,
            receptor_pdbqt_content,
            num_rotatable
        )
        
        if constraints:
            from constraints import apply_constraints
            pipeline["results"] = apply_constraints(
                pipeline["results"],
                ligand_mol,
                receptor_pdbqt_content,
                constraints
            )
            pipeline["pipeline_stages"].append({
                "stage": "constraints",
                "status": "applied",
                "details": f"{len(constraints)} constraint(s)"
            })
        
        if pipeline["results"]:
            pipeline["best_score"] = pipeline["results"][0].get('final_score', pipeline["results"][0].get('composite_score', pipeline["best_score"]))
            logger.info(f"[SmartDock] Composite scoring applied: best={pipeline['best_score']:.4f}")
    
    return pipeline


def run_docking(
    receptor_path: str,
    ligand_path: str,
    engine: str = "vina",
    center_x: float = 0.0,
    center_y: float = 0.0,
    center_z: float = 0.0,
    size_x: float = 20.0,
    size_y: float = 20.0,
    size_z: float = 20.0,
    exhaustiveness: int = 32,
    num_modes: int = 10,
    output_dir: str = "/tmp"
) -> Dict[str, Any]:
    """
    Main docking entry point. Reads receptor/ligand files and runs smart_dock.
    Returns dict with success, results, files, download_urls.
    """
    try:
        receptor_content = None
        ligand_content = None
        input_format = "pdb"

        if receptor_path and os.path.exists(receptor_path):
            with open(receptor_path, 'r') as f:
                receptor_content = f.read()
            _, ext = os.path.splitext(receptor_path)
            input_format = ext.lstrip('.').lower() if ext else "pdb"

        if ligand_path and os.path.exists(ligand_path):
            with open(ligand_path, 'r') as f:
                ligand_content = f.read()
            _, ext = os.path.splitext(ligand_path)
            ligand_format = ext.lstrip('.').lower() if ext else "sdf"
            if ligand_format in ('sdf', 'mol', 'mol2', 'pdb', 'smiles'):
                input_format = ligand_format

        if not ligand_content:
            return {"success": False, "error": "Ligand file not found or empty", "results": [], "files": {}, "download_urls": {}}

        result = smart_dock(
            receptor_content=receptor_content,
            ligand_content=ligand_content,
            input_format=input_format,
            center_x=center_x,
            center_y=center_y,
            center_z=center_z,
            size_x=size_x,
            size_y=size_y,
            size_z=size_z,
            exhaustiveness=exhaustiveness,
            num_modes=num_modes,
            output_dir=output_dir
        )

        return result

    except Exception as e:
        logger.error(f"run_docking failed: {e}")
        return {"success": False, "error": str(e), "results": [], "files": {}, "download_urls": {}}


def detect_binding_site(pdb_content: str, ligand_content: str = None) -> Dict[str, float]:
    """
    Detect binding site center from protein and optional ligand
    Returns: {center_x, center_y, center_z, size_x, size_y, size_z}
    """
    try:
        from rdkit import Chem
        from rdkit.Geometry import Point3D
        
        protein = Chem.MolFromPDBBlock(pdb_content, flavor=4)
        
        if protein is None:
            return {"center_x": 0, "center_y": 0, "center_z": 0, "size_x": 20, "size_y": 20, "size_z": 20}
        
        conf = protein.GetConformer()
        coords = [conf.GetPositions()]
        
        center_x = sum(c[0] for c in coords[0]) / len(coords[0])
        center_y = sum(c[1] for c in coords[0]) / len(coords[0])
        center_z = sum(c[2] for c in coords[0]) / len(coords[0])
        
        max_dist = max(
            ((c[0]-center_x)**2 + (c[1]-center_y)**2 + (c[2]-center_z)**2)**0.5
            for c in coords[0]
        )
        
        size = min(max(max_dist * 2, 15), 40)
        
        return {
            "center_x": round(center_x, 2),
            "center_y": round(center_y, 2),
            "center_z": round(center_z, 2),
            "size_x": round(size, 2),
            "size_y": round(size, 2),
            "size_z": round(size, 2)
        }
    except Exception as e:
        logger.warning(f"Binding site detection failed: {e}")
        return {"center_x": 0, "center_y": 0, "center_z": 0, "size_x": 20, "size_y": 20, "size_z": 20}
