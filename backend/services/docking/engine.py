import subprocess
import os
import uuid
from pathlib import Path
from typing import Optional, List, Tuple
from rdkit import Chem

class DockingService:
    def __init__(self, work_dir: Path):
        self.work_dir = work_dir
        self.vina_path = "/usr/bin/autodock_vina"  # Installed via apt
    
    def prepare_ligand_from_smiles(self, smiles: str, job_id: str) -> Optional[str]:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            mol_h = Chem.AddHs(mol)
            pdb_path = self.work_dir / job_id / "ligand.pdb"
            pdb_path.parent.mkdir(parents=True, exist_ok=True)
            Chem.MolToPDBFile(mol_h, str(pdb_path))
            
            pdbqt_path = self.work_dir / job_id / "ligand.pdbqt"
            # In production, use obabel or prepare_ligand4.py
            return str(pdbqt_path)
        except Exception as e:
            print(f"Error preparing ligand from SMILES: {e}")
            return None
    
    def prepare_receptor(self, receptor_path: str, job_id: str) -> Optional[str]:
        try:
            job_dir = self.work_dir / job_id
            job_dir.mkdir(parents=True, exist_ok=True)
            
            receptor_name = Path(receptor_path).name
            dest_path = job_dir / receptor_name
            
            if Path(receptor_path).exists():
                import shutil
                shutil.copy(receptor_path, dest_path)
            
            # Convert to PDBQT if needed
            pdbqt_path = dest_path.with_suffix('.pdbqt')
            if dest_path.suffix == '.pdb':
                # Use obabel or prepare_receptor4.py in production
                pass
            
            return str(pdbqt_path)
        except Exception as e:
            print(f"Error preparing receptor: {e}")
            return None
    
    def run_docking(
        self,
        receptor_pdbqt: str,
        ligand_pdbqt: str,
        job_id: str,
        center: Tuple[float, float, float] = (0, 0, 0),
        size: Tuple[float, float, float] = (20, 20, 20),
        exhaustiveness: int = 32,
        n_poses: int = 10
    ) -> List[dict]:
        try:
            job_dir = self.work_dir / job_id
            output_path = job_dir / "docking_results.pdbqt"
            
            cmd = [
                self.vina_path,
                "--receptor", receptor_pdbqt,
                "--ligand", ligand_pdbqt,
                "--center_x", str(center[0]),
                "--center_y", str(center[1]),
                "--center_z", str(center[2]),
                "--size_x", str(size[0]),
                "--size_y", str(size[1]),
                "--size_z", str(size[2]),
                "--exhaustiveness", str(exhaustiveness),
                "--num_modes", str(n_poses),
                "--out", str(output_path),
                "--log", str(job_dir / "vina.log")
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            
            if result.returncode != 0:
                print(f"Vina error: {result.stderr}")
                return []
            
            return self.parse_vina_log(job_dir / "vina.log")
        except subprocess.TimeoutExpired:
            print("Docking timeout")
            return []
        except Exception as e:
            print(f"Docking error: {e}")
            return []
    
    def parse_vina_log(self, log_path: Path) -> List[dict]:
        poses = []
        try:
            content = log_path.read_text()
            lines = content.split('\n')
            for line in lines:
                if line.startswith('----') or 'Ki' in line or 'kcal' in line:
                    parts = line.split()
                    if len(parts) >= 3 and parts[0].replace('.', '').replace('-', '').isdigit():
                        try:
                            pose = {
                                'rank': len(poses) + 1,
                                'binding_energy': float(parts[0]),
                                'ki': ' '.join(parts[1:]) if len(parts) > 1 else ''
                            }
                            poses.append(pose)
                        except ValueError:
                            continue
        except Exception as e:
            print(f"Error parsing log: {e}")
        return poses
