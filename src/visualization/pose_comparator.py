"""
Pose Comparator
Split-view pose comparison with RMSD calculation.
"""

from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple, Callable
import numpy as np
import logging

logger = logging.getLogger(__name__)


@dataclass
class PoseInfo:
    """Information about a docking pose"""
    id: int
    pdb_content: str
    score: float
    rmsd_to_ref: Optional[float] = None
    energy: Optional[float] = None
    atoms: Optional[List[Dict]] = None


@dataclass
class ComparisonResult:
    """Result of pose comparison"""
    pose1_id: int
    pose2_id: int
    rmsd: float
    max_displacement: float
    aligned: bool


class RMSDCalculator:
    """
    Calculates RMSD (Root Mean Square Deviation) between molecular poses.
    
    Supports:
    - Heavy atom RMSD
    - All-atom RMSD
    - Carbon-alpha RMSD
    - Ligand RMSD
    """
    
    def __init__(self):
        """Initialize RMSD calculator"""
        self.reference_pose: Optional[PoseInfo] = None
        
        logger.info("RMSDCalculator initialized")
    
    def set_reference(self, pose: PoseInfo):
        """Set reference pose for RMSD calculations"""
        self.reference_pose = pose
        logger.info(f"Reference pose set: ID {pose.id}")
    
    def calculate_rmsd(
        self,
        pose1: PoseInfo,
        pose2: PoseInfo,
        mode: str = "heavy"
    ) -> float:
        """
        Calculate RMSD between two poses.
        
        Args:
            pose1: First pose
            pose2: Second pose
            mode: RMSD calculation mode ("heavy", "all", "ca", "ligand")
            
        Returns:
            RMSD value in Angstroms
        """
        atoms1 = pose1.atoms or self._parse_pdb(pose1.pdb_content)
        atoms2 = pose2.atoms or self._parse_pdb(pose2.pdb_content)
        
        if not atoms1 or not atoms2:
            logger.warning("Could not parse atoms for RMSD calculation")
            return float('inf')
        
        # Filter atoms based on mode
        filtered1 = self._filter_atoms(atoms1, mode)
        filtered2 = self._filter_atoms(atoms2, mode)
        
        if len(filtered1) != len(filtered2):
            logger.warning("Atom counts don't match for RMSD")
            return float('inf')
        
        if not filtered1:
            return float('inf')
        
        # Calculate RMSD
        coords1 = np.array([(a['x'], a['y'], a['z']) for a in filtered1])
        coords2 = np.array([(a['x'], a['y'], a['z']) for a in filtered2])
        
        # Center the structures
        center1 = coords1.mean(axis=0)
        center2 = coords2.mean(axis=0)
        
        coords1 = coords1 - center1
        coords2 = coords2 - center2
        
        # Calculate RMSD
        diff = coords1 - coords2
        rmsd = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
        
        return float(rmsd)
    
    def calculate_to_reference(
        self,
        pose: PoseInfo,
        mode: str = "heavy"
    ) -> float:
        """Calculate RMSD to reference pose"""
        if self.reference_pose is None:
            logger.warning("No reference pose set")
            return float('inf')
        
        rmsd = self.calculate_rmsd(self.reference_pose, pose, mode)
        
        # Store in pose
        pose.rmsd_to_ref = rmsd
        
        return rmsd
    
    def _filter_atoms(
        self,
        atoms: List[Dict],
        mode: str
    ) -> List[Dict]:
        """Filter atoms based on RMSD mode"""
        if mode == "heavy":
            # Non-hydrogen atoms
            return [a for a in atoms if a.get('elem', '').upper() != 'H']
        elif mode == "all":
            return atoms
        elif mode == "ca":
            # Carbon alpha atoms
            return [a for a in atoms if a.get('name', '').upper() == 'CA']
        elif mode == "ligand":
            # Ligand atoms (usually hetero atoms)
            return [a for a in atoms if a.get('resn', '') != '']
        else:
            return atoms
    
    def _parse_pdb(self, pdb_content: str) -> List[Dict]:
        """Parse PDB content to atoms"""
        atoms = []
        
        for line in pdb_content.split('\n'):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    atom = {
                        'serial': int(line[6:11].strip()),
                        'name': line[12:16].strip(),
                        'elem': line[76:78].strip() or line[12:14].strip(),
                        'x': float(line[30:38].strip()),
                        'y': float(line[38:46].strip()),
                        'z': float(line[46:54].strip()),
                        'resn': line[17:20].strip(),
                        'resi': int(line[22:26].strip()),
                    }
                    atoms.append(atom)
                except (ValueError, IndexError):
                    continue
        
        return atoms
    
    def calculate_all_vs_all(
        self,
        poses: List[PoseInfo],
        mode: str = "heavy"
    ) -> np.ndarray:
        """Calculate RMSD matrix for all poses"""
        n = len(poses)
        rmsd_matrix = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i + 1, n):
                rmsd = self.calculate_rmsd(poses[i], poses[j], mode)
                rmsd_matrix[i, j] = rmsd
                rmsd_matrix[j, i] = rmsd
        
        return rmsd_matrix


class PoseComparator:
    """
    Professional pose comparison tool with split-view support.
    
    Features:
    - Side-by-side pose comparison
    - RMSD calculation
    - Superposition
    - Difference highlighting
    """
    
    def __init__(self):
        """Initialize pose comparator"""
        self.rmsd_calculator = RMSDCalculator()
        self.poses: List[PoseInfo] = []
        self.reference_pose: Optional[PoseInfo] = None
        
        logger.info("PoseComparator initialized")
    
    def add_pose(
        self,
        pdb_content: str,
        pose_id: int,
        score: float = 0.0
    ) -> PoseInfo:
        """Add a pose to the comparator"""
        atoms = self.rmsd_calculator._parse_pdb(pdb_content)
        
        pose = PoseInfo(
            id=pose_id,
            pdb_content=pdb_content,
            score=score,
            atoms=atoms
        )
        
        self.poses.append(pose)
        
        # Set first pose as reference if none
        if self.reference_pose is None:
            self.reference_pose = pose
            self.rmsd_calculator.set_reference(pose)
        
        return pose
    
    def set_reference(self, pose_id: int):
        """Set reference pose by ID"""
        for pose in self.poses:
            if pose.id == pose_id:
                self.reference_pose = pose
                self.rmsd_calculator.set_reference(pose)
                logger.info(f"Reference pose set to ID {pose_id}")
                break
    
    def compare_poses(
        self,
        pose1_id: int,
        pose2_id: int
    ) -> ComparisonResult:
        """Compare two poses"""
        pose1 = next((p for p in self.poses if p.id == pose1_id), None)
        pose2 = next((p for p in self.poses if p.id == pose2_id), None)
        
        if not pose1 or not pose2:
            logger.warning("Pose not found for comparison")
            return None
        
        rmsd = self.rmsd_calculator.calculate_rmsd(pose1, pose2)
        
        # Calculate max displacement
        if pose1.atoms and pose2.atoms:
            coords1 = np.array([(a['x'], a['y'], a['z']) for a in pose1.atoms])
            coords2 = np.array([(a['x'], a['y'], a['z']) for a in pose2.atoms])
            
            min_len = min(len(coords1), len(coords2))
            if min_len > 0:
                displacements = np.linalg.norm(coords1[:min_len] - coords2[:min_len], axis=1)
                max_disp = float(np.max(displacements))
            else:
                max_disp = 0.0
        else:
            max_disp = 0.0
        
        return ComparisonResult(
            pose1_id=pose1_id,
            pose2_id=pose2_id,
            rmsd=rmsd,
            max_displacement=max_disp,
            aligned=rmsd < 2.0  # Typical threshold for "aligned"
        )
    
    def get_rmsd_to_reference(
        self,
        pose_id: int,
        mode: str = "heavy"
    ) -> float:
        """Get RMSD to reference for a pose"""
        pose = next((p for p in self.poses if p.id == pose_id), None)
        
        if not pose:
            return float('inf')
        
        if self.reference_pose is None:
            return float('inf')
        
        return self.rmsd_calculator.calculate_rmsd(
            self.reference_pose,
            pose,
            mode
        )
    
    def rank_poses_by_rmsd(self) -> List[Tuple[int, float]]:
        """Rank poses by RMSD to reference (best first)"""
        if self.reference_pose is None:
            return []
        
        rmsd_list = []
        for pose in self.poses:
            if pose.id == self.reference_pose.id:
                rmsd_list.append((pose.id, 0.0))
            else:
                rmsd = self.rmsd_calculator.calculate_rmsd(
                    self.reference_pose,
                    pose
                )
                rmsd_list.append((pose.id, rmsd))
        
        return sorted(rmsd_list, key=lambda x: x[1])
    
    def get_comparison_summary(self) -> Dict:
        """Get summary of all comparisons"""
        if not self.poses:
            return {'poses': 0}
        
        # Calculate RMSD matrix
        rmsd_matrix = self.rmsd_calculator.calculate_all_vs_all(self.poses)
        
        # Find best pose
        ranked = self.rank_poses_by_rmsd()
        
        return {
            'poses': len(self.poses),
            'reference_id': self.reference_pose.id if self.reference_pose else None,
            'rmsd_matrix': rmsd_matrix.tolist(),
            'best_pose_id': ranked[0][0] if ranked else None,
            'best_rmsd': ranked[0][1] if ranked else None,
            'avg_rmsd': float(np.mean(rmsd_matrix[np.triu_indices(len(self.poses), 1)])) if len(self.poses) > 1 else 0,
        }
    
    def get_viewer_config(self) -> Dict:
        """Get configuration for split-view viewer"""
        if not self.poses:
            return {}
        
        # Prepare pose data for viewer
        pose_data = []
        for pose in self.poses:
            pose_data.append({
                'id': pose.id,
                'score': pose.score,
                'rmsd': pose.rmsd_to_ref,
            })
        
        return {
            'poses': pose_data,
            'reference_id': self.reference_pose.id if self.reference_pose else None,
            'view_mode': 'split',  # or 'overlay' or 'single'
        }


class PoseOverlay:
    """Overlay multiple poses for comparison"""
    
    def __init__(self):
        self.overlays: List[Tuple[PoseInfo, str]] = []  # pose, color
    
    def add_pose_with_color(
        self,
        pose: PoseInfo,
        color: str = "#FF00FF"
    ):
        """Add pose with color for overlay"""
        self.overlays.append((pose, color))
    
    def generate_overlay_pdb(self) -> str:
        """Generate combined PDB for overlay"""
        lines = []
        
        for pose, color in self.overlays:
            # Add comment with color info
            lines.append(f"REMARK   COLOR {color}")
            lines.append(pose.pdb_content)
        
        return '\n'.join(lines)
    
    def clear(self):
        """Clear all overlays"""
        self.overlays = []
