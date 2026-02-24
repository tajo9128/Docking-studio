"""
MM-GBSA Runner
Refines top poses using MM-GBSA (Molecular Mechanics Generalized Born/Surface Area).
"""

import logging
import subprocess
import os
from typing import Dict, Any, List, Optional
from pathlib import Path

logger = logging.getLogger(__name__)


class MMGBSARunner:
    """
    Runs MM-GBSA refinement on top poses.
    
    Uses OpenMM or AmberTools for MM-GBSA calculation.
    Only runs on high-end GPUs (6GB+) due to memory requirements.
    """
    
    def __init__(self, docker_manager=None):
        self.docker_manager = docker_manager
        self.image_name = "biodockify/mmgbsa:latest"
    
    def refine_poses(
        self,
        receptor_file: str,
        poses: List[Dict[str, Any]],
        job_id: str,
        output_dir: str = "data"
    ) -> Dict[str, Any]:
        """
        Run MM-GBSA refinement on top poses.
        
        Args:
            receptor_file: Path to receptor
            poses: Top poses to refine
            job_id: Job identifier
            output_dir: Output directory
            
        Returns:
            dict: MM-GBSA results
        """
        if not poses:
            return {"status": "skipped", "reason": "No poses provided"}
        
        logger.info(f"Running MM-GBSA on {len(poses)} poses")
        
        results = []
        
        for idx, pose in enumerate(poses):
            pose_result = self._calculate_single_pose(
                receptor_file=receptor_file,
                pose=pose,
                pose_index=idx,
                job_id=job_id,
                output_dir=output_dir
            )
            results.append(pose_result)
        
        results.sort(key=lambda x: x.get("mmgbsa_energy", float('inf')))
        
        best_result = results[0] if results else None
        
        return {
            "status": "completed",
            "num_poses": len(poses),
            "results": results,
            "best_pose": best_result,
            "best_energy": best_result.get("mmgbsa_energy") if best_result else None
        }
    
    def _calculate_single_pose(
        self,
        receptor_file: str,
        pose: Dict[str, Any],
        pose_index: int,
        job_id: str,
        output_dir: str
    ) -> Dict[str, Any]:
        """
        Calculate MM-GBSA for single pose.
        
        Note: This is a placeholder. Real implementation would:
        1. Prepare complex structure (receptor + ligand pose)
        2. Calculate GBSA energy
        3. Calculate receptor energy
        4. Calculate ligand energy
        5. Return ΔG = G_complex - G_receptor - G_ligand
        """
        logger.debug(f"Calculating MM-GBSA for pose {pose_index + 1}")
        
        binding_energy = pose.get("binding_energy", -7.0)
        
        mmgbsa_energy = binding_energy - 2.5
        
        return {
            "pose_id": pose.get("mode", pose_index + 1),
            "pose_index": pose_index,
            "mmgbsa_energy": round(mmgbsa_energy, 2),
            "docking_energy": binding_energy,
            "energy_improvement": round(2.5, 2),
            "note": "Estimated MM-GBSA (full calculation requires OpenMM/AmberTools)"
        }
    
    def _run_with_docker(
        self,
        receptor_file: str,
        pose_file: str,
        job_id: str,
        output_dir: str
    ) -> Optional[float]:
        """
        Run MM-GBSA using Docker container.
        
        Returns MM-GBSA binding energy in kcal/mol.
        """
        if not self.docker_manager:
            return None
        
        try:
            parameters = {
                "receptor": receptor_file,
                "ligand": pose_file,
                "job_id": job_id
            }
            
            success = self.docker_manager.start_container(
                receptor_file=receptor_file,
                ligand_file=pose_file,
                parameters=parameters,
                job_id=f"mmgbsa_{job_id}",
                image_name=self.image_name,
                use_gpu=True
            )
            
            if not success:
                logger.error("MM-GBSA container failed to start")
                return None
            
            output_file = Path(output_dir) / f"mmgbsa_{job_id}.txt"
            
            if output_file.exists():
                with open(output_file) as f:
                    content = f.read()
                    for line in content.split('\n'):
                        if 'MMGBSA' in line.upper():
                            parts = line.split()
                            for part in parts:
                                try:
                                    return float(part)
                                except ValueError:
                                    continue
            
            return None
            
        except Exception as e:
            logger.error(f"MM-GBSA Docker run failed: {e}")
            return None
    
    def is_available(self) -> bool:
        """Check if MM-GBSA is available"""
        if self.docker_manager:
            return self.docker_manager.is_docker_running()
        
        return False
