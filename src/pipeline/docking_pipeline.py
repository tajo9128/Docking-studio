"""
Docking Pipeline
Orchestrates the complete 6-step docking workflow.
"""

import logging
from typing import Dict, Any, List, Optional
from pathlib import Path

from src.core.engines.factory import DockingEngineFactory, EngineType
from src.core.engines.base import DockingConfig
from src.ml.rf_scoring.consensus_scorer import ConsensusScorer
from src.hardware.gpu_detector import has_gpu, get_gpu_info

logger = logging.getLogger(__name__)


class DockingPipeline:
    """
    Complete docking pipeline with intelligent hardware selection.
    
    Steps:
    1. Input Validation & Preparation
    2. Engine Selection (auto GPU/CPU)
    3. Docking Execution
    4. (Optional) GNINA Rescoring
    5. (Optional) RF Consensus Scoring
    6. Result Aggregation
    """
    
    def __init__(self, docker_manager=None):
        self.docker_manager = docker_manager
        self.factory = DockingEngineFactory(docker_manager)
        self.consensus_scorer = ConsensusScorer()
    
    def run(
        self,
        receptor_file: str,
        ligand_file: str,
        center_x: float,
        center_y: float,
        center_z: float,
        size_x: float = 22.0,
        size_y: float = 22.0,
        size_z: float = 22.0,
        exhaustiveness: int = 8,
        num_modes: int = 9,
        job_id: str = "default",
        output_dir: str = "data",
        engine_type: Optional[str] = None,
        use_rf: bool = True,
        use_gnina: bool = True
    ) -> Dict[str, Any]:
        """
        Run complete docking pipeline.
        
        Args:
            receptor_file: Path to receptor PDBQT
            ligand_file: Path to ligand PDBQT
            center_x/y/z: Grid box center
            size_x/y/z: Grid box size
            exhaustiveness: Vina exhaustiveness
            num_modes: Number of poses to generate
            job_id: Job identifier
            output_dir: Output directory
            engine_type: Specific engine or None for auto
            use_rf: Use RF consensus scoring
            use_gnina: Use GNINA rescoring
            
        Returns:
            dict: Complete results with all scores
        """
        logger.info(f"Starting docking pipeline for job {job_id}")
        
        result = {
            "status": "PENDING",
            "job_id": job_id,
            "engine_used": None,
            "gpu_info": None,
            "poses": [],
            "consensus_poses": [],
            "best_pose": None
        }
        
        result["gpu_info"] = self._get_hardware_info()
        
        config = DockingConfig(
            receptor_file=receptor_file,
            ligand_file=ligand_file,
            center_x=center_x,
            center_y=center_y,
            center_z=center_z,
            size_x=size_x,
            size_y=size_y,
            size_z=size_z,
            exhaustiveness=exhaustiveness,
            num_modes=num_modes,
            job_id=job_id,
            output_dir=output_dir
        )
        
        if engine_type is None:
            engine_type = self._select_engine(use_gnina)
        
        engine, engine_desc = self.factory.create_engine(
            engine_type=engine_type,
            prefer_cnn=use_gnina
        )
        
        result["engine_used"] = engine_desc
        logger.info(f"Using engine: {engine_desc}")
        
        try:
            docking_result = engine.run(config)
            
            result["status"] = docking_result.status
            result["poses"] = docking_result.poses
            result["logs"] = docking_result.logs
            
            if docking_result.status == "COMPLETED":
                if use_rf:
                    result["consensus_poses"] = self.consensus_scorer.score_poses(
                        poses=docking_result.poses,
                        use_rf=True,
                        receptor_path=receptor_file,
                        ligand_path=ligand_file
                    )
                    result["best_pose"] = result["consensus_poses"][0] if result["consensus_poses"] else None
                else:
                    result["consensus_poses"] = docking_result.poses
                    result["best_pose"] = result["poses"][0] if result["poses"] else None
                
                result["binding_energy"] = docking_result.binding_energy
                
                if hasattr(docking_result, 'gnina_cnn_score'):
                    result["gnina_cnn_score"] = docking_result.gnina_cnn_score
                    result["gnina_cnn_affinity"] = docking_result.gnina_cnn_affinity
                
                logger.info(f"Docking completed successfully for job {job_id}")
            else:
                logger.error(f"Docking failed for job {job_id}")
                
        except Exception as e:
            logger.error(f"Pipeline error: {e}")
            result["status"] = "ERROR"
            result["error"] = str(e)
        
        return result
    
    def _select_engine(self, use_gnina: bool) -> str:
        """Select best engine based on hardware"""
        if not has_gpu():
            return EngineType.VINA_CPU
        
        if use_gnina:
            return EngineType.GNINA
        
        return EngineType.VINA_GPU
    
    def _get_hardware_info(self) -> Dict[str, Any]:
        """Get hardware information"""
        if has_gpu():
            info = get_gpu_info()
            if info:
                return {
                    "gpu_available": True,
                    "gpu_name": info.get("name"),
                    "vram_mb": info.get("vram_mb"),
                    "cuda_version": info.get("cuda_version")
                }
        
        return {"gpu_available": False}
    
    def run_batch(
        self,
        receptor_file: str,
        ligand_files: List[str],
        **kwargs
    ) -> List[Dict[str, Any]]:
        """
        Run docking for multiple ligands.
        
        Args:
            receptor_file: Path to receptor
            ligand_files: List of ligand paths
            **kwargs: Other arguments passed to run()
            
        Returns:
            list: Results for each ligand
        """
        results = []
        
        for idx, ligand_file in enumerate(ligand_files):
            job_id = kwargs.get("job_id", "batch")
            kwargs["job_id"] = f"{job_id}_{idx}"
            
            result = self.run(
                receptor_file=receptor_file,
                ligand_file=ligand_file,
                **kwargs
            )
            results.append(result)
        
        return results
