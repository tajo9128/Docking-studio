"""
Pipeline Manager - 6-Step Docking Pipeline Orchestration
Executes the complete docking workflow with all steps.
"""

import logging
from typing import Dict, Any, List, Optional, Callable
from dataclasses import dataclass
from datetime import datetime

from src.agent_zero.brain import AgentZeroBrain, PipelinePlan, PipelineStep, DockingRequest
from src.agent_zero.resource_guard import ResourceGuard
from src.pipeline.docking_pipeline import DockingPipeline
from src.ml.rf_scoring.rf_model_service import RFModelService
from src.ml.rf_scoring.consensus_scorer import ConsensusScorer

logger = logging.getLogger(__name__)


@dataclass
class PipelineResult:
    """Result from pipeline execution"""
    status: str
    job_id: str
    steps_completed: List[str]
    poses: List[Dict[str, Any]]
    best_pose: Optional[Dict[str, Any]]
    vina_result: Optional[Dict[str, Any]] = None
    gnina_result: Optional[Dict[str, Any]] = None
    rf_results: Optional[Dict[str, Any]] = None
    consensus_results: Optional[Dict[str, Any]] = None
    ai_ranking: Optional[Dict[str, Any]] = None
    mmgbsa_result: Optional[Dict[str, Any]] = None
    export_result: Optional[Dict[str, Any]] = None
    errors: List[str] = None
    warnings: List[str] = None
    execution_time_seconds: float = 0
    
    def __post_init__(self):
        if self.errors is None:
            self.errors = []
        if self.warnings is None:
            self.warnings = []


class PipelineManager:
    """
    Orchestrates the complete 6-step docking pipeline.
    
    Steps:
    1. Vina Docking (CPU or GPU)
    2. GNINA CNN Rescoring (if GPU)
    3. RF pKd Prediction
    4. Consensus Scoring
    5. AI Pose Ranking
    6. MM-GBSA Refinement (if GPU 6GB+)
    7. Export Results
    """
    
    def __init__(self, docker_manager=None):
        self.docker_manager = docker_manager
        self.brain = AgentZeroBrain(ResourceGuard())
        self.docking_pipeline = DockingPipeline(docker_manager)
        self.rf_service = RFModelService()
        self.consensus_scorer = ConsensusScorer()
        self.progress_callback: Optional[Callable[[str, int], None]] = None
    
    def set_progress_callback(self, callback: Callable[[str, int], None]) -> None:
        """Set callback for progress updates (message, percentage)"""
        self.progress_callback = callback
    
    def _report_progress(self, message: str, percentage: int) -> None:
        """Report progress to callback or log"""
        logger.info(f"Pipeline: {message} ({percentage}%)")
        if self.progress_callback:
            self.progress_callback(message, percentage)
    
    def execute_pipeline(
        self,
        request: DockingRequest,
        user_callback: Optional[Callable[[str, str], bool]] = None
    ) -> PipelineResult:
        """
        Execute complete docking pipeline.
        
        Args:
            request: User docking request
            user_callback: Optional callback for user prompts
            
        Returns:
            PipelineResult: Complete results
        """
        import time
        start_time = time.time()
        
        result = PipelineResult(
            status="PENDING",
            job_id=request.job_id,
            steps_completed=[],
            poses=[],
            best_pose=None
        )
        
        try:
            self._report_progress("Analyzing hardware", 0)
            
            plan = self.brain.create_pipeline_plan(request, user_callback)
            result.warnings.extend(plan.warnings)
            
            self._report_progress(f"Starting Vina docking ({plan.engine_type})", 5)
            
            vina_result = self._execute_vina(request, plan)
            result.vina_result = vina_result
            result.steps_completed.append(PipelineStep.VINA_DOCKING.value)
            
            if vina_result.get("status") != "COMPLETED":
                result.status = "FAILED"
                result.errors.append("Vina docking failed")
                return result
            
            result.poses = vina_result.get("poses", [])
            
            if plan.use_gnina and PipelineStep.GNINA_RESCORING in plan.steps:
                self._report_progress("Running GNINA CNN rescoring", 25)
                gnina_result = self._execute_gnina(request, plan, result.poses)
                result.gnina_result = gnina_result
                result.steps_completed.append(PipelineStep.GNINA_RESCORING.value)
                
                if gnina_result.get("poses"):
                    result.poses = gnina_result.get("poses", result.poses)
            
            if plan.use_rf and PipelineStep.RF_PREDICTION in plan.steps:
                self._report_progress("Running RF prediction", 50)
                rf_result = self._execute_rf(request, result.poses)
                result.rf_results = rf_result
                result.steps_completed.append(PipelineStep.RF_PREDICTION.value)
            
            if plan.use_consensus and PipelineStep.CONSENSUS_SCORING in plan.steps:
                self._report_progress("Computing consensus scores", 60)
                consensus_result = self._execute_consensus(result.poses)
                result.consensus_results = consensus_result
                result.steps_completed.append(PipelineStep.CONSENSUS_SCORING.value)
                
                result.poses = consensus_result.get("scored_poses", result.poses)
            
            if plan.use_ai_ranking and PipelineStep.AI_RANKING in plan.steps:
                self._report_progress("Running AI pose ranking", 70)
                ai_result = self._execute_ai_ranking(result.poses)
                result.ai_ranking = ai_result
                result.steps_completed.append(PipelineStep.AI_RANKING.value)
            
            if plan.use_mmgbsa and PipelineStep.MMGBSA_REFINEMENT in plan.steps:
                self._report_progress("Running MM-GBSA refinement", 80)
                mmgbsa_result = self._execute_mmgbsa(request, result.poses[:5])
                result.mmgbsa_result = mmgbsa_result
                result.steps_completed.append(PipelineStep.MMGBSA_REFINEMENT.value)
            
            if PipelineStep.EXPORT_RESULTS in plan.steps:
                self._report_progress("Exporting results", 90)
                export_result = self._execute_export(request, result)
                result.export_result = export_result
                result.steps_completed.append(PipelineStep.EXPORT_RESULTS.value)
            
            result.best_pose = result.poses[0] if result.poses else None
            result.status = "COMPLETED"
            self._report_progress("Pipeline completed", 100)
            
        except Exception as e:
            logger.error(f"Pipeline error: {e}")
            result.status = "ERROR"
            result.errors.append(str(e))
        
        result.execution_time_seconds = time.time() - start_time
        return result
    
    def _execute_vina(self, request: DockingRequest, plan: PipelinePlan) -> Dict[str, Any]:
        """Execute Vina docking"""
        return self.docking_pipeline.run(
            receptor_file=request.receptor_file,
            ligand_file=request.ligand_file,
            center_x=request.center_x,
            center_y=request.center_y,
            center_z=request.center_z,
            size_x=request.size_x,
            size_y=request.size_y,
            size_z=request.size_z,
            exhaustiveness=plan.parameters.get("exhaustiveness", 8),
            num_modes=request.num_modes,
            job_id=request.job_id,
            output_dir=request.output_dir,
            engine_type=plan.engine_type,
            use_rf=False,
            use_gnina=False
        )
    
    def _execute_gnina(self, request: DockingRequest, plan: PipelinePlan, poses: List[Dict]) -> Dict[str, Any]:
        """Execute GNINA rescoring"""
        from src.core.engines.factory import get_engine, EngineType
        
        engine, _ = get_engine(
            docker_manager=self.docker_manager,
            engine_type=EngineType.GNINA
        )
        
        from src.core.engines.base import DockingConfig
        
        config = DockingConfig(
            receptor_file=request.receptor_file,
            ligand_file=request.ligand_file,
            center_x=request.center_x,
            center_y=request.center_y,
            center_z=request.center_z,
            size_x=request.size_x,
            size_y=request.size_y,
            size_z=request.size_z,
            exhaustiveness=plan.parameters.get("exhaustiveness", 8),
            num_modes=request.num_modes,
            job_id=request.job_id,
            output_dir=request.output_dir
        )
        
        return engine.run(config)
    
    def _execute_rf(self, request: DockingRequest, poses: List[Dict]) -> Dict[str, Any]:
        """Execute RF prediction"""
        predictions = []
        
        for pose in poses:
            vina_energy = pose.get("binding_energy", -7.0)
            gnina_score = pose.get("gnina_cnn_score")
            
            features = self._build_features(vina_energy, gnina_score)
            pKd = self.rf_service.predict_pkd(features)
            
            predictions.append({
                "pose_id": pose.get("mode"),
                "predicted_pKd": pKd
            })
        
        return {"predictions": predictions}
    
    def _build_features(self, vina_energy: float, gnina_score: Optional[float]) -> Any:
        """Build feature vector for RF"""
        import numpy as np
        
        features = np.zeros(54)
        
        features[0] = vina_energy
        
        if gnina_score is not None:
            features[5] = gnina_score
        
        return features
    
    def _execute_consensus(self, poses: List[Dict]) -> Dict[str, Any]:
        """Execute consensus scoring"""
        scored_poses = self.consensus_scorer.score_poses(poses, use_rf=False)
        
        return {
            "scored_poses": scored_poses,
            "best_score": scored_poses[0].get("consensus_score") if scored_poses else None
        }
    
    def _execute_ai_ranking(self, poses: List[Dict]) -> Dict[str, Any]:
        """Execute AI pose ranking (weighted ensemble)"""
        from src.agent_zero.ai.pose_ranker import PoseRanker
        
        ranker = PoseRanker()
        ranked_poses = ranker.rank_poses(poses)
        
        return {
            "ranked_poses": ranked_poses,
            "top_pose": ranked_poses[0] if ranked_poses else None
        }
    
    def _execute_mmgbsa(self, request: DockingRequest, top_poses: List[Dict]) -> Dict[str, Any]:
        """Execute MM-GBSA refinement"""
        from src.agent_zero.mmgbsa.runner import MMGBSARunner
        
        runner = MMGBSARunner(self.docker_manager)
        
        results = runner.refine_poses(
            receptor_file=request.receptor_file,
            poses=top_poses,
            job_id=request.job_id,
            output_dir=request.output_dir
        )
        
        return results
    
    def _execute_export(self, request: DockingRequest, result: PipelineResult) -> Dict[str, Any]:
        """Execute export"""
        from src.agent_zero.export.exporter import ResultExporter
        
        exporter = ResultExporter()
        
        export_paths = exporter.export_all(
            job_id=request.job_id,
            output_dir=request.output_dir,
            poses=result.poses,
            best_pose=result.best_pose,
            vina_result=result.vina_result,
            gnina_result=result.gnina_result,
            rf_result=result.rf_results,
            consensus_result=result.consensus_results,
            ai_result=result.ai_ranking,
            mmgbsa_result=result.mmgbsa_result
        )
        
        return export_paths
