"""
Agent Zero Brain - Pipeline Decision Engine
Handles GPU detection, pipeline planning, and user prompts.
"""

import logging
from typing import Dict, Any, List, Optional, Callable
from dataclasses import dataclass
from enum import Enum

from src.agent_zero.resource_guard import ResourceGuard, GPUCapabilities, VRAMLevel

logger = logging.getLogger(__name__)


class PipelineStep(Enum):
    """Pipeline step identifiers"""
    VINA_DOCKING = "vina_docking"
    GNINA_RESCORING = "gnina_rescoring"
    RF_PREDICTION = "rf_prediction"
    CONSENSUS_SCORING = "consensus_scoring"
    AI_RANKING = "ai_ranking"
    MMGBSA_REFINEMENT = "mmgbsa_refinement"
    EXPORT_RESULTS = "export_results"


@dataclass
class DockingRequest:
    """User docking request"""
    receptor_file: str
    ligand_file: str
    center_x: float
    center_y: float
    center_z: float
    size_x: float = 22.0
    size_y: float = 22.0
    size_z: float = 22.0
    exhaustiveness: int = 8
    num_modes: int = 9
    job_id: str = "default"
    output_dir: str = "data"


@dataclass
class PipelinePlan:
    """Planned execution pipeline"""
    steps: List[PipelineStep]
    engine_type: str
    use_gpu: bool
    use_gnina: bool
    use_rf: bool
    use_consensus: bool
    use_ai_ranking: bool
    use_mmgbsa: bool
    parameters: Dict[str, Any]
    warnings: List[str]


@dataclass
class UserDecision:
    """User decision on a prompt"""
    question: str
    answer: bool
    reason: str = ""


class AgentZeroBrain:
    """
    Intelligent brain for pipeline decision making.
    Checks hardware, plans pipeline, prompts user for key decisions.
    """
    
    def __init__(self, resource_guard: Optional[ResourceGuard] = None):
        self.resource_guard = resource_guard or ResourceGuard()
        self.user_decisions: Dict[str, UserDecision] = {}
        self.current_plan: Optional[PipelinePlan] = None
    
    def analyze_hardware(self) -> GPUCapabilities:
        """Analyze available hardware"""
        return self.resource_guard.get_gpu_capabilities()
    
    def create_pipeline_plan(
        self,
        request: DockingRequest,
        user_callback: Optional[Callable[[str, str], bool]] = None
    ) -> PipelinePlan:
        """
        Create pipeline plan based on hardware and user preferences.
        
        Args:
            request: User docking request
            user_callback: Optional callback for user prompts (question, details) -> bool
            
        Returns:
            PipelinePlan: Executable pipeline plan
        """
        logger.info("Analyzing hardware and creating pipeline plan")
        
        caps = self.analyze_hardware()
        
        use_gpu = False
        use_gnina = False
        use_mmgbsa = False
        warnings = []
        
        if caps.supports_vina_gpu:
            if user_callback:
                question = self._build_gpu_question(caps)
                answer = user_callback(question, self._build_gpu_details(caps))
                if answer:
                    use_gpu = True
                    self.user_decisions["use_gpu"] = UserDecision(
                        question=question, answer=True
                    )
                else:
                    warnings.append("User chose CPU mode")
                    self.user_decisions["use_gpu"] = UserDecision(
                        question=question, answer=False, reason="User declined GPU"
                    )
            else:
                use_gpu = True
        else:
            warnings.append("GPU not available - using CPU")
        
        if use_gpu and caps.supports_gnina:
            use_gnina = True
        elif not use_gpu:
            warnings.append("GNINA requires GPU - skipping CNN rescoring")
        
        if caps.supports_mmgbsa:
            if user_callback:
                question = "Enable MM-GBSA refinement for top 5 poses?"
                answer = user_callback(question, self._build_mmgbsa_details(caps))
                if answer:
                    use_mmgbsa = True
                    self.user_decisions["use_mmgbsa"] = UserDecision(
                        question=question, answer=True
                    )
                else:
                    self.user_decisions["use_mmgbsa"] = UserDecision(
                        question=question, answer=False, reason="User declined"
                    )
            else:
                use_mmgbsa = True
        elif caps.total_vram_mb > 0:
            warnings.append(f"MM-GBSA requires {VRAM_ENABLE_MMGBSA//1024}+GB GPU - skipping")
        
        steps = [PipelineStep.VINA_DOCKING]
        
        if use_gnina:
            steps.append(PipelineStep.GNINA_RESCORING)
        
        steps.extend([
            PipelineStep.RF_PREDICTION,
            PipelineStep.CONSENSUS_SCORING,
            PipelineStep.AI_RANKING
        ])
        
        if use_mmgbsa:
            steps.append(PipelineStep.MMGBSA_REFINEMENT)
        
        steps.append(PipelineStep.EXPORT_RESULTS)
        
        engine_type = "vina-gpu" if use_gpu else "vina-cpu"
        
        parameters = {
            "exhaustiveness": caps.recommended_exhaustiveness,
            "num_modes": request.num_modes,
            "cnn_batch_size": caps.recommended_cnn_batch if use_gnina else 0,
            "mmgbsa_top_poses": 5 if use_mmgbsa else 0
        }
        
        plan = PipelinePlan(
            steps=steps,
            engine_type=engine_type,
            use_gpu=use_gpu,
            use_gnina=use_gnina,
            use_rf=True,
            use_consensus=True,
            use_ai_ranking=True,
            use_mmgbsa=use_mmgbsa,
            parameters=parameters,
            warnings=warnings
        )
        
        self.current_plan = plan
        logger.info(f"Pipeline plan created: {len(steps)} steps, engine={engine_type}")
        
        return plan
    
    def _build_gpu_question(self, caps: GPUCapabilities) -> str:
        """Build GPU usage question"""
        return f"GPU detected ({caps.total_vram_mb}MB). Use GPU acceleration for faster docking?"
    
    def _build_gpu_details(self, caps: GPUCapabilities) -> str:
        """Build GPU details for user"""
        details = [
            f"GPU: {caps.total_vram_mb}MB VRAM",
            f"Level: {caps.vram_level.value.upper()}",
            f"Vina-GPU: {'✓' if caps.supports_vina_gpu else '✗'}",
            f"GNINA CNN: {'✓' if caps.supports_gnina else '✗'}",
            f"MM-GBSA: {'✓' if caps.supports_mmgbsa else '✗'}",
            f"Recommended exhaustiveness: {caps.recommended_exhaustiveness}"
        ]
        return "\n".join(details)
    
    def _build_mmgbsa_details(self, caps: GPUCapabilities) -> str:
        """Build MM-GBSA details for user"""
        return f"Requires ~2GB VRAM\nRecommended batch size: {caps.recommended_cnn_batch}\nWill refine top 5 poses"
    
    def get_status_for_ui(self) -> Dict[str, Any]:
        """Get status summary for UI display"""
        caps = self.analyze_hardware()
        
        return {
            "gpu_available": caps.supports_vina_gpu,
            "gpu_name": f"{caps.total_vram_mb}MB GPU" if caps.total_vram_mb > 0 else "None",
            "vram_level": caps.vram_level.value,
            "supports_gnina": caps.supports_gnina,
            "supports_mmgbsa": caps.supports_mmgbsa,
            "recommended_exhaustiveness": caps.recommended_exhaustiveness,
            "current_plan": {
                "steps": [s.value for s in self.current_plan.steps] if self.current_plan else [],
                "engine_type": self.current_plan.engine_type if self.current_plan else None
            } if self.current_plan else None,
            "user_decisions": {
                k: {"answer": v.answer, "reason": v.reason}
                for k, v in self.user_decisions.items()
            }
        }
    
    def reset(self) -> None:
        """Reset brain state"""
        self.user_decisions = {}
        self.current_plan = None
        logger.info("Agent Zero brain reset")


def create_default_callback() -> Callable[[str, str], bool]:
    """
    Create a default user callback that returns True (auto-accept).
    Replace this with actual UI callback in production.
    """
    def default_callback(question: str, details: str) -> bool:
        logger.info(f"User prompt: {question}")
        logger.debug(f"Details: {details}")
        return True
    
    return default_callback
