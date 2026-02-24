"""
Docking Engine Factory
Automatically selects the best engine based on hardware availability.
"""

import logging
from typing import Optional, Tuple

from .base import DockingEngine, DockingConfig, DockingResult
from .vina_cpu import VinaCPUEngine
from .vina_gpu import VinaGPUEngine
from .gnina_engine import GninaEngine

logger = logging.getLogger(__name__)


class EngineType:
    """Engine type constants"""
    VINA_CPU = "vina-cpu"
    VINA_GPU = "vina-gpu"
    GNINA = "gnina"
    GNINA_CPU = "gnina-cpu"


class DockingEngineFactory:
    """
    Factory for creating docking engines with automatic hardware detection.
    
    Usage:
        factory = DockingEngineFactory(docker_manager)
        engine = factory.create_engine()
        
        # Or with specific engine
        engine = factory.create_engine(engine_type=EngineType.GNINA)
    """
    
    def __init__(self, docker_manager=None):
        self.docker_manager = docker_manager
    
    def get_hardware_info(self) -> dict:
        """Get current hardware information"""
        from src.hardware.gpu_detector import has_gpu, get_gpu_info, get_optimal_engine_config
        
        if has_gpu():
            gpu_info = get_gpu_info()
            if gpu_info:
                return {
                    "gpu_available": True,
                    "gpu_name": gpu_info["name"],
                    "vram_mb": gpu_info["vram_mb"],
                    "recommendation": get_optimal_engine_config()
                }
        
        return {
            "gpu_available": False,
            "recommendation": {
                "engine_type": EngineType.VINA_CPU,
                "reason": "No GPU detected"
            }
        }
    
    def create_engine(
        self,
        engine_type: Optional[str] = None,
        prefer_cnn: bool = False,
        force_cpu: bool = False
    ) -> Tuple[DockingEngine, str]:
        """
        Create appropriate docking engine.
        
        Args:
            engine_type: Specific engine to use (None for auto)
            prefer_cnn: Prefer CNN-based scoring (GNINA) when available
            force_cpu: Force CPU-only mode
            
        Returns:
            tuple: (engine, description)
        """
        if engine_type:
            return self._create_specific_engine(engine_type, force_cpu)
        
        return self._create_auto_engine(prefer_cnn, force_cpu)
    
    def _create_auto_engine(self, prefer_cnn: bool, force_cpu: bool) -> Tuple[DockingEngine, str]:
        """Auto-select engine based on hardware"""
        from src.hardware.gpu_detector import has_gpu, get_gpu_info
        
        if force_cpu:
            logger.info("CPU mode forced, using Vina CPU")
            return VinaCPUEngine(self.docker_manager), "Vina CPU (forced)"
        
        if not has_gpu():
            logger.info("No GPU detected, using Vina CPU")
            return VinaCPUEngine(self.docker_manager), "Vina CPU (no GPU)"
        
        gpu_info = get_gpu_info()
        vram = gpu_info["vram_mb"] if gpu_info else 0
        
        if vram <= 4096:
            if prefer_cnn:
                logger.info(f"4GB GPU detected, using GNINA (CNN scoring)")
                return GninaEngine(self.docker_manager, force_cpu=False), f"GNINA (GPU, {vram}MB)"
            logger.info(f"4GB GPU detected, using Vina-GPU")
            return VinaGPUEngine(self.docker_manager), f"Vina-GPU ({vram}MB)"
        
        if prefer_cnn or vram >= 8192:
            logger.info(f"GPU with {vram}MB, using GNINA (CNN scoring)")
            return GninaEngine(self.docker_manager, force_cpu=False), f"GNINA (GPU, {vram}MB)"
        
        logger.info(f"GPU with {vram}MB, using Vina-GPU")
        return VinaGPUEngine(self.docker_manager), f"Vina-GPU ({vram}MB)"
    
    def _create_specific_engine(self, engine_type: str, force_cpu: bool) -> Tuple[DockingEngine, str]:
        """Create a specific engine type"""
        
        if engine_type == EngineType.VINA_CPU:
            return VinaCPUEngine(self.docker_manager), "Vina CPU"
        
        elif engine_type == EngineType.VINA_GPU:
            return VinaGPUEngine(self.docker_manager), "Vina GPU"
        
        elif engine_type == EngineType.GNINA:
            return GninaEngine(self.docker_manager, force_cpu=force_cpu), "GNINA"
        
        elif engine_type == EngineType.GNINA_CPU:
            return GninaEngine(self.docker_manager, force_cpu=True), "GNINA (CPU)"
        
        else:
            logger.warning(f"Unknown engine type: {engine_type}, defaulting to Vina CPU")
            return VinaCPUEngine(self.docker_manager), "Vina CPU (default)"
    
    def get_available_engines(self) -> list:
        """Get list of available engine types"""
        from src.hardware.gpu_detector import has_gpu
        
        engines = [
            {"type": EngineType.VINA_CPU, "name": "Vina CPU", "requires_gpu": False}
        ]
        
        if has_gpu():
            engines.append(
                {"type": EngineType.VINA_GPU, "name": "Vina GPU", "requires_gpu": True}
            )
            engines.append(
                {"type": EngineType.GNINA, "name": "GNINA (GPU)", "requires_gpu": True}
            )
        
        engines.append(
            {"type": EngineType.GNINA_CPU, "name": "GNINA CPU", "requires_gpu": False}
        )
        
        return engines


def get_engine(
    docker_manager=None,
    engine_type: Optional[str] = None,
    prefer_cnn: bool = False,
    force_cpu: bool = False
) -> Tuple[DockingEngine, str]:
    """
    Convenience function to get docking engine.
    
    Args:
        docker_manager: Docker manager instance
        engine_type: Specific engine or None for auto
        prefer_cnn: Prefer CNN scoring when GPU available
        force_cpu: Force CPU mode
        
    Returns:
        tuple: (engine_instance, description)
    """
    factory = DockingEngineFactory(docker_manager)
    return factory.create_engine(engine_type, prefer_cnn, force_cpu)
