"""
Consensus Scorer
Computes weighted consensus scores from Vina, GNINA, and RF predictions.
"""

import logging
from typing import Dict, Any, List, Optional
import numpy as np

logger = logging.getLogger(__name__)


class ConsensusScorer:
    """
    Computes consensus scores from multiple scoring methods.
    
    Consensus formula:
        consensus = 0.30 * Vina_norm + 0.40 * GNINA_norm + 0.30 * RF_norm
    
    Where each score is normalized to 0-10 scale.
    """
    
    DEFAULT_WEIGHTS = {
        "vina": 0.30,
        "gnina": 0.40,
        "rf": 0.30
    }
    
    def __init__(self, weights: Optional[Dict[str, float]] = None):
        """
        Initialize consensus scorer.
        
        Args:
            weights: Custom weights for each scoring method.
                    Must sum to 1.0. Defaults to 0.30/0.40/0.30
        """
        if weights:
            total = sum(weights.values())
            if abs(total - 1.0) > 0.01:
                raise ValueError("Weights must sum to 1.0")
            self.weights = weights
        else:
            self.weights = self.DEFAULT_WEIGHTS.copy()
    
    def normalize_vina(self, affinity: float) -> float:
        """
        Normalize Vina affinity to 0-10 scale.
        
        Vina typically ranges from -12 to 0 kcal/mol
        
        Args:
            affinity: Vina binding affinity
            
        Returns:
            float: Normalized score (0-10)
        """
        if affinity is None:
            return 5.0
        
        normalized = ((affinity + 12) / 12) * 10
        return max(0.0, min(10.0, normalized))
    
    def normalize_gnina(self, cnn_score: float) -> float:
        """
        Normalize GNINA CNN score to 0-10 scale.
        
        CNN score typically ranges from 0 to 1
        
        Args:
            cnn_score: GNINA CNN score
            
        Returns:
            float: Normalized score (0-10)
        """
        if cnn_score is None:
            return 5.0
        
        return max(0.0, min(10.0, cnn_score * 10))
    
    def normalize_rf(self, pkd: float) -> float:
        """
        Normalize RF pKd prediction.
        
        pKd typically ranges from 0 to 10 (or more)
        
        Args:
            pkd: RF predicted pKd
            
        Returns:
            float: Normalized score (0-10)
        """
        if pkd is None:
            return 5.0
        
        return max(0.0, min(10.0, pkd))
    
    def calculate_consensus(
        self,
        vina_affinity: Optional[float] = None,
        gnina_cnn_score: Optional[float] = None,
        rf_pkd: Optional[float] = None
    ) -> Dict[str, Any]:
        """
        Calculate consensus score from available predictions.
        
        Args:
            vina_affinity: Vina binding affinity (kcal/mol)
            gnina_cnn_score: GNINA CNN score (0-1)
            rf_pkd: RF predicted pKd
            
        Returns:
            dict: Consensus results with individual normalized scores
        """
        vina_norm = self.normalize_vina(vina_affinity)
        gnina_norm = self.normalize_gnina(gnina_cnn_score)
        rf_norm = self.normalize_rf(rf_pkd)
        
        available_weights = []
        weighted_sum = 0.0
        total_weight = 0.0
        
        if vina_affinity is not None:
            weighted_sum += vina_norm * self.weights["vina"]
            total_weight += self.weights["vina"]
            available_weights.append("vina")
        
        if gnina_cnn_score is not None:
            weighted_sum += gnina_norm * self.weights["gnina"]
            total_weight += self.weights["gnina"]
            available_weights.append("gnina")
        
        if rf_pkd is not None:
            weighted_sum += rf_norm * self.weights["rf"]
            total_weight += self.weights["rf"]
            available_weights.append("rf")
        
        if total_weight > 0:
            consensus = (weighted_sum / total_weight) * 10
        else:
            consensus = 5.0
        
        return {
            "consensus_score": round(consensus, 2),
            "vina_normalized": round(vina_norm, 2),
            "gnina_normalized": round(gnina_norm, 2),
            "rf_normalized": round(rf_norm, 2),
            "vina_raw": vina_affinity,
            "gnina_raw": gnina_cnn_score,
            "rf_raw": rf_pkd,
            "methods_used": available_weights
        }
    
    def score_poses(
        self,
        poses: List[Dict[str, Any]],
        use_rf: bool = True,
        receptor_path: Optional[str] = None,
        ligand_path: Optional[str] = None
    ) -> List[Dict[str, Any]]:
        """
        Score multiple poses with consensus.
        
        Args:
            poses: List of pose dictionaries
            use_rf: Whether to use RF scoring
            receptor_path: Path to receptor for RF features
            ligand_path: Path to ligand for RF features
            
        Returns:
            list: Poses with consensus scores
        """
        from .rf_model_service import RFModelService
        from ..feature_extractor import FeatureExtractor
        
        rf_service = RFModelService()
        feature_extractor = FeatureExtractor()
        
        scored_poses = []
        
        for pose in poses:
            vina_affinity = pose.get("binding_energy")
            gnina_cnn = pose.get("gnina_cnn_score")
            
            rf_pkd = None
            if use_rf and receptor_path and ligand_path:
                try:
                    features = feature_extractor.build_feature_vector(
                        vina_affinity=vina_affinity or -7.0,
                        gnina_cnn_score=gnina_cnn,
                        gnina_cnn_affinity=pose.get("gnina_cnn_affinity"),
                        receptor_path=receptor_path,
                        ligand_path=ligand_path
                    )
                    rf_pkd = rf_service.predict_pkd(features)
                except Exception as e:
                    logger.warning(f"RF prediction failed: {e}")
            
            consensus = self.calculate_consensus(
                vina_affinity=vina_affinity,
                gnina_cnn_score=gnina_cnn,
                rf_pkd=rf_pkd
            )
            
            scored_pose = pose.copy()
            scored_pose["rf_predicted_pkd"] = rf_pkd
            scored_pose["consensus_score"] = consensus["consensus_score"]
            scored_pose["consensus_details"] = consensus
            
            scored_poses.append(scored_pose)
        
        scored_poses.sort(key=lambda x: x.get("consensus_score", 0), reverse=True)
        
        return scored_poses
