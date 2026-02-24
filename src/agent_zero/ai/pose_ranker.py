"""
AI Pose Ranker
Ranks poses using weighted ensemble of multiple scoring methods.
"""

import logging
from typing import Dict, Any, List, Optional
import numpy as np

logger = logging.getLogger(__name__)


class PoseRanker:
    """
    AI-powered pose ranking using weighted ensemble.
    
    Final Score = 0.4 × Docking + 0.3 × CNN + 0.3 × RF
    
    This provides better correlation with experimental binding affinity
    than any single scoring method alone.
    """
    
    DEFAULT_WEIGHTS = {
        "docking": 0.4,
        "cnn": 0.3,
        "rf": 0.3
    }
    
    def __init__(self, weights: Optional[Dict[str, float]] = None):
        self.weights = weights or self.DEFAULT_WEIGHTS.copy()
    
    def rank_poses(
        self,
        poses: List[Dict[str, Any]],
        weights: Optional[Dict[str, float]] = None
    ) -> List[Dict[str, Any]]:
        """
        Rank poses using weighted ensemble.
        
        Args:
            poses: List of pose dictionaries
            weights: Optional custom weights
            
        Returns:
            List: Ranked poses with AI scores
        """
        if not poses:
            return []
        
        w = weights or self.weights
        
        scored_poses = []
        
        for pose in poses:
            ai_score = self._calculate_ai_score(pose, w)
            
            scored_pose = pose.copy()
            scored_pose["ai_rank_score"] = ai_score
            scored_poses.append(scored_pose)
        
        scored_poses.sort(key=lambda x: x.get("ai_rank_score", 0), reverse=True)
        
        for idx, pose in enumerate(scored_poses, 1):
            pose["ai_rank"] = idx
        
        return scored_poses
    
    def _calculate_ai_score(
        self,
        pose: Dict[str, Any],
        weights: Dict[str, float]
    ) -> float:
        """
        Calculate AI ensemble score.
        
        All components normalized to 0-10 scale.
        """
        docking_score = self._normalize_docking(pose.get("binding_energy"))
        cnn_score = self._normalize_cnn(pose.get("gnina_cnn_score"))
        rf_score = self._normalize_rf(pose.get("rf_predicted_pKd"))
        
        ai_score = (
            weights["docking"] * docking_score +
            weights["cnn"] * cnn_score +
            weights["rf"] * rf_score
        )
        
        return round(ai_score, 2)
    
    def _normalize_docking(self, binding_energy: Optional[float]) -> float:
        """Normalize Vina binding energy to 0-10"""
        if binding_energy is None:
            return 5.0
        
        normalized = ((binding_energy + 12) / 12) * 10
        return max(0.0, min(10.0, normalized))
    
    def _normalize_cnn(self, cnn_score: Optional[float]) -> float:
        """Normalize GNINA CNN score to 0-10"""
        if cnn_score is None:
            return 5.0
        
        return max(0.0, min(10.0, cnn_score * 10))
    
    def _normalize_rf(self, rf_pKd: Optional[float]) -> float:
        """Normalize RF pKd prediction to 0-10"""
        if rf_pKd is None:
            return 5.0
        
        return max(0.0, min(10.0, rf_pKd))
    
    def get_top_poses(
        self,
        poses: List[Dict[str, Any]],
        n: int = 5
    ) -> List[Dict[str, Any]]:
        """Get top N poses by AI ranking"""
        ranked = self.rank_poses(poses)
        return ranked[:n]
    
    def compare_poses(
        self,
        pose_a: Dict[str, Any],
        pose_b: Dict[str, Any]
    ) -> int:
        """
        Compare two poses.
        
        Returns:
            1 if pose_a is better
            -1 if pose_b is better
            0 if equal
        """
        score_a = self._calculate_ai_score(pose_a, self.weights)
        score_b = self._calculate_ai_score(pose_b, self.weights)
        
        if score_a > score_b:
            return 1
        elif score_b > score_a:
            return -1
        return 0
