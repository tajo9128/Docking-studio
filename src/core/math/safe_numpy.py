"""Safe Numpy Operations - Unit-safe Calculations"""

import numpy as np
from typing import List, Optional, Tuple
from ...logging_config import get_logger

logger = get_logger(__name__)

class SafeNumpy:
    """Safe numpy operations to prevent numerical errors"""
    
    @staticmethod
    def safe_sqrt(value: float) -> float:
        """Safe square root"""
        if np.isnan(value) or value < 0:
            return 0.0
        return np.sqrt(value)
    
    @staticmethod
    def safe_arccos(cos_angle: float) -> float:
        """Safe arccos with clamping"""
        clamped = max(-1.0, min(1.0, cos_angle))
        return np.arccos(clamped)
    
    @staticmethod
    def safe_mean(values: List[float]) -> float:
        """Safe mean calculation"""
        if not values:
            return 0.0
        
        clean_values = np.array([v for v in values if not (np.isnan(v) or np.isinf(v))])
        
        if clean_values.size == 0:
            return 0.0
        
        return np.mean(clean_values)
    
    @staticmethod
    def safe_std(values: List[float]) -> float:
        """Safe standard deviation calculation"""
        if not values:
            return 0.0
        
        clean_values = np.array([v for v in values if not (np.isnan(v) or np.isinf(v))])
        
        if clean_values.size < 2:
            return 0.0
        
        return np.std(clean_values)
    
    @staticmethod
    def safe_normalize(vector: List[float]) -> List[float]:
        """Safe vector normalization"""
        norm = np.linalg.norm(vector)
        
        if norm == 0 or not np.isfinite(norm):
            return [0.0] * len(vector)
        
        return (vector / norm).tolist()
    
    @staticmethod
    def safe_distance(a: List[float], b: List[float]) -> float:
        """Safe Euclidean distance calculation"""
        if len(a) != 3 or len(b) != 3:
            return 0.0
        
        dx = a[0] - b[0]
        dy = a[1] - b[1]
        dz = a[2] - b[2]
        
        return SafeNumpy.safe_sqrt(dx**2 + dy**2 + dz**2)
    
    @staticmethod
    def safe_bounding_box(coords: List[List[float]]) -> Tuple[List[float], List[float]]:
        """Calculate safe bounding box"""
        if not coords:
            return [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]
        
        coords_arr = np.array(coords)
        
        if np.any(np.isnan(coords_arr)) or np.any(np.isinf(coords_arr)):
            return [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]
        
        min_coords = np.min(coords_arr, axis=0)
        max_coords = np.max(coords_arr, axis=0)
        
        return min_coords.tolist(), max_coords.tolist()
    
    @staticmethod
    def safe_center_of_mass(coords: List[List[float]]) -> List[float]:
        """Calculate safe center of mass (geometric center)"""
        if not coords:
            return [0.0, 0.0, 0.0]
        
        coords_arr = np.array(coords)
        
        if np.any(np.isnan(coords_arr)) or np.any(np.isinf(coords_arr)):
            return [0.0, 0.0, 0.0]
        
        center = np.mean(coords_arr, axis=0)
        
        if not np.all(np.isfinite(center)):
            return [0.0, 0.0, 0.0]
        
        return center.tolist()
