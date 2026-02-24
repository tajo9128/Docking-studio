"""
RF Model Service
Loads and predicts binding affinity using Random Forest models.
"""

import logging
import os
import joblib
import numpy as np
from typing import Optional, Dict, Any

logger = logging.getLogger(__name__)

oddt_available = False

try:
    import oddt
    from oddt.scoring.functions import RFScore
    oddt_available = True
except ImportError:
    oddt_available = False


class RFModelService:
    """
    Manages Random Forest models for binding affinity prediction.
    Supports custom sklearn models and ODDT RFScore.
    """
    
    _model_cache = None
    _model_type = None
    
    def __init__(self, model_dir: str = "models"):
        self.model_dir = model_dir
        os.makedirs(self.model_dir, exist_ok=True)
    
    @classmethod
    def get_model(cls):
        """Get or load the RF model"""
        if cls._model_cache is None:
            cls._load_model()
        return cls._model_cache
    
    @classmethod
    def _load_model(cls):
        """Load model from disk or use ODDT default"""
        candidates = [
            os.path.join("models", "rf_model.pkl"),
            os.path.join("models", "rf_v3.0_production.pkl"),
            os.path.join("models", "rf_v2.1.pkl"),
        ]
        
        for path in candidates:
            if os.path.exists(path):
                try:
                    logger.info(f"Loading custom RF model from {path}")
                    cls._model_cache = joblib.load(path)
                    cls._model_type = 'sklearn'
                    return
                except Exception as e:
                    logger.warning(f"Failed to load {path}: {e}")
        
        if oddt_available:
            try:
                logger.info("Using ODDT RFScore v4")
                cls._model_cache = RFScore.rfscore(version=4)
                cls._model_type = 'oddt'
            except Exception as e:
                logger.warning(f"ODDT RFScore failed: {e}")
        
        if cls._model_cache is None:
            logger.warning("No RF model available - using fallback")
    
    def predict_pkd(self, features: np.ndarray) -> float:
        """
        Predict binding affinity (pKd) from feature vector.
        
        Args:
            features: Feature vector (54 dims) or (n_samples, 54)
            
        Returns:
            float: Predicted pKd
        """
        model = self.get_model()
        
        if model is None:
            return 6.0
        
        try:
            if self._model_type == 'sklearn':
                if features.ndim == 1:
                    features = features.reshape(1, -1)
                return float(model.predict(features)[0])
            elif self._model_type == 'oddt':
                return float(model.predict(features)[0])
            else:
                return 6.0
        except Exception as e:
            logger.error(f"RF prediction failed: {e}")
            return 6.0
    
    def predict_batch(self, features_list: list) -> list:
        """Predict for multiple samples"""
        return [self.predict_pkd(f) for f in features_list]
    
    @staticmethod
    def load_custom_model(path: str) -> bool:
        """Load a custom model from path"""
        try:
            RFModelService._model_cache = joblib.load(path)
            RFModelService._model_type = 'sklearn'
            logger.info(f"Loaded custom model from {path}")
            return True
        except Exception as e:
            logger.error(f"Failed to load model: {e}")
            return False
