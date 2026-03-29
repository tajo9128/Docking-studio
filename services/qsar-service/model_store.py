"""
Model persistence using joblib
"""

import os
import json
import uuid
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional
from datetime import datetime
import joblib
import numpy as np

logger = logging.getLogger("qsar-model-store")

MODELS_DIR = Path("/app/storage/qsar/models")
MODELS_DIR.mkdir(parents=True, exist_ok=True)


def _model_path(model_id: str) -> Path:
    return MODELS_DIR / f"{model_id}.pkl"


def _metadata_path(model_id: str) -> Path:
    return MODELS_DIR / f"{model_id}_meta.json"


def save_model(
    model: Any,
    scaler: Any,
    feature_names: List[str],
    model_type: str,
    model_name: str,
    metrics: Dict[str, float],
    activity_column: str,
    descriptor_groups: List[str],
    X_train: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Save a trained QSAR model to disk.

    Args:
        model: Trained sklearn model
        scaler: Fitted StandardScaler
        feature_names: List of descriptor names
        model_type: Algorithm name (RandomForest, PLS, etc.)
        model_name: Human-readable name
        metrics: Dict of evaluation metrics (r2, rmse, mae, cv_std)
        activity_column: Name of the activity column
        descriptor_groups: Descriptor groups used
        X_train: Training descriptor matrix (for AD)

    Returns:
        Model metadata
    """
    model_id = str(uuid.uuid4())[:8]

    model_obj = {
        "model": model,
        "scaler": scaler,
        "feature_names": feature_names,
        "X_train": X_train,
        "activity_column": activity_column,
        "descriptor_groups": descriptor_groups,
    }

    joblib.dump(model_obj, _model_path(model_id))

    metadata = {
        "model_id": model_id,
        "name": model_name,
        "model_type": model_type,
        "feature_names": feature_names,
        "n_features": len(feature_names),
        "metrics": metrics,
        "activity_column": activity_column,
        "descriptor_groups": descriptor_groups,
        "created_at": datetime.now().isoformat(),
    }

    with open(_metadata_path(model_id), "w") as f:
        json.dump(metadata, f, indent=2)

    logger.info(f"Saved model {model_id}: {model_name} ({model_type})")
    return metadata


def load_model(model_id: str) -> Optional[Dict[str, Any]]:
    """Load a model from disk"""
    try:
        model_obj = joblib.load(_model_path(model_id))
        with open(_metadata_path(model_id)) as f:
            metadata = json.load(f)
        return {"model_obj": model_obj, "metadata": metadata}
    except FileNotFoundError:
        logger.warning(f"Model {model_id} not found")
        return None
    except Exception as e:
        logger.error(f"Failed to load model {model_id}: {e}")
        return None


def list_models() -> List[Dict[str, Any]]:
    """List all saved models"""
    models = []
    for pkl_path in MODELS_DIR.glob("*_meta.json"):
        try:
            with open(pkl_path) as f:
                metadata = json.load(f)
            models.append(metadata)
        except Exception as e:
            logger.warning(f"Failed to read metadata {pkl_path}: {e}")
    return sorted(models, key=lambda m: m.get("created_at", ""), reverse=True)


def delete_model(model_id: str) -> bool:
    """Delete a model from disk"""
    try:
        pkl = _model_path(model_id)
        meta = _metadata_path(model_id)
        if pkl.exists():
            pkl.unlink()
        if meta.exists():
            meta.unlink()
        logger.info(f"Deleted model {model_id}")
        return True
    except Exception as e:
        logger.error(f"Failed to delete model {model_id}: {e}")
        return False


def get_model(model_id: str) -> Optional[Dict[str, Any]]:
    """Get model metadata"""
    try:
        with open(_metadata_path(model_id)) as f:
            return json.load(f)
    except FileNotFoundError:
        return None
    except Exception as e:
        logger.error(f"Failed to read model {model_id}: {e}")
        return None
