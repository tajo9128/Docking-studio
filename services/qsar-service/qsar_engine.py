"""
QSAR engine - model training, evaluation, and prediction
"""

import logging
import json
import uuid
import asyncio
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from sklearn.model_selection import cross_val_predict, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.cross_decomposition import PLSRegression
from sklearn.svm import SVR
from sklearn.linear_model import Ridge, Lasso
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
import plotly.graph_objects as go
import plotly.io as pio

from .descriptors import descriptors_for_smiles_list, calculate_descriptors
from .model_store import save_model, load_model
from .ad_checker import assess_applicability_domain

logger = logging.getLogger("qsar-engine")

MODEL_REGISTRY = {
    "RandomForest": RandomForestRegressor,
    "GradientBoosting": GradientBoostingRegressor,
    "PLS": PLSRegression,
    "SVR": SVR,
    "Ridge": Ridge,
    "Lasso": Lasso,
}

DEFAULT_PARAMS = {
    "RandomForest": {"n_estimators": 100, "max_depth": 10, "random_state": 42},
    "GradientBoosting": {"n_estimators": 100, "max_depth": 5, "random_state": 42},
    "SVR": {"C": 1.0, "kernel": "rbf"},
    "Ridge": {"alpha": 1.0},
    "Lasso": {"alpha": 1.0},
    "PLS": {"n_components": 2},
}


def _get_model(model_type: str, params: Optional[Dict[str, Any]] = None):
    if model_type not in MODEL_REGISTRY:
        raise ValueError(
            f"Unknown model type: {model_type}. Available: {list(MODEL_REGISTRY.keys())}"
        )

    model_cls = MODEL_REGISTRY[model_type]
    params = params or DEFAULT_PARAMS.get(model_type, {})

    if model_type == "PLS":
        n_components = params.get("n_components", 2)
        return model_cls(n_components=n_components)
    elif model_type in ("SVR",):
        return model_cls(**params)
    elif model_type in ("Ridge", "Lasso"):
        return model_cls(**params)
    else:
        return model_cls(**params)


def train_model(
    X: List[List[float]],
    y: List[float],
    feature_names: List[str],
    model_type: str,
    model_name: str,
    activity_column: str,
    descriptor_groups: List[str],
    cv_folds: int = 5,
    model_params: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Train a QSAR model with cross-validation.

    Returns:
        Dictionary with metrics and Plotly figure data
    """
    X_arr = np.array(X, dtype=float)
    y_arr = np.array(y, dtype=float)

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_arr)

    model = _get_model(model_type, model_params)

    kfold = KFold(n_splits=cv_folds, shuffle=True, random_state=42)

    try:
        if model_type == "PLS":
            y_pred_cv = cross_val_predict(model, X_scaled, y_arr, cv=kfold)
            model.fit(X_scaled, y_arr)
            y_pred_train = model.predict(X_scaled)
        elif model_type in ("SVR",):
            y_pred_cv = cross_val_predict(model, X_scaled, y_arr, cv=kfold)
            model.fit(X_scaled, y_arr)
            y_pred_train = model.predict(X_scaled)
        else:
            y_pred_cv = cross_val_predict(model, X_scaled, y_arr, cv=kfold)
            model.fit(X_scaled, y_arr)
            y_pred_train = model.predict(X_scaled)
    except Exception as e:
        raise RuntimeError(f"Training failed: {e}")

    cv_r2 = r2_score(y_arr, y_pred_cv)
    cv_rmse = np.sqrt(mean_squared_error(y_arr, y_pred_cv))
    cv_mae = mean_absolute_error(y_arr, y_pred_cv)
    train_r2 = r2_score(y_arr, y_pred_train)

    cv_scores = []
    for train_idx, val_idx in kfold.split(X_scaled):
        X_tr, X_val = X_scaled[train_idx], X_scaled[val_idx]
        y_tr, y_val = y_arr[train_idx], y_arr[val_idx]
        m = _get_model(model_type, model_params)
        m.fit(X_tr, y_tr)
        cv_scores.append(r2_score(y_val, m.predict(X_val)))

    cv_std = float(np.std(cv_scores))

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=y_arr,
            y=y_pred_cv,
            mode="markers",
            marker=dict(color="steelblue", size=8, opacity=0.7),
            name="Predicted vs Actual",
        )
    )
    min_val = float(min(y_arr.min(), y_pred_cv.min()))
    max_val = float(max(y_arr.max(), y_pred_cv.max()))
    fig.add_trace(
        go.Scatter(
            x=[min_val, max_val],
            y=[min_val, max_val],
            mode="lines",
            line=dict(color="red", dash="dash"),
            name="Perfect Prediction",
        )
    )
    fig.update_layout(
        title=dict(text=f"Actual vs Predicted ({model_type})", font=dict(size=14)),
        xaxis_title="Actual Activity",
        yaxis_title="Predicted Activity",
        width=500,
        height=400,
        margin=dict(l=60, r=20, t=40, b=60),
    )

    scatter_plot = pio.to_json(fig)

    metadata = save_model(
        model=model,
        scaler=scaler,
        feature_names=feature_names,
        model_type=model_type,
        model_name=model_name,
        metrics={
            "cv_r2": float(cv_r2),
            "cv_rmse": float(cv_rmse),
            "cv_mae": float(cv_mae),
            "train_r2": float(train_r2),
            "cv_std": cv_std,
            "cv_scores": [float(s) for s in cv_scores],
        },
        activity_column=activity_column,
        descriptor_groups=descriptor_groups,
        X_train=X_arr,
    )

    return {
        "model_id": metadata["model_id"],
        "model_name": model_name,
        "model_type": model_type,
        "metrics": {
            "cv_r2": float(cv_r2),
            "cv_rmse": float(cv_rmse),
            "cv_mae": float(cv_mae),
            "train_r2": float(train_r2),
            "cv_std": cv_std,
        },
        "cv_scores": [float(s) for s in cv_scores],
        "cv_folds": cv_folds,
        "n_compounds": len(y_arr),
        "n_features": len(feature_names),
        "scatter_plot": scatter_plot,
        "created_at": metadata["created_at"],
    }


def predict_single(
    model_id: str,
    smiles: str,
    X_train: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Predict activity for a single SMILES using a saved model.
    """
    result = load_model(model_id)
    if result is None:
        raise ValueError(f"Model {model_id} not found")

    model_obj = result["model_obj"]
    model = model_obj["model"]
    scaler = model_obj["scaler"]
    feature_names = model_obj["feature_names"]
    X_train_ref = model_obj.get("X_train")

    mol_data, _, failed = descriptors_for_smiles_list(
        [smiles], include_fingerprints=False
    )
    if failed:
        raise ValueError(f"Invalid SMILES: {smiles}")

    x_dict = mol_data[0]
    x_vec = np.array([[x_dict.get(fn, 0.0) for fn in feature_names]], dtype=float)

    x_scaled = scaler.transform(x_vec)
    predicted = float(model.predict(x_scaled)[0])

    ad_result = None
    if X_train_ref is not None:
        ad_result = assess_applicability_domain(X_train_ref, x_vec, method="leverage")

    return {
        "smiles": smiles,
        "predicted_activity": predicted,
        "ad_status": ad_result["statuses"][0] if ad_result else "unknown",
        "ad_leverage": ad_result["leverages"][0]
        if ad_result and ad_result.get("leverages")
        else None,
        "ad_warning_threshold": ad_result.get("warning_threshold"),
        "ad_danger_threshold": ad_result.get("danger_threshold"),
    }


def predict_batch(
    model_id: str,
    smiles_list: List[str],
) -> Dict[str, Any]:
    """
    Predict activity for multiple SMILES.
    """
    result = load_model(model_id)
    if result is None:
        raise ValueError(f"Model {model_id} not found")

    model_obj = result["model_obj"]
    model = model_obj["model"]
    scaler = model_obj["scaler"]
    feature_names = model_obj["feature_names"]
    X_train_ref = model_obj.get("X_train")

    predictions = []
    failed = []

    for smiles in smiles_list:
        mol_data, valid, fail = descriptors_for_smiles_list([smiles])
        if fail:
            failed.append(smiles)
            predictions.append(
                {
                    "smiles": smiles,
                    "predicted_activity": None,
                    "error": "Invalid SMILES",
                }
            )
            continue

        x_dict = mol_data[0]
        x_vec = np.array([[x_dict.get(fn, 0.0) for fn in feature_names]], dtype=float)
        x_scaled = scaler.transform(x_vec)
        pred = float(model.predict(x_scaled)[0])

        ad_result = None
        if X_train_ref is not None:
            ad_result = assess_applicability_domain(
                X_train_ref, x_vec, method="leverage"
            )

        predictions.append(
            {
                "smiles": smiles,
                "predicted_activity": pred,
                "ad_status": ad_result["statuses"][0] if ad_result else "unknown",
                "ad_leverage": ad_result["leverages"][0]
                if ad_result and ad_result.get("leverages")
                else None,
            }
        )

    valid_preds = [p for p in predictions if p.get("predicted_activity") is not None]
    n_in_domain = sum(1 for p in valid_preds if p.get("ad_status") == "in_domain")
    n_warning = sum(1 for p in valid_preds if p.get("ad_status") == "warning")
    n_out = sum(1 for p in valid_preds if p.get("ad_status") == "out_of_domain")

    return {
        "predictions": predictions,
        "n_total": len(predictions),
        "n_failed": len(failed),
        "n_in_domain": n_in_domain,
        "n_warning": n_warning,
        "n_out_of_domain": n_out,
    }
