"""
Applicability Domain assessment using leverage method
"""

import numpy as np
from typing import Dict, Any, List, Tuple, Optional


def calculate_leverage(X_train: np.ndarray, X_test: np.ndarray) -> np.ndarray:
    """
    Calculate leverage (Mahalanobis distance) for each test compound.

    Leverage = diag(X_test @ H @ X_test.T) where H = X(X'X)^-1 X'

    Args:
        X_train: Training descriptor matrix (n_train, n_features)
        X_test: Test descriptor matrix (n_test, n_features)

    Returns:
        Leverage values for each test compound
    """
    n_train, n_features = X_train.shape
    k = n_features

    try:
        XTX_inv = np.linalg.inv(X_train.T @ X_train + 1e-10 * np.eye(n_features))
    except np.linalg.LinAlgError:
        XTX_inv = np.linalg.pinv(X_train.T @ X_train)

    H = X_train @ XTX_inv @ X_train.T
    leverage_values = np.diag(X_test @ H @ X_test.T)

    return leverage_values


def get_leverage_threshold(
    n_train: int, n_features: int, warning: float = 3.0, danger: float = 2.0
) -> Tuple[float, float]:
    """
    Get leverage warning and danger thresholds.

    Warning zone: h > 3(k+1)/n
    Out of domain: h > 2*3(k+1)/n

    Args:
        n_train: Number of training compounds
        n_features: Number of descriptors
        warning: Multiplier for warning threshold
        danger: Multiplier for danger threshold

    Returns:
        (warning_threshold, danger_threshold)
    """
    if n_train <= n_features + 1:
        return 0.5, 0.9

    warning_threshold = warning * (n_features + 1) / n_train
    danger_threshold = danger * (n_features + 1) / n_train

    return float(warning_threshold), float(danger_threshold)


def assess_applicability_domain(
    X_train: np.ndarray,
    X_test: np.ndarray,
    y_train: Optional[np.ndarray] = None,
    method: str = "leverage",
) -> Dict[str, Any]:
    """
    Assess applicability domain for test compounds.

    Args:
        X_train: Training descriptor matrix
        X_test: Test descriptor matrix
        y_train: Optional training activities (for descriptor range method)
        method: "leverage" or "descriptor_range"

    Returns:
        Dictionary with AD status per compound and summary statistics
    """
    if method == "leverage":
        leverage = calculate_leverage(X_train, X_test)
        n_train, n_features = X_train.shape
        warn_thresh, danger_thresh = get_leverage_threshold(n_train, n_features)

        statuses = []
        for lev in leverage:
            if lev > danger_thresh:
                statuses.append("out_of_domain")
            elif lev > warn_thresh:
                statuses.append("warning")
            else:
                statuses.append("in_domain")

        return {
            "method": "leverage",
            "leverages": leverage.tolist(),
            "warning_threshold": warn_thresh,
            "danger_threshold": danger_thresh,
            "statuses": statuses,
            "n_in_domain": sum(1 for s in statuses if s == "in_domain"),
            "n_warning": sum(1 for s in statuses if s == "warning"),
            "n_out_of_domain": sum(1 for s in statuses if s == "out_of_domain"),
        }

    elif method == "descriptor_range":
        n_train, n_features = X_train.shape
        train_mins = X_train.min(axis=0)
        train_maxs = X_train.max(axis=0)
        train_ranges = train_maxs - train_mins
        margin = 0.05 * train_ranges

        in_domain_mask = np.all(
            (X_test >= (train_mins - margin)) & (X_test <= (train_maxs + margin)),
            axis=1,
        )

        statuses = ["in_domain" if m else "out_of_domain" for m in in_domain_mask]

        return {
            "method": "descriptor_range",
            "statuses": statuses,
            "n_in_domain": int(np.sum(in_domain_mask)),
            "n_out_of_domain": int(np.sum(~in_domain_mask)),
        }

    else:
        raise ValueError(f"Unknown AD method: {method}")


def check_single_compound_ad(
    X_train: np.ndarray,
    x_test: np.ndarray,
    method: str = "leverage",
) -> Dict[str, Any]:
    """Check AD status for a single compound"""
    result = assess_applicability_domain(X_train, x_test.reshape(1, -1), method=method)
    status = result["statuses"][0] if result["statuses"] else "unknown"
    leverage = result["leverages"][0] if result.get("leverages") else None

    return {
        "status": status,
        "leverage": leverage,
        "warning_threshold": result.get("warning_threshold"),
        "danger_threshold": result.get("danger_threshold"),
    }
