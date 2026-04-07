"""
Meta-Parameter Self-Learning System
Learns optimal docking/MD/QSAR parameters per protein family from historical results.
"""

import json
import logging
import hashlib
from pathlib import Path
from typing import Dict, Any, List, Optional
from datetime import datetime

logger = logging.getLogger(__name__)


class MetaParameterLearner:
    """
    Learns and suggests optimal parameters for docking/MD/QSAR
    based on historical experiment outcomes per target family.
    """

    # Default parameters per service
    DEFAULTS = {
        "docking": {
            "exhaustiveness": 8,
            "box_size": [10.0, 10.0, 10.0],
            "num_modes": 10,
            "scoring": "vina",
        },
        "md": {
            "steps": 50000,
            "temperature": 300.0,
            "frame_interval": 500,
            "solvent_model": "tip3p",
        },
        "qsar": {
            "model_type": "RandomForest",
            "cv_folds": 5,
            "n_estimators": 100,
        },
    }

    # Target family classification by PDB keywords
    # Order matters: more specific families must come before general ones.
    # nuclear_receptor before gpcr to prevent "estrogen receptor" matching gpcr via "receptor".
    FAMILY_KEYWORDS = {
        "kinase": ["kinase", "atp", "tyrosine", "serine", "threonine"],
        "nuclear_receptor": ["nuclear receptor", "estrogen", "androgen", "ppar", "fxr", "rar", "rxr"],
        "gpcr": ["gpcr", "adrenergic", "dopamine", "serotonin", "rhodopsin", "muscarinic"],
        "protease": ["protease", "hiv", "thrombin", "caspase", "peptidase"],
        "ion_channel": ["channel", "ion", "potassium", "sodium", "calcium"],
        "enzyme": ["enzyme", "dehydrogenase", "oxidase", "transferase", "hydrolase"],
    }

    def __init__(self, persist_dir: str = "./data/meta_params"):
        self.persist_dir = Path(persist_dir)
        self.persist_dir.mkdir(parents=True, exist_ok=True)
        self._records: Dict[str, List[Dict[str, Any]]] = {}
        self._load()

    def _load(self):
        data_file = self.persist_dir / "meta_params.json"
        if data_file.exists():
            try:
                with open(data_file, 'r') as f:
                    self._records = json.load(f)
            except Exception as e:
                logger.warning(f"Failed to load meta-params: {e}")

    def _save(self):
        data_file = self.persist_dir / "meta_params.json"
        try:
            with open(data_file, 'w') as f:
                json.dump(self._records, f, indent=2, default=str)
        except Exception as e:
            logger.warning(f"Failed to save meta-params: {e}")

    def classify_target(self, target_info: str) -> str:
        """Classify target into family based on name/description."""
        info_lower = target_info.lower()
        for family, keywords in self.FAMILY_KEYWORDS.items():
            if any(kw in info_lower for kw in keywords):
                return family
        return "unknown"

    def suggest_params(self, target_info: str, service: str = "docking", ligand_size: int = 0) -> Dict[str, Any]:
        """
        Suggest optimal parameters based on historical success for similar targets.
        """
        family = self.classify_target(target_info)
        family_records = self._records.get(family, [])
        service_records = [r for r in family_records if r.get("service") == service]

        if not service_records:
            return self.DEFAULTS.get(service, {}).copy()

        # Weight by success rate and score
        successful = [r for r in service_records if r.get("success", False)]
        if successful:
            # Find params with highest average score among successful runs
            param_scores: Dict[str, List[float]] = {}
            for r in successful:
                params_key = json.dumps(r.get("params", {}), sort_keys=True)
                score = r.get("score", 0)
                if params_key not in param_scores:
                    param_scores[params_key] = []
                param_scores[params_key].append(score)

            best_key = max(param_scores.keys(), key=lambda k: sum(param_scores[k]) / len(param_scores[k]))
            best_params = json.loads(best_key)
            return best_params

        # Fallback: use most recent params
        return service_records[-1].get("params", self.DEFAULTS.get(service, {}).copy())

    def record_outcome(self, target_info: str, service: str, params: Dict[str, Any],
                       success: bool, score: float, error: str = None):
        """Record the outcome of a simulation for future learning."""
        family = self.classify_target(target_info)

        if family not in self._records:
            self._records[family] = []

        record = {
            "service": service,
            "params": params,
            "success": success,
            "score": score,
            "error": error,
            "timestamp": datetime.now().isoformat(),
            "target": target_info,
        }

        self._records[family].append(record)

        # Keep only last 100 records per family to prevent bloat
        if len(self._records[family]) > 100:
            self._records[family] = self._records[family][-100:]

        self._save()

    def get_family_stats(self, family: str = None) -> Dict[str, Any]:
        """Get statistics for a target family."""
        if family:
            families = {family: self._records.get(family, [])}
        else:
            families = self._records

        stats = {}
        for fam, records in families.items():
            if not records:
                continue
            total = len(records)
            successful = sum(1 for r in records if r.get("success"))
            avg_score = sum(r.get("score", 0) for r in records if r.get("success")) / max(successful, 1)
            services = {}
            for r in records:
                svc = r.get("service", "unknown")
                if svc not in services:
                    services[svc] = {"total": 0, "success": 0}
                services[svc]["total"] += 1
                if r.get("success"):
                    services[svc]["success"] += 1

            stats[fam] = {
                "total_experiments": total,
                "success_rate": round(successful / max(total, 1), 3),
                "avg_score": round(avg_score, 3),
                "services": services,
            }

        return stats

    def get_param_history(self, target_info: str, service: str = "docking", n: int = 10) -> List[Dict]:
        """Get recent parameter history for a target."""
        family = self.classify_target(target_info)
        records = [r for r in self._records.get(family, []) if r.get("service") == service]
        return records[-n:]


meta_learner = MetaParameterLearner()
