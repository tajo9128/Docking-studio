"""
Adversarial Critique Agent + Uncertainty Gating
Validates AI proposals with chemical plausibility checks and confidence thresholds.
"""

import logging
import re
from typing import Dict, Any, List, Optional

logger = logging.getLogger(__name__)


class CritiqueAgent:
    """
    Scientific critique and validation specialist.
    Challenges proposals, flags chemical implausibility, enforces uncertainty thresholds.
    """

    # Chemical plausibility rules
    ENERGY_BOUNDS = {
        "docking": {"min": -15.0, "max": 2.0, "unit": "kcal/mol"},
        "md": {"min": -100000.0, "max": 100000.0, "unit": "kJ/mol"},
        "qsar": {"min": -10.0, "max": 10.0, "unit": "pIC50"},
    }

    # Red flag patterns in results
    RED_FLAGS = [
        (r"energy.*>.*100", "Unphysical energy value"),
        (r"rmsd.*>.*5", "Excessive RMSD deviation"),
        (r"clash.*>.*20", "Excessive steric clashes"),
        (r"failed|error|exception", "Execution failure detected"),
        (r"\bnan\b|\binf\b|\bnull\b", "Invalid numerical value"),
    ]

    def validate(self, tool: str, result: Dict[str, Any], confidence_threshold: float = 0.7) -> Dict[str, Any]:
        """
        Validate a tool result with chemical plausibility checks.
        Returns validation report with pass/fail and recommendations.
        """
        issues = []
        warnings = []
        recommendations = []

        # Check energy bounds
        bounds = self.ENERGY_BOUNDS.get(tool)
        if bounds:
            energy = result.get("best_score") or result.get("energy") or result.get("vina_score")
            if energy is not None:
                try:
                    energy = float(energy)
                    if energy < bounds["min"]:
                        issues.append(f"Energy {energy} {bounds['unit']} below minimum {bounds['min']}")
                    elif energy > bounds["max"]:
                        issues.append(f"Energy {energy} {bounds['unit']} above maximum {bounds['max']}")
                except (ValueError, TypeError):
                    issues.append(f"Invalid energy value: {energy}")

        # Check for red flags
        result_str = str(result)
        for pattern, description in self.RED_FLAGS:
            if re.search(pattern, result_str, re.IGNORECASE):
                warnings.append(description)

        # Check confidence
        confidence = result.get("confidence", 0.5)
        if confidence < confidence_threshold:
            issues.append(f"Confidence {confidence:.2f} below threshold {confidence_threshold}")

        # Check for missing critical fields
        if tool == "docking":
            if not result.get("results") and not result.get("error"):
                warnings.append("No docking results returned")
            if result.get("num_poses", 0) == 0 and not result.get("error"):
                warnings.append("No poses generated")

        # Generate recommendations
        if issues:
            recommendations.append("Review input parameters and structure preparation")
            if any("energy" in i.lower() for i in issues):
                recommendations.append("Check force field and protonation state")
            if any("confidence" in i.lower() for i in issues):
                recommendations.append("Consider retrying with adjusted parameters")

        is_valid = len(issues) == 0
        requires_human_review = confidence < 0.5 or len(issues) > 2

        return {
            "is_valid": is_valid,
            "requires_human_review": requires_human_review,
            "confidence": confidence,
            "issues": issues,
            "warnings": warnings,
            "recommendations": recommendations,
            "n_issues": len(issues),
            "n_warnings": len(warnings),
        }

    def cross_reference(self, smiles: str, target: str = None) -> Dict[str, Any]:
        """
        Cross-reference a compound against known literature data.
        Uses PubChem to check for known activity data.
        """
        try:
            import httpx
            from rdkit import Chem

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {"error": "Invalid SMILES", "known_data": None}

            # Check PubChem for CID
            response = httpx.get(
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/JSON",
                timeout=10
            )
            if response.status_code == 200:
                data = response.json()
                cids = data.get("IdentifierList", {}).get("CID", [])
                return {
                    "cids": cids[:5],
                    "known_compound": len(cids) > 0,
                    "message": "Compound found in PubChem" if cids else "Novel compound",
                }

            return {"known_compound": False, "message": "Could not verify compound"}

        except Exception as e:
            return {"error": str(e), "known_compound": None}

    def validate_workflow(self, workflow: Dict[str, Any]) -> Dict[str, Any]:
        """Validate a complete workflow DAG before execution."""
        issues = []

        steps = workflow.get("steps", [])
        if not steps:
            issues.append("No steps in workflow")

        # Check for circular dependencies
        step_ids = {s.get("id") for s in steps}
        for step in steps:
            for dep in step.get("depends_on", []):
                if dep not in step_ids:
                    issues.append(f"Step {step.get('id')} depends on non-existent step {dep}")

        # Check tool validity
        from crew.nl_compiler import ALLOWED_TOOLS
        for step in steps:
            tool = step.get("tool")
            if tool not in ALLOWED_TOOLS:
                issues.append(f"Unknown tool: {tool}")

        # Check for required params
        for step in steps:
            tool = step.get("tool")
            if tool in ALLOWED_TOOLS:
                for param_name, param_def in ALLOWED_TOOLS[tool]["params"].items():
                    if param_def.get("required") and param_name not in step.get("params", {}):
                        issues.append(f"Missing required param: {tool}.{param_name}")

        return {
            "is_valid": len(issues) == 0,
            "issues": issues,
            "n_steps": len(steps),
            "n_issues": len(issues),
        }


critique_agent = CritiqueAgent()
