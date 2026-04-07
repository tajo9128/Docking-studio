"""
NL-to-DAG Workflow Compiler + Self-Healing Execution Engine
Converts natural language requests into executable, safe DAG workflows.
"""

import json
import logging
import uuid
from typing import Dict, Any, List, Optional
from datetime import datetime

logger = logging.getLogger(__name__)


# Allowed tools and their parameter schemas
ALLOWED_TOOLS = {
    "docking": {
        "params": {
            "receptor_content": {"type": "str", "required": True},
            "smiles": {"type": "str", "required": True},
            "exhaustiveness": {"type": "int", "default": 8, "min": 1, "max": 64},
            "num_modes": {"type": "int", "default": 10, "min": 1, "max": 20},
            "center_x": {"type": "float", "default": 0.0},
            "center_y": {"type": "float", "default": 0.0},
            "center_z": {"type": "float", "default": 0.0},
            "size_x": {"type": "float", "default": 20.0},
            "size_y": {"type": "float", "default": 20.0},
            "size_z": {"type": "float", "default": 20.0},
        },
        "timeout": 3600,
        "max_retries": 3,
    },
    "md": {
        "params": {
            "pdb_content": {"type": "str", "required": True},
            "steps": {"type": "int", "default": 50000, "min": 1000, "max": 1000000},
            "temperature": {"type": "float", "default": 300.0, "min": 100, "max": 500},
            "solvent_model": {"type": "str", "default": "tip3p"},
        },
        "timeout": 7200,
        "max_retries": 2,
    },
    "qsar": {
        "params": {
            "smiles_list": {"type": "list", "required": True},
            "activity_column": {"type": "str", "required": True},
            "model_type": {"type": "str", "default": "RandomForest"},
            "cv_folds": {"type": "int", "default": 5, "min": 3, "max": 10},
        },
        "timeout": 1800,
        "max_retries": 2,
    },
    "admet": {
        "params": {
            "smiles": {"type": "str", "required": True},
        },
        "timeout": 300,
        "max_retries": 3,
    },
    "pharmacophore": {
        "params": {
            "smiles": {"type": "str", "required": False},
            "pdb": {"type": "str", "required": False},
        },
        "timeout": 600,
        "max_retries": 3,
    },
    "properties": {
        "params": {
            "smiles": {"type": "str", "required": True},
        },
        "timeout": 60,
        "max_retries": 3,
    },
}

# Keyword-to-tool mapping for NL parsing
TOOL_KEYWORDS = {
    "docking": ["dock", "docking", "bind", "binding", "pose", "vina", "gnina"],
    "md": ["md", "molecular dynamics", "dynamics", "simulation", "equilibration", "trajectory"],
    "qsar": ["qsar", "model", "predict", "activity", "train", "regression", "classification"],
    "admet": ["admet", "toxicity", "absorption", "permeability", "metabolism", "adme"],
    "pharmacophore": ["pharmacophore", "ph4", "features", "screen", "library"],
    "properties": ["properties", "descriptors", "mw", "logp", "tpsa", "lipinski"],
}


class NLWorkflowCompiler:
    """
    Compiles natural language requests into safe, executable DAG workflows.
    """

    def compile(self, natural_language: str, context: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        Parse natural language and build a workflow DAG.
        Uses keyword-based heuristic parsing (no LLM dependency).
        """
        nl_lower = natural_language.lower()
        steps = []
        step_id = 0

        # Detect which tools are needed
        detected_tools = []
        for tool, keywords in TOOL_KEYWORDS.items():
            if any(kw in nl_lower for kw in keywords):
                detected_tools.append(tool)

        # Define execution order
        tool_order = ["properties", "pharmacophore", "docking", "md", "qsar", "admet"]
        ordered_tools = [t for t in tool_order if t in detected_tools]

        if not ordered_tools:
            # Default: properties + docking
            ordered_tools = ["properties", "docking"]

        for tool in ordered_tools:
            step_id += 1
            step = {
                "id": f"step_{step_id}",
                "tool": tool,
                "params": self._extract_params(tool, natural_language, context),
                "depends_on": [f"step_{step_id - 1}"] if step_id > 1 else [],
                "timeout": ALLOWED_TOOLS[tool]["timeout"],
                "max_retries": ALLOWED_TOOLS[tool]["max_retries"],
            }
            steps.append(step)

        return {
            "workflow_id": f"wf_{uuid.uuid4().hex[:8]}",
            "natural_language": natural_language,
            "steps": steps,
            "created_at": datetime.now().isoformat(),
            "status": "compiled",
        }

    def _extract_params(self, tool: str, nl: str, context: Dict[str, Any] = None) -> Dict[str, Any]:
        """Extract parameters from natural language and context."""
        params = {}
        context = context or {}
        nl_lower = nl.lower()

        tool_schema = ALLOWED_TOOLS.get(tool, {})
        param_defs = tool_schema.get("params", {})

        for param_name, param_def in param_defs.items():
            # Check context first
            if param_name in context:
                params[param_name] = context[param_name]
                continue

            # Check for explicit values in NL
            default = param_def.get("default")
            if default is not None:
                params[param_name] = default

            # Extract numeric values from NL
            if param_def.get("type") in ("int", "float"):
                import re
                numbers = re.findall(r'(\d+\.?\d*)', nl_lower)
                if numbers:
                    val = float(numbers[-1]) if '.' in numbers[-1] else int(numbers[-1])
                    min_val = param_def.get("min", float('-inf'))
                    max_val = param_def.get("max", float('inf'))
                    if min_val <= val <= max_val:
                        params[param_name] = val

        return params

    def validate_and_secure(self, dag: Dict[str, Any]) -> Dict[str, Any]:
        """Validate DAG structure and enforce safety limits."""
        errors = []

        for step in dag.get("steps", []):
            tool = step.get("tool")
            if tool not in ALLOWED_TOOLS:
                errors.append(f"Unknown tool: {tool}")
                continue

            # Enforce timeouts
            step["timeout"] = min(step.get("timeout", 3600), ALLOWED_TOOLS[tool]["timeout"])
            step["max_retries"] = min(step.get("max_retries", 3), ALLOWED_TOOLS[tool]["max_retries"])

            # Validate params
            param_defs = ALLOWED_TOOLS[tool]["params"]
            for param_name, param_def in param_defs.items():
                if param_def.get("required") and param_name not in step.get("params", {}):
                    errors.append(f"Missing required param: {tool}.{param_name}")

                if param_name in step.get("params", {}):
                    val = step["params"][param_name]
                    if "min" in param_def and val < param_def["min"]:
                        step["params"][param_name] = param_def["min"]
                    if "max" in param_def and val > param_def["max"]:
                        step["params"][param_name] = param_def["max"]

        dag["validation_errors"] = errors
        dag["is_safe"] = len(errors) == 0
        return dag


class SelfHealingExecutor:
    """
    Executes workflow steps with automatic error diagnosis and recovery.
    """

    def __init__(self):
        self.execution_log: List[Dict[str, Any]] = []

    def execute_step(self, step: Dict[str, Any]) -> Dict[str, Any]:
        """Execute a single workflow step with retry and self-healing."""
        tool = step.get("tool")
        params = step.get("params", {})
        max_retries = step.get("max_retries", 3)
        timeout = step.get("timeout", 3600)

        last_error = None
        for attempt in range(max_retries + 1):
            try:
                result = self._call_tool(tool, params)
                self.execution_log.append({
                    "step_id": step.get("id"),
                    "tool": tool,
                    "attempt": attempt,
                    "status": "success",
                    "result": result,
                    "timestamp": datetime.now().isoformat(),
                })
                return {"status": "success", "result": result, "attempts": attempt + 1}

            except Exception as e:
                last_error = str(e)
                self.execution_log.append({
                    "step_id": step.get("id"),
                    "tool": tool,
                    "attempt": attempt,
                    "status": "failed",
                    "error": last_error,
                    "timestamp": datetime.now().isoformat(),
                })

                if attempt < max_retries:
                    # Self-heal: adjust parameters
                    params = self._heal_params(tool, params, last_error)
                    logger.info(f"Self-healing attempt {attempt + 1} for {tool}: {last_error}")

        return {
            "status": "failed",
            "error": last_error,
            "attempts": max_retries + 1,
            "execution_log": self.execution_log[-(max_retries + 1):],
        }

    def _call_tool(self, tool: str, params: Dict[str, Any]) -> Dict[str, Any]:
        """Call the appropriate backend tool."""
        if tool == "docking":
            from docking_engine import smart_dock
            return smart_dock(
                receptor_content=params.get("receptor_content", ""),
                ligand_content=params.get("smiles", ""),
                input_format="smiles",
                center_x=params.get("center_x", 0),
                center_y=params.get("center_y", 0),
                center_z=params.get("center_z", 0),
                size_x=params.get("size_x", 20),
                size_y=params.get("size_y", 20),
                size_z=params.get("size_z", 20),
                exhaustiveness=params.get("exhaustiveness", 8),
                num_modes=params.get("num_modes", 10),
            )
        elif tool == "properties":
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski, Crippen, rdMolDescriptors
            mol = Chem.MolFromSmiles(params.get("smiles", ""))
            if mol is None:
                raise ValueError("Invalid SMILES")
            return {
                "mw": Descriptors.MolWt(mol),
                "logp": Crippen.MolLogP(mol),
                "hbd": Lipinski.NumHDonors(mol),
                "hba": Lipinski.NumHAcceptors(mol),
                "tpsa": rdMolDescriptors.CalcTPSA(mol),
                "rotatable": rdMolDescriptors.CalcNumRotatableBonds(mol),
            }
        elif tool == "admet":
            return {"status": "admet_prediction", "smiles": params.get("smiles", "")}
        elif tool == "pharmacophore":
            from pharmacophore import get_engine
            engine = get_engine()
            if params.get("smiles"):
                return engine.generate_from_smiles(params["smiles"])
            elif params.get("pdb"):
                return engine.generate_from_pdb(params["pdb"])
            return {"error": "No input provided"}
        else:
            return {"status": f"{tool}_executed", "params": params}

    def _heal_params(self, tool: str, params: Dict[str, Any], error: str) -> Dict[str, Any]:
        """Diagnose error and adjust parameters for retry."""
        error_lower = error.lower()
        healed = params.copy()

        if tool == "docking":
            if "grid" in error_lower or "box" in error_lower:
                healed["size_x"] = params.get("size_x", 20) + 4.0
                healed["size_y"] = params.get("size_y", 20) + 4.0
                healed["size_z"] = params.get("size_z", 20) + 4.0
            elif "exhaustiveness" in error_lower or "timeout" in error_lower:
                healed["exhaustiveness"] = max(params.get("exhaustiveness", 8) // 2, 1)
            elif "clash" in error_lower or "steric" in error_lower:
                healed["size_x"] = params.get("size_x", 20) + 2.0
                healed["size_y"] = params.get("size_y", 20) + 2.0
                healed["size_z"] = params.get("size_z", 20) + 2.0
        elif tool == "md":
            if "timeout" in error_lower or "memory" in error_lower:
                healed["steps"] = max(params.get("steps", 50000) // 2, 1000)
            elif "unstable" in error_lower:
                healed["temperature"] = params.get("temperature", 300) - 50

        return healed


nl_compiler = NLWorkflowCompiler()
self_healing_executor = SelfHealingExecutor()
