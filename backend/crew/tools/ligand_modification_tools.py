"""
Ligand Modification Tools
CrewAI-compatible tools for fetching similar ligands, parsing prompts,
applying transformations, validating, and docking.
"""
import logging
import re
from typing import List, Dict, Optional, Literal

from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)


# ============ PYDANTIC SCHEMAS ============
class SimilaritySearchRequest(BaseModel):
    query_smiles: str
    database: Literal["pubchem", "chembl"] = "pubchem"
    similarity_threshold: float = Field(ge=0.8, le=0.95, default=0.85)
    max_results: int = Field(ge=5, le=50, default=20)


class ModificationPlan(BaseModel):
    mode: Literal["similarity_search", "prompt_based", "autonomous"]
    strategy: Optional[str] = None
    transformations: List[str] = Field(default_factory=list)
    property_constraints: Dict[str, float] = Field(
        default_factory=lambda: {
            "mw_max": 500.0,
            "logp_max": 5.0,
            "hbd_max": 5,
            "hba_max": 10,
            "rot_bonds_max": 10,
        }
    )
    max_variants: int = Field(ge=10, le=100, default=50)
    docking_params: Dict = Field(default_factory=lambda: {"exhaustiveness": 8})


class ModificationResult(BaseModel):
    parent_smiles: str
    modified_smiles: str
    applied_transform: str
    properties: Dict[str, float]
    docking_score: Optional[float] = None
    delta_score: Optional[float] = None
    source: Literal["database", "transformation", "autonomous"]


# ============ TOOLS ============
def fetch_similar_ligands(req: SimilaritySearchRequest) -> List[Dict]:
    """Fetch similar ligands from PubChem/ChEMBL using fingerprint similarity."""
    if req.database == "pubchem":
        return _fetch_pubchem_similar(req.query_smiles, req.similarity_threshold, req.max_results)
    return []


def _fetch_pubchem_similar(query_smiles: str, threshold: float, max_results: int) -> List[Dict]:
    """PubChem similarity search via PUG-REST API."""
    import requests

    try:
        cid_resp = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{query_smiles}/cids/JSON",
            timeout=10,
        )
        if cid_resp.status_code != 200:
            return []
        cid = cid_resp.json()["IdentifierList"]["CID"][0]

        similar_resp = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/similarity/fingerprint/JSON",
            params={"Threshold": int(threshold * 100), "MaxRecords": max_results},
            timeout=30,
        )
        if similar_resp.status_code != 200:
            return []

        similar_cids = similar_resp.json()["IdentifierList"]["CID"]

        results = []
        for cid in similar_cids[:max_results]:
            smiles_resp = requests.get(
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES/JSON",
                timeout=10,
            )
            if smiles_resp.status_code == 200:
                smiles = smiles_resp.json()["PropertyTable"]["Properties"][0]["IsomericSMILES"]
                results.append({"smiles": smiles, "cid": cid, "source": "pubchem"})

        return results
    except Exception as e:
        logger.error(f"PubChem fetch failed: {e}")
        return []


def parse_modification_prompt(prompt: str, parent_smiles: str) -> ModificationPlan:
    """Convert natural language prompt into a structured modification plan."""
    from chemistry.transformation_engine import STRATEGY_TRANSFORMS

    prompt_lower = prompt.lower()

    transformations = []
    strategy = None

    if "polar" in prompt_lower or "hydrophilic" in prompt_lower:
        transformations.extend(["add_OH", "add_NH2"])
        strategy = "increase_polarity"
    if "lipophilic" in prompt_lower or "hydrophobic" in prompt_lower:
        transformations.extend(["add_CH3", "add_F", "add_Cl"])
        strategy = "increase_lipophilicity"
    if "smaller" in prompt_lower or "reduce mw" in prompt_lower or "reduce molecular weight" in prompt_lower:
        transformations.extend(["demethylate"])
        strategy = "reduce_molecular_weight"
    if "drug-like" in prompt_lower or "lipinski" in prompt_lower:
        strategy = "improve_druglikeness"
        transformations.extend(STRATEGY_TRANSFORMS.get("improve_druglikeness", []))

    constraints = {}
    mw_match = re.search(r'(?:mw|weight|mass)\s*(?:<|less than|under|max|maximum)?\s*(\d+)', prompt_lower)
    if mw_match:
        constraints["mw_max"] = float(mw_match.group(1))

    logp_match = re.search(r'logp?\s*(?:<|less than|max|maximum)?\s*([\d.]+)', prompt_lower)
    if logp_match:
        constraints["logp_max"] = float(logp_match.group(1))

    return ModificationPlan(
        mode="prompt_based",
        strategy=strategy,
        transformations=transformations[:3],
        property_constraints={
            **ModificationPlan.model_fields["property_constraints"].default,
            **constraints,
        },
        max_variants=50,
    )


def apply_rdkit_transformations(parent_smiles_list: List[str], plan: ModificationPlan) -> List[Dict]:
    """Apply RDKit reaction SMARTS to generate variants."""
    from chemistry.transformation_engine import apply_transformations_safe

    return apply_transformations_safe(parent_smiles_list, plan.transformations, plan.max_variants)


def validate_and_filter(variants: List[Dict], constraints: Dict[str, float]) -> List[Dict]:
    """Remove invalid structures and filter by property constraints."""
    from chemistry.validator import validate_and_filter_variants

    return validate_and_filter_variants(variants, constraints)
