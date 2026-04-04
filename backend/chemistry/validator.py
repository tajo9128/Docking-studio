"""
Ligand Modification Validator
Validates molecular structures and filters by property constraints.
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from typing import List, Dict
import logging

logger = logging.getLogger(__name__)


def validate_and_filter_variants(
    variants: List[Dict],
    constraints: Dict[str, float]
) -> List[Dict]:
    """Remove invalid structures and filter by property constraints."""
    filtered = []
    seen_smiles = set()

    for v in variants:
        smi = v.get("modified_smiles")
        if not smi:
            continue

        if smi in seen_smiles:
            continue
        seen_smiles.add(smi)

        mol = Chem.MolFromSmiles(smi)
        if not mol:
            continue
        if Chem.SanitizeMol(mol, catchErrors=True) != 0:
            continue

        props = {
            "mw": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "rot_bonds": Descriptors.NumRotatableBonds(mol),
            "tpsa": rdMolDescriptors.CalcTPSA(mol)
        }

        if (props["mw"] > constraints.get("mw_max", 500) or
            props["logp"] > constraints.get("logp_max", 5.0) or
            props["hbd"] > constraints.get("hbd_max", 5) or
            props["hba"] > constraints.get("hba_max", 10) or
            props["rot_bonds"] > constraints.get("rot_bonds_max", 10)):
            continue

        v["properties"] = props
        filtered.append(v)

    return filtered
