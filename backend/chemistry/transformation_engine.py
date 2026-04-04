"""
Ligand Modification Transformation Engine
Applies safe, pre-defined RDKit reaction SMARTS to generate molecular variants.
"""
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from typing import List, Dict
import logging

logger = logging.getLogger(__name__)

# SAFE, PRE-DEFINED REACTION SMARTS (curated for drug-like chemistry)
REACTION_LIBRARY = {
    # Add functional groups
    "add_OH": rdChemReactions.ReactionFromSmarts("[#6;H3,H2:1]>>[#6:1][OH]"),
    "add_NH2": rdChemReactions.ReactionFromSmarts("[#6;H3,H2:1]>>[#6:1][NH2]"),
    "add_F": rdChemReactions.ReactionFromSmarts("[#6;H3,H2:1]>>[#6:1][F]"),
    "add_Cl": rdChemReactions.ReactionFromSmarts("[#6;H3,H2:1]>>[#6:1][Cl]"),
    "add_CH3": rdChemReactions.ReactionFromSmarts("[#6;H2:1]>>[#6:1][CH3]"),

    # Replace groups (bioisosteres)
    "replace_OH_with_NH2": rdChemReactions.ReactionFromSmarts("[#6:1][OH:2]>>[#6:1][NH2:2]"),
    "replace_H_with_F": rdChemReactions.ReactionFromSmarts("[#6:1][H:2]>>[#6:1][F:2]"),
    "replace_C_with_N": rdChemReactions.ReactionFromSmarts("[#6;H2:1]>>[#7;H1:1]"),

    # Chain modifications
    "extend_methyl": rdChemReactions.ReactionFromSmarts("[#6:1]>>[#6:1][CH3]"),
    "demethylate": rdChemReactions.ReactionFromSmarts("[#6:1][CH3:2]>>[#6:1][H:2]"),

    # Ring modifications (simple)
    "add_methyl_to_aromatic": rdChemReactions.ReactionFromSmarts("[a:1]>>[a:1][CH3]"),
}

# Strategy -> transformation mapping
STRATEGY_TRANSFORMS = {
    "increase_polarity": ["add_OH", "add_NH2"],
    "increase_lipophilicity": ["add_CH3", "add_F", "add_Cl"],
    "reduce_molecular_weight": ["demethylate"],
    "improve_druglikeness": ["add_OH", "add_F", "demethylate"],
}


def apply_transformations_safe(
    parent_smiles_list: List[str],
    transformations: List[str],
    max_variants: int = 50
) -> List[Dict]:
    """
    Apply RDKit transformations with safety checks.
    Returns list of {parent_smiles, modified_smiles, applied_transform, source}.
    """
    variants = []
    seen = set()

    for parent_smi in parent_smiles_list[:10]:
        parent_mol = Chem.MolFromSmiles(parent_smi)
        if not parent_mol:
            continue
        try:
            Chem.SanitizeMol(parent_mol)
        except Exception:
            continue

        for transform_name in transformations:
            if transform_name not in REACTION_LIBRARY:
                logger.warning(f"Unknown transform: {transform_name}")
                continue

            rxn = REACTION_LIBRARY[transform_name]
            if not rxn:
                continue

            try:
                products = rxn.RunReactants((parent_mol,))
                for prod_tuple in products[:2]:
                    if len(variants) >= max_variants:
                        break

                    prod = prod_tuple[0]
                    if not prod:
                        continue

                    if Chem.SanitizeMol(prod, catchErrors=True) != 0:
                        continue

                    mod_smi = Chem.MolToSmiles(prod, isomericSmiles=True)

                    if mod_smi in seen or mod_smi == parent_smi:
                        continue
                    seen.add(mod_smi)

                    variants.append({
                        "parent_smiles": parent_smi,
                        "modified_smiles": mod_smi,
                        "applied_transform": transform_name,
                        "source": "transformation"
                    })
            except Exception as e:
                logger.debug(f"Transform {transform_name} failed on {parent_smi}: {e}")
                continue

    return variants


def get_transforms_for_strategy(strategy: str) -> List[str]:
    """Get transformation list for a given strategy."""
    return STRATEGY_TRANSFORMS.get(strategy, [])
