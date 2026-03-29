"""
Molecular descriptor calculation using RDKit
"""

import logging
from typing import List, Dict, Any, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, AllChem

logger = logging.getLogger("qsar-descriptors")

DESCRIPTOR_GROUPS = {
    "all": [
        "MolWt",
        "MolLogP",
        "NumHDonors",
        "NumHAcceptors",
        "NumRotatableBonds",
        "TPSA",
        "NumAromaticRings",
        "NumHeteroatoms",
        "FractionCSP3",
        "NumAliphaticRings",
        "RingCount",
        "HeavyAtomCount",
        "NumSaturatedRings",
        "MolMR",
        "MaxPartialCharge",
        "MinPartialCharge",
        "MaxAbsPartialCharge",
        "MinAbsPartialCharge",
        "NumValenceElectrons",
        "NumRadicalElectrons",
        "BalabanJ",
        "Chi0",
        "Chi1",
        "Chi0n",
        "Chi1n",
        "HallKierAlpha",
        "Kappa1",
        "Kappa2",
        "Kappa3",
        "LabuteASA",
        "NOCount",
        "NOCharge",
        "NHOHCount",
        "NHOHCharge",
        "NumAromaticHeterocycles",
        "NumAromaticCarbocycles",
        "NumAliphaticHeterocycles",
        "NumAliphaticCarbocycles",
        "NumSaturatedHeterocycles",
        "NumSaturatedCarbocycles",
        "NumHAromatic",
        "NumHAliphatic",
        "NumAtoms",
    ],
    "lipinski": [
        "MolWt",
        "MolLogP",
        "NumHDonors",
        "NumHAcceptors",
        "NumRotatableBonds",
        "TPSA",
    ],
    "electronic": [
        "MolLogP",
        "TPSA",
        "MolMR",
        "NumHeteroatoms",
        "MaxPartialCharge",
        "MinPartialCharge",
        "MaxAbsPartialCharge",
        "MinAbsPartialCharge",
        "NumValenceElectrons",
        "NumRadicalElectrons",
    ],
    "structural": [
        "NumAtoms",
        "NumBonds",
        "NumRotatableBonds",
        "RingCount",
        "NumAromaticRings",
        "NumAliphaticRings",
        "NumSaturatedRings",
        "HeavyAtomCount",
        "NumAromaticHeterocycles",
        "NumAromaticCarbocycles",
        "NumAliphaticHeterocycles",
        "NumAliphaticCarbocycles",
        "NumSaturatedHeterocycles",
        "NumSaturatedCarbocycles",
        "FractionCSP3",
        "NumHeteroatoms",
    ],
    "morgan_fp": ["MorganFP"],
    "maccs_keys": ["MACCSKeys"],
}


def smiles_to_mol(smiles: str) -> Optional[Chem.Mol]:
    """Parse SMILES to RDKit mol object with error handling"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Invalid SMILES: {smiles[:50]}")
        return mol
    except Exception as e:
        logger.warning(f"SMILES parse error '{smiles[:30]}': {e}")
        return None


def get_available_descriptors() -> List[str]:
    """Return list of all available descriptor names"""
    descriptors = []
    for name in dir(Descriptors):
        if not name.startswith("_") and name[0].isupper():
            try:
                getattr(Descriptors, name)(None)
                descriptors.append(name)
            except Exception:
                pass
    return sorted(descriptors)


def calculate_descriptors(mol: Chem.Mol) -> Dict[str, float]:
    """Calculate all descriptors for a molecule"""
    result = {}
    for name in dir(Descriptors):
        if not name.startswith("_") and name[0].isupper():
            try:
                val = getattr(Descriptors, name)(mol)
                if isinstance(val, (int, float)) and not isinstance(val, bool):
                    result[name] = float(val)
            except Exception:
                pass

    for name in [
        "NumHDonors",
        "NumHAcceptors",
        "NumRotatableBonds",
        "NumAromaticRings",
        "NumHeteroatoms",
        "NumAromaticHeterocycles",
        "NumAromaticCarbocycles",
        "NumAliphaticHeterocycles",
        "NumAliphaticCarbocycles",
        "NumSaturatedHeterocycles",
        "NumSaturatedCarbocycles",
        "FractionCSP3",
        "NOCount",
        "NHOHCount",
        "NumValenceElectrons",
        "NumRadicalElectrons",
        "NumAliphaticRings",
        "RingCount",
        "HeavyAtomCount",
        "NumSaturatedRings",
    ]:
        try:
            result[name] = float(getattr(Lipinski, name)(mol))
        except Exception:
            pass

    try:
        result["TPSA"] = float(rdMolDescriptors.CalcTPSA(mol))
    except Exception:
        pass

    try:
        result["MolMR"] = float(rdMolDescriptors.CalcMolRef(mol))
    except Exception:
        pass

    result["NumAtoms"] = float(mol.GetNumAtoms())
    result["NumBonds"] = float(mol.GetNumBonds())

    return result


def calculate_fingerprints(mol: Chem.Mol, fp_type: str = "morgan") -> Dict[str, Any]:
    """Calculate molecular fingerprints"""
    if fp_type == "morgan":
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        arr = [0] * 2048
        for idx in fp.GetOnBits():
            arr[idx] = 1
        return {"MorganFP": arr}
    elif fp_type == "maccs":
        from rdkit.Chem import MACCSkeys

        fp = MACCSkeys.GenMACCSKeys(mol)
        arr = [0] * 167
        for idx in fp.GetOnBits():
            arr[idx] = 1
        return {"MACCSKeys": arr}
    return {}


def descriptors_for_smiles_list(
    smiles_list: List[str],
    descriptor_groups: Optional[List[str]] = None,
    include_fingerprints: bool = False,
) -> Tuple[List[Dict[str, float]], List[str], List[str]]:
    """
    Calculate descriptors for a list of SMILES.

    Returns:
        (descriptor_matrix, valid_smiles, failed_smiles)
    """
    if descriptor_groups is None:
        descriptor_groups = ["all"]

    descriptor_names = []
    for group in descriptor_groups:
        if group in DESCRIPTOR_GROUPS:
            descriptor_names.extend(DESCRIPTOR_GROUPS[group])
        elif group in ["morgan_fp", "maccs_keys"]:
            descriptor_names.append(group)
        else:
            descriptor_names.append(group)

    descriptor_names = list(set(descriptor_names))
    if include_fingerprints:
        descriptor_names.append("MorganFP")
        descriptor_names.append("MACCSKeys")

    all_descriptors = []
    valid_smiles = []
    failed_smiles = []

    for smiles in smiles_list:
        mol = smiles_to_mol(smiles)
        if mol is None:
            failed_smiles.append(smiles)
            continue

        desc = calculate_descriptors(mol)
        if include_fingerprints:
            desc.update(calculate_fingerprints(mol))

        filtered_desc = {}
        for name in descriptor_names:
            if name in desc:
                filtered_desc[name] = desc[name]
            elif name == "MorganFP":
                fp_dict = calculate_fingerprints(mol, "morgan")
                filtered_desc["MorganFP"] = fp_dict.get("MorganFP", [])
            elif name == "MACCSKeys":
                fp_dict = calculate_fingerprints(mol, "maccs")
                filtered_desc["MACCSKeys"] = fp_dict.get("MACCSKeys", [])

        all_descriptors.append(filtered_desc)
        valid_smiles.append(smiles)

    return all_descriptors, valid_smiles, failed_smiles


def descriptors_upload(
    csv_content: str,
    smiles_col: str,
    activity_col: str,
    descriptor_groups: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Process an uploaded CSV dataset.

    Args:
        csv_content: Raw CSV text
        smiles_col: Column name for SMILES
        activity_col: Column name for activity values
        descriptor_groups: List of descriptor group names

    Returns:
        Dataset info with descriptor matrix and statistics
    """
    import io
    import pandas as pd

    try:
        df = pd.read_csv(io.StringIO(csv_content))
    except Exception as e:
        raise ValueError(f"Failed to parse CSV: {e}")

    if smiles_col not in df.columns:
        raise ValueError(
            f"SMILES column '{smiles_col}' not found. Available: {list(df.columns)}"
        )
    if activity_col not in df.columns:
        raise ValueError(
            f"Activity column '{activity_col}' not found. Available: {list(df.columns)}"
        )

    smiles_list = df[smiles_col].astype(str).tolist()
    activities = pd.to_numeric(df[activity_col], errors="coerce").tolist()

    descriptors, valid_smiles, failed_smiles = descriptors_for_smiles_list(
        smiles_list, descriptor_groups
    )

    valid_indices = [i for i, s in enumerate(smiles_list) if s in valid_smiles]
    valid_activities = [activities[i] for i in valid_indices]

    activity_series = pd.Series(valid_activities)
    nan_count = activity_series.isna().sum()
    valid_pairs = [
        (d, a)
        for d, a in zip(descriptors, valid_activities)
        if a is not None and (a == a)
    ]

    if not valid_pairs:
        raise ValueError(
            "No valid (descriptor, activity) pairs found after filtering NaN activities."
        )

    clean_descriptors = [p[0] for p in valid_pairs]
    clean_activities = [p[1] for p in valid_pairs]

    feature_names = list(clean_descriptors[0].keys()) if clean_descriptors else []

    import numpy as np

    X = np.array([[d.get(fn, 0.0) for fn in feature_names] for d in clean_descriptors])
    y = np.array(clean_activities)

    nan_mask = np.isnan(X).any(axis=1)
    if nan_mask.any():
        logger.warning(f"Removing {nan_mask.sum()} rows with NaN values")
        X = X[~nan_mask]
        y = y[~nan_mask]

    return {
        "X": X.tolist(),
        "y": y.tolist(),
        "feature_names": feature_names,
        "n_compounds": len(y),
        "n_features": len(feature_names),
        "activity_mean": float(np.mean(y)),
        "activity_std": float(np.std(y)),
        "activity_min": float(np.min(y)),
        "activity_max": float(np.max(y)),
        "nan_count": int(nan_count),
        "failed_smiles": failed_smiles,
        "failed_count": len(failed_smiles),
    }
