"""
BioDockify Docking Studio - RDKit Descriptor Calculator
Calculates physicochemical properties and drug-likeness rules
"""

import logging
from typing import Dict, Any, List
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from pathlib import Path

class RDKitCalculator:
    """RDKit descriptor calculator"""
    
    def __init__(self):
        """Initialize RDKit calculator"""
        pass
    
    def calculate_descriptors(self, ligand_file: str, job_id: str) -> Dict[str, Any]:
        """Calculate molecular descriptors"""
        logger.info(f"Calculating descriptors for job {job_id}")
        
        try:
            # Detect file format
            file_path = Path(ligand_file)
            file_ext = file_path.suffix.lower()
            mol = None
            
            if file_ext == '.smi' or file_ext == '.smiles':
                with open(ligand_file, 'r') as f:
                    smiles = f.read().strip()
                mol = Chem.MolFromSmiles(smiles)
            elif file_ext == '.pdb':
                mol = Chem.MolFromPDBFile(ligand_file, removeHs=False)
            elif file_ext == '.mol' or file_ext == '.sdf':
                mol = Chem.MolFromMolFile(ligand_file)
            elif file_ext == '.pdbqt':  # Not directly supported
                logger.warning(f"PDBQT format not directly supported by RDKit for {job_id}")
                return {"status": "FAILED", "error": "PDBQT format requires conversion.", "descriptors": None}
            else:
                return {"status": "FAILED", "error": f"Unsupported format: {file_ext}", "descriptors": None}

            if mol is None:
                return {"status": "FAILED", "error": "Failed to parse molecule", "descriptors": None}
            
            # Calculate descriptors
            mol_weight = Descriptors.MolWt(mol)
            log_p = Descriptors.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            num_rotatable_bonds = Descriptors.NumRotatableBonds(mol)
            num_hbd = Descriptors.NumHDonors(mol)
            num_hba = Descriptors.NumHAcceptors(mol)
            num_aromatic_rings = Descriptors.NumAromaticRings(mol)
            
            # Lipinski's Rule of Five
            violations = 0
            if mol_weight > 500:
                violations += 1
            if log_p > 5:
                violations += 1
            if num_hbd > 5:
                violations += 1
            if num_hba > 10:
                violations += 1
            if num_rotatable_bonds > 10:
                violations += 1
            
            # Drug-likeness
            if violations == 0:
                drug_likeness = "High"
            elif violations == 1:
                drug_likeness = "Good"
            elif violations == 2:
                drug_likeness = "Moderate"
            else:
                drug_likeness = "Low"
            
            # Additional descriptors
            num_fused_aromatics = Descriptors.NumAromaticRings(mol)
            num_heteroatoms = Descriptors.NumHeteroatoms(mol)
            
            descriptors_data = {
                "molecular_weight": mol_weight,
                "log_p": log_p,
                "tpsa": tpsa,
                "num_rotatable_bonds": num_rotatable_bonds,
                "num_hbd": num_hbd,
                "num_hba": num_hba,
                "num_aromatic_rings": num_aromatic_rings,
                "num_heteroatoms": num_heteroatoms,
                "lipinski_violations": violations,
                "drug_likeness": drug_likeness,
                "veber_rule": "Yes" if tpsa < 140 else "No",
                "egan_rule": "Yes" if log_p < 3 and tpsa < 130 else "No"
            }
            
            logger.info(f"RDKit descriptors calculated for job {job_id}")
            
            return {
                "status": "COMPLETED",
                "descriptors": descriptors_data
            }
            
        except Exception as e:
            logger.error(f"RDKit calculation failed: {e}")
            return {
                "status": "FAILED",
                "error": str(e),
                "descriptors": None
            }
