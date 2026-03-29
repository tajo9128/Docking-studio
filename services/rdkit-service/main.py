"""
RDKit Service - Molecule processing API
"""

import os
import logging
from pathlib import Path
from typing import Optional, List, Dict, Any
from datetime import datetime
from contextlib import asynccontextmanager

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("rdkit-service")


class MolInput(BaseModel):
    smiles: Optional[str] = None
    pdb: Optional[str] = None


class ConvertRequest(BaseModel):
    input_format: str
    output_format: str
    content: str


class PropertiesRequest(BaseModel):
    smiles: str


@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("RDKit Service starting up...")
    yield
    logger.info("RDKit Service shutting down...")


app = FastAPI(
    title="RDKit Service API",
    description="RDKit molecule processing",
    version="2.0.0",
    lifespan=lifespan,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/health")
def health():
    return {
        "status": "healthy",
        "service": "rdkit-service",
        "timestamp": datetime.now().isoformat(),
    }


@app.get("/")
def root():
    return {"service": "rdkit-service", "version": "2.0.0"}


@app.post("/from_smiles")
def from_smiles(request: MolInput):
    """Convert SMILES to molecule"""
    try:
        from rdkit import Chem

        if not request.smiles:
            return {"success": False, "error": "No SMILES provided"}

        mol = Chem.MolFromSmiles(request.smiles)
        if mol is None:
            return {"success": False, "error": "Invalid SMILES"}

        return {
            "success": True,
            "smiles": request.smiles,
            "num_atoms": mol.GetNumAtoms(),
            "num_bonds": mol.GetNumBonds(),
        }

    except Exception as e:
        logger.error(f"SMILES conversion error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/from_pdb")
def from_pdb(request: MolInput):
    """Convert PDB to molecule"""
    try:
        from rdkit import Chem

        if not request.pdb:
            return {"success": False, "error": "No PDB provided"}

        mol = Chem.MolFromPDBBlock(request.pdb)
        if mol is None:
            return {"success": False, "error": "Invalid PDB"}

        return {
            "success": True,
            "num_atoms": mol.GetNumAtoms(),
            "num_bonds": mol.GetNumBonds(),
        }

    except Exception as e:
        logger.error(f"PDB conversion error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/properties")
def molecular_properties(request: PropertiesRequest):
    """Calculate molecular properties"""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski

        mol = Chem.MolFromSmiles(request.smiles)
        if mol is None:
            return {"success": False, "error": "Invalid SMILES"}

        return {
            "success": True,
            "smiles": request.smiles,
            "mw": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "num_h_donors": Lipinski.NumHDonors(mol),
            "num_h_acceptors": Lipinski.NumHAcceptors(mol),
            "num_rotatable_bonds": Lipinski.NumRotatableBonds(mol),
            "tpsa": Descriptors.TPSA(mol),
        }

    except Exception as e:
        logger.error(f"Property calculation error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/similarity")
def similarity(smiles1: str, smiles2: str):
    """Calculate Tanimoto similarity between two molecules"""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)

        if mol1 is None or mol2 is None:
            return {"success": False, "error": "Invalid SMILES"}

        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)

        from rdkit import DataStructs

        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)

        return {
            "success": True,
            "smiles1": smiles1,
            "smiles2": smiles2,
            "tanimoto_similarity": similarity,
        }

    except Exception as e:
        logger.error(f"Similarity calculation error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


class SmilesTo3DRequest(BaseModel):
    smiles: str
    name: Optional[str] = "molecule"


@app.post("/smiles-to-3d")
def smiles_to_3d(request: SmilesTo3DRequest):
    """Convert SMILES to 3D structure and save PDB file"""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        import uuid

        mol = Chem.MolFromSmiles(request.smiles)
        if mol is None:
            return {"success": False, "error": "Invalid SMILES"}

        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)

        name = request.name or "molecule"
        safe_name = "".join(c if c.isalnum() else "_" for c in name)
        pdb_id = str(uuid.uuid4())[:8]
        filename = f"{safe_name}_{pdb_id}.pdb"
        output_path = Path("/app/storage") / filename
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w") as f:
            f.write(Chem.MolToPDBBlock(mol))

        return {
            "success": True,
            "pdb_path": str(output_path),
            "num_atoms": mol.GetNumAtoms(),
            "smiles": request.smiles,
        }

    except Exception as e:
        logger.error(f"SMILES to 3D error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


class OptimizeRequest(BaseModel):
    pdb_path: str


@app.post("/optimize")
def optimize_molecule(request: OptimizeRequest):
    """Optimize 3D geometry of a molecule"""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        if not os.path.exists(request.pdb_path):
            return {"success": False, "error": "File not found"}

        mol = Chem.MolFromPDBFile(request.pdb_path)
        if mol is None:
            return {"success": False, "error": "Invalid PDB file"}

        mol = Chem.AddHs(mol)
        AllChem.MMFFOptimizeMolecule(mol)

        parts = request.pdb_path.rsplit(".", 1)
        optimized_path = f"{parts[0]}_optimized.pdb"
        with open(optimized_path, "w") as f:
            f.write(Chem.MolToPDBBlock(mol))

        return {
            "success": True,
            "output_path": optimized_path,
            "num_atoms": mol.GetNumAtoms(),
        }

    except Exception as e:
        logger.error(f"Optimization error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


class ConvertRequest(BaseModel):
    input_path: str
    output_format: str


@app.post("/convert")
def convert_format(request: ConvertRequest):
    """Convert molecule between formats"""
    try:
        from rdkit import Chem

        ext = os.path.splitext(request.input_path)[1].lower()
        if ext == ".pdb":
            mol = Chem.MolFromPDBFile(request.input_path)
        elif ext in [".mol", ".sdf"]:
            supplier = Chem.SDMolSupplier(request.input_path)
            mol = next(supplier, None)
        elif ext == ".mol2":
            mol = Chem.MolFromMol2File(request.input_path)
        else:
            return {"success": False, "error": f"Unsupported input format: {ext}"}

        if mol is None:
            return {"success": False, "error": "Failed to read molecule"}

        out_name = os.path.splitext(os.path.basename(request.input_path))[0]
        output_path = f"/app/storage/{out_name}.{request.output_format}"

        if request.output_format == "pdb":
            content = Chem.MolToPDBBlock(mol)
            with open(output_path, "w") as f:
                f.write(content)
        elif request.output_format in ["mol", "sdf"]:
            with open(output_path, "w") as f:
                f.write(Chem.MolToMolBlock(mol))
        elif request.output_format == "mol2":
            with open(output_path, "w") as f:
                f.write(Chem.MolToMol2Block(mol))
        else:
            return {
                "success": False,
                "error": f"Unsupported output format: {request.output_format}",
            }

        return {
            "success": True,
            "output_path": output_path,
            "num_atoms": mol.GetNumAtoms(),
        }

    except Exception as e:
        logger.error(f"Conversion error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


class ProcessRequest(BaseModel):
    smiles: Optional[str] = None
    pdb: Optional[str] = None
    operation: str = "properties"


@app.post("/process")
def process_molecule(request: ProcessRequest):
    """Process molecule with specified operation"""
    try:
        if request.operation == "smiles-to-3d" and request.smiles:
            r = SmilesTo3DRequest(smiles=request.smiles)
            return smiles_to_3d(r)
        elif request.operation == "optimize" and request.pdb:
            r = OptimizeRequest(pdb_path=request.pdb)
            return optimize_molecule(r)
        elif request.operation == "properties" and request.smiles:
            r = PropertiesRequest(smiles=request.smiles)
            return molecular_properties(r)
        else:
            return {"success": False, "error": "Invalid operation or missing input"}
    except Exception as e:
        logger.error(f"Process error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


class PrepareProteinRequest(BaseModel):
    pdb_content: str
    name: Optional[str] = "protein"
    remove_waters: bool = True
    add_hydrogens: bool = True
    pH: float = 7.4


@app.post("/prepare_protein")
def prepare_protein(request: PrepareProteinRequest):
    """Prepare protein: remove waters, add hydrogens, return cleaned PDB"""
    try:
        from rdkit import Chem
        import uuid

        mol = Chem.MolFromPDBBlock(request.pdb_content)
        if mol is None:
            return {"success": False, "error": "Invalid PDB content"}

        original_atoms = mol.GetNumAtoms()

        if request.remove_waters:
            waters = []
            for atom in mol.GetAtoms():
                try:
                    elem = atom.GetSymbol()
                    degree = atom.GetTotalDegree()
                    if elem == "O" and degree == 1:
                        waters.append(atom.GetIdx())
                except:
                    pass
            if waters:
                editor = Chem.RWMol(mol)
                for idx in sorted(waters, reverse=True):
                    editor.RemoveAtom(idx)
                mol = editor.GetMol()

        if request.add_hydrogens:
            mol = Chem.AddHs(mol, addCoords=True)

        safe_name = "".join(
            c if c.isalnum() else "_" for c in (request.name or "protein")
        )
        pdb_id = str(uuid.uuid4())[:8]
        filename = f"{safe_name}_prep_{pdb_id}.pdb"
        output_path = Path("/app/storage") / filename
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w") as f:
            f.write(Chem.MolToPDBBlock(mol))

        return {
            "success": True,
            "pdb_path": str(output_path),
            "original_atoms": original_atoms,
            "final_atoms": mol.GetNumAtoms(),
            "waters_removed": request.remove_waters,
            "hydrogens_added": request.add_hydrogens,
        }

    except Exception as e:
        logger.error(f"Protein preparation error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


class PrepareLigandRequest(BaseModel):
    pdb_content: str
    name: Optional[str] = "ligand"
    pH: float = 7.4


@app.post("/prepare_ligand")
def prepare_ligand(request: PrepareLigandRequest):
    """Prepare ligand: add hydrogens, return PDB and PDBQT"""
    try:
        from rdkit import Chem
        import uuid

        mol = Chem.MolFromPDBBlock(request.pdb_content)
        if mol is None:
            return {"success": False, "error": "Invalid PDB content"}

        mol = Chem.AddHs(mol, addCoords=True)

        try:
            from meeko import MoleculePreparation, PDBQTWriterLegacy

            prep = MoleculePreparation()
            mol_set_list = prep.prepare(mol)
            pdbqt_string = None
            for setup in mol_set_list:
                result = PDBQTWriterLegacy.write_string(setup)
                pdbqt_string = result[0] if isinstance(result, tuple) else result
                break

            if pdbqt_string is None:
                raise ValueError("Meeko preparation failed")

            safe_name = "".join(
                c if c.isalnum() else "_" for c in (request.name or "ligand")
            )
            pdb_id = str(uuid.uuid4())[:8]
            pdbqt_filename = f"{safe_name}_ligand_{pdb_id}.pdbqt"
            pdbqt_path = Path("/app/storage") / pdbqt_filename
            pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
            with open(pdbqt_path, "w") as f:
                f.write(pdbqt_string)

            return {
                "success": True,
                "pdbqt_path": str(pdbqt_path),
                "num_atoms": mol.GetNumAtoms(),
                "meeko_used": True,
            }
        except (ImportError, AttributeError, ValueError):
            safe_name = "".join(
                c if c.isalnum() else "_" for c in (request.name or "ligand")
            )
            pdb_id = str(uuid.uuid4())[:8]
            pdb_filename = f"{safe_name}_ligand_{pdb_id}.pdb"
            pdb_path = Path("/app/storage") / pdb_filename
            pdb_path.parent.mkdir(parents=True, exist_ok=True)
            with open(pdb_path, "w") as f:
                f.write(Chem.MolToPDBBlock(mol))

            return {
                "success": True,
                "pdb_path": str(pdb_path),
                "num_atoms": mol.GetNumAtoms(),
                "meeko_used": False,
                "message": "Meeko not available, PDB with hydrogens returned. Install meeko for PDBQT.",
            }

    except Exception as e:
        logger.error(f"Ligand preparation error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


class DetectInteractionsRequest(BaseModel):
    receptor_pdb_content: str
    ligand_pdb_content: str


@app.post("/detect_interactions")
def detect_interactions(request: DetectInteractionsRequest):
    """Detect H-bonds and hydrophobic interactions between receptor and ligand"""
    try:
        from rdkit import Chem
        from rdkit.Chem import Lipinski
        import numpy as np

        receptor = Chem.MolFromPDBBlock(request.receptor_pdb_content)
        ligand = Chem.MolFromPDBBlock(request.ligand_pdb_content)

        if receptor is None or ligand is None:
            return {"success": False, "error": "Invalid PDB content"}

        receptor_conf = receptor.GetConformer(0)
        ligand_conf = ligand.GetConformer(0)

        h_bonds = []
        hydrophobic = []

        for ligand_atom in ligand.GetAtoms():
            elem = ligand_atom.GetSymbol()
            if elem not in ["N", "O", "F", "P", "S"]:
                continue
            ligand_idx = ligand_atom.GetIdx()
            ligand_pos = np.array(ligand_conf.GetAtomPosition(ligand_idx))

            min_dist = float("inf")
            closest_receptor_atom = None

            for receptor_atom in receptor.GetAtoms():
                rec_elem = receptor_atom.GetSymbol()
                if rec_elem not in ["C", "N", "O", "F", "P", "S"]:
                    continue
                receptor_idx = receptor_atom.GetIdx()
                receptor_pos = np.array(receptor_conf.GetAtomPosition(receptor_idx))
                dist = np.linalg.norm(ligand_pos - receptor_pos)
                if dist < min_dist:
                    min_dist = dist
                    closest_receptor_atom = receptor_atom
                    closest_rec_elem = rec_elem
                    closest_rec_idx = receptor_idx
                    closest_rec_pos = receptor_pos

            if closest_receptor_atom is None:
                continue

            if elem in ["N", "O"] and closest_rec_elem == "N":
                if min_dist < 3.5:
                    h_bonds.append(
                        {
                            "ligand_atom": elem,
                            "receptor_atom": closest_rec_elem,
                            "distance_A": round(min_dist, 2),
                            "type": "H-bond",
                            "ligand_pos": ligand_pos.tolist(),
                            "receptor_pos": closest_rec_pos.tolist(),
                        }
                    )
            elif elem == "C" and closest_rec_elem == "C":
                if min_dist < 4.5:
                    hydrophobic.append(
                        {
                            "ligand_atom": elem,
                            "receptor_atom": closest_rec_elem,
                            "distance_A": round(min_dist, 2),
                            "type": "hydrophobic",
                            "ligand_pos": ligand_pos.tolist(),
                            "receptor_pos": closest_rec_pos.tolist(),
                        }
                    )

        return {
            "success": True,
            "h_bonds": h_bonds,
            "hydrophobic_contacts": hydrophobic,
            "total_h_bonds": len(h_bonds),
            "total_hydrophobic": len(hydrophobic),
        }

    except Exception as e:
        logger.error(f"Interaction detection error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8003)
