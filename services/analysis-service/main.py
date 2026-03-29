"""
Analysis Service — Insight Generator
Responsibilities:
  - Rank ligands by binding affinity, stability, ADMET
  - Compare docking poses (RMSD clustering, pose similarity)
  - Consensus scoring across multiple methods
  - Multi-ligand comparison and filtering
  - ADMET property prediction
  - Generate analysis reports (JSON/text)
"""

import os
import logging
import json
import io
import base64
from datetime import datetime
from typing import Dict, Any, Optional, List
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from contextlib import asynccontextmanager

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("analysis-service")

STORAGE_DIR = Path("/app/storage/analysis")
STORAGE_DIR.mkdir(parents=True, exist_ok=True)


class LigandRecord(BaseModel):
    ligand_id: str
    smiles: Optional[str] = None
    vina_score: Optional[float] = None
    gnina_score: Optional[float] = None
    rf_score: Optional[float] = None
    consensus_score: Optional[float] = None
    md_stability: Optional[float] = None
    md_time_ns: Optional[float] = None
    rmsd_from_crystal: Optional[float] = None
    h_bond_count: Optional[int] = None
    hydrophobic_count: Optional[int] = None
    mw: Optional[float] = None
    logp: Optional[float] = None
    tpsa: Optional[float] = None
    num_rotatable_bonds: Optional[int] = None


class RankingRequest(BaseModel):
    ligands: List[LigandRecord]
    weights: Optional[Dict[str, float]] = None


class ComparisonRequest(BaseModel):
    job_uuid: str
    ligand_ids: List[str]


app = FastAPI(
    title="Analysis Service",
    version="2.0.0",
    description="Insight Generator: ranking, consensus scoring, multi-ligand comparison",
)


@app.get("/health")
def health():
    return {
        "service": "analysis-service",
        "status": "healthy",
        "timestamp": datetime.now().isoformat(),
    }


@app.get("/")
def root():
    return {"service": "Analysis Service (Insight Generator)", "version": "2.0.0"}


@app.post("/rank")
def rank_ligands(request: RankingRequest):
    """
    Rank a list of ligands using weighted consensus scoring.
    Weights default: vina=0.4, md_stability=0.3, admet=0.3
    """
    ligands = request.ligands
    if not ligands:
        return {"ranked": [], "message": "No ligands provided"}

    weights = request.weights or {
        "vina_score": 0.40,
        "md_stability": 0.30,
        "admet": 0.30,
    }

    scored = []
    for lig in ligands:
        s = {}
        s["ligand_id"] = lig.ligand_id

        vina_norm = _normalize(-lig.vina_score if lig.vina_score else 0, -15, 0)
        s["vina_norm"] = vina_norm

        md_norm = lig.md_stability if lig.md_stability else 0.5
        s["md_norm"] = md_norm

        admet_score = _admet_score(lig)
        s["admet_norm"] = admet_score

        total = (
            weights.get("vina_score", 0.4) * vina_norm
            + weights.get("md_stability", 0.3) * md_norm
            + weights.get("admet", 0.3) * admet_score
        )
        s["consensus_score"] = round(total, 4)
        s["smiles"] = lig.smiles
        s["details"] = {
            "vina_score": lig.vina_score,
            "gnina_score": lig.gnina_score,
            "md_stability": lig.md_stability,
            "mw": lig.mw,
            "logp": lig.logp,
            "tpsa": lig.tpsa,
            "h_bond_count": lig.h_bond_count,
            "hydrophobic_count": lig.hydrophobic_count,
        }
        scored.append(s)

    scored.sort(key=lambda x: x["consensus_score"], reverse=True)
    for i, s in enumerate(scored):
        s["rank"] = i + 1

    return {
        "ranked": scored,
        "weights_used": weights,
        "count": len(scored),
        "top_ligand": scored[0]["ligand_id"] if scored else None,
    }


@app.post("/compare/poses")
def compare_poses(request: ComparisonRequest):
    """
    Compare docking poses — compute pairwise RMSD and group by similarity.
    Uses RDKit for molecule loading and FMCat RMSD.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors

        rmsd_matrix = []
        mol_list = []
        ids = request.ligand_ids

        for lid in ids:
            sdf_path = STORAGE_DIR / f"pose_{request.job_uuid}_{lid}.sdf"
            if not sdf_path.exists():
                continue
            suppl = Chem.SDMolSupplier(str(sdf_path))
            for mol in suppl:
                if mol is not None:
                    mol.SetProp("_Name", lid)
                    mol_list.append(mol)
                    break

        n = len(mol_list)
        for i in range(n):
            row = []
            for j in range(n):
                if i == j:
                    row.append(0.0)
                elif j < i:
                    row.append(rmsd_matrix[j][i])
                else:
                    rmsd = _compute_rmsd(mol_list[i], mol_list[j])
                    row.append(rmsd)
            rmsd_matrix.append(row)

        clusters = _cluster_poses(rmsd_matrix, ids[:n])

        return {
            "job_uuid": request.job_uuid,
            "ligand_ids": ids[:n],
            "rmsd_matrix": rmsd_matrix,
            "clusters": clusters,
            "n_compared": n,
        }
    except Exception as e:
        logger.error(f"Pose comparison failed: {e}")
        return {"error": str(e), "rmsd_matrix": [], "clusters": []}


@app.post("/filter/admet")
def filter_admet(ligands: List[LigandRecord]):
    """
    Filter ligands by ADMET rules:
      - MW < 500
      - LogP < 5
      - TPSA > 40
      - Num rotatable bonds < 10
    Returns pass/fail per ligand.
    """
    results = []
    for lig in ligands:
        checks = {}
        passed = True

        if lig.mw:
            checks["mw_ok"] = lig.mw < 500
            if not checks["mw_ok"]:
                passed = False
        else:
            checks["mw_ok"] = None

        if lig.logp:
            checks["logp_ok"] = lig.logp < 5
            if not checks["logp_ok"]:
                passed = False
        else:
            checks["logp_ok"] = None

        if lig.tpsa:
            checks["tpsa_ok"] = 40 < lig.tpsa < 140
            if not checks["tpsa_ok"]:
                passed = False
        else:
            checks["tpsa_ok"] = None

        if lig.num_rotatable_bonds is not None:
            checks["rotatable_ok"] = lig.num_rotatable_bonds < 10
            if not checks["rotatable_ok"]:
                passed = False
        else:
            checks["rotatable_ok"] = None

        results.append(
            {
                "ligand_id": lig.ligand_id,
                "passed": passed,
                "checks": checks,
            }
        )

    passed_ids = [r["ligand_id"] for r in results if r["passed"]]
    return {
        "total": len(results),
        "passed": len(passed_ids),
        "failed": len(results) - len(passed_ids),
        "passed_ligands": passed_ids,
        "details": results,
    }


@app.post("/consensus")
def consensus_score(ligands: List[LigandRecord]):
    """
    Compute consensus score: average of all available scoring methods,
    weighted equally. Methods: Vina, GNINA, RF-score, MD stability.
    """
    results = []
    for lig in ligands:
        scores = []
        weights = []

        if lig.vina_score is not None:
            scores.append(-lig.vina_score)
            weights.append(0.35)
        if lig.gnina_score is not None:
            scores.append(-lig.gnina_score)
            weights.append(0.25)
        if lig.rf_score is not None:
            scores.append(-lig.rf_score)
            weights.append(0.20)
        if lig.md_stability is not None:
            scores.append(lig.md_stability)
            weights.append(0.20)

        if not scores:
            consensus = None
        else:
            total_w = sum(weights)
            norm_weights = [w / total_w for w in weights]
            consensus = round(sum(s * w for s, w in zip(scores, norm_weights)), 4)

        results.append(
            {
                "ligand_id": lig.ligand_id,
                "consensus_score": consensus,
                "n_methods": len(scores),
                "methods_used": {
                    "vina": lig.vina_score,
                    "gnina": lig.gnina_score,
                    "rf_score": lig.rf_score,
                    "md_stability": lig.md_stability,
                },
            }
        )

    results.sort(key=lambda x: x["consensus_score"] or -999, reverse=True)
    return {
        "ranked": results,
        "count": len(results),
        "top_ligand": results[0]["ligand_id"] if results else None,
    }


@app.post("/report")
def generate_report(
    job_uuid: str, ligand_ids: List[str], summary: Optional[Dict] = None
):
    """
    Generate a text/JSON analysis report for a set of ligands.
    """
    report_lines = [
        f"Docking Studio Analysis Report",
        f"Job UUID: {job_uuid}",
        f"Generated: {datetime.now().isoformat()}",
        f"",
        f"Total Ligands: {len(ligand_ids)}",
        f"",
    ]

    if summary:
        report_lines.append("Summary:")
        for k, v in summary.items():
            report_lines.append(f"  {k}: {v}")

    report_lines.extend(
        [
            "",
            "---",
            "This report was generated by the Analysis Service (Insight Generator).",
            "Pipeline: Docking → MD Simulation → Consensus Ranking → ADMET Filter",
        ]
    )

    report_text = "\n".join(report_lines)

    return {
        "job_uuid": job_uuid,
        "report": report_text,
        "timestamp": datetime.now().isoformat(),
    }


@app.post("/interactions/summary")
def interactions_summary(interactions: List[Dict[str, Any]]):
    """
    Summarize protein-ligand interactions across multiple ligands.
    Groups by type (H-bond, hydrophobic, pi-pi, etc.) and counts.
    """
    type_counts: Dict[str, int] = {}
    distance_sum: Dict[str, float] = {}
    distance_count: Dict[str, int] = {}

    for interaction in interactions:
        itype = interaction.get("interaction_type", "unknown")
        dist = interaction.get("distance")

        type_counts[itype] = type_counts.get(itype, 0) + 1
        if dist is not None:
            if itype not in distance_sum:
                distance_sum[itype] = 0.0
                distance_count[itype] = 0
            distance_sum[itype] += dist
            distance_count[itype] += 1

    avg_distances = {
        k: round(distance_sum[k] / distance_count[k], 3) for k in distance_sum
    }

    return {
        "total_interactions": len(interactions),
        "by_type": [
            {"type": k, "count": v, "avg_distance_nm": avg_distances.get(k)}
            for k, v in type_counts.items()
        ],
        "most_common": max(type_counts, key=type_counts.get) if type_counts else None,
    }


class ExportTopHitsRequest(BaseModel):
    docking_results: List[Dict[str, Any]]
    top_n: int = 10
    sort_by: str = "vina_score"
    format: str = "csv"


@app.post("/export/top-hits")
def export_top_hits(request: ExportTopHitsRequest):
    """
    Export top N docking hits as CSV or JSON.
    Sort by vina_score (ascending, more negative = better).
    """
    results = request.docking_results

    if request.sort_by == "vina_score":
        results = sorted(
            results, key=lambda x: x.get("vina_score", float("inf")), reverse=False
        )
    elif request.sort_by == "gnina_score" and request.sort_by:
        results = sorted(
            results, key=lambda x: x.get("gnina_score", float("inf")), reverse=False
        )
    elif request.sort_by == "rf_score" and request.sort_by:
        results = sorted(
            results, key=lambda x: x.get("rf_score", float("inf")), reverse=False
        )

    top_hits = results[: request.top_n]

    if request.format == "csv":
        import csv
        import io

        output = io.StringIO()
        if top_hits:
            writer = csv.DictWriter(output, fieldnames=top_hits[0].keys())
            writer.writeheader()
            writer.writerows(top_hits)
        csv_content = output.getvalue()
        return {
            "format": "csv",
            "content": csv_content,
            "filename": f"top_hits_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
            "count": len(top_hits),
        }
    else:
        return {
            "format": "json",
            "results": top_hits,
            "filename": f"top_hits_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
            "count": len(top_hits),
        }


def _normalize(value: float, min_val: float, max_val: float) -> float:
    val = (value - min_val) / (max_val - min_val) if max_val != min_val else 0.5
    return max(0.0, min(1.0, val))


def _admet_score(lig: LigandRecord) -> float:
    score = 1.0
    if lig.mw and lig.mw > 500:
        score -= 0.2
    if lig.logp and lig.logp > 5:
        score -= 0.2
    if lig.tpsa:
        if lig.tpsa < 40 or lig.tpsa > 140:
            score -= 0.15
    if lig.num_rotatable_bonds is not None and lig.num_rotatable_bonds > 10:
        score -= 0.15
    return max(0.0, score)


def _compute_rmsd(mol1, mol2, num_confs: int = 10) -> float:
    try:
        from rdkit.Chem import AllChem
        import numpy as np

        c1 = AllChem.EmbedMolecule(mol1)
        c2 = AllChem.EmbedMolecule(mol2)
        if c1 < 0 or c2 < 0:
            return 2.0

        mmff_c1 = AllChem.MMFFOptimizeMolecule(mol1)
        mmff_c2 = AllChem.MMFFOptimizeMolecule(mol2)

        coords1 = mol1.GetConformer().GetPositions()
        coords2 = mol2.GetConformer().GetPositions()
        rmsd = np.sqrt(np.mean(np.sum((coords1 - coords2) ** 2, axis=1)))
        return round(float(rmsd), 3)
    except Exception:
        return 2.0


def _cluster_poses(
    rmsd_matrix: List[List[float]], ids: List[str], threshold: float = 2.0
) -> List[List[str]]:
    n = len(ids)
    clusters: List[List[int]] = []
    assigned = [False] * n

    for i in range(n):
        if assigned[i]:
            continue
        cluster = [i]
        assigned[i] = True
        for j in range(i + 1, n):
            if not assigned[j] and rmsd_matrix[i][j] < threshold:
                cluster.append(j)
                assigned[j] = True
        clusters.append(cluster)

    return [[ids[i] for i in c] for c in clusters]


@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("=" * 60)
    logger.info("Analysis Service starting up...")
    logger.info("Role: Insight Generator (ranking, consensus, ADMET, comparison)")
    logger.info("=" * 60)
    yield
    logger.info("Analysis Service shutting down...")


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8008)
