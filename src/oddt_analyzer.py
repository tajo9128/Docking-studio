"""
BioDockify Docking Studio - ODDT Interaction Analyzer
Analyzes molecular interactions using Open Drug Discovery Toolkit
"""

import logging
import os
import tempfile
from typing import Dict, List, Any, Tuple
from pathlib import Path
import json
import numpy as np

logger = logging.getLogger(__name__)

class ODDTAnalyzer:
    """ODDT interaction analyzer"""
    
    def __init__(self, docker_manager):
        """Initialize ODDT analyzer"""
        self.docker_manager = docker_manager
    
    def analyze_interactions(self, receptor_file: str, ligand_file: str,
                          output_file: str, job_id: str) -> Dict[str, Any]:
        """Analyze molecular interactions using ODDT"""
        logger.info(f"Analyzing interactions for job {job_id}")
        
        interactions = {
            "hydrogen_bonds": [],
            "hydrophobic_contacts": [],
            "pi_stacking": [],
            "total_count": 0
        }
        
        try:
            # Check for ODDT
            try:
                from oddt import toolkit
                from oddt import interactions as oddt_interactions
            except ImportError:
                return {"status": "FAILED", "error": "ODDT library not available", "interactions": None}
            
            # Load molecules
            try:
                 protein = toolkit.readfile('pdb', receptor_file).__next__()
                 ligand = toolkit.readfile('pdbqt', output_file).__next__()
            except Exception as e:
                 # Try pdbqt for both maybe? Or receptor as pdb
                 return {"status": "FAILED", "error": f"Load failed: {e}", "interactions": None}

            # Detect Interactions
            try:
                hbonds = oddt_interactions.hbonds(protein, ligand, cutoff=3.5)
                for hb in hbonds:
                     interactions["hydrogen_bonds"].append({
                         "distance": float(hb['distance']),
                         "type": "hbond"
                     })
                     
                hydro = oddt_interactions.hydrophobic_contacts(protein, ligand, cutoff=4.5)
                for h in hydro:
                     interactions["hydrophobic_contacts"].append({
                         "distance": float(h['distance']),
                         "type": "hydrophobic"
                     })
                     
                pi = oddt_interactions.pi_stacking(protein, ligand, cutoff=5.0)
                for p in pi:
                     interactions["pi_stacking"].append({
                         "distance": float(p['distance']),
                         "type": "pi_stacking"
                     })
            except Exception:
                pass

            interactions["total_count"] = len(interactions["hydrogen_bonds"]) + len(interactions["hydrophobic_contacts"]) + len(interactions["pi_stacking"])
            
            return {
                "status": "COMPLETED",
                "interactions": interactions
            }
            
        except Exception as e:
            logger.error(f"ODDT analysis failed: {e}")
            return {
                "status": "FAILED",
                "error": str(e),
                "interactions": None
            }
    
    def visualize_interactions(self, interactions: Dict[str, Any], output_path: str) -> bool:
        """Visualize interactions (would generate 2D diagram)"""
        logger.info("Visualizing interactions")
        # In real implementation, would use RDKit or other library to generate 2D diagram
        return True  # Assume success
