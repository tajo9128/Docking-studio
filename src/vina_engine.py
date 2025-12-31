"""
BioDockify Docking Studio - AutoDock Vina Engine Wrapper
Executes AutoDock Vina molecular docking simulations
"""

import logging
import subprocess
import os
from pathlib import Path
from typing import Dict, Optional, Any
import json
import re

logger = logging.getLogger(__name__)

class VinaEngine:
    """AutoDock Vina engine wrapper"""
    
    def __init__(self, docker_manager):
        """Initialize Vina engine"""
        self.docker_manager = docker_manager
    
    def run_docking(self, receptor_file: str, ligand_file: str,
                   parameters: Dict[str, Any], job_id: str) -> Dict[str, Any]:
        """Run AutoDock Vina docking"""
        logger.info(f"Running AutoDock Vina docking for job {job_id}")
        
        results = {
            "status": "PENDING",
            "binding_energy": None,
            "num_modes": 0,
            "poses": [],
            "logs": "",
            "job_id": job_id
        }
        
        try:
            # Prepare receptor and ligand paths
            receptor_path = Path(receptor_file)
            ligand_path = Path(ligand_file)
            
            # Load template config if available (Fix #21)
            config = {}
            template_path = Path(__file__).parent / "templates" / "dock_vina.conf"
            if template_path.exists():
                try:
                    import configparser
                    parser = configparser.ConfigParser()
                    parser.read(template_path)
                    for key in ['exhaustiveness', 'num_modes', 'energy_range', 'cpu']:
                        if 'DEFAULT' in parser and key in parser['DEFAULT']:
                            config[key] = parser.getint('DEFAULT', key)
                except Exception as e:
                    logger.warning(f"Failed to parse Vina template: {e}")

            # Merge config into parameters if not present
            for k, v in config.items():
                if k not in parameters:
                    parameters[k] = v

            # Start Docker container
            success = self.docker_manager.start_container(
                receptor_file=str(receptor_path),
                ligand_file=str(ligand_path),
                parameters=parameters,
                job_id=job_id
            )
            
            if not success:
                results["status"] = "FAILED"
                results["logs"] = "Failed to start Docker container or timed out"
                return results
            
            # Get logs
            logs = self.docker_manager.get_container_logs()
            results["logs"] = logs
            
            # Check for output (Fix Logic)
            # In docker_manager we mount output to /data/output, mapped to local 'data'
            output_file = Path("data") / f"output_{job_id}.pdbqt"
            
            if output_file.exists():
                logger.info(f"Output file found: {output_file}")
                parsed = self.parse_output(str(output_file))
                results.update(parsed)
                results["status"] = "COMPLETED"
            else:
                logger.error(f"Output file NOT found: {output_file}")
                results["status"] = "FAILED"
                results["logs"] += "\nError: Output file was not generated."
            
            logger.info(f"Vina docking finished for job {job_id}")
            return results
            
        except Exception as e:
            logger.error(f"Vina docking failed: {e}")
            results["status"] = "FAILED"
            results["logs"] = f"Vina error: {str(e)}"
            return results
    
    def parse_output(self, output_file: str) -> Dict[str, Any]:
        """Parse Vina output file"""
        energies = []
        poses = []
        
        try:
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            # Parse output
            for i, line in enumerate(lines):
                 if "REMARK VINA RESULT" in line or (line.startswith("MODE") and "RE" not in line): # Simple heuristics
                      pass 
                 # Look for Vina Result line: RE  -9.2  0.0 0.0 ... (This explanation varies, usually it's REMARK VINA RESULT:    -6.1      0.000      0.000)
                 # Or in standard output table: 
                 # mode |   affinity | dist from best mode
                 #      | (kcal/mol) | rmsd l.b.| rmsd u.b.
                 # -----+------------+----------+----------
                 #    1       -6.1      0.000      0.000
                 
                 # Regex for table lines:
                 #    1       -6.1      0.000      0.000
                 match = re.search(r'^\s*(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)', line)
                 if match:
                      mode = int(match.group(1))
                      energy = float(match.group(2))
                      rmsd_lb = float(match.group(3))
                      rmsd_ub = float(match.group(4))
                      
                      energies.append(energy)
                      poses.append({
                          "mode": mode,
                          "binding_energy": energy,
                          "rmsd_lb": rmsd_lb,
                          "rmsd_ub": rmsd_ub
                      })

            if not energies:
                 # Fallback for REMARK VINA RESULT lines in PDBQT
                 for line in lines:
                      if "REMARK VINA RESULT" in line:
                           # REMARK VINA RESULT:    -6.1      0.000      0.000
                           parts = line.split()
                           if len(parts) >= 4:
                                energy = float(parts[3])
                                energies.append(energy)
                                poses.append({
                                     "mode": len(energies),
                                     "binding_energy": energy, 
                                     "rmsd_lb": 0, "rmsd_ub": 0 # simplified
                                })

            binding_energy = min(energies) if energies else None
            
            return {
                "binding_energy": binding_energy,
                "num_modes": len(poses),
                "poses": poses
            }
            
        except Exception as e:
            logger.error(f"Failed to parse Vina output: {e}")
            return {"binding_energy": None, "num_modes": 0, "poses": []}
