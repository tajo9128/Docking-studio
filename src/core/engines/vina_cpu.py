"""
AutoDock Vina CPU Engine
Runs Vina docking using Docker container on CPU.
"""

import logging
import subprocess
import os
import re
from pathlib import Path
from typing import Dict, Any, List, Optional

from .base import DockingEngine, DockingConfig, DockingResult

logger = logging.getLogger(__name__)


class VinaCPUEngine(DockingEngine):
    """AutoDock Vina CPU implementation"""
    
    def __init__(self, docker_manager=None):
        self.docker_manager = docker_manager
        self.image_name = "biodockify/vina-cpu:latest"
    
    @property
    def name(self) -> str:
        return "AutoDock Vina CPU"
    
    @property
    def requires_gpu(self) -> bool:
        return False
    
    @property
    def supports_cnn_scoring(self) -> bool:
        return False
    
    def run(self, config: DockingConfig) -> DockingResult:
        """Run Vina docking on CPU"""
        logger.info(f"Running {self.name} for job {config.job_id}")
        
        result = DockingResult(
            status="PENDING",
            job_id=config.job_id
        )
        
        errors = self.validate_config(config)
        if errors:
            result.status = "FAILED"
            result.logs = "\n".join(errors)
            return result
        
        try:
            if self.docker_manager:
                return self._run_docker(config)
            else:
                return self._run_local(config)
                
        except Exception as e:
            logger.error(f"Vina CPU docking failed: {e}")
            result.status = "FAILED"
            result.logs = str(e)
            return result
    
    def _run_docker(self, config: DockingConfig) -> DockingResult:
        """Run using Docker container"""
        result = DockingResult(status="PENDING", job_id=config.job_id)
        
        parameters = {
            "center_x": config.center_x,
            "center_y": config.center_y,
            "center_z": config.center_z,
            "size_x": config.size_x,
            "size_y": config.size_y,
            "size_z": config.size_z,
            "exhaustiveness": config.exhaustiveness,
            "num_modes": config.num_modes,
            "energy_range": config.energy_range,
            "cpu": config.cpu
        }
        
        success = self.docker_manager.start_container(
            receptor_file=config.receptor_file,
            ligand_file=config.ligand_file,
            parameters=parameters,
            job_id=config.job_id,
            image_name=self.image_name
        )
        
        if not success:
            result.status = "FAILED"
            result.logs = "Failed to start Docker container"
            return result
        
        result.logs = self.docker_manager.get_container_logs()
        
        output_file = Path(config.output_dir) / f"output_{config.job_id}.pdbqt"
        if output_file.exists():
            parsed = self.parse_output(str(output_file))
            result.status = "COMPLETED"
            result.binding_energy = parsed.get("binding_energy")
            result.num_modes = parsed.get("num_modes", 0)
            result.poses = parsed.get("poses", [])
        else:
            result.status = "FAILED"
            result.logs += "\nOutput file not generated"
        
        return result
    
    def _run_local(self, config: DockingConfig) -> DockingResult:
        """Run Vina locally (requires vina installed)"""
        result = DockingResult(status="PENDING", job_id=config.job_id)
        
        config_file = self._create_config_file(config)
        
        try:
            cmd = [
                "vina",
                "--config", config_file,
                "--receptor", config.receptor_file,
                "--ligand", config.ligand_file,
                "--out", os.path.join(config.output_dir, f"output_{config.job_id}.pdbqt"),
                "--log", os.path.join(config.output_dir, f"vina_log_{config.job_id}.txt")
            ]
            
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600
            )
            
            result.logs = proc.stdout + proc.stderr
            
            if proc.returncode == 0:
                output_file = os.path.join(config.output_dir, f"output_{config.job_id}.pdbqt")
                if os.path.exists(output_file):
                    parsed = self.parse_output(output_file)
                    result.status = "COMPLETED"
                    result.binding_energy = parsed.get("binding_energy")
                    result.num_modes = parsed.get("num_modes", 0)
                    result.poses = parsed.get("poses", [])
                else:
                    result.status = "FAILED"
            else:
                result.status = "FAILED"
                
        except subprocess.TimeoutExpired:
            result.status = "TIMEOUT"
            result.logs = "Docking timed out"
        except FileNotFoundError:
            result.status = "FAILED"
            result.logs = "Vina executable not found. Install vina or use Docker mode."
        
        return result
    
    def _create_config_file(self, config: DockingConfig) -> str:
        """Create Vina config file"""
        config_path = os.path.join(config.output_dir, f"vina_config_{config.job_id}.txt")
        
        with open(config_path, 'w') as f:
            f.write(f"center_x = {config.center_x}\n")
            f.write(f"center_y = {config.center_y}\n")
            f.write(f"center_z = {config.center_z}\n")
            f.write(f"size_x = {config.size_x}\n")
            f.write(f"size_y = {config.size_y}\n")
            f.write(f"size_z = {config.size_z}\n")
            f.write(f"exhaustiveness = {config.exhaustiveness}\n")
            f.write(f"num_modes = {config.num_modes}\n")
            f.write(f"energy_range = {config.energy_range}\n")
            if config.cpu > 0:
                f.write(f"cpu = {config.cpu}\n")
        
        return config_path
    
    def parse_output(self, output_file: str) -> Dict[str, Any]:
        """Parse Vina output file"""
        energies = []
        poses = []
        
        try:
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            for line in lines:
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
                for line in lines:
                    if "REMARK VINA RESULT" in line:
                        parts = line.split()
                        if len(parts) >= 4:
                            try:
                                energy = float(parts[3])
                                energies.append(energy)
                                poses.append({
                                    "mode": len(energies),
                                    "binding_energy": energy,
                                    "rmsd_lb": 0,
                                    "rmsd_ub": 0
                                })
                            except (ValueError, IndexError):
                                pass
            
            binding_energy = min(energies) if energies else None
            
            return {
                "binding_energy": binding_energy,
                "num_modes": len(poses),
                "poses": poses
            }
            
        except Exception as e:
            logger.error(f"Failed to parse Vina output: {e}")
            return {"binding_energy": None, "num_modes": 0, "poses": []}
