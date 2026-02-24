"""
Gnina Engine
GNINA - Deep Learning CNN-based molecular docking.
Supports both GPU and CPU execution (auto-detects).
"""

import logging
import subprocess
import os
import re
from pathlib import Path
from typing import Dict, Any, Optional

from .base import DockingEngine, DockingConfig, DockingResult
from .vina_cpu import VinaCPUEngine

logger = logging.getLogger(__name__)


class GninaEngine(DockingEngine):
    """
    GNINA deep learning docking engine.
    Automatically uses GPU if available, falls back to CPU otherwise.
    """
    
    def __init__(self, docker_manager=None, force_cpu: bool = False):
        self.docker_manager = docker_manager
        self.force_cpu = force_cpu
        self.image_name = "biodockify/gnina:latest"
        self.fallback_engine = VinaCPUEngine(docker_manager)
    
    @property
    def name(self) -> str:
        if self.force_cpu or not self._has_cuda():
            return "GNINA (CPU)"
        return "GNINA (GPU)"
    
    @property
    def requires_gpu(self) -> bool:
        return self.force_cpu is False
    
    @property
    def supports_cnn_scoring(self) -> bool:
        return True
    
    def _has_cuda(self) -> bool:
        """Check if CUDA/GPU is available for GNINA"""
        try:
            result = subprocess.run(
                ["nvidia-smi"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=5
            )
            return result.returncode == 0
        except Exception:
            return False
    
    def run(self, config: DockingConfig) -> DockingResult:
        """Run GNINA docking"""
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
                result = self._run_docker(config)
            else:
                result = self._run_local(config)
            
            if result.status == "FAILED" and "GPU" in result.logs:
                logger.warning("GPU execution failed, falling back to Vina CPU")
                return self.fallback_engine.run(config)
            
            return result
            
        except Exception as e:
            logger.error(f"GNINA docking failed: {e}")
            logger.info("Falling back to Vina CPU engine")
            return self.fallback_engine.run(config)
    
    def _run_docker(self, config: DockingConfig) -> DockingResult:
        """Run using Docker container"""
        result = DockingResult(status="PENDING", job_id=config.job_id)
        
        use_gpu = not self.force_cpu and self._has_cuda()
        
        parameters = {
            "center_x": config.center_x,
            "center_y": config.center_y,
            "center_z": config.center_z,
            "size_x": config.size_x,
            "size_y": config.size_y,
            "size_z": config.size_z,
            "exhaustiveness": config.exhaustiveness,
            "num_modes": config.num_modes,
            "cnn_scoring": "rescore",
            "gpu": use_gpu,
            "cpu": 0 if use_gpu else config.cpu
        }
        
        success = self.docker_manager.start_container(
            receptor_file=config.receptor_file,
            ligand_file=config.ligand_file,
            parameters=parameters,
            job_id=config.job_id,
            image_name=self.image_name,
            use_gpu=use_gpu
        )
        
        if not success:
            result.status = "FAILED"
            result.logs = "Failed to start Docker container"
            return result
        
        result.logs = self.docker_manager.get_container_logs()
        
        output_file = Path(config.output_dir) / f"output_{config.job_id}.sdf"
        if not output_file.exists():
            output_file = Path(config.output_dir) / f"output_{config.job_id}.pdbqt"
        
        if output_file.exists():
            parsed = self.parse_output(str(output_file))
            result.status = "COMPLETED"
            result.binding_energy = parsed.get("binding_energy")
            result.num_modes = parsed.get("num_modes", 0)
            result.poses = parsed.get("poses", [])
            result.gnina_cnn_score = parsed.get("gnina_cnn_score")
            result.gnina_cnn_affinity = parsed.get("gnina_cnn_affinity")
        else:
            result.status = "FAILED"
            result.logs += "\nOutput file not generated"
        
        return result
    
    def _run_local(self, config: DockingConfig) -> DockingResult:
        """Run GNINA locally"""
        result = DockingResult(status="PENDING", job_id=config.job_id)
        
        use_gpu = not self.force_cpu and self._has_cuda()
        
        try:
            cmd = [
                "gnina",
                "-r", config.receptor_file,
                "-l", config.ligand_file,
                "--center_x", str(config.center_x),
                "--center_y", str(config.center_y),
                "--center_z", str(config.center_z),
                "--size_x", str(config.size_x),
                "--size_y", str(config.size_y),
                "--size_z", str(config.size_z),
                "--num_modes", str(config.num_modes),
                "--cnn_scoring", "rescore",
                "-o", os.path.join(config.output_dir, f"output_{config.job_id}.sdf")
            ]
            
            if not use_gpu:
                cmd.extend(["--cpu", str(config.cpu or 8)])
            
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=1200
            )
            
            result.logs = proc.stdout + proc.stderr
            
            if proc.returncode == 0:
                output_file = os.path.join(config.output_dir, f"output_{config.job_id}.sdf")
                if os.path.exists(output_file):
                    parsed = self.parse_output(output_file)
                    result.status = "COMPLETED"
                    result.binding_energy = parsed.get("binding_energy")
                    result.num_modes = parsed.get("num_modes", 0)
                    result.poses = parsed.get("poses", [])
                    result.gnina_cnn_score = parsed.get("gnina_cnn_score")
                    result.gnina_cnn_affinity = parsed.get("gnina_cnn_affinity")
                else:
                    result.status = "FAILED"
            else:
                result.status = "FAILED"
                
        except subprocess.TimeoutExpired:
            result.status = "TIMEOUT"
            result.logs = "Docking timed out"
        except FileNotFoundError:
            result.status = "FAILED"
            result.logs = "gnina executable not found"
        
        return result
    
    def parse_output(self, output_file: str) -> Dict[str, Any]:
        """Parse GNINA output (SDF or PDBQT with CNN scores)"""
        poses = []
        energies = []
        cnn_scores = []
        cnn_affinities = []
        
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            if output_file.endswith('.sdf'):
                poses = self._parse_sdf(content)
                for pose in poses:
                    if pose.get("binding_energy"):
                        energies.append(pose["binding_energy"])
                    if pose.get("gnina_cnn_score"):
                        cnn_scores.append(pose["gnina_cnn_score"])
                    if pose.get("gnina_cnn_affinity"):
                        cnn_affinities.append(pose["gnina_cnn_affinity"])
            else:
                poses = self._parse_pdbqt(content)
                for pose in poses:
                    if pose.get("binding_energy"):
                        energies.append(pose["binding_energy"])
            
            binding_energy = min(energies) if energies else None
            cnn_score = max(cnn_scores) if cnn_scores else None
            cnn_affinity = max(cnn_affinities) if cnn_affinities else None
            
            return {
                "binding_energy": binding_energy,
                "num_modes": len(poses),
                "poses": poses,
                "gnina_cnn_score": cnn_score,
                "gnina_cnn_affinity": cnn_affinity
            }
            
        except Exception as e:
            logger.error(f"Failed to parse GNINA output: {e}")
            return {"binding_energy": None, "num_modes": 0, "poses": []}
    
    def _parse_sdf(self, content: str) -> list:
        """Parse SDF format with CNN scores"""
        poses = []
        records = content.split("$$$$")
        
        for idx, record in enumerate(records):
            if not record.strip():
                continue
            
            lines = record.strip().split('\n')
            pose = {
                "mode": idx + 1,
                "binding_energy": None,
                "gnina_cnn_score": None,
                "gnina_cnn_affinity": None
            }
            
            for line in lines:
                line = line.strip()
                
                if line.startswith("> <CNNscore>"):
                    continue
                if line.startswith("> <CNNaffinity>"):
                    continue
                
                try:
                    if line.replace('.', '').replace('-', '').isdigit():
                        val = float(line)
                        if -20 <= val <= 0 and pose["binding_energy"] is None:
                            pose["binding_energy"] = val
                        elif 0 <= val <= 1:
                            pose["gnina_cnn_score"] = val
                        elif 0 <= val <= 15:
                            pose["gnina_cnn_affinity"] = val
                except ValueError:
                    continue
            
            if pose["binding_energy"] is not None:
                poses.append(pose)
        
        return poses
    
    def _parse_pdbqt(self, content: str) -> list:
        """Parse PDBQT format"""
        poses = []
        
        for line in content.split('\n'):
            if "REMARK VINA RESULT" in line:
                parts = line.split()
                try:
                    energy = float(parts[3])
                    poses.append({
                        "mode": len(poses) + 1,
                        "binding_energy": energy
                    })
                except (ValueError, IndexError):
                    pass
        
        return poses
