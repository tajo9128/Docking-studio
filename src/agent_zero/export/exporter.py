"""
Result Exporter - PDF/CSV/JSON Export
Exports docking results in multiple formats.
"""

import json
import csv
import logging
import os
from typing import Dict, Any, List, Optional
from pathlib import Path
from datetime import datetime

logger = logging.getLogger(__name__)


class ResultExporter:
    """
    Exports docking results to PDF, CSV, and JSON formats.
    """
    
    def __init__(self):
        self.export_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    def export_all(
        self,
        job_id: str,
        output_dir: str,
        poses: List[Dict[str, Any]],
        best_pose: Optional[Dict[str, Any]],
        **results
    ) -> Dict[str, str]:
        """
        Export all result formats.
        
        Returns:
            dict: Paths to exported files
        """
        exported = {}
        
        base_path = Path(output_dir) / f"results_{job_id}"
        base_path.parent.mkdir(parents=True, exist_ok=True)
        
        json_path = self.export_json(
            poses=poses,
            best_pose=best_pose,
            job_id=job_id,
            **results
        )
        if json_path:
            exported["json"] = str(json_path)
        
        csv_path = self.export_csv(
            poses=poses,
            job_id=job_id
        )
        if csv_path:
            exported["csv"] = str(csv_path)
        
        pdf_path = self.export_pdf(
            job_id=job_id,
            output_dir=str(base_path.parent),
            poses=poses,
            best_pose=best_pose,
            **results
        )
        if pdf_path:
            exported["pdf"] = str(pdf_path)
        
        return exported
    
    def export_json(
        self,
        poses: List[Dict[str, Any]],
        best_pose: Optional[Dict[str, Any]],
        job_id: str,
        **results
    ) -> Optional[Path]:
        """Export results to JSON"""
        try:
            data = {
                "job_id": job_id,
                "timestamp": self.export_timestamp,
                "total_poses": len(poses),
                "best_pose": best_pose,
                "all_poses": poses,
                "additional_results": {
                    k: v for k, v in results.items() if v is not None
                }
            }
            
            output_path = Path("data") / f"results_{job_id}.json"
            with open(output_path, 'w') as f:
                json.dump(data, f, indent=2)
            
            logger.info(f"JSON export: {output_path}")
            return output_path
            
        except Exception as e:
            logger.error(f"JSON export failed: {e}")
            return None
    
    def export_csv(self, poses: List[Dict[str, Any]], job_id: str) -> Optional[Path]:
        """Export results to CSV"""
        try:
            if not poses:
                return None
            
            output_path = Path("data") / f"results_{job_id}.csv"
            
            fieldnames = [
                "rank", "mode", "binding_energy", "rmsd_lb", "rmsd_ub",
                "gnina_cnn_score", "gnina_cnn_affinity",
                "rf_predicted_pKd", "consensus_score"
            ]
            
            with open(output_path, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                
                for idx, pose in enumerate(poses, 1):
                    row = {
                        "rank": idx,
                        "mode": pose.get("mode", idx),
                        "binding_energy": pose.get("binding_energy"),
                        "rmsd_lb": pose.get("rmsd_lb"),
                        "rmsd_ub": pose.get("rmsd_ub"),
                        "gnina_cnn_score": pose.get("gnina_cnn_score"),
                        "gnina_cnn_affinity": pose.get("gnina_cnn_affinity"),
                        "rf_predicted_pKd": pose.get("rf_predicted_pKd"),
                        "consensus_score": pose.get("consensus_score")
                    }
                    writer.writerow(row)
            
            logger.info(f"CSV export: {output_path}")
            return output_path
            
        except Exception as e:
            logger.error(f"CSV export failed: {e}")
            return None
    
    def export_pdf(
        self,
        job_id: str,
        output_dir: str,
        poses: List[Dict[str, Any]],
        best_pose: Optional[Dict[str, Any]],
        **results
    ) -> Optional[Path]:
        """Export results to PDF report"""
        try:
            from src.agent_zero.export.report import PDFReportGenerator
            
            generator = PDFReportGenerator()
            output_path = generator.generate(
                job_id=job_id,
                output_dir=output_dir,
                poses=poses,
                best_pose=best_pose,
                **results
            )
            
            logger.info(f"PDF export: {output_path}")
            return output_path
            
        except ImportError as e:
            logger.warning(f"PDF generation requires reportlab: {e}")
            return self._generate_simple_pdf(job_id, output_dir, poses, best_pose)
        except Exception as e:
            logger.error(f"PDF export failed: {e}")
            return None
    
    def _generate_simple_pdf(
        self,
        job_id: str,
        output_dir: str,
        poses: List[Dict],
        best_pose: Optional[Dict]
    ) -> Optional[Path]:
        """Generate simple text-based PDF fallback"""
        try:
            output_path = Path(output_dir) / f"results_{job_id}.txt"
            
            with open(output_path, 'w') as f:
                f.write("=" * 60 + "\n")
                f.write("BioDockify Docking Results\n")
                f.write("=" * 60 + "\n\n")
                f.write(f"Job ID: {job_id}\n")
                f.write(f"Generated: {self.export_timestamp}\n\n")
                
                f.write("BEST POSE\n")
                f.write("-" * 40 + "\n")
                if best_pose:
                    f.write(f"Mode: {best_pose.get('mode')}\n")
                    f.write(f"Binding Energy: {best_pose.get('binding_energy')} kcal/mol\n")
                    f.write(f"Consensus Score: {best_pose.get('consensus_score')}\n")
                f.write("\n")
                
                f.write("ALL POSES (Ranked)\n")
                f.write("-" * 40 + "\n")
                for idx, pose in enumerate(poses[:10], 1):
                    f.write(f"{idx}. Mode {pose.get('mode')}: {pose.get('binding_energy')} kcal/mol\n")
            
            return output_path
            
        except Exception as e:
            logger.error(f"Simple PDF export failed: {e}")
            return None
