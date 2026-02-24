"""
Reporter - Execution and Scientific Report Generation
Generates structured reports for docking jobs
"""

import os
import json
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


@dataclass
class ExecutionReport:
    """Execution summary report"""
    job_id: str
    engine_used: str
    gpu_detected: bool
    gpu_name: str = ""
    retry_count: int = 0
    runtime_seconds: float = 0.0
    cpu_usage_avg: float = 0.0
    gpu_memory_peak_mb: float = 0.0
    failures_encountered: List[str] = field(default_factory=list)
    fixes_applied: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


@dataclass
class ScientificReport:
    """Scientific results report"""
    job_id: str
    receptor_file: str
    ligand_count: int
    poses_generated: int
    
    # Scores
    vina_scores: List[float] = field(default_factory=list)
    gnina_scores: List[float] = field(default_factory=list)
    consensus_scores: List[float] = field(default_factory=list)
    
    # Best pose
    best_pose_id: int = 0
    best_vina_score: float = 0.0
    best_gnina_score: float = 0.0
    best_consensus_score: float = 0.0
    
    # Interactions
    hbond_count: int = 0
    hydrophobic_count: int = 0
    pi_stacking_count: int = 0


@dataclass
class DiagnosticReport:
    """Diagnostic report"""
    job_id: str
    failures: List[Dict] = field(default_factory=list)
    retry_history: List[Dict] = field(default_factory=list)
    resource_warnings: List[str] = field(default_factory=list)
    performance_warnings: List[str] = field(default_factory=list)


class Reporter:
    """
    Generates comprehensive reports for docking jobs.
    """
    
    def __init__(self, output_dir: str = "results"):
        """Initialize reporter"""
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        logger.info(f"Reporter initialized: {output_dir}")
    
    def generate_all_reports(
        self,
        job_id: str,
        job_data: Dict,
        execution_data: Optional[Dict] = None,
        scientific_data: Optional[Dict] = None,
        diagnostic_data: Optional[Dict] = None
    ) -> Dict[str, str]:
        """
        Generate all reports for a job.
        
        Args:
            job_id: Job identifier
            job_data: Job configuration and results
            execution_data: Execution metrics
            scientific_data: Scientific results
            diagnostic_data: Diagnostic information
            
        Returns:
            Dict of report file paths
        """
        reports = {}
        job_dir = os.path.join(self.output_dir, job_id)
        os.makedirs(job_dir, exist_ok=True)
        
        # Generate execution report
        if execution_data:
            exec_report = self._build_execution_report(job_id, execution_data)
            exec_path = self._save_json(exec_report, job_dir, "execution_report.json")
            reports["execution"] = exec_path
        
        # Generate scientific report
        if scientific_data:
            sci_report = self._build_scientific_report(job_id, scientific_data)
            sci_path = self._save_json(sci_report, job_dir, "scientific_report.json")
            reports["scientific"] = sci_path
        
        # Generate diagnostic report
        if diagnostic_data:
            diag_report = self._build_diagnostic_report(job_id, diagnostic_data)
            diag_path = self._save_json(diag_report, job_dir, "diagnostic_report.json")
            reports["diagnostic"] = diag_path
        
        # Generate summary markdown
        summary_path = self._generate_summary_markdown(
            job_id, execution_data, scientific_data, job_dir
        )
        reports["summary"] = summary_path
        
        logger.info(f"Generated {len(reports)} reports for job {job_id}")
        
        return reports
    
    def _build_execution_report(self, job_id: str, data: Dict) -> Dict:
        """Build execution report"""
        return {
            "job_id": job_id,
            "generated_at": datetime.now().isoformat(),
            "engine": data.get("engine", "unknown"),
            "gpu_detected": data.get("gpu_detected", False),
            "gpu_name": data.get("gpu_name", ""),
            "retry_count": data.get("retry_count", 0),
            "runtime_seconds": data.get("runtime_seconds", 0.0),
            "cpu_usage_avg": data.get("cpu_usage_avg", 0.0),
            "gpu_memory_peak_mb": data.get("gpu_memory_peak_mb", 0.0),
            "failures_encountered": data.get("failures", []),
            "fixes_applied": data.get("fixes", []),
            "warnings": data.get("warnings", []),
            "status": data.get("status", "unknown"),
        }
    
    def _build_scientific_report(self, job_id: str, data: Dict) -> Dict:
        """Build scientific report"""
        return {
            "job_id": job_id,
            "generated_at": datetime.now().isoformat(),
            "receptor": data.get("receptor_file", ""),
            "ligand_count": data.get("ligand_count", 0),
            "poses_generated": data.get("poses_generated", 0),
            "best_pose": {
                "pose_id": data.get("best_pose_id", 0),
                "vina_score": data.get("best_vina_score", 0.0),
                "gnina_score": data.get("best_gnina_score", 0.0),
                "consensus_score": data.get("best_consensus_score", 0.0),
            },
            "score_distribution": {
                "vina": data.get("vina_scores", []),
                "gnina": data.get("gnina_scores", []),
                "consensus": data.get("consensus_scores", []),
            },
            "interactions": {
                "hydrogen_bonds": data.get("hbond_count", 0),
                "hydrophobic": data.get("hydrophobic_count", 0),
                "pi_stacking": data.get("pi_stacking_count", 0),
            },
            "ranking": data.get("ranking", [])[:10],  # Top 10
        }
    
    def _build_diagnostic_report(self, job_id: str, data: Dict) -> Dict:
        """Build diagnostic report"""
        return {
            "job_id": job_id,
            "generated_at": datetime.now().isoformat(),
            "failures": data.get("failures", []),
            "retry_history": data.get("retry_history", []),
            "resource_warnings": data.get("resource_warnings", []),
            "performance_warnings": data.get("performance_warnings", []),
        }
    
    def _save_json(self, data: Dict, job_dir: str, filename: str) -> str:
        """Save report as JSON"""
        filepath = os.path.join(job_dir, filename)
        
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2)
        
        return filepath
    
    def _generate_summary_markdown(
        self,
        job_id: str,
        execution_data: Optional[Dict],
        scientific_data: Optional[Dict],
        job_dir: str
    ) -> str:
        """Generate summary markdown report"""
        lines = [
            f"# Docking Results - {job_id}",
            f"",
            f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "",
            "## Execution Summary",
            "",
        ]
        
        if execution_data:
            lines.extend([
                f"- **Engine:** {execution_data.get('engine', 'N/A')}",
                f"- **GPU:** {'Yes' if execution_data.get('gpu_detected') else 'No'}",
                f"- **Runtime:** {execution_data.get('runtime_seconds', 0):.1f}s",
                f"- **Retries:** {execution_data.get('retry_count', 0)}",
                f"- **Status:** {execution_data.get('status', 'unknown')}",
                "",
            ])
        
        if scientific_data:
            lines.extend([
                "## Scientific Results",
                "",
                f"- **Ligands:** {scientific_data.get('ligand_count', 0)}",
                f"- **Poses:** {scientific_data.get('poses_generated', 0)}",
                "",
                "### Best Pose",
                "",
                f"- **VINA Score:** {scientific_data.get('best_vina_score', 0):.2f} kcal/mol",
                f"- **GNINA Score:** {scientific_data.get('best_gnina_score', 0):.2f}",
                f"- **Consensus:** {scientific_data.get('best_consensus_score', 0):.2f}",
                "",
                "### Interactions",
                "",
                f"- **H-bonds:** {scientific_data.get('hbond_count', 0)}",
                f"- **Hydrophobic:** {scientific_data.get('hydrophobic_count', 0)}",
                f"- **π-stacking:** {scientific_data.get('pi_stacking_count', 0)}",
                "",
            ])
        
        lines.extend([
            "## Files",
            "",
            "- `poses/` - Docked pose PDB files",
            "- `scientific_report.json` - Detailed scores",
            "- `execution_report.json` - Execution metrics",
            "- `diagnostic_report.json` - Diagnostics",
            "",
            "*Generated by BioDockify Docking Studio*",
        ])
        
        filepath = os.path.join(job_dir, "README.md")
        with open(filepath, 'w') as f:
            f.write('\n'.join(lines))
        
        return filepath
    
    def generate_summary_html(
        self,
        job_id: str,
        data: Dict
    ) -> str:
        """Generate HTML summary for GUI display"""
        html = f"""
        <div style="font-family: Arial, sans-serif; padding: 20px;">
            <h2>Job: {job_id}</h2>
            <table style="width: 100%; border-collapse: collapse;">
                <tr>
                    <td style="padding: 8px; border: 1px solid #ddd;"><b>Engine</b></td>
                    <td style="padding: 8px; border: 1px solid #ddd;">{data.get('engine', 'N/A')}</td>
                </tr>
                <tr>
                    <td style="padding: 8px; border: 1px solid #ddd;"><b>Status</b></td>
                    <td style="padding: 8px; border: 1px solid #ddd;">{data.get('status', 'N/A')}</td>
                </tr>
                <tr>
                    <td style="padding: 8px; border: 1px solid #ddd;"><b>Runtime</b></td>
                    <td style="padding: 8px; border: 1px solid #ddd;">{data.get('runtime', 0):.1f}s</td>
                </tr>
                <tr>
                    <td style="padding: 8px; border: 1px solid #ddd;"><b>Poses</b></td>
                    <td style="padding: 8px; border: 1px solid #ddd;">{data.get('poses', 0)}</td>
                </tr>
                <tr>
                    <td style="padding: 8px; border: 1px solid #ddd;"><b>Best Score</b></td>
                    <td style="padding: 8px; border: 1px solid #ddd;">{data.get('best_score', 0):.2f}</td>
                </tr>
            </table>
        </div>
        """
        return html
