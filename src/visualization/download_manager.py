"""
Download Manager
Handles downloading and exporting docking results.
"""

import os
import shutil
import zipfile
import json
import logging
from typing import List, Dict, Optional, Callable
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class ResultFile:
    """A single result file"""
    filename: str
    path: str
    size: int
    file_type: str  # pdbqt, csv, json, pdf, png, sdf


@dataclass
class DockingResult:
    """Complete docking result package"""
    job_id: str
    receptor: str
    ligands: List[str]
    poses: List[Dict]  # List of pose data with score, pdb content
    interactions: List[Dict]
    summary: Dict
    timestamp: str
    output_dir: str


class DownloadManager:
    """
    Manages downloading and exporting of docking results.
    
    Supports:
    - Individual file download
    - CSV ranking export
    - JSON results export
    - PDF report export
    - ZIP bundle export
    - Snapshot images
    """
    
    def __init__(self, base_output_dir: str = "results"):
        """
        Initialize download manager.
        
        Args:
            base_output_dir: Base directory for storing results
        """
        self.base_output_dir = base_output_dir
        os.makedirs(base_output_dir, exist_ok=True)
        
        logger.info(f"DownloadManager initialized: {base_output_dir}")
    
    def create_job_directory(self, job_id: str) -> str:
        """Create directory for a job's results"""
        job_dir = os.path.join(self.base_output_dir, job_id)
        os.makedirs(job_dir, exist_ok=True)
        
        # Create subdirectories
        os.makedirs(os.path.join(job_dir, "poses"), exist_ok=True)
        os.makedirs(os.path.join(job_dir, "snapshots"), exist_ok=True)
        
        return job_dir
    
    def save_pose(self, job_id: str, pose_id: int, pdb_content: str, score: float) -> str:
        """Save a docking pose PDB file"""
        job_dir = self.create_job_directory(job_id)
        pose_file = os.path.join(job_dir, "poses", f"pose_{pose_id}.pdb")
        
        with open(pose_file, 'w') as f:
            f.write(pdb_content)
        
        logger.info(f"Saved pose {pose_id} to {pose_file}")
        return pose_file
    
    def save_poses(self, job_id: str, poses: List[Dict]) -> List[str]:
        """Save all docking poses"""
        pose_files = []
        for pose in poses:
            pose_file = self.save_pose(
                job_id,
                pose.get('id', 0),
                pose.get('pdb', ''),
                pose.get('score', 0.0)
            )
            pose_files.append(pose_file)
        
        return pose_files
    
    def export_csv_ranking(
        self,
        job_id: str,
        poses: List[Dict],
        output_filename: Optional[str] = None
    ) -> str:
        """
        Export poses as CSV ranking.
        
        Args:
            job_id: Job identifier
            poses: List of pose dictionaries
            output_filename: Optional custom filename
            
        Returns:
            Path to exported CSV file
        """
        import csv
        
        job_dir = self.create_job_directory(job_id)
        
        if output_filename is None:
            csv_file = os.path.join(job_dir, "ranking.csv")
        else:
            csv_file = os.path.join(job_dir, output_filename)
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Header
            writer.writerow([
                'Pose_ID', 'Score', 'RMSD', 'Binding_Energy',
                'Ligand_Name', 'Interactions', 'Notes'
            ])
            
            # Data rows
            for pose in poses:
                writer.writerow([
                    pose.get('id', ''),
                    pose.get('score', ''),
                    pose.get('rmsd', 'N/A'),
                    pose.get('energy', 'N/A'),
                    pose.get('ligand_name', ''),
                    pose.get('interaction_count', ''),
                    pose.get('notes', '')
                ])
        
        logger.info(f"Exported CSV ranking: {csv_file}")
        return csv_file
    
    def export_json_results(
        self,
        job_id: str,
        data: Dict,
        output_filename: Optional[str] = None
    ) -> str:
        """
        Export complete results as JSON.
        
        Args:
            job_id: Job identifier
            data: Result data dictionary
            output_filename: Optional custom filename
            
        Returns:
            Path to exported JSON file
        """
        job_dir = self.create_job_directory(job_id)
        
        if output_filename is None:
            json_file = os.path.join(job_dir, "results.json")
        else:
            json_file = os.path.join(job_dir, output_filename)
        
        with open(json_file, 'w') as f:
            json.dump(data, f, indent=2)
        
        logger.info(f"Exported JSON results: {json_file}")
        return json_file
    
    def save_snapshot(
        self,
        job_id: str,
        image_data: str,  # Base64 data
        filename: str = "snapshot.png"
    ) -> str:
        """
        Save a snapshot image.
        
        Args:
            job_id: Job identifier
            image_data: Base64 encoded image data
            filename: Filename for the image
            
        Returns:
            Path to saved image
        """
        import base64
        
        job_dir = self.create_job_directory(job_id)
        snapshot_file = os.path.join(job_dir, "snapshots", filename)
        
        # Extract base64 data
        if ',' in image_data:
            header, data = image_data.split(',', 1)
        else:
            data = image_data
        
        # Decode and save
        image_bytes = base64.b64decode(data)
        with open(snapshot_file, 'wb') as f:
            f.write(image_bytes)
        
        logger.info(f"Saved snapshot: {snapshot_file}")
        return snapshot_file
    
    def export_zip_bundle(
        self,
        job_id: str,
        include_poses: bool = True,
        include_csv: bool = True,
        include_json: bool = True,
        include_snapshots: bool = True,
        output_filename: Optional[str] = None
    ) -> str:
        """
        Export all results as a ZIP bundle.
        
        Args:
            job_id: Job identifier
            include_poses: Include pose PDB files
            include_csv: Include CSV ranking
            include_json: Include JSON results
            include_snapshots: Include snapshot images
            output_filename: Optional custom filename
            
        Returns:
            Path to exported ZIP file
        """
        job_dir = self.create_job_directory(job_id)
        
        if output_filename is None:
            zip_file = os.path.join(self.base_output_dir, f"{job_id}_results.zip")
        else:
            zip_file = os.path.join(self.base_output_dir, output_filename)
        
        with zipfile.ZipFile(zip_file, 'w', zipfile.ZIP_DEFLATED) as zf:
            job_dir_name = os.path.basename(job_dir)
            
            # Add poses
            if include_poses:
                poses_dir = os.path.join(job_dir, "poses")
                if os.path.exists(poses_dir):
                    for root, _, files in os.walk(poses_dir):
                        for file in files:
                            full_path = os.path.join(root, file)
                            arcname = os.path.join(
                                job_dir_name,
                                "poses",
                                os.path.basename(file)
                            )
                            zf.write(full_path, arcname)
            
            # Add CSV
            if include_csv:
                csv_file = os.path.join(job_dir, "ranking.csv")
                if os.path.exists(csv_file):
                    zf.write(csv_file, os.path.join(job_dir_name, "ranking.csv"))
            
            # Add JSON
            if include_json:
                json_file = os.path.join(job_dir, "results.json")
                if os.path.exists(json_file):
                    zf.write(json_file, os.path.join(job_dir_name, "results.json"))
            
            # Add snapshots
            if include_snapshots:
                snaps_dir = os.path.join(job_dir, "snapshots")
                if os.path.exists(snaps_dir):
                    for root, _, files in os.walk(snaps_dir):
                        for file in files:
                            full_path = os.path.join(root, file)
                            arcname = os.path.join(
                                job_dir_name,
                                "snapshots",
                                os.path.basename(file)
                            )
                            zf.write(full_path, arcname)
        
        logger.info(f"Exported ZIP bundle: {zip_file}")
        return zip_file
    
    def get_result_files(self, job_id: str) -> List[ResultFile]:
        """
        Get list of all result files for a job.
        
        Args:
            job_id: Job identifier
            
        Returns:
            List of ResultFile objects
        """
        job_dir = os.path.join(self.base_output_dir, job_id)
        
        if not os.path.exists(job_dir):
            return []
        
        files = []
        
        for root, _, filenames in os.walk(job_dir):
            for filename in filenames:
                full_path = os.path.join(root, filename)
                rel_path = os.path.relpath(full_path, job_dir)
                
                # Determine file type
                ext = os.path.splitext(filename)[1].lower()
                file_type_map = {
                    '.pdb': 'pdbqt',
                    '.pdbqt': 'pdbqt',
                    '.csv': 'csv',
                    '.json': 'json',
                    '.pdf': 'pdf',
                    '.png': 'png',
                    '.jpg': 'jpg',
                    '.sdf': 'sdf',
                }
                
                files.append(ResultFile(
                    filename=filename,
                    path=full_path,
                    size=os.path.getsize(full_path),
                    file_type=file_type_map.get(ext, 'other')
                ))
        
        return files
    
    def get_job_summary(self, job_id: str) -> Dict:
        """Get summary of a job's results"""
        files = self.get_result_files(job_id)
        
        summary = {
            'job_id': job_id,
            'total_files': len(files),
            'total_size': sum(f.size for f in files),
            'file_types': {},
        }
        
        for f in files:
            if f.file_type not in summary['file_types']:
                summary['file_types'][f.file_type] = 0
            summary['file_types'][f.file_type] += 1
        
        return summary
    
    def delete_job_results(self, job_id: str) -> bool:
        """Delete all results for a job"""
        job_dir = os.path.join(self.base_output_dir, job_id)
        
        if os.path.exists(job_dir):
            shutil.rmtree(job_dir)
            logger.info(f"Deleted job results: {job_id}")
            return True
        
        return False


class FileExporter:
    """Utility class for file format conversions"""
    
    @staticmethod
    def pdbqt_to_sdf(pdb_content: str, output_path: str) -> bool:
        """Convert PDBQT to SDF (requires RDKit)"""
        try:
            from rdkit import Chem
            
            # Parse PDBQT - simplified conversion
            # In production, would need proper conversion
            mol = Chem.MolFromPDBBlock(pdb_content)
            
            if mol:
                writer = Chem.SDWriter(output_path)
                writer.write(mol)
                writer.close()
                return True
            
            return False
            
        except ImportError:
            logger.warning("RDKit not available for SDF conversion")
            return False
        except Exception as e:
            logger.error(f"SDF conversion failed: {e}")
            return False
    
    @staticmethod
    def generate_summary_html(
        job_id: str,
        poses: List[Dict],
        interactions: List[Dict],
        output_path: str
    ) -> bool:
        """Generate HTML summary of docking results"""
        html = f"""<!DOCTYPE html>
<html>
<head>
    <title>Docking Results - {job_id}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        h1 {{ color: #2E5AAC; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #2E5AAC; color: white; }}
        tr:nth-child(even) {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
    <h1>Docking Results: {job_id}</h1>
    <h2>Pose Rankings</h2>
    <table>
        <tr>
            <th>Rank</th>
            <th>Pose ID</th>
            <th>Score (kcal/mol)</th>
            <th>Interactions</th>
        </tr>
"""
        
        for i, pose in enumerate(poses, 1):
            html += f"""        <tr>
            <td>{i}</td>
            <td>{pose.get('id', 'N/A')}</td>
            <td>{pose.get('score', 'N/A'):.2f}</td>
            <td>{pose.get('interaction_count', 0)}</td>
        </tr>
"""
        
        html += """    </table>
    <h2>Interactions Summary</h2>
    <ul>
"""
        
        # Count interaction types
        interaction_counts = {}
        for interaction in interactions:
            itype = interaction.get('type', 'unknown')
            interaction_counts[itype] = interaction_counts.get(itype, 0) + 1
        
        for itype, count in interaction_counts.items():
            html += f"        <li>{itype}: {count}</li>\n"
        
        html += """    </ul>
</body>
</html>"""
        
        try:
            with open(output_path, 'w') as f:
                f.write(html)
            return True
        except Exception as e:
            logger.error(f"HTML summary generation failed: {e}")
            return False
