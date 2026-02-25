"""
Analysis Panel - RMSD, Interactions, Job History
Integrated with FastAPI backend for Docker Desktop deployment
"""

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTableWidget, QTableWidgetItem, QGroupBox, QFrame, QComboBox,
    QLineEdit, QTextEdit, QScrollArea, QProgressBar, QHeaderView,
    QAbstractItemView, QSplitter
)
from PyQt6.QtCore import Qt, pyqtSignal, QTimer
from PyQt6.QtGui import QFont, QColor
import logging
import json
import os

logger = logging.getLogger(__name__)

API_BASE_URL = "http://localhost:8000"


class AnalysisPanel(QWidget):
    """
    Analysis panel combining RMSD, Interactions, and Job History
    Docker Desktop compatible with FastAPI backend
    """
    
    analysis_complete = pyqtSignal(dict)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.current_job = None
        self.analysis_results = {}
        self._setup_ui()
    
    def _setup_ui(self):
        """Setup the analysis panel UI"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(5, 5, 5, 5)
        
        splitter = QSplitter(Qt.Orientation.Vertical)
        
        splitter.addWidget(self._create_rmsd_group())
        splitter.addWidget(self._create_interactions_group())
        splitter.addWidget(self._create_job_history_group())
        
        splitter.setSizes([250, 300, 250])
        
        layout.addWidget(splitter)
    
    def _create_rmsd_group(self) -> QGroupBox:
        """Create RMSD calculation group"""
        group = QGroupBox("RMSD Analysis")
        layout = QVBoxLayout(group)
        
        input_layout = QHBoxLayout()
        
        self.pose1_combo = QComboBox()
        self.pose1_combo.setMinimumWidth(150)
        self.pose1_combo.addItem("Select Pose 1")
        
        self.pose2_combo = QComboBox()
        self.pose2_combo.setMinimumWidth(150)
        self.pose2_combo.addItem("Select Pose 2")
        
        self.calculate_rmsd_btn = QPushButton("Calculate RMSD")
        self.calculate_rmsd_btn.setStyleSheet("""
            QPushButton {
                background-color: #2196F3;
                color: white;
                border: none;
                padding: 8px 15px;
                border-radius: 4px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #1976D2;
            }
        """)
        self.calculate_rmsd_btn.clicked.connect(self._calculate_rmsd)
        
        input_layout.addWidget(QLabel("Pose 1:"))
        input_layout.addWidget(self.pose1_combo)
        input_layout.addWidget(QLabel("Pose 2:"))
        input_layout.addWidget(self.pose2_combo)
        input_layout.addWidget(self.calculate_rmsd_btn)
        
        layout.addLayout(input_layout)
        
        self.rmsd_result_label = QLabel("RMSD: --")
        self.rmsd_result_label.setStyleSheet("""
            QLabel {
                font-size: 18px;
                font-weight: bold;
                color: #4CAF50;
                padding: 10px;
            }
        """)
        layout.addWidget(self.rmsd_result_label)
        
        self.rmsd_info = QTextEdit()
        self.rmsd_info.setMaximumHeight(80)
        self.rmsd_info.setReadOnly(True)
        self.rmsd_info.setStyleSheet("""
            QTextEdit {
                background-color: #F5F5F5;
                border: 1px solid #E0E0E0;
                border-radius: 4px;
                font-size: 11px;
            }
        """)
        layout.addWidget(self.rmsd_info)
        
        return group
    
    def _create_interactions_group(self) -> QGroupBox:
        """Create interactions display group"""
        group = QGroupBox("Molecular Interactions")
        layout = QVBoxLayout(group)
        
        toolbar = QHBoxLayout()
        
        self.analyze_interactions_btn = QPushButton("Analyze Interactions")
        self.analyze_interactions_btn.setStyleSheet("""
            QPushButton {
                background-color: #9C27B0;
                color: white;
                border: none;
                padding: 6px 12px;
                border-radius: 4px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #7B1FA2;
            }
        """)
        self.analyze_interactions_btn.clicked.connect(self._analyze_interactions)
        
        self.show_hbonds_btn = QPushButton("H-Bonds")
        self.show_hbonds_btn.setCheckable(True)
        self.show_hbonds_btn.setChecked(True)
        
        self.show_hydrophobic_btn = QPushButton("Hydrophobic")
        self.show_hydrophobic_btn.setCheckable(True)
        self.show_hydrophobic_btn.setChecked(True)
        
        toolbar.addWidget(self.analyze_interactions_btn)
        toolbar.addWidget(self.show_hbonds_btn)
        toolbar.addWidget(self.show_hydrophobic_btn)
        toolbar.addStretch()
        
        layout.addLayout(toolbar)
        
        self.interactions_table = QTableWidget()
        self.interactions_table.setStyleSheet("""
            QTableWidget {
                background-color: white;
                border: 1px solid #CCCCCC;
            }
            QTableWidget::item {
                padding: 4px;
            }
            QHeaderView::section {
                background-color: #F5F5F5;
                font-weight: bold;
                padding: 4px;
            }
        """)
        self.interactions_table.setColumnCount(4)
        self.interactions_table.setHorizontalHeaderLabels(["Type", "Residue A", "Residue B", "Distance (Å)"])
        self.interactions_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.interactions_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.interactions_table.verticalHeader().setVisible(False)
        
        layout.addWidget(self.interactions_table)
        
        self.interaction_summary = QLabel("No interactions analyzed")
        self.interaction_summary.setStyleSheet("color: #666; font-style: italic;")
        layout.addWidget(self.interaction_summary)
        
        return group
    
    def _create_job_history_group(self) -> QGroupBox:
        """Create job history group"""
        group = QGroupBox("Job History")
        layout = QVBoxLayout(group)
        
        toolbar = QHBoxLayout()
        
        self.refresh_jobs_btn = QPushButton("Refresh")
        self.refresh_jobs_btn.setFixedWidth(80)
        self.refresh_jobs_btn.clicked.connect(self._load_jobs)
        
        self.clear_jobs_btn = QPushButton("Clear All")
        self.clear_jobs_btn.setFixedWidth(80)
        self.clear_jobs_btn.clicked.connect(self._clear_all_jobs)
        
        toolbar.addWidget(self.refresh_jobs_btn)
        toolbar.addWidget(self.clear_jobs_btn)
        toolbar.addStretch()
        
        layout.addLayout(toolbar)
        
        self.jobs_table = QTableWidget()
        self.jobs_table.setStyleSheet("""
            QTableWidget {
                background-color: white;
                border: 1px solid #CCCCCC;
            }
            QTableWidget::item {
                padding: 4px;
            }
            QHeaderView::section {
                background-color: #F5F5F5;
                font-weight: bold;
                padding: 4px;
            }
        """)
        self.jobs_table.setColumnCount(5)
        self.jobs_table.setHorizontalHeaderLabels(["Job Name", "Status", "Energy", "Confidence", "Created"])
        self.jobs_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.jobs_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.jobs_table.verticalHeader().setVisible(False)
        self.jobs_table.cellClicked.connect(self._on_job_selected)
        
        layout.addWidget(self.jobs_table)
        
        return group
    
    def _calculate_rmsd(self):
        """Calculate RMSD between two poses"""
        pose1_idx = self.pose1_combo.currentIndex()
        pose2_idx = self.pose2_combo.currentIndex()
        
        if pose1_idx <= 0 or pose2_idx <= 0:
            self.rmsd_result_label.setText("Please select both poses")
            return
        
        try:
            import requests
            
            pose1_data = self.poses[pose1_idx - 1]
            pose2_data = self.poses[pose2_idx - 1]
            
            response = requests.post(
                f"{API_BASE_URL}/rmsd",
                json={
                    "pdb1": pose1_data.get("pdb", ""),
                    "pdb2": pose2_data.get("pdb", "")
                },
                timeout=30
            )
            
            if response.status_code == 200:
                result = response.json()
                rmsd = result.get("rmsd", -1)
                
                if rmsd >= 0:
                    self.rmsd_result_label.setText(f"RMSD: {rmsd:.3f} Å")
                    
                    if rmsd < 1.0:
                        self.rmsd_result_label.setStyleSheet("""
                            QLabel {
                                font-size: 18px;
                                font-weight: bold;
                                color: #4CAF50;
                                padding: 10px;
                            }
                        """)
                        info = "Excellent! Very similar poses."
                    elif rmsd < 2.0:
                        self.rmsd_result_label.setStyleSheet("""
                            QLabel {
                                font-size: 18px;
                                font-weight: bold;
                                color: #FFC107;
                                padding: 10px;
                            }
                        """)
                        info = "Good similarity between poses."
                    else:
                        self.rmsd_result_label.setStyleSheet("""
                            QLabel {
                                font-size: 18px;
                                font-weight: bold;
                                color: #F44336;
                                padding: 10px;
                            }
                        """)
                        info = "Different poses."
                    
                    self.rmsd_info.setText(info)
                else:
                    self.rmsd_result_label.setText("RMSD calculation failed")
            else:
                self.rmsd_result_label.setText(f"Error: {response.status_code}")
                
        except Exception as e:
            logger.error(f"RMSD calculation error: {e}")
            self.rmsd_result_label.setText(f"Error: {str(e)[:30]}")
            self.rmsd_info.setText("Using local calculation...")
            self._calculate_rmsd_local(pose1_idx - 1, pose2_idx - 1)
    
    def _calculate_rmsd_local(self, idx1: int, idx2: int):
        """Fallback local RMSD calculation"""
        try:
            from src.backend.analysis import calculate_rmsd
            
            pose1_data = self.poses[idx1]
            pose2_data = self.poses[idx2]
            
            rmsd = calculate_rmsd(
                pose1_data.get("pdb", ""),
                pose2_data.get("pdb", "")
            )
            
            self.rmsd_result_label.setText(f"RMSD: {rmsd:.3f} Å")
            
        except Exception as e:
            self.rmsd_result_label.setText("RMSD unavailable")
            logger.error(f"Local RMSD error: {e}")
    
    def _analyze_interactions(self):
        """Analyze molecular interactions"""
        if not self.current_job:
            self.interaction_summary.setText("No job loaded")
            return
        
        try:
            import requests
            
            receptor = self.current_job.get("receptor_pdb", "")
            ligand = self.current_job.get("ligand_pdb", "")
            
            if not receptor or not ligand:
                self.interaction_summary.setText("Missing receptor/ligand data")
                return
            
            response = requests.post(
                f"{API_BASE_URL}/analyze/advanced",
                json={"receptor": receptor, "ligand": ligand},
                timeout=30
            )
            
            if response.status_code == 200:
                result = response.json()
                self._display_interactions(result)
            else:
                self.interaction_summary.setText(f"Error: {response.status_code}")
                
        except Exception as e:
            logger.error(f"Interaction analysis error: {e}")
            self.interaction_summary.setText(f"Error: {str(e)[:40]}")
    
    def _display_interactions(self, data: dict):
        """Display interactions in table"""
        self.interactions_table.setRowCount(0)
        
        total = 0
        
        hbonds = data.get("hbonds", [])
        if self.show_hbonds_btn.isChecked():
            for hbond in hbonds[:10]:
                row = self.interactions_table.rowCount()
                self.interactions_table.insertRow(row)
                self.interactions_table.setItem(row, 0, QTableWidgetItem("H-Bond"))
                self.interactions_table.setItem(row, 1, QTableWidgetItem(str(hbond.get("receptor_idx", "-"))))
                self.interactions_table.setItem(row, 2, QTableWidgetItem(str(hbond.get("ligand_idx", "-"))))
                self.interactions_table.setItem(row, 3, QTableWidgetItem(f"{hbond.get('distance', 0):.2f}"))
                total += 1
        
        hydrophobic = data.get("hydrophobic", [])
        if self.show_hydrophobic_btn.isChecked():
            for hydro in hydrophobic[:10]:
                row = self.interactions_table.rowCount()
                self.interactions_table.insertRow(row)
                self.interactions_table.setItem(row, 0, QTableWidgetItem("Hydrophobic"))
                self.interactions_table.setItem(row, 1, QTableWidgetItem(str(hydro.get("receptor_idx", "-"))))
                self.interactions_table.setItem(row, 2, QTableWidgetItem(str(hydro.get("ligand_idx", "-"))))
                self.interactions_table.setItem(row, 3, QTableWidgetItem(f"{hydro.get('distance', 0):.2f}"))
                total += 1
        
        self.interaction_summary.setText(
            f"Total interactions: {total} "
            f"(H-bonds: {len(hbonds)}, Hydrophobic: {len(hydrophobic)})"
        )
        
        self.analysis_results["interactions"] = data
        self.analysis_complete.emit(data)
    
    def _load_jobs(self):
        """Load job history from API"""
        try:
            import requests
            
            response = requests.get(f"{API_BASE_URL}/jobs", timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                jobs = data.get("jobs", [])
                self._display_jobs(jobs)
            else:
                logger.warning(f"Failed to load jobs: {response.status_code}")
                
        except Exception as e:
            logger.error(f"Load jobs error: {e}")
    
    def _display_jobs(self, jobs: list):
        """Display jobs in table"""
        self.jobs_table.setRowCount(0)
        
        for job in jobs:
            row = self.jobs_table.rowCount()
            self.jobs_table.insertRow(row)
            
            self.jobs_table.setItem(row, 0, QTableWidgetItem(job.get("job_name", "Unknown")))
            
            status = job.get("status", "unknown")
            status_item = QTableWidgetItem(status.upper())
            
            if status == "completed":
                status_item.setBackground(QColor("#C8E6C9"))
            elif status == "running":
                status_item.setBackground(QColor("#FFF9C4"))
            elif status == "failed":
                status_item.setBackground(QColor("#FFCDD2"))
            
            self.jobs_table.setItem(row, 1, status_item)
            
            energy = job.get("binding_energy")
            self.jobs_table.setItem(row, 2, QTableWidgetItem(f"{energy:.2f}" if energy else "-"))
            
            confidence = job.get("confidence_score")
            self.jobs_table.setItem(row, 3, QTableWidgetItem(f"{confidence:.0f}" if confidence else "-"))
            
            created = job.get("created_at", "")
            self.jobs_table.setItem(row, 4, QTableWidgetItem(created[:19] if created else "-"))
    
    def _on_job_selected(self, row: int):
        """Handle job selection"""
        job_name = self.jobs_table.item(row, 0).text()
        logger.info(f"Selected job: {job_name}")
    
    def _clear_all_jobs(self):
        """Clear all jobs"""
        try:
            import requests
            
            response = requests.get(f"{API_BASE_URL}/jobs", timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                jobs = data.get("jobs", [])
                
                for job in jobs:
                    job_uuid = job.get("job_uuid")
                    if job_uuid:
                        requests.delete(f"{API_BASE_URL}/jobs/{job_uuid}", timeout=10)
                
                self._load_jobs()
                
        except Exception as e:
            logger.error(f"Clear jobs error: {e}")
    
    def set_poses(self, poses: list):
        """Set poses for RMSD calculation"""
        self.poses = poses
        
        self.pose1_combo.clear()
        self.pose2_combo.clear()
        
        self.pose1_combo.addItem("Select Pose 1")
        self.pose2_combo.addItem("Select Pose 2")
        
        for i, pose in enumerate(poses):
            label = f"Pose {i+1}"
            if "score" in pose:
                label += f" ({pose['score']:.2f})"
            
            self.pose1_combo.addItem(label)
            self.pose2_combo.addItem(label)
    
    def set_current_job(self, job_data: dict):
        """Set current job for analysis"""
        self.current_job = job_data
        
        if "poses" in job_data:
            self.set_poses(job_data["poses"])
    
    def clear(self):
        """Clear all data"""
        self.current_job = None
        self.analysis_results = {}
        self.poses = []
        
        self.pose1_combo.clear()
        self.pose1_combo.addItem("Select Pose 1")
        self.pose2_combo.clear()
        self.pose2_combo.addItem("Select Pose 2")
        
        self.rmsd_result_label.setText("RMSD: --")
        self.rmsd_info.clear()
        self.interactions_table.setRowCount(0)
        self.interaction_summary.setText("No interactions analyzed")
