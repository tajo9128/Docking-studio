"""
Security Panel Widget
PyQt6 panel for displaying security status and running scans
"""

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTableWidget, QTableWidgetItem, QGroupBox, QFrame,
    QProgressBar, QHeaderView, QAbstractItemView, QComboBox,
    QTextEdit
)
from PyQt6.QtCore import Qt, pyqtSignal, QTimer
from PyQt6.QtGui import QColor, QFont
import logging

logger = logging.getLogger(__name__)

API_BASE_URL = "http://localhost:8000"


class SecurityPanel(QWidget):
    """
    Security monitoring panel for Agent Zero
    Displays security scan results and allows manual scans
    """
    
    scan_completed = pyqtSignal(dict)
    security_blocked = pyqtSignal(str)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.current_status = None
        self.scan_results = {}
        self._setup_ui()
        self._load_status()
    
    def _setup_ui(self):
        """Setup the security panel UI"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(5, 5, 5, 5)
        
        layout.addWidget(self._create_header_group())
        layout.addWidget(self._create_status_group())
        layout.addWidget(self._create_scans_group())
        layout.addWidget(self._create_details_group())
    
    def _create_header_group(self) -> QGroupBox:
        """Create header with title and controls"""
        group = QGroupBox("Agent Zero Security Monitor")
        group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                font-size: 14px;
            }
        """)
        
        layout = QHBoxLayout(group)
        
        self.status_icon = QLabel("🔒")
        self.status_icon.setStyleSheet("font-size: 24px;")
        
        self.status_label = QLabel("Security Status: Checking...")
        self.status_label.setStyleSheet("font-size: 14px; font-weight: bold;")
        
        self.scan_button = QPushButton("Run Scan")
        self.scan_button.setStyleSheet("""
            QPushButton {
                background-color: #2196F3;
                color: white;
                border: none;
                padding: 8px 16px;
                border-radius: 4px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #1976D2;
            }
            QPushButton:disabled {
                background-color: #BDBDBD;
            }
        """)
        self.scan_button.clicked.connect(self._run_scan)
        
        self.auto_scan_check = QComboBox()
        self.auto_scan_check.addItems(["Manual Only", "On Startup", "Daily", "Weekly"])
        self.auto_scan_check.setToolTip("Auto-scan frequency")
        
        layout.addWidget(self.status_icon)
        layout.addWidget(self.status_label, 1)
        layout.addWidget(self.auto_scan_check)
        layout.addWidget(self.scan_button)
        
        return group
    
    def _create_status_group(self) -> QGroupBox:
        """Create overall status display"""
        group = QGroupBox("Security Status")
        layout = QVBoxLayout(group)
        
        self.security_badge = QLabel("UNKNOWN")
        self.security_badge.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.security_badge.setStyleSheet("""
            QLabel {
                font-size: 24px;
                font-weight: bold;
                padding: 15px;
                border-radius: 8px;
                background-color: #9E9E9E;
                color: white;
            }
        """)
        
        self.last_scan_label = QLabel("Last scan: Never")
        self.last_scan_label.setStyleSheet("color: #666; font-size: 12px;")
        
        self.issues_count_label = QLabel("Total issues: 0")
        self.issues_count_label.setStyleSheet("font-size: 12px;")
        
        layout.addWidget(self.security_badge)
        layout.addWidget(self.last_scan_label)
        layout.addWidget(self.issues_count_label)
        
        return group
    
    def _create_scans_group(self) -> QGroupBox:
        """Create individual scan results"""
        group = QGroupBox("Scan Results")
        layout = QVBoxLayout(group)
        
        self.scans_table = QTableWidget()
        self.scans_table.setStyleSheet("""
            QTableWidget {
                background-color: white;
                border: 1px solid #E0E0E0;
            }
            QTableWidget::item {
                padding: 5px;
            }
        """)
        self.scans_table.setColumnCount(3)
        self.scans_table.setHorizontalHeaderLabels(["Scan Type", "Severity", "Issues"])
        self.scans_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.scans_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.scans_table.verticalHeader().setVisible(False)
        
        layout.addWidget(self.scans_table)
        
        return group
    
    def _create_details_group(self) -> QGroupBox:
        """Create details text area"""
        group = QGroupBox("Scan Details")
        layout = QVBoxLayout(group)
        
        self.details_text = QTextEdit()
        self.details_text.setMaximumHeight(150)
        self.details_text.setReadOnly(True)
        self.details_text.setStyleSheet("""
            QTextEdit {
                background-color: #FAFAFA;
                border: 1px solid #E0E0E0;
                font-family: monospace;
                font-size: 10px;
            }
        """)
        
        layout.addWidget(self.details_text)
        
        return group
    
    def _load_status(self):
        """Load current security status from API"""
        try:
            import requests
            
            response = requests.get(f"{API_BASE_URL}/security/status", timeout=10)
            
            if response.status_code == 200:
                self.current_status = response.json()
                self._update_display()
            else:
                self._update_display_not_scanned()
                
        except Exception as e:
            logger.warning(f"Could not load security status: {e}")
            self._update_display_not_scanned()
    
    def _update_display(self):
        """Update UI with current status"""
        if not self.current_status:
            return
        
        is_secure = self.current_status.get('is_secure', True)
        severity = self.current_status.get('overall_severity', 'UNKNOWN')
        total_issues = self.current_status.get('total_issues', 0)
        last_scan = self.current_status.get('last_scan_at', 'Never')
        
        if is_secure:
            self.status_icon.setText("🔒")
            self.security_badge.setText("SECURE")
            self.security_badge.setStyleSheet("""
                QLabel {
                    font-size: 24px;
                    font-weight: bold;
                    padding: 15px;
                    border-radius: 8px;
                    background-color: #4CAF50;
                    color: white;
                }
            """)
            self.status_label.setText("Security Status: Secure")
            self.status_label.setStyleSheet("font-size: 14px; font-weight: bold; color: #4CAF50;")
        else:
            self.status_icon.setText("⚠️")
            if severity == 'CRITICAL':
                self.security_badge.setText("CRITICAL")
                self.security_badge.setStyleSheet("""
                    QLabel {
                        font-size: 24px;
                        font-weight: bold;
                        padding: 15px;
                        border-radius: 8px;
                        background-color: #F44336;
                        color: white;
                    }
                """)
                self.status_label.setText("Security Status: CRITICAL ISSUES")
                self.status_label.setStyleSheet("font-size: 14px; font-weight: bold; color: #F44336;")
                self.security_blocked.emit("CRITICAL security issues detected")
            else:
                self.security_badge.setText("WARNING")
                self.security_badge.setStyleSheet("""
                    QLabel {
                        font-size: 24px;
                        font-weight: bold;
                        padding: 15px;
                        border-radius: 8px;
                        background-color: #FFC107;
                        color: black;
                    }
                """)
                self.status_label.setText("Security Status: Issues Found")
                self.status_label.setStyleSheet("font-size: 14px; font-weight: bold; color: #FFC107;")
        
        self.last_scan_label.setText(f"Last scan: {last_scan}")
        self.issues_count_label.setText(f"Total issues: {total_issues}")
        
        scan_results = self.current_status.get('scan_results', {}).get('scans', {})
        self._update_scans_table(scan_results)
    
    def _update_scans_table(self, scans: dict):
        """Update the scans table"""
        self.scans_table.setRowCount(0)
        
        scan_names = {
            'trivy_fs': 'Trivy (Filesystem)',
            'trivy_image': 'Trivy (Docker)',
            'bandit': 'Bandit (Python)',
            'safety': 'Safety (Dependencies)',
            'gitleaks': 'Gitleaks (Secrets)'
        }
        
        for scan_type, data in scans.items():
            row = self.scans_table.rowCount()
            self.scans_table.insertRow(row)
            
            name_item = QTableWidgetItem(scan_names.get(scan_type, scan_type))
            self.scans_table.setItem(row, 0, name_item)
            
            severity = data.get('severity', 'UNKNOWN')
            severity_item = QTableWidgetItem(severity)
            
            if severity == 'CRITICAL':
                severity_item.setBackground(QColor('#FFCDD2'))
            elif severity == 'HIGH':
                severity_item.setBackground(QColor('#FFE0B2'))
            elif severity == 'MEDIUM':
                severity_item.setBackground(QColor('#FFF9C4'))
            elif severity == 'SAFE':
                severity_item.setBackground(QColor('#C8E6C9'))
            
            self.scans_table.setItem(row, 1, severity_item)
            
            issues = data.get('issues', 0)
            issues_item = QTableWidgetItem(str(issues))
            self.scans_table.setItem(row, 2, issues_item)
    
    def _update_display_not_scanned(self):
        """Display when no scan has been performed"""
        self.status_icon.setText("❓")
        self.security_badge.setText("NOT SCANNED")
        self.security_badge.setStyleSheet("""
            QLabel {
                font-size: 24px;
                font-weight: bold;
                padding: 15px;
                border-radius: 8px;
                background-color: #9E9E9E;
                color: white;
            }
        """)
        self.status_label.setText("Security Status: Not Scanned")
        self.status_label.setStyleSheet("font-size: 14px; font-weight: bold; color: #9E9E9E;")
        self.last_scan_label.setText("Last scan: Never")
        self.issues_count_label.setText("Total issues: 0")
    
    def _run_scan(self):
        """Run security scan"""
        self.scan_button.setEnabled(False)
        self.scan_button.setText("Scanning...")
        self.status_label.setText("Security Status: Scanning...")
        self.details_text.setText("Running security scans...\n")
        
        try:
            import requests
            
            response = requests.post(f"{API_BASE_URL}/security/scan", timeout=300)
            
            if response.status_code == 200:
                result = response.json()
                self.scan_results = result
                self.current_status = result
                self._update_display()
                self.details_text.setText(f"Scan completed!\n\n{json.dumps(result, indent=2)}")
                self.scan_completed.emit(result)
                
                if not result.get('is_secure', True):
                    severity = result.get('overall_severity', 'UNKNOWN')
                    self.security_blocked.emit(f"Security scan found {severity} issues")
            else:
                self.details_text.setText(f"Scan failed: {response.status_code}")
                
        except Exception as e:
            logger.error(f"Security scan failed: {e}")
            self.details_text.setText(f"Scan error: {str(e)}")
        
        self.scan_button.setEnabled(True)
        self.scan_button.setText("Run Scan")
    
    def refresh(self):
        """Refresh security status"""
        self._load_status()
    
    def is_docking_allowed(self) -> bool:
        """Check if docking is allowed"""
        if not self.current_status:
            return True
        return self.current_status.get('is_secure', True)


import json


if __name__ == "__main__":
    from PyQt6.QtWidgets import QApplication
    app = QApplication([])
    panel = SecurityPanel()
    panel.show()
    app.exec()
