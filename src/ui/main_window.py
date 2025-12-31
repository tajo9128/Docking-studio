"""
BioDockify Docking Studio - Main Window
Main PyQt6 window for BioDockify Docking Studio
"""

from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
    QPushButton, QLabel, QStackedWidget, QFrame,
    QProgressBar, QScrollArea, QSizePolicy, QDockWidget
)
from PyQt6.QtCore import Qt, QSize, QTimer, pyqtSignal, QPoint
from PyQt6.QtGui import QIcon, QFont, QAction, QPalette, QColor
import logging

logger = logging.getLogger(__name__)

class MainWindow(QMainWindow):
    """Main application window"""
    
    # Signals
    job_started = pyqtSignal(str)
    job_completed = pyqtSignal(str)
    job_failed = pyqtSignal(str, str)
    job_cancelled = pyqtSignal(str)
    job_progress = pyqtSignal(str, int)  # job_id, percentage
    
    def __init__(self):
        """Initialize main window"""
        super().__init__()
        self.setWindowTitle("BioDockify Docking Studio v1.0.0")
        self.setMinimumSize(1280, 720)
        self.resize(1600, 900)
        
        # Central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout
        main_layout = QVBoxLayout(central_widget)
        
        # Welcome section
        welcome_label = QLabel("Welcome to BioDockify Docking Studio")
        welcome_label.setStyleSheet("font-size: 18px; font-weight: bold; color: #2196F3; padding: 10px;")
        main_layout.addWidget(welcome_label)
        
        # Job controls section
        job_controls_layout = QHBoxLayout()
        
        self.new_job_button = QPushButton("New Job")
        self.new_job_button.setStyleSheet("background-color: #4CAF50; color: white; border: none; padding: 10px; border-radius: 5px; font-weight: bold;")
        self.new_job_button.clicked.connect(self.on_new_job)
        job_controls_layout.addWidget(self.new_job_button)
        
        self.cancel_button = QPushButton("Stop Job")
        self.cancel_button.setStyleSheet("background-color: #F44336; color: white; border: none; padding: 10px; border-radius: 5px; font-weight: bold;")
        self.cancel_button.clicked.connect(self.on_cancel_job)
        self.cancel_button.setEnabled(False)
        job_controls_layout.addWidget(self.cancel_button)
        
        main_layout.addLayout(job_controls_layout)
        
        # Status section
        self.status_label = QLabel("Status: Ready")
        self.status_label.setStyleSheet("padding: 10px;")
        main_layout.addWidget(self.status_label)
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setStyleSheet("QProgressBar::chunk { background-color: #E8F5E9; } QProgressBar::text { color: white; } QProgressBar::chunk { border-radius: 4px; }")
        self.progress_bar.setValue(0)
        main_layout.addWidget(self.progress_bar)
        
        # Results section (would be implemented with ResultsWidget)
        results_placeholder = QLabel("Results will appear here after docking completes")
        results_placeholder.setStyleSheet("padding: 20px; font-style: italic; color: #666;")
        main_layout.addWidget(results_placeholder)
        
        # Add stretch to push content to top
        main_layout.addStretch()
        
        # Menu bar
        self._create_menu_bar()
    
    def _create_menu_bar(self) -> None:
        """Create menu bar"""
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu("File")
        
        new_job_action = file_menu.addAction("New Job")
        new_job_action.setShortcut("Ctrl+N")
        new_job_action.triggered.connect(self.on_new_job)
        
        open_recent_job_action = file_menu.addAction("Open Recent Job...")
        open_recent_job_action.setShortcut("Ctrl+O")
        
        file_menu.addSeparator()
        
        export_results_action = file_menu.addAction("Export Results...")
        export_results_action.setShortcut("Ctrl+E")
        export_results_action.triggered.connect(self.on_export_results)
        
        file_menu.addSeparator()
        
        exit_action = file_menu.addAction("Exit")
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        
        # Edit menu
        edit_menu = menubar.addMenu("Edit")
        
        preferences_action = edit_menu.addAction("Preferences...")
        preferences_action.setShortcut("Ctrl+,")
        preferences_action.triggered.connect(self.on_preferences)
        
        clear_history_action = edit_menu.addAction("Clear Job History")
        clear_history_action.triggered.connect(self.on_clear_history)
        
        # View menu
        view_menu = menubar.addMenu("View")
        
        show_progress_action = view_menu.addAction("Show Progress")
        show_progress_action.setCheckable(True)
        show_progress_action.triggered.connect(self.on_show_progress)
        
        hide_progress_action = view_menu.addAction("Hide Progress")
        hide_progress_action.setCheckable(False)
        hide_progress_action.triggered.connect(self.on_hide_progress)
        
        view_menu.addSeparator()
        
        show_agent_zero_action = view_menu.addAction("Show Agent Zero Status")
        show_agent_zero_action.setCheckable(True)
        show_agent_zero_action.triggered.connect(self.on_show_agent_zero)
        
        hide_agent_zero_action = view_menu.addAction("Hide Agent Zero Status")
        hide_agent_zero_action.setCheckable(False)
        hide_agent_zero_action.triggered.connect(self.on_hide_agent_zero)
        
        # Help menu
        help_menu = menubar.addMenu("Help")
        
        documentation_action = help_menu.addAction("Documentation...")
        documentation_action.setShortcut("F1")
        documentation_action.triggered.connect(self.on_documentation)
        
        troubleshooting_action = help_menu.addAction("Troubleshooting...")
        troubleshooting_action.setShortcut("F2")
        troubleshooting_action.triggered.connect(self.on_troubleshooting)
        
        faq_action = help_menu.addAction("FAQ")
        faq_action.setShortcut("F3")
        faq_action.triggered.connect(self.on_faq)
        
        help_menu.addSeparator()
        
        check_for_updates_action = help_menu.addAction("Check for Updates...")
        check_for_updates_action.triggered.connect(self.on_check_for_updates)
        
        about_action = help_menu.addAction("About BioDockify Docking Studio")
        about_action.triggered.connect(self.on_about)
        
        # Status bar
        self._create_status_bar()
    
    def _create_status_bar(self) -> None:
        """Create status bar"""
        self.statusBar = self.statusBar()
        
        # Docker status
        self.docker_status_label = QLabel("Docker: Ready")
        self.docker_status_label.setStyleSheet("color: #4CAF50; font-weight: bold; padding: 5px;")
        self.statusBar.addPermanentWidget(self.docker_status_label)
        
        # Version info
        version_label = QLabel("v1.0.0")
        version_label.setStyleSheet("padding: 5px;")
        self.statusBar.addPermanentWidget(version_label)
    
    def on_new_job(self) -> None:
        """Handle new job button click"""
        logger.info("New job button clicked")
        self.job_started.emit("new_job")
    
    def on_cancel_job(self) -> None:
        """Handle cancel job button click"""
        logger.info("Cancel job button clicked")
        self.job_cancelled.emit("current_job")
        self.cancel_button.setEnabled(False)
    
    def on_export_results(self) -> None:
        """Handle export results menu action"""
        logger.info("Export results menu action clicked")
        # In real implementation, would export to CSV, JSON, or PyMOL
    
    def on_preferences(self) -> None:
        """Handle preferences menu action"""
        logger.info("Preferences menu action clicked")
        # In real implementation, would open preferences dialog
    
    def on_clear_history(self) -> None:
        """Handle clear history menu action"""
        logger.info("Clear job history menu action clicked")
        # In real implementation, would clear job history from database
    
    def on_show_progress(self) -> None:
        """Handle show progress menu action"""
        logger.info("Show progress menu action clicked")
        # In real implementation, would show progress widget
    
    def on_hide_progress(self) -> None:
        """Handle hide progress menu action"""
        logger.info("Hide progress menu action clicked")
        # In real implementation, would hide progress widget
    
    def on_show_agent_zero(self) -> None:
        """Handle show Agent Zero menu action"""
        logger.info("Show Agent Zero status menu action clicked")
        # In real implementation, would show Agent Zero widget
    
    def on_hide_agent_zero(self) -> None:
        """Handle hide Agent Zero menu action"""
        logger.info("Hide Agent Zero status menu action clicked")
        # In real implementation, would hide Agent Zero widget
    
    def on_documentation(self) -> None:
        """Handle documentation menu action"""
        logger.info("Documentation menu action clicked")
        # In real implementation, would open documentation in browser
    
    def on_troubleshooting(self) -> None:
        """Handle troubleshooting menu action"""
        logger.info("Troubleshooting menu action clicked")
        # In real implementation, would open troubleshooting documentation in browser
    
    def on_faq(self) -> None:
        """Handle FAQ menu action"""
        logger.info("FAQ menu action clicked")
        # In real implementation, would open FAQ documentation in browser
    
    def on_check_for_updates(self) -> None:
        """Handle check for updates menu action"""
        logger.info("Check for updates menu action clicked")
        # In real implementation, would check GitHub for new releases
    
    def on_about(self) -> None:
        """Handle about menu action"""
        logger.info("About menu action clicked")
        
        from PyQt6.QtWidgets import QMessageBox
        QMessageBox.about(self, 
                         "BioDockify Docking Studio",
                         "Molecular Docking with Intelligent Self-Repair",
                         f"Version: 1.0.0\n\n"
                         "Copyright (c) 2025 BioDockify Development Team\n\n"
                         "License: Apache License 2.0\n\n"
                         "Built with: AutoDock Vina, ODDT, RDKit\n"
                         "Docker-based execution environment\n\n"
                         "Agent Zero AI system for intelligent failure detection and recovery")
    
    def update_progress(self, percentage: int, message: str = "") -> None:
        """Update progress bar"""
        self.progress_bar.setValue(percentage)
        self.status_label.setText(f"Status: {message}")
        self.job_progress.emit("current_job", percentage)
        logger.debug(f"Progress updated: {percentage}% - {message}")
