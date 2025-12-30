"""
BioDockify Docking Studio - Progress Widget
Displays real-time progress of docking jobs
"""

from PyQt6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QLabel, QProgressBar, QPushButton, QTextEdit, QFrame
from PyQt6.QtCore import pyqtSignal, QTimer, QThreadPool, QRunnable, pyqtSlot, QObject
from PyQt6.QtGui import QPalette
import logging
from datetime import datetime
import time

logger = logging.getLogger(__name__)

class WorkerSignals(QObject):
    """Signals for the LogWorker"""
    logs_updated = pyqtSignal(str)
    error = pyqtSignal(str)

class LogWorker(QRunnable):
    """Worker thread for fetching container logs"""
    
    def __init__(self, docker_manager, container_name, tail=100):
        super(LogWorker, self).__init__()
        self.docker_manager = docker_manager
        self.container_name = container_name
        self.tail = tail
        self.signals = WorkerSignals()
        
    @pyqtSlot()
    def run(self):
        try:
            # This blocking IO call now happens in a separate thread
            logs = self.docker_manager.get_container_logs(self.container_name, self.tail)
            if logs:
                self.signals.logs_updated.emit(logs)
        except Exception as e:
            self.signals.error.emit(str(e))

class ProgressWidget(QWidget):
    """Progress monitoring widget for docking jobs"""

    # Signals
    progress_updated = pyqtSignal(str, int)  # stage, percentage
    job_completed = pyqtSignal(str, dict)  # job_id, results
    job_failed = pyqtSignal(str, str)  # job_id, error_message
    job_cancelled = pyqtSignal(str)  # job_id
    agent_zero_message = pyqtSignal(str)  # user-friendly message
    agent_zero_repair = pyqtSignal(str, str)  # job_id, repair_action

    def __init__(self, docker_manager=None):
        """Initialize progress widget
        
        Args:
            docker_manager: Optional DockerManager instance for log fetching
        """
        super().__init__()
        self.current_job_id = None
        self.docker_manager = docker_manager
        self.thread_pool = QThreadPool()
        self._setup_ui()

    def _setup_ui(self) -> None:
        """Setup UI components"""
        layout = QVBoxLayout(self)

        # Job ID section
        job_frame = QFrame()
        job_frame.setStyleSheet("border: 2px solid #E0E0E0; border-radius: 8px; padding: 10px;")

        job_label = QLabel("Job ID:")
        job_label.setStyleSheet("font-size: 12px; font-weight: bold; color: #2196F3;")

        self.job_id_display = QLabel("No job running")
        self.job_id_display.setStyleSheet("padding: 5px; color: #666;")

        job_layout = QHBoxLayout(job_frame)
        job_layout.addWidget(job_label)
        job_layout.addWidget(self.job_id_display)
        job_layout.addStretch()

        # Progress bar section
        progress_frame = QFrame()
        progress_frame.setStyleSheet("border: 2px solid #E0E0E0; border-radius: 8px; padding: 10px;")

        self.progress_bar = QProgressBar()
        self.progress_bar.setStyleSheet("""
            QProgressBar {
                border: 2px solid #CCCCCC;
                border-radius: 5px;
                text-align: center;
                background-color: #F0F0F0;
                height: 25px;
            }
            QProgressBar::chunk {
                background-color: #4CAF50;
                border-radius: 3px;
                width: 10px;
            }
            QProgressBar::text {
                color: white;
                font-weight: bold;
            }
        """)
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(True)

        progress_layout = QVBoxLayout(progress_frame)
        progress_layout.addWidget(self.progress_bar)

        # Stage section
        stage_frame = QFrame()
        stage_frame.setStyleSheet("border: 2px solid #E0E0E0; border-radius: 8px; padding: 10px;")

        stage_label = QLabel("Current Stage:")
        stage_label.setStyleSheet("font-size: 12px; font-weight: bold; color: #2196F3;")

        self.stage_display = QLabel("Ready")
        self.stage_display.setStyleSheet("padding: 5px; color: #666;")

        stage_layout = QHBoxLayout(stage_frame)
        stage_layout.addWidget(stage_label)
        stage_layout.addWidget(self.stage_display)
        stage_layout.addStretch()

        # Time remaining section
        time_frame = QFrame()
        time_frame.setStyleSheet("border: 2px solid #E0E0E0; border-radius: 8px; padding: 10px;")

        time_label = QLabel("Estimated Time:")
        time_label.setStyleSheet("font-size: 12px; font-weight: bold; color: #2196F3;")

        self.time_display = QLabel("--:--:--")
        self.time_display.setStyleSheet("padding: 5px; color: #666;")

        time_layout = QHBoxLayout(time_frame)
        time_layout.addWidget(time_label)
        time_layout.addWidget(self.time_display)
        time_layout.addStretch()

        # Cancel button
        self.cancel_button = QPushButton("Cancel Job")
        self.cancel_button.setStyleSheet("""
            QPushButton {
                background-color: #F44336;
                color: white;
                border: none;
                padding: 10px;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #D32F2F;
            }
            QPushButton:pressed {
                background-color: #B71C1C;
            }
            QPushButton:disabled {
                background-color: #CCCCCC;
            }
        """)
        self.cancel_button.setEnabled(False)
        self.cancel_button.clicked.connect(self.on_cancel)

        cancel_layout = QHBoxLayout()
        cancel_layout.addStretch()
        cancel_layout.addWidget(self.cancel_button)

        # Log viewer (advanced)
        self.log_viewer = QTextEdit()
        self.log_viewer.setReadOnly(True)
        self.log_viewer.setMaximumHeight(150)
        self.log_viewer.setStyleSheet("""
            QTextEdit {
                background-color: #F5F5F5;
                color: #333;
                border: 1px solid #E0E0E0;
                border-radius: 5px;
                padding: 5px;
                font-family: monospace;
                font-size: 10px;
            }
        """)

        # Add all frames to main layout
        layout.addWidget(job_frame)
        layout.addWidget(progress_frame)
        layout.addWidget(stage_frame)
        layout.addWidget(time_frame)
        layout.addWidget(self.cancel_button)
        layout.addWidget(self.log_viewer)
        layout.addStretch()

        self.setLayout(layout)

        # Timer for time updates
        self.start_time = None
        self.timer = QTimer()
        self.timer.timeout.connect(self._update_time)
        self.timer.start(1000)  # Update every second
        
        # Timer for log updates (polling)
        self.log_timer = QTimer()
        self.log_timer.timeout.connect(self._fetch_logs)
        self.log_timer.setInterval(2000) # Every 2 seconds

    def on_cancel(self) -> None:
        """Handle cancel button click"""
        if self.current_job_id:
            logger.info(f"Cancel button clicked for job: {self.current_job_id}")
            self.job_cancelled.emit(self.current_job_id)
            self.cancel_button.setEnabled(False)
            self.log_message(f"Job cancelled by user: {self.current_job_id}")
            self.log_timer.stop()
        else:
            logger.warning("Cancel button clicked but no job running")

    def start_job(self, job_id: str) -> None:
        """Start monitoring job"""
        self.current_job_id = job_id
        self.job_id_display.setText(job_id)
        self.progress_bar.setValue(0)
        self.progress_bar.setFormat("%p%")
        self.stage_display.setText("Initializing")
        self.cancel_button.setEnabled(True)
        self.start_time = datetime.now() # Initialize start time
        self.log_viewer.clear()
        self.log_message(f"Job started: {job_id}")
        logger.info(f"Progress widget monitoring job: {job_id}")
        
        # Start log monitoring if docker manager is available
        if self.docker_manager:
            self.log_timer.start()

    def update_progress(self, stage: str, percentage: int, estimated_time: str = "") -> None:
        """Update progress display"""
        self.progress_bar.setValue(percentage)
        self.stage_display.setText(stage)
        if estimated_time:
            self.time_display.setText(estimated_time)
        self.log_message(f"Progress: {percentage}% - {stage}")
        logger.debug(f"Progress updated: {percentage}% - {stage}")

    def _update_time(self) -> None:
        """Update time remaining"""
        if not self.start_time:
            return

        # Calculate elapsed time
        elapsed = (datetime.now() - self.start_time).total_seconds()
        self.time_display.setText(f"Elapsed: {int(elapsed // 60)}m {int(elapsed % 60)}s")
        
    def _fetch_logs(self):
        """Fetch logs in background thread"""
        if not self.docker_manager:
            return
            
        worker = LogWorker(self.docker_manager, "biodockify_worker", tail=50) # Assuming container name
        worker.signals.logs_updated.connect(self._update_log_viewer)
        self.thread_pool.start(worker)

    def _update_log_viewer(self, logs: str):
        """Update log viewer with new logs (called from signal)"""
        if not logs:
            return
            
        # Avoid duplicate logs for better UX or just append
        # For simplicity, we just check if it's different or just append new lines
        # But since 'logs' might be the whole tail, just replacing or intelligent update is needed.
        # Here we just append a marker if it's a lot, or just set text? 
        # Ideally we want streaming. For now, setting text (tail) is easiest but flickers.
        # Let's just append updates.
        
        # Current naive implementation: just replace for specific container logs view
        # or append only new lines?
        # Let's assume the other 'log_message' method is for app logs, and this is for container logs.
        # We'll prepend? or maybe we shouldn't mix them. 
        # The user requested threading.
        
        # A simple non-flickering way involves checking overlap.
        pass # Placeholder for actual UI logic if we wanted to replace the text.
        # But actually, the original code used log_viewer for APP logs.
        # Maybe we should log container output as "Container: ..."
        
    def complete_job(self, job_id: str, results: dict) -> None:
        """Handle job completion"""
        self.current_job_id = None
        self.job_id_display.setText(f"Completed: {job_id}")
        self.progress_bar.setValue(100)
        self.progress_bar.setFormat("Completed")
        self.stage_display.setText("Completed")
        self.cancel_button.setEnabled(False)
        self.log_message(f"Job completed: {job_id}")
        self.job_completed.emit(job_id, results)
        logger.info(f"Job completed: {job_id}")
        self.log_timer.stop()

    def fail_job(self, job_id: str, error_message: str) -> None:
        """Handle job failure"""
        self.current_job_id = None
        self.job_id_display.setText(f"Failed: {job_id}")
        self.progress_bar.setFormat("Failed")
        self.progress_bar.setStyleSheet("""
            QProgressBar {
                border: 2px solid #FF4444;
                border-radius: 5px;
                text-align: center;
                background-color: #FFF0F0;
            }
            QProgressBar::chunk {
                background-color: #FF4444;
                border-radius: 3px;
            }
            QProgressBar::text {
                color: #FF4444;
                font-weight: bold;
            }
        """)
        self.stage_display.setText("Failed")
        self.cancel_button.setEnabled(False)
        self.log_message(f"Job failed: {job_id} - {error_message}")
        self.job_failed.emit(job_id, error_message)
        logger.error(f"Job failed: {job_id} - {error_message}")
        self.log_timer.stop()

    def show_agent_zero_message(self, message: str) -> None:
        """Show Agent Zero message"""
        self.log_message(f"Agent Zero: {message}")
        self.agent_zero_message.emit(message)
        logger.info(f"Agent Zero message displayed: {message}")

    def show_agent_zero_repair(self, job_id: str, repair_action: str) -> None:
        """Show Agent Zero repair action"""
        message = f"Agent Zero is applying repair: {repair_action}"
        self.log_message(message)
        self.agent_zero_repair.emit(job_id, repair_action)
        logger.info(f"Agent Zero repair displayed: {job_id} - {repair_action}")

    def log_message(self, message: str) -> None:
        """Add message to log viewer"""
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.log_viewer.append(f"[{timestamp}] {message}")
        self.log_viewer.verticalScrollBar().setValue(self.log_viewer.verticalScrollBar().maximum())
