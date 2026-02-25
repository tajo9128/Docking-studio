"""
Log Viewer Widget for PyQt6
Displays application logs with filtering and search
"""

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTextEdit, 
    QPushButton, QLabel, QComboBox, QLineEdit,
    QFrame, QScrollBar
)
from PyQt6.QtCore import Qt, pyqtSignal, QTimer
from PyQt6.QtGui import QTextCursor, QColor, QTextCharFormat, QFont
from datetime import datetime
import logging


class LogViewerWidget(QWidget):
    """
    Professional log viewer with:
    - Real-time log streaming
    - Log level filtering
    - Search functionality
    - Auto-scroll
    - Color-coded entries
    """
    
    log_entry_clicked = pyqtSignal(str)
    
    def __init__(self, parent=None, max_lines: int = 1000):
        super().__init__(parent)
        self.max_lines = max_lines
        self.auto_scroll = True
        self.filter_level = "ALL"
        self.search_text = ""
        self._setup_ui()
        
        # Setup timer for refreshing
        self.refresh_timer = QTimer()
        self.refresh_timer.timeout.connect(self._refresh_display)
        self.refresh_timer.start(500)
    
    def _setup_ui(self):
        """Setup the log viewer UI"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Toolbar
        toolbar = self._create_toolbar()
        layout.addWidget(toolbar)
        
        # Log display
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setStyleSheet("""
            QTextEdit {
                background-color: #1e1e1e;
                color: #d4d4d4;
                border: 1px solid #3c3c3c;
                font-family: 'Consolas', 'Monaco', monospace;
                font-size: 11px;
            }
        """)
        layout.addWidget(self.log_text)
        
        # Status bar
        status_bar = self._create_status_bar()
        layout.addWidget(status_bar)
    
    def _create_toolbar(self) -> QFrame:
        """Create toolbar with controls"""
        toolbar = QFrame()
        toolbar.setStyleSheet("""
            QFrame {
                background-color: #252526;
                border-bottom: 1px solid #3c3c3c;
                padding: 5px;
            }
        """)
        layout = QHBoxLayout(toolbar)
        layout.setContentsMargins(10, 5, 10, 5)
        
        # Title
        title = QLabel("Console Output")
        title.setStyleSheet("""
            QLabel {
                color: #cccccc;
                font-weight: bold;
                font-size: 12px;
            }
        """)
        layout.addWidget(title)
        
        layout.addStretch()
        
        # Level filter
        self.level_combo = QComboBox()
        self.level_combo.addItems(["ALL", "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
        self.level_combo.setFixedWidth(100)
        self.level_combo.setStyleSheet("""
            QComboBox {
                background-color: #3c3c3c;
                color: #cccccc;
                border: 1px solid #555;
                padding: 3px;
                border-radius: 3px;
            }
        """)
        self.level_combo.currentTextChanged.connect(self._on_level_changed)
        layout.addWidget(QLabel("Level:"))
        layout.addWidget(self.level_combo)
        
        # Search
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("Search logs...")
        self.search_input.setFixedWidth(150)
        self.search_input.setStyleSheet("""
            QLineEdit {
                background-color: #3c3c3c;
                color: #cccccc;
                border: 1px solid #555;
                padding: 3px;
                border-radius: 3px;
            }
        """)
        self.search_input.textChanged.connect(self._on_search_changed)
        layout.addWidget(self.search_input)
        
        # Clear button
        clear_btn = QPushButton("Clear")
        clear_btn.setFixedWidth(60)
        clear_btn.setStyleSheet("""
            QPushButton {
                background-color: #0e639c;
                color: white;
                border: none;
                padding: 4px 8px;
                border-radius: 3px;
            }
            QPushButton:hover {
                background-color: #1177bb;
            }
        """)
        clear_btn.clicked.connect(self.clear)
        layout.addWidget(clear_btn)
        
        # Auto-scroll toggle
        self.auto_scroll_btn = QPushButton("Auto-scroll: ON")
        self.auto_scroll_btn.setFixedWidth(110)
        self.auto_scroll_btn.setStyleSheet("""
            QPushButton {
                background-color: #4ec9b0;
                color: #1e1e1e;
                border: none;
                padding: 4px 8px;
                border-radius: 3px;
                font-weight: bold;
            }
        """)
        self.auto_scroll_btn.clicked.connect(self._toggle_auto_scroll)
        layout.addWidget(self.auto_scroll_btn)
        
        return toolbar
    
    def _create_status_bar(self) -> QFrame:
        """Create status bar"""
        status = QFrame()
        status.setStyleSheet("""
            QFrame {
                background-color: #252526;
                border-top: 1px solid #3c3c3c;
                padding: 3px;
            }
        """)
        layout = QHBoxLayout(status)
        layout.setContentsMargins(10, 2, 10, 2)
        
        self.entry_count_label = QLabel("0 entries")
        self.entry_count_label.setStyleSheet("color: #858585; font-size: 10px;")
        layout.addWidget(self.entry_count_label)
        
        layout.addStretch()
        
        self.filter_status_label = QLabel("Filter: ALL")
        self.filter_status_label.setStyleSheet("color: #858585; font-size: 10px;")
        layout.addWidget(self.filter_status_label)
        
        return status
    
    def append_log(self, level: str, message: str, timestamp: datetime = None):
        """Add a log entry"""
        if timestamp is None:
            timestamp = datetime.now()
        
        # Apply level filter
        if self.filter_level != "ALL" and level != self.filter_level:
            return
        
        # Apply search filter
        if self.search_text and self.search_text.lower() not in message.lower():
            return
        
        # Color mapping
        colors = {
            'DEBUG': '#6a9955',
            'INFO': '#d4d4d4',
            'WARNING': '#dcdcaa',
            'ERROR': '#f14c4c',
            'CRITICAL': '#ff0000'
        }
        
        color = colors.get(level.upper(), '#d4d4d4')
        
        # Format timestamp
        ts_str = timestamp.strftime('%H:%M:%S')
        
        # Create formatted entry
        cursor = self.log_text.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        
        # Set text color
        fmt = QTextCharFormat()
        fmt.setForeground(QColor(color))
        
        cursor.insertText(f"[{ts_str}] ", fmt)
        
        # Set level color
        fmt_bold = QTextCharFormat()
        fmt_bold.setForeground(QColor(color))
        fmt_bold.setFontWeight(QFont.Weight.Bold)
        
        cursor.insertText(f"{level.upper():8s} ", fmt_bold)
        
        # Normal text
        fmt_normal = QTextCharFormat()
        fmt_normal.setForeground(QColor('#d4d4d4'))
        
        cursor.insertText(f"{message}\n", fmt_normal)
        
        # Limit lines
        self._trim_lines()
        
        # Update count
        self._update_count()
        
        # Auto-scroll
        if self.auto_scroll:
            self.log_text.moveCursor(QTextCursor.MoveOperation.End)
    
    def _trim_lines(self):
        """Trim lines if exceeding max"""
        doc = self.log_text.document()
        if doc.blockCount() > self.max_lines:
            cursor = QTextCursor(doc)
            cursor.movePosition(QTextCursor.MoveOperation.Start)
            for _ in range(doc.blockCount() - self.max_lines):
                cursor.select(QTextCursor.SelectionType.BlockUnderCursor)
                cursor.removeSelectedText()
                cursor.deleteChar()
    
    def _update_count(self):
        """Update entry count"""
        count = self.log_text.document().blockCount()
        self.entry_count_label.setText(f"{count} entries")
    
    def _on_level_changed(self, level: str):
        """Handle level filter change"""
        self.filter_level = level
        self.filter_status_label.setText(f"Filter: {level}")
    
    def _on_search_changed(self, text: str):
        """Handle search text change"""
        self.search_text = text
    
    def _toggle_auto_scroll(self):
        """Toggle auto-scroll"""
        self.auto_scroll = not self.auto_scroll
        if self.auto_scroll:
            self.auto_scroll_btn.setText("Auto-scroll: ON")
            self.auto_scroll_btn.setStyleSheet("""
                QPushButton {
                    background-color: #4ec9b0;
                    color: #1e1e1e;
                    border: none;
                    padding: 4px 8px;
                    border-radius: 3px;
                    font-weight: bold;
                }
            """)
        else:
            self.auto_scroll_btn.setText("Auto-scroll: OFF")
            self.auto_scroll_btn.setStyleSheet("""
                QPushButton {
                    background-color: #3c3c3c;
                    color: #cccccc;
                    border: 1px solid #555;
                    padding: 4px 8px;
                    border-radius: 3px;
                }
            """)
    
    def _refresh_display(self):
        """Refresh display (placeholder for external updates)"""
        pass
    
    def clear(self):
        """Clear all logs"""
        self.log_text.clear()
        self._update_count()
    
    def load_from_file(self, filepath: str):
        """Load logs from file"""
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    # Parse log line
                    parts = line.split('|')
                    if len(parts) >= 3:
                        timestamp = parts[0].strip()
                        level = parts[1].strip()
                        message = '|'.join(parts[2:]).strip()
                        
                        # Convert timestamp
                        try:
                            ts = datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
                        except:
                            ts = datetime.now()
                        
                        self.append_log(level, message, ts)
        except Exception as e:
            self.append_log('ERROR', f"Failed to load log file: {e}")


class LogViewerDockWidget(QFrame):
    """
    Docked log viewer with dock controls
    Can be attached to main window
    """
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self._setup_ui()
    
    def _setup_ui(self):
        """Setup UI"""
        self.setStyleSheet("""
            QFrame {
                background-color: #252526;
                border: 1px solid #3c3c3c;
            }
        """)
        
        layout = QVBoxLayout(self)
        
        # Header
        header = QFrame()
        header.setStyleSheet("background-color: #2d2d30;")
        header_layout = QHBoxLayout(header)
        
        title = QLabel("Console / Logs")
        title.setStyleSheet("color: #cccccc; font-weight: bold;")
        header_layout.addWidget(title)
        
        header_layout.addStretch()
        
        # Buttons
        self.toggle_btn = QPushButton("−")
        self.toggle_btn.setFixedSize(30, 20)
        self.toggle_btn.setStyleSheet("""
            QPushButton {
                background-color: #3c3c3c;
                color: #cccccc;
                border: none;
            }
            QPushButton:hover {
                background-color: #505050;
            }
        """)
        header_layout.addWidget(self.toggle_btn)
        
        layout.addWidget(header)
        
        # Log viewer
        self.log_viewer = LogViewerWidget()
        layout.addWidget(self.log_viewer)
    
    def append_log(self, level: str, message: str):
        """Add log entry"""
        self.log_viewer.append_log(level, message)


if __name__ == "__main__":
    from PyQt6.QtWidgets import QApplication
    import sys
    
    app = QApplication(sys.argv)
    
    viewer = LogViewerWidget()
    viewer.resize(800, 400)
    viewer.show()
    
    # Add sample logs
    viewer.append_log('INFO', 'Docking Studio starting...')
    viewer.append_log('DEBUG', 'Loading configuration from config.yaml')
    viewer.append_log('INFO', 'Backend connected at http://localhost:8000')
    viewer.append_log('WARNING', 'GPU not detected, using CPU mode')
    viewer.append_log('INFO', 'Loading receptor: 4xyz.pdb')
    viewer.append_log('INFO', 'Docking job started: job-123')
    viewer.append_log('ERROR', 'Failed to connect to Ollama: Connection refused')
    viewer.append_log('INFO', 'Using offline AI assistant')
    
    sys.exit(app.exec())
