"""
BioDockify Docking Studio - Main Application
Handles application startup and main window initialization
"""

from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import QTimer
from ui.main_window import MainWindow
from src.config import Config
from src.utils.log_utils import setup_logging
from src.api.dependencies import check_dependencies
import sys
import logging

def main():
    """Main entry point for BioDockify Docking Studio"""
    
    # Setup logging
    setup_logging()
    
    # Create Qt Application
    app = QApplication(sys.argv)
    app.setApplicationName("BioDockify Docking Studio")
    app.setApplicationVersion(__version__)
    app.setOrganizationName("BioDockify")
    
    # Check dependencies
    deps_available = check_dependencies()
    
    # Create main window
    window = MainWindow()
    
    # Show dependency warning if needed
    if not deps_available.get("docker"):
        logging.warning("Docker Desktop not detected. Please install Docker Desktop.")
        # Note: UI will handle this gracefully
    
    # Show main window
    window.show()
    
    # Qt event loop
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
