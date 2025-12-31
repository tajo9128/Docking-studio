"""
BioDockify Docking Studio - Main Application
Handles application startup and main window initialization
"""

import sys
import os
import logging

# Early debug hook for import errors
try:
    from PyQt6.QtWidgets import QApplication
    from PyQt6.QtCore import QTimer
except ImportError as e:
    import ctypes
    import traceback
    
    # helper to format list
    def fmt(l): return "\n".join(str(x) for x in l[:10])
    
    debug_info = f"""
    CRITICAL IMPORT ERROR: {e}
    
    Python Ver: {sys.version}
    
    sys.path:
    {fmt(sys.path)}
    
    Files in local dir:
    {fmt(os.listdir('.'))}
    
    Traceback:
    {traceback.format_exc()}
    """
    ctypes.windll.user32.MessageBoxW(0, debug_info, "BioDockify Boot Error", 0x10)
    sys.exit(1)

from ui.main_window import MainWindow
from src.config import Config
from src.utils.log_utils import setup_logging
from src.api.dependencies import check_dependencies

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
