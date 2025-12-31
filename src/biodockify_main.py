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
    
    Traceback:
    {traceback.format_exc()}
    """
    try:
        ctypes.windll.user32.MessageBoxW(0, debug_info[:1000], "BioDockify Boot Error", 0x10)
    except:
        print(debug_info)
    sys.exit(1)

# Debug-wrapped imports with fallback logic
try:
    try:
        from ui.main_window import MainWindow
    except ImportError:
        # Fallback: Try loading from 'src' package if root import fails
        from src.ui.main_window import MainWindow
        
    from src.config import Config
    from src.utils.log_utils import setup_logging
except ImportError as e:
    import ctypes
    import traceback
    
    # Determine base path (frozen vs source)
    if getattr(sys, 'frozen', False):
        # Fix sys._MEIPASS typo
        base_path = getattr(sys, '_MEIPASS', getattr(sys, 'MEIPASS', None))
        if not base_path:
             base_path = os.path.dirname(sys.executable)
    else:
        base_path = os.getcwd()

    # Helper: Recursive file listing limited to depth 2 to find 'ui' folder
    def list_bundle_structure(start_path):
        structure = []
        try:
            for item in sorted(os.listdir(start_path)):
                structure.append(item)
                full_item = os.path.join(start_path, item)
                if os.path.isdir(full_item) and item in ['src', 'ui']:
                    for sub in sorted(os.listdir(full_item)):
                        structure.append(f"  \\{sub}")
        except Exception as scan_err:
            structure.append(f"Scan Error: {scan_err}")
        return "\n".join(structure[:25]) # Limit output

    sys_path_str = "\n".join(str(x) for x in sys.path[:5])

    error_msg = str(e)
    traceback_str = traceback.format_exc()
    if len(traceback_str) > 5000:
        traceback_str = traceback_str[:5000] + "\n... (truncated)"

    debug_info = f"""
    Import Failure: {error_msg}
    
    Base Path: {base_path}
    
    Bundle Structure:
    {list_bundle_structure(base_path)}
    
    sys.path:
    {sys_path_str}
    
    Traceback:
    {traceback_str}
    """
    
    # Safe message box
    try:
        safe_debug = debug_info.replace('\0', '')
        ctypes.windll.user32.MessageBoxW(0, safe_debug, "BioDockify Debug Log", 0x10)
    except Exception:
        # Fallback to file
        try:
             with open("biodockify_crash.log", "w") as f:
                  f.write(debug_info)
        except: pass
        
    sys.exit(1)

from src.utils.docker_utils import check_docker_availability
from src import __version__

def main():
    """Main entry point for BioDockify Docking Studio"""
    
    # Setup logging
    setup_logging()
    
    # Create Qt Application
    app = QApplication(sys.argv)
    app.setApplicationName("BioDockify Docking Studio")
    app.setApplicationVersion(__version__)
    app.setOrganizationName("BioDockify")
    
    
    # Create main window
    window = MainWindow()
    
    # Note: Dependency checking is now handled asynchronously by the MainWindow
    
    # Show main window
    window.show()
    
    # Qt event loop
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
