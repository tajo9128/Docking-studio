"""
BioDockify Docking Studio - Logging Utilities
Handles logging configuration and logger setup
"""

import logging
import os
import sys
from pathlib import Path
from typing import Optional

def setup_logging(log_level: str = "INFO", log_to_file: bool = True, 
               log_directory: Optional[str] = None) -> logging.Logger:
    """Setup logging configuration"""
    
    # Get log directory
    if log_directory is None:
        if os.name == "nt":  # Windows
            log_dir = Path.home() / "AppData" / "Local" / "BioDockify" / "logs"
        else:  # macOS/Linux
            log_dir = Path.home() / ".config" / "BioDockify" / "logs"
    else:
        log_dir = Path(log_directory)
    
    # Create log directory
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # Configure root logger
    root_logger = logging.getLogger()
    
    # Set log level
    numeric_level = getattr(logging, log_level.upper(), logging.INFO)
    root_logger.setLevel(numeric_level)
    
    # Remove existing handlers
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(numeric_level)
    console_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    root_logger.addHandler(console_handler)
    
    # File handler
    if log_to_file:
        log_file = log_dir / "BioDockify.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(numeric_level)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
        root_logger.addHandler(file_handler)
    
    return root_logger

def get_logger(name: str) -> logging.Logger:
    """Get or create logger with specific name"""
    return logging.getLogger(name)

def set_log_level(level: str) -> None:
    """Set logging level"""
    numeric_level = getattr(logging, level.upper(), logging.INFO)
    logging.getLogger().setLevel(numeric_level)

def log_exception(logger: logging.Logger, exception: Exception, 
                  message: str = "An exception occurred") -> None:
    """Log exception with traceback"""
    logger.error(f"{message}: {str(exception)}", exc_info=True)

def log_function_call(logger: logging.Logger, function_name: str, 
                    args: tuple, kwargs: dict) -> None:
    """Log function call details"""
    logger.debug(f"Calling function: {function_name} with args={args}, kwargs={kwargs}")

def export_logs(job_id: str, log_directory: Optional[str] = None) -> Optional[str]:
    """Export logs for a specific job"""
    logger.info(f"Exporting logs for job: {job_id}")
    
    if log_directory is None:
        if os.name == "nt":  # Windows
            log_dir = Path.home() / "AppData" / "Local" / "BioDockify" / "logs"
        else:  # macOS/Linux
            log_dir = Path.home() / ".config" / "BioDockify" / "logs"
    else:
        log_dir = Path(log_directory)
    
    try:
        # Find log files for specific job
        job_log_file = None
        for file in log_dir.glob("*.log"):
            if job_id in file.name:
                job_log_file = file
                break
        
        if not job_log_file:
            return None
        
        # Read log file content
        with open(job_log_file, 'r', encoding='utf-8') as f:
            return f.read()
    
    except Exception as e:
        logger.error(f"Failed to export logs for job {job_id}: {e}")
        return None
