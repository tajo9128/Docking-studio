"""
Structured Logging System for Docking Studio
Production-grade logging with file rotation, levels, and formatting
"""

import logging
import logging.handlers
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional
from enum import Enum


class LogLevel(Enum):
    """Log level enumeration"""
    DEBUG = "DEBUG"
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"
    CRITICAL = "CRITICAL"


class StructuredLogger:
    """
    Production logging wrapper with:
    - File rotation
    - Multiple handlers
    - Structured JSON output option
    - Colorized console output
    """
    
    def __init__(
        self,
        name: str = "docking-studio",
        log_dir: str = "logs",
        level: str = "INFO",
        max_bytes: int = 10 * 1024 * 1024,  # 10MB
        backup_count: int = 5
    ):
        self.name = name
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        
        self.logger = logging.getLogger(name)
        self.logger.setLevel(getattr(logging, level.upper()))
        self.logger.handlers = []
        
        self.max_bytes = max_bytes
        self.backup_count = backup_count
        
        self._setup_handlers()
    
    def _setup_handlers(self):
        """Setup logging handlers"""
        
        # File handler with rotation
        log_file = self.log_dir / f"{self.name}.log"
        file_handler = logging.handlers.RotatingFileHandler(
            log_file,
            maxBytes=self.max_bytes,
            backupCount=self.backup_count,
            encoding='utf-8'
        )
        file_formatter = logging.Formatter(
            '%(asctime)s | %(levelname)-8s | %(name)s | %(funcName)s:%(lineno)d | %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(file_formatter)
        self.logger.addHandler(file_handler)
        
        # Error-only file handler
        error_file = self.log_dir / f"{self.name}_errors.log"
        error_handler = logging.handlers.RotatingFileHandler(
            error_file,
            maxBytes=self.max_bytes,
            backupCount=self.backup_count,
            encoding='utf-8'
        )
        error_handler.setLevel(logging.ERROR)
        error_handler.setFormatter(file_formatter)
        self.logger.addHandler(error_handler)
        
        # Console handler with colors
        console_handler = logging.StreamHandler(sys.stdout)
        console_formatter = ColoredFormatter(
            '%(asctime)s | %(levelname)-8s | %(message)s',
            datefmt='%H:%M:%S'
        )
        console_handler.setFormatter(console_formatter)
        self.logger.addHandler(console_handler)
    
    def debug(self, message: str, **kwargs):
        """Log debug message"""
        self.logger.debug(self._format_message(message, **kwargs))
    
    def info(self, message: str, **kwargs):
        """Log info message"""
        self.logger.info(self._format_message(message, **kwargs))
    
    def warning(self, message: str, **kwargs):
        """Log warning message"""
        self.logger.warning(self._format_message(message, **kwargs))
    
    def error(self, message: str, **kwargs):
        """Log error message"""
        self.logger.error(self._format_message(message, **kwargs))
    
    def critical(self, message: str, **kwargs):
        """Log critical message"""
        self.logger.critical(self._format_message(message, **kwargs))
    
    def _format_message(self, message: str, **kwargs) -> str:
        """Format message with additional context"""
        if kwargs:
            context = " | ".join(f"{k}={v}" for k, v in kwargs.items())
            return f"{message} | {context}"
        return message
    
    def log_docking_job(self, job_id: str, status: str, **details):
        """Log docking job event with structured data"""
        self.info(f"JOB [{job_id}] {status}", **details)
    
    def log_security_event(self, event: str, severity: str, **details):
        """Log security event"""
        log_func = getattr(self, severity.lower(), self.info)
        log_func(f"SECURITY: {event}", **details)
    
    def log_ai_interaction(self, provider: str, message: str, success: bool):
        """Log AI interaction"""
        status = "SUCCESS" if success else "FAILED"
        self.info(f"AI [{provider}] {status}: {message[:100]}")
    
    def log_docker_event(self, event: str, container: str, **details):
        """Log Docker event"""
        self.info(f"DOCKER [{container}] {event}", **details)
    
    def get_log_files(self) -> list:
        """Get list of log files"""
        return list(self.log_dir.glob("*.log"))


class ColoredFormatter(logging.Formatter):
    """Colored console formatter"""
    
    COLORS = {
        'DEBUG': '\033[36m',     # Cyan
        'INFO': '\033[32m',      # Green
        'WARNING': '\033[33m',    # Yellow
        'ERROR': '\033[31m',     # Red
        'CRITICAL': '\033[35m',  # Magenta
        'RESET': '\033[0m'
    }
    
    def format(self, record):
        levelname = record.levelname
        if levelname in self.COLORS:
            record.levelname = f"{self.COLORS[levelname]}{levelname}{self.COLORS['RESET']}"
        return super().format(record)


def get_logger(name: str = "docking-studio") -> StructuredLogger:
    """Get or create logger instance"""
    if not hasattr(get_logger, "_instance"):
        get_logger._instance = StructuredLogger(name)
    return get_logger


class LogCapture:
    """Capture logs for display in UI"""
    
    def __init__(self, max_entries: int = 1000):
        self.max_entries = max_entries
        self.entries = []
        self.callbacks = []
    
    def add_entry(self, level: str, message: str, timestamp: datetime = None):
        """Add log entry"""
        if timestamp is None:
            timestamp = datetime.now()
        
        entry = {
            'timestamp': timestamp,
            'level': level,
            'message': message
        }
        
        self.entries.append(entry)
        
        if len(self.entries) > self.max_entries:
            self.entries.pop(0)
        
        for callback in self.callbacks:
            callback(entry)
    
    def register_callback(self, callback):
        """Register callback for new entries"""
        self.callbacks.append(callback)
    
    def get_recent(self, count: int = 100) -> list:
        """Get recent entries"""
        return self.entries[-count:]
    
    def clear(self):
        """Clear all entries"""
        self.entries.clear()


# Global log capture for UI
log_capture = LogCapture()


class UIHandler(logging.Handler):
    """Custom handler that sends logs to UI"""
    
    def __init__(self, log_capture: LogCapture):
        super().__init__()
        self.log_capture = log_capture
    
    def emit(self, record):
        try:
            msg = self.format(record)
            self.log_capture.add_entry(
                level=record.levelname,
                message=msg,
                timestamp=datetime.fromtimestamp(record.created)
            )
        except Exception:
            pass


def setup_logging(level: str = "INFO"):
    """Setup logging for the application"""
    logger = get_logger()
    logger.logger.setLevel(getattr(logging, level.upper()))
    
    # Add UI handler
    ui_handler = UIHandler(log_capture)
    ui_handler.setFormatter(logging.Formatter('%(levelname)s | %(message)s'))
    logger.logger.addHandler(ui_handler)
    
    return logger


if __name__ == "__main__":
    logger = get_logger()
    
    logger.info("Docking Studio starting...")
    logger.debug("Debug message")
    logger.warning("Warning message")
    logger.error("Error message")
    
    logger.log_docking_job("job-123", "STARTED", receptor="4xyz.pdb", ligand="drug.sdf")
    logger.log_security_event("scan_complete", "INFO", issues=0)
    logger.log_ai_interaction("ollama", "What is vina?", True)
    
    print(f"\nLog files: {logger.get_log_files()}")
