"""Core Utilities"""

import hashlib
import uuid
from typing import Any
from datetime import datetime
import logging

logger = logging.getLogger(__name__)

def calculate_hash(content: bytes) -> str:
    """Calculate SHA-256 hash of content"""
    return hashlib.sha256(content).hexdigest()

def generate_correlation_id() -> str:
    """Generate correlation ID for request tracing"""
    return str(uuid.uuid4())

def format_size(bytes_val: int) -> str:
    """Format bytes to human readable size"""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if bytes_val < 1024.0:
            return f"{bytes_val:.2f} {unit}"
        bytes_val /= 1024.0
    return f"{bytes_val:.2f} GB"

def get_current_time_ms() -> int:
    """Get current time in milliseconds"""
    return int(datetime.now().timestamp() * 1000)

class PerformanceTimer:
    """Context manager for timing operations"""
    
    def __init__(self, operation_name: str):
        self.operation_name = operation_name
        self.start_time = None
        self.end_time = None
    
    def __enter__(self):
        self.start_time = get_current_time_ms()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end_time = get_current_time_ms()
        duration = self.end_time - self.start_time
        
        logger.info(
            f"Operation: {self.operation_name} completed in {duration}ms",
            extra={
                "operation": self.operation_name,
                "duration_ms": duration,
            }
        )
        
        return False
