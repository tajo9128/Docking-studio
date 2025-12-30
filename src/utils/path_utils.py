"""
BioDockify Docking Studio - Path Utilities
Handles path sanitization and validation
"""

import os
import re
from pathlib import Path
from typing import Optional, Tuple, List
import logging

logger = logging.getLogger(__name__)

def sanitize_path(path: str) -> str:
    """Sanitize file path for cross-platform compatibility"""
    logger.debug(f"Sanitizing path: {path}")
    
    try:
        path_obj = Path(path)
        
        # Normalize path (remove redundant separators)
        normalized = str(path_obj.normalize())
        
        # Remove invalid characters
        sanitized = re.sub(r'[<>:"|?*]', '', normalized)
        
        # Expand user directory (~) to full path
        sanitized = os.path.expanduser(sanitized)
        
        logger.debug(f"Sanitized path: {sanitized}")
        return sanitized
    
    except Exception as e:
        logger.error(f"Failed to sanitize path {path}: {e}")
        return path

def get_safe_path(path: str, base_dir: str) -> str:
    """
    Get a safe absolute path, preventing directory traversal.
    Raises ValueError if the resulting path is outside base_dir.
    """
    # Sanitize first
    clean_path = sanitize_path(path)
    
    # Resolve absolute paths
    abs_base = os.path.abspath(base_dir)
    abs_target = os.path.abspath(os.path.join(abs_base, clean_path))
    
    # Check proper prefix
    if not abs_target.startswith(abs_base):
        logger.error(f"Path traversal attempt detected: {path} -> {abs_target} (Base: {abs_base})")
        raise ValueError("Invalid file path: Path traversal detected")
    
    return abs_target

def ensure_directory_exists(directory: str, create: bool = True) -> bool:
    """Ensure directory exists, optionally creating it"""
    logger.debug(f"Ensuring directory exists: {directory} (create: {create})")
    
    try:
        path = Path(directory)
        
        if path.exists():
            if not path.is_dir():
                logger.error(f"Path exists but is not directory: {directory}")
                return False
            
            return True
        
        elif create:
            path.mkdir(parents=True, exist_ok=True)
            logger.info(f"Created directory: {directory}")
            return True
        else:
            logger.error(f"Directory does not exist and create is False: {directory}")
            return False
    
    except Exception as e:
        logger.error(f"Failed to ensure directory exists {directory}: {e}")
        return False

def get_relative_path(absolute_path: str, base_directory: str) -> Optional[str]:
    """Get relative path from absolute path to base directory"""
    logger.debug(f"Getting relative path: {absolute_path} from {base_directory}")
    
    try:
        abs_path = Path(absolute_path).absolute()
        base_path = Path(base_directory).absolute()
        
        relative_path = abs_path.relative_to(base_path)
        logger.debug(f"Relative path: {relative_path}")
        return str(relative_path)
    
    except ValueError:
        logger.error(f"Cannot get relative path: {absolute_path} is not relative to {base_directory}")
        return None
    
    except Exception as e:
        logger.error(f"Failed to get relative path: {e}")
        return None

def get_filename_from_path(filepath: str) -> str:
    """Get filename from filepath"""
    logger.debug(f"Getting filename from: {filepath}")
    
    try:
        return Path(filepath).name
    
    except Exception as e:
        logger.error(f"Failed to get filename from {filepath}: {e}")
        return Path(filepath).name

def get_extension_from_path(filepath: str) -> str:
    """Get file extension from filepath"""
    logger.debug(f"Getting extension from: {filepath}")
    
    try:
        return Path(filepath).suffix
    
    except Exception as e:
        logger.error(f"Failed to get extension from {filepath}: {e}")
        return ""

def join_paths(*paths: str) -> str:
    """Join multiple paths using OS-appropriate separator"""
    logger.debug(f"Joining paths: {paths}")
    
    try:
        return os.path.join(*paths)
    
    except Exception as e:
        logger.error(f"Failed to join paths: {e}")
        return os.path.join(*paths) if paths else ""

def validate_path(path: str, must_exist: bool = False, 
                  must_be_file: bool = False, must_be_directory: bool = False) -> Tuple[bool, Optional[str]]:
    """Validate path meets criteria"""
    logger.debug(f"Validating path: {path} (must_exist: {must_exist}, must_be_file: {must_be_file}, must_be_directory: {must_be_directory})")
    
    try:
        path_obj = Path(path)
        
        # Check if path exists
        if must_exist and not path_obj.exists():
            return False, f"Path does not exist: {path}"
        
        # Check if path is file
        if must_be_file and not path_obj.is_file():
            return False, f"Path is not a file: {path}"
        
        # Check if path is directory
        if must_be_directory and not path_obj.is_dir():
            return False, f"Path is not a directory: {path}"
        
        # Check if path is readable
        if not os.access(path, os.R_OK):
            return False, f"Path is not readable: {path}"
        
        return True, None
    
    except Exception as e:
        logger.error(f"Failed to validate path {path}: {e}")
        return False, f"Validation error: {str(e)}"
