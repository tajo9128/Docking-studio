"""
BioDockify Docking Studio - File Utilities
Handles file validation, hashing, and path operations
"""

import hashlib
import os
from pathlib import Path
from typing import Optional, Tuple, List
import logging
import shutil
import tempfile

logger = logging.getLogger(__name__)

class FileValidationError(Exception):
    """Base class for file validation errors"""
    pass

class FileNotFoundError(FileValidationError):
    """Receptor or Ligand file not found"""
    pass

class FilePermissionError(FileValidationError):
    """No read access to file"""
    pass

class FileCorruptionError(FileValidationError):
    """File format is invalid or corrupted"""
    pass

def validate_file(filepath: str, file_type: str = "auto") -> Tuple[bool, Optional[str]]:
    """Validate file structure and content"""
    logger.info(f"Validating file: {filepath} ({file_type})")
    
    try:
        path = Path(filepath)
        
        # Check if file exists
        if not path.exists():
            return False, f"File not found: {filepath}"
        
        # Check if file is readable
        if not os.access(filepath, os.R_OK):
            return False, f"File not readable: {filepath}"
        
        # Get file size
        file_size = path.stat().st_size
        max_size = 50 * 1024 * 1024  # 50 MB
        
        if file_size > max_size:
            return False, f"File too large: {file_size / (1024*1024):.1f} MB. Maximum: {max_size / (1024*1024):.1f} MB."
        
        # Read file content
        with open(filepath, 'r') as f:
            content = f.read()
        
        # Check if file is empty
        if len(content.strip()) == 0:
            return False, "File is empty"
        
        # Check file extension
        suffix = path.suffix.lower()
        valid_extensions = {
            "pdb": [".pdb"],
            "pdbqt": [".pdbqt"],
            "sdf": [".sdf"],
            "mol2": [".mol2"]
        }
        
        if file_type == "receptor":
            if suffix not in valid_extensions["pdb"] + valid_extensions["pdbqt"]:
                return False, f"Invalid receptor file type: {suffix}. Expected: .pdb or .pdbqt"
        elif file_type == "ligand":
            if suffix not in valid_extensions["sdf"] + valid_extensions["mol2"] + valid_extensions["pdb"]:
                return False, f"Invalid ligand file type: {suffix}. Expected: .sdf, .mol2, or .pdb"
        
        # Additional validation for PDB files
        if suffix == ".pdb":
            # Check for ATOM records
            has_ator = any("ATOM  " in line for line in content.split('\n'))
            if not has_ator:
                return False, "PDB file does not contain any ATOM records"
        
        return True, None
    
    except Exception as e:
        logger.error(f"Failed to validate file {filepath}: {e}")
        return False, f"Validation error: {str(e)}"

def calculate_file_hash(filepath: str, algorithm: str = "sha256") -> str:
    """Calculate file hash for integrity verification"""
    logger.debug(f"Calculating {algorithm} hash for: {filepath}")
    
    try:
        hash_func = {
            "md5": hashlib.md5,
            "sha1": hashlib.sha1,
            "sha256": hashlib.sha256,
            "sha512": hashlib.sha512
        }.get(algorithm, hashlib.sha256)
        
        with open(filepath, 'rb') as f:
            hash_value = hash_func(f.read()).hexdigest()
        
        logger.debug(f"Calculated {algorithm} hash for {filepath}: {hash_value}")
        return hash_value
    
    except Exception as e:
        logger.error(f"Failed to calculate hash for {filepath}: {e}")
        return ""

def get_file_info(filepath: str) -> dict:
    """Get file information"""
    logger.debug(f"Getting file info for: {filepath}")
    
    try:
        path = Path(filepath)
        
        return {
            "filename": path.name,
            "filepath": str(path.absolute()),
            "size": path.stat().st_size,
            "extension": path.suffix,
            "exists": path.exists(),
            "is_readable": os.access(filepath, os.R_OK),
            "is_writable": os.access(filepath, os.W_OK),
            "is_directory": path.is_dir()
        }
    
    except Exception as e:
        logger.error(f"Failed to get file info for {filepath}: {e}")
        return {
            "filename": Path(filepath).name,
            "filepath": str(Path(filepath).absolute()),
            "size": 0,
            "extension": Path(filepath).suffix,
            "exists": False,
            "is_readable": False,
            "is_writable": False,
            "is_directory": False
        }

def read_file_content(filepath: str, encoding: str = "utf-8") -> Optional[str]:
    """Read file content safely"""
    logger.debug(f"Reading file: {filepath}")
    
    try:
        with open(filepath, 'r', encoding=encoding) as f:
            return f.read()
    
    except UnicodeDecodeError:
        logger.error(f"Unicode decode error for {filepath}")
        return None
    
    except Exception as e:
        logger.error(f"Failed to read file {filepath}: {e}")
        return None

def write_file_content(filepath: str, content: str, encoding: str = "utf-8") -> bool:
    """Write content to file safely"""
    logger.debug(f"Writing to file: {filepath}")
    
    try:
        # Create directory if not exists
        path = Path(filepath)
        path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(filepath, 'w', encoding=encoding) as f:
            f.write(content)
        
        logger.info(f"Successfully wrote to file: {filepath}")
        return True
    
    except Exception as e:
        logger.error(f"Failed to write to file {filepath}: {e}")
        return False

def create_temp_directory(prefix: str = "BioDockify") -> str:
    """Create temporary directory"""
    import tempfile
    
    try:
        temp_dir = tempfile.mkdtemp(prefix=prefix)
        logger.info(f"Created temporary directory: {temp_dir}")
        return temp_dir
    
    except Exception as e:
        logger.error(f"Failed to create temporary directory: {e}")
        return ""

def cleanup_temp_directory(temp_dir: str) -> bool:
    """Clean up temporary directory"""
    import shutil
    
    try:
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
            logger.info(f"Cleaned up temporary directory: {temp_dir}")
            return True
    
    except Exception as e:
        logger.error(f"Failed to clean up temporary directory {temp_dir}: {e}")
        return False
