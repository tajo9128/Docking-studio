"""Core Validators"""

import re
import hashlib
from typing import Optional, List, Tuple, Set
from fastapi import UploadFile, HTTPException
from math import isfinite
from ..config import settings
from ..logging_config import get_logger

logger = get_logger(__name__)

class FileValidator:
    """File upload validator"""
    
    ALLOWED_EXTENSIONS = settings.ALLOWED_FILE_TYPES
    ALLOWED_MIME_TYPES = settings.MIME_TYPES
    MAX_FILE_SIZE = settings.MAX_FILE_SIZE
    
    # PDB magic numbers (first 6 characters of valid PDB records)
    PDB_MAGIC_NUMBERS = {
        "HEADER", "TITLE", "COMPND", "SOURCE", "KEYWDS",
        "EXPDTA", "AUTHOR", "REVDAT", "SPRSDE", "JRNL",
        "REMARK", "DBREF", "SEQADV", "SEQRES", "MODRES",
        "HET", "HETNAM", "FORMUL", "HELIX", "SHEET", "TURN",
        "ATOM", "HETATM", "ANISOU", "TER", "CONECT",
        "MASTER", "END", "ENDMDL"
    }
    
    # Malicious patterns to detect
    MALICIOUS_PATTERNS = [
        re.compile(r"<script", re.IGNORECASE),
        re.compile(r"javascript:", re.IGNORECASE),
        re.compile(r"<\?php", re.IGNORECASE),
        re.compile(r"<%", re.IGNORECASE),
        re.compile(r"document\.cookie", re.IGNORECASE),
        re.compile(r"window\.location", re.IGNORECASE),
        re.compile(r"eval\s*\(", re.IGNORECASE),
        re.compile(r"alert\s*\(", re.IGNORECASE),
        re.compile(r"console\.log", re.IGNORECASE),
    ]
    
    @staticmethod
    async def validate_file(file: UploadFile, filename: str) -> Tuple[bool, Optional[str]]:
        """
        Validate uploaded file
        Returns: (is_valid, error_message)
        """
        # Read file content
        file_content = await file.read()
        file_size = len(file_content)
        
        # Reset file pointer for potential re-reading
        await file.seek(0)
        
        # Check file size
        if file_size > FileValidator.MAX_FILE_SIZE:
            return False, f"File too large: {file_size / 1024 / 1024:.2f} MB (max: {FileValidator.MAX_FILE_SIZE / 1024 / 1024:.2f} MB)"
        
        # Check file extension
        extension = filename.split('.')[-1].lower() if '.' in filename else None
        if extension not in FileValidator.ALLOWED_EXTENSIONS:
            return False, f"Unsupported file type: .{extension or 'unknown'}. Supported: {', '.join(FileValidator.ALLOWED_EXTENSIONS)}"
        
        # Check MIME type (if available)
        content_type = file.content_type
        if content_type and not content_type.startswith("chemical/"):
            logger.warning(f"Suspicious MIME type: {content_type} for file: {filename}")
        
        # Check magic numbers (for PDB files)
        if extension in ["pdb", "pdbqt"]:
            # Check for valid PDB magic numbers
            lines = file_content.split(b'\n')[:10]
            valid_magic = False
            for line in lines:
                record_name = line[:6].strip().decode('utf-8', errors='ignore').upper() if len(line) >= 6 else None
                if record_name in FileValidator.PDB_MAGIC_NUMBERS:
                    valid_magic = True
                    break
            
            if not valid_magic:
                logger.error(f"Invalid PDB file magic numbers for file: {filename}")
                return False, "Invalid PDB file format: File does not start with valid PDB records (HEADER, TITLE, ATOM, etc.)"
        
        # Check for malicious patterns
        text_content = file_content.decode('utf-8', errors='ignore')
        
        for pattern in FileValidator.MALICIOUS_PATTERNS:
            if pattern.search(text_content):
                logger.error(f"Malicious pattern detected: {pattern.pattern} in file: {filename}")
                return False, "File contains potentially malicious code. Upload rejected."
        
        return True, None
    
    @staticmethod
    def sanitize_filename(filename: str) -> str:
        """
        Sanitize filename to prevent path traversal
        """
        # Remove special characters
        sanitized = ''.join(c for c in filename if c.isalnum() or c in '._-')
        
        # Limit length
        return sanitized[:100]

class ContentTypeValidator:
    """Content type validator"""
    
    ALLOWED_CONTENT_TYPES = [
        "chemical/x-pdb",
        "chemical/x-pdbqt",
        "chemical/x-mdl-sdfile",
        "chemical/x-mdl-molfile",
        "chemical/x-mmcif",
    ]
    
    @staticmethod
    def validate_content_type(content_type: str, file_extension: str) -> bool:
        """
        Validate content type
        """
        if content_type:
            expected = ContentTypeValidator.ALLOWED_CONTENT_TYPES
            return content_type in expected
        
        return True

class AtomValidator:
    """Atom data validator"""
    
    VALID_ELEMENTS = {
        'H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I',
        'Fe', 'Mg', 'Ca', 'Mn', 'Zn', 'Cu', 'Ag', 'Au',
        'Si', 'B', 'Al', 'Ti', 'V', 'Cr', 'Ni', 'Co',
        'As', 'Se', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb',
        'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Cd',
        'In', 'Sn', 'Sb', 'Te', 'Xe', 'Cs', 'Ba',
        'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu',
        'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
        'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt',
        'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra',
        'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',
        'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs',
        'Mt', 'Rg', 'Cn', 'Fl', 'Lv',
    }
    
    @staticmethod
    def validate_element(element: str) -> bool:
        """
        Validate element symbol
        """
        return element in AtomValidator.VALID_ELEMENTS
    
    @staticmethod
    def validate_coordinates(x: float, y: float, z: float) -> bool:
        """
        Validate coordinates are reasonable
        """
        # Check for NaN or Infinity
        if not (isinstance(x, (int, float)) and isinstance(y, (int, float)) and isinstance(z, (int, float))):
            return False
        
        if not (all([isfinite(coord) for coord in [x, y, z]])):
            return False
        
        # Check for reasonable range (-1000 to 1000 Angstroms)
        if any(abs(coord) > 1000 for coord in [x, y, z]):
            return False
        
        return True

class StructureValidator:
    """Structure data validator"""
    
    MIN_ATOMS = 1
    MAX_ATOMS = settings.MAX_ATOMS
    MIN_MODELS = 1
    MAX_MODELS = 100
    
    @staticmethod
    def validate_atom_count(atom_count: int) -> bool:
        """
        Validate atom count is within reasonable range
        """
        return StructureValidator.MIN_ATOMS <= atom_count <= StructureValidator.MAX_ATOMS
    
    @staticmethod
    def validate_model_count(model_count: int) -> bool:
        """
        Validate model count is within reasonable range
        """
        return StructureValidator.MIN_MODELS <= model_count <= StructureValidator.MAX_MODELS
