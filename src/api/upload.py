"""
BioDockify Docking Studio - File Upload API Router
Handles file uploads and validation
"""

from fastapi import APIRouter, UploadFile, File, HTTPException, status
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from typing import Optional
import logging
from datetime import datetime
import hashlib
import os

from src.utils.file_utils import validate_file, calculate_file_hash, get_file_info
from src.utils.path_utils import sanitize_path, get_safe_path, ensure_directory_exists

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/upload")

class UploadResponse(BaseModel):
    """File upload response model"""
    filename: str
    original_filename: str
    file_size: int
    file_hash: str
    mime_type: Optional[str]
    upload_status: str
    validation_status: str
    validation_errors: Optional[list] = None
    timestamp: str

@router.post("/receptor", response_model=UploadResponse)
async def upload_receptor(
    file: UploadFile = File(...),
    content_type: Optional[str] = None
):
    """Upload receptor file"""
    
    # Read file content
    content = await file.read()
    
    # Base directory for uploads
    base_upload_dir = "uploads/receptors"
    
    # Ensure directory exist using os.makedirs for safety
    try:
        os.makedirs(base_upload_dir, exist_ok=True)
    except OSError as e:
         logger.error(f"Failed to create directory {base_upload_dir}: {e}")
         raise HTTPException(status_code=500, detail="Internal Server Error: Could not create upload directory")

    # Use secure path validation (Fix #2)
    try:
        # We allow the user supplied filename, but we strictly validate it lands in base_upload_dir
        file_path = get_safe_path(file.filename, base_upload_dir)
        sanitized_filename = os.path.basename(file_path) # Extracted just for response
    except ValueError as e:
        logger.warning(f"Invalid upload path: {file.filename}")
        raise HTTPException(status_code=400, detail="Invalid filename or path traversal detected")
    
    # Write content
    try:
        with open(file_path, 'wb') as f:
            f.write(content)
    except Exception as e:
        logger.error(f"Failed to write file {file_path}: {e}")
        raise HTTPException(status_code=500, detail="Failed to save uploaded file")
        
    # Validate file content
    is_valid, error_message = validate_file(file_path, "receptor")
    
    # Calculate file hash
    file_hash = calculate_file_hash(file_path)
    
    # Get file info
    file_info = get_file_info(file_path)
    
    response = UploadResponse(
        filename=sanitized_filename,
        original_filename=file.filename,
        file_size=file_info["size"],
        file_hash=file_hash,
        mime_type=content_type,
        upload_status="uploaded",
        validation_status="valid" if is_valid else "invalid",
        validation_errors=[error_message] if not is_valid else None,
        timestamp=datetime.now().isoformat()
    )
    
    if is_valid:
        logger.info(f"Receptor uploaded successfully: {sanitized_filename}")
    else:
        logger.warning(f"Receptor upload failed: {sanitized_filename} - {error_message}")
    
    return response

@router.post("/ligand", response_model=UploadResponse)
async def upload_ligand(
    file: UploadFile = File(...),
    content_type: Optional[str] = None
):
    """Upload ligand file"""
    
    content = await file.read()
    
    base_upload_dir = "uploads/ligands"
    os.makedirs(base_upload_dir, exist_ok=True)

    try:
        file_path = get_safe_path(file.filename, base_upload_dir)
        sanitized_filename = os.path.basename(file_path)
    except ValueError:
        raise HTTPException(status_code=400, detail="Invalid filename or path traversal detected")
    
    try:
        with open(file_path, 'wb') as f:
            f.write(content)
    except Exception as e:
        logger.error(f"Write error: {e}")
        raise HTTPException(status_code=500, detail="Failed to save file")

    is_valid, error_message = validate_file(file_path, "ligand")
    file_hash = calculate_file_hash(file_path)
    file_info = get_file_info(file_path)
    
    response = UploadResponse(
        filename=sanitized_filename,
        original_filename=file.filename,
        file_size=file_info["size"],
        file_hash=file_hash,
        mime_type=content_type,
        upload_status="uploaded",
        validation_status="valid" if is_valid else "invalid",
        validation_errors=[error_message] if not is_valid else None,
        timestamp=datetime.now().isoformat()
    )
    
    if is_valid:
        logger.info(f"Ligand uploaded successfully: {sanitized_filename}")
    else:
        logger.warning(f"Ligand upload failed: {sanitized_filename} - {error_message}")
    
    return response
