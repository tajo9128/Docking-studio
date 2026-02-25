"""
Agent Zero - Autonomous Docking Orchestrator & Data Custodian
Extended responsibilities:
- Conversion monitoring
- Canonical data validation
- File lifecycle management
- SQLite persistence integrity
- Storage auditing
- Recovery logic
"""

import os
import uuid
import shutil
import hashlib
from datetime import datetime
from typing import Optional, List, Dict, Any, Tuple
from pathlib import Path
from dataclasses import dataclass, field
from enum import Enum
import logging

logger = logging.getLogger(__name__)


class JobState(Enum):
    """Agent Zero job state machine."""
    IDLE = "idle"
    JOB_SUBMITTED = "job_submitted"
    VALIDATING_INPUT = "validating_input"
    CANONICALIZING = "canonicalizing"
    GRID_READY = "grid_ready"
    DOCKING = "docking"
    CNN_SCORING = "cnn_scoring"
    RF_SCORING = "rf_scoring"
    CONSENSUS = "consensus"
    PERSISTING_RESULTS = "persisting_results"
    COMPLETED = "completed"
    FAILED = "failed"


@dataclass
class ValidationResult:
    """Validation result."""
    valid: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class StorageInfo:
    """Storage information."""
    total_gb: float
    used_gb: float
    free_gb: float
    job_count: int


class AgentZeroController:
    """
    Agent Zero - Autonomous Docking Orchestrator & Data Custodian.
    
    Responsibilities:
    - Conversion monitoring (protein/ligand)
    - Canonical data validation
    - File lifecycle management
    - SQLite persistence integrity
    - Storage auditing
    - Recovery logic
    """
    
    MIN_DISK_GB = 2  # Minimum free disk space required
    MAX_LIGAND_ATOMS = 200
    MAX_FORMAL_CHARGE = 5
    
    def __init__(self, db, storage_root: str = "data/jobs"):
        self.db = db
        self.storage_root = storage_root
        self._ensure_storage_root()
    
    def _ensure_storage_root(self):
        """Ensure storage root exists."""
        Path(self.storage_root).mkdir(parents=True, exist_ok=True)
    
    # ==================== DISK SAFETY ====================
    
    def check_disk_space(self, path: str = None) -> StorageInfo:
        """Check available disk space."""
        import shutil
        path = path or self.storage_root
        stat = shutil.disk_usage(path)
        
        total = stat.total / (1024**3)
        used = stat.used / (1024**3)
        free = stat.free / (1024**3)
        
        job_count = len(list(Path(self.storage_root).glob("job_*"))) if os.path.exists(self.storage_root) else 0
        
        return StorageInfo(
            total_gb=round(total, 2),
            used_gb=round(used, 2),
            free_gb=round(free, 2),
            job_count=job_count
        )
    
    def ensure_disk_space(self, required_gb: float = None) -> bool:
        """Ensure sufficient disk space available."""
        required = required_gb or self.MIN_DISK_GB
        info = self.check_disk_space()
        
        if info.free_gb < required:
            raise RuntimeError(
                f"Insufficient disk space: {info.free_gb:.1f}GB available, "
                f"{required:.1f}GB required"
            )
        return True
    
    # ==================== RECEPTOR VALIDATION ====================
    
    def validate_receptor(self, pdbqt_path: str) -> ValidationResult:
        """
        Validate receptor PDBQT file.
        
        Checks:
        - ATOM or HETATM records exist
        - No zero-length structure
        - PDBQT contains ROOT section
        - No NaN coordinates
        """
        errors = []
        warnings = []
        metadata = {}
        
        if not os.path.exists(pdbqt_path):
            errors.append(f"Receptor file not found: {pdbqt_path}")
            return ValidationResult(valid=False, errors=errors)
        
        file_size = os.path.getsize(pdbqt_path)
        if file_size == 0:
            errors.append("Receptor file is empty")
            return ValidationResult(valid=False, errors=errors)
        
        try:
            with open(pdbqt_path, 'r') as f:
                content = f.read()
            
            # Check for ATOM/HETATM records
            atom_count = content.count("ATOM")
            hetatm_count = content.count("HETATM")
            total_atoms = atom_count + hetatm_count
            
            if total_atoms == 0:
                errors.append("No ATOM/HETATM records found in receptor")
            
            metadata["atom_count"] = total_atoms
            
            # Check ROOT section
            if "ROOT" not in content:
                warnings.append("No ROOT section found - may not be valid PDBQT")
            
            # Check for NaN
            if "nan" in content.lower():
                errors.append("Receptor contains NaN coordinates")
            
            # Check for valid coordinates
            has_valid_coords = False
            for line in content.split('\n'):
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        if not (x == 0 and y == 0 and z == 0):
                            has_valid_coords = True
                    except (ValueError, IndexError):
                        continue
            
            if not has_valid_coords and total_atoms > 0:
                warnings.append("All coordinates appear to be at origin")
            
            # Check for hydrogens
            hydrogen_count = content.count(" H ")
            metadata["hydrogen_count"] = hydrogen_count
            
        except Exception as e:
            errors.append(f"Error reading receptor: {e}")
        
        valid = len(errors) == 0
        return ValidationResult(
            valid=valid,
            errors=errors,
            warnings=warnings,
            metadata=metadata
        )
    
    # ==================== LIGAND VALIDATION ====================
    
    def validate_ligand_smiles(self, smiles: str) -> ValidationResult:
        """Validate SMILES string."""
        errors = []
        warnings = []
        
        if not smiles or not smiles.strip():
            errors.append("Empty SMILES string")
            return ValidationResult(valid=False, errors=errors)
        
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                errors.append("Invalid SMILES - could not parse")
                return ValidationResult(valid=False, errors=errors)
            
            atom_count = mol.GetNumAtoms()
            if atom_count == 0:
                errors.append("No atoms in ligand")
            elif atom_count > self.MAX_LIGAND_ATOMS:
                warnings.append(f"Large ligand: {atom_count} atoms (may be slow)")
            
            # Check formal charge
            total_charge = sum(
                atom.GetFormalCharge() for atom in mol.GetAtoms()
            )
            if abs(total_charge) > self.MAX_FORMAL_CHARGE:
                errors.append(
                    f"Extreme formal charge: {total_charge} "
                    f"(max: ±{self.MAX_FORMAL_CHARGE})"
                )
            
            return ValidationResult(
                valid=len(errors) == 0,
                errors=errors,
                warnings=warnings,
                metadata={"atom_count": atom_count}
            )
            
        except ImportError:
            warnings.append("RDKit not available - basic validation only")
            return ValidationResult(valid=True, warnings=warnings)
        except Exception as e:
            errors.append(f"Validation error: {e}")
            return ValidationResult(valid=False, errors=errors)
    
    def validate_ligand_file(self, file_path: str) -> ValidationResult:
        """Validate ligand file."""
        errors = []
        warnings = []
        metadata = {}
        
        if not os.path.exists(file_path):
            errors.append(f"Ligand file not found: {file_path}")
            return ValidationResult(valid=False, errors=errors)
        
        file_size = os.path.getsize(file_path)
        if file_size == 0:
            errors.append("Ligand file is empty")
            return ValidationResult(valid=False, errors=errors)
        
        ext = Path(file_path).suffix.lower()
        
        if ext == ".pdbqt":
            # Check PDBQT validity
            with open(file_path, 'r') as f:
                content = f.read()
            
            if "ROOT" not in content:
                warnings.append("No ROOT section - may not be valid PDBQT")
            
            atom_count = content.count("ATOM") + content.count("HETATM")
            if atom_count == 0:
                errors.append("No atoms in ligand file")
            
            metadata["atom_count"] = atom_count
        
        elif ext in [".sdf", ".mol"]:
            # Try RDKit validation
            try:
                from rdkit import Chem
                supplier = Chem.SDMolSupplier(file_path)
                mol = next(supplier)
                if mol is None:
                    errors.append("Could not parse SDF file")
                else:
                    metadata["atom_count"] = mol.GetNumAtoms()
            except ImportError:
                warnings.append("RDKit not available")
            except Exception as e:
                errors.append(f"SDF parsing error: {e}")
        
        valid = len(errors) == 0
        return ValidationResult(
            valid=valid,
            errors=errors,
            warnings=warnings,
            metadata=metadata
        )
    
    # ==================== CANONICAL DIRECTORY ====================
    
    def validate_canonical_structure(self, job_dir: str) -> ValidationResult:
        """
        Validate canonical directory structure.
        
        Expected:
        job_dir/
          canonical/
            receptor.pdbqt
            ligands/
              *.pdbqt
            ligands.sdf
            ligands.smiles
            gridbox.json
        """
        errors = []
        warnings = []
        
        canonical = os.path.join(job_dir, "canonical")
        
        if not os.path.exists(canonical):
            errors.append(f"Canonical directory not found: {canonical}")
            return ValidationResult(valid=False, errors=errors)
        
        # Check receptor
        receptor = os.path.join(canonical, "receptor.pdbqt")
        if not os.path.exists(receptor):
            errors.append("receptor.pdbqt not found")
        
        # Check ligands directory
        ligands_dir = os.path.join(canonical, "ligands")
        if not os.path.exists(ligands_dir):
            errors.append("ligands/ directory not found")
        else:
            ligand_files = list(Path(ligands_dir).glob("*.pdbqt"))
            if not ligand_files:
                errors.append("No ligand PDBQT files found")
            metadata["ligand_count"] = len(ligand_files)
        
        # Check SDF
        sdf_file = os.path.join(canonical, "ligands.sdf")
        if not os.path.exists(sdf_file):
            warnings.append("ligands.sdf not found")
        
        # Check SMILES
        smiles_file = os.path.join(canonical, "ligands.smiles")
        if not os.path.exists(smiles_file):
            warnings.append("ligands.smiles not found")
        
        # Check gridbox
        gridbox = os.path.join(canonical, "gridbox.json")
        if not os.path.exists(gridbox):
            warnings.append("gridbox.json not found")
        
        valid = len(errors) == 0
        return ValidationResult(
            valid=valid,
            errors=errors,
            warnings=warnings,
            metadata=metadata if 'metadata' in locals() else {}
        )
    
    # ==================== LOG SUPERVISION ====================
    
    def supervise_vina_log(self, log_path: str) -> ValidationResult:
        """Supervise Vina log file."""
        errors = []
        warnings = []
        
        if not os.path.exists(log_path):
            errors.append("Vina log file not found")
            return ValidationResult(valid=False, errors=errors)
        
        file_size = os.path.getsize(log_path)
        if file_size == 0:
            errors.append("Vina log is empty")
            return ValidationResult(valid=False, errors=errors)
        
        try:
            with open(log_path, 'r') as f:
                content = f.read()
            
            # Check for errors
            if "error" in content.lower() and "Affinity" not in content:
                errors.append("Vina log contains errors")
            
            # Check for success
            if "mode |" not in content:
                warnings.append("No docking results found in log")
            
            # Parse runtime
            runtime = None
            for line in content.split('\n'):
                if "Total runtime" in line:
                    try:
                        runtime = float(line.split()[2])
                    except (IndexError, ValueError):
                        pass
            
            metadata = {"runtime_seconds": runtime}
            
        except Exception as e:
            errors.append(f"Error reading log: {e}")
        
        valid = len(errors) == 0
        return ValidationResult(
            valid=valid,
            errors=errors,
            warnings=warnings,
            metadata=metadata if 'metadata' in locals() else {}
        )
    
    def supervise_gnina_log(self, log_path: str) -> ValidationResult:
        """Supervise GNINA log file."""
        errors = []
        warnings = []
        
        if not os.path.exists(log_path):
            # GNINA log may not exist in CPU mode
            warnings.append("GNINA log not found (CPU mode?)")
            return ValidationResult(valid=True, warnings=warnings)
        
        try:
            with open(log_path, 'r') as f:
                content = f.read()
            
            # Check for crashes
            if "segmentation fault" in content.lower():
                errors.append("GNINA segmentation fault detected")
            
            if "CUDA" in content and "error" in content.lower():
                warnings.append("CUDA errors in GNINA log")
            
            # Check for results
            has_affinity = "Affinity:" in content or "CNNaffinity" in content
            if not has_affinity:
                warnings.append("No CNN affinity scores found")
            
            metadata = {"has_results": has_affinity}
            
        except Exception as e:
            errors.append(f"Error reading GNINA log: {e}")
        
        valid = len(errors) == 0
        return ValidationResult(
            valid=valid,
            errors=errors,
            warnings=warnings,
            metadata=metadata if 'metadata' in locals() else {}
        )
    
    # ==================== CONSENSUS INTEGRITY ====================
    
    def validate_consensus_data(
        self,
        vina_score: Optional[float],
        rf_score: Optional[float],
        gnina_score: Optional[float],
        gpu_used: bool
    ) -> ValidationResult:
        """
        Validate consensus data completeness.
        
        GPU mode: vina + gnina + rf required
        CPU mode: vina + rf required
        """
        errors = []
        warnings = []
        
        # Check Vina (always required)
        if vina_score is None:
            errors.append("Missing Vina score")
        elif vina_score > 0:
            warnings.append("Vina score should be negative (binding energy)")
        
        # Check RF (always required)
        if rf_score is None:
            errors.append("Missing RF score")
        
        # Check GNINA (GPU only)
        if gpu_used and gnina_score is None:
            warnings.append("Missing GNINA score in GPU mode")
        
        valid = len(errors) == 0
        return ValidationResult(
            valid=valid,
            errors=errors,
            warnings=warnings
        )
    
    # ==================== FILE MANAGEMENT ====================
    
    def get_storage_path(self, user_id: int, job_uuid: str) -> str:
        """Get job storage path."""
        return os.path.join(
            self.storage_root,
            f"user_{user_id}",
            f"job_{job_uuid}"
        )
    
    def sanitize_filename(self, filename: str) -> str:
        """Sanitize filename for safety."""
        # Remove path separators
        filename = filename.replace("/", "_").replace("\\", "_")
        # Remove dangerous characters
        for char in ["..", "$", "`", "|", ";", "&", "<", ">"]:
            filename = filename.replace(char, "_")
        # Limit length
        if len(filename) > 255:
            name, ext = os.path.splitext(filename)
            filename = name[:255-len(ext)] + ext
        return filename
    
    def compute_file_hash(self, file_path: str) -> str:
        """Compute SHA256 hash of file."""
        sha256 = hashlib.sha256()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                sha256.update(chunk)
        return sha256.hexdigest()
    
    # ==================== JOB STATE MANAGEMENT ====================
    
    def create_job_with_state(
        self,
        user_id: int,
        job_name: str,
        compute_mode: str = "auto",
        gpu_available: bool = False,
        gpu_name: str = None
    ) -> str:
        """Create job and set initial state."""
        # Check disk space first
        self.ensure_disk_space()
        
        # Create job in database
        job_uuid = self.db.create_job(
            user_id=user_id,
            job_name=job_name,
            compute_mode=compute_mode,
            gpu_used=gpu_available,
            gpu_name=gpu_name
        )
        
        # Create storage directory
        storage_path = self.get_storage_path(user_id, job_uuid)
        os.makedirs(storage_path, exist_ok=True)
        
        # Update state
        self.db.update_job_status(job_uuid, JobState.JOB_SUBMITTED.value)
        
        return job_uuid
    
    def transition_state(self, job_uuid: str, new_state: JobState):
        """Transition job to new state."""
        state_str = new_state.value
        logger.info(f"Job {job_uuid}: {state_str}")
        self.db.update_job_status(job_uuid, state_str)
        
        # Also update progress based on state
        progress_map = {
            JobState.JOB_SUBMITTED: 5,
            JobState.VALIDATING_INPUT: 10,
            JobState.CANONICALIZING: 20,
            JobState.GRID_READY: 30,
            JobState.DOCKING: 50,
            JobState.CNN_SCORING: 70,
            JobState.RF_SCORING: 80,
            JobState.CONSENSUS: 90,
            JobState.PERSISTING_RESULTS: 95,
            JobState.COMPLETED: 100,
            JobState.FAILED: 0,
        }
        
        if new_state in progress_map:
            try:
                self.db.update_progress(job_uuid, progress_map[new_state])
            except:
                pass
    
    def mark_job_failed(self, job_uuid: str, error: str):
        """Mark job as failed with error."""
        logger.error(f"Job {job_uuid} failed: {error}")
        self.db.update_job_status(job_uuid, JobState.FAILED.value, error_message=error)
        
        try:
            self.db.update_progress(job_uuid, 0)
        except:
            pass
    
    # ==================== STORAGE AUDITING ====================
    
    def audit_storage(self) -> Dict[str, Any]:
        """Audit storage and check for issues."""
        storage_info = self.check_disk_space()
        
        issues = []
        
        # Check disk space
        if storage_info.free_gb < self.MIN_DISK_GB:
            issues.append(f"Low disk space: {storage_info.free_gb:.1f}GB")
        
        # Check for orphaned directories
        if not os.path.exists(self.storage_root):
            return {"status": "ok", "issues": issues}
        
        # List jobs with missing files
        job_dirs = list(Path(self.storage_root).glob("*/job_*"))
        
        incomplete_jobs = []
        for job_dir in job_dirs:
            canonical = job_dir / "canonical"
            if not canonical.exists():
                incomplete_jobs.append(str(job_dir.name))
        
        if incomplete_jobs:
            warnings = [f"Incomplete jobs: {len(incomplete_jobs)}"]
        
        return {
            "status": "warning" if issues else "ok",
            "storage": {
                "total_gb": storage_info.total_gb,
                "used_gb": storage_info.used_gb,
                "free_gb": storage_info.free_gb,
                "job_count": storage_info.job_count
            },
            "issues": issues,
            "incomplete_jobs": incomplete_jobs if 'incomplete_jobs' in locals() else []
        }
    
    # ==================== RESULTS PERSISTENCE ====================
    
    def persist_results(
        self,
        job_uuid: str,
        results: List[Dict[str, Any]],
        grid_config: Dict[str, Any]
    ):
        """Persist docking results to database."""
        # Add grid metadata
        self.db.add_grid_metadata(
            job_uuid=job_uuid,
            center_x=grid_config.get("center_x", 0),
            center_y=grid_config.get("center_y", 0),
            center_z=grid_config.get("center_z", 0),
            size_x=grid_config.get("size_x", 22),
            size_y=grid_config.get("size_y", 22),
            size_z=grid_config.get("size_z", 22),
            exhaustiveness=grid_config.get("exhaustiveness", 8),
            seed=grid_config.get("seed", 123456)
        )
        
        # Add results
        self.db.add_results_batch(job_uuid, results)
        
        # Transition to completed
        self.transition_state(job_uuid, JobState.COMPLETED)


# ==================== HELPER FUNCTIONS ====================

def create_agent_zero(db, storage_root: str = "data/jobs") -> AgentZeroController:
    """Create Agent Zero controller instance."""
    return AgentZeroController(db, storage_root)
