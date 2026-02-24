"""
Surface Engine
Molecular surface generation and rendering for 3D visualization.
"""

from enum import Enum
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple, Callable
import numpy as np
import logging

logger = logging.getLogger(__name__)


class SurfaceType(Enum):
    """Types of molecular surfaces"""
    VDW = "vdw"              # Van der Waals surface
    SAS = "sas"              # Solvent Accessible Surface
    SES = "ses"              # Solvent Excluded Surface
    MOLECULAR = "molecular"  # Molecular surface (probe accessible)


class ColoringMode(Enum):
    """Surface coloring modes"""
    ELEMENT = "element"
    CHARGE = "charge"
    HYDROPHOBICITY = "hydrophobicity"
    B_FACTOR = "b_factor"
    RESIDUE = "residue"
    ELECTROSTATIC = "electrostatic"
    CUSTOM = "custom"


@dataclass
class SurfaceConfig:
    """Configuration for surface rendering"""
    surface_type: SurfaceType = SurfaceType.VDW
    coloring_mode: ColoringMode = ColoringMode.ELEMENT
    opacity: float = 0.4
    probe_radius: float = 1.4  # Water molecule radius
    resolution: float = 1.0     # Mesh resolution
    show_contours: bool = False
    contour_level: float = 0.0


@dataclass
class SurfaceData:
    """Surface mesh data"""
    vertices: np.ndarray
    faces: np.ndarray
    colors: Optional[np.ndarray] = None
    normals: Optional[np.ndarray] = None


class SurfaceEngine:
    """
    Surface generation engine for molecular surfaces.
    
    Provides:
    - VDW surface calculation
    - Solvent Accessible Surface (SAS)
    - Solvent Excluded Surface (SES)
    - Multiple coloring schemes
    - Transparency control
    """
    
    def __init__(self, resolution: float = 1.0):
        """
        Initialize surface engine.
        
        Args:
            resolution: Mesh resolution (higher = finer mesh)
        """
        self.resolution = resolution
        self.current_surface: Optional[SurfaceData] = None
        self.config = SurfaceConfig()
        
        logger.info(f"SurfaceEngine initialized (resolution: {resolution})")
    
    def generate_surface(
        self,
        atoms: List[Dict],
        surface_type: SurfaceType = SurfaceType.VDW,
        probe_radius: float = 1.4,
    ) -> SurfaceData:
        """
        Generate molecular surface.
        
        Args:
            atoms: List of atom dictionaries with x, y, z, elem
            surface_type: Type of surface to generate
            probe_radius: Radius of solvent probe
            
        Returns:
            SurfaceData with vertices and faces
        """
        if not atoms:
            logger.warning("No atoms provided for surface generation")
            return self._empty_surface()
        
        # Extract coordinates and radii
        coords = np.array([(a['x'], a['y'], a['z']) for a in atoms])
        radii = self._get_vdw_radii(atoms)
        
        if surface_type == SurfaceType.VDW:
            surface = self._generate_vdw_surface(coords, radii)
        elif surface_type == SurfaceType.SAS:
            surface = self._generate_sas_surface(coords, radii, probe_radius)
        elif surface_type == SurfaceType.SES:
            surface = self._generate_ses_surface(coords, radii, probe_radius)
        else:
            surface = self._generate_vdw_surface(coords, radii)
        
        # Generate colors based on current coloring mode
        surface.colors = self._generate_colors(atoms, coords)
        
        self.current_surface = surface
        logger.info(f"Generated {surface_type.value} surface with {len(surface.vertices)} vertices")
        
        return surface
    
    def _get_vdw_radii(self, atoms: List[Dict]) -> np.ndarray:
        """Get Van der Waals radii for atoms"""
        # Standard VDW radii (in Angstroms)
        vdw_radii = {
            'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
            'S': 1.80, 'P': 1.80, 'F': 1.47, 'CL': 1.75,
            'BR': 1.85, 'I': 1.98, 'FE': 1.80, 'ZN': 1.39,
            'MG': 1.73, 'CA': 1.97, 'MN': 2.05, 'CO': 2.00,
        }
        
        radii = np.array([
            vdw_radii.get(atom.get('elem', 'C').upper(), 1.70)
            for atom in atoms
        ])
        
        return radii
    
    def _generate_vdw_surface(
        self,
        coords: np.ndarray,
        radii: np.ndarray
    ) -> SurfaceData:
        """Generate Van der Waals surface using sphere union"""
        # Simplified implementation - creates union of spheres
        # For production, would use MSMS or similar algorithm
        
        vertices = []
        faces = []
        
        # Sample points on each sphere
        n_theta = int(20 * self.resolution)
        n_phi = int(10 * self.resolution)
        
        for i, (center, radius) in enumerate(zip(coords, radii)):
            # Generate sphere vertices
            sphere_verts = self._generate_sphere_points(center, radius, n_theta, n_phi)
            
            # Check if vertex is not inside any other atom
            for vert in sphere_verts:
                inside = False
                for j, (other_center, other_radius) in enumerate(zip(coords, radii)):
                    if i != j:
                        dist = np.linalg.norm(vert - other_center)
                        if dist < other_radius - 0.1:
                            inside = True
                            break
                
                if not inside:
                    vertices.append(vert)
        
        if not vertices:
            return self._empty_surface()
        
        vertices = np.array(vertices)
        
        # Generate faces using convex hull (simplified)
        try:
            from scipy.spatial import ConvexHull
            hull = ConvexHull(vertices)
            faces = hull.simplices
        except Exception as e:
            logger.warning(f"Face generation failed: {e}")
            faces = np.array([[0, 1, 2]])  # Placeholder
        
        return SurfaceData(
            vertices=vertices,
            faces=faces,
            normals=self._calculate_normals(vertices, faces)
        )
    
    def _generate_sas_surface(
        self,
        coords: np.ndarray,
        radii: np.ndarray,
        probe_radius: float
    ) -> SurfaceData:
        """Generate Solvent Accessible Surface"""
        # Add probe radius to atom radii
        expanded_radii = radii + probe_radius
        return self._generate_vdw_surface(coords, expanded_radii)
    
    def _generate_ses_surface(
        self,
        coords: np.ndarray,
        radii: np.ndarray,
        probe_radius: float
    ) -> SurfaceData:
        """Generate Solvent Excluded Surface"""
        # Simplified - similar to SAS for now
        return self._generate_sas_surface(coords, radii, probe_radius)
    
    def _generate_sphere_points(
        self,
        center: np.ndarray,
        radius: float,
        n_theta: int,
        n_phi: int
    ) -> np.ndarray:
        """Generate points on a sphere surface"""
        points = []
        
        for i in range(n_theta):
            theta = np.pi * i / (n_theta - 1)
            for j in range(n_phi):
                phi = 2 * np.pi * j / (n_phi - 1)
                
                x = center[0] + radius * np.sin(theta) * np.cos(phi)
                y = center[1] + radius * np.sin(theta) * np.sin(phi)
                z = center[2] + radius * np.cos(theta)
                
                points.append(np.array([x, y, z]))
        
        return np.array(points)
    
    def _generate_colors(
        self,
        atoms: List[Dict],
        coords: np.ndarray
    ) -> np.ndarray:
        """Generate colors for surface vertices"""
        n_verts = len(coords)
        colors = np.zeros((n_verts, 3))
        
        if self.config.coloring_mode == ColoringMode.ELEMENT:
            # Element-based coloring (CPK colors)
            element_colors = {
                'H': (1.0, 1.0, 1.0),      # White
                'C': (0.5, 0.5, 0.5),      # Gray
                'N': (0.2, 0.2, 0.9),      # Blue
                'O': (0.9, 0.2, 0.2),      # Red
                'S': (1.0, 0.8, 0.2),      # Yellow
                'P': (1.0, 0.5, 0.0),      # Orange
                'FE': (0.8, 0.4, 0.0),     # Brown
                'ZN': (0.4, 0.4, 0.4),     # Slate
            }
            
            # Assign colors based on nearest atom
            for i, vert in enumerate(coords):
                distances = np.linalg.norm(coords - vert, axis=1)
                nearest_idx = np.argmin(distances)
                elem = atoms[nearest_idx].get('elem', 'C').upper()
                colors[i] = element_colors.get(elem, (0.5, 0.5, 0.5))
        
        elif self.config.coloring_mode == ColoringMode.CHARGE:
            # Charge-based coloring (red = negative, blue = positive)
            for i, vert in enumerate(coords):
                distances = np.linalg.norm(coords - vert, axis=1)
                nearest_idx = np.argmin(distances)
                elem = atoms[nearest_idx].get('elem', 'C').upper()
                
                # Simplified charge assignment
                if elem in ['O', 'S']:  # Negative
                    colors[i] = (1.0, 0.3, 0.3)
                elif elem in ['N']:  # Can be positive
                    colors[i] = (0.3, 0.3, 1.0)
                else:
                    colors[i] = (0.8, 0.8, 0.8)
        
        elif self.config.coloring_mode == ColoringMode.HYDROPHOBICITY:
            # Hydrophobicity coloring (Kyte-Doolittle scale)
            hydrophobic_colors = {
                'ILE': (1.0, 0.0, 0.0),   # Very hydrophobic (red)
                'VAL': (0.8, 0.2, 0.2),
                'LEU': (0.6, 0.2, 0.2),
                'PHE': (0.4, 0.4, 0.2),
                'CYS': (0.3, 0.5, 0.3),
                'MET': (0.3, 0.5, 0.4),
                'ALA': (0.4, 0.6, 0.4),
                'GLY': (0.5, 0.7, 0.5),
                'THR': (0.4, 0.7, 0.6),
                'SER': (0.3, 0.7, 0.7),
                'TRP': (0.3, 0.3, 0.5),
                'TYR': (0.3, 0.3, 0.5),
                'PRO': (0.3, 0.6, 0.6),
                'HIS': (0.5, 0.3, 0.6),
                'GLU': (0.0, 0.0, 1.0),   # Hydrophilic (blue)
                'GLN': (0.2, 0.2, 0.8),
                'ASP': (0.0, 0.2, 0.9),
                'ASN': (0.2, 0.3, 0.7),
                'LYS': (0.0, 0.3, 1.0),
                'ARG': (0.0, 0.3, 0.9),
            }
            
            for i, vert in enumerate(coords):
                distances = np.linalg.norm(coords - vert, axis=1)
                nearest_idx = np.argmin(distances)
                resn = atoms[nearest_idx].get('resn', 'GLY')
                colors[i] = hydrophobic_colors.get(resn, (0.5, 0.5, 0.5))
        
        else:
            # Default to white
            colors[:] = 0.9
        
        return colors
    
    def _calculate_normals(
        self,
        vertices: np.ndarray,
        faces: np.ndarray
    ) -> np.ndarray:
        """Calculate vertex normals"""
        normals = np.zeros_like(vertices)
        
        for face in faces:
            if len(face) >= 3:
                v0, v1, v2 = vertices[face[:3]]
                normal = np.cross(v1 - v0, v2 - v0)
                normals[face[0]] += normal
                normals[face[1]] += normal
                normals[face[2]] += normal
        
        # Normalize
        norms = np.linalg.norm(normals, axis=1, keepdims=True)
        norms[norms == 0] = 1
        normals = normals / norms
        
        return normals
    
    def _empty_surface(self) -> SurfaceData:
        """Return empty surface"""
        return SurfaceData(
            vertices=np.array([]),
            faces=np.array([]),
            colors=None,
            normals=None
        )
    
    def get_3dmol_config(self) -> Dict:
        """Get configuration for 3Dmol.js rendering"""
        return {
            'type': self.config.surface_type.value,
            'opacity': self.config.opacity,
            'colorScheme': self.config.coloring_mode.value,
        }
    
    def set_coloring_mode(self, mode: ColoringMode):
        """Set surface coloring mode"""
        self.config.coloring_mode = mode
    
    def set_opacity(self, opacity: float):
        """Set surface opacity (0-1)"""
        self.config.opacity = max(0.0, min(1.0, opacity))
    
    def get_vertex_colors_hex(self) -> List[str]:
        """Get colors as hex strings for 3Dmol.js"""
        if self.current_surface is None or self.current_surface.colors is None:
            return []
        
        hex_colors = []
        for color in self.current_surface.colors:
            r, g, b = int(color[0] * 255), int(color[1] * 255), int(color[2] * 255)
            hex_colors.append(f'#{r:02x}{g:02x}{b:02x}')
        
        return hex_colors


class SurfaceExporter:
    """Export surfaces to various formats"""
    
    @staticmethod
    def export_obj(surface: SurfaceData, filename: str):
        """Export surface to OBJ format"""
        with open(filename, 'w') as f:
            # Write vertices
            for v in surface.vertices:
                f.write(f"v {v[0]:.4f} {v[1]:.4f} {v[2]:.4f}\n")
            
            # Write faces (OBJ is 1-indexed)
            for face in surface.faces:
                f.write(f"f {face[0]+1} {face[1]+1} {face[2]+1}\n")
        
        logger.info(f"Surface exported to {filename}")
    
    @staticmethod
    def export_stl(surface: SurfaceData, filename: str):
        """Export surface to STL format"""
        from stl import mesh
        
        triangles = []
        for face in surface.faces:
            if len(face) >= 3:
                v0 = surface.vertices[face[0]]
                v1 = surface.vertices[face[1]]
                v2 = surface.vertices[face[2]]
                
                # Calculate normal
                normal = np.cross(v1 - v0, v2 - v0)
                
                triangles.append([v0, v1, v2])
        
        if triangles:
            stl_mesh = mesh.Mesh(np.array(triangles))
            stl_mesh.save(filename)
            logger.info(f"Surface exported to {filename}")
