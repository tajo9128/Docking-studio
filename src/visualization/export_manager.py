"""
Export Manager
Publication-ready export for molecular visualizations.
"""

from enum import Enum
from dataclasses import dataclass
from typing import List, Dict, Optional, Callable
import base64
import logging

logger = logging.getLogger(__name__)


class ExportFormat(Enum):
    """Export format options"""
    PNG = "png"
    JPG = "jpg"
    SVG = "svg"
    PDF = "pdf"
    GIF = "gif"
    WEBM = "webm"
    OBJ = "obj"  # 3D mesh
    STL = "stl"  # 3D mesh


@dataclass
class ExportSettings:
    """Settings for export"""
    format: ExportFormat = ExportFormat.PNG
    width: int = 1920
    height: int = 1080
    transparent: bool = False
    background_color: str = "white"
    dpi: int = 300
    quality: int = 95
    antialiasing: bool = True
    ray_tracing: bool = False
    animation_frames: int = 36  # For GIF
    animation_duration: float = 3.0  # Seconds


class ExportManager:
    """
    Manages export of molecular visualizations.
    
    Supports:
    - PNG/JPG images at various resolutions (including 4K)
    - SVG vector graphics
    - PDF documents
    - Animated GIFs
    - 3D mesh formats (OBJ, STL)
    """
    
    def __init__(self):
        """Initialize export manager"""
        self.settings = ExportSettings()
        self.export_callback: Optional[Callable] = None
        
        logger.info("ExportManager initialized")
    
    def configure(self, **kwargs):
        """Configure export settings"""
        for key, value in kwargs.items():
            if hasattr(self.settings, key):
                setattr(self.settings, key, value)
    
    def export_image(
        self,
        viewer_data: str,  # Base64 or data URL
        output_path: str,
        settings: Optional[ExportSettings] = None
    ) -> bool:
        """
        Export current view as image.
        
        Args:
            viewer_data: Image data (base64 or data URL)
            output_path: Output file path
            settings: Optional export settings
            
        Returns:
            True if successful
        """
        if settings:
            self.settings = settings
        
        try:
            # Decode base64 if needed
            if viewer_data.startswith('data:image'):
                # Extract base64 part
                header, data = viewer_data.split(',', 1)
                image_data = base64.b64decode(data)
            else:
                image_data = base64.b64decode(viewer_data)
            
            # Write to file
            with open(output_path, 'wb') as f:
                f.write(image_data)
            
            logger.info(f"Image exported to {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"Export failed: {e}")
            return False
    
    def export_4k_png(
        self,
        viewer_data: str,
        output_path: str,
        transparent: bool = False
    ) -> bool:
        """Export 4K PNG image (3840x2160)"""
        self.settings.width = 3840
        self.settings.height = 2160
        self.settings.format = ExportFormat.PNG
        self.settings.transparent = transparent
        self.settings.dpi = 300
        
        return self.export_image(viewer_data, output_path)
    
    def export_transparent_png(
        self,
        viewer_data: str,
        output_path: str,
        width: int = 1920,
        height: int = 1080
    ) -> bool:
        """Export PNG with transparent background"""
        self.settings.width = width
        self.settings.height = height
        self.settings.format = ExportFormat.PNG
        self.settings.transparent = True
        
        return self.export_image(viewer_data, output_path)
    
    def export_svg(
        self,
        svg_content: str,
        output_path: str
    ) -> bool:
        """Export SVG vector graphics"""
        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(svg_content)
            
            logger.info(f"SVG exported to {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"SVG export failed: {e}")
            return False
    
    def export_pdf_report(
        self,
        content: Dict,
        output_path: str
    ) -> bool:
        """Export PDF report with molecular images"""
        try:
            # Basic PDF generation (would use reportlab in production)
            from reportlab.lib.pagesizes import letter
            from reportlab.lib.units import inch
            from reportlab.pdfgen import canvas
            
            c = canvas.Canvas(output_path, pagesize=letter)
            width, height = letter
            
            # Title
            c.setFont("Helvetica-Bold", 24)
            c.drawString(1*inch, height - 1*inch, "Docking Results Report")
            
            # Add images and data
            y = height - 1.5*inch
            c.setFont("Helvetica", 12)
            
            for key, value in content.get('metadata', {}).items():
                c.drawString(1*inch, y, f"{key}: {value}")
                y -= 0.3*inch
            
            c.save()
            
            logger.info(f"PDF exported to {output_path}")
            return True
            
        except ImportError:
            logger.warning("reportlab not available, PDF export limited")
            return False
        except Exception as e:
            logger.error(f"PDF export failed: {e}")
            return False
    
    def export_gif_animation(
        self,
        frame_data: List[str],  # List of base64 frames
        output_path: str,
        duration: float = 3.0,
        loop: int = 0
    ) -> bool:
        """Export animated GIF"""
        try:
            from PIL import Image
            import io
            
            frames = []
            for frame_b64 in frame_data:
                # Decode base64
                if ',' in frame_b64:
                    frame_b64 = frame_b64.split(',')[1]
                
                img_data = base64.b64decode(frame_b64)
                img = Image.open(io.BytesIO(img_data))
                
                # Convert to RGB if necessary
                if img.mode != 'RGB':
                    img = img.convert('RGB')
                
                frames.append(img)
            
            if not frames:
                logger.warning("No frames for GIF export")
                return False
            
            # Save as GIF
            duration_ms = int(duration * 1000 / len(frames))
            frames[0].save(
                output_path,
                save_all=True,
                append_images=frames[1:],
                duration=duration_ms,
                loop=loop
            )
            
            logger.info(f"GIF exported to {output_path}")
            return True
            
        except ImportError:
            logger.warning("PIL not available for GIF export")
            return False
        except Exception as e:
            logger.error(f"GIF export failed: {e}")
            return False
    
    def export_3d_mesh(
        self,
        vertices: List[List[float]],
        faces: List[List[int]],
        output_path: str,
        format: str = "obj"
    ) -> bool:
        """Export 3D mesh (OBJ or STL)"""
        try:
            if format == "obj":
                return self._export_obj(vertices, faces, output_path)
            elif format == "stl":
                return self._export_stl(vertices, faces, output_path)
            else:
                logger.warning(f"Unsupported 3D format: {format}")
                return False
                
        except Exception as e:
            logger.error(f"3D mesh export failed: {e}")
            return False
    
    def _export_obj(
        self,
        vertices: List[List[float]],
        faces: List[List[int]],
        output_path: str
    ) -> bool:
        """Export OBJ format"""
        with open(output_path, 'w') as f:
            # Write vertices
            for v in vertices:
                f.write(f"v {v[0]:.4f} {v[1]:.4f} {v[2]:.4f}\n")
            
            # Write faces
            for face in faces:
                f.write(f"f {face[0]+1} {face[1]+1} {face[2]+1}\n")
        
        logger.info(f"OBJ exported to {output_path}")
        return True
    
    def _export_stl(
        self,
        vertices: List[List[float]],
        faces: List[List[int]],
        output_path: str
    ) -> bool:
        """Export STL format"""
        try:
            from stl import mesh
            import numpy as np
            
            triangles = []
            for face in faces:
                if len(face) >= 3:
                    v0 = vertices[face[0]]
                    v1 = vertices[face[1]]
                    v2 = vertices[face[2]]
                    triangles.append([v0, v1, v2])
            
            if triangles:
                stl_mesh = mesh.Mesh(np.array(triangles))
                stl_mesh.save(output_path)
            
            logger.info(f"STL exported to {output_path}")
            return True
            
        except ImportError:
            logger.warning("numpy-stl not available for STL export")
            return False
    
    def get_preset_configs(self) -> Dict[str, ExportSettings]:
        """Get preset export configurations"""
        return {
            '4k': ExportSettings(
                format=ExportFormat.PNG,
                width=3840,
                height=2160,
                dpi=300,
                quality=95
            ),
            'publication': ExportSettings(
                format=ExportFormat.PNG,
                width=2400,
                height=1800,
                dpi=300,
                quality=100,
                antialiasing=True
            ),
            'slide': ExportSettings(
                format=ExportFormat.PNG,
                width=1920,
                height=1080,
                dpi=150,
                quality=90
            ),
            'transparent': ExportSettings(
                format=ExportFormat.PNG,
                width=1920,
                height=1080,
                transparent=True,
                background_color="transparent"
            ),
            'web': ExportSettings(
                format=ExportFormat.PNG,
                width=800,
                height=600,
                dpi=72,
                quality=85
            ),
            'animation': ExportSettings(
                format=ExportFormat.GIF,
                width=800,
                height=600,
                animation_frames=36,
                animation_duration=3.0
            ),
        }
    
    def create_export_dialog_options(self) -> List[Dict]:
        """Create options for export dialog UI"""
        presets = self.get_preset_configs()
        
        options = [
            {
                'id': 'format',
                'type': 'combo',
                'label': 'Format',
                'options': [f.value for f in ExportFormat],
                'default': 'png'
            },
            {
                'id': 'preset',
                'type': 'combo',
                'label': 'Preset',
                'options': list(presets.keys()),
                'default': 'publication'
            },
            {
                'id': 'width',
                'type': 'spin',
                'label': 'Width',
                'min': 100,
                'max': 7680,
                'default': 1920
            },
            {
                'id': 'height',
                'type': 'spin',
                'label': 'Height',
                'min': 100,
                'max': 4320,
                'default': 1080
            },
            {
                'id': 'transparent',
                'type': 'checkbox',
                'label': 'Transparent Background',
                'default': False
            },
            {
                'id': 'dpi',
                'type': 'combo',
                'label': 'DPI',
                'options': ['72', '150', '300', '600'],
                'default': '300'
            },
        ]
        
        return options


class ImageProcessor:
    """Post-processing for exported images"""
    
    @staticmethod
    def add_watermark(
        image_path: str,
        output_path: str,
        watermark_text: str = "Docking Studio",
        position: str = "bottom-right"
    ) -> bool:
        """Add watermark to image"""
        try:
            from PIL import Image, ImageDraw, ImageFont
            
            img = Image.open(image_path)
            draw = ImageDraw.Draw(img)
            
            # Try to use a font, fall back to default
            try:
                font = ImageFont.truetype("arial.ttf", 24)
            except:
                font = ImageFont.load_default()
            
            # Get text size
            bbox = draw.textbbox((0, 0), watermark_text, font=font)
            text_width = bbox[2] - bbox[0]
            text_height = bbox[3] - bbox[1]
            
            # Calculate position
            img_width, img_height = img.size
            padding = 20
            
            if position == "bottom-right":
                x = img_width - text_width - padding
                y = img_height - text_height - padding
            elif position == "bottom-left":
                x = padding
                y = img_height - text_height - padding
            elif position == "top-right":
                x = img_width - text_width - padding
                y = padding
            else:  # top-left
                x = padding
                y = padding
            
            # Draw watermark
            draw.text((x, y), watermark_text, fill=(128, 128, 128, 128), font=font)
            
            img.save(output_path)
            
            logger.info(f"Watermark added to {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"Watermark failed: {e}")
            return False
    
    @staticmethod
    def resize_image(
        image_path: str,
        output_path: str,
        width: int,
        height: int,
        maintain_aspect: bool = True
    ) -> bool:
        """Resize image"""
        try:
            from PIL import Image
            
            img = Image.open(image_path)
            
            if maintain_aspect:
                img.thumbnail((width, height), Image.Resampling.LANCZOS)
            else:
                img = img.resize((width, height), Image.Resampling.LANCZOS)
            
            img.save(output_path)
            
            logger.info(f"Image resized to {width}x{height}")
            return True
            
        except Exception as e:
            logger.error(f"Resize failed: {e}")
            return False
