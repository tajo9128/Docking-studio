"""
Scene Animator
Pose animation and rotation for molecular visualizations.
"""

from enum import Enum
from dataclasses import dataclass
from typing import List, Dict, Optional, Callable, Tuple
import numpy as np
import logging

logger = logging.getLogger(__name__)


class AnimationMode(Enum):
    """Animation modes"""
    ROTATION = "rotation"
    POSE_SEQUENCE = "pose_sequence"
    ZOOM = "zoom"
    PAN = "pan"
    MORPH = "morph"
    ROCK = "rock"


class RotationAxis(Enum):
    """Rotation axes"""
    X = "x"
    Y = "y"
    Z = "z"
    AUTO = "auto"


@dataclass
class AnimationFrame:
    """A single animation frame"""
    index: int
    rotation: Tuple[float, float, float]  # x, y, z degrees
    zoom: float
    center: Tuple[float, float, float]
    duration: float = 0.0  # Duration in ms


@dataclass
class AnimationSequence:
    """Complete animation sequence"""
    name: str
    frames: List[AnimationFrame]
    duration: float  # Total duration in seconds
    loop: bool = True
    fps: int = 30


class SceneAnimator:
    """
    Manages animations for molecular scenes.
    
    Features:
    - Auto-rotation
    - Pose sequence animation
    - Smooth transitions
    - Export to GIF
    """
    
    def __init__(self):
        """Initialize scene animator"""
        self.is_animating = False
        self.current_frame = 0
        self.total_frames = 0
        
        self.animation_speed = 1.0
        self.rotation_axis = RotationAxis.Y
        self.rotation_speed = 0.5  # Degrees per frame
        
        self.frames: List[str] = []  # Base64 frame data
        self.animation_complete_callback: Optional[Callable] = None
        
        logger.info("SceneAnimator initialized")
    
    def start_rotation(
        self,
        axis: RotationAxis = RotationAxis.Y,
        speed: float = 0.5,
        callback: Optional[Callable] = None
    ):
        """Start continuous rotation animation"""
        self.rotation_axis = axis
        self.rotation_speed = speed
        self.animation_complete_callback = callback
        
        self.is_animating = True
        self.current_frame = 0
        
        logger.info(f"Started rotation animation: axis={axis.value}, speed={speed}")
    
    def stop_animation(self):
        """Stop current animation"""
        self.is_animating = False
        logger.info("Animation stopped")
    
    def pause_animation(self):
        """Pause animation"""
        self.is_animating = False
        logger.info("Animation paused")
    
    def resume_animation(self):
        """Resume animation"""
        if self.current_frame > 0:
            self.is_animating = True
            logger.info("Animation resumed")
    
    def get_next_frame_params(self) -> Dict:
        """Get parameters for next animation frame"""
        if not self.is_animating:
            return {
                'rotation': (0, 0, 0),
                'zoom': 1.0,
                'center': (0, 0, 0),
                'frame': self.current_frame
            }
        
        # Calculate rotation based on axis
        rotation = self._calculate_rotation()
        
        # Update frame counter
        self.current_frame = (self.current_frame + 1) % 360
        
        return {
            'rotation': rotation,
            'zoom': 1.0,
            'center': (0, 0, 0),
            'frame': self.current_frame
        }
    
    def _calculate_rotation(self) -> Tuple[float, float, float]:
        """Calculate rotation angles for current frame"""
        angle = self.current_frame * self.rotation_speed
        
        if self.rotation_axis == RotationAxis.X:
            return (angle, 0, 0)
        elif self.rotation_axis == RotationAxis.Y:
            return (0, angle, 0)
        elif self.rotation_axis == RotationAxis.Z:
            return (0, 0, angle)
        else:  # AUTO - rotate around all axes slightly
            return (angle * 0.3, angle, angle * 0.5)
    
    def create_pose_animation(
        self,
        poses: List[str],  # List of PDB content
        durations: Optional[List[float]] = None,
        transition_type: str = "fade"
    ) -> AnimationSequence:
        """
        Create animation from pose sequence.
        
        Args:
            poses: List of PDB content strings
            durations: Optional list of durations per pose
            transition_type: How to transition between poses
            
        Returns:
            AnimationSequence
        """
        if not poses:
            logger.warning("No poses provided for animation")
            return AnimationSequence(name="empty", frames=[], duration=0)
        
        n_frames = len(poses)
        
        # Default durations (1 second per pose)
        if durations is None:
            durations = [1.0] * n_frames
        
        total_duration = sum(durations)
        fps = 30
        total_frames = int(total_duration * fps)
        
        frames = []
        
        for i in range(total_frames):
            # Determine which pose we're on
            time = i / fps
            cumsum = 0
            pose_idx = 0
            
            for j, dur in enumerate(durations):
                cumsum += dur
                if time < cumsum:
                    pose_idx = j
                    break
            
            # Create frame
            frame = AnimationFrame(
                index=i,
                rotation=(0, 0, 0),  # Static rotation for pose sequence
                zoom=1.0,
                center=(0, 0, 0),
                duration=1000 / fps  # ms
            )
            frames.append(frame)
        
        sequence = AnimationSequence(
            name="pose_sequence",
            frames=frames,
            duration=total_duration,
            loop=True,
            fps=fps
        )
        
        logger.info(f"Created pose animation: {n_frames} poses, {total_duration}s")
        
        return sequence
    
    def create_zoom_animation(
        self,
        start_zoom: float = 1.0,
        end_zoom: float = 2.0,
        duration: float = 2.0,
        pause: float = 1.0
    ) -> AnimationSequence:
        """Create zoom in/out animation"""
        fps = 30
        total_frames = int(duration * fps)
        
        frames = []
        
        for i in range(total_frames):
            t = i / total_frames
            
            # Ease in-out
            if t < 0.5:
                progress = 2 * t * t
            else:
                progress = 1 - pow(-2 * t + 2, 2) / 2
            
            zoom = start_zoom + (end_zoom - start_zoom) * progress
            
            frame = AnimationFrame(
                index=i,
                rotation=(0, 0, 0),
                zoom=zoom,
                center=(0, 0, 0)
            )
            frames.append(frame)
        
        # Add pause frames
        pause_frames = int(pause * fps)
        for i in range(pause_frames):
            frame = AnimationFrame(
                index=total_frames + i,
                rotation=(0, 0, 0),
                zoom=end_zoom,
                center=(0, 0, 0)
            )
            frames.append(frame)
        
        total_duration = duration + pause * 2
        
        return AnimationSequence(
            name="zoom",
            frames=frames,
            duration=total_duration,
            loop=True,
            fps=fps
        )
    
    def create_rock_animation(
        self,
        angle: float = 30,
        duration: float = 2.0
    ) -> AnimationSequence:
        """Create rock animation (back and forth)"""
        fps = 30
        total_frames = int(duration * fps)
        
        frames = []
        
        for i in range(total_frames):
            # Sine wave for back and forth
            t = i / total_frames
            angle_rad = angle * np.sin(2 * np.pi * t)
            
            frame = AnimationFrame(
                index=i,
                rotation=(angle_rad, 0, 0),
                zoom=1.0,
                center=(0, 0, 0)
            )
            frames.append(frame)
        
        return AnimationSequence(
            name="rock",
            frames=frames,
            duration=duration,
            loop=True,
            fps=fps
        )
    
    def generate_frames_for_export(
        self,
        viewer,
        sequence: AnimationSequence,
        max_frames: int = 36
    ) -> List[str]:
        """
        Generate frames for animation export.
        
        Args:
            viewer: The molecular viewer instance
            sequence: AnimationSequence to render
            max_frames: Maximum frames to capture
            
        Returns:
            List of base64 frame data
        """
        frames = []
        frame_indices = np.linspace(
            0, 
            len(sequence.frames) - 1, 
            min(max_frames, len(sequence.frames))
        ).astype(int)
        
        logger.info(f"Generating {len(frame_indices)} frames for export")
        
        for idx in frame_indices:
            # Apply frame parameters
            frame = sequence.frames[idx]
            
            # Rotate viewer
            if viewer:
                viewer.page().runJavaScript(
                    f'viewer.rotate({frame.rotation[0]}, "x");'
                    f'viewer.rotate({frame.rotation[1]}, "y");'
                    f'viewer.rotate({frame.rotation[2]}, "z");'
                    f'viewer.zoom({frame.zoom});'
                )
            
            # Capture frame
            # Note: Would need to capture from viewer
            # This is a placeholder
            
        self.frames = frames
        return frames
    
    def export_gif(
        self,
        frames: List[str],
        output_path: str,
        duration: float = 3.0
    ) -> bool:
        """Export frames as GIF"""
        try:
            from ..export_manager import ExportManager
            export_mgr = ExportManager()
            return export_mgr.export_gif_animation(frames, output_path, duration)
        except Exception as e:
            logger.error(f"GIF export failed: {e}")
            return False
    
    def get_javascript_config(self, sequence: AnimationSequence) -> str:
        """Get JavaScript configuration for web-based animation"""
        frames_data = []
        
        for frame in sequence.frames:
            frames_data.append({
                'r': frame.rotation,
                'z': frame.zoom,
                'c': frame.center
            })
        
        import json
        return json.dumps({
            'frames': frames_data,
            'fps': sequence.fps,
            'loop': sequence.loop,
            'duration': sequence.duration * 1000
        })
    
    def set_speed(self, speed: float):
        """Set animation speed multiplier"""
        self.animation_speed = max(0.1, min(10.0, speed))
        logger.info(f"Animation speed set to {self.animation_speed}x")
    
    def get_preset_animations(self) -> Dict[str, Dict]:
        """Get preset animation configurations"""
        return {
            'spin_y': {
                'mode': AnimationMode.ROTATION,
                'axis': RotationAxis.Y,
                'speed': 0.5,
            },
            'spin_xyz': {
                'mode': AnimationMode.ROTATION,
                'axis': RotationAxis.AUTO,
                'speed': 0.3,
            },
            'rock': {
                'mode': AnimationMode.ROCK,
                'angle': 30,
                'duration': 2.0,
            },
            'pulse': {
                'mode': AnimationMode.ZOOM,
                'start_zoom': 1.0,
                'end_zoom': 1.5,
                'duration': 1.0,
            },
        }


class TrajectoryPlayer:
    """Play back molecular dynamics trajectories"""
    
    def __init__(self):
        """Initialize trajectory player"""
        self.frames: List[str] = []
        self.current_frame = 0
        self.is_playing = False
        self.fps = 30
        
        logger.info("TrajectoryPlayer initialized")
    
    def load_trajectory(
        self,
        trajectory_file: str,
        format: str = "dcd"
    ) -> bool:
        """Load molecular dynamics trajectory"""
        logger.info(f"Loading trajectory from {trajectory_file}")
        
        # Placeholder - would need MDAnalysis or similar
        # For now, just log
        self.frames = []
        
        return False
    
    def play(self):
        """Start playback"""
        self.is_playing = True
        logger.info("Trajectory playback started")
    
    def pause(self):
        """Pause playback"""
        self.is_playing = False
        logger.info("Trajectory playback paused")
    
    def stop(self):
        """Stop and reset"""
        self.is_playing = False
        self.current_frame = 0
        logger.info("Trajectory stopped")
    
    def seek(self, frame: int):
        """Seek to specific frame"""
        if 0 <= frame < len(self.frames):
            self.current_frame = frame
    
    def get_frame(self) -> Optional[str]:
        """Get current frame data"""
        if 0 <= self.current_frame < len(self.frames):
            return self.frames[self.current_frame]
        return None
    
    def next_frame(self):
        """Advance to next frame"""
        if self.current_frame < len(self.frames) - 1:
            self.current_frame += 1
    
    def prev_frame(self):
        """Go to previous frame"""
        if self.current_frame > 0:
            self.current_frame -= 1
    
    def get_progress(self) -> float:
        """Get playback progress (0-1)"""
        if not self.frames:
            return 0
        return self.current_frame / len(self.frames)
