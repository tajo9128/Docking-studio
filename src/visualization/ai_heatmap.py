"""
AI Heatmap Generator
Visualizes AI confidence scores and predictions on molecular structures.
"""

from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
import numpy as np
import logging

logger = logging.getLogger(__name__)


@dataclass
class AIInsight:
    """AI prediction/insight for a region"""
    residue: str
    residue_index: int
    confidence: float  # 0-1
    contribution: float  # Contribution to binding
    prediction: str  # e.g., "favorable", "unfavorable", "neutral"
    details: Dict


@dataclass
class HeatmapData:
    """Heatmap visualization data"""
    residues: List[str]
    indices: List[int]
    values: List[float]  # Color values
    colors: List[str]  # Hex colors
    annotations: List[str]  # Text annotations


class AIHeatmapGenerator:
    """
    Generates AI-based heatmaps for molecular visualization.
    
    Provides:
    - Per-residue confidence scores
    - Binding contribution heatmaps
    - Key residue highlighting
    - Integration with Agent Zero predictions
    """
    
    # Heatmap color schemes
    CONFIDENCE_COLORS = {
        'high': '#00FF00',      # Green - high confidence
        'medium': '#FFFF00',    # Yellow - medium
        'low': '#FF6600',       # Orange - low
        'neutral': '#FFFFFF',   # White - neutral
    }
    
    CONTRIBUTION_COLORS = {
        'favorable': '#0088FF',   # Blue
        'unfavorable': '#FF4444', # Red
        'neutral': '#888888',     # Gray
    }
    
    def __init__(self):
        """Initialize AI heatmap generator"""
        self.current_heatmap: Optional[HeatmapData] = None
        self.insights: List[AIInsight] = []
        
        logger.info("AIHeatmapGenerator initialized")
    
    def generate_from_predictions(
        self,
        predictions: List[Dict],
        residues: List[str]
    ) -> HeatmapData:
        """
        Generate heatmap from AI predictions.
        
        Args:
            predictions: List of prediction dictionaries
            residues: List of residue identifiers
            
        Returns:
            HeatmapData for visualization
        """
        if not predictions or not residues:
            return self._empty_heatmap()
        
        values = []
        colors = []
        annotations = []
        
        # Build prediction lookup
        pred_lookup = {p.get('residue', ''): p for p in predictions}
        
        for res in residues:
            pred = pred_lookup.get(res, {})
            
            # Get confidence
            confidence = pred.get('confidence', 0.5)
            values.append(confidence)
            
            # Generate color
            color = self._confidence_to_color(confidence)
            colors.append(color)
            
            # Annotation
            contribution = pred.get('contribution', 0.0)
            annotations.append(f"{res}\n{contribution:.2f}")
        
        heatmap = HeatmapData(
            residues=residues,
            indices=list(range(len(residues))),
            values=values,
            colors=colors,
            annotations=annotations
        )
        
        self.current_heatmap = heatmap
        
        logger.info(f"Generated heatmap for {len(residues)} residues")
        
        return heatmap
    
    def generate_from_interaction_scores(
        self,
        interactions: List[Dict],
        residues: List[str]
    ) -> HeatmapData:
        """
        Generate heatmap from interaction analysis.
        
        Args:
            interactions: List of interaction dictionaries
            residues: List of residue identifiers
            
        Returns:
            HeatmapData for visualization
        """
        # Count interactions per residue
        interaction_counts: Dict[str, int] = {res: 0 for res in residues}
        
        for interaction in interactions:
            res = interaction.get('residue1', '')
            if res in interaction_counts:
                interaction_counts[res] += 1
        
        # Normalize to 0-1
        max_count = max(interaction_counts.values()) if interaction_counts else 1
        
        values = []
        colors = []
        annotations = []
        
        for res in residues:
            count = interaction_counts.get(res, 0)
            value = count / max_count if max_count > 0 else 0
            values.append(value)
            
            color = self._contribution_to_color(value)
            colors.append(color)
            
            annotations.append(f"{res}\n{count} interactions")
        
        heatmap = HeatmapData(
            residues=residues,
            indices=list(range(len(residues))),
            values=values,
            colors=colors,
            annotations=annotations
        )
        
        self.current_heatmap = heatmap
        
        return heatmap
    
    def generate_from_ai_scores(
        self,
        ai_scores: Dict[str, float],
        residues: List[str]
    ) -> HeatmapData:
        """
        Generate heatmap from AI scoring results.
        
        Args:
            ai_scores: Dictionary of residue -> score
            residues: List of residue identifiers
            
        Returns:
            HeatmapData for visualization
        """
        # Normalize scores
        all_scores = list(ai_scores.values())
        if all_scores:
            min_score = min(all_scores)
            max_score = max(all_scores)
            score_range = max_score - min_score if max_score != min_score else 1
        else:
            min_score = max_score = 0
            score_range = 1
        
        values = []
        colors = []
        annotations = []
        
        for res in residues:
            score = ai_scores.get(res, 0.5)
            normalized = (score - min_score) / score_range
            values.append(normalized)
            
            color = self._score_to_color(normalized)
            colors.append(color)
            
            annotations.append(f"{res}\n{score:.3f}")
        
        heatmap = HeatmapData(
            residues=residues,
            indices=list(range(len(residues))),
            values=values,
            colors=colors,
            annotations=annotations
        )
        
        self.current_heatmap = heatmap
        
        return heatmap
    
    def highlight_key_residues(
        self,
        residues: List[str],
        threshold: float = 0.7,
        top_n: int = 10
    ) -> List[str]:
        """
        Identify key residues based on current heatmap.
        
        Args:
            residues: List of residue identifiers
            threshold: Minimum value for key residue
            top_n: Number of top residues to return
            
        Returns:
            List of key residue identifiers
        """
        if self.current_heatmap is None:
            return []
        
        # Sort by value
        indexed_values = list(zip(
            self.current_heatmap.residues,
            self.current_heatmap.values
        ))
        
        # Filter and sort
        key_residues = [
            res for res, val in indexed_values
            if val >= threshold
        ]
        
        # Add top N
        sorted_by_value = sorted(indexed_values, key=lambda x: x[1], reverse=True)
        for res, val in sorted_by_value[:top_n]:
            if res not in key_residues:
                key_residues.append(res)
        
        return key_residues[:top_n]
    
    def get_residue_colors_3dmol(self) -> Dict[str, str]:
        """Get residue colors for 3Dmol.js highlighting"""
        if self.current_heatmap is None:
            return {}
        
        return {
            res: color 
            for res, color in zip(
                self.current_heatmap.residues,
                self.current_heatmap.colors
            )
        }
    
    def _confidence_to_color(self, confidence: float) -> str:
        """Convert confidence to hex color"""
        if confidence >= 0.8:
            return self.CONFIDENCE_COLORS['high']
        elif confidence >= 0.5:
            return self.CONFIDENCE_COLORS['medium']
        elif confidence >= 0.3:
            return self.CONFIDENCE_COLORS['low']
        else:
            return self.CONFIDENCE_COLORS['neutral']
    
    def _contribution_to_color(self, contribution: float) -> str:
        """Convert contribution to hex color"""
        # Contribution: -1 to 1
        if contribution > 0.2:
            # Favorable - blue gradient
            t = min(1.0, contribution)
            r, g, b = 0, int(136 * t), int(255 * t)
        elif contribution < -0.2:
            # Unfavorable - red gradient
            t = min(1.0, -contribution)
            r, g, b = int(255 * t), int(68 * t), int(68 * t)
        else:
            # Neutral
            r, g, b = 136, 136, 136
        
        return f'#{r:02x}{g:02x}{b:02x}'
    
    def _score_to_color(self, score: float) -> str:
        """Convert normalized score to hex color"""
        # Score: 0 to 1
        # Red -> Yellow -> Green
        if score < 0.5:
            t = score * 2
            r = 255
            g = int(255 * t)
            b = 0
        else:
            t = (score - 0.5) * 2
            r = int(255 * (1 - t))
            g = 255
            b = 0
        
        return f'#{r:02x}{g:02x}{b:02x}'
    
    def _empty_heatmap(self) -> HeatmapData:
        """Return empty heatmap"""
        return HeatmapData(
            residues=[],
            indices=[],
            values=[],
            colors=[],
            annotations=[]
        )
    
    def create_legend(self) -> Dict:
        """Create legend data for UI"""
        return {
            'title': 'AI Confidence',
            'type': 'gradient',
            'colors': [
                {'position': 0.0, 'color': '#FF6600', 'label': 'Low (<0.3)'},
                {'position': 0.5, 'color': '#FFFF00', 'label': 'Medium (0.3-0.5)'},
                {'position': 1.0, 'color': '#00FF00', 'label': 'High (>0.8)'},
            ],
        }
    
    def export_to_csv(self, filename: str):
        """Export heatmap data to CSV"""
        if self.current_heatmap is None:
            return
        
        import csv
        
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Residue', 'Index', 'Value', 'Color', 'Annotation'])
            
            for i in range(len(self.current_heatmap.residues)):
                writer.writerow([
                    self.current_heatmap.residues[i],
                    self.current_heatmap.indices[i],
                    self.current_heatmap.values[i],
                    self.current_heatmap.colors[i],
                    self.current_heatmap.annotations[i],
                ])
        
        logger.info(f"Heatmap exported to {filename}")


class BindingHotspotDetector:
    """Detect binding hotspots from AI predictions"""
    
    def __init__(self):
        self.threshold = 0.6
    
    def detect_hotspots(
        self,
        predictions: List[Dict],
        cluster_distance: float = 5.0
    ) -> List[Dict]:
        """Detect binding hotspots"""
        hotspots = []
        
        # Filter high-confidence predictions
        high_conf = [p for p in predictions if p.get('confidence', 0) > self.threshold]
        
        # Cluster nearby residues
        clusters = []
        for pred in high_conf:
            pos = pred.get('position', [0, 0, 0])
            
            # Check existing clusters
            added = False
            for cluster in clusters:
                center = cluster['center']
                dist = np.sqrt(
                    (pos[0] - center[0])**2 +
                    (pos[1] - center[1])**2 +
                    (pos[2] - center[2])**2
                )
                
                if dist < cluster_distance:
                    cluster['positions'].append(pos)
                    cluster['residues'].append(pred.get('residue', ''))
                    cluster['confidence'].append(pred.get('confidence', 0))
                    # Update center
                    cluster['center'] = [
                        sum(p[0] for p in cluster['positions']) / len(cluster['positions']),
                        sum(p[1] for p in cluster['positions']) / len(cluster['positions']),
                        sum(p[2] for p in cluster['positions']) / len(cluster['positions']),
                    ]
                    added = True
                    break
            
            if not added:
                clusters.append({
                    'positions': [pos],
                    'residues': [pred.get('residue', '')],
                    'confidence': [pred.get('confidence', 0)],
                    'center': pos,
                })
        
        # Create hotspot dictionaries
        for cluster in clusters:
            if len(cluster['residues']) >= 2:  # Minimum 2 residues
                hotspots.append({
                    'center': cluster['center'],
                    'residues': cluster['residues'],
                    'avg_confidence': sum(cluster['confidence']) / len(cluster['confidence']),
                    'count': len(cluster['residues']),
                })
        
        return sorted(hotspots, key=lambda x: x['avg_confidence'], reverse=True)
