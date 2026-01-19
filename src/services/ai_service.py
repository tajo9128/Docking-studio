"""
BioDockify AI Analysis Service
Automated molecular analysis insights and natural language processing
"""

from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
import logging
import json

logger = logging.getLogger(__name__)


@dataclass
class AnalysisInsight:
    """A single analysis insight."""
    category: str
    title: str
    description: str
    confidence: float
    severity: str  # 'info', 'warning', 'critical'
    recommendations: List[str] = field(default_factory=list)


@dataclass
class AIAnalysisResult:
    """Complete AI analysis result."""
    summary: str
    insights: List[AnalysisInsight] = field(default_factory=list)
    binding_analysis: Dict[str, Any] = field(default_factory=dict)
    druglikeness_summary: str = ""
    recommendations: List[str] = field(default_factory=list)
    confidence_score: float = 0.0


class AIAnalysisService:
    """
    AI-powered molecular analysis service.
    
    Provides automated insights, natural language summaries,
    and recommendations based on docking results and molecular properties.
    """
    
    def __init__(self):
        """Initialize the AI analysis service."""
        self.templates = self._load_templates()
    
    def _load_templates(self) -> Dict[str, str]:
        """Load analysis templates."""
        return {
            'binding_excellent': (
                "The compound shows excellent binding affinity with a score of {score:.2f} kcal/mol. "
                "This is in the top tier of drug candidates and warrants further investigation."
            ),
            'binding_good': (
                "The compound demonstrates good binding affinity ({score:.2f} kcal/mol), "
                "suggesting potential as a lead compound."
            ),
            'binding_moderate': (
                "The binding affinity is moderate ({score:.2f} kcal/mol). "
                "Structural optimization may improve binding."
            ),
            'binding_weak': (
                "The binding is relatively weak ({score:.2f} kcal/mol). "
                "Consider alternative scaffolds or significant modifications."
            ),
            'interaction_summary': (
                "The complex shows {hbonds} hydrogen bonds, {hydrophobic} hydrophobic contacts, "
                "and {pistacking} π-stacking interactions."
            ),
            'druglikeness_pass': (
                "The compound passes drug-likeness filters with a score of {score:.0%}. "
                "It is suitable for further development."
            ),
            'druglikeness_fail': (
                "The compound has drug-likeness concerns ({issues}). "
                "Consider structural modifications to improve properties."
            ),
            'admet_warning': (
                "ADMET analysis indicates potential issues: {warnings}. "
                "These should be addressed in lead optimization."
            ),
        }
    
    def analyze_docking_result(self, docking_result: Dict) -> AIAnalysisResult:
        """
        Analyze a docking result and generate insights.
        
        Args:
            docking_result: Dictionary containing docking scores and poses
            
        Returns:
            AIAnalysisResult with comprehensive analysis
        """
        insights = []
        recommendations = []
        
        # Analyze binding score
        score = docking_result.get('score', 0)
        binding_insight = self._analyze_binding_score(score)
        insights.append(binding_insight)
        
        # Analyze interactions
        interactions = docking_result.get('interactions', {})
        interaction_insight = self._analyze_interactions(interactions)
        if interaction_insight:
            insights.append(interaction_insight)
        
        # Analyze pose quality
        poses = docking_result.get('poses', [])
        pose_insight = self._analyze_poses(poses)
        if pose_insight:
            insights.append(pose_insight)
        
        # Generate summary
        summary = self._generate_summary(score, interactions, insights)
        
        # Generate recommendations
        recommendations = self._generate_recommendations(insights)
        
        # Calculate confidence
        confidence = self._calculate_confidence(docking_result)
        
        return AIAnalysisResult(
            summary=summary,
            insights=insights,
            binding_analysis={
                'score': score,
                'category': self._categorize_score(score),
                'interactions': interactions
            },
            recommendations=recommendations,
            confidence_score=confidence
        )
    
    def _analyze_binding_score(self, score: float) -> AnalysisInsight:
        """Analyze the binding score."""
        if score <= -10:
            return AnalysisInsight(
                category='binding',
                title='Excellent Binding Affinity',
                description=self.templates['binding_excellent'].format(score=score),
                confidence=0.9,
                severity='info',
                recommendations=['Proceed to experimental validation']
            )
        elif score <= -8:
            return AnalysisInsight(
                category='binding',
                title='Good Binding Affinity',
                description=self.templates['binding_good'].format(score=score),
                confidence=0.8,
                severity='info',
                recommendations=['Consider lead optimization']
            )
        elif score <= -6:
            return AnalysisInsight(
                category='binding',
                title='Moderate Binding Affinity',
                description=self.templates['binding_moderate'].format(score=score),
                confidence=0.7,
                severity='warning',
                recommendations=['Explore structural modifications', 'Check for better poses']
            )
        else:
            return AnalysisInsight(
                category='binding',
                title='Weak Binding Affinity',
                description=self.templates['binding_weak'].format(score=score),
                confidence=0.6,
                severity='critical',
                recommendations=['Consider alternative compounds', 'Review docking parameters']
            )
    
    def _analyze_interactions(self, interactions: Dict) -> Optional[AnalysisInsight]:
        """Analyze molecular interactions."""
        hbonds = interactions.get('hydrogen_bonds', 0)
        hydrophobic = interactions.get('hydrophobic', 0)
        pistacking = interactions.get('pi_stacking', 0)
        salt_bridges = interactions.get('salt_bridges', 0)
        
        total = hbonds + hydrophobic + pistacking + salt_bridges
        
        if total == 0:
            return AnalysisInsight(
                category='interactions',
                title='No Significant Interactions',
                description='No significant protein-ligand interactions detected.',
                confidence=0.5,
                severity='warning',
                recommendations=['Review binding pose', 'Check ligand positioning']
            )
        
        desc = self.templates['interaction_summary'].format(
            hbonds=hbonds, hydrophobic=hydrophobic, pistacking=pistacking
        )
        
        # Evaluate quality
        if hbonds >= 3 and total >= 6:
            severity = 'info'
            title = 'Strong Interaction Network'
        elif hbonds >= 1 and total >= 3:
            severity = 'info'
            title = 'Good Interaction Profile'
        else:
            severity = 'warning'
            title = 'Limited Interactions'
        
        return AnalysisInsight(
            category='interactions',
            title=title,
            description=desc,
            confidence=0.75,
            severity=severity,
            recommendations=[]
        )
    
    def _analyze_poses(self, poses: List) -> Optional[AnalysisInsight]:
        """Analyze pose diversity and quality."""
        if not poses:
            return None
        
        num_poses = len(poses)
        
        if num_poses >= 9:
            return AnalysisInsight(
                category='poses',
                title='Diverse Binding Modes',
                description=f'{num_poses} distinct poses found, indicating binding flexibility.',
                confidence=0.7,
                severity='info',
                recommendations=['Compare top poses for consistency']
            )
        elif num_poses >= 3:
            return AnalysisInsight(
                category='poses',
                title='Multiple Binding Modes',
                description=f'{num_poses} poses identified.',
                confidence=0.8,
                severity='info',
                recommendations=[]
            )
        else:
            return AnalysisInsight(
                category='poses',
                title='Limited Pose Diversity',
                description=f'Only {num_poses} pose(s) found.',
                confidence=0.6,
                severity='warning',
                recommendations=['Consider increasing exhaustiveness']
            )
    
    def _categorize_score(self, score: float) -> str:
        """Categorize binding score."""
        if score <= -10:
            return 'excellent'
        elif score <= -8:
            return 'good'
        elif score <= -6:
            return 'moderate'
        else:
            return 'weak'
    
    def _generate_summary(self, score: float, interactions: Dict, 
                          insights: List[AnalysisInsight]) -> str:
        """Generate a natural language summary."""
        category = self._categorize_score(score)
        
        summary_parts = []
        
        # Binding summary
        if category == 'excellent':
            summary_parts.append(f"Excellent binding affinity ({score:.2f} kcal/mol).")
        elif category == 'good':
            summary_parts.append(f"Good binding ({score:.2f} kcal/mol).")
        elif category == 'moderate':
            summary_parts.append(f"Moderate binding ({score:.2f} kcal/mol).")
        else:
            summary_parts.append(f"Weak binding ({score:.2f} kcal/mol).")
        
        # Interaction summary
        hbonds = interactions.get('hydrogen_bonds', 0)
        if hbonds > 0:
            summary_parts.append(f"{hbonds} hydrogen bond(s) detected.")
        
        # Warnings
        warnings = [i for i in insights if i.severity in ['warning', 'critical']]
        if warnings:
            summary_parts.append(f"{len(warnings)} concern(s) noted.")
        
        return " ".join(summary_parts)
    
    def _generate_recommendations(self, insights: List[AnalysisInsight]) -> List[str]:
        """Generate overall recommendations."""
        recommendations = []
        
        for insight in insights:
            recommendations.extend(insight.recommendations)
        
        # Remove duplicates while preserving order
        seen = set()
        unique = []
        for r in recommendations:
            if r not in seen:
                seen.add(r)
                unique.append(r)
        
        return unique[:5]  # Top 5 recommendations
    
    def _calculate_confidence(self, result: Dict) -> float:
        """Calculate overall analysis confidence."""
        score = result.get('score', 0)
        interactions = result.get('interactions', {})
        
        confidence = 0.5
        
        # Better score = higher confidence
        if score <= -8:
            confidence += 0.2
        elif score <= -6:
            confidence += 0.1
        
        # More interactions = higher confidence
        total_interactions = sum(interactions.values()) if interactions else 0
        if total_interactions >= 5:
            confidence += 0.2
        elif total_interactions >= 2:
            confidence += 0.1
        
        return min(confidence, 1.0)
    
    def analyze_druglikeness(self, druglikeness_result: Dict) -> str:
        """Generate natural language summary of drug-likeness."""
        overall = druglikeness_result.get('overall_druglikeness', 'unknown')
        score = druglikeness_result.get('score', 0)
        
        lipinski = druglikeness_result.get('lipinski', {})
        violations = lipinski.get('violations', 0)
        
        if overall == 'excellent':
            return f"Excellent drug-like properties (score: {score:.0%}). Fully compliant with Lipinski's rules."
        elif overall == 'good':
            return f"Good drug-like properties (score: {score:.0%}). Minor optimization may be beneficial."
        elif overall == 'moderate':
            return f"Moderate drug-likeness (score: {score:.0%}). {violations} Lipinski violation(s) noted."
        else:
            return f"Poor drug-likeness (score: {score:.0%}). Significant property optimization required."
    
    def analyze_admet(self, admet_result: Dict) -> str:
        """Generate natural language summary of ADMET predictions."""
        warnings = admet_result.get('warnings', [])
        overall = admet_result.get('overall_score', 0)
        
        if not warnings:
            return f"ADMET profile is favorable (score: {overall:.0%}). No major concerns identified."
        elif len(warnings) == 1:
            return f"ADMET profile is acceptable (score: {overall:.0%}). Note: {warnings[0]}"
        else:
            issues = "; ".join(warnings[:3])
            return f"ADMET concerns identified (score: {overall:.0%}): {issues}"
    
    def generate_report(self, docking_result: Dict, druglikeness: Optional[Dict] = None,
                        admet: Optional[Dict] = None) -> str:
        """Generate a complete analysis report."""
        report_parts = ["# Molecular Analysis Report\n"]
        
        # Docking analysis
        analysis = self.analyze_docking_result(docking_result)
        report_parts.append(f"## Binding Analysis\n{analysis.summary}\n")
        
        # Insights
        if analysis.insights:
            report_parts.append("## Key Insights\n")
            for insight in analysis.insights:
                icon = "✅" if insight.severity == 'info' else "⚠️" if insight.severity == 'warning' else "❌"
                report_parts.append(f"- {icon} **{insight.title}**: {insight.description}\n")
        
        # Drug-likeness
        if druglikeness:
            summary = self.analyze_druglikeness(druglikeness)
            report_parts.append(f"\n## Drug-likeness\n{summary}\n")
        
        # ADMET
        if admet:
            summary = self.analyze_admet(admet)
            report_parts.append(f"\n## ADMET Profile\n{summary}\n")
        
        # Recommendations
        if analysis.recommendations:
            report_parts.append("\n## Recommendations\n")
            for i, rec in enumerate(analysis.recommendations, 1):
                report_parts.append(f"{i}. {rec}\n")
        
        return "".join(report_parts)
