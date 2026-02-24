"""
PDF Report Generator
Generates professional PDF reports for docking results.
"""

import logging
from typing import Dict, Any, List, Optional
from pathlib import Path
from datetime import datetime

logger = logging.getLogger(__name__)

try:
    from reportlab.lib.pagesizes import letter, A4
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import inch
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
    from reportlab.lib import colors
    REPORTLAB_AVAILABLE = True
except ImportError:
    REPORTLAB_AVAILABLE = False


class PDFReportGenerator:
    """
    Generates professional PDF reports for docking results.
    """
    
    def __init__(self):
        self.styles = getSampleStyleSheet() if REPORTLAB_AVAILABLE else None
    
    def generate(
        self,
        job_id: str,
        output_dir: str,
        poses: List[Dict[str, Any]],
        best_pose: Optional[Dict[str, Any]],
        **results
    ) -> Optional[Path]:
        """
        Generate PDF report.
        
        Args:
            job_id: Job identifier
            output_dir: Output directory
            poses: All poses
            best_pose: Best pose
            **results: Additional results (vina, gnina, rf, etc.)
            
        Returns:
            Path to generated PDF or None
        """
        if not REPORTLAB_AVAILABLE:
            logger.warning("reportlab not installed - cannot generate PDF")
            return None
        
        try:
            output_path = Path(output_dir) / f"results_{job_id}.pdf"
            
            doc = SimpleDocTemplate(
                str(output_path),
                pagesize=A4,
                rightMargin=72,
                leftMargin=72,
                topMargin=72,
                bottomMargin=18
            )
            
            story = []
            
            story.extend(self._create_header(job_id))
            
            story.extend(self._create_summary(poses, best_pose, **results))
            
            story.extend(self._create_top_poses_table(poses[:10]))
            
            if results.get("rf_results"):
                story.extend(self._create_rf_section(results["rf_results"]))
            
            if results.get("mmgbsa_result"):
                story.extend(self._create_mmgbsa_section(results["mmgbsa_result"]))
            
            story.extend(self._create_footer())
            
            doc.build(story)
            
            logger.info(f"PDF report generated: {output_path}")
            return output_path
            
        except Exception as e:
            logger.error(f"PDF generation failed: {e}")
            return None
    
    def _create_header(self, job_id: str) -> List:
        """Create report header"""
        story = []
        
        title_style = ParagraphStyle(
            'CustomTitle',
            parent=self.styles['Heading1'],
            fontSize=24,
            textColor=colors.HexColor('#1976D2'),
            spaceAfter=30
        )
        
        story.append(Paragraph("BioDockify Docking Results", title_style))
        story.append(Spacer(1, 12))
        
        info = [
            f"Job ID: {job_id}",
            f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"Report Type: Molecular Docking Analysis"
        ]
        
        for item in info:
            story.append(Paragraph(item, self.styles['Normal']))
            story.append(Spacer(1, 6))
        
        story.append(Spacer(1, 20))
        story.append(Paragraph("=" * 60, self.styles['Normal']))
        story.append(Spacer(1, 20))
        
        return story
    
    def _create_summary(
        self,
        poses: List[Dict],
        best_pose: Optional[Dict],
        **results
    ) -> List:
        """Create summary section"""
        story = []
        
        story.append(Paragraph("Executive Summary", self.styles['Heading2']))
        story.append(Spacer(1, 12))
        
        summary_data = [
            ["Metric", "Value"],
            ["Total Poses", str(len(poses))],
        ]
        
        if best_pose:
            summary_data.append([
                "Best Binding Energy",
                f"{best_pose.get('binding_energy', 'N/A')} kcal/mol"
            ])
            summary_data.append([
                "Best Consensus Score",
                f"{best_pose.get('consensus_score', 'N/A')}"
            ])
        
        if results.get("gnina_result"):
            summary_data.append([
                "GNINA CNN Score",
                f"{best_pose.get('gnina_cnn_score', 'N/A')}" if best_pose else "N/A"
            ])
        
        if results.get("rf_results"):
            summary_data.append([
                "RF Predicted pKd",
                f"{best_pose.get('rf_predicted_pKd', 'N/A')}" if best_pose else "N/A"
            ])
        
        table = Table(summary_data, colWidths=[2.5*inch, 3*inch])
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 12),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.black)
        ]))
        
        story.append(table)
        story.append(Spacer(1, 30))
        
        return story
    
    def _create_top_poses_table(self, poses: List[Dict]) -> List:
        """Create top poses table"""
        story = []
        
        story.append(Paragraph("Top 10 Poses", self.styles['Heading2']))
        story.append(Spacer(1, 12))
        
        table_data = [["Rank", "Mode", "Energy", "CNN Score", "Consensus"]]
        
        for idx, pose in enumerate(poses, 1):
            table_data.append([
                str(idx),
                str(pose.get('mode', '-')),
                f"{pose.get('binding_energy', '-')}",
                f"{pose.get('gnina_cnn_score', '-')}",
                f"{pose.get('consensus_score', '-')}"
            ])
        
        table = Table(table_data, colWidths=[0.8*inch, 0.8*inch, 1.5*inch, 1.5*inch, 1.5*inch])
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1976D2')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 10),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.white),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F5F5F5')])
        ]))
        
        story.append(table)
        story.append(Spacer(1, 30))
        
        return story
    
    def _create_rf_section(self, rf_results: Dict) -> List:
        """Create RF prediction section"""
        story = []
        
        story.append(Paragraph("RF Predictions", self.styles['Heading2']))
        story.append(Spacer(1, 12))
        
        predictions = rf_results.get("predictions", [])
        
        table_data = [["Pose", "Predicted pKd"]]
        
        for pred in predictions[:5]:
            table_data.append([
                str(pred.get("pose_id", "-")),
                f"{pred.get('predicted_pKd', '-')}"
            ])
        
        table = Table(table_data, colWidths=[2*inch, 2*inch])
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#388E3C')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey)
        ]))
        
        story.append(table)
        story.append(Spacer(1, 20))
        
        return story
    
    def _create_mmgbsa_section(self, mmgbsa_result: Dict) -> List:
        """Create MM-GBSA section"""
        story = []
        
        story.append(Paragraph("MM-GBSA Refinement", self.styles['Heading2']))
        story.append(Spacer(1, 12))
        
        story.append(Paragraph(
            "MM-GBSA (Molecular Mechanics Generalized Born/Surface Area) refinement provides more accurate binding energy estimates by considering solvation effects.",
            self.styles['Normal']
        ))
        story.append(Spacer(1, 12))
        
        results = mmgbsa_result.get("results", [])
        
        if results:
            table_data = [["Pose", "MM-GBSA Energy (kcal/mol)"]]
            
            for r in results[:5]:
                table_data.append([
                    str(r.get("pose_id", "-")),
                    f"{r.get('mmgbsa_energy', '-')}"
                ])
            
            table = Table(table_data, colWidths=[2*inch, 2.5*inch])
            table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#7B1FA2')),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.grey)
            ]))
            
            story.append(table)
        
        story.append(Spacer(1, 20))
        
        return story
    
    def _create_footer(self) -> List:
        """Create report footer"""
        story = []
        
        story.append(Spacer(1, 20))
        story.append(Paragraph("=" * 60, self.styles['Normal']))
        
        footer_style = ParagraphStyle(
            'Footer',
            parent=self.styles['Normal'],
            fontSize=8,
            textColor=colors.grey
        )
        
        story.append(Paragraph(
            "Generated by BioDockify Docking Studio - Agent Zero AI System",
            footer_style
        ))
        
        return story
