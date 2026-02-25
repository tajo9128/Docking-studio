"""
BioDockify Docking Studio - Results Widget
Displays docking results and analysis
"""

from PyQt6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QTabWidget, QTableWidget, QTableWidgetItem, QHeaderView, QTextEdit, QFrame, QComboBox, QAbstractItemView
from PyQt6.QtCore import pyqtSignal, Qt
from PyQt6.QtGui import QPalette
import logging
import json

logger = logging.getLogger(__name__)

class ResultsWidget(QWidget):
    """Results display widget for docking jobs"""

    # Signals
    results_exported = pyqtSignal(str, str)  # format, content
    new_job_requested = pyqtSignal()  # User wants to start new job
    view_3d_requested = pyqtSignal(object)  # 3D viewer widget

    def __init__(self):
        """Initialize results widget"""
        super().__init__()
        self.current_results = None
        self._setup_ui()

    def _setup_ui(self) -> None:
        """Setup UI components"""
        layout = QVBoxLayout(self)

        # Tabs for different views
        self.tabs = QTabWidget()
        self.tabs.setStyleSheet("""
            QTabWidget::pane {
                border: 1px solid #CCCCCC;
                border-top: 2px solid #E0E0E0;
                background-color: white;
            }
            QTabBar::tab {
                background-color: #F0F0F0;
                color: #333;
                padding: 10px 20px;
                border-top-left-radius: 4px;
                border-top-right-radius: 4px;
                font-weight: bold;
                font-size: 13px;
            }
            QTabBar::tab:selected {
                background-color: #2196F3;
                color: white;
                border-bottom: 2px solid #2196F3;
            }
        """)

        # Summary tab
        self.summary_tab = QWidget()
        self._setup_summary_tab()
        self.tabs.addTab(self.summary_tab, "Summary")

        # Interactions tab
        self.interactions_tab = QWidget()
        self._setup_interactions_tab()
        self.tabs.addTab(self.interactions_tab, "Interactions")

        # Descriptors tab
        self.descriptors_tab = QWidget()
        self._setup_descriptors_tab()
        self.tabs.addTab(self.descriptors_tab, "Descriptors")

        # Agent Zero tab
        self.agent_zero_tab = QWidget()
        self._setup_agent_zero_tab()
        self.tabs.addTab(self.agent_zero_tab, "Agent Zero")

        # Add tabs to layout
        layout.addWidget(self.tabs)

        # Export section
        export_frame = QFrame()
        export_frame.setStyleSheet("border: 2px solid #E0E0E0; border-radius: 8px; padding: 10px;")
        export_frame.setMaximumHeight(80)

        export_label = QLabel("Export Results:")
        export_label.setStyleSheet("font-size: 12px; font-weight: bold; color: #2196F3;")

        export_layout = QHBoxLayout(export_frame)

        self.csv_button = QPushButton("Export CSV")
        self.csv_button.setStyleSheet("""
            QPushButton {
                background-color: #2196F3;
                color: white;
                border: none;
                padding: 8px 15px;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #0B7CD5;
            }
            QPushButton:pressed {
                background-color: #0A5F8D;
            }
        """)
        self.csv_button.clicked.connect(lambda: self.export_results("csv"))

        self.json_button = QPushButton("Export JSON")
        self.json_button.setStyleSheet("""
            QPushButton {
                background-color: #2196F3;
                color: white;
                border: none;
                padding: 8px 15px;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #0B7CD5;
            }
            QPushButton:pressed {
                background-color: #0A5F8D;
            }
        """)
        self.json_button.clicked.connect(lambda: self.export_results("json"))

        self.pymol_button = QPushButton("Export PyMOL")
        self.pymol_button.setStyleSheet("""
            QPushButton {
                background-color: #2196F3;
                color: white;
                border: none;
                padding: 8px 15px;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #0B7CD5;
            }
            QPushButton:pressed {
                background-color: #0A5F8D;
            }
        """)
        self.pymol_button.clicked.connect(lambda: self.export_results("pymol"))

        self.view_3d_button = QPushButton("View in 3D")
        self.view_3d_button.setStyleSheet("""
            QPushButton {
                background-color: #9C27B0;
                color: white;
                border: none;
                padding: 8px 15px;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #7B1FA2;
            }
            QPushButton:pressed {
                background-color: #6A1B9A;
            }
        """)
        self.view_3d_button.clicked.connect(self.on_view_in_3d)

        export_layout.addWidget(export_label)
        export_layout.addWidget(self.csv_button)
        export_layout.addWidget(self.json_button)
        export_layout.addWidget(self.pymol_button)
        export_layout.addWidget(self.view_3d_button)
        export_layout.addStretch()

        # New job button
        self.new_job_button = QPushButton("New Job")
        self.new_job_button.setStyleSheet("""
            QPushButton {
                background-color: #4CAF50;
                color: white;
                border: none;
                padding: 10px 20px;
                border-radius: 5px;
                font-weight: bold;
                font-size: 14px;
            }
            QPushButton:hover {
                background-color: #45A049;
            }
            QPushButton:pressed {
                background-color: #3D8B40;
            }
        """)
        self.new_job_button.clicked.connect(self.on_new_job)

        new_job_layout = QHBoxLayout()
        new_job_layout.addStretch()
        new_job_layout.addWidget(self.new_job_button)

        # Add frames to layout
        layout.addWidget(export_frame)
        layout.addWidget(new_job_layout)

        self.setLayout(layout)

    def _setup_summary_tab(self) -> None:
        """Setup summary tab"""
        layout = QVBoxLayout(self.summary_tab)

        # Binding energy section
        energy_frame = QFrame()
        energy_frame.setStyleSheet("border: 2px solid #E0E0E0; border-radius: 8px; padding: 10px;")

        energy_label = QLabel("Binding Energy")
        energy_label.setStyleSheet("font-size: 14px; font-weight: bold; color: #2196F3;")

        self.energy_display = QLabel("-- kcal/mol")
        self.energy_display.setStyleSheet("font-size: 24px; font-weight: bold; color: #4CAF50; padding: 10px;")

        energy_layout = QHBoxLayout(energy_frame)
        energy_layout.addWidget(energy_label)
        energy_layout.addWidget(self.energy_display)
        energy_layout.addStretch()

        # Confidence score section
        confidence_frame = QFrame()
        confidence_frame.setStyleSheet("border: 2px solid #E0E0E0; border-radius: 8px; padding: 10px;")

        confidence_label = QLabel("Confidence Score")
        confidence_label.setStyleSheet("font-size: 14px; font-weight: bold; color: #2196F3;")

        self.confidence_display = QLabel("-- / 100")
        self.confidence_display.setStyleSheet("font-size: 24px; font-weight: bold; color: #4CAF50; padding: 10px;")

        self.confidence_level = QLabel("--")
        self.confidence_level.setStyleSheet("font-size: 16px; color: #666; padding: 5px;")

        confidence_layout = QVBoxLayout(confidence_frame)
        confidence_layout.addWidget(self.confidence_display)
        confidence_layout.addWidget(self.confidence_level)

        confidence_header = QHBoxLayout(confidence_frame)
        confidence_header.addWidget(confidence_label)
        confidence_header.addStretch()

        # Job info section
        job_info_frame = QFrame()
        job_info_frame.setStyleSheet("border: 2px solid #E0E0E0; border-radius: 8px; padding: 10px;")

        job_info_label = QLabel("Job Information")
        job_info_label.setStyleSheet("font-size: 14px; font-weight: bold; color: #2196F3;")

        self.job_id_display = QLabel("No results available")
        self.job_id_display.setStyleSheet("color: #666;")

        self.num_modes_display = QLabel("N/A")
        self.num_modes_display.setStyleSheet("color: #666;")

        job_info_layout = QVBoxLayout(job_info_frame)
        job_info_layout.addWidget(job_info_label)
        job_info_layout.addWidget(self.job_id_display)
        job_info_layout.addWidget(self.num_modes_display)

        layout.addWidget(energy_frame)
        layout.addWidget(confidence_frame)
        layout.addWidget(job_info_frame)
        layout.addStretch()

    def _setup_interactions_tab(self) -> None:
        """Setup interactions tab"""
        layout = QVBoxLayout(self.interactions_tab)

        self.interactions_table = QTableWidget()
        self.interactions_table.setStyleSheet("""
            QTableWidget {
                background-color: white;
                border: 1px solid #CCCCCC;
                gridline-color: #E0E0E0;
            }
            QTableWidget::item {
                padding: 5px;
                border: none;
            }
            QTableWidget::item:selected {
                background-color: #E3F2FD;
                color: white;
            }
            QHeaderView::section {
                background-color: #F5F5F5;
                color: #333;
                font-weight: bold;
                padding: 5px;
                border: none;
                border-right: 1px solid #E0E0E0;
                border-bottom: 1px solid #E0E0E0;
            }
        """)
        self.interactions_table.setColumnCount(4)
        self.interactions_table.setHorizontalHeaderLabels(["Type", "Atoms", "Distance (Å)", "Strength"])
        self.interactions_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.interactions_table.verticalHeader().setVisible(False)
        self.interactions_table.setEditTriggers(QAbstractItemView.NoEditTriggers)

        layout.addWidget(self.interactions_table)

        self.no_results_label = QLabel("No results available")
        self.no_results_label.setStyleSheet("color: #666; font-style: italic;")
        self.no_results_label.setVisible(False)

        layout.addWidget(self.no_results_label)
        layout.addStretch()

    def _setup_descriptors_tab(self) -> None:
        """Setup descriptors tab"""
        layout = QVBoxLayout(self.descriptors_tab)

        self.descriptors_table = QTableWidget()
        self.descriptors_table.setStyleSheet("""
            QTableWidget {
                background-color: white;
                border: 1px solid #CCCCCC;
                gridline-color: #E0E0E0;
            }
            QTableWidget::item {
                padding: 5px;
                border: none;
            }
            QTableWidget::item:selected {
                background-color: #E3F2FD;
                color: white;
            }
            QHeaderView::section {
                background-color: #F5F5F5;
                color: #333;
                font-weight: bold;
                padding: 5px;
                border: none;
                border-right: 1px solid #E0E0E0;
                border-bottom: 1px solid #E0E0E0;
            }
        """)
        self.descriptors_table.setColumnCount(2)
        self.descriptors_table.setHorizontalHeaderLabels(["Descriptor", "Value"])
        self.descriptors_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.descriptors_table.verticalHeader().setVisible(False)
        self.descriptors_table.setEditTriggers(QAbstractItemView.NoEditTriggers)

        layout.addWidget(self.descriptors_table)

        self.no_results_label_descriptors = QLabel("No results available")
        self.no_results_label_descriptors.setStyleSheet("color: #666; font-style: italic;")
        self.no_results_label_descriptors.setVisible(False)

        layout.addWidget(self.no_results_label_descriptors)
        layout.addStretch()

    def _setup_agent_zero_tab(self) -> None:
        """Setup Agent Zero tab"""
        layout = QVBoxLayout(self.agent_zero_tab)

        # Confidence explanation
        confidence_explanation = QLabel("""
        <h3 style="color: #2196F3; margin: 0 0 10px 0;">Confidence Score Explanation</h3>
        <p style="margin: 0 0 10px 0;">
        BioDockify Docking Studio provides a confidence score (0-100) to indicate the reliability of docking results.
        </p>
        <ul style="margin: 0 0 10px 20px;">
            <li style="margin: 5px 0;"><strong style="color: #4CAF50;">High (80-100)</strong>: Results highly reliable. No repairs needed.</li>
            <li style="margin: 5px 0;"><strong style="color: #FFC107;">Medium (60-79)</strong>: Results reliable. Some repairs performed and validated.</li>
            <li style="margin: 5px 0;"><strong style="color: #F44336;">Low (<60)</strong>: Results may be less reliable. Multiple repairs or graceful degradation occurred.</li>
        </ul>
        <p style="margin: 0 0 10px 0;">
        Confidence scores are adjusted based on the number and severity of repairs performed during the docking process.
        </p>
        """)
        confidence_explanation.setWordWrap(True)
        confidence_explanation.setTextFormat(Qt.RichText)
        confidence_explanation.setStyleSheet("background-color: #F5F5F5; color: #333; padding: 15px; border-radius: 8px;")

        layout.addWidget(confidence_explanation)

        # Failures and repairs section
        failures_frame = QFrame()
        failures_frame.setStyleSheet("border: 2px solid #E0E0E0; border-radius: 8px; padding: 10px;")

        failures_label = QLabel("Failures and Repairs")
        failures_label.setStyleSheet("font-size: 14px; font-weight: bold; color: #2196F3;")

        self.failures_text = QTextEdit()
        self.failures_text.setReadOnly(True)
        self.failures_text.setMaximumHeight(200)
        self.failures_text.setStyleSheet("""
            QTextEdit {
                background-color: #F5F5F5;
                color: #333;
                border: 1px solid #E0E0E0;
                border-radius: 5px;
                padding: 5px;
                font-family: monospace;
                font-size: 10px;
            }
        """)

        failures_layout = QVBoxLayout(failures_frame)
        failures_layout.addWidget(failures_label)
        failures_layout.addWidget(self.failures_text)

        layout.addWidget(failures_frame)
        layout.addStretch()

    def display_results(self, results: dict) -> None:
        """Display docking results"""
        self.current_results = results

        # Update summary tab
        self.energy_display.setText(f"{results.get('binding_energy', '--:--')} kcal/mol")
        self.confidence_display.setText(f"{results.get('confidence_score', '--')} / 100")
        self.job_id_display.setText(results.get('job_id', 'N/A'))
        self.num_modes_display.setText(str(results.get('num_modes', 0)))

        # Update confidence level
        confidence_score = results.get('confidence_score', 0)
        if confidence_score >= 80:
            confidence_level = "High"
            self.confidence_display.setStyleSheet("font-size: 24px; font-weight: bold; color: #4CAF50; padding: 10px;")
        elif confidence_score >= 60:
            confidence_level = "Medium"
            self.confidence_display.setStyleSheet("font-size: 24px; font-weight: bold; color: #FFC107; padding: 10px;")
        else:
            confidence_level = "Low"
            self.confidence_display.setStyleSheet("font-size: 24px; font-weight: bold; color: #F44336; padding: 10px;")

        self.confidence_level.setText(f"Confidence Level: {confidence_level}")

        # Update interactions tab
        interactions = results.get('interactions', {})
        self.interactions_table.setRowCount(0)

        if interactions and interactions.get('total_count', 0) > 0:
            self.no_results_label.setVisible(False)
            self.interactions_table.setVisible(True)

            row = 0

            # Hydrogen bonds
            for hbond in interactions.get('hydrogen_bonds', []):
                self.interactions_table.insertRow(row)
                self.interactions_table.setItem(row, 0, QTableWidgetItem("Hydrogen Bond"))
                self.interactions_table.setItem(row, 1, QTableWidgetItem(f"{hbond['atom_a']} - {hbond['atom_b']}"))
                self.interactions_table.setItem(row, 2, QTableWidgetItem(f"{hbond['distance']}"))
                self.interactions_table.setItem(row, 3, QTableWidgetItem(hbond['strength']))
                row += 1

            # Hydrophobic contacts
            for hydro in interactions.get('hydrophobic_contacts', []):
                self.interactions_table.insertRow(row)
                self.interactions_table.setItem(row, 0, QTableWidgetItem("Hydrophobic"))
                self.interactions_table.setItem(row, 1, QTableWidgetItem(f"{hydro['atom_a']} - {hydro['atom_b']}"))
                self.interactions_table.setItem(row, 2, QTableWidgetItem(f"{hydro['distance']}"))
                self.interactions_table.setItem(row, 3, QTableWidgetItem("-"))
                row += 1

            # Pi-stacking
            for pi_stack in interactions.get('pi_stacking', []):
                self.interactions_table.insertRow(row)
                self.interactions_table.setItem(row, 0, QTableWidgetItem("Pi-Stacking"))
                self.interactions_table.setItem(row, 1, QTableWidgetItem(f"{pi_stack['atom_a']} - {pi_stack['atom_b']}"))
                self.interactions_table.setItem(row, 2, QTableWidgetItem(f"{pi_stack['distance']}"))
                self.interactions_table.setItem(row, 3, QTableWidgetItem("-"))
                row += 1

            # Halogen bonds
            for halogen in interactions.get('halogen_bonds', []):
                self.interactions_table.insertRow(row)
                self.interactions_table.setItem(row, 0, QTableWidgetItem("Halogen Bond"))
                self.interactions_table.setItem(row, 1, QTableWidgetItem(f"{halogen['atom_a']} - {halogen['atom_b']}"))
                self.interactions_table.setItem(row, 2, QTableWidgetItem(f"{halogen['distance']}"))
                self.interactions_table.setItem(row, 3, QTableWidgetItem(halogen['strength']))
                row += 1

            # Salt bridges
            for salt_bridge in interactions.get('salt_bridges', []):
                self.interactions_table.insertRow(row)
                self.interactions_table.setItem(row, 0, QTableWidgetItem("Salt Bridge"))
                self.interactions_table.setItem(row, 1, QTableWidgetItem(f"{salt_bridge['atom_a']} - {salt_bridge['atom_b']}"))
                self.interactions_table.setItem(row, 2, QTableWidgetItem(f"{salt_bridge['distance']}"))
                self.interactions_table.setItem(row, 3, QTableWidgetItem("-"))
                row += 1

            # Cation-pi
            for cation_pi in interactions.get('cation_pi_interactions', []):
                self.interactions_table.insertRow(row)
                self.interactions_table.setItem(row, 0, QTableWidgetItem("Cation-π"))
                self.interactions_table.setItem(row, 1, QTableWidgetItem(f"{cation_pi['atom_a']} - {cation_pi['atom_b']}"))
                self.interactions_table.setItem(row, 2, QTableWidgetItem(f"{cation_pi['distance']}"))
                self.interactions_table.setItem(row, 3, QTableWidgetItem("-"))
                row += 1

            # Metal coordination
            for metal_coord in interactions.get('metal_coordination', []):
                self.interactions_table.insertRow(row)
                self.interactions_table.setItem(row, 0, QTableWidgetItem("Metal Coordination"))
                self.interactions_table.setItem(row, 1, QTableWidgetItem(f"{metal_coord['atom_a']} - {metal_coord['atom_b']}"))
                self.interactions_table.setItem(row, 2, QTableWidgetItem(f"{metal_coord['distance']}"))
                self.interactions_table.setItem(row, 3, QTableWidgetItem("-"))
                row += 1

            self.interactions_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
            self.interactions_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
            self.interactions_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
            self.interactions_table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeToContents)

        else:
            self.no_results_label.setVisible(True)
            self.interactions_table.setVisible(False)

        # Update descriptors tab
        descriptors = results.get('descriptors', {})
        self.descriptors_table.setRowCount(0)

        if descriptors:
            self.descriptors_table.setVisible(True)
            self.no_results_label_descriptors.setVisible(False)

            row = 0

            for key, value in descriptors.items():
                self.descriptors_table.insertRow(row)
                self.descriptors_table.setItem(row, 0, QTableWidgetItem(key))
                self.descriptors_table.setItem(row, 1, QTableWidgetItem(str(value)))
                row += 1

            self.descriptors_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
            self.descriptors_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)

        else:
            self.no_results_label_descriptors.setVisible(True)
            self.descriptors_table.setVisible(False)

        # Update Agent Zero tab
        failures_detected = results.get('failures_detected', [])
        repairs_attempted = results.get('repairs_attempted', [])

        agent_zero_text = ""

        if failures_detected:
            agent_zero_text += "<h4 style='color: #F44336; margin: 0 0 10px 0;'>Failures Detected:</h4>"
            agent_zero_text += "<ul>"
            for failure in failures_detected:
                agent_zero_text += f"<li style='margin: 5px 0;'>{failure}</li>"
            agent_zero_text += "</ul>"

        if repairs_attempted:
            agent_zero_text += "<h4 style='color: #4CAF50; margin: 10px 0 0;'>Repairs Attempted:</h4>"
            agent_zero_text += "<ul>"
            for repair in repairs_attempted:
                agent_zero_text += f"<li style='margin: 5px 0;'>{repair}</li>"
            agent_zero_text += "</ul>"

        if not failures_detected and not repairs_attempted:
            agent_zero_text += "<p style='color: #666; font-style: italic;'>No failures or repairs occurred. Results are from first-pass success.</p>"

        self.failures_text.setText(agent_zero_text)

        logger.info(f"Results displayed for job: {results.get('job_id', 'unknown')}")

    def export_results(self, format: str) -> None:
        """Export results to file"""
        if not self.current_results:
            logger.warning("No results to export")
            return

        try:
            if format == "csv":
                content = self._generate_csv()
            elif format == "json":
                content = json.dumps(self.current_results, indent=2)
            elif format == "pymol":
                content = self._generate_pymol()
            else:
                logger.error(f"Unknown export format: {format}")
                return

            self.results_exported.emit(format, content)
            logger.info(f"Results exported as {format}")

        except Exception as e:
            logger.error(f"Failed to export results: {e}")

    def _generate_csv(self) -> str:
        """Generate CSV content"""
        lines = []
        lines.append("Metric,Value")

        # Summary data
        lines.append(f"Binding Energy (kcal/mol),{self.current_results.get('binding_energy', '')}")
        lines.append(f"Confidence Score,{self.current_results.get('confidence_score', '')}")

        # Interactions data
        interactions = self.current_results.get('interactions', {})
        if interactions:
            lines.append(f"Total Interactions,{interactions.get('total_count', 0)}")

            lines.append("Interaction Type,Count")
            lines.append(f"Hydrogen Bonds,{len(interactions.get('hydrogen_bonds', []))}")
            lines.append(f"Hydrophobic Contacts,{len(interactions.get('hydrophobic_contacts', []))}")
            lines.append(f"Pi-Stacking,{len(interactions.get('pi_stacking', []))}")
            lines.append(f"Halogen Bonds,{len(interactions.get('halogen_bonds', []))}")
            lines.append(f"Salt Bridges,{len(interactions.get('salt_bridges', []))}")
            lines.append(f"Cation-π,{len(interactions.get('cation_pi_interactions', []))}")
            lines.append(f"Metal Coordination,{len(interactions.get('metal_coordination', []))}")

        # Descriptors data
        descriptors = self.current_results.get('descriptors', {})
        if descriptors:
            lines.append("")
            lines.append("Descriptor,Value")

            for key, value in descriptors.items():
                lines.append(f"{key},{value}")

        return "\n".join(lines)

    def _generate_pymol(self) -> str:
        """Generate PyMOL content"""
        # In real implementation, would generate PyMOL script to load ligand and visualize interactions
        pymol_script = """
# PyMOL script to visualize BioDockify Docking Studio results
# Load this script in PyMOL to visualize docking results

# Load ligand (modify path as needed)
load your_ligand_file.pdb

# Set background to white
bg_color white

# Show binding interactions
# Add visualization commands here based on BioDockify interactions
"""

        return pymol_script

    def on_view_in_3d(self) -> None:
        """Open 3D viewer with colorful visualization of docking results"""
        if not self.current_results:
            logger.warning("No results available for 3D visualization")
            return
        
        try:
            from src.visualization import AdvancedMolecularViewer
            
            receptor_pdb = self.current_results.get('receptor_pdb')
            ligand_pdb = self.current_results.get('ligand_pdb')
            poses = self.current_results.get('poses', [])
            interactions = self.current_results.get('interactions', {})
            
            viewer = AdvancedMolecularViewer()
            viewer.setWindowTitle("Docking Results - 3D Visualization")
            viewer.resize(1200, 800)
            
            if receptor_pdb:
                viewer.load_receptor(receptor_pdb)
                logger.info("Receptor loaded in 3D viewer")
            
            if ligand_pdb:
                viewer.load_ligand(ligand_pdb, pose_id=0)
                viewer.color_ligand("element")
                logger.info("Ligand loaded in 3D viewer with CPK coloring")
            
            for i, pose_data in enumerate(poses[:5]):
                if isinstance(pose_data, dict):
                    pose_pdb = pose_data.get('pdb')
                    score = pose_data.get('score', 0.0)
                else:
                    pose_pdb = pose_data
                    score = 0.0
                
                if pose_pdb:
                    pose_colors = ['cyan', 'magenta', 'green', 'orange', 'pink']
                    viewer.add_pose(pose_pdb, i, score)
            
            viewer.apply_professional_style()
            
            if interactions and interactions.get('total_count', 0) > 0:
                interaction_list = []
                
                for hbond in interactions.get('hydrogen_bonds', []):
                    interaction_list.append({
                        'type': 'hbond',
                        'atom1_coords': hbond.get('atom1_coords', [0, 0, 0]),
                        'atom2_coords': hbond.get('atom2_coords', [0, 0, 0])
                    })
                
                for hydro in interactions.get('hydrophobic_contacts', []):
                    interaction_list.append({
                        'type': 'hydrophobic',
                        'atom1_coords': hydro.get('atom1_coords', [0, 0, 0]),
                        'atom2_coords': hydro.get('atom2_coords', [0, 0, 0])
                    })
                
                if interaction_list:
                    viewer.set_interactions(interaction_list)
                    viewer.btn_interactions.setChecked(True)
            
            viewer.show()
            
            self.view_3d_requested.emit(viewer)
            logger.info("3D visualization window opened")
            
        except Exception as e:
            logger.error(f"Failed to open 3D viewer: {e}")

    def export_colorful_image(self, filename: str = "docking_result.png") -> bool:
        """Export colorful 3D image directly from results"""
        if not self.current_results:
            logger.warning("No results available for image export")
            return False
        
        try:
            from src.visualization import AdvancedMolecularViewer
            import base64
            import os
            
            receptor_pdb = self.current_results.get('receptor_pdb')
            ligand_pdb = self.current_results.get('ligand_pdb')
            
            viewer = AdvancedMolecularViewer()
            
            if receptor_pdb:
                viewer.load_receptor(receptor_pdb)
            
            if ligand_pdb:
                viewer.load_ligand(ligand_pdb, pose_id=0)
                viewer.color_ligand("element")
            
            viewer.apply_professional_style()
            
            def save_image(data_uri):
                if data_uri and data_uri.startswith('data:image/png;base64,'):
                    data = base64.b64decode(data_uri.split(',')[1])
                    with open(filename, 'wb') as f:
                        f.write(data)
                    logger.info(f"Colorful image exported to {filename}")
            
            viewer.take_screenshot(resolution=2, callback=save_image)
            return True
            
        except Exception as e:
            logger.error(f"Failed to export colorful image: {e}")
            return False

    def on_new_job(self) -> None:
        """Handle new job button click"""
        logger.info("New job button clicked")
        self.new_job_requested.emit()

    def clear_results(self) -> None:
        """Clear results display"""
        self.current_results = None
        self.energy_display.setText("-- kcal/mol")
        self.confidence_display.setText("-- / 100")
        self.confidence_level.setText("--")
        self.job_id_display.setText("No results available")
        self.num_modes_display.setText("N/A")
        self.interactions_table.setRowCount(0)
        self.descriptors_table.setRowCount(0)
        self.failures_text.clear()
        logger.info("Results cleared")
