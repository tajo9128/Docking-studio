"""
BioDockify Docking Studio - Professional Main Window
International-quality UI implementation
"""

from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QStackedWidget, QFrame,
    QProgressBar, QScrollArea, QSizePolicy, QDockWidget,
    QToolBar, QStatusBar, QMenu, QMenuBar,
    QSplitter, QTabWidget, QSpinBox, QDoubleSpinBox, QGridLayout
)
from PyQt6.QtCore import Qt, QSize, QTimer, pyqtSignal, QPoint, QPropertyAnimation, QEasingCurve
from PyQt6.QtGui import QIcon, QFont, QAction, QPalette, QColor, QPixmap, QPainter, QLinearGradient
from PyQt6.QtSvgWidgets import QSvgWidget
import logging

from src.ui.theme import DesignTokens, Stylesheet
from src.ui.chat_widget import ChatWidget
from src.services.llm_chat_service import ChatService, ChatConfig, ModelProvider

logger = logging.getLogger(__name__)

class MainWindow(QMainWindow):
    """Professional main application window"""
    
    # Signals
    job_started = pyqtSignal(str)
    job_completed = pyqtSignal(str)
    job_failed = pyqtSignal(str, str)
    job_cancelled = pyqtSignal(str)
    job_progress = pyqtSignal(str, int)  # job_id, percentage
    
    def __init__(self):
        """Initialize main window with professional styling"""
        super().__init__()
        
        # Setup professional styling
        self._apply_professional_theme()
        
        # Window setup
        self.setWindowTitle("BioDockify Docking Studio")
        self.setMinimumSize(1400, 900)
        self.resize(1600, 1000)
        self.setWindowIcon(QIcon("src/ui/styles/icon.ico"))
        
        # Central widget with modern layout
        self._setup_central_widget()
        
        # Professional toolbar
        self._setup_toolbar()
        
        # Professional menu bar
        self._create_menu_bar()
        
        # Professional status bar
        self._create_status_bar()
        
        # Setup chat dock widget
        self._setup_chat_dock()
        
        logger.info("Professional main window initialized")
    
    def _apply_professional_theme(self):
        """Apply professional international styling"""
        # Set application stylesheet
        self.setStyleSheet(Stylesheet.get_application_stylesheet())
        
        # Set modern font
        font = QFont()
        font.setFamily(DesignTokens.Typography.PRIMARY)
        font.setPointSize(DesignTokens.Typography.BODY_L)
        self.setFont(font)
        
        # Set modern palette
        palette = QPalette()
        palette.setColor(QPalette.ColorRole.Window, QColor(DesignTokens.Colors.WHITE))
        palette.setColor(QPalette.ColorRole.WindowText, QColor(DesignTokens.Colors.TEXT_PRIMARY))
        palette.setColor(QPalette.ColorRole.Base, QColor(DesignTokens.Colors.WHITE))
        palette.setColor(QPalette.ColorRole.Text, QColor(DesignTokens.Colors.TEXT_PRIMARY))
        palette.setColor(QPalette.ColorRole.Button, QColor(DesignTokens.Colors.PRIMARY))
        palette.setColor(QPalette.ColorRole.ButtonText, QColor(DesignTokens.Colors.WHITE))
        self.setPalette(palette)
    
    def _setup_central_widget(self):
        """Setup central widget with professional layout"""
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout with proper spacing
        main_layout = QVBoxLayout(central_widget)
        main_layout.setSpacing(DesignTokens.Spacing.GAP_LG)
        main_layout.setContentsMargins(DesignTokens.Spacing.PADDING_XL, 
                                     DesignTokens.Spacing.PADDING_XL,
                                     DesignTokens.Spacing.PADDING_XL,
                                     DesignTokens.Spacing.PADDING_XL)
        
        # Professional header section
        header_widget = self._create_header()
        main_layout.addWidget(header_widget)
        
        # Content area with splitter
        content_splitter = QSplitter(Qt.Orientation.Horizontal)
        content_splitter.setHandleWidth(1)
        
        # Left sidebar (navigation)
        sidebar_widget = self._create_sidebar()
        content_splitter.addWidget(sidebar_widget)
        
        # Main content area
        main_content_widget = self._create_main_content()
        content_splitter.addWidget(main_content_widget)
        
        # Set splitter ratios
        content_splitter.setStretchFactor(0, 1)
        content_splitter.setStretchFactor(1, 4)
        
        main_layout.addWidget(content_splitter, 1)
    
    def _create_header(self) -> QWidget:
        """Create professional header"""
        header = QWidget()
        header.setProperty("card", "true")
        header.setMaximumHeight(100)
        
        header_layout = QHBoxLayout(header)
        header_layout.setSpacing(DesignTokens.Spacing.GAP_MD)
        header_layout.setContentsMargins(DesignTokens.Spacing.PADDING_LG,
                                       DesignTokens.Spacing.PADDING_LG,
                                       DesignTokens.Spacing.PADDING_LG,
                                       DesignTokens.Spacing.PADDING_LG)
        
        # Logo and title
        logo_section = QWidget()
        logo_layout = QHBoxLayout(logo_section)
        logo_layout.setSpacing(DesignTokens.Spacing.GAP_SM)
        
        # App logo (placeholder - use actual logo)
        logo_label = QLabel("🧬")
        logo_label.setStyleSheet(f"""
            font-size: 48px;
            background: linear-gradient(135deg, {DesignTokens.Colors.PRIMARY}, {DesignTokens.Colors.SECONDARY});
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            padding: 8px;
            border-radius: {DesignTokens.BorderRadius.MD}px;
        """)
        
        # Title and version
        title_section = QWidget()
        title_layout = QVBoxLayout(title_section)
        title_layout.setSpacing(DesignTokens.Spacing.GAP_XS)
        
        title_label = QLabel("BioDockify Docking Studio")
        title_label.setProperty("heading", "true")
        title_label.setStyleSheet(f"""
            QLabel[heading="true"] {{
                font-size: {DesignTokens.Typography.DISPLAY_M}px;
                font-weight: 700;
                color: {DesignTokens.Colors.TEXT_PRIMARY};
            }}
        """)
        
        version_label = QLabel("Version 1.2.2")
        version_label.setProperty("caption", "true")
        version_label.setStyleSheet(f"""
            QLabel[caption="true"] {{
                font-size: {DesignTokens.Typography.CAPTION_M}px;
                color: {DesignTokens.Colors.TEXT_TERTIARY};
                font-weight: 500;
                padding: 4px 12px;
                background: {DesignTokens.Colors.GRAY_100};
                border-radius: {DesignTokens.BorderRadius.FULL}px;
            }}
        """)
        
        title_layout.addWidget(title_label)
        title_layout.addWidget(version_label)
        
        logo_layout.addWidget(logo_label)
        logo_layout.addWidget(title_section)
        logo_layout.addStretch()
        
        # Action buttons
        actions_section = QWidget()
        actions_layout = QHBoxLayout(actions_section)
        actions_layout.setSpacing(DesignTokens.Spacing.GAP_SM)
        
        new_job_button = QPushButton("New Job")
        new_job_button.setProperty("icon", "➕")
        new_job_button.setStyleSheet(f"""
            QPushButton {{
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                    stop:0 {DesignTokens.Colors.PRIMARY},
                    stop:1 {DesignTokens.Colors.PRIMARY_HOVER});
                color: {DesignTokens.Colors.WHITE};
                border: none;
                border-radius: {DesignTokens.BorderRadius.MD}px;
                padding: {DesignTokens.Spacing.PADDING_SM}px {DesignTokens.Spacing.PADDING_LG}px;
                font-size: {DesignTokens.Typography.BODY_L}px;
                font-weight: 600;
                min-height: {DesignTokens.Spacing.XL}px;
            }}
            QPushButton:hover {{
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                    stop:0 {DesignTokens.Colors.PRIMARY_HOVER},
                    stop:1 {DesignTokens.Colors.PRIMARY_ACTIVE});
                box-shadow: {DesignTokens.Shadows.MD};
            }}
        """)
        new_job_button.clicked.connect(self.on_new_job)
        
        open_recent_button = QPushButton("Open Recent")
        open_recent_button.setProperty("secondary", "true")
        open_recent_button.clicked.connect(self.on_open_recent)
        
        settings_button = QPushButton("⚙")
        settings_button.setMaximumSize(50, 50)
        settings_button.setProperty("secondary", "true")
        settings_button.clicked.connect(self.on_preferences)
        
        actions_layout.addWidget(new_job_button)
        actions_layout.addWidget(open_recent_button)
        actions_layout.addWidget(settings_button)
        
        # Status indicator
        status_widget = self._create_status_indicator()
        
        header_layout.addWidget(logo_section)
        header_layout.addWidget(actions_section)
        header_layout.addWidget(status_widget)
        
        return header
    
    def _create_status_indicator(self) -> QWidget:
        """Create professional status indicator"""
        status_widget = QWidget()
        status_widget.setFixedWidth(200)
        
        status_layout = QVBoxLayout(status_widget)
        status_layout.setSpacing(DesignTokens.Spacing.GAP_XS)
        
        status_label = QLabel("System Status")
        status_label.setProperty("caption", "true")
        
        # Status badge
        status_badge = QLabel("Ready")
        status_badge.setStyleSheet(f"""
            QLabel {{
                background: {DesignTokens.Colors.SUCCESS_BG};
                color: {DesignTokens.Colors.SUCCESS};
                padding: 8px 16px;
                border-radius: {DesignTokens.BorderRadius.FULL}px;
                font-size: {DesignTokens.Typography.CAPTION_L}px;
                font-weight: 600;
                text-align: center;
                border: 1px solid {DesignTokens.Colors.SUCCESS};
            }}
        """)
        
        status_layout.addWidget(status_label, 0, Qt.AlignmentFlag.AlignRight)
        status_layout.addWidget(status_badge, 0, Qt.AlignmentFlag.AlignRight)
        
        return status_widget
    
    def _create_sidebar(self) -> QWidget:
        """Create professional sidebar navigation"""
        sidebar = QWidget()
        sidebar.setMinimumWidth(250)
        sidebar.setMaximumWidth(250)
        sidebar.setStyleSheet(f"""
            QWidget {{
                background: {DesignTokens.Colors.BACKGROUND_SECONDARY};
                border-right: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
            }}
        """)
        
        sidebar_layout = QVBoxLayout(sidebar)
        sidebar_layout.setSpacing(DesignTokens.Spacing.GAP_SM)
        sidebar_layout.setContentsMargins(DesignTokens.Spacing.PADDING_MD,
                                       DesignTokens.Spacing.PADDING_MD,
                                       DesignTokens.Spacing.PADDING_MD,
                                       DesignTokens.Spacing.PADDING_MD)
        
        # Navigation sections
        nav_sections = [
            ("Jobs", "📊"),
            ("Molecules", "🧬"),
            ("Results", "📈"),
            ("Analysis", "🔬"),
            ("Tools", "🛠"),
            ("Settings", "⚙"),
        ]
        
        for title, icon in nav_sections:
            nav_button = self._create_nav_button(icon, title)
            sidebar_layout.addWidget(nav_button)
        
        sidebar_layout.addStretch()
        
        # Docker status
        docker_status = self._create_docker_status()
        sidebar_layout.addWidget(docker_status)
        
        return sidebar
    
    def _create_nav_button(self, icon: str, title: str) -> QPushButton:
        """Create professional navigation button"""
        button = QPushButton(f"{icon}  {title}")
        button.setProperty("nav", "true")
        button.setStyleSheet(f"""
            QPushButton[nav="true"] {{
                background: transparent;
                border: none;
                text-align: left;
                padding: {DesignTokens.Spacing.PADDING_MD}px;
                font-size: {DesignTokens.Typography.BODY_L}px;
                color: {DesignTokens.Colors.TEXT_SECONDARY};
                border-radius: {DesignTokens.BorderRadius.MD}px;
                spacing: {DesignTokens.Spacing.GAP_SM}px;
            }}
            QPushButton[nav="true"]:hover {{
                background: {DesignTokens.Colors.WHITE};
                color: {DesignTokens.Colors.PRIMARY};
                box-shadow: {DesignTokens.Shadows.XS};
            }}
            QPushButton[nav="true"]:pressed {{
                background: {DesignTokens.Colors.GRAY_100};
            }}
        """)
        button.setCursor(Qt.CursorShape.PointingHandCursor)
        return button
    
    def _create_docker_status(self) -> QFrame:
        """Create Docker status indicator"""
        docker_frame = QFrame()
        docker_frame.setProperty("card", "true")
        docker_frame.setStyleSheet(f"""
            QFrame[card="true"] {{
                background: {DesignTokens.Colors.WHITE};
                border: 1px solid {DesignTokens.Colors.SUCCESS};
            }}
        """)
        
        docker_layout = QVBoxLayout(docker_frame)
        docker_layout.setSpacing(DesignTokens.Spacing.GAP_SM)
        docker_layout.setContentsMargins(DesignTokens.Spacing.PADDING_MD,
                                       DesignTokens.Spacing.PADDING_MD,
                                       DesignTokens.Spacing.PADDING_MD,
                                       DesignTokens.Spacing.PADDING_MD)
        
        title = QLabel("Docker Status")
        title.setProperty("subheading", "true")
        
        status = QLabel("✓ Running")
        status.setStyleSheet(f"""
            QLabel {{
                color: {DesignTokens.Colors.SUCCESS};
                font-size: {DesignTokens.Typography.BODY_L}px;
                font-weight: 600;
            }}
        """)
        
        docker_layout.addWidget(title)
        docker_layout.addWidget(status)
        
        return docker_frame
    
    def _create_main_content(self) -> QWidget:
        """Create main content area with docking form"""
        content = QWidget()
        
        content_layout = QVBoxLayout(content)
        content_layout.setSpacing(DesignTokens.Spacing.GAP_MD)
        
        content_layout.addWidget(self._create_docking_form(), 1)
        
        return content
    
    def _create_docking_form(self) -> QWidget:
        """Create docking form similar to CloudVina BatchDockingPage"""
        form_widget = QWidget()
        form_layout = QVBoxLayout(form_widget)
        form_layout.setSpacing(DesignTokens.Spacing.GAP_LG)
        
        header = QLabel("Batch Docking Page")
        header.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.DISPLAY_L}px;
                font-weight: 700;
                color: {DesignTokens.Colors.TEXT_PRIMARY};
                padding: {DesignTokens.Spacing.PADDING_MD}px 0;
            }}
        """)
        form_layout.addWidget(header)
        
        desc = QLabel("Deploy massive virtual screening campaigns using our consensus docking engine.")
        desc.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.BODY_L}px;
                color: {DesignTokens.Colors.TEXT_SECONDARY};
                margin-bottom: {DesignTokens.Spacing.PADDING_LG}px;
            }}
        """)
        form_layout.addWidget(desc)
        
        upload_splitter = QSplitter(Qt.Orientation.Horizontal)
        upload_splitter.setHandleWidth(16)
        
        receptor_card = self._create_receptor_upload_card()
        upload_splitter.addWidget(receptor_card)
        
        ligand_card = self._create_ligand_upload_card()
        upload_splitter.addWidget(ligand_card)
        
        form_layout.addWidget(upload_splitter)
        
        engine_card = self._create_engine_selection_card()
        form_layout.addWidget(engine_card)
        
        options_splitter = QSplitter(Qt.Orientation.Horizontal)
        options_splitter.setHandleWidth(16)
        
        grid_card = self._create_grid_config_card()
        options_splitter.addWidget(grid_card)
        
        gpu_card = self._create_gpu_status_card()
        options_splitter.addWidget(gpu_card)
        
        form_layout.addWidget(options_splitter)
        
        submit_layout = QHBoxLayout()
        submit_layout.addStretch()
        
        self.start_docking_btn = QPushButton("Start Experiment")
        self.start_docking_btn.setStyleSheet(f"""
            QPushButton {{
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                    stop:0 {DesignTokens.Colors.PRIMARY},
                    stop:1 {DesignTokens.Colors.SECONDARY});
                color: {DesignTokens.Colors.WHITE};
                border: none;
                border-radius: {DesignTokens.BorderRadius.LG}px;
                padding: {DesignTokens.Spacing.PADDING_LG}px {DesignTokens.Spacing.PADDING_XL}px;
                font-size: {DesignTokens.Typography.BODY_XL}px;
                font-weight: 700;
                min-height: 60px;
                box-shadow: {DesignTokens.Shadows.LG};
            }}
            QPushButton:hover {{
                box-shadow: {DesignTokens.Shadows.XL};
                transform: translateY(-2px);
            }}
            QPushButton:disabled {{
                background: {DesignTokens.Colors.GRAY_300};
                box-shadow: none;
            }}
        """)
        self.start_docking_btn.clicked.connect(self.on_start_docking)
        submit_layout.addWidget(self.start_docking_btn)
        
        form_layout.addLayout(submit_layout)
        
        return form_widget
    
    def _create_receptor_upload_card(self) -> QWidget:
        """Create receptor upload card (Step 1)"""
        card = QFrame()
        card.setStyleSheet(f"""
            QFrame {{
                background: {DesignTokens.Colors.WHITE};
                border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                border-radius: {DesignTokens.BorderRadius.XL}px;
                padding: {DesignTokens.Spacing.PADDING_LG}px;
            }}
        """)
        
        layout = QVBoxLayout(card)
        layout.setSpacing(DesignTokens.Spacing.GAP_MD)
        
        header = QWidget()
        header_layout = QHBoxLayout(header)
        header_layout.setSpacing(DesignTokens.Spacing.GAP_MD)
        
        step_number = QLabel("1")
        step_number.setStyleSheet(f"""
            QLabel {{
                background: {DesignTokens.Colors.PRIMARY_10};
                color: {DesignTokens.Colors.PRIMARY};
                border-radius: {DesignTokens.BorderRadius.LG}px;
                padding: 8px 14px;
                font-size: {DesignTokens.Typography.BODY_L}px;
                font-weight: 700;
            }}
        """)
        
        title_layout = QVBoxLayout()
        title_layout.setSpacing(2)
        
        title = QLabel("Target Receptor")
        title.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.BODY_XL}px;
                font-weight: 700;
                color: {DesignTokens.Colors.TEXT_PRIMARY};
            }}
        """)
        
        subtitle = QLabel("Protein structure (.pdb, .pdbqt)")
        subtitle.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.CAPTION_L}px;
                color: {DesignTokens.Colors.TEXT_TERTIARY};
            }}
        """)
        
        title_layout.addWidget(title)
        title_layout.addWidget(subtitle)
        
        header_layout.addWidget(step_number)
        header_layout.addWidget(title_layout, 1)
        
        layout.addWidget(header)
        
        drop_zone = QFrame()
        drop_zone.setMinimumHeight(150)
        drop_zone.setStyleSheet(f"""
            QFrame {{
                border: 2px dashed {DesignTokens.Colors.BORDER_LIGHT};
                border-radius: {DesignTokens.BorderRadius.LG}px;
                background: {DesignTokens.Colors.BACKGROUND_TERTIARY};
            }}
            QFrame:hover {{
                border-color: {DesignTokens.Colors.PRIMARY};
                background: {DesignTokens.Colors.PRIMARY_5};
            }}
        """)
        
        drop_layout = QVBoxLayout(drop_zone)
        drop_layout.setSpacing(DesignTokens.Spacing.GAP_SM)
        
        upload_icon = QLabel("📁")
        upload_icon.setStyleSheet("font-size: 32px;")
        upload_icon.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        upload_text = QLabel("Upload Receptor")
        upload_text.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.BODY_L}px;
                font-weight: 600;
                color: {DesignTokens.Colors.TEXT_PRIMARY};
            }}
        """)
        upload_text.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        upload_hint = QLabel(".pdb, .mmcif, .mol2")
        upload_hint.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.CAPTION_M}px;
                color: {DesignTokens.Colors.TEXT_TERTIARY};
            }}
        """)
        upload_hint.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        drop_layout.addWidget(upload_icon)
        drop_layout.addWidget(upload_text)
        drop_layout.addWidget(upload_hint)
        
        layout.addWidget(drop_zone)
        
        self.receptor_file_label = QLabel("No file selected")
        self.receptor_file_label.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.CAPTION_M}px;
                color: {DesignTokens.Colors.TEXT_TERTIARY};
                padding: {DesignTokens.Spacing.PADDING_SM}px;
                background: {DesignTokens.Colors.GRAY_100};
                border-radius: {DesignTokens.BorderRadius.MD}px;
            }}
        """)
        self.receptor_file_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.receptor_file_label)
        
        select_receptor_btn = QPushButton("Select Receptor File")
        select_receptor_btn.setStyleSheet(f"""
            QPushButton {{
                background: {DesignTokens.Colors.GRAY_100};
                color: {DesignTokens.Colors.TEXT_PRIMARY};
                border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                border-radius: {DesignTokens.BorderRadius.MD}px;
                padding: {DesignTokens.Spacing.PADDING_MD}px;
                font-weight: 600;
            }}
            QPushButton:hover {{
                background: {DesignTokens.Colors.GRAY_200};
            }}
        """)
        select_receptor_btn.clicked.connect(self.on_select_receptor)
        layout.addWidget(select_receptor_btn)
        
        layout.addStretch()
        
        return card
    
    def _create_ligand_upload_card(self) -> QWidget:
        """Create ligand upload card (Step 2)"""
        card = QFrame()
        card.setStyleSheet(f"""
            QFrame {{
                background: {DesignTokens.Colors.WHITE};
                border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                border-radius: {DesignTokens.BorderRadius.XL}px;
                padding: {DesignTokens.Spacing.PADDING_LG}px;
            }}
        """)
        
        layout = QVBoxLayout(card)
        layout.setSpacing(DesignTokens.Spacing.GAP_MD)
        
        header = QWidget()
        header_layout = QHBoxLayout(header)
        header_layout.setSpacing(DesignTokens.Spacing.GAP_MD)
        
        step_number = QLabel("2")
        step_number.setStyleSheet(f"""
            QLabel {{
                background: {DesignTokens.Colors.SUCCESS_10};
                color: {DesignTokens.Colors.SUCCESS};
                border-radius: {DesignTokens.BorderRadius.LG}px;
                padding: 8px 14px;
                font-size: {DesignTokens.Typography.BODY_L}px;
                font-weight: 700;
            }}
        """)
        
        title_layout = QVBoxLayout()
        title_layout.setSpacing(2)
        
        title = QLabel("Ligand Library")
        title.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.BODY_XL}px;
                font-weight: 700;
                color: {DesignTokens.Colors.TEXT_PRIMARY};
            }}
        """)
        
        subtitle = QLabel("Small molecules, Drugs")
        subtitle.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.CAPTION_L}px;
                color: {DesignTokens.Colors.TEXT_TERTIARY};
            }}
        """)
        
        title_layout.addWidget(title)
        title_layout.addWidget(subtitle)
        
        header_layout.addWidget(step_number)
        header_layout.addWidget(title_layout, 1)
        
        layout.addWidget(header)
        
        toggle_layout = QHBoxLayout()
        toggle_layout.setSpacing(4)
        
        self.files_toggle = QPushButton("Files")
        self.files_toggle.setCheckable(True)
        self.files_toggle.setChecked(True)
        self.files_toggle.setStyleSheet(f"""
            QPushButton {{
                background: {DesignTokens.Colors.PRIMARY};
                color: {DesignTokens.Colors.WHITE};
                border: none;
                border-radius: {DesignTokens.BorderRadius.MD}px;
                padding: 8px 16px;
                font-weight: 600;
            }}
            QPushButton:checked {{
                background: {DesignTokens.Colors.PRIMARY};
            }}
            QPushButton:!checked {{
                background: {DesignTokens.Colors.GRAY_100};
                color: {DesignTokens.Colors.TEXT_SECONDARY};
            }}
        """)
        self.files_toggle.clicked.connect(lambda: self.on_upload_mode_changed('files'))
        
        self.csv_toggle = QPushButton("CSV")
        self.csv_toggle.setCheckable(True)
        self.csv_toggle.setStyleSheet(f"""
            QPushButton {{
                background: {DesignTokens.Colors.GRAY_100};
                color: {DesignTokens.Colors.TEXT_SECONDARY};
                border: none;
                border-radius: {DesignTokens.BorderRadius.MD}px;
                padding: 8px 16px;
                font-weight: 600;
            }}
            QPushButton:checked {{
                background: {DesignTokens.Colors.PRIMARY};
                color: {DesignTokens.Colors.WHITE};
            }}
        """)
        self.csv_toggle.clicked.connect(lambda: self.on_upload_mode_changed('csv'))
        
        toggle_layout.addWidget(self.files_toggle)
        toggle_layout.addWidget(self.csv_toggle)
        toggle_layout.addStretch()
        
        layout.addLayout(toggle_layout)
        
        drop_zone = QFrame()
        drop_zone.setMinimumHeight(150)
        drop_zone.setStyleSheet(f"""
            QFrame {{
                border: 2px dashed {DesignTokens.Colors.BORDER_LIGHT};
                border-radius: {DesignTokens.BorderRadius.LG}px;
                background: {DesignTokens.Colors.BACKGROUND_TERTIARY};
            }}
            QFrame:hover {{
                border-color: {DesignTokens.Colors.SUCCESS};
                background: {DesignTokens.Colors.SUCCESS_5};
            }}
        """)
        
        drop_layout = QVBoxLayout(drop_zone)
        drop_layout.setSpacing(DesignTokens.Spacing.GAP_SM)
        
        upload_icon = QLabel("🧪")
        upload_icon.setStyleSheet("font-size: 32px;")
        upload_icon.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        upload_text = QLabel("Upload Ligands")
        upload_text.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.BODY_L}px;
                font-weight: 600;
                color: {DesignTokens.Colors.TEXT_PRIMARY};
            }}
        """)
        upload_text.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        upload_hint = QLabel(".sdf, .mol2, .smi")
        upload_hint.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.CAPTION_M}px;
                color: {DesignTokens.Colors.TEXT_TERTIARY};
            }}
        """)
        upload_hint.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        drop_layout.addWidget(upload_icon)
        drop_layout.addWidget(upload_text)
        drop_layout.addWidget(upload_hint)
        
        layout.addWidget(drop_zone)
        
        self.ligand_count_label = QLabel("0 files prepared")
        self.ligand_count_label.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.CAPTION_M}px;
                color: {DesignTokens.Colors.TEXT_TERTIARY};
                padding: {DesignTokens.Spacing.PADDING_SM}px;
                background: {DesignTokens.Colors.GRAY_100};
                border-radius: {DesignTokens.BorderRadius.MD}px;
            }}
        """)
        self.ligand_count_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.ligand_count_label)
        
        select_ligand_btn = QPushButton("Select Ligand Files")
        select_ligand_btn.setStyleSheet(f"""
            QPushButton {{
                background: {DesignTokens.Colors.GRAY_100};
                color: {DesignTokens.Colors.TEXT_PRIMARY};
                border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                border-radius: {DesignTokens.BorderRadius.MD}px;
                padding: {DesignTokens.Spacing.PADDING_MD}px;
                font-weight: 600;
            }}
            QPushButton:hover {{
                background: {DesignTokens.Colors.GRAY_200};
            }}
        """)
        select_ligand_btn.clicked.connect(self.on_select_ligands)
        layout.addWidget(select_ligand_btn)
        
        layout.addStretch()
        
        return card
    
    def _create_engine_selection_card(self) -> QWidget:
        """Create engine selection card (Step 3) - Tri-Score Protocol"""
        card = QFrame()
        card.setStyleSheet(f"""
            QFrame {{
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                    stop:0 {DesignTokens.Colors.PRIMARY},
                    stop:1 {DesignTokens.Colors.SECONDARY});
                border-radius: {DesignTokens.BorderRadius.XL}px;
                padding: {DesignTokens.Spacing.PADDING_LG}px;
            }}
        """)
        
        layout = QVBoxLayout(card)
        layout.setSpacing(DesignTokens.Spacing.GAP_MD)
        
        header = QWidget()
        header_layout = QHBoxLayout(header)
        header_layout.setSpacing(DesignTokens.Spacing.GAP_MD)
        
        step_number = QLabel("3")
        step_number.setStyleSheet(f"""
            QLabel {{
                background: rgba(255, 255, 255, 0.2);
                color: white;
                border-radius: {DesignTokens.BorderRadius.LG}px;
                padding: 8px 14px;
                font-size: {DesignTokens.Typography.BODY_L}px;
                font-weight: 700;
                border: 1px solid rgba(255, 255, 255, 0.3);
            }}
        """)
        
        title_layout = QVBoxLayout()
        title_layout.setSpacing(4)
        
        title = QLabel("BioDockify Tri-Score Protocol")
        title.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.BODY_XL}px;
                font-weight: 700;
                color: white;
            }}
        """)
        
        subtitle = QLabel("International Standard")
        subtitle.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.CAPTION_M}px;
                color: rgba(255, 255, 255, 0.8);
                background: rgba(255, 255, 255, 0.2);
                padding: 4px 12px;
                border-radius: {DesignTokens.BorderRadius.FULL}px;
            }}
        """)
        
        title_layout.addWidget(title)
        title_layout.addWidget(subtitle)
        
        header_layout.addWidget(step_number)
        header_layout.addWidget(title_layout, 1)
        
        layout.addWidget(header)
        
        options_frame = QFrame()
        options_frame.setStyleSheet(f"""
            QFrame {{
                background: rgba(0, 0, 0, 0.2);
                border-radius: {DesignTokens.BorderRadius.LG}px;
                padding: {DesignTokens.Spacing.PADDING_MD}px;
                border: 1px solid rgba(255, 255, 255, 0.1);
            }}
        """)
        options_layout = QVBoxLayout(options_frame)
        options_layout.setSpacing(DesignTokens.Spacing.GAP_MD)
        
        engine_options = [
            ("Vina", "Physics-based scoring", True),
            ("GNINA", "Deep learning CNN", False),
            ("RF-Score", "Random Forest ML", False),
        ]
        
        self.engine_checkboxes = {}
        for engine, desc, checked in engine_options:
            option_widget = QWidget()
            option_layout = QHBoxLayout(option_widget)
            option_layout.setSpacing(DesignTokens.Spacing.GAP_MD)
            
            checkbox = QPushButton(engine)
            checkbox.setCheckable(True)
            checkbox.setChecked(checked)
            checkbox.setMinimumWidth(100)
            checkbox.setStyleSheet(f"""
                QPushButton {{
                    background: rgba(255, 255, 255, 0.1);
                    color: white;
                    border: 1px solid rgba(255, 255, 255, 0.3);
                    border-radius: {DesignTokens.BorderRadius.MD}px;
                    padding: 8px 16px;
                    font-weight: 600;
                }}
                QPushButton:checked {{
                    background: {DesignTokens.Colors.SUCCESS};
                    border-color: {DesignTokens.Colors.SUCCESS};
                }}
            """)
            self.engine_checkboxes[engine] = checkbox
            
            desc_label = QLabel(desc)
            desc_label.setStyleSheet(f"""
                QLabel {{
                    color: rgba(255, 255, 255, 0.8);
                    font-size: {DesignTokens.Typography.CAPTION_L}px;
                }}
            """)
            
            option_layout.addWidget(checkbox)
            option_layout.addWidget(desc_label, 1)
            
            options_layout.addWidget(option_widget)
        
        layout.addWidget(options_frame)
        
        info_label = QLabel("Powered by Vina (Physics) + Gnina (Deep Learning) + RF-Score (ML)")
        info_label.setStyleSheet(f"""
            QLabel {{
                color: rgba(255, 255, 255, 0.7);
                font-size: {DesignTokens.Typography.CAPTION_L}px;
                font-style: italic;
            }}
        """)
        layout.addWidget(info_label)
        
        return card
    
    def _create_grid_config_card(self) -> QWidget:
        """Create grid configuration card"""
        card = QFrame()
        card.setStyleSheet(f"""
            QFrame {{
                background: {DesignTokens.Colors.WHITE};
                border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                border-radius: {DesignTokens.BorderRadius.XL}px;
                padding: {DesignTokens.Spacing.PADDING_LG}px;
            }}
        """)
        
        layout = QVBoxLayout(card)
        layout.setSpacing(DesignTokens.Spacing.GAP_MD)
        
        header = QWidget()
        header_layout = QHBoxLayout(header)
        header_layout.setSpacing(DesignTokens.Spacing.GAP_MD)
        
        title = QLabel("Grid Box Configuration")
        title.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.BODY_XL}px;
                font-weight: 700;
                color: {DesignTokens.Colors.TEXT_PRIMARY};
            }}
        """)
        
        header_layout.addWidget(title)
        header_layout.addStretch()
        
        layout.addWidget(header)
        
        grid_layout = QGridLayout()
        grid_layout.setSpacing(DesignTokens.Spacing.GAP_MD)
        
        center_labels = ["Center X", "Center Y", "Center Z"]
        size_labels = ["Size X", "Size Y", "Size Z"]
        defaults = [0, 0, 0, 20, 20, 20]
        
        self.center_inputs = []
        self.size_inputs = []
        
        for i, (label, default) in enumerate(zip(center_labels + size_labels, defaults)):
            lbl = QLabel(label)
            lbl.setStyleSheet(f"""
                QLabel {{
                    font-size: {DesignTokens.Typography.CAPTION_L}px;
                    font-weight: 600;
                    color: {DesignTokens.Colors.TEXT_SECONDARY};
                }}
            """)
            
            from PyQt6.QtWidgets import QSpinBox, QDoubleSpinBox
            if i < 3:
                inp = QDoubleSpinBox()
                inp.setRange(-100, 100)
            else:
                inp = QSpinBox()
                inp.setRange(1, 100)
            inp.setValue(default)
            inp.setStyleSheet(f"""
                QSpinBox, QDoubleSpinBox {{
                    border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                    border-radius: {DesignTokens.BorderRadius.MD}px;
                    padding: {DesignTokens.Spacing.PADDING_SM}px {DesignTokens.Spacing.PADDING_MD}px;
                    background: {DesignTokens.Colors.GRAY_50};
                    font-weight: 600;
                }}
                QSpinBox:focus, QDoubleSpinBox:focus {{
                    border-color: {DesignTokens.Colors.PRIMARY};
                }}
            """)
            
            if i < 3:
                self.center_inputs.append(inp)
                grid_layout.addWidget(lbl, 0, i)
                grid_layout.addWidget(inp, 1, i)
            else:
                self.size_inputs.append(inp)
                grid_layout.addWidget(lbl, 2, i - 3)
                grid_layout.addWidget(inp, 3, i - 3)
        
        layout.addLayout(grid_layout)
        
        exhaustiveness_label = QLabel("Exhaustiveness")
        exhaustiveness_label.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.CAPTION_L}px;
                font-weight: 600;
                color: {DesignTokens.Colors.TEXT_SECONDARY};
            }}
        """)
        layout.addWidget(exhaustiveness_label)
        
        self.exhaustiveness_slider = QSpinBox()
        self.exhaustiveness_slider.setRange(1, 32)
        self.exhaustiveness_slider.setValue(8)
        self.exhaustiveness_slider.setStyleSheet(f"""
            QSpinBox {{
                border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                border-radius: {DesignTokens.BorderRadius.MD}px;
                padding: {DesignTokens.Spacing.PADDING_SM}px {DesignTokens.Spacing.PADDING_MD}px;
                background: {DesignTokens.Colors.GRAY_50};
                font-weight: 600;
                font-size: {DesignTokens.Typography.BODY_L}px;
            }}
        """)
        layout.addWidget(self.exhaustiveness_slider)
        
        batch_label = QLabel("Batch Size (ligands per run)")
        batch_label.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.CAPTION_L}px;
                font-weight: 600;
                color: {DesignTokens.Colors.TEXT_SECONDARY};
                margin-top: {DesignTokens.Spacing.PADDING_MD}px;
            }}
        """)
        layout.addWidget(batch_label)
        
        batch_info = QLabel("Process 1-10 ligands at a time for memory efficiency")
        batch_info.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.CAPTION_M}px;
                color: {DesignTokens.Colors.TEXT_TERTIARY};
                margin-bottom: {DesignTokens.Spacing.PADDING_SM}px;
            }}
        """)
        layout.addWidget(batch_info)
        
        self.batch_size_slider = QSpinBox()
        self.batch_size_slider.setRange(1, 10)
        self.batch_size_slider.setValue(5)
        self.batch_size_slider.setStyleSheet(f"""
            QSpinBox {{
                border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                border-radius: {DesignTokens.BorderRadius.MD}px;
                padding: {DesignTokens.Spacing.PADDING_SM}px {DesignTokens.Spacing.PADDING_MD}px;
                background: {DesignTokens.Colors.GRAY_50};
                font-weight: 600;
                font-size: {DesignTokens.Typography.BODY_L}px;
            }}
        """)
        layout.addWidget(self.batch_size_slider)
        
        batch_buttons_layout = QHBoxLayout()
        batch_buttons_layout.setSpacing(DesignTokens.Spacing.GAP_SM)
        
        for size in [1, 3, 5, 10]:
            btn = QPushButton(str(size))
            btn.setFixedWidth(50)
            btn.setStyleSheet(f"""
                QPushButton {{
                    background: {DesignTokens.Colors.GRAY_100};
                    border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                    border-radius: {DesignTokens.BorderRadius.MD}px;
                    padding: {DesignTokens.Spacing.PADDING_SM}px;
                    font-weight: 600;
                    color: {DesignTokens.Colors.TEXT_SECONDARY};
                }}
                QPushButton:hover {{
                    background: {DesignTokens.Colors.PRIMARY_10};
                    border-color: {DesignTokens.Colors.PRIMARY};
                    color: {DesignTokens.Colors.PRIMARY};
                }}
            """)
            btn.clicked.connect(lambda checked, s=size: self.batch_size_slider.setValue(s))
            batch_buttons_layout.addWidget(btn)
        
        batch_buttons_layout.addStretch()
        layout.addLayout(batch_buttons_layout)
        
        return card
    
    def _create_gpu_status_card(self) -> QWidget:
        """Create GPU status indicator card"""
        card = QFrame()
        card.setStyleSheet(f"""
            QFrame {{
                background: {DesignTokens.Colors.WHITE};
                border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                border-radius: {DesignTokens.BorderRadius.XL}px;
                padding: {DesignTokens.Spacing.PADDING_LG}px;
            }}
        """)
        
        layout = QVBoxLayout(card)
        layout.setSpacing(DesignTokens.Spacing.GAP_MD)
        
        title = QLabel("Hardware Status")
        title.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.BODY_XL}px;
                font-weight: 700;
                color: {DesignTokens.Colors.TEXT_PRIMARY};
                margin-bottom: {DesignTokens.Spacing.PADDING_SM}px;
            }}
        """)
        layout.addWidget(title)
        
        gpu_widget = QWidget()
        gpu_layout = QHBoxLayout(gpu_widget)
        gpu_layout.setSpacing(DesignTokens.Spacing.GAP_MD)
        
        gpu_icon = QLabel("🖥️")
        gpu_icon.setStyleSheet("font-size: 24px;")
        
        gpu_info_layout = QVBoxLayout()
        gpu_info_layout.setSpacing(2)
        
        self.gpu_name_label = QLabel("Detecting GPU...")
        self.gpu_name_label.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.BODY_L}px;
                font-weight: 700;
                color: {DesignTokens.Colors.TEXT_PRIMARY};
            }}
        """)
        
        self.gpu_status_label = QLabel("Checking...")
        self.gpu_status_label.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.CAPTION_M}px;
                color: {DesignTokens.Colors.TEXT_SECONDARY};
            }}
        """)
        
        gpu_info_layout.addWidget(self.gpu_name_label)
        gpu_info_layout.addWidget(self.gpu_status_label)
        
        gpu_layout.addWidget(gpu_icon)
        gpu_layout.addWidget(gpu_info_layout, 1)
        
        layout.addWidget(gpu_widget)
        
        compute_mode_widget = QWidget()
        compute_mode_layout = QHBoxLayout(compute_mode_widget)
        compute_mode_layout.setSpacing(DesignTokens.Spacing.GAP_MD)
        
        compute_mode_label = QLabel("Compute Mode:")
        compute_mode_label.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.CAPTION_M}px;
                color: {DesignTokens.Colors.TEXT_SECONDARY};
            }}
        """)
        
        self.compute_mode_combo = QComboBox()
        self.compute_mode_combo.addItems(["AUTO (Adaptive)", "GPU Only", "CPU Only"])
        self.compute_mode_combo.setCurrentIndex(0)
        self.compute_mode_combo.setToolTip("Select compute mode. AUTO lets the system choose best available.")
        self.compute_mode_combo.setStyleSheet(f"""
            QComboBox {{
                padding: 6px 12px;
                border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                border-radius: {DesignTokens.BorderRadius.MD}px;
                background: {DesignTokens.Colors.BACKGROUND_LIGHT};
                min-width: 150px;
            }}
            QComboBox::drop-down {{
                border: none;
            }}
        """)
        
        compute_mode_layout.addWidget(compute_mode_label)
        compute_mode_layout.addWidget(self.compute_mode_combo)
        compute_mode_layout.addStretch()
        
        layout.addWidget(compute_mode_widget)
        
        self._detect_gpu()
        
        return card
    
    def _detect_gpu(self):
        """Detect GPU and update status"""
        try:
            import subprocess
            result = subprocess.run(
                ['nvidia-smi', '--query-gpu=name,memory.total', '--format=csv,noheader'],
                capture_output=True, text=True, timeout=5
            )
            if result.returncode == 0:
                lines = result.stdout.strip().split('\n')
                if lines:
                    name, memory = lines[0].split(',')
                    self.gpu_name_label.setText(name.strip())
                    self.gpu_status_label.setText(f"VRAM: {memory.strip()} | Mode: GPU")
                    self.gpu_status_label.setStyleSheet(f"""
                        QLabel {{
                            font-size: {DesignTokens.Typography.CAPTION_M}px;
                            color: {DesignTokens.Colors.SUCCESS};
                            font-weight: 600;
                        }}
                    """)
                    return
        except Exception:
            pass
        
        self.gpu_name_label.setText("No GPU Detected")
        self.gpu_status_label.setText("Mode: CPU (slower)")
        self.gpu_status_label.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.CAPTION_M}px;
                color: {DesignTokens.Colors.WARNING};
                font-weight: 600;
            }}
        """)
    
    def on_select_receptor(self):
        """Handle receptor file selection"""
        from PyQt6.QtWidgets import QFileDialog
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Receptor File",
            "",
            "Protein Files (*.pdb *.pdbqt *.mol2 *.mmcif);;All Files (*)"
        )
        if file_path:
            import os
            self.receptor_file_label.setText(os.path.basename(file_path))
            logger.info(f"Receptor selected: {file_path}")
    
    def on_select_ligands(self):
        """Handle ligand file selection"""
        from PyQt6.QtWidgets import QFileDialog
        file_paths, _ = QFileDialog.getOpenFileNames(
            self,
            "Select Ligand Files",
            "",
            "Ligand Files (*.sdf *.mol2 *.pdbqt *.smi *.smiles);;All Files (*)"
        )
        if file_paths:
            self.ligand_count_label.setText(f"{len(file_paths)} files prepared")
            logger.info(f"Ligands selected: {len(file_paths)} files")
    
    def on_upload_mode_changed(self, mode: str):
        """Handle upload mode toggle"""
        if mode == 'files':
            self.files_toggle.setChecked(True)
            self.csv_toggle.setChecked(False)
        else:
            self.files_toggle.setChecked(False)
            self.csv_toggle.setChecked(True)
        logger.info(f"Upload mode changed to: {mode}")
    
    def on_start_docking(self):
        """Handle start docking button click"""
        batch_size = getattr(self, 'batch_size_slider', None)
        batch_size_value = batch_size.value() if batch_size else 5
        
        logger.info(f"Starting docking experiment with batch size: {batch_size_value}")
        self.start_docking_btn.setEnabled(False)
        self.start_docking_btn.setText(f"Running (Batch: {batch_size_value})...")
        
        exhaustiveness = getattr(self, 'exhaustiveness_slider', None)
        exhaustiveness_value = exhaustiveness.value() if exhaustiveness else 8
        
        grid_params = {}
        if hasattr(self, 'center_inputs') and hasattr(self, 'size_inputs'):
            grid_params = {
                'center_x': self.center_inputs[0].value(),
                'center_y': self.center_inputs[1].value(),
                'center_z': self.center_inputs[2].value(),
                'size_x': self.size_inputs[0].value(),
                'size_y': self.size_inputs[1].value(),
                'size_z': self.size_inputs[2].value(),
            }
        
        logger.info(f"Docking config: batch_size={batch_size_value}, exhaustiveness={exhaustiveness_value}, grid={grid_params}")
    
    def _setup_toolbar(self):
        """Create professional toolbar"""
        toolbar = QToolBar()
        toolbar.setMovable(False)
        toolbar.setFloatable(False)
        toolbar.setStyleSheet(f"""
            QToolBar {{
                background: {DesignTokens.Colors.WHITE};
                border: none;
                spacing: {DesignTokens.Spacing.GAP_SM}px;
                padding: {DesignTokens.Spacing.PADDING_SM}px {DesignTokens.Spacing.PADDING_MD}px;
            }}
        """)
        
        # Toolbar actions
        actions = [
            ("New Job", "➕", self.on_new_job),
            ("Open", "📂", self.on_open_recent),
            ("Save", "💾", self.on_save),
            (None, None, None),  # Separator
            ("Undo", "↩", self.on_undo),
            ("Redo", "↪", self.on_redo),
            (None, None, None),  # Separator
            ("Zoom In", "🔍+", self.on_zoom_in),
            ("Zoom Out", "🔍-", self.on_zoom_out),
            (None, None, None),  # Separator
            ("Settings", "⚙", self.on_preferences),
        ]
        
        for name, icon, handler in actions:
            if name is None:
                toolbar.addSeparator()
            else:
                action = QAction(icon + " " + name, self)
                action.triggered.connect(handler)
                toolbar.addAction(action)
        
        self.addToolBar(toolbar)
    
    def _create_menu_bar(self) -> None:
        """Create professional menu bar with international structure"""
        menubar = self.menuBar()
        
        # File Menu
        file_menu = menubar.addMenu("&File")
        
        new_action = file_menu.addAction("&New Job")
        new_action.setShortcut("Ctrl+N")
        new_action.triggered.connect(self.on_new_job)
        
        file_menu.addSeparator()
        
        open_action = file_menu.addAction("&Open Recent...")
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(self.on_open_recent)
        
        file_menu.addSeparator()
        
        import_action = file_menu.addAction("&Import Molecules...")
        import_action.setShortcut("Ctrl+I")
        import_action.triggered.connect(self.on_import)
        
        file_menu.addSeparator()
        
        save_action = file_menu.addAction("&Save")
        save_action.setShortcut("Ctrl+S")
        save_action.triggered.connect(self.on_save)
        
        save_as_action = file_menu.addAction("Save &As...")
        save_as_action.setShortcut("Ctrl+Shift+S")
        save_as_action.triggered.connect(self.on_save_as)
        
        file_menu.addSeparator()
        
        export_action = file_menu.addAction("&Export Results...")
        export_action.setShortcut("Ctrl+E")
        export_action.triggered.connect(self.on_export_results)
        
        file_menu.addSeparator()
        
        exit_action = file_menu.addAction("E&xit")
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        
        # Edit Menu
        edit_menu = menubar.addMenu("&Edit")
        
        undo_action = edit_menu.addAction("&Undo")
        undo_action.setShortcut("Ctrl+Z")
        undo_action.triggered.connect(self.on_undo)
        
        redo_action = edit_menu.addAction("&Redo")
        redo_action.setShortcut("Ctrl+Y")
        redo_action.triggered.connect(self.on_redo)
        
        edit_menu.addSeparator()
        
        cut_action = edit_menu.addAction("Cu&t")
        cut_action.setShortcut("Ctrl+X")
        cut_action.triggered.connect(self.on_cut)
        
        copy_action = edit_menu.addAction("&Copy")
        copy_action.setShortcut("Ctrl+C")
        copy_action.triggered.connect(self.on_copy)
        
        paste_action = edit_menu.addAction("&Paste")
        paste_action.setShortcut("Ctrl+V")
        paste_action.triggered.connect(self.on_paste)
        
        edit_menu.addSeparator()
        
        preferences_action = edit_menu.addAction("Pre&ferences...")
        preferences_action.setShortcut("Ctrl+,")
        preferences_action.triggered.connect(self.on_preferences)
        
        # View Menu
        view_menu = menubar.addMenu("&View")
        
        fullscreen_action = view_menu.addAction("&Full Screen")
        fullscreen_action.setShortcut("F11")
        fullscreen_action.setCheckable(True)
        fullscreen_action.triggered.connect(self.on_fullscreen)
        
        view_menu.addSeparator()
        
        reset_layout_action = view_menu.addAction("Reset &Layout")
        reset_layout_action.setShortcut("Ctrl+Shift+R")
        reset_layout_action.triggered.connect(self.on_reset_layout)
        
        view_menu.addSeparator()
        
        chat_action = view_menu.addAction("&AI Assistant")
        chat_action.setShortcut("Ctrl+Shift+A")
        chat_action.setCheckable(True)
        chat_action.triggered.connect(self.toggle_chat_dock)
        
        job_monitor_action = view_menu.addAction("&Job Monitor")
        job_monitor_action.setShortcut("Ctrl+Shift+J")
        job_monitor_action.setCheckable(True)
        job_monitor_action.triggered.connect(self.toggle_job_monitor)
        
        # Tools Menu
        tools_menu = menubar.addMenu("&Tools")
        
        docker_action = tools_menu.addAction("&Docker Desktop")
        docker_action.triggered.connect(self.on_docker)
        
        tools_menu.addSeparator()
        
        validator_action = tools_menu.addAction("Molecule &Validator...")
        validator_action.triggered.connect(self.on_validator)
        
        converter_action = tools_menu.addAction("File &Converter...")
        converter_action.triggered.connect(self.on_converter)
        
        # Help Menu
        help_menu = menubar.addMenu("&Help")
        
        documentation_action = help_menu.addAction("&Documentation")
        documentation_action.setShortcut("F1")
        documentation_action.triggered.connect(self.on_documentation)
        
        tutorials_action = help_menu.addAction("&Tutorials")
        tutorials_action.triggered.connect(self.on_tutorials)
        
        help_menu.addSeparator()
        
        check_updates_action = help_menu.addAction("Check for &Updates...")
        check_updates_action.triggered.connect(self.on_check_updates)
        
        help_menu.addSeparator()
        
        about_action = help_menu.addAction("&About BioDockify...")
        about_action.triggered.connect(self.on_about)
        
        # Apply professional menu styling
        menubar.setStyleSheet(f"""
            QMenuBar {{
                background: {DesignTokens.Colors.WHITE};
                border-bottom: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                padding: 0 {DesignTokens.Spacing.PADDING_MD}px;
            }}
            
            QMenuBar::item {{
                background: transparent;
                padding: {DesignTokens.Spacing.PADDING_SM}px {DesignTokens.Spacing.PADDING_MD}px;
                font-size: {DesignTokens.Typography.BODY_L}px;
                color: {DesignTokens.Colors.TEXT_PRIMARY};
                border-radius: {DesignTokens.BorderRadius.SM}px;
                spacing: 2px;
            }}
            
            QMenuBar::item:selected {{
                background: {DesignTokens.Colors.PRIMARY_10};
                color: {DesignTokens.Colors.PRIMARY};
            }}
            
            QMenuBar::item:pressed {{
                background: {DesignTokens.Colors.PRIMARY};
                color: {DesignTokens.Colors.WHITE};
            }}
        """)
    
    def _create_status_bar(self) -> None:
        """Create professional status bar"""
        status_bar = QStatusBar()
        status_bar.setStyleSheet(f"""
            QStatusBar {{
                background: {DesignTokens.Colors.BACKGROUND_SIDEBAR};
                color: {DesignTokens.Colors.TEXT_INVERTED};
                border: none;
                padding: {DesignTokens.Spacing.PADDING_SM}px {DesignTokens.Spacing.PADDING_MD}px;
                font-size: {DesignTokens.Typography.CAPTION_L}px;
            }}
        """)
        
        self.setStatusBar(status_bar)
        
        # Status sections
        # Job status
        self.job_status_label = QLabel("Ready")
        job_status_widget = self._create_status_badge("Job Status", "Ready")
        status_bar.addPermanentWidget(job_status_widget)
        
        # Docker status
        docker_status_widget = self._create_status_badge("Docker", "✓ Running")
        status_bar.addPermanentWidget(docker_status_widget)
        
        # Version info
        version_label = QLabel("v1.2.2")
        version_label.setStyleSheet(f"""
            QLabel {{
                color: {DesignTokens.Colors.TEXT_INVERTED};
                font-weight: 500;
                opacity: 0.8;
            }}
        """)
        status_bar.addPermanentWidget(version_label)
    
    def _setup_chat_dock(self):
        """Setup dockable chat widget"""
        self.chat_widget = ChatWidget()
        
        # Initialize chat service
        try:
            from src.config_settings import Settings
            settings = Settings()
            app_settings = settings.get()
            chat_config = ChatConfig(
                provider=ModelProvider(app_settings.chat_provider.provider),
                model=app_settings.chat_provider.model,
                api_base=app_settings.chat_provider.api_base,
                api_key=app_settings.chat_provider.api_key,
                temperature=app_settings.chat_provider.temperature,
                max_tokens=app_settings.chat_provider.max_tokens
            )
            self.chat_service = ChatService(config=chat_config)
            self.chat_widget.chat_service = self.chat_service
        except Exception as e:
            logger.warning(f"Could not initialize chat service: {e}")
            self.chat_service = None
        
        # Create dock widget
        self.chat_dock = QDockWidget("AI Assistant", self)
        self.chat_dock.setWidget(self.chat_widget)
        self.chat_dock.setAllowedAreas(Qt.DockWidgetArea.RightDockWidgetArea | Qt.DockWidgetArea.LeftDockWidgetArea)
        self.chat_dock.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetClosable | QDockWidget.DockWidgetFeature.DockWidgetMovable)
        self.chat_dock.setMinimumWidth(350)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.chat_dock)
        
        # Hide by default
        self.chat_dock.hide()
        
        # Setup job monitor dock
        self._setup_job_monitor_dock()
        
        logger.info("Chat dock widget initialized")
    
    def _setup_job_monitor_dock(self):
        """Setup dockable job monitor widget"""
        from PyQt6.QtWidgets import QTableWidget, QTableWidgetItem, QHeaderView
        
        job_monitor = QWidget()
        job_layout = QVBoxLayout(job_monitor)
        job_layout.setContentsMargins(4, 4, 4, 4)
        
        # Header
        header_label = QLabel("Active Jobs")
        header_label.setStyleSheet(f"font-weight: bold; color: {DesignTokens.Colors.TEXT_PRIMARY}; padding: 4px;")
        job_layout.addWidget(header_label)
        
        # Jobs table
        self.jobs_table = QTableWidget(0, 4)
        self.jobs_table.setHorizontalHeaderLabels(["Job ID", "Status", "Progress", "Engine"])
        self.jobs_table.setStyleSheet("""
            QTableWidget {
                border: none;
                background: transparent;
            }
            QTableWidget::item {
                padding: 8px;
            }
        """)
        self.jobs_table.horizontalHeader().setStretchLastSection(True)
        self.jobs_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.jobs_table.setShowGrid(False)
        job_layout.addWidget(self.jobs_table, 1)
        
        # Create dock widget
        self.job_monitor_dock = QDockWidget("Job Monitor", self)
        self.job_monitor_dock.setWidget(job_monitor)
        self.job_monitor_dock.setAllowedAreas(Qt.DockWidgetArea.RightDockWidgetArea | Qt.DockWidgetArea.BottomDockWidgetArea)
        self.job_monitor_dock.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetClosable | QDockWidget.DockWidgetFeature.DockWidgetMovable)
        self.job_monitor_dock.setMinimumHeight(200)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.job_monitor_dock)
        
        # Hide by default
        self.job_monitor_dock.hide()
        
        # Connect signals
        self.job_started.connect(self._on_job_started)
        self.job_completed.connect(self._on_job_completed)
        self.job_failed.connect(self._on_job_failed)
        self.job_progress.connect(self._on_job_progress)
        
        logger.info("Job monitor dock widget initialized")
    
    def _on_job_started(self, job_id: str):
        """Handle job started"""
        row = self.jobs_table.rowCount()
        self.jobs_table.insertRow(row)
        self.jobs_table.setItem(row, 0, QTableWidgetItem(job_id))
        self.jobs_table.setItem(row, 1, QTableWidgetItem("Running"))
        self.jobs_table.setItem(row, 2, QTableWidgetItem("0%"))
        self.jobs_table.setItem(row, 3, QTableWidgetItem(""))
        logger.info(f"Job {job_id} started - added to monitor")
    
    def _on_job_completed(self, job_id: str):
        """Handle job completed"""
        for row in range(self.jobs_table.rowCount()):
            item = self.jobs_table.item(row, 0)
            if item and item.text() == job_id:
                self.jobs_table.setItem(row, 1, QTableWidgetItem("✓ Completed"))
                self.jobs_table.setItem(row, 2, QTableWidgetItem("100%"))
                break
        logger.info(f"Job {job_id} completed")
    
    def _on_job_failed(self, job_id: str, error: str):
        """Handle job failed"""
        for row in range(self.jobs_table.rowCount()):
            item = self.jobs_table.item(row, 0)
            if item and item.text() == job_id:
                self.jobs_table.setItem(row, 1, QTableWidgetItem("✗ Failed"))
                self.jobs_table.item(row, 1).setBackground(QColor(DesignTokens.Colors.ERROR))
                break
        logger.error(f"Job {job_id} failed: {error}")
    
    def _on_job_progress(self, job_id: str, percentage: int):
        """Handle job progress update"""
        for row in range(self.jobs_table.rowCount()):
            item = self.jobs_table.item(row, 0)
            if item and item.text() == job_id:
                self.jobs_table.setItem(row, 2, QTableWidgetItem(f"{percentage}%"))
                break
    
    def _create_status_badge(self, label: str, status: str) -> QWidget:
        """Create professional status badge for status bar"""
        badge = QWidget()
        badge_layout = QHBoxLayout(badge)
        badge_layout.setSpacing(DesignTokens.Spacing.GAP_XS)
        
        label_widget = QLabel(label)
        label_widget.setStyleSheet(f"""
            QLabel {{
                color: rgba(255, 255, 255, 0.7);
                font-size: {DesignTokens.Typography.CAPTION_L}px;
            }}
        """)
        
        status_widget = QLabel(status)
        status_widget.setStyleSheet(f"""
            QLabel {{
                background: rgba(255, 255, 255, 0.15);
                padding: 3px 8px;
                border-radius: {DesignTokens.BorderRadius.SM}px;
                font-size: {DesignTokens.Typography.CAPTION_L}px;
                font-weight: 500;
            }}
        """)
        
        badge_layout.addWidget(label_widget)
        badge_layout.addWidget(status_widget)
        
        return badge
    
    # Menu Actions
    def on_new_job(self):
        """Handle new job button click"""
        logger.info("New job button clicked")
        self.job_started.emit("new_job")
    
    def on_open_recent(self):
        """Handle open recent menu action"""
        logger.info("Open recent menu action clicked")
    
    def on_save(self):
        """Handle save menu action"""
        logger.info("Save menu action clicked")
    
    def on_save_as(self):
        """Handle save as menu action"""
        logger.info("Save as menu action clicked")
    
    def on_import(self):
        """Handle import menu action"""
        logger.info("Import menu action clicked")
    
    def on_export_results(self):
        """Handle export results menu action"""
        logger.info("Export results menu action clicked")
    
    def on_undo(self):
        """Handle undo menu action"""
        logger.info("Undo menu action clicked")
    
    def on_redo(self):
        """Handle redo menu action"""
        logger.info("Redo menu action clicked")
    
    def on_cut(self):
        """Handle cut menu action"""
        logger.info("Cut menu action clicked")
    
    def on_copy(self):
        """Handle copy menu action"""
        logger.info("Copy menu action clicked")
    
    def on_paste(self):
        """Handle paste menu action"""
        logger.info("Paste menu action clicked")
    
    def on_preferences(self):
        """Handle preferences menu action"""
        logger.info("Preferences menu action clicked")
        from src.ui.settings_dialog import show_settings_dialog
        show_settings_dialog(self)
    
    def on_fullscreen(self, checked):
        """Handle fullscreen toggle"""
        if checked:
            self.showFullScreen()
        else:
            self.showNormal()
    
    def on_reset_layout(self):
        """Handle reset layout action"""
        logger.info("Reset layout menu action clicked")
    
    def toggle_chat_dock(self, checked):
        """Toggle AI Assistant dock visibility"""
        if hasattr(self, 'chat_dock'):
            self.chat_dock.setVisible(checked)
            logger.info(f"Chat dock {'shown' if checked else 'hidden'}")
    
    def toggle_job_monitor(self, checked):
        """Toggle Job Monitor dock visibility"""
        if hasattr(self, 'job_monitor_dock'):
            self.job_monitor_dock.setVisible(checked)
            logger.info(f"Job monitor {'shown' if checked else 'hidden'}")
        else:
            logger.info("Job monitor not yet implemented")
    
    def on_docker(self):
        """Handle Docker menu action"""
        logger.info("Docker menu action clicked")
    
    def on_validator(self):
        """Handle validator menu action"""
        logger.info("Validator menu action clicked")
    
    def on_converter(self):
        """Handle converter menu action"""
        logger.info("Converter menu action clicked")
    
    def on_documentation(self):
        """Handle documentation menu action"""
        logger.info("Documentation menu action clicked")
        # In real implementation, would open documentation in browser
    
    def on_tutorials(self):
        """Handle tutorials menu action"""
        logger.info("Tutorials menu action clicked")
    
    def on_check_updates(self):
        """Handle check for updates menu action"""
        logger.info("Check for updates menu action clicked")
        # In real implementation, would check GitHub for new releases
    
    def on_about(self):
        """Handle about menu action"""
        from PyQt6.QtWidgets import QMessageBox
        QMessageBox.about(
            self,
            "BioDockify Docking Studio",
            """
            <h2 style='margin: 0 0 10px 0; color: #2E5AAC;'>BioDockify Docking Studio</h2>
            <p style='margin: 0; color: #6B7280;'>Version 1.2.2</p>
            <hr style='border: none; border-top: 1px solid #E5E7EB; margin: 20px 0;'>
            <p style='margin: 0; color: #4B5563;'>
            Molecular Docking with Intelligent Self-Repair
            </p>
            <p style='margin: 20px 0 10px 0; color: #9CA3AF; font-size: 12px;'>
            Copyright © 2025 BioDockify Development Team<br>
            License: Apache License 2.0<br><br>
            <strong>Built with:</strong><br>
            • AutoDock Vina<br>
            • ODDT (Open Drug Discovery Toolkit)<br>
            • RDKit (Chemistry Development Kit)<br>
            • Docker-based execution environment<br><br>
            <strong>Agent Zero™ AI System</strong> for intelligent 
            failure detection and recovery
            </p>
            """
        )
    
    def on_zoom_in(self):
        """Handle zoom in"""
        logger.info("Zoom in")
    
    def on_zoom_out(self):
        """Handle zoom out"""
        logger.info("Zoom out")
    
    def on_cancel_job(self):
        """Handle cancel job button click"""
        logger.info("Cancel job button clicked")
        self.job_cancelled.emit("current_job")
    
    def update_progress(self, percentage: int, message: str = "") -> None:
        """Update progress bar with smooth animation"""
        self.job_progress.emit("current_job", percentage)
        logger.debug(f"Progress updated: {percentage}% - {message}")
