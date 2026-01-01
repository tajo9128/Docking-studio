"""
BioDockify Docking Studio - Professional Main Window
International-quality UI implementation
"""

from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QStackedWidget, QFrame,
    QProgressBar, QScrollArea, QSizePolicy, QDockWidget,
    QToolBar, QAction, QStatusBar, QMenu, QMenuBar,
    QSplitter, QTabWidget
)
from PyQt6.QtCore import Qt, QSize, QTimer, pyqtSignal, QPoint, QPropertyAnimation, QEasingCurve
from PyQt6.QtGui import QIcon, QFont, QAction, QPalette, QColor, QPixmap, QPainter, QLinearGradient
from PyQt6.QtSvgWidgets import QSvgWidget
import logging

from src.ui.theme import DesignTokens, Stylesheet

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
        
        version_label = QLabel("Version 1.0.73")
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
        """Create main content area"""
        content = QWidget()
        
        content_layout = QVBoxLayout(content)
        content_layout.setSpacing(DesignTokens.Spacing.GAP_MD)
        
        # Breadcrumb
        breadcrumb = QLabel("Dashboard → Jobs")
        breadcrumb.setProperty("caption", "true")
        content_layout.addWidget(breadcrumb)
        
        # Page title
        page_title = QLabel("Active Jobs")
        page_title.setProperty("heading", "true")
        page_title.setStyleSheet(f"""
            QLabel[heading="true"] {{
                font-size: {DesignTokens.Typography.DISPLAY_L}px;
                font-weight: 700;
                color: {DesignTokens.Colors.TEXT_PRIMARY};
            }}
        """)
        content_layout.addWidget(page_title)
        
        # Content cards grid
        cards_widget = self._create_job_cards()
        content_layout.addWidget(cards_widget, 1)
        
        return content
    
    def _create_job_cards(self) -> QWidget:
        """Create professional job cards grid"""
        cards_widget = QWidget()
        cards_layout = QVBoxLayout(cards_widget)
        
        # Stats row
        stats_widget = self._create_stats_row()
        cards_layout.addWidget(stats_widget)
        
        # Recent jobs
        recent_label = QLabel("Recent Jobs")
        recent_label.setProperty("heading", "true")
        cards_layout.addWidget(recent_label)
        
        # Scrollable jobs list
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setStyleSheet(f"""
            QScrollArea {{
                border: none;
                background: {DesignTokens.Colors.BACKGROUND_TERTIARY};
                border-radius: {DesignTokens.BorderRadius.LG}px;
            }}
        """)
        
        jobs_container = QWidget()
        jobs_layout = QVBoxLayout(jobs_container)
        jobs_layout.setSpacing(DesignTokens.Spacing.GAP_MD)
        jobs_layout.setContentsMargins(DesignTokens.Spacing.PADDING_MD,
                                       DesignTokens.Spacing.PADDING_MD,
                                       DesignTokens.Spacing.PADDING_MD,
                                       DesignTokens.Spacing.PADDING_MD)
        
        # Add sample job cards
        for i in range(5):
            job_card = self._create_job_card(f"Job-{i+1}")
            jobs_layout.addWidget(job_card)
        
        jobs_layout.addStretch()
        scroll.setWidget(jobs_container)
        cards_layout.addWidget(scroll, 1)
        
        return cards_widget
    
    def _create_job_card(self, job_id: str) -> QFrame:
        """Create professional job card"""
        card = QFrame()
        card.setProperty("card", "true")
        card.setMinimumHeight(120)
        card.setCursor(Qt.CursorShape.PointingHandCursor)
        
        card_layout = QVBoxLayout(card)
        card_layout.setSpacing(DesignTokens.Spacing.GAP_SM)
        
        # Card header
        header_layout = QHBoxLayout()
        job_title = QLabel(f"{job_id}: Protein Docking")
        job_title.setProperty("heading", "true")
        
        status_badge = QLabel("Running")
        status_badge.setStyleSheet(f"""
            QLabel {{
                background: {DesignTokens.Colors.INFO_BG};
                color: {DesignTokens.Colors.INFO};
                padding: 4px 12px;
                border-radius: {DesignTokens.BorderRadius.FULL}px;
                font-size: {DesignTokens.Typography.CAPTION_M}px;
                font-weight: 600;
            }}
        """)
        
        header_layout.addWidget(job_title)
        header_layout.addStretch()
        header_layout.addWidget(status_badge)
        
        # Card content
        content = QLabel(\"\"\"
            <div style="color: #6B7280; font-size: 13px;">
            <strong>Receptor:</strong> 3CL7.pdbqt &nbsp;&nbsp;
            <strong>Ligand:</strong> N3.pdbqt<br>
            <strong>Progress:</strong> 67% &nbsp;&nbsp;
            <strong>ETA:</strong> 12m 34s
            </div>
        \"\"\")
        content.setTextFormat(Qt.TextFormat.RichText)
        content.setProperty("caption", "true")
        
        # Progress bar
        progress = QProgressBar()
        progress.setValue(67)
        progress.setStyleSheet(f"""
            QProgressBar {{
                background: {DesignTokens.Colors.GRAY_200};
                border: none;
                border-radius: {DesignTokens.BorderRadius.FULL}px;
                height: 6px;
            }}
            QProgressBar::chunk {{
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                        stop:0 {DesignTokens.Colors.PRIMARY},
                        stop:1 {DesignTokens.Colors.SECONDARY});
                border-radius: {DesignTokens.BorderRadius.FULL}px;
            }}
        """)
        
        card_layout.addLayout(header_layout)
        card_layout.addWidget(content)
        card_layout.addWidget(progress)
        
        return card
    
    def _create_stats_row(self) -> QWidget:
        """Create professional statistics row"""
        stats_widget = QWidget()
        stats_layout = QHBoxLayout(stats_widget)
        stats_layout.setSpacing(DesignTokens.Spacing.GAP_MD)
        
        stats = [
            ("Total Jobs", "127", DesignTokens.Colors.PRIMARY),
            ("Completed", "98", DesignTokens.Colors.SUCCESS),
            ("Running", "5", DesignTokens.Colors.INFO),
            ("Failed", "24", DesignTokens.Colors.ERROR),
        ]
        
        for title, value, color in stats:
            stat_card = QWidget()
            stat_card.setProperty("card", "true")
            stat_layout_card = QVBoxLayout(stat_card)
            stat_layout_card.setSpacing(DesignTokens.Spacing.GAP_XS)
            
            stat_title = QLabel(title)
            stat_title.setProperty("caption", "true")
            
            stat_value = QLabel(value)
            stat_value.setStyleSheet(f"""
                QLabel {{
                    font-size: {DesignTokens.Typography.DISPLAY_XL}px;
                    font-weight: 700;
                    color: {color};
                }}
            """)
            
            stat_layout_card.addWidget(stat_title)
            stat_layout_card.addWidget(stat_value)
            
            stats_layout.addWidget(stat_card)
        
        return stats_widget
    
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
        version_label = QLabel("v1.0.73")
        version_label.setStyleSheet(f"""
            QLabel {{
                color: {DesignTokens.Colors.TEXT_INVERTED};
                font-weight: 500;
                opacity: 0.8;
            }}
        """)
        status_bar.addPermanentWidget(version_label)
    
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
        # In real implementation, would open preferences dialog
    
    def on_fullscreen(self, checked):
        """Handle fullscreen toggle"""
        if checked:
            self.showFullScreen()
        else:
            self.showNormal()
    
    def on_reset_layout(self):
        """Handle reset layout action"""
        logger.info("Reset layout menu action clicked")
    
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
            <p style='margin: 0; color: #6B7280;'>Version 1.0.73</p>
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
