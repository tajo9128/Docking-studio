"""
BioDockify Docking Studio - Professional Main Window
International-quality UI implementation
"""

from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QStackedWidget, QFrame,
    QProgressBar, QScrollArea, QSizePolicy, QDockWidget,
    QToolBar, QStatusBar, QMenu, QMenuBar,
    QSplitter, QTabWidget, QSpinBox, QDoubleSpinBox, QGridLayout,
    QComboBox
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
        
        # Check backend connection
        self._check_backend_connection()
        
        logger.info("Professional main window initialized")
    
    def _check_backend_connection(self):
        """Check backend connection and update status"""
        try:
            from src.backend_connection import get_connection, BackendStatus
            
            conn = get_connection()
            info = conn.check_health()
            
            if hasattr(self, 'log_viewer'):
                if info.status == BackendStatus.CONNECTED:
                    self.log_viewer.append_log('INFO', f"Backend connected: v{info.version}")
                    self.log_viewer.append_log('INFO', f"AI Provider: {info.ollama_provider}")
                else:
                    self.log_viewer.append_log('WARNING', "Backend not connected")
                    self.log_viewer.append_log('INFO', "Run: docker compose up -d")
            
            if hasattr(self, 'status_label'):
                if info.status == BackendStatus.CONNECTED:
                    self.status_label.setText(f"Backend: Connected (v{info.version})")
                else:
                    self.status_label.setText("Backend: Disconnected")
                    
        except Exception as e:
            logger.warning(f"Backend check failed: {e}")
    
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
        """Create professional header with quick actions"""
        header = QWidget()
        header.setStyleSheet(f"""
            QWidget {{
                background: {DesignTokens.Colors.WHITE};
                border-bottom: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
            }}
        """)
        
        header_layout = QHBoxLayout(header)
        header_layout.setSpacing(16)
        header_layout.setContentsMargins(24, 16, 24, 16)
        
        # Left: Breadcrumb / Page title
        breadcrumb = QWidget()
        breadcrumb_layout = QVBoxLayout(breadcrumb)
        breadcrumb_layout.setSpacing(2)
        breadcrumb_layout.setContentsMargins(0, 0, 0, 0)
        
        page_title = QLabel("Molecular Docking")
        page_title.setStyleSheet(f"""
            QLabel {{
                font-size: 20px;
                font-weight: 700;
                color: {DesignTokens.Colors.TEXT_PRIMARY};
            }}
        """)
        
        subtitle = QLabel("Virtual screening and binding affinity prediction")
        subtitle.setStyleSheet(f"""
            QLabel {{
                font-size: 13px;
                color: {DesignTokens.Colors.TEXT_SECONDARY};
            }}
        """)
        
        breadcrumb_layout.addWidget(page_title)
        breadcrumb_layout.addWidget(subtitle)
        
        header_layout.addWidget(breadcrumb, 1)
        
        # Right: Quick action buttons
        actions_layout = QHBoxLayout()
        actions_layout.setSpacing(8)
        
        # Quick action buttons with icons
        quick_actions = [
            ("Import", "import", False),
            ("Export", "export", False),
            ("Help", "help", False),
        ]
        
        for label, icon_name, primary in quick_actions:
            btn = self._create_action_button(label, icon_name, primary)
            actions_layout.addWidget(btn)
        
        header_layout.addLayout(actions_layout)
        
        return header
    
    def _create_action_button(self, label: str, icon_name: str, primary: bool = False) -> QPushButton:
        """Create header action button"""
        icon_map = {
            "import": "📥",
            "export": "📤",
            "help": "❓",
            "settings": "⚙",
            "refresh": "🔄",
        }
        icon = icon_map.get(icon_name, "•")
        
        if primary:
            btn = QPushButton(f"{icon} {label}")
            btn.setStyleSheet(f"""
                QPushButton {{
                    background: {DesignTokens.Colors.PRIMARY};
                    color: {DesignTokens.Colors.WHITE};
                    border: none;
                    border-radius: 8px;
                    padding: 10px 20px;
                    font-size: 13px;
                    font-weight: 600;
                }}
                QPushButton:hover {{
                    background: {DesignTokens.Colors.PRIMARY_HOVER};
                }}
            """)
        else:
            btn = QPushButton(f"{icon} {label}")
            btn.setStyleSheet(f"""
                QPushButton {{
                    background: {DesignTokens.Colors.GRAY_100};
                    color: {DesignTokens.Colors.TEXT_PRIMARY};
                    border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                    border-radius: 8px;
                    padding: 10px 16px;
                    font-size: 13px;
                    font-weight: 500;
                }}
                QPushButton:hover {{
                    background: {DesignTokens.Colors.GRAY_200};
                    border-color: {DesignTokens.Colors.GRAY_400};
                }}
            """)
        
        btn.setCursor(Qt.CursorShape.PointingHandCursor)
        return btn
    
    def _create_status_indicator(self) -> QWidget:
        """Create professional status indicator with security status"""
        status_widget = QWidget()
        status_widget.setFixedWidth(220)
        
        status_layout = QVBoxLayout(status_widget)
        status_layout.setSpacing(DesignTokens.Spacing.GAP_XS)
        
        system_label = QLabel("System Status")
        system_label.setProperty("caption", "true")
        
        self.system_status_badge = QLabel("Ready")
        self.system_status_badge.setStyleSheet(f"""
            QLabel {{
                background: {DesignTokens.Colors.SUCCESS_BG};
                color: {DesignTokens.Colors.SUCCESS};
                padding: 6px 12px;
                border-radius: {DesignTokens.BorderRadius.FULL}px;
                font-size: {DesignTokens.Typography.CAPTION_L}px;
                font-weight: 600;
                text-align: center;
                border: 1px solid {DesignTokens.Colors.SUCCESS};
            }}
        """)
        
        security_label = QLabel("Security")
        security_label.setProperty("caption", "true")
        
        self.security_status_badge = QLabel("Checking...")
        self.security_status_badge.setStyleSheet(f"""
            QLabel {{
                background: {DesignTokens.Colors.GRAY_100};
                color: {DesignTokens.Colors.TEXT_SECONDARY};
                padding: 6px 12px;
                border-radius: {DesignTokens.BorderRadius.FULL}px;
                font-size: {DesignTokens.Typography.CAPTION_L}px;
                font-weight: 600;
                text-align: center;
            }}
        """)
        
        status_layout.addWidget(system_label, 0, Qt.AlignmentFlag.AlignRight)
        status_layout.addWidget(self.system_status_badge, 0, Qt.AlignmentFlag.AlignRight)
        status_layout.addWidget(security_label, 0, Qt.AlignmentFlag.AlignRight)
        status_layout.addWidget(self.security_status_badge, 0, Qt.AlignmentFlag.AlignRight)
        
        self._check_security_status()
        
        self.security_timer = QTimer(self)
        self.security_timer.timeout.connect(self._check_security_status)
        self.security_timer.start(30000)
        
        return status_widget
    
    def _check_security_status(self):
        """Check security status from backend"""
        try:
            from src.backend_connection import get_connection
            conn = get_connection()
            sec_data = conn.get("/security/status")
            
            if sec_data:
                is_secure = sec_data.get('is_secure', True)
                severity = sec_data.get('overall_severity', 'UNKNOWN')
                total_issues = sec_data.get('total_issues', 0)
                
                if is_secure and total_issues == 0:
                    self.security_status_badge.setText("Secure")
                    self.security_status_badge.setStyleSheet(f"""
                        QLabel {{
                            background: {DesignTokens.Colors.SUCCESS_BG};
                            color: {DesignTokens.Colors.SUCCESS};
                            padding: 6px 12px;
                            border-radius: {DesignTokens.BorderRadius.FULL}px;
                            font-size: {DesignTokens.Typography.CAPTION_L}px;
                            font-weight: 600;
                            text-align: center;
                            border: 1px solid {DesignTokens.Colors.SUCCESS};
                        }}
                    """)
                elif severity in ['HIGH', 'CRITICAL']:
                    self.security_status_badge.setText(f"{severity}: {total_issues}")
                    self.security_status_badge.setStyleSheet(f"""
                        QLabel {{
                            background: #ffebee;
                            color: #c62828;
                            padding: 6px 12px;
                            border-radius: {DesignTokens.BorderRadius.FULL}px;
                            font-size: {DesignTokens.Typography.CAPTION_L}px;
                            font-weight: 600;
                            text-align: center;
                            border: 1px solid #c62828;
                        }}
                    """)
                else:
                    self.security_status_badge.setText(f"Issues: {total_issues}")
                    self.security_status_badge.setStyleSheet(f"""
                        QLabel {{
                            background: #fff3e0;
                            color: #e65100;
                            padding: 6px 12px;
                            border-radius: {DesignTokens.BorderRadius.FULL}px;
                            font-size: {DesignTokens.Typography.CAPTION_L}px;
                            font-weight: 600;
                            text-align: center;
                            border: 1px solid #e65100;
                        }}
                    """)
        except Exception as e:
            logger.debug(f"Security status check failed: {e}")
            self.security_status_badge.setText("N/A")
            self.security_status_badge.setStyleSheet(f"""
                QLabel {{
                    background: {DesignTokens.Colors.GRAY_100};
                    color: {DesignTokens.Colors.TEXT_TERTIARY};
                    padding: 6px 12px;
                    border-radius: {DesignTokens.BorderRadius.FULL}px;
                    font-size: {DesignTokens.Typography.CAPTION_L}px;
                    font-weight: 600;
                    text-align: center;
                }}
            """)
    
    def _create_sidebar(self) -> QWidget:
        """Create professional sidebar navigation with organized sections"""
        sidebar = QWidget()
        sidebar.setMinimumWidth(240)
        sidebar.setMaximumWidth(240)
        sidebar.setStyleSheet(f"""
            QWidget {{
                background: {DesignTokens.Colors.BACKGROUND_SIDEBAR};
            }}
        """)
        
        sidebar_layout = QVBoxLayout(sidebar)
        sidebar_layout.setSpacing(0)
        sidebar_layout.setContentsMargins(0, 0, 0, 0)
        
        # Logo section
        logo_section = QWidget()
        logo_section.setStyleSheet(f"""
            QWidget {{
                background: {DesignTokens.Colors.BACKGROUND_DARK};
                padding: 20px;
            }}
        """)
        logo_layout = QVBoxLayout(logo_section)
        logo_layout.setSpacing(8)
        
        logo_label = QLabel("Docking Studio")
        logo_label.setStyleSheet(f"""
            QLabel {{
                color: {DesignTokens.Colors.WHITE};
                font-size: 18px;
                font-weight: 700;
            }}
        """)
        
        version_label = QLabel("Version 1.2.3")
        version_label.setStyleSheet(f"""
            QLabel {{
                color: {DesignTokens.Colors.GRAY_500};
                font-size: 11px;
            }}
        """)
        
        logo_layout.addWidget(logo_label)
        logo_layout.addWidget(version_label)
        sidebar_layout.addWidget(logo_section)
        
        # Navigation sections
        nav_container = QWidget()
        nav_container.setStyleSheet(f"""
            QWidget {{
                background: {DesignTokens.Colors.BACKGROUND_SIDEBAR};
            }}
        """)
        nav_layout = QVBoxLayout(nav_container)
        nav_layout.setSpacing(4)
        nav_layout.setContentsMargins(12, 16, 12, 16)
        
        # Section: Main
        section_label = self._create_section_label("MAIN")
        nav_layout.addWidget(section_label)
        
        main_navs = [
            ("Dashboard", "home", True),
            ("New Docking", "flask", False),
            ("Job Queue", "list", False),
        ]
        
        for title, icon, active in main_navs:
            nav_btn = self._create_sidebar_button(title, icon, active)
            nav_layout.addWidget(nav_btn)
        
        # Section: Analysis
        nav_layout.addSpacing(16)
        section_label = self._create_section_label("ANALYSIS")
        nav_layout.addWidget(section_label)
        
        analysis_navs = [
            ("Results", "chart", False),
            ("RMSD Analysis", "target", False),
            ("Interactions", "link", False),
        ]
        
        for title, icon, active in analysis_navs:
            nav_btn = self._create_sidebar_button(title, icon, active)
            nav_layout.addWidget(nav_btn)
        
        # Section: Tools
        nav_layout.addSpacing(16)
        section_label = self._create_section_label("TOOLS")
        nav_layout.addWidget(section_label)
        
        tools_navs = [
            ("File Converter", "convert", False),
            ("Molecule Viewer", "eye", False),
            ("AI Assistant", "bot", False),
        ]
        
        for title, icon, active in tools_navs:
            nav_btn = self._create_sidebar_button(title, icon, active)
            nav_layout.addWidget(nav_btn)
        
        nav_layout.addStretch()
        sidebar_layout.addWidget(nav_container)
        
        # Bottom status section
        bottom_section = QWidget()
        bottom_section.setStyleSheet(f"""
            QWidget {{
                background: {DesignTokens.Colors.BACKGROUND_DARK};
                padding: 16px;
            }}
        """)
        bottom_layout = QVBoxLayout(bottom_section)
        bottom_layout.setSpacing(8)
        
        # Settings button
        settings_btn = QPushButton("Settings")
        settings_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        settings_btn.setStyleSheet(f"""
            QPushButton {{
                background: transparent;
                color: {DesignTokens.Colors.GRAY_400};
                border: 1px solid {DesignTokens.Colors.GRAY_700};
                border-radius: 6px;
                padding: 8px 16px;
                text-align: left;
                font-size: 13px;
            }}
            QPushButton:hover {{
                background: {DesignTokens.Colors.GRAY_800};
                color: {DesignTokens.Colors.WHITE};
                border-color: {DesignTokens.Colors.GRAY_600};
            }}
        """)
        settings_btn.clicked.connect(self.on_preferences)
        bottom_layout.addWidget(settings_btn)
        
        # Docker status
        docker_status = self._create_docker_status_compact()
        bottom_layout.addWidget(docker_status)
        
        sidebar_layout.addWidget(bottom_section)
        
        return sidebar
    
    def _create_section_label(self, text: str) -> QLabel:
        """Create sidebar section label"""
        label = QLabel(text)
        label.setStyleSheet(f"""
            QLabel {{
                color: {DesignTokens.Colors.GRAY_600};
                font-size: 11px;
                font-weight: 600;
                letter-spacing: 1px;
                padding: 4px 8px;
            }}
        """)
        return label
    
    def _create_sidebar_button(self, title: str, icon_name: str, active: bool = False) -> QPushButton:
        """Create sidebar navigation button"""
        icon_map = {
            "home": "🏠",
            "flask": "🧪",
            "list": "📋",
            "chart": "📊",
            "target": "🎯",
            "link": "🔗",
            "convert": "🔄",
            "eye": "👁",
            "bot": "🤖",
        }
        icon = icon_map.get(icon_name, "•")
        
        bg_color = DesignTokens.Colors.PRIMARY if active else "transparent"
        text_color = DesignTokens.Colors.WHITE if active else DesignTokens.Colors.GRAY_300
        hover_bg = DesignTokens.Colors.GRAY_800
        
        btn = QPushButton(f"  {icon}  {title}")
        btn.setCursor(Qt.CursorShape.PointingHandCursor)
        btn.setStyleSheet(f"""
            QPushButton {{
                background: {bg_color};
                color: {text_color};
                border: none;
                border-radius: 8px;
                padding: 10px 12px;
                text-align: left;
                font-size: 13px;
                font-weight: 500;
            }}
            QPushButton:hover {{
                background: {hover_bg};
                color: {DesignTokens.Colors.WHITE};
            }}
        """)
        return btn
    
    def _create_docker_status_compact(self) -> QWidget:
        """Create compact Docker status indicator"""
        status_widget = QFrame()
        status_widget.setStyleSheet(f"""
            QFrame {{
                background: {DesignTokens.Colors.GRAY_900};
                border-radius: 8px;
                padding: 12px;
            }}
        """)
        
        layout = QHBoxLayout(status_widget)
        layout.setSpacing(8)
        
        indicator = QLabel("●")
        indicator.setStyleSheet(f"""
            QLabel {{
                color: {DesignTokens.Colors.SUCCESS};
                font-size: 10px;
            }}
        """)
        
        status_text = QLabel("Docker Running")
        status_text.setStyleSheet(f"""
            QLabel {{
                color: {DesignTokens.Colors.GRAY_400};
                font-size: 12px;
            }}
        """)
        
        layout.addWidget(indicator)
        layout.addWidget(status_text)
        layout.addStretch()
        
        return status_widget
    
    def _create_main_content(self) -> QWidget:
        """Create main content area with viewer and docking form"""
        content = QWidget()
        
        content_layout = QVBoxLayout(content)
        content_layout.setSpacing(16)
        content_layout.setContentsMargins(0, 0, 0, 0)
        
        # Main splitter: Docking Form | Molecular Viewer
        main_splitter = QSplitter(Qt.Orientation.Horizontal)
        main_splitter.setHandleWidth(1)
        
        # Left: Docking Form
        docking_form = self._create_docking_form()
        main_splitter.addWidget(docking_form)
        
        # Right: Molecular Viewer
        viewer_panel = self._create_viewer_panel()
        main_splitter.addWidget(viewer_panel)
        
        main_splitter.setStretchFactor(0, 3)
        main_splitter.setStretchFactor(1, 2)
        
        content_layout.addWidget(main_splitter, 1)
        
        return content
    
    def _create_viewer_panel(self) -> QWidget:
        """Create molecular viewer panel with all controls"""
        viewer_container = QWidget()
        viewer_container.setStyleSheet(f"""
            QWidget {{
                background: {DesignTokens.Colors.WHITE};
                border-left: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
            }}
        """)
        
        layout = QVBoxLayout(viewer_container)
        layout.setSpacing(0)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Viewer Header
        header = QWidget()
        header.setStyleSheet(f"""
            QWidget {{
                background: {DesignTokens.Colors.GRAY_50};
                border-bottom: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                padding: 12px 16px;
            }}
        """)
        header_layout = QHBoxLayout(header)
        header_layout.setSpacing(12)
        
        title = QLabel("3D Molecular Viewer")
        title.setStyleSheet(f"""
            QLabel {{
                font-size: 14px;
                font-weight: 700;
                color: {DesignTokens.Colors.TEXT_PRIMARY};
            }}
        """)
        header_layout.addWidget(title)
        
        header_layout.addStretch()
        
        # View control buttons
        view_btns = [
            ("⟳", "Auto Rotate", self.on_viewer_rotate),
            ("⊕", "Zoom Fit", self.on_viewer_zoom_fit),
            ("◎", "Center", self.on_viewer_center),
            ("↺", "Reset", self.on_viewer_reset),
        ]
        
        for icon, tooltip, handler in view_btns:
            btn = QPushButton(icon)
            btn.setToolTip(tooltip)
            btn.setFixedSize(32, 32)
            btn.setCursor(Qt.CursorShape.PointingHandCursor)
            btn.setStyleSheet(f"""
                QPushButton {{
                    background: {DesignTokens.Colors.WHITE};
                    border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                    border-radius: 6px;
                    font-size: 14px;
                }}
                QPushButton:hover {{
                    background: {DesignTokens.Colors.PRIMARY_10};
                    border-color: {DesignTokens.Colors.PRIMARY};
                }}
            """)
            btn.clicked.connect(handler)
            header_layout.addWidget(btn)
        
        layout.addWidget(header)
        
        # Style toolbar
        style_toolbar = QWidget()
        style_toolbar.setStyleSheet(f"""
            QWidget {{
                background: {DesignTokens.Colors.WHITE};
                padding: 8px 16px;
                border-bottom: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
            }}
        """)
        style_layout = QHBoxLayout(style_toolbar)
        style_layout.setSpacing(8)
        
        style_label = QLabel("Style:")
        style_label.setStyleSheet(f"""
            QLabel {{
                font-size: 12px;
                color: {DesignTokens.Colors.TEXT_SECONDARY};
                font-weight: 600;
            }}
        """)
        style_layout.addWidget(style_label)
        
        styles = [
            ("stick", "Stick"),
            ("ball", "Ball & Stick"),
            ("cartoon", "Cartoon"),
            ("surface", "Surface"),
            ("sphere", "Spacefill"),
        ]
        
        for style_id, label in styles:
            btn = QPushButton(label)
            btn.setCursor(Qt.CursorShape.PointingHandCursor)
            btn.setStyleSheet(f"""
                QPushButton {{
                    background: {DesignTokens.Colors.GRAY_100};
                    border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                    border-radius: 4px;
                    padding: 4px 12px;
                    font-size: 11px;
                    font-weight: 500;
                    color: {DesignTokens.Colors.TEXT_PRIMARY};
                }}
                QPushButton:hover {{
                    background: {DesignTokens.Colors.PRIMARY_10};
                    border-color: {DesignTokens.Colors.PRIMARY};
                }}
            """)
            btn.clicked.connect(lambda checked, s=style_id: self.on_viewer_style(s))
            style_layout.addWidget(btn)
        
        style_layout.addStretch()
        
        # Action buttons
        actions_label = QLabel("Actions:")
        actions_label.setStyleSheet(f"""
            QLabel {{
                font-size: 12px;
                color: {DesignTokens.Colors.TEXT_SECONDARY};
                font-weight: 600;
            }}
        """)
        style_layout.addWidget(actions_label)
        
        action_btns = [
            ("📷", "Screenshot", self.on_viewer_screenshot),
            ("🔗", "Interactions", self.on_viewer_interactions),
            ("⊞", "Overlay", self.on_viewer_overlay),
            ("⛶", "Fullscreen", self.on_viewer_fullscreen),
        ]
        
        for icon, tooltip, handler in action_btns:
            btn = QPushButton(icon)
            btn.setToolTip(tooltip)
            btn.setFixedSize(32, 32)
            btn.setCursor(Qt.CursorShape.PointingHandCursor)
            btn.setStyleSheet(f"""
                QPushButton {{
                    background: {DesignTokens.Colors.WHITE};
                    border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                    border-radius: 6px;
                    font-size: 14px;
                }}
                QPushButton:hover {{
                    background: {DesignTokens.Colors.PRIMARY_10};
                    border-color: {DesignTokens.Colors.PRIMARY};
                }}
            """)
            btn.clicked.connect(handler)
            style_layout.addWidget(btn)
        
        layout.addWidget(style_toolbar)
        
        # Viewer area (placeholder - would be AdvancedMolecularViewer)
        viewer_frame = QFrame()
        viewer_frame.setStyleSheet(f"""
            QFrame {{
                background: {DesignTokens.Colors.GRAY_900};
            }}
        """)
        viewer_layout = QVBoxLayout(viewer_frame)
        
        # Placeholder label for viewer
        placeholder = QWidget()
        placeholder_layout = QVBoxLayout(placeholder)
        placeholder_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        placeholder_label = QLabel("🧬")
        placeholder_label.setStyleSheet("""
            QLabel {
                font-size: 72px;
            }
        """)
        placeholder_layout.addWidget(placeholder_label)
        
        placeholder_text = QLabel("No structure loaded")
        placeholder_text.setStyleSheet(f"""
            QLabel {{
                font-size: 16px;
                color: {DesignTokens.Colors.GRAY_500};
                margin-top: 16px;
            }}
        """)
        placeholder_layout.addWidget(placeholder_text)
        
        placeholder_subtext = QLabel("Run a docking job to visualize results")
        placeholder_subtext.setStyleSheet(f"""
            QLabel {{
                font-size: 13px;
                color: {DesignTokens.Colors.GRAY_600};
            }}
        """)
        placeholder_layout.addWidget(placeholder_subtext)
        
        viewer_layout.addWidget(placeholder)
        
        # Store reference for later use
        self.viewer_placeholder = placeholder
        self.viewer_frame = viewer_frame
        
        layout.addWidget(viewer_frame, 1)
        
        # Status bar for viewer
        viewer_status = QWidget()
        viewer_status.setStyleSheet(f"""
            QWidget {{
                background: {DesignTokens.Colors.GRAY_50};
                border-top: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                padding: 8px 16px;
            }}
        """)
        status_layout = QHBoxLayout(viewer_status)
        status_layout.setSpacing(16)
        
        self.viewer_info = QLabel("Ready")
        self.viewer_info.setStyleSheet(f"""
            QLabel {{
                font-size: 12px;
                color: {DesignTokens.Colors.TEXT_SECONDARY};
            }}
        """)
        status_layout.addWidget(self.viewer_info)
        
        status_layout.addStretch()
        
        self.viewer_coords = QLabel("X: --  Y: --  Z: --")
        self.viewer_coords.setStyleSheet(f"""
            QLabel {{
                font-size: 11px;
                color: {DesignTokens.Colors.GRAY_500};
                font-family: monospace;
            }}
        """)
        status_layout.addWidget(self.viewer_coords)
        
        layout.addWidget(viewer_status)
        
        return viewer_container
    
    # Viewer control handlers
    def on_viewer_rotate(self):
        """Toggle auto-rotate"""
        logger.info("Viewer: Toggle auto-rotate")
        self.viewer_info.setText("Auto-rotate toggled")
    
    def on_viewer_zoom_fit(self):
        """Zoom to fit structure"""
        logger.info("Viewer: Zoom to fit")
        self.viewer_info.setText("Zoomed to fit")
    
    def on_viewer_center(self):
        """Center view on structure"""
        logger.info("Viewer: Center view")
        self.viewer_info.setText("View centered")
    
    def on_viewer_reset(self):
        """Reset viewer"""
        logger.info("Viewer: Reset view")
        self.viewer_info.setText("View reset")
    
    def on_viewer_style(self, style: str):
        """Change molecular style"""
        logger.info(f"Viewer: Style changed to {style}")
        self.viewer_info.setText(f"Style: {style}")
    
    def on_viewer_screenshot(self):
        """Take screenshot"""
        from PyQt6.QtWidgets import QFileDialog
        from datetime import datetime
        
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Screenshot",
            f"molecule_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png",
            "PNG Files (*.png)"
        )
        
        if file_path:
            logger.info(f"Viewer: Screenshot saved to {file_path}")
            self.viewer_info.setText("Screenshot saved")
            if hasattr(self, 'log_viewer'):
                self.log_viewer.append_log('INFO', f"Screenshot saved: {file_path}")
    
    def on_viewer_interactions(self):
        """Toggle interactions display"""
        logger.info("Viewer: Toggle interactions")
        self.viewer_info.setText("Interactions toggled")
    
    def on_viewer_overlay(self):
        """Toggle overlay mode"""
        logger.info("Viewer: Toggle overlay")
        self.viewer_info.setText("Overlay mode toggled")
    
    def on_viewer_fullscreen(self):
        """Toggle fullscreen"""
        logger.info("Viewer: Toggle fullscreen")
        if self.viewer_frame.isFullScreen():
            self.viewer_frame.showNormal()
        else:
            self.viewer_frame.showFullScreen()
    
    def _create_docking_form(self) -> QWidget:
        """Create docking form with tabbed interface"""
        form_widget = QWidget()
        form_layout = QVBoxLayout(form_widget)
        form_layout.setSpacing(16)
        form_layout.setContentsMargins(0, 0, 0, 0)
        
        # Create tab widget
        tab_widget = QTabWidget()
        tab_widget.setStyleSheet(f"""
            QTabWidget::pane {{
                border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                border-radius: 12px;
                background: {DesignTokens.Colors.WHITE};
            }}
            QTabBar::tab {{
                background: {DesignTokens.Colors.GRAY_100};
                color: {DesignTokens.Colors.TEXT_SECONDARY};
                padding: 12px 24px;
                font-size: 14px;
                font-weight: 600;
                border-top-left-radius: 10px;
                border-top-right-radius: 10px;
                margin-right: 4px;
            }}
            QTabBar::tab:selected {{
                background: {DesignTokens.Colors.WHITE};
                color: {DesignTokens.Colors.PRIMARY};
                border-bottom: 2px solid {DesignTokens.Colors.PRIMARY};
            }}
            QTabBar::tab:hover:!selected {{
                background: {DesignTokens.Colors.GRAY_200};
            }}
        """)
        
        # Tab 1: Input Files
        input_tab = QWidget()
        input_layout = QVBoxLayout(input_tab)
        input_layout.setSpacing(16)
        
        upload_splitter = QSplitter(Qt.Orientation.Horizontal)
        upload_splitter.setHandleWidth(16)
        
        receptor_card = self._create_receptor_upload_card()
        upload_splitter.addWidget(receptor_card)
        
        ligand_card = self._create_ligand_upload_card()
        upload_splitter.addWidget(ligand_card)
        
        input_layout.addWidget(upload_splitter)
        input_layout.addStretch()
        
        tab_widget.addTab(input_tab, "📁 Input Files")
        
        # Tab 2: Docking Settings
        settings_tab = QWidget()
        settings_layout = QVBoxLayout(settings_tab)
        settings_layout.setSpacing(16)
        
        engine_card = self._create_engine_selection_card()
        settings_layout.addWidget(engine_card)
        
        grid_card = self._create_grid_config_card()
        settings_layout.addWidget(grid_card)
        
        settings_layout.addStretch()
        
        tab_widget.addTab(settings_tab, "⚙ Docking Settings")
        
        # Tab 3: Hardware & Status
        hw_tab = QWidget()
        hw_layout = QVBoxLayout(hw_tab)
        hw_layout.setSpacing(16)
        
        gpu_card = self._create_gpu_status_card()
        hw_layout.addWidget(gpu_card)
        
        hw_layout.addStretch()
        
        tab_widget.addTab(hw_tab, "🖥 Hardware")
        
        # Tab 4: Viewer (inline preview)
        viewer_tab = QWidget()
        viewer_tab_layout = QVBoxLayout(viewer_tab)
        viewer_tab_layout.setSpacing(8)
        
        # Mini viewer controls
        mini_controls = QWidget()
        mini_controls.setStyleSheet(f"""
            QWidget {{
                background: {DesignTokens.Colors.GRAY_50};
                padding: 12px;
                border-bottom: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
            }}
        """)
        mini_layout = QHBoxLayout(mini_controls)
        mini_layout.setSpacing(8)
        
        mini_layout.addWidget(QLabel("Quick View:"))
        
        quick_btns = [
            ("Stick", "stick"),
            ("Cartoon", "cartoon"),
            ("Surface", "surface"),
        ]
        
        for label, style in quick_btns:
            btn = QPushButton(label)
            btn.setStyleSheet(f"""
                QPushButton {{
                    background: {DesignTokens.Colors.WHITE};
                    border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                    border-radius: 4px;
                    padding: 6px 12px;
                    font-size: 12px;
                }}
                QPushButton:hover {{
                    background: {DesignTokens.Colors.PRIMARY_10};
                }}
            """)
            btn.clicked.connect(lambda checked, s=style: self.on_viewer_style(s))
            mini_layout.addWidget(btn)
        
        mini_layout.addStretch()
        
        capture_btn = QPushButton("📷 Capture")
        capture_btn.setStyleSheet(f"""
            QPushButton {{
                background: {DesignTokens.Colors.WHITE};
                border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                border-radius: 4px;
                padding: 6px 12px;
                font-size: 12px;
            }}
            QPushButton:hover {{
                background: {DesignTokens.Colors.PRIMARY_10};
            }}
        """)
        capture_btn.clicked.connect(self.on_viewer_screenshot)
        mini_layout.addWidget(capture_btn)
        
        viewer_tab_layout.addWidget(mini_controls)
        
        # Placeholder for viewer
        placeholder_frame = QFrame()
        placeholder_frame.setStyleSheet(f"""
            QFrame {{
                background: {DesignTokens.Colors.GRAY_100};
            }}
        """)
        placeholder_layout = QVBoxLayout(placeholder_frame)
        placeholder_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        placeholder_icon = QLabel("🧬")
        placeholder_icon.setStyleSheet("font-size: 64px;")
        placeholder_layout.addWidget(placeholder_icon)
        
        placeholder_text = QLabel("Load a structure to preview")
        placeholder_text.setStyleSheet(f"""
            QLabel {{
                font-size: 14px;
                color: {DesignTokens.Colors.TEXT_SECONDARY};
                margin-top: 12px;
            }}
        """)
        placeholder_layout.addWidget(placeholder_text)
        
        viewer_tab_layout.addWidget(placeholder_frame, 1)
        
        tab_widget.addTab(viewer_tab, "👁 Viewer")
        
        # Add tab widget to form
        form_layout.addWidget(tab_widget, 1)
        
        # Submit section
        submit_frame = QWidget()
        submit_frame.setStyleSheet(f"""
            QWidget {{
                background: {DesignTokens.Colors.GRAY_50};
                border-top: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                padding: 16px 24px;
            }}
        """)
        submit_layout = QHBoxLayout(submit_frame)
        submit_layout.setSpacing(12)
        
        # Status info
        status_info = QLabel("Ready to start docking")
        status_info.setStyleSheet(f"""
            QLabel {{
                color: {DesignTokens.Colors.TEXT_SECONDARY};
                font-size: 13px;
            }}
        """)
        submit_layout.addWidget(status_info, 1)
        
        self.start_docking_btn = QPushButton("▶ Start Docking")
        self.start_docking_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self.start_docking_btn.setStyleSheet(f"""
            QPushButton {{
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                    stop:0 {DesignTokens.Colors.PRIMARY},
                    stop:1 {DesignTokens.Colors.SECONDARY});
                color: {DesignTokens.Colors.WHITE};
                border: none;
                border-radius: 10px;
                padding: 14px 32px;
                font-size: 15px;
                font-weight: 700;
                min-width: 180px;
            }}
            QPushButton:hover {{
                box-shadow: 0 4px 12px rgba(46, 90, 172, 0.4);
            }}
            QPushButton:disabled {{
                background: {DesignTokens.Colors.GRAY_300};
                box-shadow: none;
            }}
        """)
        self.start_docking_btn.clicked.connect(self.on_start_docking)
        submit_layout.addWidget(self.start_docking_btn)
        
        form_layout.addLayout(submit_layout)
        
        self.progress_container = QFrame()
        self.progress_container.setStyleSheet(f"""
            QFrame {{
                background: {DesignTokens.Colors.WHITE};
                border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
                border-radius: {DesignTokens.BorderRadius.LG}px;
                padding: {DesignTokens.Spacing.PADDING_MD}px;
            }}
        """)
        self.progress_container.setVisible(False)
        progress_layout = QVBoxLayout(self.progress_container)
        progress_layout.setSpacing(DesignTokens.Spacing.GAP_SM)
        
        self.progress_label = QLabel("Initializing...")
        self.progress_label.setStyleSheet(f"""
            QLabel {{
                font-size: {DesignTokens.Typography.BODY_L}px;
                font-weight: 600;
                color: {DesignTokens.Colors.TEXT_PRIMARY};
            }}
        """)
        progress_layout.addWidget(self.progress_label)
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setStyleSheet(f"""
            QProgressBar {{
                border: none;
                border-radius: {DesignTokens.BorderRadius.MD}px;
                height: 12px;
                background: {DesignTokens.Colors.GRAY_200};
                text-align: center;
            }}
            QProgressBar::chunk {{
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                    stop:0 {DesignTokens.Colors.PRIMARY},
                    stop:1 {DesignTokens.Colors.SECONDARY});
                border-radius: {DesignTokens.BorderRadius.MD}px;
            }}
        """)
        progress_layout.addWidget(self.progress_bar)
        
        self.cancel_docking_btn = QPushButton("Cancel")
        self.cancel_docking_btn.setStyleSheet(f"""
            QPushButton {{
                background: {DesignTokens.Colors.ERROR};
                color: {DesignTokens.Colors.WHITE};
                border: none;
                border-radius: {DesignTokens.BorderRadius.MD}px;
                padding: {DesignTokens.Spacing.PADDING_SM}px {DesignTokens.Spacing.PADDING_MD}px;
                font-weight: 600;
            }}
            QPushButton:hover {{
                background: #d32f2f;
            }}
        """)
        self.cancel_docking_btn.clicked.connect(self.on_cancel_docking)
        self.cancel_docking_btn.setVisible(False)
        progress_layout.addWidget(self.cancel_docking_btn, 0, Qt.AlignmentFlag.AlignRight)
        
        form_layout.addWidget(self.progress_container)
        
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
        
        self.gpu_timer = QTimer(self)
        self.gpu_timer.timeout.connect(self._detect_gpu)
        self.gpu_timer.start(10000)
        
        return card
    
    def _detect_gpu(self):
        """Detect GPU and update status - uses backend GPU endpoint"""
        try:
            from src.backend_connection import get_connection
            conn = get_connection()
            gpu_data = conn.get_gpu_status()
            
            if gpu_data and gpu_data.get('available') and gpu_data.get('gpus'):
                gpu = gpu_data['gpus'][0]
                self.gpu_name_label.setText(gpu.get('name', 'GPU'))
                
                util = gpu.get('utilization', 0)
                mem_used = gpu.get('memory_used', 0)
                mem_total = gpu.get('memory_total', 0)
                temp = gpu.get('temperature', 0)
                
                self.gpu_status_label.setText(f"Util: {util}% | VRAM: {mem_used}/{mem_total}MB | {temp}°C")
                self.gpu_status_label.setStyleSheet(f"""
                    QLabel {{
                        font-size: {DesignTokens.Typography.CAPTION_M}px;
                        color: {DesignTokens.Colors.SUCCESS};
                        font-weight: 600;
                    }}
                """)
                return
        except Exception as e:
            pass
        
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
        
        self.progress_container.setVisible(True)
        self.cancel_docking_btn.setVisible(True)
        self.progress_bar.setValue(0)
        self.progress_label.setText("Initializing docking job...")
        
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
        
        import uuid
        job_id = str(uuid.uuid4())
        
        try:
            from src.backend_connection import DockingWorker
            
            self.docking_worker = DockingWorker(
                job_id=job_id,
                progress_callback=self.on_docking_progress
            )
            self.docking_worker.start_docking(total_ligands=batch_size_value * 10)
            
            if hasattr(self, 'log_viewer'):
                self.log_viewer.append_log('INFO', f"Docking job started: {job_id[:8]}...")
                
        except Exception as e:
            logger.error(f"Failed to start docking: {e}")
            self.progress_label.setText(f"Error: {str(e)}")
            self.start_docking_btn.setEnabled(True)
            self.start_docking_btn.setText("Start Experiment")
            self.progress_container.setVisible(False)
    
    def on_docking_progress(self, progress: int, status: str, message: str):
        """Handle docking progress updates from SSE stream"""
        self.progress_bar.setValue(progress)
        self.progress_label.setText(message)
        
        if hasattr(self, 'log_viewer'):
            self.log_viewer.append_log('DEBUG', f"Progress: {progress}% - {message}")
        
        if status == 'completed':
            self.on_docking_complete()
        elif status == 'cancelled':
            self.on_docking_cancelled()
        elif status == 'error':
            self.on_docking_error(message)
    
    def on_docking_complete(self):
        """Handle docking completion"""
        logger.info("Docking completed successfully")
        self.progress_label.setText("Docking completed!")
        self.progress_bar.setValue(100)
        
        self.start_docking_btn.setEnabled(True)
        self.start_docking_btn.setText("Start Experiment")
        self.cancel_docking_btn.setVisible(False)
        
        if hasattr(self, 'log_viewer'):
            self.log_viewer.append_log('INFO', "Docking job completed successfully")
    
    def on_docking_cancelled(self):
        """Handle docking cancellation"""
        logger.info("Docking cancelled")
        self.progress_label.setText("Docking cancelled")
        
        self.start_docking_btn.setEnabled(True)
        self.start_docking_btn.setText("Start Experiment")
        self.cancel_docking_btn.setVisible(False)
        
        if hasattr(self, 'log_viewer'):
            self.log_viewer.append_log('WARNING', "Docking job cancelled by user")
    
    def on_docking_error(self, error_message: str):
        """Handle docking error"""
        logger.error(f"Docking error: {error_message}")
        self.progress_label.setText(f"Error: {error_message}")
        
        self.start_docking_btn.setEnabled(True)
        self.start_docking_btn.setText("Start Experiment")
        self.cancel_docking_btn.setVisible(False)
        
        if hasattr(self, 'log_viewer'):
            self.log_viewer.append_log('ERROR', f"Docking failed: {error_message}")
    
    def on_cancel_docking(self):
        """Handle cancel button click"""
        if hasattr(self, 'docking_worker') and self.docking_worker:
            self.docking_worker.cancel()
            logger.info("Docking job cancellation requested")
    
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
        """Create professional status bar with comprehensive info"""
        status_bar = QStatusBar()
        status_bar.setStyleSheet(f"""
            QStatusBar {{
                background: {DesignTokens.Colors.BACKGROUND_DARK};
                color: {DesignTokens.Colors.GRAY_300};
                border: none;
                padding: 6px 16px;
                font-size: 12px;
            }}
        """)
        
        self.setStatusBar(status_bar)
        
        # Left side: Connection status
        self.backend_status_label = QLabel("● Backend: Checking...")
        self.backend_status_label.setStyleSheet(f"""
            QLabel {{
                color: {DesignTokens.Colors.GRAY_400};
                padding: 4px 8px;
            }}
        """)
        status_bar.addWidget(self.backend_status_label)
        
        # Center: Job status
        self.job_status_label = QLabel("Ready")
        self.job_status_label.setStyleSheet(f"""
            QLabel {{
                color: {DesignTokens.Colors.GRAY_300};
                padding: 4px 12px;
            }}
        """)
        status_bar.addWidget(self.job_status_label)
        
        # Right side: Additional info
        info_container = QWidget()
        info_layout = QHBoxLayout(info_container)
        info_layout.setSpacing(16)
        info_layout.setContentsMargins(0, 0, 0, 0)
        
        # GPU status
        self.gpu_status_mini = QLabel("GPU: --")
        self.gpu_status_mini.setStyleSheet(f"""
            QLabel {{
                color: {DesignTokens.Colors.GRAY_500};
            }}
        """)
        info_layout.addWidget(self.gpu_status_mini)
        
        # AI status
        self.ai_status_mini = QLabel("AI: --")
        self.ai_status_mini.setStyleSheet(f"""
            QLabel {{
                color: {DesignTokens.Colors.GRAY_500};
            }}
        """)
        info_layout.addWidget(self.ai_status_mini)
        
        # Docker status
        self.docker_status_mini = QLabel("Docker: ●")
        self.docker_status_mini.setStyleSheet(f"""
            QLabel {{
                color: {DesignTokens.Colors.SUCCESS};
            }}
        """)
        info_layout.addWidget(self.docker_status_mini)
        
        # Version
        version_label = QLabel("v1.2.3")
        version_label.setStyleSheet(f"""
            QLabel {{
                color: {DesignTokens.Colors.GRAY_600};
                font-weight: 500;
            }}
        """)
        info_layout.addWidget(version_label)
        
        status_bar.addPermanentWidget(info_container)
        
        # Start periodic status updates
        self._update_status_bar()
        self.status_timer = QTimer(self)
        self.status_timer.timeout.connect(self._update_status_bar)
        self.status_timer.start(5000)
    
    def _update_status_bar(self):
        """Update status bar information"""
        try:
            from src.backend_connection import get_connection, BackendStatus
            
            conn = get_connection()
            info = conn.check_health()
            
            if info.status == BackendStatus.CONNECTED:
                self.backend_status_label.setText("● Backend: Connected")
                self.backend_status_label.setStyleSheet(f"""
                    QLabel {{
                        color: {DesignTokens.Colors.SUCCESS};
                        padding: 4px 8px;
                    }}
                """)
                
                if info.ollama_available:
                    self.ai_status_mini.setText(f"AI: {info.ollama_provider or 'Ollama'}")
                    self.ai_status_mini.setStyleSheet(f"""
                        QLabel {{ color: {DesignTokens.Colors.SUCCESS}; }}
                    """)
                else:
                    self.ai_status_mini.setText("AI: Offline")
                    self.ai_status_mini.setStyleSheet(f"""
                        QLabel {{ color: {DesignTokens.Colors.GRAY_500}; }}
                    """)
                
                self.docker_status_mini.setText("Docker: ● Running")
                self.docker_status_mini.setStyleSheet(f"""
                    QLabel {{ color: {DesignTokens.Colors.SUCCESS}; }}
                """)
            else:
                self.backend_status_label.setText("● Backend: Disconnected")
                self.backend_status_label.setStyleSheet(f"""
                    QLabel {{
                        color: {DesignTokens.Colors.ERROR};
                        padding: 4px 8px;
                    }}
                """)
                self.ai_status_mini.setText("AI: N/A")
                self.docker_status_mini.setText("Docker: ○ Stopped")
                self.docker_status_mini.setStyleSheet(f"""
                    QLabel {{ color: {DesignTokens.Colors.GRAY_500}; }}
                """)
        except Exception as e:
            logger.debug(f"Status bar update failed: {e}")
    
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
        
        # Setup log viewer dock
        self._setup_log_dock()
        
        # Setup job monitor dock
        self._setup_job_monitor_dock()
        
        logger.info("Chat dock widget initialized")
    
    def _setup_log_dock(self):
        """Setup dockable log viewer widget"""
        try:
            from src.ui.log_viewer import LogViewerWidget
            
            self.log_viewer = LogViewerWidget()
            
            self.log_dock = QDockWidget("Console / Logs", self)
            self.log_dock.setWidget(self.log_viewer)
            self.log_dock.setAllowedAreas(Qt.DockWidgetArea.BottomDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)
            self.log_dock.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetClosable | QDockWidget.DockWidgetFeature.DockWidgetMovable | QDockWidget.DockWidgetFeature.DockWidgetFloatable)
            self.log_dock.setMinimumHeight(150)
            self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.log_dock)
            
            self.log_viewer.append_log('INFO', 'Docking Studio initialized')
            self.log_viewer.append_log('INFO', 'Ready for docking jobs')
            
            logger.info("Log viewer dock initialized")
            
        except ImportError as e:
            logger.warning(f"Log viewer not available: {e}")
    
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
        
        from PyQt6.QtWidgets import QFileDialog
        from datetime import datetime
        
        file_path, selected_filter = QFileDialog.getSaveFileName(
            self,
            "Export Docking Report",
            f"docking_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
            "PDF Files (*.pdf);;HTML Files (*.html);;CSV Files (*.csv)"
        )
        
        if not file_path:
            return
        
        try:
            if file_path.endswith('.pdf'):
                self._export_pdf_report(file_path)
            elif file_path.endswith('.html'):
                self._export_html_report(file_path)
            elif file_path.endswith('.csv'):
                self._export_csv_report(file_path)
            
            if hasattr(self, 'log_viewer'):
                self.log_viewer.append_log('INFO', f"Report exported: {file_path}")
                
        except Exception as e:
            logger.error(f"Export failed: {e}")
            if hasattr(self, 'log_viewer'):
                self.log_viewer.append_log('ERROR', f"Export failed: {e}")
    
    def _export_pdf_report(self, file_path: str):
        """Export docking results as PDF report"""
        try:
            from reportlab.lib.pagesizes import letter
            from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
            from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
            from reportlab.lib import colors
            
            doc = SimpleDocTemplate(file_path, pagesize=letter)
            styles = getSampleStyleSheet()
            story = []
            
            title_style = ParagraphStyle(
                'CustomTitle',
                parent=styles['Title'],
                fontSize=24,
                spaceAfter=30,
            )
            story.append(Paragraph("Docking Studio - Scientific Report", title_style))
            story.append(Spacer(1, 12))
            
            from datetime import datetime
            story.append(Paragraph(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", styles['Normal']))
            story.append(Spacer(1, 20))
            
            story.append(Paragraph("Job Configuration", styles['Heading2']))
            config_data = []
            
            if hasattr(self, 'batch_size_slider'):
                config_data.append(['Batch Size', str(self.batch_size_slider.value())])
            if hasattr(self, 'exhaustiveness_slider'):
                config_data.append(['Exhaustiveness', str(self.exhaustiveness_slider.value())])
            
            if hasattr(self, 'center_inputs') and hasattr(self, 'size_inputs'):
                config_data.append(['Grid Center', f"({self.center_inputs[0].value()}, {self.center_inputs[1].value()}, {self.center_inputs[2].value()})"])
                config_data.append(['Grid Size', f"({self.size_inputs[0].value()}, {self.size_inputs[1].value()}, {self.size_inputs[2].value()})"])
            
            if config_data:
                t = Table(config_data, colWidths=[150, 200])
                t.setStyle(TableStyle([
                    ('BACKGROUND', (0, 0), (0, -1), colors.lightgrey),
                    ('TEXTCOLOR', (0, 0), (-1, -1), colors.black),
                    ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                    ('FONTNAME', (0, 0), (-1, -1), 'Helvetica'),
                    ('FONTSIZE', (0, 0), (-1, -1), 10),
                    ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
                    ('GRID', (0, 0), (-1, -1), 1, colors.grey),
                ]))
                story.append(t)
            
            story.append(Spacer(1, 20))
            story.append(Paragraph("Results Summary", styles['Heading2']))
            story.append(Paragraph("No docking results available yet. Run a docking job to generate results.", styles['Normal']))
            
            story.append(Spacer(1, 20))
            story.append(Paragraph("Security Status", styles['Heading2']))
            if hasattr(self, 'security_status_badge'):
                story.append(Paragraph(f"Status: {self.security_status_badge.text()}", styles['Normal']))
            
            story.append(Spacer(1, 20))
            story.append(Paragraph("GPU Status", styles['Heading2']))
            if hasattr(self, 'gpu_name_label'):
                story.append(Paragraph(f"GPU: {self.gpu_name_label.text()}", styles['Normal']))
                story.append(Paragraph(f"Status: {self.gpu_status_label.text()}", styles['Normal']))
            
            doc.build(story)
            logger.info(f"PDF report exported: {file_path}")
            
        except ImportError:
            logger.warning("reportlab not installed, falling back to HTML export")
            html_path = file_path.replace('.pdf', '.html')
            self._export_html_report(html_path)
    
    def _export_html_report(self, file_path: str):
        """Export docking results as HTML report"""
        from datetime import datetime
        
        html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>Docking Studio Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; background: #f5f5f5; }}
        .container {{ background: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        h1 {{ color: #1976d2; border-bottom: 2px solid #1976d2; padding-bottom: 10px; }}
        h2 {{ color: #424242; margin-top: 30px; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background: #1976d2; color: white; }}
        tr:nth-child(even) {{ background: #f9f9f9; }}
        .badge {{ display: inline-block; padding: 4px 12px; border-radius: 12px; font-size: 12px; }}
        .success {{ background: #e8f5e9; color: #2e7d32; }}
        .warning {{ background: #fff3e0; color: #e65100; }}
        .error {{ background: #ffebee; color: #c62828; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Docking Studio - Scientific Report</h1>
        <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        
        <h2>Job Configuration</h2>
        <table>
            <tr><th>Parameter</th><th>Value</th></tr>
"""
        
        if hasattr(self, 'batch_size_slider'):
            html_content += f"<tr><td>Batch Size</td><td>{self.batch_size_slider.value()}</td></tr>\n"
        if hasattr(self, 'exhaustiveness_slider'):
            html_content += f"<tr><td>Exhaustiveness</td><td>{self.exhaustiveness_slider.value()}</td></tr>\n"
        if hasattr(self, 'center_inputs') and hasattr(self, 'size_inputs'):
            html_content += f"<tr><td>Grid Center</td><td>({self.center_inputs[0].value()}, {self.center_inputs[1].value()}, {self.center_inputs[2].value()})</td></tr>\n"
            html_content += f"<tr><td>Grid Size</td><td>({self.size_inputs[0].value()}, {self.size_inputs[1].value()}, {self.size_inputs[2].value()})</td></tr>\n"
        
        html_content += """        </table>
        
        <h2>Results Summary</h2>
        <p>No docking results available yet. Run a docking job to generate results.</p>
        
        <h2>System Status</h2>
        <table>
            <tr><th>Component</th><th>Status</th></tr>
"""
        
        if hasattr(self, 'security_status_badge'):
            html_content += f"<tr><td>Security</td><td>{self.security_status_badge.text()}</td></tr>\n"
        if hasattr(self, 'gpu_name_label'):
            html_content += f"<tr><td>GPU</td><td>{self.gpu_name_label.text()} - {self.gpu_status_label.text()}</td></tr>\n"
        
        html_content += """        </table>
        
        <p style="margin-top: 40px; color: #666; font-size: 12px;">
            Generated by BioDockify Docking Studio
        </p>
    </div>
</body>
</html>"""
        
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        logger.info(f"HTML report exported: {file_path}")
    
    def _export_csv_report(self, file_path: str):
        """Export docking results as CSV"""
        import csv
        
        with open(file_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['Parameter', 'Value'])
            
            if hasattr(self, 'batch_size_slider'):
                writer.writerow(['Batch Size', self.batch_size_slider.value()])
            if hasattr(self, 'exhaustiveness_slider'):
                writer.writerow(['Exhaustiveness', self.exhaustiveness_slider.value()])
            if hasattr(self, 'center_inputs') and hasattr(self, 'size_inputs'):
                writer.writerow(['Grid Center X', self.center_inputs[0].value()])
                writer.writerow(['Grid Center Y', self.center_inputs[1].value()])
                writer.writerow(['Grid Center Z', self.center_inputs[2].value()])
                writer.writerow(['Grid Size X', self.size_inputs[0].value()])
                writer.writerow(['Grid Size Y', self.size_inputs[1].value()])
                writer.writerow(['Grid Size Z', self.size_inputs[2].value()])
            if hasattr(self, 'gpu_name_label'):
                writer.writerow(['GPU', self.gpu_name_label.text()])
                writer.writerow(['GPU Status', self.gpu_status_label.text()])
            if hasattr(self, 'security_status_badge'):
                writer.writerow(['Security', self.security_status_badge.text()])
        
        logger.info(f"CSV report exported: {file_path}")
    
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
