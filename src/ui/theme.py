"""
BioDockify Docking Studio - Professional Design Theme
International software design system
"""

from PyQt6.QtGui import QColor, QPalette, QFont
from PyQt6.QtCore import Qt

class DesignTokens:
    """Design tokens for consistent styling"""
    
    # Color Palette - Professional International Style
    class Colors:
        # Primary Brand Colors (Scientific Blue)
        PRIMARY = "#2E5AAC"
        PRIMARY_HOVER = "#254A8F"
        PRIMARY_ACTIVE = "#1E3A70"
        PRIMARY_LIGHT = "#E8F0FF"
        PRIMARY_10 = "rgba(46, 90, 172, 0.1)"
        
        # Secondary Colors (Teal for science/chemistry)
        SECONDARY = "#00B8A3"
        SECONDARY_HOVER = "#00A092"
        SECONDARY_ACTIVE = "#008E81"
        
        # Semantic Colors
        SUCCESS = "#00C853"
        SUCCESS_BG = "#E8F5E9"
        WARNING = "#FFAB00"
        WARNING_BG = "#FFF3E0"
        ERROR = "#FF1744"
        ERROR_BG = "#FFEBEE"
        INFO = "#2196F3"
        INFO_BG = "#E3F2FD"
        
        # Neutral Colors
        WHITE = "#FFFFFF"
        GRAY_50 = "#FAFAFA"
        GRAY_100 = "#F5F5F5"
        GRAY_200 = "#EEEEEE"
        GRAY_300 = "#E0E0E0"
        GRAY_400 = "#BDBDBD"
        GRAY_500 = "#9E9E9E"
        GRAY_600 = "#757575"
        GRAY_700 = "#616161"
        GRAY_800 = "#424242"
        GRAY_900 = "#212121"
        BLACK = "#000000"
        
        # Background Colors
        BACKGROUND_PRIMARY = "#FFFFFF"
        BACKGROUND_SECONDARY = "#F5F7FA"
        BACKGROUND_TERTIARY = "#F0F2F5"
        BACKGROUND_CARD = "#FFFFFF"
        BACKGROUND_SIDEBAR = "#1A1F2C"
        BACKGROUND_DARK = "#12141D"
        
        # Text Colors
        TEXT_PRIMARY = "#1F2937"
        TEXT_SECONDARY = "#6B7280"
        TEXT_TERTIARY = "#9CA3AF"
        TEXT_INVERTED = "#FFFFFF"
        TEXT_DISABLED = "#D1D5DB"
        TEXT_PLACEHOLDER = "#9CA3AF"
        
        # Border Colors
        BORDER_LIGHT = "#E5E7EB"
        BORDER_DEFAULT = "#D1D5DB"
        BORDER_FOCUS = "#2E5AAC"
        BORDER_ERROR = "#EF4444"
        BORDER_SUCCESS = "#10B981"
    
    # Typography
    class Typography:
        # Font families
        PRIMARY = "SF Pro Display, -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif"
        MONO = "SF Mono, Monaco, 'Cascadia Code', 'Roboto Mono', monospace"
        
        # Font sizes (in points)
        DISPLAY_2XL = 32
        DISPLAY_XL = 28
        DISPLAY_L = 24
        DISPLAY_M = 20
        HEADING_XL = 18
        HEADING_L = 16
        HEADING_M = 14
        HEADING_S = 12
        BODY_L = 14
        BODY_M = 13
        BODY_S = 12
        CAPTION_L = 11
        CAPTION_M = 10
        
        # Font weights
        THIN = "QFont::Thin"
        EXTRA_LIGHT = "QFont::ExtraLight"
        LIGHT = "QFont::Light"
        REGULAR = "QFont::Normal"
        MEDIUM = "QFont::Medium"
        SEMI_BOLD = "QFont::DemiBold"
        BOLD = "QFont::Bold"
        EXTRA_BOLD = "QFont::ExtraBold"
        
        # Line heights
        LINE_HEIGHT_TIGHT = 1.2
        LINE_HEIGHT_NORMAL = 1.5
        LINE_HEIGHT_RELAXED = 1.75
    
    # Spacing
    class Spacing:
        # Spacing scale (4px base unit)
        XS = 4
        SM = 8
        MD = 16
        LG = 24
        XL = 32
        XXL = 48
        XXXL = 64
        
        # Specific spacing
        PADDING_XS = 8
        PADDING_SM = 12
        PADDING_MD = 16
        PADDING_LG = 24
        PADDING_XL = 32
        
        MARGIN_XS = 4
        MARGIN_SM = 8
        MARGIN_MD = 16
        MARGIN_LG = 24
        MARGIN_XL = 32
        
        GAP_XS = 4
        GAP_SM = 8
        GAP_MD = 16
        GAP_LG = 24
        GAP_XL = 32
    
    # Border Radius
    class BorderRadius:
        XS = 4
        SM = 6
        MD = 8
        LG = 12
        XL = 16
        FULL = 9999
    
    # Shadows
    class Shadows:
        NONE = "none"
        XS = "0 1px 2px 0 rgba(0, 0, 0, 0.05)"
        SM = "0 1px 3px 0 rgba(0, 0, 0, 0.1)"
        MD = "0 4px 6px -1px rgba(0, 0, 0, 0.1)"
        LG = "0 10px 15px -3px rgba(0, 0, 0, 0.1)"
        XL = "0 20px 25px -5px rgba(0, 0, 0, 0.1)"
        XXL = "0 25px 50px -12px rgba(0, 0, 0, 0.25)"
        INNER = "inset 0 2px 4px 0 rgba(0, 0, 0, 0.06)"
    
    # Transitions
    class Transitions:
        FAST = "150ms ease-in-out"
        NORMAL = "200ms ease-in-out"
        SLOW = "300ms ease-in-out"
        SLOWER = "500ms ease-in-out"
        EASE_OUT = "cubic-bezier(0, 0, 0.2, 1)"
        EASE_IN_OUT = "cubic-bezier(0.4, 0, 0.2, 1)"

class Stylesheet:
    """Ready-to-use stylesheet strings"""
    
    @staticmethod
    def get_application_stylesheet() -> str:
        """Get complete application stylesheet"""
        return f"""
        /* Global Reset */
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        /* Application */
        QMainWindow {{
            background-color: {DesignTokens.Colors.BACKGROUND_PRIMARY};
        }}
        
        QWidget {{
            font-family: {DesignTokens.Typography.PRIMARY};
            color: {DesignTokens.Colors.TEXT_PRIMARY};
            background-color: {DesignTokens.Colors.BACKGROUND_PRIMARY};
        }}
        
        /* Buttons - Primary */
        QPushButton {{
            background-color: {DesignTokens.Colors.PRIMARY};
            color: {DesignTokens.Colors.WHITE};
            border: none;
            border-radius: {DesignTokens.BorderRadius.MD}px;
            padding: {DesignTokens.Spacing.PADDING_SM}px {DesignTokens.Spacing.PADDING_LG}px;
            font-size: {DesignTokens.Typography.BODY_L}px;
            font-weight: 600;
            min-height: {DesignTokens.Spacing.XL}px;
            letter-spacing: 0.5px;
            text-transform: none;
        }}
        
        QPushButton:hover {{
            background-color: {DesignTokens.Colors.PRIMARY_HOVER};
            box-shadow: {DesignTokens.Shadows.SM};
        }}
        
        QPushButton:pressed {{
            background-color: {DesignTokens.Colors.PRIMARY_ACTIVE};
            box-shadow: {DesignTokens.Shadows.INNER};
        }}
        
        QPushButton:disabled {{
            background-color: {DesignTokens.Colors.GRAY_300};
            color: {DesignTokens.Colors.TEXT_DISABLED};
            opacity: 0.6;
        }}
        
        /* Buttons - Secondary */
        QPushButton[secondary="true"] {{
            background-color: {DesignTokens.Colors.WHITE};
            color: {DesignTokens.Colors.PRIMARY};
            border: 1px solid {DesignTokens.Colors.BORDER_DEFAULT};
            border-radius: {DesignTokens.BorderRadius.MD}px;
            padding: {DesignTokens.Spacing.PADDING_SM}px {DesignTokens.Spacing.PADDING_LG}px;
            font-size: {DesignTokens.Typography.BODY_L}px;
            font-weight: 500;
            min-height: {DesignTokens.Spacing.XL}px;
        }}
        
        QPushButton[secondary="true"]:hover {{
            background-color: {DesignTokens.Colors.PRIMARY_10};
            border-color: {DesignTokens.Colors.PRIMARY};
        }}
        
        /* Buttons - Success */
        QPushButton[success="true"] {{
            background-color: {DesignTokens.Colors.SUCCESS};
            color: {DesignTokens.Colors.WHITE};
        }}
        
        QPushButton[success="true"]:hover {{
            background-color: #00A664;
        }}
        
        /* Buttons - Danger */
        QPushButton[danger="true"] {{
            background-color: {DesignTokens.Colors.ERROR};
            color: {DesignTokens.Colors.WHITE};
        }}
        
        QPushButton[danger="true"]:hover {{
            background-color: #DC2626;
        }}
        
        /* Cards */
        QFrame[card="true"] {{
            background-color: {DesignTokens.Colors.BACKGROUND_CARD};
            border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
            border-radius: {DesignTokens.BorderRadius.LG}px;
            padding: {DesignTokens.Spacing.PADDING_LG}px;
        }}
        
        QFrame[card="true"]:hover {{
            border-color: {DesignTokens.Colors.BORDER_DEFAULT};
            box-shadow: {DesignTokens.Shadows.MD};
        }}
        
        /* Progress Bars */
        QProgressBar {{
            background-color: {DesignTokens.Colors.GRAY_200};
            border: none;
            border-radius: {DesignTokens.BorderRadius.FULL}px;
            height: 8px;
            text-align: center;
        }}
        
        QProgressBar::chunk {{
            background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                stop:0 {DesignTokens.Colors.PRIMARY},
                stop:1 {DesignTokens.Colors.SECONDARY});
            border-radius: {DesignTokens.BorderRadius.FULL}px;
        }}
        
        /* Labels */
        QLabel {{
            font-size: {DesignTokens.Typography.BODY_L}px;
            color: {DesignTokens.Colors.TEXT_PRIMARY};
        }}
        
        QLabel[heading="true"] {{
            font-size: {DesignTokens.Typography.HEADING_L}px;
            font-weight: 600;
            color: {DesignTokens.Colors.TEXT_PRIMARY};
        }}
        
        QLabel[subheading="true"] {{
            font-size: {DesignTokens.Typography.BODY_L}px;
            font-weight: 500;
            color: {DesignTokens.Colors.TEXT_SECONDARY};
        }}
        
        QLabel[caption="true"] {{
            font-size: {DesignTokens.Typography.CAPTION_L}px;
            color: {DesignTokens.Colors.TEXT_TERTIARY};
        }}
        
        /* Text Inputs */
        QLineEdit, QTextEdit {{
            background-color: {DesignTokens.Colors.WHITE};
            border: 1px solid {DesignTokens.Colors.BORDER_DEFAULT};
            border-radius: {DesignTokens.BorderRadius.MD}px;
            padding: {DesignTokens.Spacing.PADDING_SM}px;
            font-size: {DesignTokens.Typography.BODY_L}px;
            color: {DesignTokens.Colors.TEXT_PRIMARY};
            min-height: {DesignTokens.Spacing.XL}px;
        }}
        
        QLineEdit:focus, QTextEdit:focus {{
            border: 1px solid {DesignTokens.Colors.PRIMARY};
            outline: none;
        }}
        
        /* Spin Boxes */
        QSpinBox, QDoubleSpinBox {{
            background-color: {DesignTokens.Colors.WHITE};
            border: 1px solid {DesignTokens.Colors.BORDER_DEFAULT};
            border-radius: {DesignTokens.BorderRadius.MD}px;
            padding: {DesignTokens.Spacing.PADDING_SM}px;
            font-size: {DesignTokens.Typography.BODY_L}px;
            min-height: {DesignTokens.Spacing.XL}px;
        }}
        
        /* Tabs */
        QTabWidget::pane {{
            border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
            border-radius: {DesignTokens.BorderRadius.MD}px;
            background-color: {DesignTokens.Colors.BACKGROUND_CARD};
            top: -1px;
        }}
        
        QTabBar::tab {{
            background-color: {DesignTokens.Colors.BACKGROUND_SECONDARY};
            color: {DesignTokens.Colors.TEXT_SECONDARY};
            padding: {DesignTokens.Spacing.PADDING_MD}px {DesignTokens.Spacing.PADDING_XL}px;
            font-size: {DesignTokens.Typography.BODY_L}px;
            font-weight: 500;
            border-top-left-radius: {DesignTokens.BorderRadius.MD}px;
            border-top-right-radius: {DesignTokens.BorderRadius.MD}px;
            margin-right: {DesignTokens.Spacing.XS}px;
        }}
        
        QTabBar::tab:selected {{
            background-color: {DesignTokens.Colors.WHITE};
            color: {DesignTokens.Colors.PRIMARY};
            border-bottom: 2px solid {DesignTokens.Colors.PRIMARY};
        }}
        
        QTabBar::tab:hover:!selected {{
            background-color: {DesignTokens.Colors.GRAY_100};
        }}
        
        /* Tables */
        QTableView {{
            background-color: {DesignTokens.Colors.WHITE};
            border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
            border-radius: {DesignTokens.BorderRadius.MD}px;
            gridline-color: {DesignTokens.Colors.BORDER_LIGHT};
            selection-background-color: {DesignTokens.Colors.PRIMARY_10};
        }}
        
        QTableView::item:selected {{
            background-color: {DesignTokens.Colors.PRIMARY_10};
            color: {DesignTokens.Colors.PRIMARY};
        }}
        
        QHeaderView::section {{
            background-color: {DesignTokens.Colors.BACKGROUND_SECONDARY};
            color: {DesignTokens.Colors.TEXT_SECONDARY};
            padding: {DesignTokens.Spacing.PADDING_SM}px {DesignTokens.Spacing.PADDING_MD}px;
            font-weight: 600;
            font-size: {DesignTokens.Typography.CAPTION_L}px;
            border: none;
            border-bottom: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
            border-right: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
        }}
        
        /* Scroll Bars */
        QScrollBar:vertical {{
            background: {DesignTokens.Colors.GRAY_100};
            width: 12px;
            border-radius: 6px;
            margin: 0px;
        }}
        
        QScrollBar::handle:vertical {{
            background: {DesignTokens.Colors.GRAY_400};
            border-radius: 6px;
            min-height: 30px;
            margin: 2px;
        }}
        
        QScrollBar::handle:vertical:hover {{
            background: {DesignTokens.Colors.GRAY_500};
        }}
        
        QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
            height: 0px;
        }}
        
        QScrollBar:horizontal {{
            background: {DesignTokens.Colors.GRAY_100};
            height: 12px;
            border-radius: 6px;
            margin: 0px;
        }}
        
        QScrollBar::handle:horizontal {{
            background: {DesignTokens.Colors.GRAY_400};
            border-radius: 6px;
            min-width: 30px;
            margin: 2px;
        }}
        
        /* Group Boxes */
        QGroupBox {{
            font-size: {DesignTokens.Typography.HEADING_M}px;
            font-weight: 600;
            color: {DesignTokens.Colors.TEXT_PRIMARY};
            border: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
            border-radius: {DesignTokens.BorderRadius.MD}px;
            margin-top: {DesignTokens.Spacing.MARGIN_LG}px;
            padding: {DesignTokens.Spacing.PADDING_MD}px;
        }}
        
        QGroupBox::title {{
            subcontrol-origin: margin;
            subcontrol-position: top left;
            left: {DesignTokens.Spacing.PADDING_SM}px;
            padding: 0 {DesignTokens.Spacing.PADDING_XS}px 0 {DesignTokens.Spacing.PADDING_XS}px;
        }}
        
        /* Menu Bar */
        QMenuBar {{
            background-color: {DesignTokens.Colors.WHITE};
            border-bottom: 1px solid {DesignTokens.Colors.BORDER_LIGHT};
            padding: 0 {DesignTokens.Spacing.PADDING_MD}px;
        }}
        
        QMenuBar::item {{
            background-color: transparent;
            padding: {DesignTokens.Spacing.PADDING_SM}px {DesignTokens.Spacing.PADDING_MD}px;
            font-size: {DesignTokens.Typography.BODY_L}px;
            color: {DesignTokens.Colors.TEXT_PRIMARY};
            border-radius: {DesignTokens.BorderRadius.SM}px;
            margin: {DesignTokens.Spacing.MARGIN_XS}px;
        }}
        
        QMenuBar::item:selected {{
            background-color: {DesignTokens.Colors.PRIMARY_10};
            color: {DesignTokens.Colors.PRIMARY};
        }}
        
        /* Status Bar */
        QStatusBar {{
            background-color: {DesignTokens.Colors.BACKGROUND_SIDEBAR};
            color: {DesignTokens.Colors.TEXT_INVERTED};
            border: none;
            padding: {DesignTokens.Spacing.PADDING_SM}px {DesignTokens.Spacing.PADDING_MD}px;
            font-size: {DesignTokens.Typography.CAPTION_L}px;
        }}
        
        /* Tool Tips */
        QToolTip {{
            background-color: {DesignTokens.Colors.GRAY_900};
            color: {DesignTokens.Colors.WHITE};
            border: none;
            border-radius: {DesignTokens.BorderRadius.MD}px;
            padding: {DesignTokens.Spacing.PADDING_SM}px;
            font-size: {DesignTokens.Typography.CAPTION_L}px;
        }}
        """
