"""
BioDockify Docking Studio - Main Window Stylesheet
Professional dark theme styles for main window
"""

MAIN_WINDOW_STYLE = """
/* Main Window */
QMainWindow {
    background-color: #F5F5F5;
    color: #FFFFFF;
}

/* Central Widget */
QWidget#centralwidget {
    background-color: #F5F5F5;
    color: #FFFFFF;
}

/* Status Bar */
QStatusBar {
    background-color: #E0E0E0;
    color: #333333;
    border-top: 1px solid #CCCCCC;
    padding: 5px;
}

/* Menu Bar */
QMenuBar {
    background-color: #2196F3;
    color: #FFFFFF;
    border-bottom: 2px solid #0B7CD5;
    padding: 5px;
}

QMenuBar::item {
    background-color: transparent;
    padding: 10px 20px;
}

QMenuBar::item:selected {
    background-color: #0B7CD5;
}

QMenuBar::item:pressed {
    background-color: #0A5F8D;
}

QMenu {
    background-color: #FFFFFF;
    color: #333333;
    border: 1px solid #E0E0E0;
    padding: 5px;
}

QMenu::item {
    background-color: transparent;
    padding: 8px 25px;
}

QMenu::item:selected {
    background-color: #2196F3;
    color: #FFFFFF;
}

QMenu::item:pressed {
    background-color: #0B7CD5;
    color: #FFFFFF;
}

/* Buttons */
QPushButton {
    background-color: #2196F3;
    color: white;
    border: none;
    padding: 10px 20px;
    border-radius: 5px;
    font-weight: bold;
    font-size: 14px;
}

QPushButton:hover {
    background-color: #0B7CD5;
}

QPushButton:pressed {
    background-color: #0A5F8D;
}

QPushButton:disabled {
    background-color: #CCCCCC;
    color: #888888;
}

/* Frame Borders */
QFrame {
    border: 2px solid #E0E0E0;
    border-radius: 8px;
    padding: 10px;
    background-color: #FFFFFF;
}

/* Labels */
QLabel {
    color: #333333;
    font-size: 12px;
}

/* Progress Bars */
QProgressBar {
    border: 2px solid #E0E0E0;
    border-radius: 5px;
    text-align: center;
    background-color: #F5F5F5;
    height: 25px;
}

QProgressBar::chunk {
    background-color: #4CAF50;
    border-radius: 3px;
}

QProgressBar::text {
    color: #333333;
    font-weight: bold;
}

/* Text Edit */
QTextEdit {
    background-color: #FFFFFF;
    color: #333333;
    border: 1px solid #CCCCCC;
    border-radius: 5px;
    padding: 5px;
    font-family: monospace;
    font-size: 11px;
}

/* Tables */
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
    color: #333333;
    font-weight: bold;
    padding: 5px;
    border: none;
    border-right: 1px solid #E0E0E0;
    border-bottom: 1px solid #E0E0E0;
}

/* Scroll Bars */
QScrollBar:vertical {
    background-color: #F5F5F5;
    width: 12px;
    border-radius: 6px;
}

QScrollBar::handle:vertical {
    background-color: #CCCCCC;
    border-radius: 3px;
    min-height: 30px;
}

QScrollBar::add-line:vertical {
    background: none;
}

QScrollBar::sub-line:vertical {
    background: none;
}

/* Group Boxes */
QGroupBox {
    border: 2px solid #E0E0E0;
    border-radius: 8px;
    padding: 10px;
    margin-top: 10px;
    font-weight: bold;
    color: #333333;
}

QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top left;
    padding: 0 5px 5px 5px;
    border-radius: 4px;
    border: none;
    background-color: #F5F5F5;
    margin: -10px;
}

/* Spin Boxes */
QSpinBox, QDoubleSpinBox {
    padding: 5px;
    border: 1px solid #CCCCCC;
    border-radius: 5px;
    background-color: #FFFFFF;
    color: #333333;
}

QSpinBox::up-button, QSpinBox::down-button,
QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {
    background-color: #F5F5F5;
    border: none;
    width: 20px;
    subcontrol-position: right;
    subcontrol-origin: margin;
}

QSpinBox::up-button:hover, QSpinBox::down-button:hover,
QDoubleSpinBox::up-button:hover, QDoubleSpinBox::down-button:hover {
    background-color: #E0E0E0;
}

/* Check Boxes */
QCheckBox {
    padding: 5px;
    background-color: #F5F5F5;
    color: #333333;
}

QCheckBox::indicator {
    width: 18px;
    height: 18px;
    border-radius: 3px;
}

QCheckBox::indicator:checked {
    background-color: #4CAF50;
    border: 1px solid #4CAF50;
}

QCheckBox::indicator:unchecked {
    background-color: #FFFFFF;
    border: 1px solid #CCCCCC;
}

/* Combo Boxes */
QComboBox {
    padding: 5px;
    border: 1px solid #CCCCCC;
    border-radius: 5px;
    background-color: #FFFFFF;
    color: #333333;
}

QComboBox QAbstractItemView {
    background-color: #FFFFFF;
    border: none;
    selection-background-color: #2196F3;
    color: white;
}
"""
