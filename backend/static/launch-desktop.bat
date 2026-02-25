@echo off
echo ============================================
echo   Docking Studio - Desktop Launcher
echo ============================================
echo.
echo Starting Docking Studio Desktop App...
echo.

REM Check if Python is available
python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python is not installed or not in PATH
    echo Please install Python 3.9+ from https://python.org
    pause
    exit /b 1
)

REM Check if pip is available
pip --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: pip is not installed
    pause
    exit /b 1
)

REM Install dependencies if needed
echo Installing dependencies...
pip install -q PyQt6 pyqtgraph 2>nul

REM Check if repo exists, if not clone
if not exist "Docking-studio" (
    echo Cloning Docking Studio repository...
    git clone https://github.com/tajo9128/Docking-studio.git
    cd Docking-studio
    pip install -q -r requirements.txt
) else (
    cd Docking-studio
)

REM Start the desktop app
echo Starting Desktop App...
python -m src.biodockify_main

if errorlevel 1 (
    echo.
    echo ERROR: Failed to start application
    pause
)
