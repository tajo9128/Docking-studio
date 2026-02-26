@echo off
setlocal enabledelayedexpansion

echo ============================================
echo   Docking Studio - Desktop Launcher
echo ============================================
echo.

REM Get current directory
set "SCRIPT_DIR=%~dp0"
set "SCRIPT_DIR=%SCRIPT_DIR:~0,-1%"
cd /d "%SCRIPT_DIR%"

REM Check if Python is available
python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python is not installed or not in PATH
    echo Please install Python 3.9+ from https://python.org
    pause
    exit /b 1
)

echo Python found:
python --version
echo.

REM Install core dependencies only (skip problematic ones)
echo Installing core dependencies...
pip install -q PyQt6 pyqtgraph requests 2>nul

REM Start the desktop app from current directory
echo Starting Desktop App...
echo.

python -m src.biodockify_main

if errorlevel 1 (
    echo.
    echo ERROR: Failed to start application
    pause
)
