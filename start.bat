@echo off
REM ============================================================
REM BioDockify Studio AI v4.0.0 - Localhost Launcher (Windows)
REM No Docker required. Runs everything in one process.
REM ============================================================

echo.
echo ============================================================
echo   BioDockify Studio AI v4.0.0 - Localhost Mode
echo ============================================================
echo.

REM Check Python
python --version >nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo [ERROR] Python is not installed or not in PATH.
    echo Install Python 3.11+ from https://python.org
    pause
    exit /b 1
)

echo [OK] Python found
echo.
echo Starting BioDockify Studio AI...
echo Frontend: bundled (React build)
echo Backend:  all 160+ API routes loaded
echo.
echo Press Ctrl+C to stop the server.
echo.

python app/launcher.py
