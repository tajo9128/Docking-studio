@echo off
title Docking Studio Launcher
color 0a

echo.
echo ============================================
echo   Docking Studio - Smart Launcher
echo ============================================
echo.

REM Check if Docker is running
echo [1/4] Checking Docker...
docker info >nul 2>&1
if errorlevel 1 (
    echo    Docker not running. Starting Docker Desktop...
    start "" "C:\Program Files\Docker\Docker\Docker Desktop.exe"
    echo    Waiting for Docker to start...
    :wait_docker
    docker info >nul 2>&1
    if errorlevel 1 (
        timeout /t 2 >nul
        goto wait_docker
    )
    echo    Docker started!
) else (
    echo    Docker is running!
)

REM Check if backend container is running
echo [2/4] Checking Backend Container...
docker ps --filter "name=docking-studio" --format "{{.Names}}" | findstr "docking-studio" >nul
if errorlevel 1 (
    echo    Starting backend container...
    docker run -d -p 8000:8000 --name docking-studio tajo9128/docking-studio:latest
) else (
    echo    Backend container already running!
)

REM Wait for API to be ready
echo [3/4] Waiting for API...
:wait_api
curl -s http://localhost:8000/health >nul 2>&1
if errorlevel 1 (
    timeout /t 1 >nul
    goto wait_api
)
echo    API is ready!

REM Check Python
echo [4/4] Checking Python...
python --version >nul 2>&1
if errorlevel 1 (
    echo.
    echo ERROR: Python is not installed!
    echo Please install Python 3.9+ from https://python.org
    pause
    exit /b 1
)

echo.
echo ============================================
echo   All Ready! Starting Docking Studio...
echo ============================================
echo.
echo Access points:
echo   - Web Dashboard: http://localhost:8000
echo   - API Docs:       http://localhost:8000/docs
echo.

REM Start the desktop app
cd /d "%~dp0Docking-studio"
python -m src.biodockify_main

if errorlevel 1 (
    echo.
    echo ERROR: Failed to start application
    pause
)
