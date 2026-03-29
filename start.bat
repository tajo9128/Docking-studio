@echo off
rem ============================================================
rem Docking Studio v2.0 - Easy Start Script for Windows
rem ============================================================

setlocal enabledelayedexpansion

set "RED=[31m"
set "GREEN=[32m"
set "YELLOW=[33m"
set "CYAN=[36m"
set "BOLD=[1m"
set "RESET=[0m"

echo.
echo %CYAN%============================================================%RESET%
echo %CYAN%  Docking Studio v2.0  -  AI Drug Discovery%RESET%
echo %CYAN%============================================================%RESET%
echo.

rem Check Docker
echo %BOLD%Checking Docker...%RESET%
docker info >nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo %RED%  ✗ Docker is not running.%RESET%
    echo.
    echo   Please start Docker Desktop and wait for it to be ready.
    echo.
    pause
    exit /b 1
)
echo %GREEN%  ✓ Docker is running%RESET%

rem Start
echo.
echo %BOLD%Starting Docking Studio v2.0...%RESET%
echo %YELLOW%  (First time: downloads all services - may take 10-15 min)%RESET%
echo.

docker compose up -d --build
if %ERRORLEVEL% neq 0 (
    echo %RED%✗ Build failed.%RESET%
    echo.
    echo   Check logs with: docker compose logs
    pause
    exit /b 1
)

rem Wait for health
echo.
echo %BOLD%Waiting for gateway to be ready...%RESET%

set /a "waited=0"
set "max_wait=120"

:wait_loop
    curl -sf http://localhost:8000 >nul 2>&1
    if %ERRORLEVEL% equ 0 goto :ready

    timeout /t 3 /nobreak >nul
    set /a "waited+=3"
    if !waited! lss %max_wait% (
        powershell -Command "Write-Host -NoNewline '.'"
        goto :wait_loop
    )

:ready
echo.
echo.
echo %GREEN%============================================================%RESET%
echo %GREEN%  Docking Studio v2.0 is ready!%RESET%
echo %GREEN%============================================================%RESET%
echo.
echo   🌐 Open your browser and go to:
echo   %CYAN%   http://localhost:8000%RESET%
echo.
echo -----------------------------------------------------------
echo   Stop:        docker compose down
echo   View logs:   docker compose logs -f
echo   Restart:     start.bat
echo -----------------------------------------------------------
echo.
echo Opening browser...
start http://localhost:8000
pause
