@echo off
REM BioDockify - Docker Status Check Script
REM Checks Docker Desktop installation and running status

setlocal enabledelayedexpansion

echo Checking Docker Desktop installation...
reg query "HKLM\SOFTWARE\Docker\Docker Desktop" >nul 2>&1
if %ERRORLEVEL% == 0 (
    echo [OK] Docker Desktop is installed
    reg query "HKLM\SOFTWARE\Docker\Docker Desktop" /v InstallLocation
) else (
    reg query "HKLM\SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\Docker Desktop" >nul 2>&1
    if %ERRORLEVEL% == 0 (
        echo [OK] Docker Desktop is installed (legacy)
    ) else (
        echo [WARNING] Docker Desktop is NOT installed
        echo.
        echo Please install Docker Desktop from:
        echo https://www.docker.com/products/docker-desktop/
        echo.
        pause
        exit /b 1
    )
)

echo.
echo Checking if Docker daemon is running...
docker ps >nul 2>&1
if %ERRORLEVEL% == 0 (
    echo [OK] Docker daemon is RUNNING
    docker version
) else (
    echo [WARNING] Docker daemon is NOT running
    echo.
    echo Please start Docker Desktop from:
    echo Start Menu ^> Docker Desktop
    echo.
)

pause
