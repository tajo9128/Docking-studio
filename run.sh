#!/bin/bash

# Docking Studio Launcher for Linux/Mac

echo "============================================"
echo "  Docking Studio - Smart Launcher"
echo "============================================"
echo ""

# Check Docker
echo "[1/4] Checking Docker..."
if ! command -v docker &> /dev/null; then
    echo "ERROR: Docker is not installed!"
    echo "Please install Docker from https://docker.com"
    exit 1
fi

# Start Docker if not running
docker info > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "Starting Docker..."
    # Try to start Docker
    if command -v open &> /dev/null; then
        open -a Docker
    elif command -v systemctl &> /dev/null; then
        sudo systemctl start docker
    fi
    
    echo "Waiting for Docker..."
    while ! docker info > /dev/null 2>&1; do
        sleep 2
    done
    echo "Docker started!"
else
    echo "Docker is running!"
fi

# Start backend container
echo "[2/4] Starting Backend Container..."
if docker ps | grep -q docking-studio; then
    echo "Backend already running!"
else
    docker run -d -p 8000:8000 --name docking-studio tajo9128/docking-studio:latest
    echo "Backend started!"
fi

# Wait for API
echo "[3/4] Waiting for API..."
while ! curl -s http://localhost:8000/health > /dev/null 2>&1; do
    sleep 1
done
echo "API is ready!"

# Check Python
echo "[4/4] Checking Python..."
if ! command -v python &> /dev/null; then
    echo "ERROR: Python is not installed!"
    exit 1
fi

echo ""
echo "============================================"
echo "  All Ready! Starting Docking Studio..."
echo "============================================"
echo ""
echo "Access points:"
echo "  - Web Dashboard: http://localhost:8000"
echo "  - API Docs:       http://localhost:8000/docs"
echo ""

# Start desktop app
cd "$(dirname "$0")/Docking-studio"
python -m src.biodockify_main
