#!/bin/bash
# =========================================================
# Docking Studio v2.0 - Easy Start Script for Students
# =========================================================

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
RESET='\033[0m'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

print_banner() {
    echo ""
    echo -e "${CYAN}============================================================${RESET}"
    echo -e "${CYAN}  🧬 Docking Studio v2.0${RESET}  -  AI Drug Discovery"
    echo -e "${CYAN}============================================================${RESET}"
    echo ""
}

check_docker() {
    echo -e "${BOLD}Checking Docker...${RESET}"
    if ! command -v docker &> /dev/null; then
        echo -e "${RED}✗ Docker is not installed.${RESET}"
        echo ""
        echo "  Please install Docker Desktop:"
        echo "  → https://www.docker.com/products/docker-desktop"
        echo ""
        exit 1
    fi

    if ! docker info &> /dev/null; then
        echo -e "${RED}✗ Docker is not running.${RESET}"
        echo ""
        echo "  Please start Docker Desktop and wait for it to be ready."
        echo ""
        exit 1
    fi

    echo -e "${GREEN}✓ Docker is installed and running${RESET}"
}

start_services() {
    echo ""
    echo -e "${BOLD}Starting Docking Studio v2.0...${RESET}"
    echo -e "${YELLOW}  (First time: downloads all services - may take 10-15 minutes)${RESET}"
    echo ""

    docker compose up -d --build

    echo ""
    echo -e "${BOLD}Waiting for gateway to be ready...${RESET}"

    MAX_WAIT=120
    WAITED=0
    while [ $WAITED -lt $MAX_WAIT ]; do
        if curl -sf http://localhost:8000 &>/dev/null; then
            echo ""
            echo -e "${GREEN}============================================================${RESET}"
            echo -e "${GREEN}  ✅ Docking Studio v2.0 is ready!${RESET}"
            echo -e "${GREEN}============================================================${RESET}"
            echo ""
    echo -e "  🌐 ${BOLD}Open your browser and go to:${RESET}"
    echo -e "     → ${CYAN}http://localhost:8000${RESET}"
            echo ""
            return 0
        fi

        printf "."
        sleep 3
        WAITED=$((WAITED + 3))
    done

    echo ""
    echo -e "${RED}✗ Gateway failed to start within $MAX_WAIT seconds.${RESET}"
    echo ""
    echo "  Check logs with:"
    echo "  → docker compose logs gateway"
    echo "  → docker compose logs api-backend"
    echo ""
    return 1
}

open_browser() {
    case "$(uname -s)" in
        Darwin*)
            sleep 2 && open http://localhost:8000 &>/dev/null || true
            ;;
        MINGW*|MSYS*|CYGWIN*)
            sleep 2 && start http://localhost:8000 &>/dev/null || true
            ;;
        Linux*)
            sleep 2 && xdg-open http://localhost:3000 &>/dev/null || true
            ;;
    esac
}

show_status() {
    echo ""
    echo -e "${BOLD}Container Status:${RESET}"
    docker compose ps --format "table {{.Name}}\t{{.Status}}\t{{.Ports}}" 2>/dev/null || echo "  (docker compose not available)"
    echo ""
    echo -e "  🌐 Open:      ${CYAN}http://localhost:8000${RESET}"
    echo -e "  Stop:        ${CYAN}docker compose down${RESET}"
    echo -e "  View logs:   ${CYAN}docker compose logs -f${RESET}"
    echo -e "  Restart:     ${CYAN}./start.sh${RESET}"
    echo ""
}

# =========================================================
# Main
# =========================================================

print_banner
check_docker

if start_services; then
    show_status
    open_browser
else
    show_status
    exit 1
fi
