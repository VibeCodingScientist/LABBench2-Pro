#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Deploy LABBench2-Pro to a fresh Ubuntu VPS
#
# Usage: scp this file to the server, then run:
#   chmod +x deploy.sh && ./deploy.sh
#
# Prerequisites on server: nothing (this script installs everything)
###############################################################################

echo "=== LABBench2-Pro Server Deployment ==="

# Install Docker if missing
if ! command -v docker &>/dev/null; then
    echo "Installing Docker..."
    curl -fsSL https://get.docker.com | sh
    sudo usermod -aG docker "$USER"
    echo "Docker installed. You may need to log out and back in for group changes."
fi

# Install Python 3.11+ if missing
if ! command -v python3.11 &>/dev/null && ! python3 --version 2>&1 | grep -qE "3\.1[1-9]"; then
    echo "Installing Python 3.11..."
    sudo apt-get update
    sudo apt-get install -y python3.11 python3.11-venv python3.11-dev python3-pip
fi

# Install psql client for schema setup
if ! command -v psql &>/dev/null; then
    echo "Installing PostgreSQL client..."
    sudo apt-get update
    sudo apt-get install -y postgresql-client
fi

# Install tmux for persistent sessions
if ! command -v tmux &>/dev/null; then
    echo "Installing tmux..."
    sudo apt-get install -y tmux
fi

# Clone repo (or pull if already cloned)
REPO_DIR="$HOME/LABBench2-Pro"
if [ -d "$REPO_DIR" ]; then
    echo "Updating existing repo..."
    cd "$REPO_DIR"
    git pull
else
    echo "Cloning repo..."
    git clone https://github.com/VibeCodingScientist/LABBench2-Pro.git "$REPO_DIR"
    cd "$REPO_DIR"
fi

# Create venv and install deps
echo "Setting up Python environment..."
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -e .

# Prompt for .env if missing
if [ ! -f .env ]; then
    echo ""
    echo "=== .env file not found ==="
    echo "Copy .env.example to .env and fill in your keys:"
    echo "  cp .env.example .env"
    echo "  nano .env"
    echo ""
    echo "Required:"
    echo "  ANTHROPIC_API_KEY=sk-ant-..."
    echo "  HF_TOKEN=hf_... (get from https://huggingface.co/settings/tokens)"
    echo ""
    cp .env.example .env
    echo "Created .env from template. Edit it, then run:"
    echo "  ./run_all.sh --model claude-opus-4.6"
    exit 0
fi

# Start services
echo "Starting Docker services..."
docker compose up -d

echo ""
echo "=== Deployment complete ==="
echo ""
echo "To run the full pipeline in a tmux session:"
echo "  tmux new -s labbench"
echo "  source .venv/bin/activate"
echo "  ./run_all.sh --model claude-opus-4.6"
echo ""
echo "Detach with Ctrl+B then D. Reattach with: tmux attach -t labbench"
