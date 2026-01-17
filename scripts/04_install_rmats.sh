#!/bin/bash

# Install rMATS-turbo from source
# Run from project root: bash scripts/04_install_rmats.sh
# Note: Dependencies should be installed via environment.yml first

set -e

echo "Installing rMATS-turbo..."

cd ~
if [ ! -d "rmats-turbo" ]; then
    git clone https://github.com/Xinglab/rmats-turbo.git
fi

cd rmats-turbo
git pull || true

# Build
echo "Building rMATS..."
./build_rmats

# Test
echo "Testing installation..."
./rmats.py --version

# Symlink ins PATH
mkdir -p ~/.local/bin
ln -sf ~/rmats-turbo/rmats.py ~/.local/bin/rmats.py

echo ""
echo "rMATS installation complete!"
echo "Location: ~/rmats-turbo/rmats.py"
echo "Symlink: ~/.local/bin/rmats.py"
