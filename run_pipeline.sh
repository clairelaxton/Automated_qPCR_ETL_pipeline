#!/bin/bash
# qPCR Pipeline Runner Script
# This script sets up a Python 3.14 virtual environment, installs dependencies,
# and runs the qPCR pipeline.

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Check if python3 is installed
echo -e "${YELLOW}Checking for Python 3.14...${NC}"
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}Error: python3 is not installed or not in PATH.${NC}"
    echo "Please install Python 3.14 and ensure it's available as 'python3'"
    exit 1
fi

# Verify Python version
PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
echo -e "${GREEN}Found Python ${PYTHON_VERSION}${NC}"

# Create virtual environment directory name
VENV_DIR="venv"

# Remove existing venv if it exists
if [ -d "$VENV_DIR" ]; then
    echo -e "${YELLOW}Removing existing virtual environment...${NC}"
    rm -rf "$VENV_DIR"
fi

# Create new virtual environment with python3
echo -e "${YELLOW}Creating virtual environment with Python 3.14...${NC}"
python3 -m venv "$VENV_DIR"

# Activate virtual environment
echo -e "${YELLOW}Activating virtual environment...${NC}"
source "$VENV_DIR/bin/activate"

# Upgrade pip
echo -e "${YELLOW}Upgrading pip...${NC}"
pip install --upgrade pip --quiet

# Install dependencies
echo -e "${YELLOW}Installing dependencies...${NC}"
if [ -f "requirements.txt" ]; then
    pip install -r requirements.txt --quiet
else
    echo -e "${RED}Error: requirements.txt not found!${NC}"
    exit 1
fi

# Run the pipeline
echo -e "${GREEN}Running qPCR pipeline...${NC}"
echo ""

# Pass any command-line arguments to the Python script
python qpcr_pipeline_main.py "$@"

# Deactivate virtual environment
deactivate

echo ""
echo -e "${GREEN}Pipeline execution completed!${NC}"
