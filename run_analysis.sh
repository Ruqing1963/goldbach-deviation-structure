#!/bin/bash
# Run Goldbach deviation analysis

echo "Goldbach Deviation Structure Analysis"
echo "======================================"

# Check Python
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 not found"
    exit 1
fi

# Install dependencies
echo "Installing dependencies..."
pip install -r code/requirements.txt

# Run analysis
echo "Running analysis..."
cd code
python3 goldbach_analysis.py

echo "Done!"
