#!/bin/bash
# Quick start script for running checkpoint tests
# Usage: ./tests/run_tests.sh

set -e

echo "========================================"
echo "THG Checkpoint Tests - Quick Start"
echo "========================================"
echo ""

# Check if we're in the project root
if [ ! -f "generate_data-base/generate_db.py" ]; then
    echo "Error: Must run from project root directory"
    exit 1
fi

# Check if pytest is installed
if ! python3 -c "import pytest" 2>/dev/null; then
    echo "Installing test dependencies..."
    pip install -r tests/requirements.txt
    echo ""
fi

echo "Running checkpoint tests..."
echo ""

# Run tests with verbose output
python3 -m pytest tests/test_checkpoint.py -v --tb=short

echo ""
echo "========================================"
echo "Test run complete!"
echo "========================================"
