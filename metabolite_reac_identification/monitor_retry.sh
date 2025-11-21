#!/bin/bash
# Monitor the retry process

LOG_FILE="retry_run_nov3.log"
ANNOTATION_FILE="reports/met_annotation.tsv"

echo "========================================"
echo "Retry Process Monitor"
echo "========================================"
echo ""

# Check if process is running
if ps aux | grep -q "[r]etry_failures_only.py"; then
    echo "Status: RUNNING ✓"
    PID=$(ps aux | grep "[r]etry_failures_only.py" | awk '{print $2}')
    echo "PID: $PID"
else
    echo "Status: NOT RUNNING"
fi

echo ""
echo "Current annotations:"
if [ -f "$ANNOTATION_FILE" ]; then
    COUNT=$(wc -l < "$ANNOTATION_FILE")
    echo "  Total: $COUNT metabolites"
else
    echo "  File not found"
fi

echo ""
echo "Log file size:"
if [ -f "$LOG_FILE" ]; then
    ls -lh "$LOG_FILE" | awk '{print "  " $5}'
    echo "  Lines: $(wc -l < $LOG_FILE)"
else
    echo "  No log file yet"
fi

echo ""
echo "Recent activity (last 10 lines):"
if [ -f "$LOG_FILE" ]; then
    tail -10 "$LOG_FILE"
else
    echo "  No log file yet"
fi

echo ""
echo "========================================"
