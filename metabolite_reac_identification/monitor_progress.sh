#!/bin/bash

# Monitor script for metabolite_reac_identification.py with failure analysis

REPORTS_FILE="reports/met_annotation.tsv"
LOG_FILE="metabolite_run.log"
TOTAL_METABOLITES=4169

echo "=== Metabolite Identification Progress Monitor ==="
echo ""

# Check if process is running (matches both old and new script names)
if ps aux | grep -q "[m]etabolite_reac_identification"; then
    echo "✓ Script is RUNNING"
    PID=$(ps aux | grep "[m]etabolite_reac_identification" | grep -v grep | awk '{print $2}')
    echo "  PID: $PID"
else
    echo "✗ Script is NOT RUNNING"
fi

echo ""
echo "--- Progress ---"

# Count completed metabolites
if [ -f "$REPORTS_FILE" ]; then
    COMPLETED=$(wc -l < "$REPORTS_FILE")
    PERCENTAGE=$(awk "BEGIN {printf \"%.2f\", ($COMPLETED/$TOTAL_METABOLITES)*100}")
    REMAINING=$((TOTAL_METABOLITES - COMPLETED))
    echo "Metabolites annotated: $COMPLETED / $TOTAL_METABOLITES ($PERCENTAGE%)"
    echo "Remaining: $REMAINING"
    
    if [ $COMPLETED -gt 0 ]; then
        echo ""
        echo "Last 3 successful annotations:"
        tail -3 "$REPORTS_FILE" | cut -f1 | sed 's/^/  ✓ /'
    fi
else
    echo "No results file found yet"
fi

echo ""
echo "--- Failure Analysis ---"
if [ -f "$LOG_FILE" ]; then
    # Count different types of errors
    PUBCHEM_BUSY=$(grep -c "PubChem busy" "$LOG_FILE" 2>/dev/null || echo 0)
    FAILED_RETRIES=$(grep -c "Failed after .* retries" "$LOG_FILE" 2>/dev/null || echo 0)
    GET_CIDS_FAILED=$(grep -c "get_cids did not work" "$LOG_FILE" 2>/dev/null || echo 0)
    COMPOUND_ERRORS=$(grep -c "Error getting compound" "$LOG_FILE" 2>/dev/null || echo 0)
    DETAIL_ERRORS=$(grep -c "Error getting detailed compound info" "$LOG_FILE" 2>/dev/null || echo 0)
    
    echo "Error occurrences in log:"
    echo "  PubChem busy warnings:        $PUBCHEM_BUSY"
    echo "  Failed after max retries:     $FAILED_RETRIES"
    echo "  get_cids failures:            $GET_CIDS_FAILED"
    echo "  Compound retrieval errors:    $COMPOUND_ERRORS"
    echo "  Detailed info errors:         $DETAIL_ERRORS"
    
    # Show last few failures
    echo ""
    echo "Recent failures:"
    grep "Failed after\|get_cids extra did not work\|Error getting" "$LOG_FILE" | tail -5 | sed 's/^/  /'
else
    echo "No log file found"
fi

echo ""
echo "--- Recent Log Activity (last 5 lines) ---"
if [ -f "$LOG_FILE" ]; then
    tail -5 "$LOG_FILE" | grep -v "^$" | sed 's/^/  /'
else
    echo "No log file found"
fi

echo ""
echo "--- Commands ---"
echo "Run this monitor: ./monitor_progress.sh"
echo "Watch live log:   tail -f $LOG_FILE"
echo "Check results:    wc -l $REPORTS_FILE"
echo "Check failures:   wc -l reports/met_annotation_failures.tsv"
if ps aux | grep -q "[m]etabolite_reac_identification"; then
    PID=$(ps aux | grep "[m]etabolite_reac_identification" | grep -v grep | awk '{print $2}')
    echo "Stop script:      kill $PID"
fi
