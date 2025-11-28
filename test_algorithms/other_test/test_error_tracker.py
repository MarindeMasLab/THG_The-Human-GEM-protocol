#!/usr/bin/env python3
"""
Test script for ErrorTracker functionality.

This validates that the comprehensive error reporting system works correctly.
"""

import sys
import os

# Add project root to path
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

from functions.error_tracker import ErrorTracker

print("="*80)
print("TESTING ERROR TRACKER")
print("="*80)
print()

# Create tracker
tracker = ErrorTracker()

# Simulate various errors
print("Simulating errors...")

# Gene annotation errors
tracker.add_error('gene_annotation', 'Gene BRCA1 not found in Ensembl',
                  context={'gene_id': 'BRCA1', 'reaction': 'R00001'})
tracker.add_error('gene_annotation', 'Gene TP53 not found in Ensembl',
                  context={'gene_id': 'TP53', 'reaction': 'R00002'})
tracker.add_error('gene_annotation', 'Gene BRCA1 not found in Ensembl',
                  context={'gene_id': 'BRCA1', 'reaction': 'R00003'})
tracker.add_success('gene_annotation')
tracker.add_success('gene_annotation')

# Mass balance errors
tracker.add_error('mass_balance', 'Cannot balance equation - missing formula',
                  context={'reaction_id': 'R12345', 'equation': 'A + B -> C'})
tracker.add_success('mass_balance')
tracker.add_success('mass_balance')
tracker.add_success('mass_balance')

# BioCyc query errors
tracker.add_error('biocyc_query', 'Timeout querying BioCyc for gene location',
                  context={'gene': 'GENE1'}, level='warning')
tracker.add_success('biocyc_query')
tracker.add_success('biocyc_query')

# Ensembl API error
tracker.add_error('ensembl_api', 'Batch annotation fetch failed: Connection timeout',
                  context={'gene_count': 100}, level='error')
tracker.add_success('ensembl_api')

print("Errors simulated successfully!")
print()

# Generate summary
print("="*80)
print("GENERATING ERROR SUMMARY")
print("="*80)
print()

summary = tracker.generate_summary(verbose=False)
print(summary)

# Save reports
print("\n" + "="*80)
print("SAVING ERROR REPORTS")
print("="*80)

tracker.save_report('logs/test_error_report.txt', verbose=True)
tracker.export_json('logs/test_error_report.json')

print("\n✅ Test completed successfully!")
print(f"   Text report: logs/test_error_report.txt")
print(f"   JSON report: logs/test_error_report.json")
