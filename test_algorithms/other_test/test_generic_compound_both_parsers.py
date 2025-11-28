#!/usr/bin/env python3
"""
Test script to verify generic compound handling in both REST API and HTML parsers.

This tests whether both getCompParamFromRestAPI and getCompParam correctly handle 
generic compounds like C15602 (Quinone) and C15603 (Hydroquinone) by extracting 
alternative specific compounds from the COMMENT field and using their formulas.
"""

import sys
import os

# Add project root to path
project_root = os.path.join(os.path.dirname(__file__), '..', '..')
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from functions.function_bm_gdb import (
    getCompParamFromRestAPI, 
    getCompParam,
    getHtml,
    parse_kegg_flat_file
)
import urllib.request
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
LOGGER = logging.getLogger(__name__)

def test_rest_api_parser(compound_id):
    """Test REST API parser (getCompParamFromRestAPI)."""
    LOGGER.info(f"\n{'='*70}")
    LOGGER.info(f"Testing REST API parser for: {compound_id}")
    LOGGER.info(f"{'='*70}")
    
    try:
        # Fetch compound data from KEGG REST API
        url = f"https://rest.kegg.jp/get/{compound_id}"
        response = urllib.request.urlopen(url, timeout=10).read()
        text = response.decode("utf-8") if isinstance(response, bytes) else response
        
        # Parse with REST API parser
        EF = []  # Empty extra formulas list
        specialCompounds = "/tmp/special_compounds_test.txt"
        
        result = getCompParamFromRestAPI(text, compound_id, 10, EF, specialCompounds, None)
        
        # Extract formula
        formula = result[1]  # (formula1, formula2)
        
        LOGGER.info(f"Formula from REST API parser: {formula}")
        
        has_formula = bool(formula and len(formula) > 0 and formula[0] and str(formula[0]).strip())
        
        if has_formula:
            LOGGER.info(f"✓ REST API parser SUCCESS: {compound_id} has formula: {formula[0]}")
        else:
            LOGGER.error(f"✗ REST API parser FAILED: {compound_id} has no formula!")
        
        return has_formula
        
    except Exception as e:
        LOGGER.error(f"Error testing REST API parser for {compound_id}: {e}", exc_info=True)
        return False

def test_html_parser(compound_id):
    """Test HTML parser (getCompParam)."""
    LOGGER.info(f"\n{'='*70}")
    LOGGER.info(f"Testing HTML parser for: {compound_id}")
    LOGGER.info(f"{'='*70}")
    
    try:
        # Fetch compound data from KEGG HTML page
        url = f"http://www.genome.jp/dbget-bin/www_bget?cpd:{compound_id}"
        page = getHtml(url, 10)
        
        # Parse with HTML parser
        EF = []  # Empty extra formulas list
        specialCompounds = "/tmp/special_compounds_test.txt"
        
        result = getCompParam(page, compound_id, 10, EF, specialCompounds, None)
        
        # Extract formula
        formula = result[1]  # (formula1, formula2)
        
        LOGGER.info(f"Formula from HTML parser: {formula}")
        
        has_formula = bool(formula and len(formula) > 0 and formula[0] and str(formula[0]).strip())
        
        if has_formula:
            LOGGER.info(f"✓ HTML parser SUCCESS: {compound_id} has formula: {formula[0]}")
        else:
            LOGGER.error(f"✗ HTML parser FAILED: {compound_id} has no formula!")
        
        return has_formula
        
    except Exception as e:
        LOGGER.error(f"Error testing HTML parser for {compound_id}: {e}", exc_info=True)
        return False

def main():
    """Test both parsers with generic and specific compounds."""
    print("\n" + "="*70)
    print("GENERIC COMPOUND HANDLING TEST - BOTH PARSERS")
    print("="*70)
    
    test_compounds = [
        ("C15602", "Quinone (generic) - should resolve to C00472"),
        ("C15603", "Hydroquinone (generic) - should resolve to C00530"),
        ("C00472", "p-Benzoquinone (specific) - should have formula C6H4O2"),
        ("C00530", "Hydroquinone (specific) - should have formula C6H6O2"),
    ]
    
    results = {
        'rest_api': {},
        'html': {}
    }
    
    for cpd_id, description in test_compounds:
        print(f"\n{'='*70}")
        print(f"Testing: {cpd_id} - {description}")
        print(f"{'='*70}")
        
        # Test REST API parser
        results['rest_api'][cpd_id] = test_rest_api_parser(cpd_id)
        
        # Test HTML parser
        results['html'][cpd_id] = test_html_parser(cpd_id)
    
    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    
    all_passed = True
    
    print("\nREST API Parser Results:")
    for cpd_id, description in test_compounds:
        status = "✓ PASS" if results['rest_api'][cpd_id] else "✗ FAIL"
        print(f"  {status}: {cpd_id} - {description}")
        if not results['rest_api'][cpd_id]:
            all_passed = False
    
    print("\nHTML Parser Results:")
    for cpd_id, description in test_compounds:
        status = "✓ PASS" if results['html'][cpd_id] else "✗ FAIL"
        print(f"  {status}: {cpd_id} - {description}")
        if not results['html'][cpd_id]:
            all_passed = False
    
    print("="*70)
    if all_passed:
        print("✓ ALL TESTS PASSED - BOTH PARSERS WORKING")
        return 0
    else:
        print("✗ SOME TESTS FAILED")
        return 1

if __name__ == "__main__":
    sys.exit(main())
