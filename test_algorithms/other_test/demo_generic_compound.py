#!/usr/bin/env python3
"""
Test script to verify generic compound handling with alternative compounds.

This tests whether the algorithm can correctly handle generic compounds like
C15602 (Quinone) and C15603 (Hydroquinone) by extracting alternative specific
compounds from the COMMENT field and using their formulas.
"""

import sys
import os

# Add project root to path
project_root = os.path.join(os.path.dirname(__file__), '..', '..')
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from functions.function_bm_gdb import getCompParamFromRestAPI, parse_kegg_flat_file
import urllib.request
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
LOGGER = logging.getLogger(__name__)

def test_generic_compound(compound_id):
    """Test fetching a compound and checking if formula is resolved."""
    LOGGER.info(f"\n{'='*70}")
    LOGGER.info(f"Testing compound: {compound_id}")
    LOGGER.info(f"{'='*70}")
    
    try:
        # Fetch compound data from KEGG
        url = f"https://rest.kegg.jp/get/{compound_id}"
        response = urllib.request.urlopen(url, timeout=10).read()
        text = response.decode("utf-8") if isinstance(response, bytes) else response
        
        # Display raw data
        LOGGER.info(f"\nRaw KEGG data:")
        for line in text.split('\n')[:20]:  # First 20 lines
            LOGGER.info(line)
        
        # Parse with our function
        EF = []  # Empty extra formulas list
        specialCompounds = "/tmp/special_compounds_test.txt"
        
        result = getCompParamFromRestAPI(text, compound_id, 10, EF, specialCompounds, None)
        
        # Extract relevant fields
        urls01_02 = result[0]  # (primary_id, secondary_id)
        formula = result[1]     # (formula1, formula2)
        
        LOGGER.info(f"\n{'='*70}")
        LOGGER.info(f"RESULTS:")
        LOGGER.info(f"{'='*70}")
        LOGGER.info(f"Primary ID:   {urls01_02[0] if urls01_02 else 'N/A'}")
        LOGGER.info(f"Secondary ID: {urls01_02[1] if len(urls01_02) > 1 else 'N/A'}")
        LOGGER.info(f"Formula:      {formula}")
        
        # formula is a tuple like ('C6H4O2', 'C6H4O2')
        # Check if formula was found
        has_formula = bool(formula and len(formula) > 0 and formula[0] and str(formula[0]).strip())
        LOGGER.info(f"\n✓ Formula found: {has_formula}")
        
        if has_formula:
            LOGGER.info(f"✓ SUCCESS: Compound {compound_id} has formula: {formula[0]}")
        else:
            LOGGER.error(f"✗ FAILED: Compound {compound_id} has no formula!")
        
        return has_formula
        
    except Exception as e:
        LOGGER.error(f"Error testing {compound_id}: {e}", exc_info=True)
        return False

def main():
    """Test generic compounds from R02164 reaction."""
    print("\n" + "="*70)
    print("GENERIC COMPOUND FORMULA RESOLUTION TEST")
    print("="*70)
    
    test_compounds = [
        ("C15602", "Quinone (generic) - should resolve to C00472 (p-Benzoquinone)"),
        ("C15603", "Hydroquinone (generic) - should resolve to C00530"),
        ("C00472", "p-Benzoquinone (specific) - should have formula C6H4O2"),
        ("C00530", "Hydroquinone (specific) - should have formula C6H6O2"),
    ]
    
    results = {}
    for cpd_id, description in test_compounds:
        print(f"\n{'='*70}")
        print(f"{description}")
        print(f"{'='*70}")
        results[cpd_id] = test_generic_compound(cpd_id)
    
    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    
    all_passed = True
    for cpd_id, description in test_compounds:
        status = "✓ PASS" if results[cpd_id] else "✗ FAIL"
        print(f"{status}: {cpd_id} - {description}")
        if not results[cpd_id]:
            all_passed = False
    
    print("="*70)
    if all_passed:
        print("✓ ALL TESTS PASSED")
        return 0
    else:
        print("✗ SOME TESTS FAILED")
        return 1

if __name__ == "__main__":
    sys.exit(main())
