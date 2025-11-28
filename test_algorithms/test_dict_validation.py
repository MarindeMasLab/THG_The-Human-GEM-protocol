#!/usr/bin/env python3
"""
Test script to verify dictionary access validation in equations_bm_gdb.py

This tests that the WrapRxnSubsProdParam and RxnParam2Eq functions properly
validate nested dictionary accesses and raise informative KeyErrors when
dictionaries are inconsistent.
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from functions.equations_bm_gdb import WrapRxnSubsProdParam, RxnParam2Eq


class MockMetabolite:
    """Mock metabolite object for testing"""
    def __init__(self, id1, formula1, formula2):
        self.ID1 = id1
        self.Formula1 = formula1
        self.Formula2 = formula2


class MockReaction:
    """Mock reaction object for testing"""
    def __init__(self, substrates, products):
        self._substrates = substrates
        self._products = products
    
    def Substrate(self):
        return self._substrates
    
    def Product(self):
        return self._products


def test_valid_access():
    """Test that valid dictionary access works correctly"""
    print("Test 1: Valid dictionary access...")
    
    # Setup valid dictionaries
    MetList = {
        'C00001': MockMetabolite('C00001', 'H2O', 'H2O'),
        'C00002': MockMetabolite('C00002', 'ATP', 'ATP'),
        'C00008': MockMetabolite('C00008', 'ADP', 'ADP'),
    }
    
    MetEquiv = {
        'G00001': 'C00001',  # Glycan -> Compound mapping
    }
    
    # Create mock reaction: G00001 + C00002 -> C00008
    substrates = [[1, '', 'G00001'], [1, '', 'C00002']]
    products = [[1, '', 'C00008']]
    reaction = MockReaction(substrates, products)
    
    try:
        DS, DP = WrapRxnSubsProdParam(reaction, MetList, MetEquiv)
        print("  ✓ WrapRxnSubsProdParam succeeded with valid data")
        print(f"    Substrates: {DS}")
        print(f"    Products: {DP}")
    except Exception as e:
        print(f"  ✗ Unexpected error: {e}")
        return False
    
    return True


def test_missing_equiv_in_metlist():
    """Test that missing MetEquiv mapping in MetList raises proper error"""
    print("\nTest 2: MetEquiv points to ID not in MetList...")
    
    # Setup INVALID dictionaries - MetEquiv references ID not in MetList
    MetList = {
        'C00001': MockMetabolite('C00001', 'H2O', 'H2O'),
        'C00002': MockMetabolite('C00002', 'ATP', 'ATP'),
        # Missing C00999!
    }
    
    MetEquiv = {
        'G00001': 'C00999',  # Points to non-existent metabolite!
    }
    
    # Create mock reaction: G00001 + C00002 -> C00001
    substrates = [[1, '', 'G00001'], [1, '', 'C00002']]
    products = [[1, '', 'C00001']]
    reaction = MockReaction(substrates, products)
    
    try:
        DS, DP = WrapRxnSubsProdParam(reaction, MetList, MetEquiv)
        print("  ✗ Should have raised KeyError but didn't!")
        return False
    except KeyError as e:
        if 'C00999' in str(e) and 'MetEquiv not found in MetList' in str(e):
            print(f"  ✓ Correctly raised KeyError: {e}")
            return True
        else:
            print(f"  ✗ KeyError raised but wrong message: {e}")
            return False
    except Exception as e:
        print(f"  ✗ Wrong exception type: {type(e).__name__}: {e}")
        return False


def test_missing_metabolite_in_metlist():
    """Test that direct metabolite ID not in MetList raises proper error"""
    print("\nTest 3: Direct metabolite ID not in MetList...")
    
    # Setup INVALID dictionaries - metabolite not in MetList
    MetList = {
        'C00001': MockMetabolite('C00001', 'H2O', 'H2O'),
        # Missing C00002!
    }
    
    MetEquiv = {}  # No glycan mappings
    
    # Create mock reaction: C00001 -> C00002
    substrates = [[1, '', 'C00001']]
    products = [[1, '', 'C00002']]  # This is missing from MetList!
    reaction = MockReaction(substrates, products)
    
    try:
        DS, DP = WrapRxnSubsProdParam(reaction, MetList, MetEquiv)
        print("  ✗ Should have raised KeyError but didn't!")
        return False
    except KeyError as e:
        if 'C00002' in str(e) and 'not found in MetList' in str(e):
            print(f"  ✓ Correctly raised KeyError: {e}")
            return True
        else:
            print(f"  ✗ KeyError raised but wrong message: {e}")
            return False
    except Exception as e:
        print(f"  ✗ Wrong exception type: {type(e).__name__}: {e}")
        return False


if __name__ == '__main__':
    print("=" * 70)
    print("Dictionary Access Validation Tests")
    print("=" * 70)
    
    results = []
    results.append(test_valid_access())
    results.append(test_missing_equiv_in_metlist())
    results.append(test_missing_metabolite_in_metlist())
    
    print("\n" + "=" * 70)
    print(f"Results: {sum(results)}/{len(results)} tests passed")
    print("=" * 70)
    
    if all(results):
        print("\n✓ All dictionary validation tests passed!")
        sys.exit(0)
    else:
        print("\n✗ Some tests failed!")
        sys.exit(1)
