#!/usr/bin/env python3
"""
Test script to verify the lambda pickle bug fix.

This script tests that:
1. Reaction objects with GPR2 methods can be pickled with dill
2. Unpickling works without NameError
3. The GPR2 accessor returns correct data after sanitization
"""

import dill
import sys
from types import MethodType


# Define the fixed accessor function (copy from pattern_generate_database.py)
def _rxn_gpr2_accessor(self, compartment=None):
    """Module-level GPR2 accessor to avoid lambda pickle issues."""
    try:
        if hasattr(self, '_gpr2_data'):
            if compartment is None:
                # Return entire dictionary (used when called as GPR2()[compartment])
                return self._gpr2_data
            else:
                # Return specific compartment (direct call with parameter)
                return self._gpr2_data.get(compartment, "")
        return {} if compartment is None else ""
    except Exception:
        return {} if compartment is None else ""


def sanitize_gpr2_accessor(rxn):
    """Sanitize GPR2 accessor for a reaction after unpickling."""
    try:
        # Check if GPR2 accessor exists and works
        if hasattr(rxn, 'GPR2') and callable(rxn.GPR2):
            try:
                # Try calling it to see if it works
                _ = rxn.GPR2()
            except NameError:
                # Bad closure - rebind to stable method
                if hasattr(rxn, '_gpr2_data'):
                    rxn.GPR2 = MethodType(_rxn_gpr2_accessor, rxn)
        elif hasattr(rxn, '_gpr2_data'):
            # No GPR2 method but has data - bind it
            rxn.GPR2 = MethodType(_rxn_gpr2_accessor, rxn)
    except Exception as e:
        print(f"Warning: Could not sanitize GPR2 for reaction: {e}")


class MockReaction:
    """Mock reaction class for testing"""
    def __init__(self, reaction_id):
        self.ID = reaction_id
        self.Subcel = [{"c": "gene1 or gene2", "m": "gene3 and gene4"}]


def test_gpr2_pickle():
    """Test that GPR2 method can be pickled and unpickled"""
    print("=" * 70)
    print("Testing Lambda Pickle Fix for GPR2 Method")
    print("=" * 70)
    
    # Create a mock reaction with compartment-specific GPRs
    rxn = MockReaction("R00001_c")
    rxn._gpr2_data = {"c": "gene1 or gene2", "m": "gene3 and gene4"}
    rxn.GPR2 = MethodType(_rxn_gpr2_accessor, rxn)
    
    print("\n1. Testing GPR2 accessor before pickling...")
    print(f"   Reaction ID: {rxn.ID}")
    print(f"   GPR2() returns: {rxn.GPR2()}")
    print(f"   GPR2()['c'] returns: {rxn.GPR2()['c']}")
    print(f"   GPR2()['m'] returns: {rxn.GPR2()['m']}")
    
    # Test pickling
    print("\n2. Pickling reaction object with dill...")
    try:
        pickled = dill.dumps(rxn)
        print("   ✓ Pickling successful")
    except Exception as e:
        print(f"   ✗ Pickling failed: {e}")
        return False
    
    # Test unpickling
    print("\n3. Unpickling reaction object...")
    try:
        rxn_restored = dill.loads(pickled)
        print("   ✓ Unpickling successful")
    except Exception as e:
        print(f"   ✗ Unpickling failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Sanitize the unpickled object (mimics what generate_db.py does)
    print("\n4. Sanitizing unpickled object...")
    sanitize_gpr2_accessor(rxn_restored)
    print("   ✓ Sanitization complete")
    
    # Test that unpickled object works correctly
    print("\n5. Testing GPR2 accessor after unpickling and sanitization...")
    try:
        result_dict = rxn_restored.GPR2()
        print(f"   GPR2() returns: {result_dict}")
        
        result_c = rxn_restored.GPR2()['c']
        print(f"   GPR2()['c'] returns: {result_c}")
        
        result_m = rxn_restored.GPR2()['m']
        print(f"   GPR2()['m'] returns: {result_m}")
        
        # Verify results match original
        assert result_dict == {"c": "gene1 or gene2", "m": "gene3 and gene4"}
        assert result_c == "gene1 or gene2"
        assert result_m == "gene3 and gene4"
        
        print("   ✓ All assertions passed")
    except Exception as e:
        print(f"   ✗ GPR2 accessor failed after unpickling: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    print("\n" + "=" * 70)
    print("✓ ALL TESTS PASSED - Lambda pickle bug is fixed!")
    print("=" * 70)
    return True


if __name__ == "__main__":
    success = test_gpr2_pickle()
    sys.exit(0 if success else 1)
