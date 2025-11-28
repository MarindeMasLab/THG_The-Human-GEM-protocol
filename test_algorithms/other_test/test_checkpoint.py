"""
Unit tests for checkpoint save/load functionality.

These tests verify that:
1. Checkpoint files can be saved successfully
2. Checkpoint files can be loaded without errors
3. Data integrity is maintained (reactions, metabolites match)
4. Reaction methods work after reload (no NameError from lambdas)
5. Compartmentalized reactions preserve their data
6. Failed location reactions are properly filtered

Run with: pytest tests/test_checkpoint.py -v
"""

import os
import sys
import tempfile
import pytest
import dill
from types import MethodType
from copy import deepcopy

# Add project to path
project_root = os.path.join(os.path.dirname(__file__), '..', '..')
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from functions.class_generate_database import reaction, compound, pathway
sys.path.insert(0, os.path.join(project_root, 'generate_data-base'))
from generate_db import sanitize_loaded_reactions, _rxn_substrate, _rxn_product


class TestReactionPickle:
    """Test that reaction objects with MethodType methods pickle correctly."""
    
    def test_reaction_pickle_no_lambda_error(self):
        """Test that reactions with MethodType methods pickle without NameError."""
        # Create a reaction object
        rxn = reaction("http://example.com/R00001", 0.5, "R00001", "TestPathway", 0)
        
        # Simulate substrate and product data
        substrate_data = [
            (1.0, "C00001", "C00001"),
            (1.0, "C00002", "C00002")
        ]
        product_data = [
            (1.0, "C00003", "C00003"),
            (2.0, "C00004", "C00004")
        ]
        
        # Bind using the safe method (how generate_db.py does it)
        rxn.subs = substrate_data
        rxn.prods = product_data
        rxn.Substrate = MethodType(_rxn_substrate, rxn)
        rxn.Product = MethodType(_rxn_product, rxn)
        rxn.SetSubstrate = MethodType(_rxn_substrate, rxn)
        rxn.SetProduct = MethodType(_rxn_product, rxn)
        
        # Pickle and unpickle
        with tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix='.pkl') as f:
            temp_file = f.name
            dill.dump(rxn, f)
        
        try:
            # Load in a new context (simulating checkpoint resume)
            with open(temp_file, 'rb') as f:
                loaded_rxn = dill.load(f)
            
            # Test that methods work (no NameError)
            subs = loaded_rxn.Substrate()
            prods = loaded_rxn.Product()
            
            # Verify data integrity
            assert subs == substrate_data, "Substrate data not preserved"
            assert prods == product_data, "Product data not preserved"
            assert len(subs) == 2, "Wrong number of substrates"
            assert len(prods) == 2, "Wrong number of products"
            
            # Verify stoichiometry preserved
            assert subs[0][0] == 1.0
            assert prods[1][0] == 2.0
            
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)
    
    def test_multiple_reactions_pickle(self):
        """Test that multiple reactions can be pickled together."""
        RxnList = {}
        
        # Create multiple reactions
        for i in range(5):
            rxn_id = f"R{i:05d}"
            rxn = reaction(f"http://example.com/{rxn_id}", 0.5, rxn_id, "TestPath", 0)
            rxn.subs = [(1.0, f"C{i:05d}", f"C{i:05d}")]
            rxn.prods = [(1.0, f"C{i+10:05d}", f"C{i+10:05d}")]
            rxn.Substrate = MethodType(_rxn_substrate, rxn)
            rxn.Product = MethodType(_rxn_product, rxn)
            RxnList[rxn_id] = rxn
        
        # Pickle the entire dictionary
        with tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix='.pkl') as f:
            temp_file = f.name
            dill.dump(RxnList, f)
        
        try:
            # Load and verify
            with open(temp_file, 'rb') as f:
                loaded_list = dill.load(f)
            
            assert len(loaded_list) == 5
            
            # Test each reaction
            for i in range(5):
                rxn_id = f"R{i:05d}"
                assert rxn_id in loaded_list
                
                rxn = loaded_list[rxn_id]
                subs = rxn.Substrate()
                prods = rxn.Product()
                
                assert len(subs) == 1
                assert len(prods) == 1
                assert subs[0][1] == f"C{i:05d}"
                assert prods[0][1] == f"C{i+10:05d}"
                
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)


class TestCheckpointDataStructure:
    """Test the complete checkpoint data structure."""
    
    def test_checkpoint_full_roundtrip(self):
        """Test that complete checkpoint data can be saved and loaded."""
        # Create test data structures matching generate_db.py
        PathList = {}
        RxnList = {}
        MetList = {}
        GPRList = {}
        MetEquiv = {}
        PathNameRxn = {"TestPath": "R00001 R00002 "}
        RxnEquiv = {}
        PathIdent = ["path1", "path2", "path3"]
        RxnIdent = []
        MetIdent = []
        GPRIdent = []
        
        # Add some test reactions
        for i in range(3):
            rxn_id = f"R{i:05d}"
            rxn = reaction(f"http://example.com/{rxn_id}", 0.5, rxn_id, "TestPath", 0)
            rxn.subs = [(1.0, f"C{i:05d}", f"C{i:05d}")]
            rxn.prods = [(1.0, f"C{i+100:05d}", f"C{i+100:05d}")]
            rxn.Substrate = MethodType(_rxn_substrate, rxn)
            rxn.Product = MethodType(_rxn_product, rxn)
            RxnList[rxn_id] = rxn
            RxnIdent.append(rxn_id)
        
        # Add compartmentalized reactions
        RxnList_CL = {}
        for base_id in ["R00001", "R00002"]:
            for comp in ["c", "n", "m"]:
                rxn_id = f"{base_id}_{comp}"
                rxn = reaction(f"http://example.com/{base_id}", 0.5, rxn_id, "TestPath", 0)
                rxn.subs = [(1.0, f"C00001_{comp}", f"C00001_{comp}")]
                rxn.prods = [(1.0, f"C00002_{comp}", f"C00002_{comp}")]
                rxn.Substrate = MethodType(_rxn_substrate, rxn)
                rxn.Product = MethodType(_rxn_product, rxn)
                RxnList_CL[rxn_id] = rxn
        
        MetList_CL = {}
        RxnIdent_CL = list(RxnList_CL.keys())
        MetIdent_CL = []
        Compartment_CL = {"c": "cytosol", "n": "nucleus", "m": "mitochondrion"}
        
        # Create checkpoint data (matching generate_db.py structure)
        checkpoint_data = {
            "last_completed_pathway": 5,
            "PathList": PathList,
            "RxnList": RxnList,
            "MetList": MetList,
            "GPRList": GPRList,
            "MetEquiv": MetEquiv,
            "PathNameRxn": PathNameRxn,
            "RxnEquiv": RxnEquiv,
            "PathIdent": PathIdent,
            "RxnIdent": RxnIdent,
            "MetIdent": MetIdent,
            "GPRIdent": GPRIdent,
            "RxnList_CL": RxnList_CL,
            "MetList_CL": MetList_CL,
            "RxnIdent_CL": RxnIdent_CL,
            "MetIdent_CL": MetIdent_CL,
            "Compartment_CL": Compartment_CL,
            "failed_location_reactions": set(),
        }
        
        # Save checkpoint
        with tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix='.pkl') as f:
            temp_file = f.name
            dill.dump(checkpoint_data, f)
        
        try:
            # Load checkpoint (simulating resume)
            with open(temp_file, 'rb') as f:
                loaded_data = dill.load(f)
            
            # Verify all keys present
            required_keys = [
                "last_completed_pathway", "PathList", "RxnList", "MetList",
                "GPRList", "MetEquiv", "PathNameRxn", "RxnEquiv",
                "PathIdent", "RxnIdent", "MetIdent", "GPRIdent",
                "RxnList_CL", "MetList_CL", "RxnIdent_CL", "MetIdent_CL",
                "Compartment_CL", "failed_location_reactions"
            ]
            for key in required_keys:
                assert key in loaded_data, f"Missing key: {key}"
            
            # Verify scalar values
            assert loaded_data["last_completed_pathway"] == 5
            assert loaded_data["PathNameRxn"] == {"TestPath": "R00001 R00002 "}
            assert len(loaded_data["PathIdent"]) == 3
            
            # Verify reactions
            assert len(loaded_data["RxnList"]) == 3
            assert "R00000" in loaded_data["RxnList"]
            
            # Test reaction methods work after reload
            loaded_rxn = loaded_data["RxnList"]["R00000"]
            subs = loaded_rxn.Substrate()
            prods = loaded_rxn.Product()
            assert len(subs) == 1
            assert subs[0][1] == "C00000"
            assert len(prods) == 1
            assert prods[0][1] == "C00100"
            
            # Verify compartmentalized reactions
            assert len(loaded_data["RxnList_CL"]) == 6  # 2 reactions × 3 compartments
            assert "R00001_c" in loaded_data["RxnList_CL"]
            assert "R00001_n" in loaded_data["RxnList_CL"]
            assert "R00002_m" in loaded_data["RxnList_CL"]
            
            # Test compartmentalized reaction methods
            comp_rxn = loaded_data["RxnList_CL"]["R00001_c"]
            comp_subs = comp_rxn.Substrate()
            assert comp_subs[0][1] == "C00001_c"
            
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)
    
    def test_checkpoint_pathway_index_tracking(self):
        """Test that pathway index is correctly saved and loaded."""
        for test_index in [0, 5, 42, 100]:
            checkpoint_data = {
                "last_completed_pathway": test_index,
                "RxnList": {},
                "MetList": {},
            }
            
            with tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix='.pkl') as f:
                temp_file = f.name
                dill.dump(checkpoint_data, f)
            
            try:
                with open(temp_file, 'rb') as f:
                    loaded = dill.load(f)
                
                assert loaded["last_completed_pathway"] == test_index
                
                # Simulate resume logic from generate_db.py
                start_pathway_index = loaded.get("last_completed_pathway", 0) + 1
                assert start_pathway_index == test_index + 1
                
            finally:
                if os.path.exists(temp_file):
                    os.unlink(temp_file)


class TestSanitizeLoadedReactions:
    """Test the sanitization function for backward compatibility."""
    
    def test_sanitize_adds_missing_methods(self):
        """Test that sanitization adds missing Substrate/Product methods."""
        # Create a reaction without proper method binding (simulating old checkpoint)
        rxn = reaction("http://example.com/R00001", 0.5, "R00001", "TestPath", 0)
        rxn.subs = [(1.0, "C00001", "C00001")]
        rxn.prods = [(1.0, "C00002", "C00002")]
        # Don't bind methods - simulate old checkpoint
        
        RxnList = {"R00001": rxn}
        
        # Sanitize (should add missing methods)
        sanitize_loaded_reactions(RxnList, name="RxnList")
        
        # Verify methods were added
        assert hasattr(RxnList["R00001"], "Substrate")
        assert callable(RxnList["R00001"].Substrate)
        assert hasattr(RxnList["R00001"], "Product")
        assert callable(RxnList["R00001"].Product)
        
        # Verify they work
        subs = RxnList["R00001"].Substrate()
        prods = RxnList["R00001"].Product()
        assert len(subs) == 1
        assert len(prods) == 1
    
    def test_sanitize_preserves_existing_methods(self):
        """Test that sanitization doesn't break existing methods."""
        # Create properly bound reaction
        rxn = reaction("http://example.com/R00001", 0.5, "R00001", "TestPath", 0)
        rxn.subs = [(1.0, "C00001", "C00001")]
        rxn.prods = [(1.0, "C00002", "C00002")]
        rxn.Substrate = MethodType(_rxn_substrate, rxn)
        rxn.Product = MethodType(_rxn_product, rxn)
        
        RxnList = {"R00001": rxn}
        
        # Get data before sanitization
        subs_before = rxn.Substrate()
        prods_before = rxn.Product()
        
        # Sanitize
        sanitize_loaded_reactions(RxnList, name="RxnList")
        
        # Verify data unchanged
        subs_after = RxnList["R00001"].Substrate()
        prods_after = RxnList["R00001"].Product()
        
        assert subs_after == subs_before
        assert prods_after == prods_before
    
    def test_sanitize_handles_empty_dict(self):
        """Test that sanitization handles empty reaction dict."""
        RxnList = {}
        # Should not raise exception
        sanitize_loaded_reactions(RxnList, name="RxnList")
        assert RxnList == {}
    
    def test_sanitize_handles_multiple_reactions(self):
        """Test sanitization works on dict with multiple reactions."""
        RxnList = {}
        
        # Create several reactions without proper binding
        for i in range(10):
            rxn_id = f"R{i:05d}"
            rxn = reaction(f"http://example.com/{rxn_id}", 0.5, rxn_id, "Test", 0)
            rxn.subs = [(1.0, f"C{i:05d}", f"C{i:05d}")]
            rxn.prods = [(1.0, f"C{i+10:05d}", f"C{i+10:05d}")]
            RxnList[rxn_id] = rxn
        
        # Sanitize all at once
        sanitize_loaded_reactions(RxnList, name="RxnList")
        
        # Verify all have working methods
        for i in range(10):
            rxn_id = f"R{i:05d}"
            subs = RxnList[rxn_id].Substrate()
            prods = RxnList[rxn_id].Product()
            assert len(subs) == 1
            assert len(prods) == 1


class TestCompartmentalizedReactions:
    """Test compartmentalized reactions in checkpoints."""
    
    def test_compartmentalized_reactions_preserve_suffix(self):
        """Test that _c, _n, _m suffixes are preserved."""
        RxnList_CL = {}
        
        # Create reactions with different compartment suffixes
        compartments = ["c", "n", "m", "er", "x", "r", "v", "g"]
        for comp in compartments:
            rxn_id = f"R00001_{comp}"
            rxn = reaction("http://example.com/R00001", 0.5, rxn_id, "Test", 0)
            rxn.subs = [(1.0, f"C00001_{comp}", f"C00001_{comp}")]
            rxn.prods = [(1.0, f"C00002_{comp}", f"C00002_{comp}")]
            rxn.Substrate = MethodType(_rxn_substrate, rxn)
            rxn.Product = MethodType(_rxn_product, rxn)
            RxnList_CL[rxn_id] = rxn
        
        # Pickle and unpickle
        with tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix='.pkl') as f:
            temp_file = f.name
            dill.dump(RxnList_CL, f)
        
        try:
            with open(temp_file, 'rb') as f:
                loaded = dill.load(f)
            
            # Verify all compartments preserved
            assert len(loaded) == len(compartments)
            for comp in compartments:
                rxn_id = f"R00001_{comp}"
                assert rxn_id in loaded, f"Missing {rxn_id}"
                
                # Verify metabolite IDs have correct compartment
                rxn = loaded[rxn_id]
                subs = rxn.Substrate()
                prods = rxn.Product()
                assert subs[0][1] == f"C00001_{comp}"
                assert prods[0][1] == f"C00002_{comp}"
                
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)
    
    def test_compartmentalized_reactions_data_integrity(self):
        """Test that compartmentalized reactions maintain data integrity."""
        RxnList_CL = {}
        
        # Create same reaction in multiple compartments with different stoichiometry
        test_data = {
            "c": [(1.0, "C00001_c", "C00001_c"), (2.0, "C00002_c", "C00002_c")],
            "n": [(1.5, "C00001_n", "C00001_n")],
            "m": [(3.0, "C00001_m", "C00001_m"), (1.0, "C00002_m", "C00002_m"), (1.0, "C00003_m", "C00003_m")],
        }
        
        for comp, subs_data in test_data.items():
            rxn_id = f"R00042_{comp}"
            rxn = reaction("http://example.com/R00042", 0.5, rxn_id, "Test", 0)
            rxn.subs = subs_data
            rxn.prods = [(1.0, f"C00099_{comp}", f"C00099_{comp}")]
            rxn.Substrate = MethodType(_rxn_substrate, rxn)
            rxn.Product = MethodType(_rxn_product, rxn)
            RxnList_CL[rxn_id] = rxn
        
        # Pickle and verify
        with tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix='.pkl') as f:
            temp_file = f.name
            dill.dump(RxnList_CL, f)
        
        try:
            with open(temp_file, 'rb') as f:
                loaded = dill.load(f)
            
            # Verify exact stoichiometry and metabolite counts
            for comp, expected_subs in test_data.items():
                rxn_id = f"R00042_{comp}"
                actual_subs = loaded[rxn_id].Substrate()
                assert len(actual_subs) == len(expected_subs)
                assert actual_subs == expected_subs
                
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)


class TestFailedLocationReactions:
    """Test filtering of reactions with failed location detection."""
    
    def test_failed_location_reactions_are_filtered(self):
        """Test that reactions with failed location detection are excluded from checkpoint."""
        RxnList_CL = {
            "R00001_c": "reaction_obj_1",
            "R00001_n": "reaction_obj_1n",
            "R00002_c": "reaction_obj_2",
            "R00003_c": "reaction_obj_3",  # This one failed
            "R00003_m": "reaction_obj_3m",  # This one also failed (same base)
            "R00004_c": "reaction_obj_4",
        }
        failed_location_reactions = {"R00003"}
        
        # Filter (simulating checkpoint save logic from generate_db.py line 1385-1390)
        RxnList_CL_filtered = {
            rxn_key: rxn_obj 
            for rxn_key, rxn_obj in RxnList_CL.items()
            if not any(failed_rxn_id in rxn_key for failed_rxn_id in failed_location_reactions)
        }
        
        # Verify R00003_c and R00003_m were excluded
        assert len(RxnList_CL_filtered) == 4
        assert "R00001_c" in RxnList_CL_filtered
        assert "R00001_n" in RxnList_CL_filtered
        assert "R00002_c" in RxnList_CL_filtered
        assert "R00004_c" in RxnList_CL_filtered
        assert "R00003_c" not in RxnList_CL_filtered
        assert "R00003_m" not in RxnList_CL_filtered
    
    def test_failed_location_reactions_set_persists(self):
        """Test that failed_location_reactions set is saved and loaded."""
        failed_reactions = {"R00005", "R00123", "R00456"}
        
        checkpoint_data = {
            "last_completed_pathway": 10,
            "RxnList": {},
            "failed_location_reactions": failed_reactions,
        }
        
        with tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix='.pkl') as f:
            temp_file = f.name
            dill.dump(checkpoint_data, f)
        
        try:
            with open(temp_file, 'rb') as f:
                loaded = dill.load(f)
            
            assert "failed_location_reactions" in loaded
            assert loaded["failed_location_reactions"] == failed_reactions
            assert len(loaded["failed_location_reactions"]) == 3
            assert "R00123" in loaded["failed_location_reactions"]
            
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)


class TestOperatorAddReplacement:
    """Test that operator.add works correctly for list concatenation."""
    
    def test_operator_add_concatenates_lists(self):
        """Test operator.add works like lambda x, y: x + y."""
        import operator
        from functools import reduce
        
        lists = [[1, 2], [3, 4], [5, 6]]
        result = reduce(operator.add, lists, [])
        
        assert result == [1, 2, 3, 4, 5, 6]
    
    def test_operator_add_with_empty_lists(self):
        """Test operator.add handles empty lists correctly."""
        import operator
        from functools import reduce
        
        lists = [[], [1], [], [2, 3], []]
        result = reduce(operator.add, lists, [])
        
        assert result == [1, 2, 3]
    
    def test_operator_add_with_single_list(self):
        """Test operator.add with single list."""
        import operator
        from functools import reduce
        
        lists = [[1, 2, 3]]
        result = reduce(operator.add, lists, [])
        
        assert result == [1, 2, 3]
    
    def test_operator_add_with_nested_structures(self):
        """Test operator.add with complex nested structures."""
        import operator
        from functools import reduce
        
        # Simulating reaction query results (lists of objects)
        lists = [
            [{"id": "R1", "name": "reaction1"}],
            [{"id": "R2", "name": "reaction2"}, {"id": "R3", "name": "reaction3"}],
            [{"id": "R4", "name": "reaction4"}],
        ]
        result = reduce(operator.add, lists, [])
        
        assert len(result) == 4
        assert result[0]["id"] == "R1"
        assert result[2]["id"] == "R3"


class TestCheckpointBackwardCompatibility:
    """Test backward compatibility with old checkpoint formats."""
    
    def test_load_checkpoint_without_failed_locations(self):
        """Test loading old checkpoint that doesn't have failed_location_reactions."""
        # Create checkpoint without the new field
        checkpoint_data = {
            "last_completed_pathway": 3,
            "RxnList": {},
            "MetList": {},
            "PathIdent": ["path1"],
        }
        
        with tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix='.pkl') as f:
            temp_file = f.name
            dill.dump(checkpoint_data, f)
        
        try:
            with open(temp_file, 'rb') as f:
                loaded = dill.load(f)
            
            # Simulate generate_db.py loading logic
            failed_location_reactions = loaded.get("failed_location_reactions", set())
            
            # Should default to empty set
            assert isinstance(failed_location_reactions, set)
            assert len(failed_location_reactions) == 0
            
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v", "--tb=short"])
