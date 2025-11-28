"""
Unit tests for reaction identification module

Tests verify:
1. Jaccard similarity calculation
2. KEGG metabolite gathering
3. Metabolite ID replacement
4. Jaccard-based reaction matching
5. Process jaccard results
6. Edge cases and error handling

Run with: pytest tests/test_reaction_identification.py -v
"""

import os
import sys
import pytest
import pandas as pd
from unittest.mock import Mock, patch, mock_open
from io import StringIO

# Add project to path
project_root = os.path.join(os.path.dirname(__file__), '..', '..')
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from functions.function_reac_identification import (
    jaccard,
    gather_kegg_metabolites,
    replace_met_id_by_met_kegg,
    execute_jaccard,
    process_jaccard
)


class TestJaccardSimilarity:
    """Test Jaccard similarity calculation."""
    
    def test_identical_lists(self):
        """Test Jaccard similarity of identical lists."""
        list1 = ['A', 'B', 'C']
        list2 = ['A', 'B', 'C']
        
        similarity = jaccard(list1, list2)
        
        assert similarity == 1.0
    
    def test_disjoint_lists(self):
        """Test Jaccard similarity of completely different lists."""
        list1 = ['A', 'B', 'C']
        list2 = ['D', 'E', 'F']
        
        similarity = jaccard(list1, list2)
        
        assert similarity == 0.0
    
    def test_partial_overlap(self):
        """Test Jaccard similarity with partial overlap."""
        list1 = ['A', 'B', 'C', 'D']
        list2 = ['C', 'D', 'E', 'F']
        
        similarity = jaccard(list1, list2)
        
        # Intersection: {C, D} = 2
        # Union: len(list1) + len(list2) - intersection = 4 + 4 - 2 = 6
        # Jaccard = 2/6 = 0.333...
        assert abs(similarity - 0.333) < 0.01
    
    def test_empty_intersection(self):
        """Test Jaccard similarity with no overlap."""
        list1 = ['X']
        list2 = ['A', 'B']
        
        # Intersection = 0, Union = 3, Jaccard = 0/3 = 0
        similarity = jaccard(list1, list2)
        assert similarity == 0.0
    
    def test_single_element_match(self):
        """Test Jaccard with single matching element."""
        list1 = ['A']
        list2 = ['A']
        
        similarity = jaccard(list1, list2)
        
        assert similarity == 1.0
    
    def test_duplicates_in_lists(self):
        """Test that duplicates affect calculation (uses list length, not sets)."""
        list1 = ['A', 'A', 'B']
        list2 = ['A', 'B', 'B']
        
        similarity = jaccard(list1, list2)
        
        # Intersection: {A, B} = 2 items
        # Union: len(list1) + len(list2) - intersection = 3 + 3 - 2 = 4
        # Jaccard = 2/4 = 0.5
        assert similarity == 0.5


class TestGatherKeggMetabolites:
    """Test gathering KEGG metabolites from SBML."""
    
    @patch('builtins.open', new_callable=mock_open, read_data='''
        <species metaid="meta_M_glc_D_c">
            <rdf:Description rdf:about="#M_glc_D_c">
                <bqbiol:is rdf:resource="https://identifiers.org/kegg.compound/C00031"/>
            </rdf:Description>
        </species>
        <species metaid="meta_M_atp_c">
            <rdf:Description rdf:about="#M_atp_c">
                <bqbiol:is rdf:resource="https://identifiers.org/kegg.compound/C00002"/>
            </rdf:Description>
        </species>
    ''')
    def test_gather_basic_metabolites(self, mock_file):
        """Test gathering basic KEGG metabolites."""
        result = gather_kegg_metabolites('dummy_file.xml')
        
        assert isinstance(result, dict)
        assert 'M_glc_D_c' in result
        assert 'M_atp_c' in result
        assert 'C00031' in result['M_glc_D_c']
        assert 'C00002' in result['M_atp_c']
    
    @patch('builtins.open', new_callable=mock_open, read_data='<root></root>')
    def test_gather_no_metabolites(self, mock_file):
        """Test when no metabolites are found."""
        result = gather_kegg_metabolites('empty_file.xml')
        
        assert result == {}
    
    @patch('builtins.open', new_callable=mock_open, read_data='''
        <species metaid="meta_M_glc_D_c">
            <rdf:Description rdf:about="#M_glc_D_c">
                <bqbiol:is rdf:resource="https://identifiers.org/kegg.compound/C00031"/>
            </rdf:Description>
        </species>
        <species metaid="meta_M_glc_D_c">
            <rdf:Description rdf:about="#M_glc_D_c">
                <bqbiol:is rdf:resource="https://identifiers.org/kegg.compound/C00031"/>
            </rdf:Description>
        </species>
    ''')
    def test_gather_removes_duplicates(self, mock_file):
        """Test that duplicate metabolites are handled."""
        result = gather_kegg_metabolites('duplicate_file.xml')
        
        # Should only have one entry
        assert len([k for k in result.keys() if 'glc_D' in k]) == 1


class TestReplaceMetIdByKegg:
    """Test metabolite ID replacement."""
    
    def test_replace_basic_ids(self):
        """Test replacing metabolite IDs."""
        df = pd.DataFrame({
            0: ['R1'],
            5: ['M_glc_D_c,M_atp_c'],
            6: ['M_g6p_c,M_adp_c']
        })
        
        dictionary = {
            'M_glc_D_c': 'C00031c',
            'M_atp_c': 'C00002c',
            'M_g6p_c': 'C00092c',
            'M_adp_c': 'C00008c'
        }
        
        result = replace_met_id_by_met_kegg(df, dictionary)
        
        assert 'C00031c' in result[5].iloc[0]
        assert 'C00002c' in result[5].iloc[0]
    
    def test_replace_filters_na(self):
        """Test that NaN values are filtered out."""
        df = pd.DataFrame({
            0: ['R1', 'R2'],
            5: ['M_glc_D_c', None],
            6: ['M_g6p_c', 'M_adp_c']
        })
        
        dictionary = {'M_glc_D_c': 'C00031c'}
        
        result = replace_met_id_by_met_kegg(df, dictionary)
        
        # Should only have one row (NaN filtered)
        assert len(result) == 1


class TestExecuteJaccard:
    """Test Jaccard-based reaction matching."""
    
    def test_exact_match_forward(self):
        """Test exact match in forward direction."""
        b = pd.DataFrame({
            0: ['R1'],
            5: ['C00031c,C00002c'],
            6: ['C00092c,C00008c']
        })
        
        c2 = pd.DataFrame({
            1: ['MAR1'],
            5: ['C00031c,C00002c'],
            6: ['C00092c,C00008c']
        })
        
        result = execute_jaccard(b, c2)
        
        assert len(result) == 1
        assert result[6].iloc[0] == 2.0  # Perfect match
    
    def test_exact_match_reverse(self):
        """Test exact match in reverse direction."""
        b = pd.DataFrame({
            0: ['R1'],
            5: ['C00031c,C00002c'],
            6: ['C00092c,C00008c']
        })
        
        c2 = pd.DataFrame({
            1: ['MAR1'],
            5: ['C00092c,C00008c'],  # Reversed
            6: ['C00031c,C00002c']   # Reversed
        })
        
        result = execute_jaccard(b, c2)
        
        assert len(result) == 1
        assert result[6].iloc[0] == 2.0  # Perfect match (reversible)
    
    def test_no_match(self):
        """Test when reactions don't match."""
        b = pd.DataFrame({
            0: ['R1'],
            5: ['C00031c'],
            6: ['C00092c']
        })
        
        c2 = pd.DataFrame({
            1: ['MAR1'],
            5: ['C00001c'],
            6: ['C00002c']
        })
        
        result = execute_jaccard(b, c2)
        
        # Should return empty DataFrame with correct columns
        assert len(result) == 0
        assert len(result.columns) == 7
    
    def test_multiple_reactions(self):
        """Test matching with multiple reactions."""
        b = pd.DataFrame({
            0: ['R1', 'R2'],
            5: ['C00031c', 'C00001c'],
            6: ['C00092c', 'C00002c']
        })
        
        c2 = pd.DataFrame({
            1: ['MAR1', 'MAR2'],
            5: ['C00031c', 'C00001c'],
            6: ['C00092c', 'C00002c']
        })
        
        result = execute_jaccard(b, c2)
        
        # Should find 2 matches
        assert len(result) == 2


class TestProcessJaccard:
    """Test Jaccard result processing."""
    
    def test_process_basic_results(self):
        """Test basic processing of Jaccard results."""
        left_model = pd.DataFrame({
            0: ['R1'],
            1: ['R00001'],
            2: ['c1'],
            3: ['c1,c2'],
            4: ['true'],
            5: ['C00031c'],
            6: ['C00092c'],
            7: ['1'],
            8: ['1,1'],
            9: ['2.7.1.1'],
            10: ['gene1'],
            11: ['MAR1']
        })
        
        right_model = pd.DataFrame({
            13: [0],
            14: ['MAR1'],
            15: ['c1'],
            16: ['c1,c2'],
            17: ['true'],
            18: ['C00031c'],
            19: ['C00092c'],
            20: ['1'],
            21: ['1,1'],
            22: ['2.7.1.1'],
            23: ['gene1'],
            24: ['MAR00001']
        })
        
        jaccard_result = pd.DataFrame({
            0: ['R1'],
            5: ['C00031c'],
            6: ['C00092c'],
            14: ['MAR1'],
            18: ['C00031c'],
            19: ['C00092c'],
            12: [2.0]
        })
        
        result = process_jaccard(left_model, right_model, jaccard_result)
        
        assert isinstance(result, dict)
        assert len(result) > 0


class TestEdgeCases:
    """Test edge cases and error handling."""
    
    def test_jaccard_with_single_elements(self):
        """Test Jaccard with minimal lists."""
        similarity = jaccard(['A'], ['B'])
        
        assert similarity == 0.0
    
    def test_execute_jaccard_empty_dataframes(self):
        """Test Jaccard with empty DataFrames."""
        b = pd.DataFrame({0: [], 5: [], 6: []})
        c2 = pd.DataFrame({1: [], 5: [], 6: []})
        
        result = execute_jaccard(b, c2)
        
        assert len(result) == 0
    
    def test_jaccard_order_independence(self):
        """Test that Jaccard is order-independent."""
        list1 = ['A', 'B', 'C']
        list2_ordered = ['A', 'B', 'C']
        list2_shuffled = ['C', 'A', 'B']
        
        sim1 = jaccard(list1, list2_ordered)
        sim2 = jaccard(list1, list2_shuffled)
        
        assert sim1 == sim2 == 1.0


class TestReactionComparison:
    """Test reaction comparison logic."""
    
    def test_bidirectional_matching(self):
        """Test that reactions match in both directions."""
        # A + B -> C + D should match C + D -> A + B (reverse)
        b = pd.DataFrame({
            0: ['R1'],
            5: ['A,B'],
            6: ['C,D']
        })
        
        c2_forward = pd.DataFrame({
            1: ['MAR1'],
            5: ['A,B'],
            6: ['C,D']
        })
        
        c2_reverse = pd.DataFrame({
            1: ['MAR2'],
            5: ['C,D'],
            6: ['A,B']
        })
        
        result_forward = execute_jaccard(b, c2_forward)
        result_reverse = execute_jaccard(b, c2_reverse)
        
        assert len(result_forward) == 1
        assert len(result_reverse) == 1


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
