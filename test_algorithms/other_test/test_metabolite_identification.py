"""
Unit tests for metabolite identification module

Tests verify:
1. Formula similarity computation
2. PubChem API interaction (mocked)
3. Metabolite identification logic
4. Formula parsing (atom function)
5. Annotation processing
6. Rate limiting
7. Retry logic
8. Error handling

Run with: pytest tests/test_metabolite_identification.py -v
"""

import os
import sys
import pytest
from unittest.mock import Mock, patch, MagicMock, mock_open
import tempfile
import time

# Add project to path
project_root = os.path.join(os.path.dirname(__file__), '..', '..')
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from functions.function_metabolite_identification import (
    formula_similarity,
    identify_metabolite,
    atom,
    gather_metabolites,
    remove_null_value,
    _rate_limit,
    setup_proxy
)


class TestFormulaSimilarity:
    """Test formula similarity calculation."""
    
    def test_identical_formulas(self):
        """Test similarity of identical formulas."""
        formula1 = "C6H12O6"
        formula2 = "C6H12O6"
        
        similarity = formula_similarity(formula1, formula2)
        
        assert similarity == 1.0
    
    def test_similar_formulas(self):
        """Test similarity of similar formulas."""
        formula1 = "C6H12O6"
        formula2 = "C6H14O6"
        
        similarity = formula_similarity(formula1, formula2)
        
        # After H normalization, formulas become very similar
        # (H12 and H14 both become H)
        assert similarity >= 0.9
    
    def test_different_formulas(self):
        """Test similarity of different formulas."""
        formula1 = "C6H12O6"  # glucose
        formula2 = "C2H5OH"   # ethanol
        
        similarity = formula_similarity(formula1, formula2)
        
        # After normalization (H becomes H), similarity is moderate
        # Original: C6HO6 vs C2HO
        assert 0.5 <= similarity < 0.8
    
    def test_formula_with_charge(self):
        """Test that charges are stripped before comparison."""
        formula1 = "C6H12O6-2"
        formula2 = "C6H12O6"
        
        similarity = formula_similarity(formula1, formula2)
        
        # Charges should be ignored
        assert similarity == 1.0
    
    def test_formula_with_hydrogen_numbers(self):
        """Test that hydrogen numbers are normalized."""
        formula1 = "C6H12O6"
        formula2 = "C6HO6"  # H without number
        
        # Should handle H normalization
        similarity = formula_similarity(formula1, formula2)
        assert isinstance(similarity, float)


class TestAtomParsing:
    """Test chemical formula atom parsing."""
    
    def test_simple_formula_carbon_hydrogen(self):
        """Test parsing simple C and H formula."""
        result = atom("C6H12")
        
        assert "C6" in result
        assert "H12" in result
    
    def test_formula_with_oxygen(self):
        """Test formula with oxygen."""
        result = atom("C6H12O6")
        
        assert "C6" in result
        assert "H12" in result
        assert "O6" in result
    
    def test_formula_with_nitrogen(self):
        """Test formula with nitrogen."""
        result = atom("C2H5NO2")
        
        assert "N" in result
    
    def test_formula_with_phosphorus(self):
        """Test formula with phosphorus."""
        result = atom("C10H16N5O13P3")
        
        assert "P3" in result
    
    def test_formula_with_sulfur(self):
        """Test formula with sulfur."""
        result = atom("C3H7NO2S")
        
        assert "S" in result
    
    def test_single_carbon(self):
        """Test that single C is represented as C1."""
        result = atom("CH4")
        
        assert "C" in result
        assert "H" in result
    
    def test_empty_formula(self):
        """Test handling of empty formula."""
        result = atom("")
        
        # Should return empty string or handle gracefully
        assert result == "" or result is None


class TestIdentifyMetabolite:
    """Test metabolite identification with PubChem."""
    
    @patch('functions.function_metabolite_identification.pcp')
    @patch('functions.function_metabolite_identification.time.sleep')
    def test_successful_identification(self, mock_sleep, mock_pcp):
        """Test successful metabolite identification."""
        # Mock PubChem responses
        mock_pcp.get_cids.return_value = [5793]  # Glucose CID
        
        mock_compound = Mock()
        mock_compound.cid = 5793
        mock_compound.molecular_formula = "C6H12O6"
        mock_compound.synonyms = ['Glucose', 'D-Glucose', 'Dextrose', 'C00031']
        mock_compound.inchi = "InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1"
        mock_compound.inchikey = "WQZGKKKJIJFFOK-GASJEMHNSA-N"
        
        mock_pcp.Compound.from_cid.return_value = mock_compound
        
        result = identify_metabolite("Glucose", "C6H12O6", "glc_D")
        
        assert result is not None
        assert "Glucose" in result or "Dextrose" in result
        assert "C6H12O6" in result
        assert "5793" in result
    
    @patch('functions.function_metabolite_identification.pcp')
    @patch('functions.function_metabolite_identification.time.sleep')
    def test_formula_mismatch_rejection(self, mock_sleep, mock_pcp):
        """Test that mismatched formulas are rejected."""
        # Mock PubChem to return wrong formula
        mock_pcp.get_cids.return_value = [123]
        
        mock_compound = Mock()
        mock_compound.molecular_formula = "C2H5OH"  # Wrong formula
        mock_pcp.Compound.from_cid.return_value = mock_compound
        
        result = identify_metabolite("Glucose", "C6H12O6", "glc_D", threshold=0.9)
        
        # Should reject due to formula mismatch
        assert result is None
    
    @patch('functions.function_metabolite_identification.pcp')
    @patch('functions.function_metabolite_identification.time.sleep')
    def test_not_found_in_pubchem(self, mock_sleep, mock_pcp):
        """Test handling when compound not found."""
        mock_pcp.get_cids.return_value = None
        
        result = identify_metabolite("UnknownCompound", "C100H200", "unk")
        
        assert result is None
    
    @patch('functions.function_metabolite_identification.pcp')
    @patch('functions.function_metabolite_identification.time.sleep')
    def test_retry_on_server_busy(self, mock_sleep, mock_pcp):
        """Test retry logic when PubChem server is busy."""
        # First call raises 503 error, second succeeds
        mock_pcp.get_cids.side_effect = [
            Exception("503 Server Busy"),
            [5793]
        ]
        
        mock_compound = Mock()
        mock_compound.cid = 5793
        mock_compound.molecular_formula = "C6H12O6"
        mock_compound.synonyms = ['Glucose']
        mock_compound.inchi = "test_inchi"
        mock_compound.inchikey = "test_key"
        mock_pcp.Compound.from_cid.return_value = mock_compound
        
        result = identify_metabolite("Glucose", "C6H12O6", "glc_D", max_retries=3)
        
        # Should succeed after retry
        assert mock_pcp.get_cids.call_count == 2
    
    @patch('functions.function_metabolite_identification.pcp')
    @patch('functions.function_metabolite_identification.time.sleep')
    def test_max_retries_exceeded(self, mock_sleep, mock_pcp):
        """Test that max retries limit is respected."""
        # Always raise 503 error
        mock_pcp.get_cids.side_effect = Exception("503 ServerBusy")
        
        result = identify_metabolite("Glucose", "C6H12O6", "glc_D", max_retries=2)
        
        assert result is None
        assert mock_pcp.get_cids.call_count == 2
    
    @patch('functions.function_metabolite_identification.pcp')
    @patch('functions.function_metabolite_identification.time.sleep')
    def test_kegg_id_extraction(self, mock_sleep, mock_pcp):
        """Test extraction of KEGG ID from synonyms."""
        mock_pcp.get_cids.return_value = [5793]
        
        mock_compound = Mock()
        mock_compound.cid = 5793
        mock_compound.molecular_formula = "C6H12O6"
        mock_compound.synonyms = ['Glucose', 'C00031', 'D-Glucose']  # KEGG ID
        mock_compound.inchi = "test"
        mock_compound.inchikey = "test"
        mock_pcp.Compound.from_cid.return_value = mock_compound
        
        result = identify_metabolite("Glucose", "C6H12O6", "glc_D")
        
        assert result is not None
        assert "C00031" in result
    
    @patch('functions.function_metabolite_identification.pcp')
    @patch('functions.function_metabolite_identification.time.sleep')
    def test_chebi_id_extraction(self, mock_sleep, mock_pcp):
        """Test extraction of ChEBI ID from synonyms."""
        mock_pcp.get_cids.return_value = [5793]
        
        mock_compound = Mock()
        mock_compound.cid = 5793
        mock_compound.molecular_formula = "C6H12O6"
        mock_compound.synonyms = ['Glucose', 'CHEBI:17234']
        mock_compound.inchi = "test"
        mock_compound.inchikey = "test"
        mock_pcp.Compound.from_cid.return_value = mock_compound
        
        result = identify_metabolite("Glucose", "C6H12O6", "glc_D")
        
        assert result is not None
        assert "CHEBI:17234" in result
    
    @patch('functions.function_metabolite_identification.pcp')
    @patch('functions.function_metabolite_identification.time.sleep')
    def test_lipidmaps_id_extraction(self, mock_sleep, mock_pcp):
        """Test extraction of LIPIDMAPS ID from synonyms."""
        mock_pcp.get_cids.return_value = [123]
        
        mock_compound = Mock()
        mock_compound.cid = 123
        mock_compound.molecular_formula = "C18H32O2"
        mock_compound.synonyms = ['Linoleic acid', 'LMFA01030120']
        mock_compound.inchi = "test"
        mock_compound.inchikey = "test"
        mock_pcp.Compound.from_cid.return_value = mock_compound
        
        result = identify_metabolite("Linoleic acid", "C18H32O2", "lnl")
        
        assert result is not None
        assert "LMFA01030120" in result


class TestGatherMetabolites:
    """Test gathering metabolites from COBRA model."""
    
    @patch('cobra.Model')
    def test_gather_unique_metabolites(self, mock_model):
        """Test gathering unique metabolites from model."""
        # Create mock metabolites
        met1 = Mock()
        met1.name = "Glucose"
        met1.formula = "C6H12O6"
        met1.annotation = {}
        met1.id = "glc_D_c"
        
        met2 = Mock()
        met2.name = "ATP"
        met2.formula = "C10H16N5O13P3"
        met2.annotation = {}
        met2.id = "atp_c"
        
        mock_model.metabolites = [met1, met2]
        
        result = gather_metabolites(mock_model)
        
        assert len(result) == 2
        assert ("Glucose", "C6H12O6", {}, "glc_D_") in result
        assert ("ATP", "C10H16N5O13P3", {}, "atp_") in result
    
    @patch('cobra.Model')
    def test_gather_removes_duplicates(self, mock_model):
        """Test that duplicate metabolites are removed."""
        met1 = Mock()
        met1.name = "Glucose"
        met1.formula = "C6H12O6"
        met1.annotation = {}
        met1.id = "glc_D_c"
        
        met2 = Mock()
        met2.name = "Glucose"
        met2.formula = "C6H12O6"
        met2.annotation = {}
        met2.id = "glc_D_e"  # Different compartment
        
        mock_model.metabolites = [met1, met2]
        
        result = gather_metabolites(mock_model)
        
        # Should only have one unique entry
        assert len(result) == 1


class TestRemoveNullValue:
    """Test null value removal from dictionaries."""
    
    def test_remove_empty_dict(self):
        """Test removing empty nested dictionaries."""
        input_dict = {
            'a': {'b': {}},
            'c': {'d': 'value'}
        }
        
        result = remove_null_value(input_dict)
        
        assert 'a' not in result  # Empty nested dict removed
        assert 'c' in result
        assert result['c'] == {'d': 'value'}
    
    def test_remove_deeply_nested_empty(self):
        """Test removing deeply nested empty values."""
        input_dict = {
            'a': {'b': {'c': {}}},
            'd': {'e': {'f': 'value'}}
        }
        
        result = remove_null_value(input_dict)
        
        assert 'a' not in result
        assert 'd' in result
    
    def test_keep_non_empty_values(self):
        """Test that non-empty values are kept."""
        input_dict = {
            'a': 'value1',
            'b': {'c': 'value2'}
        }
        
        result = remove_null_value(input_dict)
        
        assert result == input_dict


class TestRateLimiting:
    """Test rate limiting functionality."""
    
    @patch('functions.function_metabolite_identification.time')
    def test_rate_limit_enforces_delay(self, mock_time):
        """Test that rate limiting enforces minimum delay."""
        mock_time.time.side_effect = [0, 0.05, 0.1]  # Simulate fast requests
        mock_time.sleep = Mock()
        
        _rate_limit()
        
        # Should have called sleep to enforce delay
        assert mock_time.sleep.called


class TestProxySetup:
    """Test proxy configuration."""
    
    def test_proxy_disabled_by_env_var(self, monkeypatch):
        """Test proxy disabled via environment variable."""
        monkeypatch.setenv("DISABLE_PROXY", "1")
        
        result = setup_proxy()
        
        assert result is None
    
    def test_proxy_from_env_vars(self, monkeypatch):
        """Test proxy configuration from environment variables."""
        monkeypatch.setenv("HTTP_PROXY", "http://proxy:8080")
        monkeypatch.setenv("HTTPS_PROXY", "http://proxy:8080")
        
        result = setup_proxy()
        
        assert result is not None
        assert result['http'] == "http://proxy:8080"
        assert result['https'] == "http://proxy:8080"
    
    def test_no_proxy_configured(self, monkeypatch):
        """Test when no proxy is configured."""
        monkeypatch.delenv("DISABLE_PROXY", raising=False)
        monkeypatch.delenv("HTTP_PROXY", raising=False)
        monkeypatch.delenv("HTTPS_PROXY", raising=False)
        
        result = setup_proxy()
        
        assert result is None


class TestAnnotationProcessing:
    """Test annotation file processing."""
    
    @patch('functions.function_metabolite_identification.pd.read_csv')
    def test_process_annotation_file(self, mock_read_csv):
        """Test processing annotation file."""
        # Create mock DataFrame with proper indexing support
        import pandas as pd
        
        # Create a real DataFrame with scalar values (not lists)
        mock_data = {
            0: 'Glucose',
            1: 'Synonym',
            2: 'C6H12O6',
            3: 'LMFA01',
            4: 'C00031',
            5: 'CHEBI:17234',
            6: '5793',
            7: '',
            8: 'WQZGKKKJIJFFOK',
            9: 'InChI=test',
            10: '',
            11: 'glc_D'
        }
        mock_df = pd.DataFrame([mock_data])
        mock_read_csv.return_value = mock_df
        
        # This should not raise an exception
        from functions.function_metabolite_identification import process_annotation
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("test\tdata\n")
            temp_file = f.name
        
        try:
            result = process_annotation(temp_file)
            assert isinstance(result, dict)
            # Should have one entry for the metabolite
            assert 'glc_D' in result
            assert 'kegg.compound' in result['glc_D']
        finally:
            os.unlink(temp_file)


class TestEdgeCases:
    """Test edge cases and error conditions."""
    
    def test_formula_similarity_with_empty_strings(self):
        """Test formula similarity with empty strings."""
        similarity = formula_similarity("", "")
        
        assert isinstance(similarity, float)
    
    def test_atom_with_complex_formula(self):
        """Test atom parsing with complex formula."""
        result = atom("C21H27N7O14P2")
        
        assert "C21" in result
        assert "N7" in result
        assert "P2" in result
    
    @patch('functions.function_metabolite_identification.pcp')
    @patch('functions.function_metabolite_identification.time.sleep')
    def test_identify_metabolite_with_parentheses_in_name(self, mock_sleep, mock_pcp):
        """Test identification with parentheses in metabolite name."""
        # First call fails, second (without parens) succeeds
        mock_pcp.get_cids.side_effect = [None, [123]]
        
        mock_compound = Mock()
        mock_compound.cid = 123
        mock_compound.molecular_formula = "C6H12O6"
        mock_compound.synonyms = ['Test']
        mock_compound.inchi = "test"
        mock_compound.inchikey = "test"
        mock_pcp.Compound.from_cid.return_value = mock_compound
        
        result = identify_metabolite("Glucose (D)", "C6H12O6", "glc_D")
        
        # Should try without parentheses
        assert mock_pcp.get_cids.call_count == 2


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
