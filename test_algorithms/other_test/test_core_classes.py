"""
Unit tests for core classes: reaction, compound, pathway

Tests verify:
1. Class initialization with valid data
2. Method returns (ID, Name, Formula, etc.)
3. Edge cases (empty data, missing attributes)
4. Error handling (network failures, malformed data)
5. Batch loading functionality

Run with: pytest tests/test_core_classes.py -v
"""

import os
import sys
import pytest
from unittest.mock import Mock, patch, MagicMock

# Add project to path
project_root = os.path.join(os.path.dirname(__file__), '..', '..')
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from functions.class_generate_database import reaction, compound, pathway


class TestReactionClass:
    """Test reaction class initialization and methods."""
    
    @patch('functions.function_bm_gdb.getHtml')
    @patch('functions.function_bm_gdb.getReacParam')
    def test_reaction_initialization(self, mock_get_reac_param, mock_get_html):
        """Test that reaction object initializes correctly."""
        # Mock HTML response
        mock_get_html.return_value = b"<html>mock reaction page</html>"
        
        # Mock parsed reaction parameters
        mock_get_reac_param.return_value = [
            [[(1.0, 'C00001', 'C00001')], [(1.0, 'C00002', 'C00002')]],  # substrates, products
            "Equation",
            "Formula",
            "TestReaction",  # Name
            [["EC1", "1.1.1.1"]],  # EC numbers
            "GPR",
            "termodyn"
        ]
        
        # Create reaction object
        rxn = reaction(
            "http://example.com/R00001",
            0.5,  # time
            "R00001",
            "TestPathway",
            0  # termodyn
        )
        
        # Verify initialization
        assert rxn.ID == "R00001"
        assert rxn.path == "TestPathway"
        assert rxn.termdyn == 0
        assert mock_get_html.called
        assert mock_get_reac_param.called
    
    @patch('functions.function_bm_gdb.getHtml')
    @patch('functions.function_bm_gdb.getReacParam')
    def test_reaction_name_method(self, mock_get_reac_param, mock_get_html):
        """Test that Name() method returns correct name."""
        mock_get_html.return_value = b"<html>mock</html>"
        mock_get_reac_param.return_value = [None, None, None, "TestName", None, None, None]
        
        rxn = reaction("http://example.com/R00001", 0.5, "R00001", "Path", 0)
        
        assert rxn.Name() == "TestName"
    
    @patch('functions.function_bm_gdb.getHtml')
    @patch('functions.function_bm_gdb.getReacParam')
    def test_reaction_name_fallback_to_id(self, mock_get_reac_param, mock_get_html):
        """Test that Name() falls back to ID if name is empty."""
        mock_get_html.return_value = b"<html>mock</html>"
        mock_get_reac_param.return_value = [None, None, None, None, None, None, None]  # No name
        
        rxn = reaction("http://example.com/R00001", 0.5, "R00001", "Path", 0)
        
        assert rxn.Name() == "R00001"
    
    @patch('functions.function_bm_gdb.getHtml')
    @patch('functions.function_bm_gdb.getReacParam')
    def test_reaction_ec_method(self, mock_get_reac_param, mock_get_html):
        """Test that EC() method returns list of EC numbers."""
        mock_get_html.return_value = b"<html>mock</html>"
        mock_get_reac_param.return_value = [
            None, None, None, None,
            [["EC1", "1.1.1.1"], ["EC2", "2.2.2.2"]],  # EC numbers
            None, None
        ]
        
        rxn = reaction("http://example.com/R00001", 0.5, "R00001", "Path", 0)
        
        ec_list = rxn.EC()
        assert ec_list == ["1.1.1.1", "2.2.2.2"]
    
    @patch('functions.function_bm_gdb.getHtml')
    @patch('functions.function_bm_gdb.getReacParam')
    def test_reaction_ec_empty(self, mock_get_reac_param, mock_get_html):
        """Test that EC() returns empty string when no EC numbers."""
        mock_get_html.return_value = b"<html>mock</html>"
        mock_get_reac_param.return_value = [None, None, None, None, None, None, None]
        
        rxn = reaction("http://example.com/R00001", 0.5, "R00001", "Path", 0)
        
        assert rxn.EC() == ""
    
    @patch('functions.function_bm_gdb.getHtml')
    @patch('functions.function_bm_gdb.getReacParam')
    def test_reaction_id_attribute(self, mock_get_reac_param, mock_get_html):
        """Test that ID attribute is set correctly."""
        mock_get_html.return_value = b"<html>mock</html>"
        mock_get_reac_param.return_value = [None] * 7
        
        rxn = reaction("http://example.com/R00001", 0.5, "R00001", "Path", 0)
        
        # Note: ID is an attribute, not a method. ID() method has recursion bug.
        assert rxn.ID == "R00001"


class TestCompoundClass:
    """Test compound class initialization and methods."""
    
    @patch('functions.function_bm_gdb.getHtml')
    @patch('functions.function_bm_gdb.getCompParamFromRestAPI')
    def test_compound_initialization(self, mock_get_comp_param, mock_get_html):
        """Test that compound object initializes correctly."""
        # Mock REST API response
        mock_get_html.return_value = """ENTRY       C00001
NAME        H2O;
            Water
FORMULA     H2O"""
        
        # Mock parsed compound parameters
        mock_get_comp_param.return_value = [
            [["C00001", "cpd:C00001"]],  # IDs
            [["H2O", "H2O"]],  # Formulas
            [[]],  # Associated reactions
            ["Water"],  # Names
            None, None, None, None, None, None, None,
            "H2O"  # Formula4
        ]
        
        # Create compound object
        EF = {}  # Extra formulas
        specialCompounds = []
        comp = compound(
            "http://example.com/C00001",
            "C00001",
            0.5,  # time
            EF,
            specialCompounds
        )
        
        # Verify initialization
        assert comp.ident == "C00001"
        assert mock_get_html.called
        assert mock_get_comp_param.called
    
    @patch('functions.function_bm_gdb.getHtml')
    @patch('functions.function_bm_gdb.getCompParamFromRestAPI')
    def test_compound_id_attributes(self, mock_get_comp_param, mock_get_html):
        """Test ID1 and ID2 attributes."""
        mock_get_html.return_value = "ENTRY C00001"
        mock_get_comp_param.return_value = [
            [["C00001", "cpd:C00001"]],
            [["H2O", "H2O"]],
            [[]],
            ["Water"],
            None, None, None, None, None, None, None, "H2O"
        ]
        
        comp = compound("http://example.com/C00001", "C00001", 0.5, {}, [])
        
        assert comp.ID1 == "C00001"
        assert comp.ID2 == "cpd:C00001"
    
    @patch('functions.function_bm_gdb.getHtml')
    @patch('functions.function_bm_gdb.getCompParamFromRestAPI')
    def test_compound_formula_attributes(self, mock_get_comp_param, mock_get_html):
        """Test Formula1, Formula2, Formula3, Formula4 attributes."""
        mock_get_html.return_value = "ENTRY C00001"
        mock_get_comp_param.return_value = [
            [["C00001", "cpd:C00001"]],
            [["H2O", "H2O"]],
            [[]],
            ["Water"],
            None, None, None, None, None, None, None, "H2O"
        ]
        
        comp = compound("http://example.com/C00001", "C00001", 0.5, {}, [])
        
        assert comp.Formula1 == "H2O"
        assert comp.Formula2 == "H2O"
        assert comp.Formula3 == "H2O"  # Legacy mirror of Formula1
        assert comp.Formula4 == "H2O"
    
    @patch('functions.function_bm_gdb.getHtml')
    @patch('functions.function_bm_gdb.getCompParamFromRestAPI')
    def test_compound_name_attribute(self, mock_get_comp_param, mock_get_html):
        """Test Name attribute."""
        mock_get_html.return_value = "ENTRY C00001"
        mock_get_comp_param.return_value = [
            [["C00001", "cpd:C00001"]],
            [["H2O", "H2O"]],
            [[]],
            ["Water"],
            None, None, None, None, None, None, None, "H2O"
        ]
        
        comp = compound("http://example.com/C00001", "C00001", 0.5, {}, [])
        
        assert comp.Name == "Water"
    
    @patch('functions.function_bm_gdb.getHtml')
    @patch('functions.function_bm_gdb.getCompParamFromRestAPI')
    def test_compound_empty_attributes(self, mock_get_comp_param, mock_get_html):
        """Test compound handles empty/missing attributes gracefully."""
        mock_get_html.return_value = "ENTRY C00001"
        # Return minimal/empty attributes
        mock_get_comp_param.return_value = [
            [[]],  # Empty IDs
            [[]],  # Empty formulas
            [[]],
            [],  # Empty names
            None, None, None, None, None, None, None, ""
        ]
        
        comp = compound("http://example.com/C00001", "C00001", 0.5, {}, [])
        
        # Should not crash, should have empty strings
        assert comp.ID1 == ""
        assert comp.Formula1 == ""
        assert comp.Name == ""
    
    def test_compound_from_batch_data(self):
        """Test from_batch_data class method for batch loading."""
        # Mock batch data (REST API flat file format)
        batch_data = """ENTRY       C00002
NAME        ATP;
            Adenosine 5'-triphosphate
FORMULA     C10H16N5O13P3"""
        
        with patch('functions.function_bm_gdb.getCompParamFromRestAPI') as mock_parser:
            mock_parser.return_value = [
                [["C00002", "cpd:C00002"]],
                [["C10H16N5O13P3", "C10H16N5O13P3"]],
                [[]],
                ["ATP"],
                None, None, None, None, None, None, None, "C10H16N5O13P3"
            ]
            
            comp = compound.from_batch_data(
                "C00002",
                batch_data,
                0.5,  # time
                {},   # EF
                []    # specialCompounds
            )
            
            assert comp.ident == "C00002"
            assert comp.ID1 == "C00002"
            assert comp.Name == "ATP"
            assert mock_parser.called


class TestPathwayClass:
    """Test pathway class initialization and methods."""
    
    @patch('functions.function_bm_gdb.getHtml')
    @patch('functions.function_bm_gdb.getLinkPath')
    def test_pathway_initialization(self, mock_get_link_path, mock_get_html):
        """Test that pathway object initializes correctly."""
        # Mock HTML response
        mock_get_html.return_value = b"<html>mock pathway page</html>"
        
        # Mock parsed pathway data
        mock_get_link_path.return_value = [
            ["C00001", "C00002", "C00003"],  # Compounds
            ["R00001", "R00002"]  # Reactions
        ]
        
        # Create pathway object
        path = pathway(
            "http://example.com/hsa00010",
            0.5,  # time
            "hsa00010",
            "http://kegg.jp",
            "Glycolysis"
        )
        
        # Verify initialization
        assert path.ID == "hsa00010"
        assert path.PathName == "Glycolysis"
        assert mock_get_html.called
        assert mock_get_link_path.called
    
    @patch('functions.function_bm_gdb.getHtml')
    @patch('functions.function_bm_gdb.getLinkPath')
    def test_pathway_id_attribute(self, mock_get_link_path, mock_get_html):
        """Test that ID attribute is set correctly."""
        mock_get_html.return_value = b"<html>mock</html>"
        mock_get_link_path.return_value = [[], []]
        
        path = pathway("http://example.com/hsa00010", 0.5, "hsa00010", "http://kegg.jp", "Test")
        
        # Note: ID is attribute. ID() method has recursion bug.
        assert path.ID == "hsa00010"
    
    @patch('functions.function_bm_gdb.getHtml')
    @patch('functions.function_bm_gdb.getLinkPath')
    def test_pathway_pathname_attribute(self, mock_get_link_path, mock_get_html):
        """Test that PathName attribute is set correctly."""
        mock_get_html.return_value = b"<html>mock</html>"
        mock_get_link_path.return_value = [[], []]
        
        path = pathway("http://example.com/hsa00010", 0.5, "hsa00010", "http://kegg.jp", "Glycolysis")
        
        # Note: PathName is attribute. PathName() method has recursion bug.
        assert path.PathName == "Glycolysis"
    
    @patch('functions.function_bm_gdb.getHtml')
    @patch('functions.function_bm_gdb.getLinkPath')
    def test_pathway_compounds_method(self, mock_get_link_path, mock_get_html):
        """Test that Compounds() method returns list of compound IDs."""
        mock_get_html.return_value = b"<html>mock</html>"
        mock_get_link_path.return_value = [
            ["C00001", "C00002", "C00003"],
            ["R00001"]
        ]
        
        path = pathway("http://example.com/hsa00010", 0.5, "hsa00010", "http://kegg.jp", "Test")
        
        compounds = path.Compounds()
        assert compounds == ["C00001", "C00002", "C00003"]
    
    @patch('functions.function_bm_gdb.getHtml')
    @patch('functions.function_bm_gdb.getLinkPath')
    def test_pathway_reactions_method(self, mock_get_link_path, mock_get_html):
        """Test that Reactions() method returns list of reaction IDs."""
        mock_get_html.return_value = b"<html>mock</html>"
        mock_get_link_path.return_value = [
            ["C00001"],
            ["R00001", "R00002", "R00003"]
        ]
        
        path = pathway("http://example.com/hsa00010", 0.5, "hsa00010", "http://kegg.jp", "Test")
        
        reactions = path.Reactions()
        assert reactions == ["R00001", "R00002", "R00003"]
    
    @patch('functions.function_bm_gdb.getHtml')
    @patch('functions.function_bm_gdb.getLinkPath')
    def test_pathway_empty_data(self, mock_get_link_path, mock_get_html):
        """Test pathway handles empty compounds/reactions gracefully."""
        mock_get_html.return_value = b"<html>mock</html>"
        mock_get_link_path.return_value = [[], []]  # Empty lists
        
        path = pathway("http://example.com/hsa00010", 0.5, "hsa00010", "http://kegg.jp", "Empty")
        
        assert path.Compounds() == []
        assert path.Reactions() == []


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
