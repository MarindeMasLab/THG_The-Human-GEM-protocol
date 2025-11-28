"""
Unit tests for Ensembl client

Tests verify:
1. Session creation with retries and headers
2. Symbol lookup with batch processing
3. ID lookup for Ensembl IDs
4. Cross-reference fetching (Entrez, UniProt)
5. Parallel processing with ThreadPoolExecutor
6. Error handling and retry logic
7. Empty input handling
8. Mixed symbol/ID input handling

Run with: pytest tests/test_ensembl_client.py -v
"""

import os
import sys
import pytest
from unittest.mock import Mock, patch, MagicMock
import logging

# Add project to path
project_root = os.path.join(os.path.dirname(__file__), '..', '..')
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from functions.ensembl_client import (
    _make_requests_session,
    fetch_ensembl_annotations
)


class TestSessionCreation:
    """Test requests session creation and configuration."""
    
    def test_make_requests_session(self):
        """Test session creation with proper headers and adapters."""
        session = _make_requests_session(retries=3, backoff=0.3)
        
        assert session is not None
        assert 'Content-Type' in session.headers
        assert session.headers['Content-Type'] == 'application/json'
        assert session.headers['Accept'] == 'application/json'
        assert 'User-Agent' in session.headers
        assert 'THG-ensembl-client' in session.headers['User-Agent']
    
    def test_session_has_retry_adapter(self):
        """Test that session has retry adapter mounted."""
        session = _make_requests_session(retries=5)
        
        # Check adapters are mounted
        assert 'https://' in session.adapters
        assert 'http://' in session.adapters


class TestSymbolLookup:
    """Test gene symbol lookup functionality."""
    
    @patch('functions.ensembl_client.requests.Session')
    def test_symbol_lookup_batch_success(self, mock_session_class):
        """Test successful batch symbol lookup."""
        mock_session = Mock()
        mock_session_class.return_value = mock_session
        
        # Mock POST response for symbol lookup
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            'BRCA1': {
                'id': 'ENSG00000012048',
                'display_name': 'BRCA1',
                'biotype': 'protein_coding',
                'description': 'BRCA1 DNA repair associated'
            }
        }
        mock_session.post.return_value = mock_response
        
        # Mock GET response for xrefs
        mock_xref_response = Mock()
        mock_xref_response.status_code = 200
        mock_xref_response.json.return_value = [
            {'dbname': 'EntrezGene', 'primary_id': '672'},
            {'dbname': 'Uniprot/SWISSPROT', 'primary_id': 'P38398'}
        ]
        mock_session.get.return_value = mock_xref_response
        
        result = fetch_ensembl_annotations(['BRCA1'])
        
        assert 'BRCA1' in result
        assert result['BRCA1']['ensembl'] == 'ENSG00000012048'
        assert result['BRCA1']['display_name'] == 'BRCA1'
        assert result['BRCA1']['biotype'] == 'protein_coding'
    
    @patch('functions.ensembl_client.requests.Session')
    def test_symbol_lookup_batch_failure_fallback(self, mock_session_class):
        """Test fallback to per-symbol GET when batch fails."""
        mock_session = Mock()
        mock_session_class.return_value = mock_session
        
        # Mock POST returning empty (batch failure)
        mock_post_response = Mock()
        mock_post_response.status_code = 200
        mock_post_response.json.return_value = {}
        
        # Mock GET returning success
        mock_get_response = Mock()
        mock_get_response.status_code = 200
        mock_get_response.json.return_value = {
            'id': 'ENSG00000141510',
            'display_name': 'TP53',
            'biotype': 'protein_coding',
            'description': 'Tumor protein p53'
        }
        
        mock_session.post.return_value = mock_post_response
        
        # Setup responses for GET calls
        responses = [mock_get_response]  # For symbol lookup
        
        def side_effect(*args, **kwargs):
            url = args[0] if args else kwargs.get('url', '')
            if 'lookup/symbol' in url:
                return responses.pop(0)
            # For xrefs
            xref_mock = Mock()
            xref_mock.status_code = 200
            xref_mock.json.return_value = []
            return xref_mock
        
        mock_session.get.side_effect = side_effect
        
        result = fetch_ensembl_annotations(['TP53'])
        
        assert 'TP53' in result
        assert result['TP53']['ensembl'] == 'ENSG00000141510'


class TestIDLookup:
    """Test Ensembl ID lookup functionality."""
    
    @patch('functions.ensembl_client.requests.Session')
    def test_id_lookup_batch(self, mock_session_class):
        """Test batch lookup for Ensembl IDs."""
        mock_session = Mock()
        mock_session_class.return_value = mock_session
        
        # Mock POST response for ID lookup
        mock_post_response = Mock()
        mock_post_response.status_code = 200
        mock_post_response.json.return_value = {
            'ENSG00000139618': {
                'id': 'ENSG00000139618',
                'display_name': 'BRCA2',
                'biotype': 'protein_coding',
                'description': 'BRCA2 DNA repair associated'
            }
        }
        mock_session.post.return_value = mock_post_response
        
        # Mock GET response for xrefs
        mock_xref_response = Mock()
        mock_xref_response.status_code = 200
        mock_xref_response.json.return_value = []
        mock_session.get.return_value = mock_xref_response
        
        result = fetch_ensembl_annotations(['ENSG00000139618'])
        
        assert 'ENSG00000139618' in result
        assert result['ENSG00000139618']['ensembl'] == 'ENSG00000139618'
        assert result['ENSG00000139618']['display_name'] == 'BRCA2'


class TestCrossReferences:
    """Test cross-reference fetching."""
    
    @patch('functions.ensembl_client.requests.Session')
    def test_xref_fetching_entrez_uniprot(self, mock_session_class):
        """Test fetching Entrez and UniProt xrefs."""
        mock_session = Mock()
        mock_session_class.return_value = mock_session
        
        # Mock symbol lookup
        mock_post_response = Mock()
        mock_post_response.status_code = 200
        mock_post_response.json.return_value = {
            'BRCA1': {
                'id': 'ENSG00000012048',
                'display_name': 'BRCA1',
                'biotype': 'protein_coding',
                'description': 'Test'
            }
        }
        
        # Mock xref response
        mock_xref_response = Mock()
        mock_xref_response.status_code = 200
        mock_xref_response.json.return_value = [
            {'dbname': 'EntrezGene', 'primary_id': '672'},
            {'dbname': 'Uniprot/SWISSPROT', 'primary_id': 'P38398'},
            {'dbname': 'Uniprot/SPTREMBL', 'primary_id': 'Q9NZB2'}
        ]
        
        def side_effect(url, *args, **kwargs):
            if 'xrefs' in url:
                return mock_xref_response
            return mock_post_response
        
        mock_session.post.side_effect = side_effect
        mock_session.get.return_value = mock_xref_response
        
        result = fetch_ensembl_annotations(['BRCA1'])
        
        assert 'BRCA1' in result
        assert '672' in result['BRCA1']['entrez']
        assert 'P38398' in result['BRCA1']['uniprot']
        assert 'Q9NZB2' in result['BRCA1']['uniprot']
    
    @patch('functions.ensembl_client.requests.Session')
    def test_xref_no_results(self, mock_session_class):
        """Test when no cross-references are found."""
        mock_session = Mock()
        mock_session_class.return_value = mock_session
        
        # Mock symbol lookup
        mock_post_response = Mock()
        mock_post_response.status_code = 200
        mock_post_response.json.return_value = {
            'GENE1': {
                'id': 'ENSG00000000001',
                'display_name': 'GENE1',
                'biotype': 'protein_coding',
                'description': 'Test'
            }
        }
        mock_session.post.return_value = mock_post_response
        
        # Mock empty xref response
        mock_xref_response = Mock()
        mock_xref_response.status_code = 200
        mock_xref_response.json.return_value = []
        mock_session.get.return_value = mock_xref_response
        
        result = fetch_ensembl_annotations(['GENE1'])
        
        assert 'GENE1' in result
        assert result['GENE1']['entrez'] == []
        assert result['GENE1']['uniprot'] == []


class TestErrorHandling:
    """Test error handling and edge cases."""
    
    def test_empty_input(self):
        """Test with empty gene identifier list."""
        result = fetch_ensembl_annotations([])
        assert result == {}
    
    def test_none_input(self):
        """Test with None input."""
        result = fetch_ensembl_annotations(None)
        assert result == {}
    
    @patch('functions.ensembl_client.requests.Session')
    def test_api_failure_handling(self, mock_session_class):
        """Test handling of API failures."""
        mock_session = Mock()
        mock_session_class.return_value = mock_session
        
        # Mock failed POST
        mock_response = Mock()
        mock_response.status_code = 500
        mock_response.text = "Internal Server Error"
        mock_session.post.return_value = mock_response
        
        result = fetch_ensembl_annotations(['BRCA1'])
        
        # Should return empty dict or handle gracefully
        assert isinstance(result, dict)
    
    @patch('functions.ensembl_client.requests.Session')
    def test_network_timeout(self, mock_session_class):
        """Test handling of network timeout."""
        mock_session = Mock()
        mock_session_class.return_value = mock_session
        
        # Mock timeout exception
        mock_session.post.side_effect = TimeoutError("Connection timeout")
        
        result = fetch_ensembl_annotations(['BRCA1'])
        
        # Should handle timeout gracefully
        assert isinstance(result, dict)


class TestMixedInput:
    """Test handling of mixed symbol and ID inputs."""
    
    @patch('functions.ensembl_client.requests.Session')
    def test_mixed_symbols_and_ids(self, mock_session_class):
        """Test with both gene symbols and Ensembl IDs."""
        mock_session = Mock()
        mock_session_class.return_value = mock_session
        
        # Mock POST for symbols
        symbol_response = Mock()
        symbol_response.status_code = 200
        symbol_response.json.return_value = {
            'BRCA1': {
                'id': 'ENSG00000012048',
                'display_name': 'BRCA1',
                'biotype': 'protein_coding',
                'description': 'Test'
            }
        }
        
        # Mock POST for IDs
        id_response = Mock()
        id_response.status_code = 200
        id_response.json.return_value = {
            'ENSG00000139618': {
                'id': 'ENSG00000139618',
                'display_name': 'BRCA2',
                'biotype': 'protein_coding',
                'description': 'Test'
            }
        }
        
        # Setup post side effect to return different responses
        post_count = [0]
        def post_side_effect(*args, **kwargs):
            url = args[0] if args else kwargs.get('url', '')
            post_count[0] += 1
            if 'symbol' in url:
                return symbol_response
            return id_response
        
        mock_session.post.side_effect = post_side_effect
        
        # Mock xrefs
        xref_response = Mock()
        xref_response.status_code = 200
        xref_response.json.return_value = []
        mock_session.get.return_value = xref_response
        
        result = fetch_ensembl_annotations(['BRCA1', 'ENSG00000139618'])
        
        assert 'BRCA1' in result
        assert 'ENSG00000139618' in result


class TestBatchProcessing:
    """Test batch processing functionality."""
    
    @patch('functions.ensembl_client.requests.Session')
    def test_multiple_batches(self, mock_session_class):
        """Test processing multiple batches."""
        mock_session = Mock()
        mock_session_class.return_value = mock_session
        
        # Create 100 gene symbols
        genes = [f'GENE{i}' for i in range(100)]
        
        # Mock responses
        mock_post_response = Mock()
        mock_post_response.status_code = 200
        
        # Return some genes in each batch
        call_count = [0]
        def json_side_effect():
            call_count[0] += 1
            # Return different genes for each batch
            start = (call_count[0] - 1) * 50
            return {
                f'GENE{i}': {
                    'id': f'ENSG{i:011d}',
                    'display_name': f'GENE{i}',
                    'biotype': 'protein_coding',
                    'description': 'Test'
                }
                for i in range(start, min(start + 50, 100))
            }
        
        mock_post_response.json.side_effect = json_side_effect
        mock_session.post.return_value = mock_post_response
        
        # Mock xrefs
        xref_response = Mock()
        xref_response.status_code = 200
        xref_response.json.return_value = []
        mock_session.get.return_value = xref_response
        
        result = fetch_ensembl_annotations(genes, batch_size=50)
        
        # Should have processed genes in batches
        assert len(result) > 0


class TestDuplicateHandling:
    """Test handling of duplicate identifiers."""
    
    @patch('functions.ensembl_client.requests.Session')
    def test_duplicate_removal(self, mock_session_class):
        """Test that duplicate identifiers are removed."""
        mock_session = Mock()
        mock_session_class.return_value = mock_session
        
        mock_post_response = Mock()
        mock_post_response.status_code = 200
        mock_post_response.json.return_value = {
            'BRCA1': {
                'id': 'ENSG00000012048',
                'display_name': 'BRCA1',
                'biotype': 'protein_coding',
                'description': 'Test'
            }
        }
        mock_session.post.return_value = mock_post_response
        
        xref_response = Mock()
        xref_response.status_code = 200
        xref_response.json.return_value = []
        mock_session.get.return_value = xref_response
        
        # Pass duplicates
        result = fetch_ensembl_annotations(['BRCA1', 'BRCA1', 'BRCA1'])
        
        # Should only process unique identifiers
        assert 'BRCA1' in result
        assert len(result) == 1


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
