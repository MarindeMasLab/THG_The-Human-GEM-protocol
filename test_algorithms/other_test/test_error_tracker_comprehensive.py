"""
Unit tests for ErrorTracker class

Tests verify:
1. Error addition to correct categories
2. Success/failure counting
3. Summary report generation
4. JSON export
5. Top error identification
6. Success rate calculation
7. Edge cases (empty tracker, unknown categories)

Run with: pytest tests/test_error_tracker.py -v
"""

import os
import sys
import pytest
import json
from datetime import datetime
from io import StringIO

# Add project to path
project_root = os.path.join(os.path.dirname(__file__), '..', '..')
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from functions.error_tracker import ErrorTracker


class TestErrorTrackerBasics:
    """Test basic ErrorTracker functionality."""
    
    def test_initialization(self):
        """Test that ErrorTracker initializes correctly."""
        tracker = ErrorTracker()
        
        assert len(tracker.errors) == 0
        assert len(tracker.counters) == 0
        assert tracker.start_time is not None
        assert len(ErrorTracker.CATEGORIES) == 11
    
    def test_add_error_basic(self):
        """Test adding a basic error."""
        tracker = ErrorTracker()
        
        tracker.add_error('gene_annotation', 'Failed to fetch gene BRCA1')
        
        assert tracker.get_error_count('gene_annotation') == 1
        assert len(tracker.errors['gene_annotation']) == 1
        assert tracker.errors['gene_annotation'][0]['message'] == 'Failed to fetch gene BRCA1'
    
    def test_add_error_with_context(self):
        """Test adding error with context."""
        tracker = ErrorTracker()
        
        context = {'gene_id': 'BRCA1', 'reaction_id': 'R00001'}
        tracker.add_error('gene_annotation', 'API timeout', context=context)
        
        error = tracker.errors['gene_annotation'][0]
        assert error['context'] == context
        assert error['message'] == 'API timeout'
    
    def test_add_error_with_level(self):
        """Test adding errors with different severity levels."""
        tracker = ErrorTracker()
        
        tracker.add_error('mass_balance', 'Cannot balance', level='error')
        tracker.add_error('mass_balance', 'Minor issue', level='warning')
        tracker.add_error('mass_balance', 'Critical failure', level='critical')
        
        assert tracker.get_error_count('mass_balance') == 3
        assert tracker.errors['mass_balance'][0]['level'] == 'error'
        assert tracker.errors['mass_balance'][1]['level'] == 'warning'
        assert tracker.errors['mass_balance'][2]['level'] == 'critical'
    
    def test_add_error_unknown_category(self):
        """Test that unknown category defaults to 'other'."""
        tracker = ErrorTracker()
        
        tracker.add_error('unknown_category', 'Some error')
        
        assert tracker.get_error_count('other') == 1
        assert tracker.get_error_count('unknown_category') == 0
    
    def test_add_success(self):
        """Test adding success counts."""
        tracker = ErrorTracker()
        
        tracker.add_success('gene_annotation')
        tracker.add_success('gene_annotation')
        tracker.add_success('mass_balance')
        
        assert tracker.success_counters['gene_annotation']['success'] == 2
        assert tracker.success_counters['mass_balance']['success'] == 1


class TestErrorTrackerCounting:
    """Test error counting and statistics."""
    
    def test_get_error_count_by_category(self):
        """Test getting error count for specific category."""
        tracker = ErrorTracker()
        
        tracker.add_error('gene_annotation', 'Error 1')
        tracker.add_error('gene_annotation', 'Error 2')
        tracker.add_error('mass_balance', 'Error 3')
        
        assert tracker.get_error_count('gene_annotation') == 2
        assert tracker.get_error_count('mass_balance') == 1
        assert tracker.get_error_count('biocyc_query') == 0
    
    def test_get_total_error_count(self):
        """Test getting total error count across all categories."""
        tracker = ErrorTracker()
        
        tracker.add_error('gene_annotation', 'Error 1')
        tracker.add_error('mass_balance', 'Error 2')
        tracker.add_error('kegg_api', 'Error 3')
        
        assert tracker.get_error_count() == 3
    
    def test_get_success_rate(self):
        """Test calculating success rate."""
        tracker = ErrorTracker()
        
        # 8 successes, 2 failures = 80% success rate
        for _ in range(8):
            tracker.add_success('gene_annotation')
        for _ in range(2):
            tracker.add_error('gene_annotation', 'Failed')
        
        rate = tracker.get_success_rate('gene_annotation')
        assert rate == 80.0
    
    def test_get_success_rate_no_operations(self):
        """Test success rate when no operations recorded."""
        tracker = ErrorTracker()
        
        rate = tracker.get_success_rate('gene_annotation')
        assert rate == 0.0  # Returns 0.0 if no operations
    
    def test_get_success_rate_all_failures(self):
        """Test success rate with only failures."""
        tracker = ErrorTracker()
        
        tracker.add_error('mass_balance', 'Failed')
        tracker.add_error('mass_balance', 'Failed')
        
        rate = tracker.get_success_rate('mass_balance')
        assert rate == 0.0


class TestErrorTrackerReporting:
    """Test report generation."""
    
    def test_get_top_errors(self):
        """Test getting most common errors."""
        tracker = ErrorTracker()
        
        tracker.add_error('gene_annotation', 'Timeout')
        tracker.add_error('gene_annotation', 'Timeout')
        tracker.add_error('gene_annotation', 'Timeout')
        tracker.add_error('gene_annotation', 'Not found')
        tracker.add_error('gene_annotation', 'Not found')
        tracker.add_error('gene_annotation', 'Invalid ID')
        
        top_errors = tracker.get_top_errors('gene_annotation', n=2)
        
        assert len(top_errors) == 2
        assert top_errors[0][0] == 'Timeout'
        assert top_errors[0][1] == 3
        assert top_errors[1][0] == 'Not found'
        assert top_errors[1][1] == 2
    
    def test_generate_summary(self):
        """Test summary generation."""
        tracker = ErrorTracker()
        
        tracker.add_error('gene_annotation', 'Error 1')
        tracker.add_error('mass_balance', 'Error 2')
        tracker.add_success('gene_annotation')
        
        summary = tracker.generate_summary()
        
        assert 'Total errors:' in summary
        assert 'Gene Annotation Errors' in summary
        assert 'Duration:' in summary
        assert '=' * 80 in summary
    
    def test_generate_text_report(self):
        """Test text report generation."""
        tracker = ErrorTracker()
        
        tracker.add_error('gene_annotation', 'Failed gene fetch', 
                         context={'gene_id': 'BRCA1'})
        tracker.add_success('gene_annotation')
        
        # generate_summary returns the text report
        report = tracker.generate_summary()
        
        assert 'DATABASE GENERATION ERROR REPORT' in report
        assert 'Gene Annotation Errors' in report
        assert 'Failed gene fetch' in report
        assert 'Duration' in report
    
    def test_generate_json_export(self):
        """Test JSON export generation."""
        tracker = ErrorTracker()
        
        tracker.add_error('mass_balance', 'Cannot balance equation',
                         context={'reaction_id': 'R00001'})
        
        # Note: export_json saves to file, we'll test save_report instead
        # Just test that we can access the data structure
        assert 'mass_balance' in tracker.errors
        assert len(tracker.errors['mass_balance']) == 1


class TestErrorTrackerEdgeCases:
    """Test edge cases and error handling."""
    
    def test_empty_tracker_summary(self):
        """Test summary generation on empty tracker."""
        tracker = ErrorTracker()
        
        summary = tracker.generate_summary()
        
        assert 'Total errors: 0' in summary or 'No errors recorded' in summary
    
    def test_empty_tracker_text_report(self):
        """Test text report on empty tracker."""
        tracker = ErrorTracker()
        
        # generate_summary returns the text report
        report = tracker.generate_summary()
        
        assert 'No errors recorded' in report or 'Total errors: 0' in report
    
    def test_multiple_categories(self):
        """Test tracking errors across multiple categories."""
        tracker = ErrorTracker()
        
        categories = ['gene_annotation', 'mass_balance', 'kegg_api', 'biocyc_query']
        for cat in categories:
            tracker.add_error(cat, f'Error in {cat}')
        
        assert tracker.get_error_count() == 4
        for cat in categories:
            assert tracker.get_error_count(cat) == 1
    
    def test_error_timestamp(self):
        """Test that errors have timestamps."""
        tracker = ErrorTracker()
        
        tracker.add_error('gene_annotation', 'Test error')
        
        error = tracker.errors['gene_annotation'][0]
        assert 'timestamp' in error
        # Should be ISO format datetime
        datetime.fromisoformat(error['timestamp'])  # Will raise if invalid
    
    def test_counter_increments(self):
        """Test that error counters increment correctly."""
        tracker = ErrorTracker()
        
        tracker.add_error('mass_balance', 'Same error')
        tracker.add_error('mass_balance', 'Same error')
        tracker.add_error('mass_balance', 'Different error')
        
        assert tracker.counters['mass_balance']['Same error'] == 2
        assert tracker.counters['mass_balance']['Different error'] == 1


class TestErrorTrackerRecommendations:
    """Test error recommendation system."""
    
    def test_recommendations_in_summary(self):
        """Test that recommendations appear in summary."""
        tracker = ErrorTracker()
        
        tracker.add_error('gene_annotation', 'Timeout')
        
        summary = tracker.generate_summary()
        
        # Should have RECOMMENDATIONS section
        assert 'RECOMMENDATIONS' in summary


class TestErrorTrackerSaving:
    """Test saving reports to files."""
    
    def test_save_report(self, tmp_path):
        """Test saving text report to file."""
        tracker = ErrorTracker()
        tracker.add_error('gene_annotation', 'Test error')
        
        report_file = tmp_path / "error_report.txt"
        tracker.save_report(str(report_file))
        
        assert report_file.exists()
        content = report_file.read_text()
        assert 'Test error' in content
    
    def test_export_json(self, tmp_path):
        """Test exporting JSON data to file."""
        tracker = ErrorTracker()
        tracker.add_error('mass_balance', 'Test error')
        
        json_file = tmp_path / "error_report.json"
        tracker.export_json(str(json_file))
        
        assert json_file.exists()
        data = json.loads(json_file.read_text())
        assert 'categories' in data
    
    def test_clear(self):
        """Test clearing tracker data."""
        tracker = ErrorTracker()
        tracker.add_error('gene_annotation', 'Error 1')
        tracker.add_success('gene_annotation')
        
        assert tracker.get_error_count() > 0
        
        tracker.clear()
        
        assert tracker.get_error_count() == 0
        assert len(tracker.errors) == 0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
