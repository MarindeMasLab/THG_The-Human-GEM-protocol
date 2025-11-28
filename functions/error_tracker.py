"""
Comprehensive Error Tracking System for Database Generation

This module provides centralized error tracking and reporting capabilities
for the metabolic model database generation pipeline.

Features:
- Categorized error collection
- Detailed error context (reactions, genes, pathways)
- Summary statistics and reports
- Top error identification
- Success rate tracking per category

Usage:
    tracker = ErrorTracker()
    tracker.add_error('gene_annotation', 'Failed to fetch gene BRCA1', 
                      context={'gene_id': 'BRCA1', 'reaction': 'R00001'})
    summary = tracker.generate_summary()
"""

import logging
from collections import defaultdict, Counter
from datetime import datetime
from typing import Dict, List, Optional, Any
import json


class ErrorTracker:
    """Centralized error tracking system for database generation pipeline."""
    
    # Define error categories
    CATEGORIES = {
        'gene_annotation': 'Gene Annotation Errors',
        'biocyc_query': 'BioCyc Location Query Errors',
        'mass_balance': 'Mass Balance Validation Errors',
        'compartmentalization': 'Compartmentalization Errors',
        'reaction_processing': 'Reaction Processing Errors',
        'compound_processing': 'Compound Processing Errors',
        'kegg_api': 'KEGG API Errors',
        'ensembl_api': 'Ensembl API Errors',
        'data_validation': 'Data Validation Errors',
        'checkpoint': 'Checkpoint Save/Load Errors',
        'other': 'Other Errors'
    }
    
    def __init__(self):
        """Initialize error tracker."""
        self.errors = defaultdict(list)  # category -> list of error dicts
        self.counters = defaultdict(Counter)  # category -> Counter of error messages
        self.success_counters = defaultdict(lambda: {'success': 0, 'failure': 0})
        self.start_time = datetime.now()
        self.logger = logging.getLogger(__name__)
        
    def add_error(self, category: str, message: str, context: Optional[Dict[str, Any]] = None,
                  level: str = 'error'):
        """Add an error to the tracker.
        
        Args:
            category: Error category (must be in CATEGORIES)
            message: Error message
            context: Optional dict with additional context (reaction_id, gene_id, etc.)
            level: Severity level ('error', 'warning', 'critical')
        """
        if category not in self.CATEGORIES:
            self.logger.warning(f"Unknown error category '{category}', using 'other'")
            category = 'other'
        
        error_entry = {
            'timestamp': datetime.now().isoformat(),
            'message': message,
            'context': context or {},
            'level': level
        }
        
        self.errors[category].append(error_entry)
        self.counters[category][message] += 1
        self.success_counters[category]['failure'] += 1
        
        # Log the error
        log_func = getattr(self.logger, level.lower(), self.logger.error)
        context_str = f" | Context: {context}" if context else ""
        log_func(f"[{self.CATEGORIES[category]}] {message}{context_str}")
    
    def add_success(self, category: str):
        """Increment success counter for a category.
        
        Args:
            category: Error category
        """
        if category not in self.CATEGORIES:
            category = 'other'
        self.success_counters[category]['success'] += 1
    
    def get_error_count(self, category: Optional[str] = None) -> int:
        """Get total error count for a category or all categories.
        
        Args:
            category: Specific category or None for total
            
        Returns:
            Error count
        """
        if category:
            return len(self.errors.get(category, []))
        return sum(len(errors) for errors in self.errors.values())
    
    def get_success_rate(self, category: str) -> float:
        """Calculate success rate for a category.
        
        Args:
            category: Error category
            
        Returns:
            Success rate as percentage (0-100)
        """
        stats = self.success_counters[category]
        total = stats['success'] + stats['failure']
        if total == 0:
            return 0.0
        return (stats['success'] / total) * 100
    
    def get_top_errors(self, category: str, n: int = 5) -> List[tuple]:
        """Get top N most common errors in a category.
        
        Args:
            category: Error category
            n: Number of top errors to return
            
        Returns:
            List of (error_message, count) tuples
        """
        return self.counters[category].most_common(n)
    
    def generate_summary(self, verbose: bool = False) -> str:
        """Generate comprehensive error summary report.
        
        Args:
            verbose: If True, include detailed error listings
            
        Returns:
            Formatted summary string
        """
        lines = []
        lines.append("=" * 80)
        lines.append("DATABASE GENERATION ERROR REPORT")
        lines.append("=" * 80)
        lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append(f"Duration: {self._format_duration()}")
        lines.append("")
        
        # Overall statistics
        total_errors = self.get_error_count()
        total_operations = sum(
            stats['success'] + stats['failure'] 
            for stats in self.success_counters.values()
        )
        
        lines.append("OVERALL STATISTICS")
        lines.append("-" * 80)
        lines.append(f"Total operations: {total_operations:,}")
        lines.append(f"Total errors: {total_errors:,}")
        if total_operations > 0:
            overall_success_rate = ((total_operations - total_errors) / total_operations) * 100
            lines.append(f"Overall success rate: {overall_success_rate:.2f}%")
        lines.append("")
        
        # Category-specific statistics
        lines.append("ERRORS BY CATEGORY")
        lines.append("-" * 80)
        
        categories_with_errors = [
            cat for cat in self.CATEGORIES.keys() 
            if self.get_error_count(cat) > 0
        ]
        
        if not categories_with_errors:
            lines.append("✅ No errors recorded!")
        else:
            for category in categories_with_errors:
                count = self.get_error_count(category)
                success_rate = self.get_success_rate(category)
                stats = self.success_counters[category]
                total_ops = stats['success'] + stats['failure']
                
                lines.append(f"\n{self.CATEGORIES[category]}:")
                lines.append(f"  Total errors: {count}")
                lines.append(f"  Operations: {stats['success']} succeeded, {stats['failure']} failed")
                lines.append(f"  Success rate: {success_rate:.2f}%")
                
                # Show top errors
                top_errors = self.get_top_errors(category, n=3)
                if top_errors:
                    lines.append(f"  Top errors:")
                    for msg, error_count in top_errors:
                        # Truncate long messages
                        display_msg = msg[:70] + "..." if len(msg) > 70 else msg
                        lines.append(f"    - {display_msg} ({error_count}x)")
        
        lines.append("")
        lines.append("=" * 80)
        
        # Verbose mode: detailed error listings
        if verbose and categories_with_errors:
            lines.append("\nDETAILED ERROR LISTINGS")
            lines.append("=" * 80)
            
            for category in categories_with_errors:
                lines.append(f"\n{self.CATEGORIES[category]}:")
                lines.append("-" * 80)
                
                errors = self.errors[category][:50]  # Limit to first 50
                for i, error in enumerate(errors, 1):
                    lines.append(f"\n{i}. {error['message']}")
                    if error['context']:
                        lines.append(f"   Context: {json.dumps(error['context'], indent=11)[11:]}")
                    lines.append(f"   Time: {error['timestamp']}")
                    lines.append(f"   Level: {error['level'].upper()}")
                
                if len(self.errors[category]) > 50:
                    remaining = len(self.errors[category]) - 50
                    lines.append(f"\n... and {remaining} more errors")
        
        # Recommendations
        lines.append("\nRECOMMENDATIONS")
        lines.append("-" * 80)
        lines.extend(self._generate_recommendations())
        
        lines.append("")
        lines.append("=" * 80)
        
        return "\n".join(lines)
    
    def _format_duration(self) -> str:
        """Format elapsed time since tracking started."""
        duration = datetime.now() - self.start_time
        hours, remainder = divmod(int(duration.total_seconds()), 3600)
        minutes, seconds = divmod(remainder, 60)
        
        if hours > 0:
            return f"{hours}h {minutes}m {seconds}s"
        elif minutes > 0:
            return f"{minutes}m {seconds}s"
        else:
            return f"{seconds}s"
    
    def _generate_recommendations(self) -> List[str]:
        """Generate recommendations based on error patterns."""
        recommendations = []
        
        # Check gene annotation errors
        if self.get_error_count('gene_annotation') > 0:
            rate = self.get_success_rate('gene_annotation')
            if rate < 50:
                recommendations.append(
                    "⚠️  High gene annotation failure rate. Check Ensembl API connectivity."
                )
            elif rate < 80:
                recommendations.append(
                    "⚡ Moderate gene annotation failures. Some genes may not be in Ensembl database."
                )
        
        # Check BioCyc query errors
        if self.get_error_count('biocyc_query') > 0:
            rate = self.get_success_rate('biocyc_query')
            if rate < 50:
                recommendations.append(
                    "⚠️  High BioCyc query failure rate. Check internet connectivity or BioCyc service status."
                )
        
        # Check mass balance errors
        if self.get_error_count('mass_balance') > 0:
            rate = self.get_success_rate('mass_balance')
            if rate < 70:
                recommendations.append(
                    "⚠️  Many reactions failed mass balance. Check compound formula quality in source data."
                )
        
        # Check compartmentalization errors
        if self.get_error_count('compartmentalization') > 0:
            recommendations.append(
                "⚡ Some reactions have ambiguous compartments. Review gene location assignments."
            )
        
        # Overall health check
        total_errors = self.get_error_count()
        if total_errors == 0:
            recommendations.append("✅ No errors detected! Pipeline completed successfully.")
        elif total_errors < 10:
            recommendations.append("✅ Very few errors. Pipeline is healthy.")
        elif total_errors < 100:
            recommendations.append("⚡ Minor issues detected. Review error details above.")
        else:
            recommendations.append("⚠️  Significant issues detected. Investigate error patterns.")
        
        return recommendations
    
    def save_report(self, filepath: str, verbose: bool = True):
        """Save error report to file.
        
        Args:
            filepath: Output file path
            verbose: Include detailed error listings
        """
        summary = self.generate_summary(verbose=verbose)
        with open(filepath, 'w') as f:
            f.write(summary)
        self.logger.info(f"Error report saved to: {filepath}")
    
    def export_json(self, filepath: str):
        """Export errors to JSON format.
        
        Args:
            filepath: Output JSON file path
        """
        data = {
            'generated_at': datetime.now().isoformat(),
            'duration_seconds': (datetime.now() - self.start_time).total_seconds(),
            'total_errors': self.get_error_count(),
            'categories': {}
        }
        
        for category in self.CATEGORIES.keys():
            if self.get_error_count(category) > 0:
                data['categories'][category] = {
                    'name': self.CATEGORIES[category],
                    'error_count': self.get_error_count(category),
                    'success_rate': self.get_success_rate(category),
                    'errors': self.errors[category],
                    'top_errors': [
                        {'message': msg, 'count': count}
                        for msg, count in self.get_top_errors(category, n=10)
                    ]
                }
        
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2)
        self.logger.info(f"Error data exported to: {filepath}")
    
    def clear(self):
        """Clear all tracked errors and reset counters."""
        self.errors.clear()
        self.counters.clear()
        self.success_counters.clear()
        self.start_time = datetime.now()
