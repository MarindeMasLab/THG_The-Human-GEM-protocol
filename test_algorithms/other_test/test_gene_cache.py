#!/usr/bin/env python3
"""
Test script for gene location caching functionality.
"""

import sys
import os

# Add project root to path
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.append(current_dir)

from functions.class_generate_database import gene

def test_gene_cache():
    """Test gene location caching."""
    
    print("Testing Gene Location Cache")
    print("=" * 80)
    
    # Clear any existing cache
    gene.clear_cache()
    print("\n1. Cache cleared")
    print(f"   Stats: {gene.get_cache_stats()}")
    
    # Test cache miss (first lookup)
    print("\n2. Testing cache MISS (first lookup)")
    result = gene.get_locations("BRCA1", "HS12345")
    print(f"   Result: {result}")
    print(f"   Stats: {gene.get_cache_stats()}")
    assert result is None, "First lookup should return None (cache miss)"
    
    # Manually cache a location
    print("\n3. Manually caching location for BRCA1")
    gene.cache_location("BRCA1", ["nucleus", "cytosol"], "HS12345")
    print(f"   Stats: {gene.get_cache_stats()}")
    
    # Test cache hit (second lookup)
    print("\n4. Testing cache HIT (second lookup)")
    result = gene.get_locations("BRCA1", "HS12345")
    print(f"   Result: {result}")
    print(f"   Stats: {gene.get_cache_stats()}")
    assert result == ["nucleus", "cytosol"], "Second lookup should return cached data"
    
    # Test another gene (cache miss)
    print("\n5. Testing another gene (cache MISS)")
    result = gene.get_locations("TP53", "HS67890")
    print(f"   Result: {result}")
    assert result is None, "New gene should be cache miss"
    
    # Cache second gene
    gene.cache_location("TP53", ["nucleus"], "HS67890")
    
    # Test third lookup on first gene (cache hit)
    print("\n6. Testing BRCA1 again (cache HIT)")
    result = gene.get_locations("BRCA1")
    print(f"   Result: {result}")
    assert result == ["nucleus", "cytosol"], "Should still be cached"
    
    # Final stats
    print("\n7. Final cache statistics:")
    stats = gene.get_cache_stats()
    for key, value in stats.items():
        print(f"   {key}: {value}")
    
    # Test save/load
    print("\n8. Testing cache persistence")
    cache_file = "/tmp/test_gene_cache.pkl"
    gene.save_cache(cache_file)
    print(f"   Saved cache to {cache_file}")
    
    # Clear and verify empty
    gene.clear_cache()
    result = gene.get_locations("BRCA1")
    assert result is None, "Cache should be empty after clear"
    print("   Cache cleared")
    
    # Load and verify
    success = gene.load_cache(cache_file)
    assert success, "Cache load should succeed"
    result = gene.get_locations("BRCA1")
    assert result == ["nucleus", "cytosol"], "Cache should be restored"
    print(f"   Cache loaded successfully")
    print(f"   BRCA1 location: {result}")
    
    # Clean up
    os.remove(cache_file)
    
    print("\n" + "=" * 80)
    print("✅ All tests passed!")
    print("=" * 80)

if __name__ == "__main__":
    test_gene_cache()
