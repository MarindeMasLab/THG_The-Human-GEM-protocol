"""
Shared Metabolite Matching Utilities
=====================================

These functions provide a unified approach for matching metabolites across
different data sources (KEGG, ChEBI, Rhea, etc.) using multiple identifiers.

Used by: merge_metabolic_networks, rhea_extension, and other modules.

This is a lightweight module with minimal dependencies to avoid import issues.
"""

import re
import logging

LOGGER = logging.getLogger(__name__)

# Standard attributes used for metabolite matching (priority order)
METABOLITE_MATCH_ATTRIBUTES = ['kegg.compound', 'chebi', 'inchikey', 'inchi', 'pubchem.compound', 'lipidmaps']


def build_metabolite_attribute_index(compounds: dict, attribute_getters: dict = None) -> dict:
    """
    Build an index of metabolite attributes for fast lookup.
    
    This function scans a collection of compounds and builds reverse lookup
    dictionaries for each attribute type, enabling O(1) matching.
    
    Args:
        compounds: Dict of compound_id -> compound_object
                  Works with both raw compound objects (from generate_db) 
                  and cobra.Metabolite objects (from COBRA models)
        attribute_getters: Optional dict mapping attribute_name -> function(compound) -> value
                          If None, uses default getters for common attributes
    
    Returns:
        Dict of {attribute_name: {attribute_value: compound_id}}
        
    Example:
        >>> index = build_metabolite_attribute_index(MetList)
        >>> index['chebi']['15377']  # Returns 'C00001' (water)
    """
    if attribute_getters is None:
        # Default getters that work with both raw compound objects and COBRA metabolites
        attribute_getters = {
            'chebi': lambda c: _get_attr_normalized(c, ['CheBI', 'chebi_id'], prefix='CHEBI:'),
            'inchi': lambda c: _get_attr(c, ['inchi']),
            'inchikey': lambda c: _get_attr(c, ['inchikey']),
            'pubchem.compound': lambda c: _get_attr(c, ['PubChem', 'CID']),
            'lipidmaps': lambda c: _get_attr(c, ['LIPIDMAPS']),
            'kegg.compound': lambda c: _get_kegg_id(c),
            'name': lambda c: _get_attr_lower(c, ['Name', 'name']),
        }
    
    index = {attr: {} for attr in attribute_getters.keys()}
    
    for compound_id, compound in compounds.items():
        for attr_name, getter in attribute_getters.items():
            try:
                value = getter(compound)
                if value and str(value).strip():
                    normalized = str(value).strip()
                    if normalized and normalized.lower() not in ('none', 'nan', ''):
                        index[attr_name][normalized] = compound_id
            except (AttributeError, KeyError, TypeError):
                continue
                
    return index


def find_equivalent_metabolite(chebi_id: str = None, inchi: str = None, 
                                inchikey: str = None, pubchem: str = None,
                                lipidmaps: str = None, name: str = None,
                                attribute_index: dict = None) -> str:
    """
    Find an existing metabolite using multiple identifiers (priority-based).
    
    Checks identifiers in order of reliability:
    1. ChEBI (curated, high confidence)
    2. InChIKey (structure-based, unique)
    3. InChI (structure-based)
    4. PubChem
    5. LIPID MAPS
    6. Name (exact match only, lowest confidence)
    
    Args:
        chebi_id: ChEBI identifier (with or without 'CHEBI:' prefix)
        inchi: InChI string
        inchikey: InChIKey
        pubchem: PubChem CID
        lipidmaps: LIPID MAPS ID
        name: Compound name (exact match)
        attribute_index: Index built by build_metabolite_attribute_index()
        
    Returns:
        Compound ID if found, None otherwise
        
    Example:
        >>> kegg_id = find_equivalent_metabolite(chebi_id='15377', attribute_index=index)
        >>> print(kegg_id)  # 'C00001' (water)
    """
    if attribute_index is None:
        return None
        
    # 1. Check ChEBI
    if chebi_id:
        chebi_norm = str(chebi_id).replace('CHEBI:', '').strip()
        if chebi_norm in attribute_index.get('chebi', {}):
            return attribute_index['chebi'][chebi_norm]
            
    # 2. Check InChIKey
    if inchikey and inchikey in attribute_index.get('inchikey', {}):
        return attribute_index['inchikey'][inchikey]
        
    # 3. Check InChI
    if inchi and inchi in attribute_index.get('inchi', {}):
        return attribute_index['inchi'][inchi]
        
    # 4. Check PubChem
    if pubchem:
        pubchem_str = str(pubchem).strip()
        if pubchem_str in attribute_index.get('pubchem.compound', {}):
            return attribute_index['pubchem.compound'][pubchem_str]
            
    # 5. Check LIPID MAPS
    if lipidmaps and lipidmaps in attribute_index.get('lipidmaps', {}):
        return attribute_index['lipidmaps'][lipidmaps]
        
    # 6. Check name (exact match, lowest priority)
    if name:
        name_norm = name.lower().strip()
        if name_norm in attribute_index.get('name', {}):
            return attribute_index['name'][name_norm]
            
    return None


# =============================================================================
# Helper functions for attribute extraction
# =============================================================================

def _get_attr(compound, attr_names: list) -> str:
    """Get first available attribute from compound."""
    for attr in attr_names:
        if hasattr(compound, attr):
            val = getattr(compound, attr)
            if val:
                return str(val)
        # Also check annotation dict (for COBRA metabolites)
        if hasattr(compound, 'annotation') and attr in compound.annotation:
            val = compound.annotation[attr]
            if isinstance(val, list):
                val = val[0] if val else None
            if val:
                return str(val)
    return None


def _get_attr_normalized(compound, attr_names: list, prefix: str = '') -> str:
    """Get attribute and remove prefix if present."""
    val = _get_attr(compound, attr_names)
    if val and prefix:
        return val.replace(prefix, '').strip()
    return val


def _get_attr_lower(compound, attr_names: list) -> str:
    """Get attribute and convert to lowercase."""
    val = _get_attr(compound, attr_names)
    if val:
        return val.lower().strip()
    return None


def _get_kegg_id(compound) -> str:
    """Extract KEGG compound ID from various sources."""
    # Direct ID attributes
    for attr in ['ID1', 'ID2', 'ident']:
        if hasattr(compound, attr):
            val = getattr(compound, attr)
            if val and isinstance(val, str) and re.match(r'^[CG]\d{5}$', val):
                return val
    # From annotation (COBRA)
    if hasattr(compound, 'annotation') and 'kegg.compound' in compound.annotation:
        val = compound.annotation['kegg.compound']
        if isinstance(val, list):
            val = val[0] if val else None
        if val:
            return val
    return None
