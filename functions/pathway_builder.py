#!/usr/bin/env python3
"""
Generic Pathway Builder Module

This module provides generic functions for adding any pathway/compartment
to metabolic models. It's designed to be pathway-agnostic and work with
any compartment configuration from config.json.

Usage:
    from functions.pathway_builder import add_compartment, create_compartment_metabolites
    
    # Add any compartment
    model = add_compartment(model, 'gc', 'glycocalyx')
    
    # Create metabolites for any pathway
    create_compartment_metabolites(model, id_database, 'gc', config, 'glycocalyx_specific')
"""

import json
import copy
import os
import re
import sys


# ============================================================================
# Metabolite Matching Functions (Online/Offline Hybrid)
# ============================================================================

def find_metabolite_by_formula_in_model(model, formula, compartment=None):
    """
    Find a metabolite in the model by its chemical formula.
    
    Args:
        model: Model dictionary
        formula: Chemical formula string (e.g., "C6H9N3O2")
        compartment: Optional compartment to search in. If None, searches all compartments.
        
    Returns:
        Metabolite dict if found, None otherwise
    """
    if not formula:
        return None
        
    for met in model.get('metabolites', []):
        met_formula = met.get('formula', '')
        met_compartment = met.get('compartment', '')
        
        # Match formula exactly
        if met_formula == formula:
            if compartment is None or met_compartment == compartment:
                return met
                
    return None


def get_metabolite_id_base(met_id):
    """
    Extract the base ID from a metabolite ID (remove compartment suffix).
    
    Args:
        met_id: Metabolite ID (e.g., "MAM02125c")
        
    Returns:
        Base ID without compartment (e.g., "MAM02125")
    """
    # Match pattern like MAM#####X where X is compartment letter(s)
    match = re.match(r'^([A-Z]+\d+)', met_id)
    if match:
        return match.group(1)
    return met_id


def create_metabolite_in_new_compartment(model, source_met, target_compartment):
    """
    Create a copy of a metabolite in a new compartment with the same formula.
    
    Args:
        model: Model dictionary
        source_met: Source metabolite dictionary to copy from
        target_compartment: Target compartment abbreviation
        
    Returns:
        New metabolite ID in target compartment
    """
    # Get base ID and create new ID with target compartment
    base_id = get_metabolite_id_base(source_met['id'])
    new_id = f"{base_id}{target_compartment}"
    
    # Check if it already exists
    existing = next((m for m in model['metabolites'] if m['id'] == new_id), None)
    if existing:
        return new_id
    
    # Create new metabolite with same formula and annotations
    new_met = {
        'id': new_id,
        'name': source_met.get('name', ''),
        'compartment': target_compartment,
        'formula': source_met.get('formula', ''),  # Critical: preserve formula!
        'charge': source_met.get('charge', 0),
        'annotation': copy.deepcopy(source_met.get('annotation', {}))
    }
    
    model['metabolites'].append(new_met)
    print(f"  ✓ Created {new_id} ({new_met['name']}) in [{target_compartment}] with formula {new_met['formula']}")
    
    return new_id


# Utility Functions
def get_next_metabolite_id(model):
    """Get the next available metabolite ID in the MAM##### format."""
    max_id = 0
    for met in model['metabolites']:
        if met['id'].startswith('MAM') and len(met['id']) >= 8:
            try:
                num = int(met['id'][3:8])
                if num > max_id:
                    max_id = num
            except ValueError:
                continue
    return f"MAM{max_id + 1:05d}"


def get_next_reaction_id(model):
    """Get the next available reaction ID in the MAR##### format."""
    max_id = 0
    for rxn in model['reactions']:
        if rxn['id'].startswith('MAR') and len(rxn['id']) >= 8:
            try:
                num = int(rxn['id'][3:8])
                if num > max_id:
                    max_id = num
            except ValueError:
                continue
    return f"MAR{max_id + 1:05d}"


def find_metabolite_by_annotation(model, metabolite_name, id_database, compartment='c'):
    """
    Find a metabolite by matching its annotations against the ID database.
    Returns: metabolite dict if found, None otherwise
    """
    if metabolite_name not in id_database:
        return None
    
    known_ids = id_database[metabolite_name]
    
    # Priority weights for different annotation types
    priority_weights = {
        'inchikey': 10,
        'inchi': 9,
        'chebi': 8,
        'kegg.compound': 7,
        'bigg.metabolite': 6,
        'pubchem.compound': 5,
        'vmhmetabolite': 4,
        'metanetx.chemical': 3,
    }
    
    best_match = None
    best_score = 0
    best_match_count = 0
    
    for met in model['metabolites']:
        if met['compartment'] != compartment:
            continue
        
        if 'annotation' not in met or not met['annotation']:
            continue
        
        match_count = 0
        score = 0
        
        for key, value in met['annotation'].items():
            if key in known_ids and value == known_ids[key]:
                match_count += 1
                score += priority_weights.get(key, 1)
        
        if match_count >= 2 and score > best_score:
            best_match = met
            best_score = score
            best_match_count = match_count
    
    return best_match


def find_metabolite_robust(model, metabolite_name, id_database, compartment='c', auto_create=False):
    """
    Find a metabolite using multiple strategies:
    1. Exact name match (case-sensitive)
    2. Case-insensitive exact match
    3. Annotation-based matching
    4. Case-insensitive substring match
    5. Formula-based matching (finds in other compartments and creates copy)
    
    Args:
        model: Model dictionary
        metabolite_name: Name of metabolite to find
        id_database: Metabolite ID database
        compartment: Target compartment abbreviation
        auto_create: If True, automatically creates metabolite in target compartment
                     if found in another compartment (preserves formula)
    
    Returns: metabolite dict if found, None otherwise
    """
    # Strategy 1: Exact name match (case-sensitive)
    for met in model['metabolites']:
        if met['name'] == metabolite_name and met['compartment'] == compartment:
            return met
    
    # Strategy 2: Case-insensitive exact match
    for met in model['metabolites']:
        if met['name'].lower() == metabolite_name.lower() and met['compartment'] == compartment:
            return met
    
    # Strategy 3: Annotation-based matching
    met = find_metabolite_by_annotation(model, metabolite_name, id_database, compartment)
    if met:
        return met
    
    # Strategy 4: Case-insensitive substring match
    for met in model['metabolites']:
        if metabolite_name.lower() in met['name'].lower() and met['compartment'] == compartment:
            return met
    
    # Strategy 5: Formula-based matching (new!) - find in other compartments
    if auto_create:
        # Try to find the formula for this metabolite
        formula = None
        source_met = None
        
        # 5a: Try to find metabolite with matching name in ANY compartment
        for met in model['metabolites']:
            if met['name'].lower() == metabolite_name.lower():
                source_met = met
                if met.get('formula'):
                    formula = met['formula']
                    break
        
        # 5b: If not found by name, check the database for this metabolite
        if not formula and metabolite_name in id_database:
            db_entry = id_database[metabolite_name]
            if 'id' in db_entry:
                # Get the base ID (remove compartment suffix if present)
                base_id = get_metabolite_id_base(db_entry['id'])
                # Search for any metabolite with this base ID in ANY compartment
                for met in model['metabolites']:
                    met_base = get_metabolite_id_base(met['id'])
                    if met_base == base_id:
                        source_met = met
                        if met.get('formula'):
                            formula = met['formula']
                            break
        
        if formula:
            # Check if we already have this formula in target compartment
            existing = find_metabolite_by_formula_in_model(model, formula, compartment)
            if existing:
                return existing
            
            # Create new metabolite in target compartment with same formula
            if source_met:
                new_id = create_metabolite_in_new_compartment(model, source_met, compartment)
                # Return the newly created metabolite
                return next((m for m in model['metabolites'] if m['id'] == new_id), None)
    
    return None


def add_compartment(model, compartment_abbrev, compartment_name):
    """
    Add a compartment to the model (generic version).
    
    Args:
        model: The metabolic model dictionary
        compartment_abbrev: Abbreviation for the compartment (e.g., 'ck', 'gc', 'ecm')
        compartment_name: Full name of the compartment (e.g., 'cytoskeleton', 'glycocalyx')
    
    Returns:
        tuple: (model, compartment_abbrev, already_existed)
    """
    already_existed = False
    
    if compartment_abbrev in model['compartments']:
        print(f"\nWARNING: Compartment [{compartment_abbrev}] '{compartment_name}' already exists.")
        already_existed = True
    else:
        print(f"\nAdding compartment [{compartment_abbrev}] '{compartment_name}'...")
        model['compartments'][compartment_abbrev] = compartment_name
        print(f"Compartment [{compartment_abbrev}] added.")
    
    return model, compartment_abbrev, already_existed


def create_compartment_metabolites(model, id_database, compartment_abbrev, config, pathway_key=None, abbrev_map=None):
    """
    Create metabolites for a pathway in a specified compartment (generic version).
    
    Args:
        model: The metabolic model dictionary
        id_database: Database of metabolite IDs for annotation-based matching
        compartment_abbrev: Target compartment abbreviation (e.g., 'ck', 'gc')
        config: Configuration dictionary
        pathway_key: Key in config['metabolites'] for pathway-specific metabolites
                    (e.g., 'cytoskeleton_specific', 'glycocalyx_specific')
                    If None, will try to infer from compartment name
        abbrev_map: Optional dict mapping config abbreviations to model abbreviations
                   (e.g., {'er': 'r', 'gl': 'gc'})
    
    Returns:
        int: Number of metabolites created
    """
    print(f"\nCreating metabolites in [{compartment_abbrev}]...")
    
    new_metabolites = []
    
    # Helper function to create metabolite in target compartment
    def create_compartment_metabolite(metabolite_name):
        """Create a compartment version of a metabolite using annotation-based matching."""
        # Use robust matching to find the metabolite in cytosol
        base_met = find_metabolite_robust(model, metabolite_name, id_database, 'c')
        
        if not base_met:
            print(f"  ⚠ Warning: Could not find metabolite '{metabolite_name}' in cytosol")
            return None
        
        # Check if target compartment version already exists
        target_id = base_met['id'][:-1] + compartment_abbrev
        for met in model['metabolites']:
            if met['id'] == target_id:
                print(f"  ✓ Metabolite {target_id} ({base_met['name']}) already exists in [{compartment_abbrev}]")
                return target_id
        
        # Create new metabolite in target compartment
        new_met = copy.deepcopy(base_met)
        new_met['id'] = target_id
        new_met['compartment'] = compartment_abbrev
        new_metabolites.append(new_met)
        print(f"  ✓ Created {target_id} ({base_met['name']}) in [{compartment_abbrev}]")
        return target_id
    
    # Create shared metabolites from targets list
    # ONLY if they're actually needed in this compartment (check reactions)
    if 'targets' in config.get('metabolites', {}):
        # Parse reactions to find which target metabolites are needed in this compartment
        needed_targets = set()
        for reaction in config.get('reactions', []):
            equation = reaction.get('equation', '')
            # Extract metabolites in this compartment from equation
            # Pattern: metabolite_name[compartment]
            import re
            pattern = r'([A-Za-z0-9\-\+\(\),\. ]+)\[' + re.escape(compartment_abbrev) + r'\]'
            matches = re.findall(pattern, equation)
            for match in matches:
                met_name = match.strip()
                # Check if this is a target metabolite
                if met_name in config['metabolites']['targets']:
                    needed_targets.add(met_name)
        
        # Only create target metabolites that are actually used in reactions in this compartment
        for met_name in needed_targets:
            create_compartment_metabolite(met_name)
    
    # Create pathway-specific metabolites
    if pathway_key:
        # Handle nested keys (e.g., ('metabolites', 'Syndecan1_HS_specific'))
        if isinstance(pathway_key, tuple):
            parent_key, child_key = pathway_key
            if parent_key in config and child_key in config[parent_key]:
                pathway_metabolites = config[parent_key][child_key]
            else:
                pathway_metabolites = None
        elif pathway_key in config.get('metabolites', {}):
            pathway_metabolites = config['metabolites'][pathway_key]
        elif pathway_key in config:
            pathway_metabolites = config[pathway_key]
        else:
            pathway_metabolites = None
        
        if pathway_metabolites and isinstance(pathway_metabolites, dict):
            print(f"  ✓ Loaded {len(pathway_metabolites)} pathway-specific metabolite definitions from config")
        
            # Filter metabolites for this compartment
            compartment_metabolites = {}
            for met_name, met_data in pathway_metabolites.items():
                # Get compartment from metabolite data
                config_compartment = met_data.get('compartment', compartment_abbrev)
                
                # Map config compartment to model compartment if mapping provided
                model_compartment = abbrev_map.get(config_compartment, config_compartment) if abbrev_map else config_compartment
                
                # Only include if it matches current compartment
                if model_compartment == compartment_abbrev:
                    compartment_metabolites[met_name] = met_data
            
            if not compartment_metabolites:
                print(f"  ℹ No pathway-specific metabolites for compartment [{compartment_abbrev}]")
            else:
                print(f"  ✓ Found {len(compartment_metabolites)} metabolites for compartment [{compartment_abbrev}]")
            
            base_id = int(get_next_metabolite_id(model)[3:8])
            for idx, (met_name, met_data) in enumerate(compartment_metabolites.items()):
                # Use the compartment specified in the metabolite data, mapped to model abbreviation
                config_comp = met_data.get('compartment', compartment_abbrev)
                met_compartment = abbrev_map.get(config_comp, config_comp) if abbrev_map else config_comp
                
                # Generate new ID
                new_id = f"MAM{base_id + idx:05d}{met_compartment}"
                
                # Check if already exists
                exists = False
                for met in model['metabolites']:
                    if met['id'] == new_id or (met['name'] == met_name and met['compartment'] == met_compartment):
                        exists = True
                        break
                
                if exists:
                    print(f"  ✓ Metabolite {new_id} already exists")
                    continue
                
                # Create new metabolite
                new_met = {
                    'id': new_id,
                    'name': met_name,
                    'compartment': met_compartment,
                    'charge': met_data.get('charge', 0),
                    'formula': met_data.get('formula', ''),
                }
                
                # Add annotations if present
                annotation = {}
                
                # Check for annotation sub-dict first (new format)
                if 'annotation' in met_data and isinstance(met_data['annotation'], dict):
                    annotation.update(met_data['annotation'])
                
                # Also check for direct annotation fields (old format)
                for key in ['bigg.metabolite', 'chebi', 'hmdb', 'inchi', 'inchikey', 
                           'kegg.compound', 'metanetx.chemical', 'pubchem.compound', 
                           'uniprot', 'vmhmetabolite', 'hgnc', 'ensembl']:
                    if key in met_data and met_data[key]:
                        annotation[key] = met_data[key]
                
                if annotation:
                    new_met['annotation'] = annotation
                
                # Add SBO term if present
                if 'sbo' in met_data:
                    new_met['annotation'] = new_met.get('annotation', {})
                    new_met['annotation']['sbo'] = met_data['sbo']
                
                new_metabolites.append(new_met)
                print(f"  ✓ Prepared {met_name} ({new_id}) in [{met_compartment}]")
    
    # Add all new metabolites to model
    if new_metabolites:
        model['metabolites'].extend(new_metabolites)
        print(f"Added {len(new_metabolites)} new metabolites to the model.")
    
    return len(new_metabolites)


def parse_universal_reaction(equation, model, id_database):
    """
    Parse universal reaction equation format.
    
    Format: "n A[comp1] + m B[comp2] --> k C[comp3] + l D[comp4]"
    
    Args:
        equation: Reaction equation string
        model: Model dictionary
        id_database: Metabolite ID database
    
    Returns:
        dict: {'metabolites': {met_id: stoich}} or None if parsing fails
    """
    # Split by arrow to get reactants and products
    if '-->' in equation:
        arrow = '-->'
        reversible = False
    elif '<=>' in equation:
        arrow = '<=>'
        reversible = True
    else:
        print(f"  ⚠ Warning: No valid arrow found in equation: {equation}")
        return None
    
    parts = equation.split(arrow)
    if len(parts) != 2:
        print(f"  ⚠ Warning: Invalid equation format: {equation}")
        return None
    
    reactants_str, products_str = parts
    
    # Parse metabolites
    metabolites = {}
    
    def parse_side(side_str, sign):
        """Parse one side of the equation."""
        # Split by ' + ' (with spaces) to separate terms
        # This preserves metabolites like H+ which contain + in their name
        terms = [t.strip() for t in side_str.split(' + ')]
        
        for term in terms:
            if not term:
                continue
                
            # Pattern for each term: [optional_number] metabolite_name[compartment]
            # Matches: "2 ATP[c]" or "ATP[c]" or "0.5 H+[m]" or "9 H2O[gl]"
            # Two alternatives: with stoichiometry or without
            pattern = r'^(\d+(?:\.\d+)?)\s+(.+?)\[(\w+)\]$|^(.+?)\[(\w+)\]$'
            match = re.match(pattern, term)
            
            if not match:
                continue
            
            groups = match.groups()
            if groups[0]:  # Has stoichiometry
                stoich_str = groups[0]
                met_name = groups[1]
                compartment = groups[2]
            else:  # No stoichiometry
                stoich_str = None
                met_name = groups[3]
                compartment = groups[4]
            
            # Parse stoichiometry (default to 1.0 if not provided)
            stoich = float(stoich_str) if stoich_str else 1.0
            
            # Clean metabolite name
            met_name = met_name.strip()
            compartment = compartment.strip()
            
            if not met_name or not compartment:
                continue
            
            # Find metabolite in model (with auto-creation if found in other compartments)
            met = find_metabolite_robust(model, met_name, id_database, compartment, auto_create=True)
            if not met:
                print(f"  ⚠ Warning: Metabolite '{met_name}' not found in [{compartment}]")
                continue
            
            met_id = met['id']
            metabolites[met_id] = metabolites.get(met_id, 0.0) + (sign * stoich)
    
    # Parse reactants (negative stoichiometry)
    parse_side(reactants_str, -1.0)
    
    # Parse products (positive stoichiometry)
    parse_side(products_str, 1.0)
    
    if not metabolites:
        print(f"  ⚠ Warning: No metabolites parsed from equation: {equation}")
        return None
    
    return {'metabolites': metabolites, 'reversible': reversible}


def create_compartment_reactions(model, id_database, compartment_abbrev, config, abbrev_map=None):
    """
    Create reactions for a pathway (generic version).
    
    Args:
        model: The metabolic model dictionary
        id_database: Database of metabolite IDs
        compartment_abbrev: Target compartment abbreviation
        config: Configuration dictionary containing reactions
        abbrev_map: Optional dict mapping config abbreviations to model abbreviations
                   (e.g., {'er': 'r', 'gl': 'gc'})
    
    Returns:
        int: Number of reactions created
    """
    print(f"\nCreating reactions for [{compartment_abbrev}]...")
    
    if 'reactions' not in config:
        print("  ⚠ Warning: No reactions found in config")
        return 0
    
    reactions = config['reactions']
    new_reactions = []
    
    print(f"\n  Block 1: Creating substrate transport reactions...")
    print(f"  Block 2: Creating macromolecule synthesis reactions...")
    print(f"  Block 3: Creating protein modification reactions...")
    print(f"  Block 4: Creating energy dissipation reactions...")
    print(f"  Block 5: Creating turnover and export reactions...")
    
    for rxn_config in reactions:
        # Check if reaction already exists
        rxn_id = rxn_config.get('id', get_next_reaction_id(model))
        
        exists = False
        for rxn in model['reactions']:
            if rxn['id'] == rxn_id:
                exists = True
                break
        
        if exists:
            continue
        
        # Get equation and substitute compartment abbreviations if needed
        equation = rxn_config['equation']
        if abbrev_map:
            # Replace config compartment abbreviations with model abbreviations
            for config_abbrev, model_abbrev in abbrev_map.items():
                if config_abbrev != model_abbrev:
                    equation = equation.replace(f'[{config_abbrev}]', f'[{model_abbrev}]')
        
        # Parse universal reaction format
        reaction_data = parse_universal_reaction(equation, model, id_database)
        
        if not reaction_data:
            print(f"  ⚠ Warning: Could not parse reaction: {rxn_id}")
            continue
        
        # Create reaction object
        new_rxn = {
            'id': rxn_id,
            'name': rxn_config.get('name', ''),
            'metabolites': reaction_data['metabolites'],
            'lower_bound': rxn_config.get('lower_bound', 0.0),
            'upper_bound': rxn_config.get('upper_bound', 1000.0),
            'gene_reaction_rule': rxn_config.get('gpr', ''),
            'subsystem': rxn_config.get('subsystem', ''),
        }
        
        # Add optional fields
        if 'notes' in rxn_config or 'description' in rxn_config:
            new_rxn['notes'] = {
                'description': rxn_config.get('description', rxn_config.get('notes', ''))
            }
        
        # Add annotations
        annotation = {}
        if 'ec' in rxn_config and rxn_config['ec']:
            annotation['ec-code'] = rxn_config['ec']
        if 'sbo' in rxn_config:
            annotation['sbo'] = rxn_config['sbo']
        
        if annotation:
            new_rxn['annotation'] = annotation
        
        new_reactions.append(new_rxn)
    
    # Add all new reactions to model
    if new_reactions:
        model['reactions'].extend(new_reactions)
        print(f"Added {len(new_reactions)} new reactions to the model.")
    
    return len(new_reactions)


def check_pathway_exists(model, compartment_abbrev, pathway_name=None):
    """
    Check if a pathway/compartment already exists in the model.
    
    Args:
        model: The metabolic model dictionary
        compartment_abbrev: Compartment abbreviation to check (e.g., 'ck', 'gc')
        pathway_name: Optional pathway name for more detailed reporting
    
    Returns:
        dict: Information about pathway existence
            {
                'has_compartment': bool,
                'compartment': str or None,
                'metabolite_count': int,
                'reaction_count': int,
                'metabolites': list,
                'reactions': list
            }
    """
    result = {
        'has_compartment': False,
        'compartment': None,
        'metabolite_count': 0,
        'reaction_count': 0,
        'metabolites': [],
        'reactions': []
    }
    
    # Check for compartment
    if compartment_abbrev in model.get('compartments', {}):
        result['has_compartment'] = True
        result['compartment'] = model['compartments'][compartment_abbrev]
    
    # Count metabolites in compartment
    for met in model.get('metabolites', []):
        if met.get('compartment') == compartment_abbrev:
            result['metabolites'].append(met['id'])
            result['metabolite_count'] += 1
    
    # Count reactions involving compartment
    for rxn in model.get('reactions', []):
        involves_compartment = False
        for met_id, coeff in rxn.get('metabolites', {}).items():
            # Check if metabolite is in the compartment
            for met in model['metabolites']:
                if met['id'] == met_id and met['compartment'] == compartment_abbrev:
                    involves_compartment = True
                    break
            if involves_compartment:
                break
        
        if involves_compartment:
            result['reactions'].append(rxn['id'])
            result['reaction_count'] += 1
    
    return result


# Backward compatibility: import old function names
def add_cytoskeleton_compartment(model, compartment_abbrev='ck', compartment_name='cytoskeleton'):
    """Backward compatibility wrapper for add_compartment."""
    return add_compartment(model, compartment_abbrev, compartment_name)


def create_cytoskeleton_metabolites(model, id_database, compartment_abbrev='ck', config=None):
    """Backward compatibility wrapper for create_compartment_metabolites."""
    if config is None:
        from .config import load_config
        config = load_config()
    return create_compartment_metabolites(model, id_database, compartment_abbrev, config, 'cytoskeleton_specific')


def create_cytoskeleton_reactions(model, id_database, compartment_abbrev='ck', config=None):
    """Backward compatibility wrapper for create_compartment_reactions."""
    if config is None:
        from .config import load_config
        config = load_config()
    return create_compartment_reactions(model, id_database, compartment_abbrev, config)


def check_cytoskeleton_exists(model):
    """Backward compatibility wrapper for check_pathway_exists."""
    return check_pathway_exists(model, 'ck', 'cytoskeleton')
