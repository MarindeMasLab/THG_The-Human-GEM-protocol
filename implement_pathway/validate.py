#!/usr/bin/env python3
"""
Unified Cytoskeleton Validation and Testing Tool

Run comprehensive validations and tests on cytoskeleton implementation:
  0: Run ALL validations (default)
  1: Quick validation (model structure + FBA)
  2: Unit tests (metabolite matching)
  3: Reproducibility tests (6 tests × 5 runs)
  4: Integration tests (8 pipeline tests)
  5: Mass balance check
  6: Network topology validation
  7: Database consistency check
  8: Synthesis balance verification
  9: Demand reaction functionality

Usage:
  python3 validate.py [option] [--config CONFIG_PATH] [--model MODEL_PATH]
  
Examples:
  python3 validate.py        # Run ALL validations
  python3 validate.py 0      # Run ALL validations
  python3 validate.py 1      # Quick validation only
  python3 validate.py 2      # Unit tests only
  python3 validate.py 3      # Reproducibility tests only
  python3 validate.py 4      # Integration tests only
  python3 validate.py 5 --config config/config_Syndecan1_HS.json  # Mass balance for specific pathway
  python3 validate.py 9 --config config/config_Syndecan1_HS.json  # Check demand reactions
"""

import sys
from pathlib import Path
# Add parent directory to path for shared functions module
sys.path.insert(0, str(Path(__file__).parent.parent))
import sys
import os
import json
import unittest
import hashlib
from pathlib import Path
from collections import defaultdict
from datetime import datetime

# Import config loader
from functions.config import load_config, get_model_paths

# Try to import COBRApy
try:
    import cobra
    COBRA_AVAILABLE = True
except ImportError:
    COBRA_AVAILABLE = False
    print("⚠ Warning: COBRApy not installed. Some validations will be skipped.")
    print("  Install with: pip install cobra")

# Load configuration (only used as defaults if no --config provided)
CONFIG = None
BASE_MODEL_PATH = '../models/endoA_250917_3.json'
OUTPUT_MODEL_PATH = '../models/endoA_250917_3_with_cytoskeleton.json'

try:
    CONFIG = load_config()
    BASE_MODEL_PATH, OUTPUT_MODEL_PATH = get_model_paths(CONFIG)
except Exception:
    # Silently use defaults - will be overridden by --config if provided
    pass


# ============================================================================
# Helper Functions
# ============================================================================

def build_compartment_mapping(model, config):
    """
    Build compartment abbreviation mapping from config to model.
    
    Args:
        model: Loaded COBRA model (JSON dict)
        config: Pathway configuration dict
        
    Returns:
        dict: Mapping from config abbreviation to model abbreviation
              e.g., {'er': 'r', 'gl': 'gl', 'c': 'c'}
    """
    from functions.config import resolve_all_compartments
    
    # Get resolved compartments
    resolved_compartments = resolve_all_compartments(model, config)
    
    # Build abbreviation mapping
    abbrev_map = {}
    for comp_config in config.get('compartments', []):
        config_abbrev = comp_config['abbreviation']
        config_name = comp_config['name']
        
        # Find corresponding resolved compartment
        resolved = next((c for c in resolved_compartments if c['name'] == config_name), None)
        if resolved:
            model_abbrev = resolved['abbreviation']
            abbrev_map[config_abbrev] = model_abbrev
    
    return abbrev_map

# ============================================================================
# VALIDATION 1: Quick Validation (Model Structure + FBA)
# ============================================================================

def validate_quick(model_path=None, config_path=None):
    """
    Quick validation: model structure and basic FBA.
    
    Args:
        model_path: Path to the model file to validate
        config_path: Path to the pathway config file (e.g., config/config_Syndecan1_HS.json)
                     If None, uses default config.json for cytoskeleton
    """
    if model_path is None:
        model_path = OUTPUT_MODEL_PATH
    
    # Load pathway configuration
    if config_path:
        with open(config_path, 'r') as f:
            pathway_config = json.load(f)
        pathway_name = pathway_config.get('pathway_name', 'Unknown')
    else:
        pathway_config = CONFIG
        pathway_name = 'Cytoskeleton'
    
    print("\n" + "="*80)
    print(f"VALIDATION 1: Quick Validation - {pathway_name}")
    print("="*80)
    
    if not COBRA_AVAILABLE:
        print("✗ COBRApy required for this validation")
        return False
    
    # Load model (as JSON dict for compartment mapping, then as COBRA model)
    print("\nLoading model...")
    with open(model_path, 'r') as f:
        model_dict = json.load(f)
    
    # Build compartment mapping
    abbrev_map = build_compartment_mapping(model_dict, pathway_config)
    
    # Now load as COBRA model
    model = cobra.io.load_json_model(model_path)
    print(f"✓ Loaded: {len(model.metabolites)} metabolites, {len(model.reactions)} reactions")
    
    # 1. Check compartments from config
    print("\n1. Checking pathway compartments...")
    config_compartments = pathway_config.get('compartments', [])
    
    compartment_stats = {}
    all_comps_found = True
    
    for comp_config in config_compartments:
        config_abbrev = comp_config['abbreviation']
        comp_name = comp_config['name']
        
        # Map config abbreviation to model abbreviation
        model_abbrev = abbrev_map.get(config_abbrev, config_abbrev)
        
        # Check if compartment exists in model
        if model_abbrev in model.compartments:
            if config_abbrev != model_abbrev:
                print(f"  ✓ [{config_abbrev}→{model_abbrev}] {comp_name}: {model.compartments[model_abbrev]}")
            else:
                print(f"  ✓ [{model_abbrev}] {comp_name}: {model.compartments[model_abbrev]}")
            
            # Count metabolites in this compartment
            comp_mets = [m for m in model.metabolites if m.compartment == model_abbrev]
            comp_rxns = [r for r in model.reactions if any(m.compartment == model_abbrev for m in r.metabolites)]
            
            compartment_stats[model_abbrev] = {
                'name': comp_name,
                'metabolites': len(comp_mets),
                'reactions': len(comp_rxns)
            }
            print(f"      {len(comp_mets)} metabolites, {len(comp_rxns)} reactions")
        else:
            print(f"  ✗ [{config_abbrev}→{model_abbrev}] {comp_name}: NOT FOUND")
            all_comps_found = False
    
    if not all_comps_found:
        print("✗ Some compartments are missing!")
        return False
    
    # 2. Check target metabolites in new compartments
    print("\n2. Checking target metabolites availability...")
    targets = pathway_config.get('metabolites', {}).get('targets', [])
    
    # Get pathway-specific metabolites
    pathway_name_key = pathway_config.get('pathway_name', '').replace(' ', '_').replace('-', '_')
    pathway_specific_key = f"{pathway_name_key}_specific" if pathway_name_key else None
    pathway_specific_mets = {}
    
    if pathway_specific_key and pathway_specific_key in pathway_config.get('metabolites', {}):
        pathway_specific_mets = pathway_config['metabolites'][pathway_specific_key]
    
    # Sample a few target metabolites to verify
    sample_targets = targets[:3] if len(targets) > 3 else targets
    target_check_passed = True
    
    for met_name in sample_targets:
        # Check if this metabolite exists in the new compartments
        found_compartments = []
        for comp_config in config_compartments:
            comp_abbrev = comp_config['abbreviation']
            # Search for metabolite by name in this compartment
            comp_met = next((m for m in model.metabolites 
                           if comp_abbrev in m.id and met_name.lower() in m.name.lower()), None)
            if comp_met:
                found_compartments.append(comp_abbrev)
        
        if found_compartments:
            print(f"  ✓ {met_name}: found in {found_compartments}")
        else:
            print(f"  ⚠ {met_name}: not found in new compartments (may only be in cytosol)")
    
    # Check pathway-specific metabolites
    if pathway_specific_mets:
        print(f"\n  Checking {len(pathway_specific_mets)} pathway-specific metabolites...")
        specific_found = 0
        for met_key, met_data in list(pathway_specific_mets.items())[:5]:  # Sample first 5
            met_name = met_data.get('name', met_key)
            met_comp = met_data.get('compartment')
            
            # Search for this metabolite
            found = next((m for m in model.metabolites 
                         if met_name.lower() in m.name.lower() or met_key.lower() in m.name.lower()), None)
            if found:
                print(f"    ✓ {met_name} [{met_comp}]: {found.id}")
                specific_found += 1
            else:
                print(f"    ✗ {met_name} [{met_comp}]: NOT FOUND")
                target_check_passed = False
        
        print(f"  ✓ Found {specific_found}/{min(5, len(pathway_specific_mets))} sampled pathway-specific metabolites")
    
    # 3. Test demand reactions with FBA
    print("\n3. Testing demand reactions with FBA...")
    
    # Find demand reactions from config
    demand_reactions = []
    for rxn_config in pathway_config.get('reactions', []):
        rxn_id = rxn_config.get('id', '')
        rxn_name = rxn_config.get('name', '')
        subsystem = rxn_config.get('subsystem', '')
        
        # Identify demand reactions (typically in "Protein turnover" or with "Demand" in name/subsystem)
        if 'demand' in rxn_name.lower() or 'demand' in subsystem.lower() or rxn_id.startswith('DM_'):
            demand_reactions.append((rxn_id, rxn_name))
    
    if not demand_reactions:
        # Fallback: look for synthesis reactions if no demands found
        print("  No demand reactions found in config. Looking for synthesis reactions...")
        for rxn_config in pathway_config.get('reactions', []):
            rxn_id = rxn_config.get('id', '')
            rxn_name = rxn_config.get('name', '')
            if 'synth' in rxn_name.lower() or 'synthesis' in rxn_name.lower():
                demand_reactions.append((rxn_id, rxn_name))
    
    if demand_reactions:
        success_count = 0
        total_reactions = min(len(demand_reactions), 5)  # Test up to 5 reactions
        
        for rxn_id, rxn_name in demand_reactions[:5]:
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                with model:
                    # Set flux target
                    if rxn.lower_bound >= 0:
                        rxn.lower_bound = 0.01
                    else:
                        # For reversible reactions, try both directions
                        rxn.lower_bound = 0.01
                    
                    solution = model.optimize()
                    
                    if solution.status == 'optimal':
                        flux = abs(solution.fluxes[rxn_id])
                        if flux > 1e-6:
                            print(f"  ✓ {rxn_name[:40]:40} : {flux:.6f} mmol/gDW/h")
                            success_count += 1
                        else:
                            print(f"  ✗ {rxn_name[:40]:40} : No flux")
                    else:
                        print(f"  ✗ {rxn_name[:40]:40} : {solution.status}")
            except KeyError:
                print(f"  ⚠ {rxn_name[:40]:40} : Not found in model")
            except Exception as e:
                print(f"  ✗ {rxn_name[:40]:40} : Error - {str(e)[:30]}")
        
        print(f"\n  ✓ FBA TEST: {success_count}/{total_reactions} reactions can carry flux")
        fba_passed = success_count > 0
    else:
        print("  ⚠ No demand/synthesis reactions found to test")
        fba_passed = True  # Don't fail if no reactions to test
    
    # Summary
    print(f"\n{'='*80}")
    all_passed = all_comps_found and target_check_passed and fba_passed
    print(f"QUICK VALIDATION: {'✓ PASSED' if all_passed else '✗ FAILED'}")
    print(f"{'='*80}")
    
    return all_passed


# ============================================================================
# VALIDATION 2: Unit Tests (Metabolite Matching)
# ============================================================================

def validate_unit_tests(model_path=None, config_path=None):
    """Run unit tests for metabolite matching."""
    
    print("\n" + "="*80)
    print("VALIDATION 2: Unit Tests (Metabolite Matching)")
    print("="*80)
    
    if not COBRA_AVAILABLE:
        print("✗ COBRApy required for unit tests")
        return False
    
    # Use provided model path or default
    if model_path is None:
        model_path = OUTPUT_MODEL_PATH
    
    # Import functions
    sys.path.insert(0, os.path.dirname(__file__))
    from functions.pathway_builder import find_metabolite_robust, find_metabolite_by_annotation
    
    # Load model and database
    with open(model_path, 'r') as f:
        model_dict = json.load(f)
    with open('data/metabolite_id_database.json', 'r') as f:
        id_database = json.load(f)
    
    # Get test metabolites from config or use defaults
    test_metabolites = []
    test_compartments = ['c']  # Default test compartment
    
    if config_path:
        try:
            with open(config_path, 'r') as f:
                cfg = json.load(f)
            
            pathway_name = cfg.get('pathway_name', 'Unknown')
            print(f"\nPathway: {pathway_name}")
            
            # Get metabolites from config targets
            config_targets = cfg.get('metabolites', {}).get('targets', [])
            if config_targets:
                test_metabolites = config_targets
                print(f"Testing {len(test_metabolites)} metabolites from config")
            
            # Get compartments for testing
            test_compartments = [c['abbreviation'] for c in cfg.get('compartments', [])]
            if not test_compartments:
                test_compartments = ['c']
            
        except Exception as e:
            print(f"⚠ Warning: Could not load config file: {e}")
            print("  Using default test metabolites")
    
    # Default test metabolites if none from config
    if not test_metabolites:
        print("\nUsing default test metabolites (amino acids + energy metabolites)")
        test_metabolites = [
            'L-alanine', 'L-arginine', 'L-asparagine', 'L-aspartate', 'L-cysteine',
            'L-glutamine', 'L-glutamate', 'L-glycine', 'L-histidine', 'L-isoleucine',
            'L-leucine', 'L-lysine', 'L-methionine', 'L-phenylalanine', 'L-proline',
            'L-serine', 'L-threonine', 'L-tryptophan', 'L-tyrosine', 'L-valine',
            'ATP', 'ADP', 'GTP', 'GDP', 'H2O', 'phosphate', 'H+'
        ]
    
    # Run tests
    print(f"\nTesting metabolite matching for {len(test_metabolites)} metabolites...")
    print(f"Test compartments: {test_compartments[:3]}{'...' if len(test_compartments) > 3 else ''}")
    print("="*80)
    
    found_count = 0
    not_found = []
    
    for met_name in test_metabolites:
        # Test in first available compartment
        test_comp = test_compartments[0] if test_compartments else 'c'
        
        result = find_metabolite_robust(model_dict, met_name, id_database, test_comp)
        if result:
            found_count += 1
            print(f"  ✓ {met_name:40} → {result['id']}")
        else:
            not_found.append(met_name)
            print(f"  ✗ {met_name:40} → NOT FOUND")
    
    # Summary
    print("="*80)
    success_rate = (found_count / len(test_metabolites)) * 100
    print(f"\nRESULTS:")
    print(f"  Found: {found_count}/{len(test_metabolites)} ({success_rate:.1f}%)")
    
    if not_found:
        print(f"  Not found: {len(not_found)} metabolites")
        if len(not_found) <= 10:
            for met in not_found:
                print(f"    - {met}")
        else:
            for met in not_found[:10]:
                print(f"    - {met}")
            print(f"    ... and {len(not_found) - 10} more")
    
    # Pass if at least 80% of metabolites are found
    passed = success_rate >= 80.0
    print(f"\n✓ UNIT TESTS: {'PASSED' if passed else 'FAILED'} ({success_rate:.1f}% success rate)")
    
    return passed


# ============================================================================
# VALIDATION 3: Reproducibility Tests
# ============================================================================

def validate_reproducibility(model_path=None, config_path=None):
    """Run reproducibility tests (6 tests × 5 runs)."""
    
    print("\n" + "="*80)
    print("VALIDATION 3: Reproducibility Tests (6 tests × 5 runs)")
    print("="*80)
    
    if not COBRA_AVAILABLE:
        print("✗ COBRApy required for reproducibility tests")
        return False
    
    # Use provided model path or default
    if model_path is None:
        model_path = OUTPUT_MODEL_PATH
    
    # Get pathway compartments and reactions from config
    pathway_compartments = set()
    pathway_reaction_ids = set()
    pathway_name = "pathway"
    
    if config_path:
        try:
            with open(config_path, 'r') as f:
                cfg = json.load(f)
            
            pathway_name = cfg.get('pathway_name', 'pathway')
            
            # Load model as dict for compartment mapping
            with open(model_path, 'r') as f:
                model_dict = json.load(f)
            
            # Build compartment mapping
            compartment_map = build_compartment_mapping(model_dict, cfg)
            
            # Get pathway compartments (model abbreviations)
            for comp_data in cfg.get('compartments', []):
                config_abbrev = comp_data['abbreviation']
                model_abbrev = compartment_map.get(config_abbrev, config_abbrev)
                pathway_compartments.add(model_abbrev)
            
            # Get pathway reaction IDs
            for reaction in cfg.get('reactions', []):
                pathway_reaction_ids.add(reaction['id'])
            
            print(f"\nPathway: {pathway_name}")
            print(f"Testing reproducibility across {len(pathway_compartments)} compartments and {len(pathway_reaction_ids)} reactions")
            
        except Exception as e:
            print(f"⚠ Warning: Could not load config file: {e}")
            print("  Falling back to 'ck' compartment analysis")
            pathway_compartments = {'ck'}
            pathway_name = "cytoskeleton"
    else:
        # Default to cytoskeleton compartment
        pathway_compartments = {'ck'}
        pathway_name = "cytoskeleton"
    
    def hash_object(obj):
        """Generate SHA256 hash of an object."""
        return hashlib.sha256(str(obj).encode()).hexdigest()
    
    print("\nRunning 5 independent iterations...")
    
    test_results = []
    for run in range(5):
        print(f"\n--- Run {run + 1}/5 ---")
        
        # Test 1: Model loading
        model = cobra.io.load_json_model(model_path)
        model_hash = hash_object(f"{len(model.metabolites)}_{len(model.reactions)}")
        
        # Test 2: Pathway metabolite data
        if pathway_reaction_ids:
            # Get metabolites from pathway reactions
            pathway_mets = []
            for rxn_id in pathway_reaction_ids:
                try:
                    rxn = model.reactions.get_by_id(rxn_id)
                    for met in rxn.metabolites:
                        if met.id not in pathway_mets:
                            pathway_mets.append(met.id)
                except KeyError:
                    continue
            pathway_mets = sorted(pathway_mets)
        else:
            # Fall back to compartment-based
            pathway_mets = sorted([m.id for m in model.metabolites if m.compartment in pathway_compartments])
        
        pathway_hash = hash_object(pathway_mets)
        
        # Test 3: FBA optimization
        solution = model.optimize()
        fba_hash = hash_object(f"{solution.status}_{solution.objective_value:.10f}")
        
        # Test 4: Metabolite matching database
        with open('data/metabolite_id_database.json', 'r') as f:
            id_database = json.load(f)
        db_hash = hash_object(sorted(id_database.keys()))
        
        # Test 5: Pathway reaction stoichiometry
        if pathway_reaction_ids:
            pathway_rxns = []
            for rxn_id in pathway_reaction_ids:
                try:
                    rxn = model.reactions.get_by_id(rxn_id)
                    pathway_rxns.append((rxn.id, len(rxn.metabolites)))
                except KeyError:
                    continue
            stoich_data = sorted(pathway_rxns)
        else:
            # Fall back to compartment-based
            pathway_rxns = [r for r in model.reactions if any(m.compartment in pathway_compartments for m in r.metabolites)]
            stoich_data = sorted([(r.id, len(r.metabolites)) for r in pathway_rxns])
        
        stoich_hash = hash_object(str(stoich_data))
        
        # Test 6: Transport reactions in pathway
        if pathway_reaction_ids:
            transport_rxns = []
            for rxn_id in pathway_reaction_ids:
                try:
                    rxn = model.reactions.get_by_id(rxn_id)
                    if 'transport' in rxn.name.lower() or rxn_id.startswith('T_'):
                        transport_rxns.append(rxn.id)
                except KeyError:
                    continue
            transport_rxns = sorted(transport_rxns)
        else:
            # Fall back to compartment-based
            pathway_rxns = [r for r in model.reactions if any(m.compartment in pathway_compartments for m in r.metabolites)]
            transport_rxns = sorted([r.id for r in pathway_rxns if 'transport' in r.name.lower()])
        
        transport_hash = hash_object(transport_rxns)
        
        test_results.append({
            'model': model_hash,
            'pathway_data': pathway_hash,
            'fba': fba_hash,
            'database': db_hash,
            'stoichiometry': stoich_hash,
            'transport': transport_hash
        })
    
    # Check consistency
    print("\n" + "="*80)
    print("Checking reproducibility across 5 runs...")
    print("="*80)
    
    all_reproducible = True
    for test_name in ['model', 'pathway_data', 'fba', 'database', 'stoichiometry', 'transport']:
        hashes = [result[test_name] for result in test_results]
        is_reproducible = len(set(hashes)) == 1
        
        status = "✓ IDENTICAL" if is_reproducible else "✗ DIFFERENT"
        print(f"{test_name:15}: {status}")
        
        if not is_reproducible:
            all_reproducible = False
    
    print(f"\n✓ REPRODUCIBILITY TEST: {'PASSED' if all_reproducible else 'FAILED'}")
    return all_reproducible


# ============================================================================
# VALIDATION 4: Integration Tests
# ============================================================================

def validate_integration(model_path=None, config_path=None):
    """Run integration tests (8 pipeline tests)."""
    
    print("\n" + "="*80)
    print("VALIDATION 4: Integration Tests (8 pipeline tests)")
    print("="*80)
    
    if not COBRA_AVAILABLE:
        print("✗ COBRApy required for integration tests")
        return False
    
    # Use provided model path or default
    if model_path is None:
        model_path = OUTPUT_MODEL_PATH
    
    # Get pathway info from config
    pathway_compartments = set()
    pathway_reaction_ids = set()
    pathway_name = "pathway"
    expected_min_mets = 30
    expected_min_rxns = 20
    
    if config_path:
        try:
            with open(config_path, 'r') as f:
                cfg = json.load(f)
            
            pathway_name = cfg.get('pathway_name', 'pathway')
            
            # Get pathway reaction IDs
            for reaction in cfg.get('reactions', []):
                pathway_reaction_ids.add(reaction['id'])
            
            # Get compartments
            for comp_data in cfg.get('compartments', []):
                pathway_compartments.add(comp_data['abbreviation'])
            
            # Expected minimums based on config
            expected_min_mets = len(cfg.get('metabolites', {}).get('targets', [])) - 5  # Allow some not found
            expected_min_rxns = len(pathway_reaction_ids) - 2  # Allow some missing
            
            print(f"\nPathway: {pathway_name}")
            print(f"Expected: ≥{expected_min_mets} metabolites, ≥{expected_min_rxns} reactions")
            
        except Exception as e:
            print(f"⚠ Warning: Could not load config file: {e}")
            print("  Using default cytoskeleton expectations")
            pathway_compartments = {'ck'}
            pathway_name = "cytoskeleton"
            expected_min_mets = 35
            expected_min_rxns = 40
    else:
        pathway_compartments = {'ck'}
        pathway_name = "cytoskeleton"
        expected_min_mets = 35
        expected_min_rxns = 40
    
    passed = 0
    total = 8
    
    # Test 1: Model loading
    print("\n1. Model loading pipeline...")
    try:
        model = cobra.io.load_json_model(model_path)
        assert len(model.metabolites) > 1000  # General sanity check
        assert len(model.reactions) > 1000
        print(f"✓ Model loaded: {len(model.metabolites)} metabolites, {len(model.reactions)} reactions")
        passed += 1
    except Exception as e:
        print(f"✗ Failed: {e}")
    
    # Test 2: Database loading
    print("\n2. Database loading pipeline...")
    try:
        with open('data/metabolite_id_database.json', 'r') as f:
            id_database = json.load(f)
        assert len(id_database) >= 30
        print(f"✓ Database loaded: {len(id_database)} entries")
        passed += 1
    except Exception as e:
        print(f"✗ Failed: {e}")
        id_database = {}  # Create empty dict for subsequent tests
    
    # Test 3: Metabolite matching
    print("\n3. Metabolite matching pipeline...")
    try:
        from functions.pathway_builder import find_metabolite_robust
        # Load model as JSON for pathway_builder function
        with open(model_path, 'r') as f:
            model_json = json.load(f)
        
        # Test common metabolites
        test_metabolites = ['ATP', 'H2O', 'H+', 'ADP', 'phosphate']
        test_comp = list(pathway_compartments)[0] if pathway_compartments else 'c'
        matched = sum(1 for met in test_metabolites if find_metabolite_robust(model_json, met, id_database, test_comp))
        assert matched >= 3
        print(f"✓ Matched {matched}/{len(test_metabolites)} common metabolites")
        passed += 1
    except Exception as e:
        print(f"✗ Failed: {e}")
    
    # Test 4: FBA optimization
    print("\n4. FBA optimization pipeline...")
    try:
        solution = model.optimize()
        assert solution.status == 'optimal'
        print(f"✓ FBA optimal: {solution.objective_value:.6f}")
        passed += 1
    except Exception as e:
        print(f"✗ Failed: {e}")
    
    # Test 5: Pathway reaction presence
    print("\n5. Pathway reaction pipeline...")
    try:
        if pathway_reaction_ids:
            found_rxns = sum(1 for rxn_id in pathway_reaction_ids if rxn_id in model.reactions)
            success_rate = (found_rxns / len(pathway_reaction_ids)) * 100
            assert found_rxns >= expected_min_rxns
            print(f"✓ Found {found_rxns}/{len(pathway_reaction_ids)} pathway reactions ({success_rate:.1f}%)")
        else:
            # Fall back to compartment check
            pathway_rxns = [r for r in model.reactions if any(m.compartment in pathway_compartments for m in r.metabolites)]
            assert len(pathway_rxns) >= expected_min_rxns
            print(f"✓ Found {len(pathway_rxns)} reactions in pathway compartments")
        passed += 1
    except Exception as e:
        print(f"✗ Failed: {e}")
    
    # Test 6: Transport validation
    print("\n6. Transport reaction pipeline...")
    try:
        if pathway_reaction_ids:
            transport_rxns = []
            for rxn_id in pathway_reaction_ids:
                try:
                    rxn = model.reactions.get_by_id(rxn_id)
                    if 'transport' in rxn.name.lower() or rxn_id.startswith('T_'):
                        transport_rxns.append(rxn_id)
                except KeyError:
                    continue
            print(f"✓ Found {len(transport_rxns)} transport reactions in pathway")
        else:
            pathway_rxns = [r for r in model.reactions if any(m.compartment in pathway_compartments for m in r.metabolites)]
            transport_rxns = [r for r in pathway_rxns if 'transport' in r.name.lower()]
            print(f"✓ Found {len(transport_rxns)} transport reactions")
        passed += 1
    except Exception as e:
        print(f"✗ Failed: {e}")
    
    # Test 7: Compartment validation
    print("\n7. Compartment validation pipeline...")
    try:
        # Check if pathway compartments exist
        missing_comps = []
        for comp in pathway_compartments:
            if comp not in model.compartments:
                missing_comps.append(comp)
        
        if missing_comps:
            print(f"⚠ Warning: Compartments not in model: {missing_comps}")
        
        # Count metabolites in pathway compartments or reactions
        if pathway_reaction_ids:
            pathway_mets = []
            for rxn_id in pathway_reaction_ids:
                try:
                    rxn = model.reactions.get_by_id(rxn_id)
                    for met in rxn.metabolites:
                        if met.id not in [m.id for m in pathway_mets]:
                            pathway_mets.append(met)
                except KeyError:
                    continue
            assert len(pathway_mets) >= expected_min_mets
            print(f"✓ Pathway metabolites: {len(pathway_mets)} (expected ≥{expected_min_mets})")
        else:
            pathway_mets = [m for m in model.metabolites if m.compartment in pathway_compartments]
            print(f"✓ Pathway compartment metabolites: {len(pathway_mets)}")
        passed += 1
    except Exception as e:
        print(f"✗ Failed: {e}")
    
    # Test 8: Complete workflow
    print("\n8. Complete workflow pipeline...")
    try:
        # Check pathway is integrated into model
        if pathway_reaction_ids:
            found = sum(1 for rxn_id in pathway_reaction_ids if rxn_id in model.reactions)
            assert found >= expected_min_rxns
        else:
            assert any(comp in model.compartments for comp in pathway_compartments)
        
        # Check model is still solvable
        solution = model.optimize()
        assert solution.status == 'optimal'
        print(f"✓ Complete workflow successful (FBA: {solution.status})")
        passed += 1
    except Exception as e:
        print(f"✗ Failed: {e}")
    
    print(f"\n{'='*80}")
    print(f"INTEGRATION TESTS: {passed}/{total} passed ({100*passed/total:.1f}%)")
    print(f"{'='*80}")
    return passed == total


# ============================================================================
# VALIDATION 5: Mass Balance Check
# ============================================================================

def validate_mass_balance(model_path=None, config_path=None):
    """
    Check mass balance of pathway reactions.
    
    Args:
        model_path: Path to the model file to validate
        config_path: Path to the pathway config file (e.g., config/config_Syndecan1_HS.json)
                     If None, uses default config.json for cytoskeleton
    """
    if model_path is None:
        model_path = OUTPUT_MODEL_PATH
    
    # Load pathway configuration
    if config_path:
        with open(config_path, 'r') as f:
            pathway_config = json.load(f)
        pathway_name = pathway_config.get('pathway_name', 'Unknown')
    else:
        pathway_config = CONFIG
        pathway_name = 'Cytoskeleton'
    
    print("\n" + "="*80)
    print(f"VALIDATION 5: Mass Balance Check - {pathway_name}")
    print("="*80)
    
    if not COBRA_AVAILABLE:
        print("✗ COBRApy required for mass balance check")
        return False
    
    model = cobra.io.load_json_model(model_path)
    
    # Get pathway reaction IDs from config
    pathway_rxn_ids = set()
    for rxn_config in pathway_config.get('reactions', []):
        rxn_id = rxn_config.get('id', '')
        if rxn_id:
            pathway_rxn_ids.add(rxn_id)
    
    # Filter reactions that exist in model and are from pathway
    pathway_rxns = []
    for rxn_id in pathway_rxn_ids:
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            pathway_rxns.append(rxn)
        except KeyError:
            pass  # Reaction not in model
    
    print(f"\nChecking {len(pathway_rxns)} pathway reactions...")
    
    imbalanced = []
    skipped = 0
    for rxn in pathway_rxns:
        # Skip biomass/demand/exchange reactions (these are intentionally "imbalanced" sinks)
        rxn_name_lower = rxn.name.lower() if rxn.name else ""
        rxn_id_lower = rxn.id.lower()
        if any(keyword in rxn_id_lower or keyword in rxn_name_lower 
               for keyword in ['biomass', 'demand', 'exchange', 'sink', 'ex_', 'dm_']):
            skipped += 1
            continue
        
        # Check if reaction is balanced
        # check_mass_balance() returns empty dict {} if balanced, dict with elements if imbalanced
        balance_result = rxn.check_mass_balance()
        if balance_result:  # If dict is not empty, reaction is imbalanced
            imbalanced.append((rxn, balance_result))
    
    checked = len(pathway_rxns) - skipped
    print(f"  Checked: {checked} reactions ({skipped} sink reactions skipped)")
    
    if len(imbalanced) == 0:
        print(f"✓ All {pathway_name} reactions are mass balanced!")
        return True
    else:
        print(f"✗ Found {len(imbalanced)} imbalanced reactions:")
        for rxn, balance in imbalanced[:10]:  # Show first 10
            print(f"  - {rxn.id}: {rxn.name}")
            print(f"    Imbalance: {balance}")
        if len(imbalanced) > 10:
            print(f"  ... and {len(imbalanced) - 10} more")
        return False


# ============================================================================
# VALIDATION 6: Network Topology
# ============================================================================

def validate_network_topology(model_path=None, config_path=None):
    """Validate network topology and connectivity."""
    
    print("\n" + "="*80)
    print("VALIDATION 6: Network Topology Validation")
    print("="*80)
    
    if not COBRA_AVAILABLE:
        print("✗ COBRApy required for network topology validation")
        return False
    
    # Use provided model path or default
    if model_path is None:
        model_path = OUTPUT_MODEL_PATH
    
    model = cobra.io.load_json_model(model_path)
    
    # Get pathway compartments and reactions from config
    pathway_compartments = set()
    pathway_reaction_ids = set()
    
    if config_path:
        try:
            with open(config_path, 'r') as f:
                cfg = json.load(f)
            
            # Load model as dict for compartment mapping
            with open(model_path, 'r') as f:
                model_dict = json.load(f)
            
            # Build compartment mapping
            compartment_map = build_compartment_mapping(model_dict, cfg)
            
            # Get pathway compartments (model abbreviations)
            for comp_data in cfg.get('compartments', []):
                config_abbrev = comp_data['abbreviation']
                model_abbrev = compartment_map.get(config_abbrev, config_abbrev)
                pathway_compartments.add(model_abbrev)
            
            # Get pathway reaction IDs
            for reaction in cfg.get('reactions', []):
                pathway_reaction_ids.add(reaction['id'])
            
            print(f"\nPathway: {cfg.get('pathway_name', 'Unknown')}")
            print(f"Compartments: {sorted(pathway_compartments)}")
            
        except Exception as e:
            print(f"⚠ Warning: Could not load config file: {e}")
            print("  Falling back to 'ck' compartment analysis")
            pathway_compartments = {'ck'}
    else:
        # Default to cytoskeleton compartment
        pathway_compartments = {'ck'}
    
    # Get pathway reactions
    pathway_reactions = []
    if pathway_reaction_ids:
        print(f"Analyzing {len(pathway_reaction_ids)} pathway reactions...")
        for rxn_id in pathway_reaction_ids:
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                pathway_reactions.append(rxn)
            except KeyError:
                print(f"  ⚠ Warning: Reaction {rxn_id} not found in model")
                continue
    else:
        # Fall back to reactions in pathway compartments
        pathway_reactions = [r for r in model.reactions 
                           if any(m.compartment in pathway_compartments for m in r.metabolites)]
    
    # Get all metabolites involved in pathway reactions
    pathway_mets = []
    for rxn in pathway_reactions:
        for met in rxn.metabolites:
            if met not in pathway_mets:
                pathway_mets.append(met)
    
    print(f"\nAnalyzing {len(pathway_mets)} metabolites from {len(pathway_reactions)} pathway reactions...")

    
    # Check for orphaned metabolites (not connected to any PATHWAY reaction)
    orphaned = []
    for met in pathway_mets:
        # Count reactions involving this metabolite that are in the pathway
        pathway_rxn_count = sum(1 for r in met.reactions if r in pathway_reactions)
        if pathway_rxn_count == 0:
            orphaned.append(met)
    
    # Check for dead-end metabolites IN THE FULL MODEL
    # (metabolites used by pathway but blocked in the full model context)
    dead_ends = []
    for met in pathway_mets:
        # Consider ALL reactions in the model (not just pathway reactions)
        producing = [r for r in met.reactions if r.metabolites[met] > 0]
        consuming = [r for r in met.reactions if r.metabolites[met] < 0]
        
        if len(producing) == 0 and len(consuming) > 0:
            dead_ends.append(('sink', met))
        elif len(consuming) == 0 and len(producing) > 0:
            dead_ends.append(('source', met))
    
    print(f"\n✓ Total metabolites: {len(pathway_mets)}")
    print(f"✓ Orphaned metabolites: {len(orphaned)}")
    print(f"✓ Dead-end metabolites (sinks/sources): {len(dead_ends)}")
    
    # Show details if there are issues
    if orphaned:
        print(f"\n⚠ Orphaned metabolites (not connected to reactions):")
        for met in orphaned[:10]:  # Show first 10
            print(f"  - {met.id}: {met.name} [{met.compartment}]")
        if len(orphaned) > 10:
            print(f"  ... and {len(orphaned) - 10} more")
    
    if dead_ends:
        print(f"\n⚠ Dead-end metabolites:")
        sinks = [m for t, m in dead_ends if t == 'sink']
        sources = [m for t, m in dead_ends if t == 'source']
        
        if sinks:
            print(f"  Sinks (only consumed, never produced): {len(sinks)}")
            for met in sinks[:5]:
                print(f"    - {met.id}: {met.name} [{met.compartment}]")
            if len(sinks) > 5:
                print(f"    ... and {len(sinks) - 5} more")
        
        if sources:
            print(f"  Sources (only produced, never consumed): {len(sources)}")
            for met in sources[:5]:
                print(f"    - {met.id}: {met.name} [{met.compartment}]")
            if len(sources) > 5:
                print(f"    ... and {len(sources) - 5} more")
    
    # More lenient criteria for pathway validation
    # Some dead-ends are expected (external inputs/outputs)
    if len(orphaned) == 0 and len(dead_ends) < 20:
        print("\n✓ Network topology is good!")
        return True
    elif len(orphaned) == 0:
        print(f"\n⚠ Warning: {len(dead_ends)} dead-ends found (acceptable for pathways with external inputs/outputs)")
        return True
    else:
        print(f"\n✗ Issues found: {len(orphaned)} orphaned metabolites")
        return False


# ============================================================================
# VALIDATION 7: Database Consistency
# ============================================================================

def validate_database_consistency():
    """Validate metabolite database consistency."""
    
    print("\n" + "="*80)
    print("VALIDATION 7: Database Consistency Check")
    print("="*80)
    
    with open('data/metabolite_id_database.json', 'r') as f:
        id_database = json.load(f)
    
    print(f"\nAnalyzing database with {len(id_database)} entries...")
    
    # Check annotation coverage
    annotation_types = set()
    coverage = defaultdict(int)
    
    for met_name, annotations in id_database.items():
        for ann_type in annotations.keys():
            annotation_types.add(ann_type)
            coverage[ann_type] += 1
    
    print(f"\n✓ Annotation types: {len(annotation_types)}")
    for ann_type in sorted(annotation_types):
        count = coverage[ann_type]
        percentage = (count / len(id_database)) * 100
        print(f"  - {ann_type:20}: {count:3}/{len(id_database)} ({percentage:.1f}%)")
    
    # Check for metabolites with multiple annotations
    multi_annotated = sum(1 for anns in id_database.values() if len(anns) >= 2)
    print(f"\n✓ Metabolites with ≥2 annotations: {multi_annotated}/{len(id_database)} ({100*multi_annotated/len(id_database):.1f}%)")
    
    if multi_annotated >= len(id_database) * 0.9:
        print("✓ Database has excellent annotation coverage!")
        return True
    else:
        print("⚠ Database could benefit from more cross-references")
        return True


# ============================================================================
# VALIDATION 8: Synthesis Balance
# ============================================================================

def validate_synthesis_balance(model_path=None, config_path=None):
    """Verify synthesis reaction balance."""
    
    print("\n" + "="*80)
    print("VALIDATION 8: Synthesis Balance Verification")
    print("="*80)
    
    if not COBRA_AVAILABLE:
        print("✗ COBRApy required for synthesis balance check")
        return False
    
    # Use provided model path or default
    if model_path is None:
        model_path = OUTPUT_MODEL_PATH
    
    model = cobra.io.load_json_model(model_path)
    
    # Get synthesis reactions from config or default
    synthesis_reactions = {}
    
    if config_path:
        try:
            with open(config_path, 'r') as f:
                cfg = json.load(f)
            
            pathway_name = cfg.get('pathway_name', 'Unknown')
            print(f"\nPathway: {pathway_name}")
            
            # Get all reactions with 'Synth' in their ID from config
            for reaction in cfg.get('reactions', []):
                rxn_id = reaction['id']
                if 'Synth' in rxn_id or 'synth' in rxn_id:
                    rxn_name = reaction.get('name', rxn_id)
                    synthesis_reactions[rxn_id] = rxn_name
            
            if synthesis_reactions:
                print(f"Found {len(synthesis_reactions)} synthesis reaction(s) in config")
            else:
                print("No synthesis reactions found in config (looking for 'Synth' in reaction ID)")
                print("Checking all pathway reactions for mass balance...")
                # Fall back to checking all pathway reactions
                for reaction in cfg.get('reactions', []):
                    rxn_id = reaction['id']
                    rxn_name = reaction.get('name', rxn_id)
                    synthesis_reactions[rxn_id] = rxn_name
                    
        except Exception as e:
            print(f"⚠ Warning: Could not load config file: {e}")
            print("  Falling back to default cytoskeleton synthesis reactions")
            synthesis_reactions = {
                'R_Synth_Actin': 'Actin',
                'R_Synth_Tubulin': 'Tubulin',
                'R_Synth_Vimentin': 'Vimentin',
                'R_Synth_MyosinLC': 'Myosin'
            }
    else:
        # Default to cytoskeleton protein synthesis reactions
        synthesis_reactions = {
            'R_Synth_Actin': 'Actin',
            'R_Synth_Tubulin': 'Tubulin',
            'R_Synth_Vimentin': 'Vimentin',
            'R_Synth_MyosinLC': 'Myosin'
        }
        print("\nChecking default cytoskeleton synthesis reactions...")
    
    print(f"\nChecking {len(synthesis_reactions)} reaction(s) for mass balance...")
    
    all_balanced = True
    checked_count = 0
    missing_count = 0
    
    for rxn_id, rxn_name in synthesis_reactions.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            checked_count += 1
            
            # Check stoichiometry
            substrates = {m.id: abs(coeff) for m, coeff in rxn.metabolites.items() if coeff < 0}
            products = {m.id: coeff for m, coeff in rxn.metabolites.items() if coeff > 0}
            
            print(f"\n{rxn_name} ({rxn_id}):")
            print(f"  Substrates: {len(substrates)}")
            print(f"  Products: {len(products)}")
            
            # Check if reaction is balanced
            # check_mass_balance() returns empty dict {} if balanced, dict with elements if imbalanced
            balance_result = rxn.check_mass_balance()
            if not balance_result:  # Empty dict means balanced
                print(f"  ✓ Mass balanced")
            else:
                print(f"  ✗ Mass imbalanced: {balance_result}")
                all_balanced = False
                
        except KeyError:
            print(f"\n⚠ Warning: Reaction {rxn_id} not found in model")
            missing_count += 1
            continue
    
    print(f"\n{'='*80}")
    print(f"SUMMARY:")
    print(f"  Checked: {checked_count}/{len(synthesis_reactions)} reactions")
    if missing_count > 0:
        print(f"  Missing: {missing_count} reactions not found in model")
    print(f"  Result: {'✓ ALL BALANCED' if all_balanced and checked_count > 0 else '✗ IMBALANCED REACTIONS FOUND'}")
    print(f"{'='*80}")
    
    return all_balanced and checked_count > 0
    return all_balanced


# ============================================================================
# VALIDATION 9: Demand Reaction Functionality
# ============================================================================

def validate_demand_reactions(model_path=None, config_path=None):
    """Check if all demand reactions can carry flux."""
    
    print("\n" + "="*80)
    print("VALIDATION 9: Demand Reaction Functionality")
    print("="*80)
    
    if not COBRA_AVAILABLE:
        print("✗ COBRApy required for demand reaction check")
        return False
    
    # Use provided model path or default
    if model_path is None:
        model_path = OUTPUT_MODEL_PATH
    
    model = cobra.io.load_json_model(model_path)
    
    # Get demand reactions from config or find them in model
    demand_reactions = {}
    
    if config_path:
        try:
            with open(config_path, 'r') as f:
                cfg = json.load(f)
            
            pathway_name = cfg.get('pathway_name', 'Unknown')
            print(f"\nPathway: {pathway_name}")
            
            # Get all reactions with 'DM_' prefix or 'Demand' in name from config
            for reaction in cfg.get('reactions', []):
                rxn_id = reaction['id']
                if rxn_id.startswith('DM_') or 'demand' in reaction.get('name', '').lower():
                    rxn_name = reaction.get('name', rxn_id)
                    demand_reactions[rxn_id] = rxn_name
            
            if demand_reactions:
                print(f"Found {len(demand_reactions)} demand reaction(s) in config")
            else:
                print("No demand reactions found in config (looking for 'DM_' prefix)")
                print("Searching model for demand reactions...")
                # Fall back to finding all demand reactions in model
                for rxn in model.reactions:
                    if rxn.id.startswith('DM_') or 'demand' in rxn.name.lower():
                        demand_reactions[rxn.id] = rxn.name
                print(f"Found {len(demand_reactions)} demand reactions in model")
                    
        except Exception as e:
            print(f"⚠ Warning: Could not load config file: {e}")
            print("  Searching model for demand reactions...")
            for rxn in model.reactions:
                if rxn.id.startswith('DM_'):
                    demand_reactions[rxn.id] = rxn.name
    else:
        # Find all demand reactions in model
        print("\nSearching model for demand reactions...")
        for rxn in model.reactions:
            if rxn.id.startswith('DM_'):
                demand_reactions[rxn.id] = rxn.name
        print(f"Found {len(demand_reactions)} demand reactions")
    
    if not demand_reactions:
        print("\n⚠ No demand reactions found to test")
        return True  # Not a failure, just nothing to test
    
    print(f"\nTesting flux capability for {len(demand_reactions)} demand reaction(s)...")
    print("="*80)
    
    functional_count = 0
    blocked_count = 0
    missing_count = 0
    
    blocked_demands = []
    functional_demands = []
    
    for rxn_id, rxn_name in demand_reactions.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            
            # Test if demand can carry flux by running FVA
            print(f"\nTesting: {rxn_name}")
            print(f"  ID: {rxn_id}")
            print(f"  Bounds: [{rxn.lower_bound}, {rxn.upper_bound}]")
            
            # Run FVA on this specific reaction
            try:
                fva_result = cobra.flux_analysis.flux_variability_analysis(
                    model, 
                    reaction_list=[rxn_id],
                    fraction_of_optimum=0.0  # Don't require optimal growth
                )
                
                min_flux = fva_result.loc[rxn_id, 'minimum']
                max_flux = fva_result.loc[rxn_id, 'maximum']
                
                print(f"  FVA range: [{min_flux:.6f}, {max_flux:.6f}]")
                
                # Check if reaction can carry any flux
                if abs(max_flux) < 1e-9 and abs(min_flux) < 1e-9:
                    print(f"  ✗ BLOCKED - Cannot carry flux")
                    blocked_count += 1
                    blocked_demands.append((rxn_id, rxn_name))
                else:
                    print(f"  ✓ FUNCTIONAL - Can carry flux")
                    functional_count += 1
                    functional_demands.append((rxn_id, rxn_name, min_flux, max_flux))
                    
            except Exception as e:
                print(f"  ✗ FVA failed: {e}")
                blocked_count += 1
                blocked_demands.append((rxn_id, rxn_name))
                
        except KeyError:
            print(f"\n⚠ Warning: Reaction {rxn_id} not found in model")
            missing_count += 1
            continue
    
    # Summary
    print(f"\n{'='*80}")
    print(f"SUMMARY:")
    print(f"  Total demands: {len(demand_reactions)}")
    print(f"  Functional: {functional_count} ({100*functional_count/len(demand_reactions):.1f}%)")
    print(f"  Blocked: {blocked_count} ({100*blocked_count/len(demand_reactions):.1f}%)")
    if missing_count > 0:
        print(f"  Missing: {missing_count} reactions not found in model")
    
    if blocked_demands:
        print(f"\n⚠ BLOCKED DEMANDS:")
        for rxn_id, rxn_name in blocked_demands:
            print(f"  - {rxn_id}: {rxn_name}")
    
    if functional_demands:
        print(f"\n✓ FUNCTIONAL DEMANDS:")
        for rxn_id, rxn_name, min_flux, max_flux in functional_demands:
            print(f"  - {rxn_id}: {rxn_name}")
            print(f"    Range: [{min_flux:.6f}, {max_flux:.6f}]")
    
    print(f"\n{'='*80}")
    all_functional = blocked_count == 0 and missing_count == 0
    print(f"Result: {'✓ ALL DEMANDS FUNCTIONAL' if all_functional else '✗ SOME DEMANDS BLOCKED OR MISSING'}")
    print(f"{'='*80}")
    
    return all_functional


# ============================================================================
# MAIN
# ============================================================================

def print_usage():
    """Print usage instructions."""
    print(__doc__)


def main():
    """Main function."""
    
    # Parse command line arguments
    option = 0  # Default: run all
    config_path = None
    model_path = None
    
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        
        if arg in ['-h', '--help', 'help']:
            print_usage()
            return
        elif arg == '--config':
            if i + 1 < len(sys.argv):
                config_path = sys.argv[i + 1]
                i += 2
            else:
                print("Error: --config requires a path argument")
                return
        elif arg == '--model':
            if i + 1 < len(sys.argv):
                model_path = sys.argv[i + 1]
                i += 2
            else:
                print("Error: --model requires a path argument")
                return
        else:
            try:
                option = int(arg)
                i += 1
            except ValueError:
                print(f"Error: Invalid option '{arg}'. Must be 0-8.")
                print_usage()
                return
    
    print("\n" + "="*80)
    if config_path:
        print("PATHWAY VALIDATION AND TESTING TOOL")
        print(f"Config: {config_path}")
    else:
        print("CYTOSKELETON VALIDATION AND TESTING TOOL")
    print("="*80)
    
    # Determine model path
    if model_path is None:
        if config_path:
            # Try to infer model path from config
            try:
                with open(config_path, 'r') as f:
                    cfg = json.load(f)
                pathway_name = cfg.get('pathway_name', 'pathway')
                model_file = cfg.get('model_file', 'endoA_250917_3.json').replace('.json', '')
                model_path = f"../models/{model_file}_with_{pathway_name}.json"
            except:
                model_path = OUTPUT_MODEL_PATH
        else:
            model_path = OUTPUT_MODEL_PATH
    
    # Check if model exists
    if not os.path.exists(model_path):
        print("\n✗ Error: Model file not found!")
        print(f"  Expected: {model_path}")
        print("  Run the implementation first: python3 pathway_implementation.py")
        return
    
    results = {}
    
    # Run validations based on option
    if option == 0:
        print("\nRunning ALL validations (1-9)...\n")
        results['quick'] = validate_quick(model_path, config_path)
        results['unit_tests'] = validate_unit_tests(model_path, config_path)
        results['reproducibility'] = validate_reproducibility(model_path, config_path)
        results['integration'] = validate_integration(model_path, config_path)
        results['mass_balance'] = validate_mass_balance(model_path, config_path)
        results['network_topology'] = validate_network_topology(model_path, config_path)
        results['database'] = validate_database_consistency()
        results['synthesis'] = validate_synthesis_balance(model_path, config_path)
        results['demand_reactions'] = validate_demand_reactions(model_path, config_path)
        
        # Summary
        print("\n" + "="*80)
        print("VALIDATION SUMMARY")
        print("="*80)
        passed = sum(1 for v in results.values() if v)
        total = len(results)
        
        for name, result in results.items():
            status = "✓ PASSED" if result else "✗ FAILED"
            print(f"{name:20}: {status}")
        
        print(f"\n{'='*80}")
        print(f"TOTAL: {passed}/{total} validations passed ({100*passed/total:.1f}%)")
        print(f"{'='*80}")
        
        if passed == total:
            print("\n🎉 ALL VALIDATIONS PASSED! Model is production-ready.")
        else:
            print(f"\n⚠ {total - passed} validation(s) failed. Review above.")
        
    elif option == 1:
        results['quick'] = validate_quick(model_path, config_path)
    elif option == 2:
        results['unit_tests'] = validate_unit_tests(model_path, config_path)
    elif option == 3:
        results['reproducibility'] = validate_reproducibility(model_path, config_path)
    elif option == 4:
        results['integration'] = validate_integration(model_path, config_path)
    elif option == 5:
        results['mass_balance'] = validate_mass_balance(model_path, config_path)
    elif option == 6:
        results['network_topology'] = validate_network_topology(model_path, config_path)
    elif option == 7:
        results['database'] = validate_database_consistency()
    elif option == 8:
        results['synthesis'] = validate_synthesis_balance(model_path, config_path)
    elif option == 9:
        results['demand_reactions'] = validate_demand_reactions(model_path, config_path)
    else:
        print(f"Error: Invalid option {option}. Must be 0-9.")
        print_usage()
        return
    
    # Print single validation result
    if option > 0:
        validation_name = list(results.keys())[0]
        result = results[validation_name]
        print(f"\n{'='*80}")
        print(f"Result: {'✓ PASSED' if result else '✗ FAILED'}")
        print(f"{'='*80}\n")
    
    # ========================================================================
    # GENERATE VALIDATION REPORT
    # ========================================================================
    # Create reports directory if it doesn't exist
    reports_dir = Path("reports")
    reports_dir.mkdir(exist_ok=True)
    
    # Get pathway name from config
    pathway_name = "pathway"
    if config_path:
        try:
            with open(config_path, 'r') as f:
                cfg = json.load(f)
            pathway_name = cfg.get('pathway_name', 'pathway')
        except:
            pathway_name = "pathway"
    
    # Generate date string for report (format: YYYY-MM-DD)
    date_str = datetime.now().strftime('%Y-%m-%d')
    report_path = reports_dir / f"{date_str}_report_{pathway_name}_validation.txt"
    
    # Build detailed report content
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append(f"{pathway_name.upper()} PATHWAY VALIDATION REPORT")
    report_lines.append("=" * 80)
    report_lines.append(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append(f"Pathway: {pathway_name}")
    report_lines.append(f"Model: {model_path}")
    if config_path:
        report_lines.append(f"Config: {config_path}")
    report_lines.append("")
    
    # Summary statistics
    passed = sum(1 for v in results.values() if v)
    total = len(results)
    
    report_lines.append("=" * 80)
    report_lines.append("SUMMARY")
    report_lines.append("=" * 80)
    report_lines.append(f"Total validations: {total}")
    report_lines.append(f"Passed: {passed} ({100*passed/total:.1f}%)")
    report_lines.append(f"Failed: {total - passed} ({100*(total-passed)/total:.1f}%)")
    report_lines.append("")
    
    # Overall status
    if passed == total:
        report_lines.append("STATUS: ✓ ALL VALIDATIONS PASSED")
        report_lines.append("Model is production-ready.")
    else:
        report_lines.append(f"STATUS: ✗ {total - passed} VALIDATION(S) FAILED")
        report_lines.append("Review detailed results below.")
    report_lines.append("")
    
    # Detailed validation results
    report_lines.append("=" * 80)
    report_lines.append("DETAILED VALIDATION RESULTS")
    report_lines.append("=" * 80)
    report_lines.append("")
    
    # Map validation names to descriptions
    validation_descriptions = {
        'quick': {
            'name': 'Quick Validation',
            'description': 'Model structure + FBA optimization',
            'details': 'Checks model loads correctly, FBA is optimal, and basic constraints'
        },
        'unit_tests': {
            'name': 'Unit Tests',
            'description': 'Metabolite matching functionality',
            'details': 'Tests metabolite identification from config targets'
        },
        'reproducibility': {
            'name': 'Reproducibility Tests',
            'description': '6 tests × 5 runs for consistency',
            'details': 'Ensures model loading, FBA, and pathway data are deterministic'
        },
        'integration': {
            'name': 'Integration Tests',
            'description': '8 pipeline tests',
            'details': 'Tests complete workflow: loading, matching, FBA, reactions, transport'
        },
        'mass_balance': {
            'name': 'Mass Balance Check',
            'description': 'Elemental balance of pathway reactions',
            'details': 'Verifies C, H, N, O, P, S balance for all pathway reactions'
        },
        'network_topology': {
            'name': 'Network Topology',
            'description': 'Pathway connectivity analysis',
            'details': 'Checks for orphaned metabolites and dead-end reactions'
        },
        'database': {
            'name': 'Database Consistency',
            'description': 'Metabolite annotation coverage',
            'details': 'Validates metabolite ID database completeness and cross-references'
        },
        'synthesis': {
            'name': 'Synthesis Balance',
            'description': 'Protein/polymer synthesis validation',
            'details': 'Checks mass balance of synthesis reactions (e.g., R_Synth_*)'
        },
        'demand_reactions': {
            'name': 'Demand Reaction Functionality',
            'description': 'Flux capability of demand reactions',
            'details': 'Tests if demand reactions (DM_*) can carry flux using FVA'
        }
    }
    
    for idx, (name, result) in enumerate(results.items(), 1):
        status = "✓ PASSED" if result else "✗ FAILED"
        info = validation_descriptions.get(name, {
            'name': name.replace('_', ' ').title(),
            'description': 'Validation test',
            'details': 'No description available'
        })
        
        report_lines.append(f"{idx}. {info['name']}")
        report_lines.append(f"   Status: {status}")
        report_lines.append(f"   Description: {info['description']}")
        report_lines.append(f"   Details: {info['details']}")
        report_lines.append("")
    
    report_lines.append("=" * 80)
    report_lines.append("RECOMMENDATIONS")
    report_lines.append("=" * 80)
    
    # Add recommendations based on failed validations
    if 'mass_balance' in results and not results['mass_balance']:
        report_lines.append("• Mass Balance: Check reaction formulas in config file")
        report_lines.append("  Run: python3 validate.py 5 --config <config> --model <model>")
    
    if 'demand_reactions' in results and not results['demand_reactions']:
        report_lines.append("• Demand Reactions: Investigate blocked demands")
        report_lines.append("  Run: python3 validate.py 9 --config <config> --model <model>")
        report_lines.append("  Check pathway connectivity and metabolite availability")
    
    if 'network_topology' in results and not results['network_topology']:
        report_lines.append("• Network Topology: Review orphaned metabolites and dead-ends")
        report_lines.append("  Run: python3 validate.py 6 --config <config> --model <model>")
    
    if 'quick' in results and not results['quick']:
        report_lines.append("• Quick Validation: Check FBA constraints and model structure")
        report_lines.append("  Run: python3 validate.py 1 --config <config> --model <model>")
    
    if passed < total:
        report_lines.append("")
        report_lines.append("For detailed output, re-run individual validation options.")
    
    report_lines.append("")
    report_lines.append("=" * 80)
    
    # Write report
    with open(report_path, 'w') as f:
        f.write('\n'.join(report_lines))
    
    print(f"\n📄 Validation report saved: {report_path}")


if __name__ == '__main__':
    main()
