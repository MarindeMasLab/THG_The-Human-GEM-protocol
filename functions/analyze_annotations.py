#!/usr/bin/env python3
"""
Step 1: Analyze model annotations and extract metabolite information.
Reads target metabolites from config.json and extracts their complete annotation data.
"""

import json
from collections import defaultdict
from pathlib import Path

def analyze_model_annotations(model_path):
    """Analyze all annotation fields used in the model."""
    print("=" * 70)
    print("STEP 1: ANALYZING MODEL ANNOTATION FIELDS")
    print("=" * 70)
    
    with open(model_path, 'r') as f:
        model = json.load(f)
    
    # Count annotation fields
    annotation_fields = defaultdict(int)
    metabolites_with_annotations = 0
    
    print(f"\nAnalyzing {len(model['metabolites'])} metabolites...")
    
    for met in model['metabolites']:
        if 'annotation' in met and met['annotation']:
            metabolites_with_annotations += 1
            for key in met['annotation'].keys():
                annotation_fields[key] += 1
    
    print(f"\nMetabolites with annotations: {metabolites_with_annotations} / {len(model['metabolites'])}")
    print(f"Unique annotation field types found: {len(annotation_fields)}")
    
    print("\nAnnotation fields (sorted by usage frequency):")
    print("-" * 70)
    for field, count in sorted(annotation_fields.items(), key=lambda x: x[1], reverse=True):
        pct = (count / len(model['metabolites'])) * 100
        print(f"  {field:30s}: {count:6d} metabolites ({pct:5.1f}%)")
    
    return annotation_fields

def extract_metabolite_annotations(model_path, targets, verbose=True):
    """
    Extract complete annotation data for target metabolites.
    
    Args:
        model_path: Path to model JSON file
        targets: List of metabolite names to search for
        verbose: If True, print progress information
        
    Returns:
        dict: {metabolite_name: {annotation_data}}
    """
    if verbose:
        print("\n" + "=" * 70)
        print("EXTRACTING METABOLITE ANNOTATIONS")
        print("=" * 70)
    
    with open(model_path, 'r') as f:
        model = json.load(f)
    
    metabolite_database = {}
    not_found = []
    
    for target in targets:
        if verbose:
            print(f"\n>>> Searching for: {target}")
        
        found = False
        for met in model['metabolites']:
            # Try exact match first, then case-insensitive substring
            if (target.lower() == met['name'].lower() or 
                (target.lower() in met['name'].lower() and met['compartment'] == 'c')):
                
                if verbose:
                    print(f"  ✓ Found: {met['name']} [{met['id']}]")
                
                # Extract all annotation data
                metabolite_data = {}
                
                if 'annotation' in met and met['annotation']:
                    for key, value in met['annotation'].items():
                        # Handle list values (take first element)
                        if isinstance(value, list) and len(value) > 0:
                            metabolite_data[key] = value[0] if len(value) == 1 else value
                        else:
                            metabolite_data[key] = value
                
                # Add InChI and InChIKey if available
                if 'inchi' in met:
                    metabolite_data['inchi'] = met['inchi']
                if 'inchikey' in met:
                    metabolite_data['inchikey'] = met['inchikey']
                
                metabolite_database[target] = metabolite_data
                found = True
                break
        
        if not found:
            if verbose:
                print(f"  ✗ Not found in model")
            not_found.append(target)
    
    if verbose:
        print(f"\n" + "=" * 70)
        print(f"EXTRACTION SUMMARY")
        print(f"=" * 70)
        print(f"  Found: {len(metabolite_database)}/{len(targets)}")
        if not_found:
            print(f"  Not found: {', '.join(not_found)}")
    
    return metabolite_database, not_found

def main():
    """Main function - can be called standalone or from other scripts."""
    try:
        from functions.config import load_config, get_model_paths
    except ModuleNotFoundError:
        # Direct execution - adjust sys.path
        import sys
        from pathlib import Path
        sys.path.insert(0, str(Path(__file__).parent.parent))
        from functions.config import load_config, get_model_paths
    
    try:
        # Load config
        config = load_config()
        model_path, _ = get_model_paths(config)
        targets = config.get('metabolites', {}).get('targets', [])
        
        if not targets:
            print("⚠ Warning: No metabolite targets found in config.json")
            print("Using default sample targets...")
            targets = ['ATP', 'ADP', 'GTP', 'GDP', 'L-alanine', 'phosphate', 'H2O']
        
    except Exception as e:
        print(f"⚠ Warning: Could not load config: {e}")
        print("Using default values...")
        model_path = 'model/endoA_250917_3.json'
        targets = ['ATP', 'ADP', 'GTP', 'GDP', 'L-alanine', 'phosphate', 'H2O']
    
    # Analyze annotation fields
    print("=" * 70)
    print("STEP 1: METABOLITE ANNOTATION ANALYSIS")
    print("=" * 70)
    annotation_fields = analyze_model_annotations(model_path)
    
    # Extract metabolite annotations
    metabolite_db, not_found = extract_metabolite_annotations(model_path, targets, verbose=True)
    
    print("\n" + "=" * 70)
    print("Analysis complete.")
    print(f"Extracted {len(metabolite_db)} metabolites.")
    print("=" * 70)
    
    return metabolite_db, annotation_fields

if __name__ == '__main__':
    main()
