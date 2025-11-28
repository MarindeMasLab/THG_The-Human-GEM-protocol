#!/usr/bin/env python3
"""
Main script for Pathway Implementation into Genome-Scale Metabolic Models
Orchestrates the entire process from annotation analysis to model modification

Usage:
    python3 pathway_implementation.py <model_name>
    
    where <model_name> is the base name of your model file (without .json extension)
    Example: python3 pathway_implementation.py endoA_250917_3
"""

import sys
from pathlib import Path
# Add parent directory to path for shared functions module
sys.path.insert(0, str(Path(__file__).parent.parent))
import os
import sys
import json
import argparse
from pathlib import Path
from datetime import datetime

# Check Python version
if sys.version_info < (3, 6):
    print("❌ ERROR: This script requires Python 3.6 or higher")
    print(f"   Current version: Python {sys.version_info.major}.{sys.version_info.minor}")
    sys.exit(1)

# Import the individual modules from functions/
from functions import analyze_annotations
from functions import build_id_database
from functions import pathway_builder
from functions.config import load_config, get_model_paths


def check_python_requirements():
    """Check Python version and display system information."""
    print("System Information:")
    print(f"  Python version: {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}")
    print(f"  Platform: {sys.platform}")
    print(f"  Working directory: {os.getcwd()}")
    
    # Check minimum version
    if sys.version_info < (3, 6):
        print("\n❌ ERROR: Python 3.6 or higher is required")
        return False
    
    print("  ✓ Python version OK")
    return True


def check_and_create_directories():
    """Ensure required directories exist."""
    base_path = Path(__file__).parent
    reports_dir = base_path / "reports"
    
    if not reports_dir.exists():
        print(f"Creating directory: {reports_dir}")
        reports_dir.mkdir(parents=True, exist_ok=True)
        print("✓ Directory created")
    else:
        print(f"✓ Directory exists: {reports_dir}")
    
    print()
    return reports_dir


def main(model_name=None, config_file=None):
    """Main execution pipeline."""
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Implement pathways into a genome-scale metabolic model',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  python3 pathway_implementation.py                           # Use config.json
  python3 pathway_implementation.py -c config_glycocalyx.json # Use specific config file
  python3 pathway_implementation.py --config configs/cytoskeleton.json  # Use config from subfolder

The config file specifies the base model, output path, compartments, metabolites, and reactions.
        '''
    )
    parser.add_argument('-c', '--config', dest='config_file', 
                       help='Path to configuration file (default: config.json)')
    parser.add_argument('model_name', nargs='?', 
                       help='[DEPRECATED] Base name of the model file. Use --config instead.')
    
    args = parser.parse_args()
    
    # Use provided config_file or from command line
    if config_file is None:
        config_file = args.config_file
    
    # Use provided model_name or from command line (legacy support)
    if model_name is None:
        model_name = args.model_name
    
    # Start timing and create timestamp
    start_time = datetime.now()
    timestamp = start_time.strftime("%Y%m%d_%H%M%S")
    
    print("=" * 70)
    print("PATHWAY IMPLEMENTATION - MAIN PIPELINE")
    print("=" * 70)
    print(f"\nExecution started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    base_path = Path(__file__).parent
    
    # Determine model paths
    if model_name is None:
        # Use config file (default: config.json)
        try:
            config_path = config_file if config_file else "config.json"
            print(f"Loading configuration from: {config_path}")
            config = load_config(config_path)
            model_path_str, output_path_str = get_model_paths(config)
            model_path = Path(model_path_str)
            output_path = Path(output_path_str)
            print(f"  ✓ Configuration loaded successfully")
        except Exception as e:
            print(f"❌ ERROR: Could not load configuration file: {e}")
            print("\nPlease either:")
            print("  1. Ensure config.json exists in the project root, OR")
            print("  2. Specify a config file with -c/--config option")
            print("\nExample: python3 pathway_implementation.py --config config_glycocalyx.json")
            sys.exit(1)
    else:
        # Use command-line argument (legacy mode)
        print("WARNING: Using legacy model_name argument is deprecated.")
        print("         Please use --config option instead.")
        model_path = base_path / "model" / f"{model_name}.json"
        output_path = base_path / "model" / f"{model_name}_with_pathways.json"
        print(f"Using command-line model name: {model_name}")
        # Load config for other settings
        config = load_config(config_file if config_file else "config.json")
    
    print(f"Model: {model_path.name}\n")
    
    # Check if model exists
    if not model_path.exists():
        print(f"❌ ERROR: Model file not found: {model_path}")
        print(f"\nPlease ensure your model file exists at: {model_path}")
        sys.exit(1)
    
    print(f"Input model: {model_path}")
    print(f"Output model: {output_path}")
    print()
    
    # Check Python requirements
    if not check_python_requirements():
        sys.exit(1)
    
    print("\nThis script will execute the following steps:")
    print("  1. Analyze model annotations")
    print("  2. Build metabolite ID database")
    print(f"  3. Implement pathways defined in {config_file if config_file else 'config.json'}")
    print()
    
    # User confirmation
    response = input("Proceed with implementation? (yes/no): ").strip().lower()
    if response not in ['yes', 'y']:
        print("Implementation cancelled by user.")
        sys.exit(0)
    
    print("\n" + "=" * 70)
    
    # Initialize report data
    report_data = {
        "timestamp": start_time.strftime('%Y-%m-%d %H:%M:%S'),
        "input_model": str(model_path),
        "output_model": str(output_path),
        "steps": []
    }
    
    # ========================================================================
    # STEP 0: Check and create directories
    # ========================================================================
    print("\nSTEP 0: Checking directories...")
    print("-" * 70)
    reports_dir = check_and_create_directories()
    
    # Get database path from config
    db_path_from_config = config.get('metabolites', {}).get('database_file', 'data/metabolite_id_database.json')
    id_database_path = base_path / db_path_from_config
    
    # Generate date string for report (format: YYYY-MM-DD)
    date_str = datetime.now().strftime('%Y-%m-%d')
    
    # Get pathway name from config for report filename
    pathway_name = config.get('pathway_name', 'pathway')
    
    report_path = reports_dir / f"{date_str}_report_{pathway_name}_implementation.txt"
    
    # ========================================================================
    # STEP 1: Analyze model annotations
    # ========================================================================
    print("\nSTEP 1: Analyzing model annotations...")
    print("-" * 70)
    
    # Skip annotation analysis if database already exists
    if id_database_path.exists():
        print(f"Metabolite ID database already exists at: {id_database_path}")
        print("✓ Skipping annotation analysis (using existing database)")
        annotation_fields = {}
        report_data["steps"].append({
            "step": 1,
            "name": "Analyze model annotations",
            "status": "skipped",
            "reason": "Database already exists"
        })
    else:
        try:
            annotation_fields = analyze_annotations.analyze_model_annotations(str(model_path))
            
            # Store in report
            report_data["steps"].append({
                "step": 1,
                "name": "Analyze model annotations",
                "status": "completed",
                "annotation_fields": dict(annotation_fields)
            })
            
            print("\n✓ Step 1 complete: Annotation analysis finished")
        except Exception as e:
            print(f"❌ ERROR in Step 1: {e}")
            report_data["steps"].append({
                "step": 1,
                "name": "Analyze model annotations",
                "status": "failed",
                "error": str(e)
            })
            sys.exit(1)
    
    # ========================================================================
    # STEP 2: Build metabolite ID database
    # ========================================================================
    print("\n" + "=" * 70)
    print("\nSTEP 2: Building metabolite ID database...")
    print("-" * 70)
    try:
        # Check if database already exists
        if id_database_path.exists():
            print(f"Database already exists at: {id_database_path}")
            print("Loading existing database...")
            with open(id_database_path, 'r') as f:
                id_db = json.load(f)
            print(f"✓ Loaded database with {len(id_db)} metabolites")
        else:
            print("Database not found. Building new database...")
            # Use integrated pipeline: analyze_annotations -> build_id_database
            id_db = build_id_database.main()
            
            # Verify database was created
            if id_db is None:
                raise Exception("Database creation returned None")
        
        report_data["steps"].append({
            "step": 2,
            "name": "Build metabolite ID database",
            "status": "completed",
            "database_path": str(id_database_path),
            "metabolite_count": len(id_db)
        })
        
        print("✓ Step 2 complete: Metabolite ID database ready")
    except Exception as e:
        print(f"❌ ERROR in Step 2: {e}")
        report_data["steps"].append({
            "step": 2,
            "name": "Build metabolite ID database",
            "status": "failed",
            "error": str(e)
        })
        sys.exit(1)
    
    # ========================================================================
    # STEP 3: Implement pathways from config
    # ========================================================================
    print("\n" + "=" * 70)
    print("\nSTEP 3: Implementing pathways...")
    print("-" * 70)
    
    original_met_count = 0
    original_rxn_count = 0
    met_count = 0
    rxn_count = 0
    
    try:
        # Load model
        print(f"Loading model from {model_path}...")
        with open(model_path, 'r') as f:
            model = json.load(f)
        original_met_count = len(model['metabolites'])
        original_rxn_count = len(model['reactions'])
        print(f"✓ Model loaded: {original_met_count} metabolites, {original_rxn_count} reactions")
        
        # Check for existing pathway components
        print("\nChecking for existing pathway components...")
        
        # Determine pathway identifier from config
        pathway_name = config.get('pathway_name', '')
        # Get first compartment abbreviation as a proxy for pathway-specific compartment
        first_comp_abbrev = config.get('compartments', [{}])[0].get('abbreviation', 'ck') if config.get('compartments') else 'ck'
        
        # Skip the check - we'll handle duplicates during creation
        print("✓ Skipping duplicate check (will be handled during metabolite/reaction creation)")
        
        # Load ID database
        print(f"\nLoading metabolite ID database from {id_database_path}...")
        with open(id_database_path, 'r') as f:
            id_database = json.load(f)
        print(f"✓ Loaded {len(id_database)} metabolites from ID database")
        
        # Resolve compartments from config
        from functions.config import resolve_all_compartments
        print("\nResolving compartments from config...")
        resolved_compartments = resolve_all_compartments(model, config)
        
        # Add compartments and build abbreviation mapping
        compartment_map = {}  # Maps config name to actual abbreviation used
        abbrev_map = {}  # Maps config abbreviation to model abbreviation
        
        for comp_config in config.get('compartments', []):
            config_abbrev = comp_config['abbreviation']
            config_name = comp_config['name']
            
            # Find corresponding resolved compartment
            resolved = next((c for c in resolved_compartments if c['name'] == config_name), None)
            if resolved:
                model_abbrev = resolved['abbreviation']
                abbrev_map[config_abbrev] = model_abbrev
                compartment_map[config_name] = model_abbrev
                
                if resolved['exists']:
                    print(f"  ✓ Using existing compartment '{config_name}' [{config_abbrev}→{model_abbrev}]")
                else:
                    print(f"  + Creating compartment '{config_name}' as [{model_abbrev}]")
                    model, _, _ = pathway_builder.add_compartment(
                        model, 
                        model_abbrev, 
                        config_name
                    )
        
        # Process compartments in two phases:
        # Phase 1: Create all metabolites across all compartments
        # Phase 2: Create all reactions (which may reference metabolites from multiple compartments)
        total_met_count = 0
        total_rxn_count = 0
        
        print("\n" + "="*70)
        print("PHASE 1: Creating metabolites in all compartments")
        print("="*70)
        
        for comp in resolved_compartments:
            comp_name = comp['name']
            comp_abbrev = comp['abbreviation']
            
            # Determine which metabolite key to use for this compartment
            # Try pathway-specific keys first
            met_key = None
            pathway_name = config.get('pathway_name', '').replace(' ', '_').replace('-', '_')
            
            # Try several possible keys in order of specificity
            possible_keys = [
                f"{pathway_name}_specific",  # e.g., "Syndecan1_HS_specific"
                'pathway_specific',  # Generic pathway-specific metabolites
                'cytoskeleton_specific',  # Legacy cytoskeleton
                'glycocalyx_specific',  # Legacy glycocalyx
                f"{comp_name}_specific",  # Compartment-specific
            ]
            
            # First check if any of these keys exist at top level
            for key in possible_keys:
                if key in config:
                    met_key = key
                    break
            
            # If not found at top level, check inside 'metabolites' dict
            if met_key is None and 'metabolites' in config and isinstance(config['metabolites'], dict):
                for key in possible_keys:
                    if key in config['metabolites']:
                        met_key = ('metabolites', key)  # Tuple indicating nested location
                        break
            
            # Fall back to generic "metabolites" if nothing else works
            if met_key is None and 'metabolites' in config:
                met_key = 'metabolites'
            
            # Create metabolites for this compartment
            print(f"\nCreating metabolites for [{comp_abbrev}] {comp_name}...")
            met_count = pathway_builder.create_compartment_metabolites(
                model, id_database, comp_abbrev, config, met_key, abbrev_map
            )
            print(f"✓ Created {met_count} metabolites in [{comp_abbrev}] compartment")
            total_met_count += met_count
        
        print("\n" + "="*70)
        print("PHASE 2: Creating reactions across all compartments")
        print("="*70)
        
        for comp in resolved_compartments:
            comp_name = comp['name']
            comp_abbrev = comp['abbreviation']
            
            # Create reactions for this compartment
            print(f"\nCreating reactions for [{comp_abbrev}] {comp_name}...")
            rxn_count = pathway_builder.create_compartment_reactions(
                model, id_database, comp_abbrev, config, abbrev_map
            )
            print(f"✓ Created {rxn_count} reactions involving [{comp_abbrev}]")
            total_rxn_count += rxn_count
        
        met_count = total_met_count
        rxn_count = total_rxn_count
        
        # Save modified model
        print(f"\nSaving modified model to {output_path}...")
        with open(output_path, 'w') as f:
            json.dump(model, f, indent=2)
        print("✓ Model saved successfully")
        
        # Store in report
        report_data["steps"].append({
            "step": 3,
            "name": "Implement pathways",
            "status": "completed",
            "original_metabolites": original_met_count,
            "original_reactions": original_rxn_count,
            "final_metabolites": len(model['metabolites']),
            "final_reactions": len(model['reactions']),
            "new_metabolites": met_count,
            "new_reactions": rxn_count
        })
        
        print("\n✓ Step 3 complete: Pathway implementation finished")
        
    except Exception as e:
        print(f"❌ ERROR in Step 3: {e}")
        import traceback
        traceback.print_exc()
        report_data["steps"].append({
            "step": 3,
            "name": "Implement pathways",
            "status": "failed",
            "error": str(e)
        })
        sys.exit(1)
    
    # ========================================================================
    # GENERATE REPORT
    # ========================================================================
    end_time = datetime.now()
    duration = end_time - start_time
    
    report_data["completion_time"] = end_time.strftime('%Y-%m-%d %H:%M:%S')
    report_data["duration_seconds"] = duration.total_seconds()
    report_data["summary"] = {
        "compartments_added": f"Pathways defined in {config_file if config_file else 'config.json'}",
        "original_metabolites": original_met_count,
        "original_reactions": original_rxn_count,
        "final_metabolites": len(model['metabolites']),
        "final_reactions": len(model['reactions']),
        "new_metabolites": met_count,
        "new_reactions": rxn_count
    }
    
    # Write report to file
    with open(report_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("PATHWAY IMPLEMENTATION REPORT\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Execution Date: {report_data['timestamp']}\n")
        f.write(f"Completion Date: {report_data['completion_time']}\n")
        f.write(f"Duration: {duration.total_seconds():.2f} seconds\n\n")
        
        f.write("-" * 70 + "\n")
        f.write("INPUT/OUTPUT FILES\n")
        f.write("-" * 70 + "\n")
        f.write(f"Input model:  {report_data['input_model']}\n")
        f.write(f"Output model: {report_data['output_model']}\n")
        f.write(f"ID database:  {id_database_path}\n\n")
        
        f.write("-" * 70 + "\n")
        f.write("EXECUTION STEPS\n")
        f.write("-" * 70 + "\n\n")
        
        for step_info in report_data['steps']:
            f.write(f"Step {step_info['step']}: {step_info['name']}\n")
            f.write(f"  Status: {step_info['status']}\n")
            
            if step_info['step'] == 1 and 'annotation_fields' in step_info:
                f.write(f"  Annotation fields found: {len(step_info['annotation_fields'])}\n")
                f.write("  Top annotation databases:\n")
                sorted_fields = sorted(step_info['annotation_fields'].items(), 
                                     key=lambda x: x[1], reverse=True)
                for field, count in sorted_fields[:5]:
                    f.write(f"    - {field}: {count} metabolites\n")
            
            elif step_info['step'] == 2 and 'metabolite_count' in step_info:
                f.write(f"  Metabolites in ID database: {step_info['metabolite_count']}\n")
            
            elif step_info['step'] == 3 and 'new_metabolites' in step_info:
                f.write(f"  Original metabolites: {step_info['original_metabolites']}\n")
                f.write(f"  Original reactions:   {step_info['original_reactions']}\n")
                f.write(f"  Final metabolites:    {step_info['final_metabolites']}\n")
                f.write(f"  Final reactions:      {step_info['final_reactions']}\n")
                f.write(f"  New metabolites:      {step_info['new_metabolites']}\n")
                f.write(f"  New reactions:        {step_info['new_reactions']}\n")
            
            f.write("\n")
        
        f.write("-" * 70 + "\n")
        f.write("SUMMARY\n")
        f.write("-" * 70 + "\n")
        f.write(f"Pathways added: {report_data['summary']['compartments_added']}\n")
        f.write(f"Total metabolites: {report_data['summary']['final_metabolites']} "
                f"(+{report_data['summary']['new_metabolites']})\n")
        f.write(f"Total reactions:   {report_data['summary']['final_reactions']} "
                f"(+{report_data['summary']['new_reactions']})\n")
        f.write("\n" + "=" * 70 + "\n")
        f.write("Implementation completed successfully!\n")
        f.write("=" * 70 + "\n")
    
    # ========================================================================
    # FINAL SUMMARY TO CONSOLE
    # ========================================================================
    print("\n" + "=" * 70)
    print("IMPLEMENTATION COMPLETE!")
    print("=" * 70)
    print(f"\nModified model saved to:")
    print(f"  {output_path}")
    print(f"\nMetabolite ID database saved to:")
    print(f"  {id_database_path}")
    print(f"\nReport saved to:")
    print(f"  {report_path}")
    
    # Build compartment summary string
    comp_summary = ', '.join([f"[{c['abbreviation']}] {c['name']}" for c in resolved_compartments])
    
    print("\nSummary:")
    print(f"  - New compartments: {comp_summary}")
    print(f"  - Total metabolites: {len(model['metabolites'])} (+{met_count})")
    print(f"  - Total reactions: {len(model['reactions'])} (+{rxn_count})")
    print(f"  - Execution time: {duration.total_seconds():.2f} seconds")
    
    print("\n" + "=" * 70)
    print("All steps completed successfully! ✓")
    print("=" * 70)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nProcess interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\n\n❌ FATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
