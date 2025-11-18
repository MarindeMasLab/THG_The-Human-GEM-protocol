#!/usr/bin/env python3
"""
Example: Complete workflow with isolated inputs/outputs

This example demonstrates the full pathway implementation workflow with:
- Input files in examples/inputs/
- Output files in examples/outputs/ (data, reports, figures, models)

Run from the implement_pathway directory:
    python3 examples/run_example.py
"""

import sys
from pathlib import Path

# Add parent directory to path for functions module
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from functions import pathway_builder, analyze_annotations, build_id_database
from functions.config import load_config
import json
import os
from datetime import datetime

def main():
    """Run complete example workflow."""
    
    print("=" * 80)
    print("EXAMPLE: COMPLETE PATHWAY IMPLEMENTATION WORKFLOW")
    print("=" * 80)
    print()
    
    # Define paths
    example_dir = Path(__file__).parent
    inputs_dir = example_dir / "inputs"
    outputs_dir = example_dir / "outputs"
    
    # Input files
    config_path = inputs_dir / "config_example.json"
    model_path = inputs_dir / "endoA_250917_3.json"
    
    # Output directories
    data_dir = outputs_dir / "data"
    reports_dir = outputs_dir / "reports"
    figures_dir = outputs_dir / "figures"
    
    # Output files
    output_model_path = outputs_dir / "endoA_250917_3_with_Example_Pathway.json"
    database_path = data_dir / "metabolite_id_database.json"
    
    # Ensure output directories exist
    for dir_path in [data_dir, reports_dir, figures_dir]:
        dir_path.mkdir(parents=True, exist_ok=True)
    
    print("📁 Directory structure:")
    print(f"  Inputs:  {inputs_dir}")
    print(f"  Outputs: {outputs_dir}")
    print()
    
    # ========================================================================
    # STEP 1: Load configuration
    # ========================================================================
    print("STEP 1: Loading configuration")
    print("-" * 70)
    
    with open(config_path, 'r') as f:
        config = json.load(f)
    
    pathway_name = config.get('pathway_name', 'Example_Pathway')
    print(f"✓ Loaded config for: {pathway_name}")
    print(f"  Config file: {config_path}")
    print()
    
    # ========================================================================
    # STEP 2: Build metabolite ID database
    # ========================================================================
    print("STEP 2: Building metabolite ID database")
    print("-" * 70)
    
    if database_path.exists():
        print(f"✓ Database already exists: {database_path}")
        with open(database_path, 'r') as f:
            id_database = json.load(f)
    else:
        # Extract metabolite annotations from model
        targets = config['metabolites']['targets']
        print(f"  Analyzing {len(targets)} target metabolites...")
        
        metabolite_annotations, not_found = analyze_annotations.extract_metabolite_annotations(
            str(model_path), targets, verbose=False
        )
        
        print(f"  ✓ Found annotations for {len(metabolite_annotations)} metabolites")
        if not_found:
            print(f"  ⚠ {len(not_found)} metabolites not found: {', '.join(not_found)}")
        
        # Build database
        id_database = build_id_database.build_database_from_annotations(metabolite_annotations)
        
        # Save database
        with open(database_path, 'w') as f:
            json.dump(id_database, f, indent=2)
        
        print(f"  ✓ Database saved: {database_path}")
    
    print(f"  Database contains {len(id_database)} metabolites")
    print()
    
    # ========================================================================
    # STEP 3: Load base model
    # ========================================================================
    print("STEP 3: Loading base model")
    print("-" * 70)
    
    with open(model_path, 'r') as f:
        model = json.load(f)
    
    initial_mets = len(model['metabolites'])
    initial_rxns = len(model['reactions'])
    
    print(f"  ✓ Base model loaded")
    print(f"    Metabolites: {initial_mets}")
    print(f"    Reactions: {initial_rxns}")
    print()
    
    # ========================================================================
    # STEP 4: Add compartment
    # ========================================================================
    print("STEP 4: Adding compartment")
    print("-" * 70)
    
    compartments = config['compartments']
    comp_abbrev = compartments[0]['abbreviation']
    comp_name = compartments[0]['name']
    
    model, abbrev_used, already_exists = pathway_builder.add_compartment(
        model, comp_abbrev, comp_name
    )
    
    if already_exists:
        print(f"  ✓ Compartment '{comp_name}' [{abbrev_used}] already exists")
    else:
        print(f"  ✓ Added compartment '{comp_name}' [{abbrev_used}]")
    print()
    
    # ========================================================================
    # STEP 5: Create metabolites
    # ========================================================================
    print("STEP 5: Creating metabolites")
    print("-" * 70)
    
    met_count = pathway_builder.create_compartment_metabolites(
        model, id_database, abbrev_used, config, 'cytoskeleton_specific'
    )
    
    print(f"  ✓ Created {met_count} metabolites in {comp_name}")
    print()
    
    # ========================================================================
    # STEP 6: Create reactions
    # ========================================================================
    print("STEP 6: Creating reactions")
    print("-" * 70)
    
    rxn_count = pathway_builder.create_compartment_reactions(
        model, id_database, abbrev_used, config
    )
    
    print(f"  ✓ Created {rxn_count} reactions")
    print()
    
    # ========================================================================
    # STEP 7: Save modified model
    # ========================================================================
    print("STEP 7: Saving modified model")
    print("-" * 70)
    
    with open(output_model_path, 'w') as f:
        json.dump(model, f, indent=2)
    
    final_mets = len(model['metabolites'])
    final_rxns = len(model['reactions'])
    
    print(f"  ✓ Model saved: {output_model_path}")
    print(f"    Metabolites: {initial_mets} → {final_mets} (+{final_mets - initial_mets})")
    print(f"    Reactions: {initial_rxns} → {final_rxns} (+{final_rxns - initial_rxns})")
    print()
    
    # ========================================================================
    # STEP 8: Generate report
    # ========================================================================
    print("STEP 8: Generating report")
    print("-" * 70)
    
    date_str = datetime.now().strftime('%Y-%m-%d')
    report_path = reports_dir / f"{date_str}_example_implementation_report.txt"
    
    with open(report_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("EXAMPLE PATHWAY IMPLEMENTATION REPORT\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Pathway: {pathway_name}\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("INPUT FILES\n")
        f.write("-" * 80 + "\n")
        f.write(f"Config: {config_path}\n")
        f.write(f"Model:  {model_path}\n\n")
        f.write("OUTPUT FILES\n")
        f.write("-" * 80 + "\n")
        f.write(f"Database: {database_path}\n")
        f.write(f"Model:    {output_model_path}\n")
        f.write(f"Report:   {report_path}\n\n")
        f.write("RESULTS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Compartment added: {comp_name} [{abbrev_used}]\n")
        f.write(f"Metabolites added: {met_count}\n")
        f.write(f"Reactions added:   {rxn_count}\n\n")
        f.write(f"Total metabolites: {initial_mets} → {final_mets}\n")
        f.write(f"Total reactions:   {initial_rxns} → {final_rxns}\n\n")
        f.write("=" * 80 + "\n")
    
    print(f"  ✓ Report saved: {report_path}")
    print()
    
    # ========================================================================
    # SUMMARY
    # ========================================================================
    print("=" * 80)
    print("✅ EXAMPLE WORKFLOW COMPLETE")
    print("=" * 80)
    print()
    print("Generated files:")
    print(f"  📊 Database:  {database_path}")
    print(f"  🧬 Model:     {output_model_path}")
    print(f"  📄 Report:    {report_path}")
    print()
    print("Next steps:")
    print("  1. Validate: python3 validate.py 0 --config examples/inputs/config_example.json --model examples/outputs/endoA_250917_3_with_Example_Pathway.json")
    print("  2. Visualize: python3 visualize.py 0 --config examples/inputs/config_example.json --model examples/outputs/endoA_250917_3_with_Example_Pathway.json")
    print()
    print(f"All outputs saved to: {outputs_dir}")
    print()

if __name__ == '__main__':
    main()
