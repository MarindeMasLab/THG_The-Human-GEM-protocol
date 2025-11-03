#!/usr/bin/env python3
"""
Two-pass metabolite annotation: 
1. First pass annotates all metabolites
2. Second pass retries failures (which may be temporary API errors)
"""

import cobra.io
import os
import sys

# Setup paths
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.join(current_dir, "..")

if project_root not in sys.path:
    sys.path.append(project_root)

from functions.function_metabolite_identification import *
from functions.function_reac_identification import *
from functions.function_annotate_cobra_model import *

# Load model
model = os.path.join(project_root, "models", "Human-GEM_2022-06-21.xml")
try:
    cobra_model = cobra.io.read_sbml_model(model)
except:
    # Try alternative model name
    model = os.path.join(project_root, "models", "Human-GEM1_19.xml")
    cobra_model = cobra.io.read_sbml_model(model)

print("\n" + "="*70)
print("METABOLITE ANNOTATION - INTELLIGENT MULTI-PASS SYSTEM")
print("="*70 + "\n")

# Configuration
MAX_RETRY_ROUNDS = 10  # Maximum number of retry rounds
STOP_THRESHOLD = 0.10  # Stop if remaining failures < 10% of potential achievable

print(f"Configuration:")
print(f"  Max retry rounds: {MAX_RETRY_ROUNDS}")
print(f"  Stop threshold: {STOP_THRESHOLD*100:.0f}% of potential achievable")
print(f"{'='*70}\n")

# Gather metabolites
metabolites = gather_metabolites(cobra_model)
total_metabolites = len(metabolites)

print(f"Total metabolites to process: {total_metabolites}\n")

# ===== FIRST PASS =====
print("="*70)
print("STARTING FIRST PASS")
print("="*70 + "\n")

annotated, unannotated = generate_met_annotation(metabolites)

print(f"\nFirst pass complete!")
print(f"Annotated: {len(annotated)}, Unannotated: {len(unannotated)}")

# ===== INTELLIGENT MULTI-PASS RETRY SYSTEM =====
# Use configuration from above

if len(unannotated) > 0:
    total_annotated = len(annotated)
    to_retry = unannotated
    retry_round = 1
    
    while retry_round <= MAX_RETRY_ROUNDS and len(to_retry) > 0:
        print("\n" + "="*70)
        print(f"STARTING RETRY ROUND {retry_round}")
        print("="*70 + "\n")
        print(f"Attempting to annotate {len(to_retry)} failed metabolites...")
        
        annotation_file = global_met_annotation_file()
        retry_success = 0
        retry_fail = []
        
        with open(annotation_file, 'a') as f:
            for idx, (name, formula, annotation, iden) in enumerate(to_retry, 1):
                if idx % 50 == 0:
                    print(f"  Progress: {idx}/{len(to_retry)} - New successes this round: {retry_success}")
                
                try:
                    result = identify_metabolite(name, formula, iden)
                    if result is not None:
                        f.write(result)
                        f.flush()
                        retry_success += 1
                        total_annotated += 1
                    else:
                        retry_fail.append((name, formula, annotation, iden))
                except Exception as e:
                    retry_fail.append((name, formula, annotation, iden))
                    continue
        
        print(f"\n{'='*70}")
        print(f"RETRY ROUND {retry_round} COMPLETE")
        print(f"{'='*70}")
        print(f"New successes this round: {retry_success}")
        print(f"Still failed: {len(retry_fail)}")
        print(f"Total annotated so far: {total_annotated} / {total_metabolites} ({(total_annotated/total_metabolites)*100:.1f}%)")
        
        # Calculate if we should continue
        if retry_success == 0:
            print(f"\nNo new successes in this round. Stopping retries.")
            print(f"Remaining {len(retry_fail)} failures are likely not in PubChem database.")
            break
        
        # Calculate potential maximum (total - those likely not in DB)
        # Estimate: if we got X successes this round, there might be more API errors
        potential_max = total_annotated + len(retry_fail)
        remaining_ratio = len(retry_fail) / potential_max if potential_max > 0 else 0
        
        print(f"\nEstimated maximum achievable: {potential_max} metabolites")
        print(f"Remaining failures: {len(retry_fail)} ({remaining_ratio*100:.1f}% of potential max)")
        
        if remaining_ratio <= STOP_THRESHOLD:
            print(f"\n✓ Remaining failures ({remaining_ratio*100:.1f}%) are below threshold ({STOP_THRESHOLD*100:.0f}%)")
            print(f"  These are likely metabolites not in PubChem database.")
            print(f"  Stopping retry process.")
            break
        
        # Prepare for next round
        to_retry = retry_fail
        retry_round += 1
        
        if retry_round <= MAX_RETRY_ROUNDS and len(to_retry) > 0:
            print(f"\n→ Will attempt round {retry_round} with {len(to_retry)} remaining failures...")
    
    if retry_round > MAX_RETRY_ROUNDS:
        print(f"\n⚠ Reached maximum retry limit ({MAX_RETRY_ROUNDS} rounds)")
    
    print(f"\n{'='*70}")
    print("ALL RETRY ROUNDS COMPLETE")
    print(f"{'='*70}")
    print(f"Total retry rounds executed: {retry_round - 1}")
    print(f"Final annotated count: {total_annotated} / {total_metabolites}")
    print(f"Final success rate: {(total_annotated/total_metabolites)*100:.1f}%")
    print(f"Remaining unannotated: {len(retry_fail) if retry_round > 1 else len(unannotated)}")
    print(f"{'='*70}\n")
else:
    print("\nNo failures to retry!")

print("Annotation process complete!")
print("Continuing with reaction identification...\n")

# Continue with rest of the pipeline
met_annotation = process_annotation()

database = os.path.join(project_root, "models", "Human Database.xml")

# Reactions
reac = process_reac(model, '[A-Z]+[0-9]+[a-z]+[0-9]*', 'MAM02040', 'MAM02039')
reac = replace_met_id_by_met_kegg(reac, gather_kegg_metabolites(model))
reac_y = process_reac(database, '([A-Z][0-9]+_?[a-z]+[0-9]*)', 'C00080', 'C00001')
jaccard = execute_jaccard(reac_y, reac)
reac_annotation = process_jaccard(reac_y, reac, jaccard)

# SBML
annotate_cobra_model(cobra_model, met_annotation, reac_annotation)

print("\nAll done!")
