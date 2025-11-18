#!/usr/bin/env python3
"""
Helper script to add a new reaction to config.json

Usage:
    python3 add_reaction.py

The script will prompt for all necessary information and add the reaction to config.json
"""

import json
import sys

def add_reaction_interactive():
    """Interactive helper to add a new reaction."""
    
    print("=" * 70)
    print("ADD NEW REACTION TO CONFIG.JSON")
    print("=" * 70)
    print("\nThis helper will guide you through adding a new reaction.")
    print("Press Ctrl+C at any time to cancel.\n")
    
    try:
        # Basic information
        rxn_id = input("Reaction ID (e.g., R_MyReaction): ").strip()
        if not rxn_id:
            print("✗ Reaction ID is required")
            return False
        
        name = input("Reaction name: ").strip()
        if not name:
            print("✗ Reaction name is required")
            return False
        
        description = input("Description (optional): ").strip()
        
        # Equation
        print("\nEquation format examples:")
        print("  - Forward: 1 A[c] + 2 B[c] --> 1 C[c]")
        print("  - Reversible: 1 A[c] + 2 B[c] <=> 1 C[c]")
        print("  - Sink: 1 A[c] -->")
        print("  Note: Use metabolite names as defined in config.json")
        print("        Include compartment in brackets [c], [ck], etc.")
        
        equation = input("\nEquation: ").strip()
        if not equation:
            print("✗ Equation is required")
            return False
        
        # Check if reversible
        is_reversible = '<=>' in equation
        
        # Optional fields
        subsystem = input("\nSubsystem (optional): ").strip()
        ec = input("EC number (optional): ").strip()
        gpr = input("Gene-protein-reaction rule (optional): ").strip()
        sbo = input("SBO term (default: SBO:0000176 for biochemical reaction): ").strip()
        if not sbo:
            sbo = "SBO:0000176"
        
        # Bounds
        print("\nReaction bounds:")
        if is_reversible:
            lb_default = "-1000.0"
            print("  (Reversible reaction detected, default: -1000 to 1000)")
        elif '-->' in equation and equation.strip().endswith('-->'):
            lb_default = "0.0"
            print("  (Sink reaction detected, default: 0 to 1000)")
        else:
            lb_default = "0.0"
            print("  (Forward reaction, default: 0 to 1000)")
        
        lb_input = input(f"Lower bound (default: {lb_default}): ").strip()
        lb = float(lb_input) if lb_input else float(lb_default)
        
        ub_input = input("Upper bound (default: 1000.0): ").strip()
        ub = float(ub_input) if ub_input else 1000.0
        
        # Build reaction dict
        new_reaction = {
            "id": rxn_id,
            "name": name,
            "description": description,
            "equation": equation,
            "subsystem": subsystem,
            "ec": ec,
            "gpr": gpr,
            "lower_bound": lb,
            "upper_bound": ub,
            "sbo": sbo
        }
        
        # Show summary
        print("\n" + "=" * 70)
        print("REACTION SUMMARY")
        print("=" * 70)
        for key, value in new_reaction.items():
            if value or key in ['lower_bound', 'upper_bound']:
                print(f"  {key}: {value}")
        
        # Confirm
        confirm = input("\nAdd this reaction to config.json? (yes/no): ").strip().lower()
        if confirm not in ['yes', 'y']:
            print("✗ Cancelled")
            return False
        
        # Load config
        with open('config.json', 'r') as f:
            config = json.load(f)
        
        # Check if reactions exist
        if 'reactions' not in config:
            config['reactions'] = []
        
        # Check for duplicate ID
        if any(r['id'] == rxn_id for r in config['reactions']):
            print(f"\n✗ ERROR: Reaction with ID '{rxn_id}' already exists")
            overwrite = input("Overwrite existing reaction? (yes/no): ").strip().lower()
            if overwrite not in ['yes', 'y']:
                return False
            # Remove old reaction
            config['reactions'] = [r for r in config['reactions'] if r['id'] != rxn_id]
        
        # Add reaction
        config['reactions'].append(new_reaction)
        
        # Save
        with open('config.json', 'w') as f:
            json.dump(config, f, indent=2)
        
        print(f"\n✓ Reaction '{rxn_id}' added successfully!")
        print(f"✓ Total reactions in config: {len(config['reactions'])}")
        
        return True
        
    except KeyboardInterrupt:
        print("\n\n✗ Cancelled by user")
        return False
    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        return False

def add_reaction_from_args(args):
    """Add reaction from command-line arguments."""
    # TODO: Implement command-line argument parsing
    print("Command-line mode not yet implemented. Use interactive mode.")
    return False

if __name__ == '__main__':
    if len(sys.argv) > 1:
        success = add_reaction_from_args(sys.argv[1:])
    else:
        success = add_reaction_interactive()
    
    sys.exit(0 if success else 1)
