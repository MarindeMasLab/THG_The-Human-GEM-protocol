#!/usr/bin/env python3
"""Test script to verify compartment names are set correctly in COBRA models."""

import cobra
import tempfile

# Create a simple model
model = cobra.Model("test_model")

# Add metabolites with different compartments
met1 = cobra.Metabolite("M_glc_c", name="Glucose", compartment="c")
met2 = cobra.Metabolite("M_atp_m", name="ATP", compartment="m")
met3 = cobra.Metabolite("M_nadh_n", name="NADH", compartment="n")

model.add_metabolites([met1, met2, met3])

print("Before setting compartment names:")
print(f"  model.compartments = {model.compartments}")

# Now set compartment names (this is what our fix does)
location_dict = {
    "cytosol": "c",
    "mitochondria": "m",
    "nucleus": "n"
}

# Reverse the dict: compartment_id -> compartment_name
compartment_id_to_name = {v: k for k, v in location_dict.items()}

# Set compartment names
model.compartments = {comp_id: compartment_id_to_name.get(comp_id, comp_id) 
                      for comp_id in model.compartments.keys()}

print("\nAfter setting compartment names:")
print(f"  model.compartments = {model.compartments}")

# Write to SBML and check
with tempfile.NamedTemporaryFile(mode='w', suffix='.xml', delete=False) as f:
    output_file = f.name
    cobra.io.write_sbml_model(model, output_file)
    print(f"\nWrote model to: {output_file}")
    
# Read the SBML and check compartments
with open(output_file, 'r') as f:
    content = f.read()
    if '<listOfCompartments>' in content:
        start = content.index('<listOfCompartments>')
        end = content.index('</listOfCompartments>') + len('</listOfCompartments>')
        compartments_section = content[start:end]
        print("\nCompartments in SBML:")
        print(compartments_section)
        
        # Check if names are present
        if 'name="cytosol"' in compartments_section:
            print("\n✅ SUCCESS: Compartment names are correctly set!")
        else:
            print("\n❌ FAILED: Compartment names are missing!")
