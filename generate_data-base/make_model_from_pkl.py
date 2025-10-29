#!/usr/bin/python
import sys
import os
import dill
import cobra

# Add the parent directory to the Python path so 'functions' module can be imported
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.join(current_dir, "..")

if project_root not in sys.path:
    sys.path.append(project_root)

# Import the cobra_reconstruction function from generate_db
from generate_db import cobra_reconstruction

# Try both pickle files
pickle_files = [
    os.path.join(project_root, "files", "pre_sbml_pos_comp.pk"),
    os.path.join(project_root, "files", "pre_sbml_raw.pk"),
]

for pickle_file in pickle_files:
    print(f"\n{'='*80}")
    print(f"Checking: {os.path.basename(pickle_file)}")
    print("=" * 80)

    if not os.path.exists(pickle_file):
        print(f"File not found: {pickle_file}")
        continue

    with open(pickle_file, "rb") as f:
        data = dill.load(f)

    print(f"\nData type: {type(data)}")
    print(f"Keys: {list(data.keys())}")

    # Check which has actual data
    print(f"\nData counts:")
    for key in ["mets", "mets_cl", "reactions", "reactions_cl", "genes"]:
        if key in data:
            count = len(data[key]) if isinstance(data[key], dict) else "N/A"
            print(f"  {key}: {count}")

    # Use the file that has compartmentalized data
    if data.get("reactions_cl") and len(data["reactions_cl"]) > 0:
        print(f"\n✓ This file has compartmentalized data!")
        selected_pickle = pickle_file
        break
    elif data.get("reactions") and len(data["reactions"]) > 0:
        print(f"\n⚠ This file has reactions but they're not compartmentalized yet")
        # We'll use this if it's the best we have
        selected_pickle = pickle_file

print(f"\n{'='*80}")
print(f"Using pickle file: {os.path.basename(selected_pickle)}")
print("=" * 80)

# Load the selected pickle data
with open(selected_pickle, "rb") as f:
    data = dill.load(f)

print("\n" + "=" * 80)
print("PICKLE DATA STRUCTURE:")
print("=" * 80)
print(f"\nData type: {type(data)}")
print(f"\nKeys in data dictionary: {list(data.keys())}")
print("\nData structure details:")
for key, value in data.items():
    print(f"  - {key}: {type(value)}")
    if isinstance(value, dict):
        print(f"    Number of items: {len(value)}")
        if value:
            first_key = list(value.keys())[0]
            print(f"    Sample key: {first_key}")
            print(f"    Sample value type: {type(value[first_key])}")
    elif isinstance(value, str):
        print(
            f"    Value: {value[:100]}..."
            if len(value) > 100
            else f"    Value: {value}"
        )

print("\n" + "=" * 80)
print("GENERATING MODEL FROM PICKLE DATA...")
print("=" * 80)

# Extract data from the pickle
ModName = data.get("name", "Human_Database")
ModID = data.get("id", "Human_Database")
MetList_CL = data.get("mets_cl", {})
RxnList_CL = data.get("reactions_cl", {})
GeneList = data.get("genes", {})
PathNameRxn = data.get("pathways", {})
LocVar = data.get("loc", {})
MetEquiv = data.get("met_equiv", {})
MetList = data.get("mets", {})

# Check if we have compartmentalized data, if not, we need to use the raw data
if len(RxnList_CL) == 0 and len(data.get("reactions", {})) > 0:
    print("\n⚠ WARNING: Compartmentalized reactions are empty!")
    print("The pickle file appears to be incomplete.")
    print(
        "This suggests the generate_db.py script didn't finish the compartmentalization step."
    )
    print("\nYou need to either:")
    print("1. Re-run generate_db.py to completion, or")
    print("2. Check the log file: logs/generate_db.log for errors")

    print(f"\nRaw data available:")
    print(f"  - reactions: {len(data.get('reactions', {}))}")
    print(f"  - mets: {len(data.get('mets', {}))}")
    print(f"  - genes: {len(data.get('genes', {}))}")

    import sys

    sys.exit(1)

print(f"\nModel Name: {ModName}")
print(f"Model ID: {ModID}")
print(f"Number of metabolites (compartmentalized): {len(MetList_CL)}")
print(f"Number of reactions (compartmentalized): {len(RxnList_CL)}")
print(f"Number of genes: {len(GeneList)}")
print(f"Number of pathways: {len(PathNameRxn)}")
print(f"Number of compartments: {len(LocVar)}")

# Generate the COBRA model
print("\nBuilding COBRA model...")
model = cobra_reconstruction(
    ModName,
    ModID,
    MetList_CL,
    RxnList_CL,
    GeneList,
    PathNameRxn,
    LocVar,
    MetEquiv,
    MetList,
)

# Define output path
Output = os.path.join(project_root, "models", "Human_Database_from_pickle.xml")
print(f"\nWriting model to: {Output}")
cobra.io.write_sbml_model(model, Output)

print("\n" + "=" * 80)
print("MODEL GENERATION COMPLETE!")
print("=" * 80)
print(f"\nModel summary:")
print(f"  Metabolites: {len(model.metabolites)}")
print(f"  Reactions: {len(model.reactions)}")
print(f"  Genes: {len(model.genes)}")
print(f"  Groups (pathways): {len(model.groups)}")
print(f"\nModel saved to: {Output}")
print("=" * 80)
