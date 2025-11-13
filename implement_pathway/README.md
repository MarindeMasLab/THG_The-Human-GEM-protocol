# Pathway Implementation for Genome-Scale Metabolic Models

A comprehensive, **pathway-agnostic pipeline** for adding any metabolic pathway (cytoskeleton, glycocalyx, ECM, etc.) to genome-scale metabolic models (GEMs). This implementation uses a fully configurable approach with metabolite annotation-based matching and universal reaction format.

## Repository Structure

This module (`implement_pathway/`) is part of a larger repository structure:

```
/
├── models/          # Shared metabolic models (parent directory)
├── functions/       # Shared core implementation modules (parent directory)
└── implement_pathway/  # This module - pathway implementation tools
```

**Key Design**: The `models/` and `functions/` directories are at the parent level to allow multiple pathway modules to share resources. Scripts in `implement_pathway/` import from the parent `functions/` module and reference models in the parent `models/` directory.

## 📋 Table of Contents

- [Repository Structure](#repository-structure)
- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Validation & Visualization](#validation--visualization)
- [Examples](#examples)
- [Extending to Other Networks](#extending-to-other-networks)
- [Troubleshooting](#troubleshooting)

---

## 🔬 Overview

This project provides a **generic, pathway-agnostic framework** for implementing metabolic pathways into genome-scale models. Originally developed for cytoskeletal pathways (actin, tubulin, vimentin, myosin), it now supports **any pathway/compartment** through configuration files.

### Pathway Features
- **Protein synthesis reactions** (amino acid polymerization)
- **Post-translational modifications** (glycosylation, phosphorylation)
- **Energy dissipation** (ATP/GTP hydrolysis)
- **Transport reactions** (inter-compartment exchange)
- **Demand reactions** (protein turnover, flux testing)

### Implementation Results (Cytoskeleton Example)

- **Base Model**: endoA_250917_3.json (16,670 metabolites, 19,300 reactions)
- **Modified Model**: +38 metabolites, +46 reactions
- **Cytoskeleton Compartment**: Added as `[ck]`
- **Metabolite Database**: 30/30 metabolites matched (100%)
- **Validation**: ✅ 9/9 tests passed
- **Demand Reactions**: ✅ 4/4 functional (100% flux capability)
- **Visualizations**: ✅ 3 types generated (D3.js, Plotly, GraphViz)

> **📊 Validation reports available in `reports/` folder**
> 
> **🎨 Interactive visualizations in `figures/` folder** - Open HTML files in browser

---

## ✨ Features

### Core Capabilities
- ✅ **Generic pathway builder** - Add ANY pathway/compartment (not just cytoskeleton!)
- ✅ **Annotation-based metabolite matching** - Robust identification using multiple databases
- ✅ **Universal reaction format** - Human-readable equations that work for any network
- ✅ **Fully configurable** - All metabolites and reactions defined in JSON config files
- ✅ **Multi-compartment support** - Easy to add new compartments
- ✅ **Automatic ID generation** - Smart metabolite and reaction ID assignment
- ✅ **Comprehensive validation** - 9 validation tests including FBA, mass balance, topology
- ✅ **Rich visualizations** - D3.js interactive, Plotly with FVA, GraphViz static images
- ✅ **Detailed reporting** - Pathway-specific implementation and validation reports

### Validation Suite (9 Tests)
1. **Quick Validation** - Model structure + FBA optimization
2. **Unit Tests** - Metabolite matching functionality
3. **Reproducibility Tests** - Deterministic results across runs
4. **Integration Tests** - Complete pipeline verification
5. **Mass Balance Check** - Elemental balance validation
6. **Network Topology** - Connectivity analysis (orphaned metabolites, dead-ends)
7. **Database Consistency** - Annotation coverage verification
8. **Synthesis Balance** - Protein/polymer synthesis validation
9. **Demand Reactions** - Flux capability testing with FVA

### Visualization Types
1. **D3.js Interactive HTML** - Force-directed network with zoom/pan
2. **Plotly Interactive** - Network with FVA analysis and flux ranges
3. **GraphViz Static** - Publication-quality PNG/SVG/PDF images

---

## 💻 Installation

### Directory Structure Requirements

This module requires the following parent-level directories:
```
/
├── models/          # Shared metabolic models
├── functions/       # Shared core implementation modules
└── implement_pathway/  # This module
```

Scripts automatically add the parent directory to Python's path:
```python
sys.path.insert(0, str(Path(__file__).parent.parent))
from functions import pathway_builder
```

Config files use paths like `models/endoA_250917_3.json` which are resolved relative to the project root by `functions/config.py`.

### Requirements
- Python 3.6+
- pip package manager

### Install Dependencies

Navigate to the `implement_pathway/` directory and install dependencies:

```bash
cd implement_pathway
pip install -r requirements.txt
```

Required packages:
```
cobra>=0.26.0
pandas>=1.3.0
numpy>=1.21.0
matplotlib>=3.4.0
plotly>=5.0.0
networkx>=2.6.0
graphviz>=0.16.0
```

---

## 🚀 Quick Start

**Note**: All commands should be run from the `implement_pathway/` directory.

### 1. Run the Complete Pipeline

The main pipeline script supports any pathway configuration:

```bash
cd implement_pathway

# Using a specific config file (recommended)
python3 pathway_implementation.py --config config/config_cytoskeleton.json

# Or use shorthand
python3 pathway_implementation.py -c config/config_Syndecan1_HS.json
```

This will:
1. Load the specified configuration file
2. Analyze model annotations (from `../models/`)
3. Build metabolite ID database
4. Add compartments defined in config
5. Create metabolites for each compartment
6. Create reactions for each compartment
7. Save modified model to `../models/` with pathway-specific name
8. Generate implementation report in `reports/`

### 2. Validate the Implementation (All 9 Tests)

```bash
# Run all validations with specific pathway
python3 validate.py 0 --config config/config_cytoskeleton.json --model ../models/endoA_250917_3_with_Cytoskeleton.json

# Run a specific validation (1-9)
python3 validate.py 5 --config config/config_cytoskeleton.json --model ../models/endoA_250917_3_with_Cytoskeleton.json
```

Available validation options:
- `0`: Run ALL validations (generates detailed report)
- `1`: Quick validation (FBA + structure)
- `2`: Unit tests (metabolite matching)
- `3`: Reproducibility tests
- `4`: Integration tests
- `5`: Mass balance check
- `6`: Network topology
- `7`: Database consistency
- `8`: Synthesis balance
- `9`: Demand reaction functionality (FVA)

Validation reports saved to `reports/YYYY-MM-DD_report_{pathway}_validation.txt`

### 3. Visualize the Network

```bash
# Generate all visualizations
python3 visualize.py 0 --config config/config_cytoskeleton.json --model ../models/endoA_250917_3_with_Cytoskeleton.json

# Generate specific visualization type
python3 visualize.py 1 --config config/config_cytoskeleton.json --model ../models/endoA_250917_3_with_Cytoskeleton.json  # D3.js only
python3 visualize.py 2 --config config/config_cytoskeleton.json --model ../models/endoA_250917_3_with_Cytoskeleton.json  # Plotly only
python3 visualize.py 3 --config config/config_cytoskeleton.json --model ../models/endoA_250917_3_with_Cytoskeleton.json  # GraphViz only
```

This generates (option 0):
1. **D3.js interactive network** - `figures/{Pathway}_network_d3.html`
2. **Plotly with FVA analysis** - `figures/{Pathway}_network_plotly.html`
3. **GraphViz static images** - `figures/{Pathway}_network.png/svg/pdf`

All visualizations saved to `figures/` directory with pathway-specific names.

**Visualization Features:**
- **GraphViz**: Optimized for widescreen monitors with horizontal layout (left-to-right)
- **Reaction Categories**: 10 distinct categories with color coding (Synthesis, Demand, AA Transport, etc.)
- **Dynamic Legend**: Shows only categories present in the network
- **Edge Labels**: Display stoichiometric coefficients (numbers on edges)
- **Flexibility**: Bold borders = flexible reactions, Dashed = blocked reactions

### 4. Implement Multiple Pathways Sequentially

```bash
# First pathway: Cytoskeleton
python3 pathway_implementation.py -c config/config_cytoskeleton.json

# Second pathway: Syndecan1_HS (builds on cytoskeleton model)
# Update config to use previous output as input
python3 pathway_implementation.py -c config/config_Syndecan1_HS.json
```

---

## ⚙️ Configuration

All pathway configurations are stored in JSON files in the `config/` folder. Each config file specifies compartments, metabolites, and reactions for a pathway.

### Configuration File Structure

```json
{
  "pathway_name": "Cytoskeleton",
  "model": {
    "base": "models/endoA_250917_3.json",
    "output": "models/endoA_250917_3_with_cytoskeleton.json"
  },
  "visualizations": {
    "output_directory": "figures"
  },
  "compartments": [
    {
      "name": "cytoskeleton",
      "abbreviation": "ck"
    }
  ],
  "metabolites": {
    "database_file": "data/metabolite_id_database.json",
    "targets": ["ATP", "ADP", "GTP", "GDP", "L-alanine", ...],
    "cytoskeleton_specific": {
      "actin protein (bare)": {
        "charge": -18,
        "formula": "C1946H3018N528O610S32",
        "uniprot": "P60709",
        "sbo": "SBO:0000247"
      }
    }
  },
  "reactions": [...]
}
```

### Model Paths

Configuration files specify model paths relative to the project root:

```json
{
  "model": {
    "base": "models/endoA_250917_3.json",
    "output": "models/endoA_250917_3_with_cytoskeleton.json"
  }
}
```

**Note**: Paths use `models/` (not `../models/`) as they are resolved relative to the project root by `functions/config.py`.

**For sequential pathways**, chain the output of one as the input to the next:

```json
{
  "model": {
    "base": "models/endoA_250917_3_with_cytoskeleton.json",
    "output": "models/endoA_250917_3_with_cytoskeleton_and_glycocalyx.json"
  }
}
```

### Compartments

```json
{
  "compartments": [
    {
      "name": "cytoskeleton",
      "abbreviation": "ck"
    }
  ]
}
```

### Metabolites

**Standard metabolites** (found via annotations in base model):
```json
{
  "metabolites": {
    "targets": ["ATP", "ADP", "GTP", "GDP", "L-alanine", "L-arginine", ...]
  }
}
```

**Pathway-specific metabolites** (proteins, modified species not in base model):
```json
{
  "metabolites": {
    "cytoskeleton_specific": {
      "actin protein (bare)": {
        "charge": -18,
        "formula": "C1946H3018N528O610S32",
        "uniprot": "P60709",
        "sbo": "SBO:0000247"
      }
    }
  }
}
```

### Reactions

Reactions use a **universal, human-readable format**:

```json
{
  "reactions": [
    {
      "id": "R_Synth_Actin",
      "name": "Synthesis of actin protein in cytoskeleton",
      "description": "Polymerization of 395 amino acids into actin protein",
      "equation": "22 L-alanine[ck] + 18 L-arginine[ck] + ... + 1576 ATP[ck] --> 1 actin protein (bare)[ck] + 1576 ADP[ck] + ...",
      "subsystem": "Cytoskeleton synthesis",
      "ec": "6.1.-.-",
      "gpr": "ENSG00000075624",
      "lower_bound": 0.0,
      "upper_bound": 1000.0,
      "sbo": "SBO:0000176"
    }
  ]
}
```

**Equation syntax:**
```
[stoich] metabolite[compartment] + [stoich] metabolite[compartment] --> [stoich] metabolite[compartment]
```

**Arrow types:**
- `-->` for irreversible (forward) reactions
- `<=>` for reversible reactions

**Key fields:**
- **id** (required): Unique identifier
- **name** (required): Human-readable name
- **equation** (required): Reaction equation with stoichiometry and compartments
- **subsystem** (optional): Pathway or subsystem
- **lower_bound** (required): Minimum flux (negative for reversible)
- **upper_bound** (required): Maximum flux

### Cytoskeleton Reaction Categories

| Category | Count | Description |
|----------|-------|-------------|
| Protein Synthesis | 4 | Actin, tubulin, vimentin, myosin |
| Protein Modification | 6 | Glycosylation, phosphorylation |
| Energy Dissipation | 2 | ATP/GTP hydrolysis |
| Demand Reactions | 4 | Protein turnover |
| Transport | 30 | Cytosol ↔ cytoskeleton |
| **Total** | **46** | |

---

## 💻 Usage

### Command-Line Usage

The `pathway_implementation.py` script supports flexible configuration file selection:

```bash
# Change to the implement_pathway directory first
cd implement_pathway

# Use a specific config file from the config folder
python3 pathway_implementation.py --config config/config_glycocalyx.json

# Short form
python3 pathway_implementation.py -c config/config_cytoskeleton.json
```

### Sequential Pathway Implementation

To implement multiple pathways sequentially (e.g., cytoskeleton, then glycocalyx):

```bash
# Step 1: Implement cytoskeleton (creates model with cytoskeleton)
python3 pathway_implementation.py -c config/config_cytoskeleton.json

# Step 2: Implement glycocalyx on top of cytoskeleton model
# Note: config_glycocalyx.json should have base model pointing to 
#       models/endoA_250917_3_with_cytoskeleton.json
python3 pathway_implementation.py -c config/config_glycocalyx.json

# Step 3: Implement another pathway on top of combined model
python3 pathway_implementation.py -c config/config_another_pathway.json
```

### Available Config Files

Config files are organized in the `config/` folder:

- `config/config_cytoskeleton.json` - Cytoskeleton pathway implementation
  - Base: `models/endoA_250917_3.json`
  - Output: `models/endoA_250917_3_with_cytoskeleton.json`
  - Compartments: [ck]
  - Reactions: 46 cytoskeleton reactions

- `config/config_glycocalyx.json` - Glycocalyx pathway implementation
  - Base: `models/endoA_250917_3_with_cytoskeleton.json`
  - Output: `models/endoA_250917_3_with_cytoskeleton_and_glycocalyx.json`
  - Compartments: [gl, g, e]
  - Reactions: 33 glycocalyx reactions

### Programmatic Usage (Python API)

**Important**: When using the API programmatically, ensure you're importing from the `functions` module in the parent directory:

```python
# Add parent directory to path if needed
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from functions import pathway_builder
from functions.config import load_config
import json

# Load config from specific file
config = load_config('implement_pathway/config/config_cytoskeleton.json')

# Load model from parent models directory
with open('../models/endoA_250917_3.json', 'r') as f:
    model = json.load(f)

# Add cytoskeleton compartment
model, ck_abbrev, already_exists = pathway_builder.add_compartment(
    model, 'ck', 'cytoskeleton'
)

# Load metabolite database
with open('implement_pathway/data/metabolite_id_database.json', 'r') as f:
    id_database = json.load(f)

# Add metabolites
met_count = pathway_builder.create_compartment_metabolites(
    model, id_database, ck_abbrev, config, 'cytoskeleton_specific'
)

# Add reactions
rxn_count = pathway_builder.create_compartment_reactions(
    model, id_database, ck_abbrev, config
)

print(f"Added {met_count} metabolites and {rxn_count} reactions")

# Save to parent models directory
with open('../models/endoA_250917_3_with_cytoskeleton.json', 'w') as f:
    json.dump(model, f, indent=2)
```

### Advanced: Custom Metabolite Database

```python
# Add parent directory to path for imports
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from functions import analyze_annotations, build_id_database
from functions.config import load_config

# Load config
config = load_config('implement_pathway/config/config_cytoskeleton.json')

# Extract annotations from model in parent directory
model_path = config['model']['base']  # e.g., 'models/endoA_250917_3.json'
targets = config['metabolites']['targets']

metabolite_annotations, not_found = analyze_annotations.extract_metabolite_annotations(
    model_path, targets, verbose=True
)

# Build and save database
database = build_id_database.build_database_from_annotations(metabolite_annotations)
output_path = 'implement_pathway/data/metabolite_id_database.json'
build_id_database.save_database(database, output_path)

print(f"Database saved with {len(database)} metabolites")
```

### Validation with COBRApy

```python
import cobra

# Load model from parent models directory
model = cobra.io.load_json_model('../models/endoA_250917_3_with_Cytoskeleton.json')

# Run FBA
solution = model.optimize()
print(f"Growth rate: {solution.objective_value}")

# Check cytoskeleton reactions
ck_reactions = [r for r in model.reactions if any('ck' in m.id for m in r.metabolites)]
print(f"Cytoskeleton reactions: {len(ck_reactions)}")

# Test demand reactions
for rxn in model.reactions:
    if rxn.id.startswith('DM_'):
        print(f"{rxn.id}: bounds [{rxn.lower_bound}, {rxn.upper_bound}]")
```

---

## 📁 Project Structure

```
/
├── models/ ⭐                                # Metabolic models (shared parent directory)
│   ├── endoA_250917_3.json                 # Base GEM model
│   ├── endoA_250917_3_with_Cytoskeleton.json         # With cytoskeleton
│   └── endoA_250917_3_with_Syndecan1_HS.json        # With Syndecan1_HS
│
├── functions/ ⭐                             # Core implementation modules (shared parent directory)
│   ├── __init__.py                         # Package initialization
│   ├── config.py                           # Configuration management
│   ├── pathway_builder.py                  # Generic pathway builder
│   ├── add_reaction.py                     # Reaction utilities
│   ├── analyze_annotations.py              # Annotation extraction
│   ├── build_id_database.py                # Database builder
│   └── ... (other utility functions)
│
└── implement_pathway/                       # This repository/module ⭐
    ├── README.md                           # This file - complete documentation
    ├── requirements.txt                    # Python dependencies
    │
    ├── pathway_implementation.py ⭐         # Main implementation script
    ├── validate.py ⭐                       # Validation suite (9 tests)
    ├── visualize.py ⭐                      # Visualization generator (3 types)
    │
    ├── config/                             # Configuration files (pathway definitions)
    │   ├── config_cytoskeleton.json ⭐     # Cytoskeleton pathway config
    │   ├── config_Syndecan1_HS.json       # Syndecan1/Heparan Sulfate config
    │   └── config_glycocalyx.json         # Glycocalyx pathway config (template)
    │
    ├── data/                               # Generated data files
    │   └── metabolite_id_database.json    # Auto-generated metabolite database
    │
    ├── reports/ ⭐                          # Generated validation reports
    │   ├── YYYY-MM-DD_report_Cytoskeleton_implementation.txt
    │   ├── YYYY-MM-DD_report_Cytoskeleton_validation.txt
    │   ├── YYYY-MM-DD_report_Syndecan1_HS_implementation.txt
    │   └── YYYY-MM-DD_report_Syndecan1_HS_validation.txt
    │
    ├── figures/ ⭐                          # Generated visualizations
    │   ├── Cytoskeleton_network_d3.html   # D3.js interactive
    │   ├── Cytoskeleton_network_plotly.html # Plotly with FVA
    │   ├── Cytoskeleton_network.png       # GraphViz PNG
    │   ├── Cytoskeleton_network.svg       # GraphViz SVG
    │   ├── Cytoskeleton_network.pdf       # GraphViz PDF
    │   └── ... (other pathway visualizations)
    │
    ├── examples/                           # Usage examples
    │   ├── README.md                      # Examples documentation
    │   ├── basic_integration.py           # Simple integration example
    │   └── custom_validation.py           # Custom validation example
    │
    ├── docs/                               # Additional documentation
    └── dev/                                # Development files and archives

```

⭐ = Most critical files for using the pipeline

**Key Points:**
- The `implement_pathway/` subdirectory contains this module's main scripts and configuration
- `functions/` and `models/` are in the parent directory and shared across modules
- Configuration files use relative paths (`models/` not `../models/`) resolved from project root
- Each pathway has its own config file in `config/`
- Reports and figures are pathway-specific (include pathway name)
- Models are saved with pathway names in filename

---

## 📚 Examples

### Example 1: Complete Automated Workflow

Run the complete pipeline (implementation → validation → visualization) with a single script:

```bash
python examples/complete_workflow.py
```

This automated script will:
1. Run `pathway_implementation.py` (add pathways)
2. Run `validate.py` (verify implementation)
3. Run `visualize.py` (generate network visualizations)
4. Verify all outputs were created successfully

Main pipeline outputs saved to:
- `model/endoA_250917_3_with_Cytoskeleton.json` (main model)
- `reports/YYYY-MM-DD_report_Cytoskeleton_implementation.txt`
- `reports/YYYY-MM-DD_report_Cytoskeleton_validation.txt`
- `figures/Cytoskeleton_network_*.html/png/svg/pdf`

### Example 1b: Manual Step-by-Step

Or run each step manually:

```bash
# Step 1: Implement pathway
python3 pathway_implementation.py -c config/config_cytoskeleton.json

# Step 2: Validate (all 9 tests)
python3 validate.py 0 --config config/config_cytoskeleton.json --model model/endoA_250917_3_with_Cytoskeleton.json

# Step 3: Generate visualizations (all 3 types)
python3 visualize.py 0 --config config/config_cytoskeleton.json --model model/endoA_250917_3_with_Cytoskeleton.json
```

### Example 2: Programmatic Integration

See `examples/basic_integration.py` for usage. Here's a generic pathway integration workflow:

```python
#!/usr/bin/env python3
"""Add any pathway to a model programmatically"""

from functions.config import load_config
from functions import pathway_builder
import json

def main():
    # Load configuration
    config = load_config('config/config_cytoskeleton.json')
    
    # Load model
    with open(config['model']['base'], 'r') as f:
        model = json.load(f)
    
    # Load metabolite database (pre-generated)
    with open('data/metabolite_id_database.json', 'r') as f:
        id_database = json.load(f)
    
    # Check if pathway already exists
    pathway_check = pathway_builder.check_pathway_exists(model, 'ck')
    if pathway_check['exists']:
        print("⚠ Pathway already exists!")
        return
    
    # Add compartment
    compartments = config['compartments']
    abbrev = compartments[0]['abbreviation']
    name = compartments[0]['name']
    
    model, abbrev_used, exists = pathway_builder.add_compartment(
        model, abbrev, name
    )
    
    # Add metabolites (shared + pathway-specific)
    met_count = pathway_builder.create_compartment_metabolites(
        model, id_database, abbrev_used, config, 'cytoskeleton_specific'
    )
    
    # Add reactions
    rxn_count = pathway_builder.create_compartment_reactions(
        model, id_database, abbrev_used, config
    )
    
    # Save
    output_path = config['model']['output']
    with open(output_path, 'w') as f:
        json.dump(model, f, indent=2)
    
    print(f"✓ Complete: {met_count} metabolites, {rxn_count} reactions")
    print(f"✓ Saved to: {output_path}")

if __name__ == '__main__':
    main()
```

### Example 3: Custom Validation

```python
#!/usr/bin/env python3
"""Custom validation with COBRApy"""

import cobra

# Load the enhanced model from parent directory
model = cobra.io.load_json_model('../models/endoA_250917_3_with_Cytoskeleton.json')

print(f"Model: {len(model.metabolites)} metabolites, {len(model.reactions)} reactions")

# Test 1: Basic FBA
solution = model.optimize()
print(f"\nTest 1 - FBA: {solution.status}")
print(f"  Objective value: {solution.objective_value:.6f}")

# Test 2: Protein synthesis reactions
print("\nTest 2 - Protein Synthesis:")
protein_reactions = {
    'MAR24020': 'Actin',
    'MAR24021': 'Tubulin',
    'MAR24022': 'Vimentin',
    'MAR24023': 'Myosin'
}

for rxn_id, protein_name in protein_reactions.items():
    rxn = model.reactions.get_by_id(rxn_id)
    rxn.lower_bound = 0.01  # Force synthesis
    
    solution = model.optimize()
    flux = solution.fluxes[rxn_id]
    
    status = "✓" if solution.status == 'optimal' and flux >= 0.01 else "✗"
    print(f"  {status} {protein_name}: flux={flux:.6f}")
    
    rxn.lower_bound = 0  # Reset

# Test 3: Mass balance
print("\nTest 3 - Mass Balance:")
from cobra.manipulation import check_mass_balance

ck_reactions = [r for r in model.reactions if 'ck' in r.id.lower()]
imbalanced = []

for rxn in ck_reactions:
    balance = check_mass_balance(rxn)
    if any(abs(v) > 1e-6 for v in balance.values()):
        imbalanced.append(rxn.id)

print(f"  Balanced: {len(ck_reactions) - len(imbalanced)}/{len(ck_reactions)}")
if imbalanced:
    print(f"  ✗ Imbalanced reactions: {', '.join(imbalanced)}")
else:
    print(f"  ✓ All cytoskeleton reactions are mass balanced")

print("\n✅ Validation complete")
```

For more examples, see `examples/basic_integration.py` and `examples/custom_validation.py`.

---

## 🌐 Extending to Other Networks

This framework is **pathway-agnostic** and can be easily extended to any network (glycocalyx, ECM, nucleus, etc.):

### 1. Create New Config File

Copy an existing config and modify:
```bash
cp config/config_cytoskeleton.json config/config_mynetwork.json
```

### 2. Update Pathway Name and Model Paths

```json
{
  "pathway_name": "MyNetwork",
  "model": {
    "base": "model/endoA_250917_3.json",
    "output": "model/endoA_250917_3_with_MyNetwork.json"
  }
}
```

### 3. Define Compartments

```json
{
  "compartments": [
    {
      "name": "mynetwork",
      "abbreviation": "mn",
      "existing": false
    }
  ]
}
```

### 4. Add Network-Specific Metabolites

```json
{
  "metabolites": {
    "targets": [
      "ATP", "ADP", "H2O", "Pi"
    ],
    "mynetwork_specific": {
      "my_protein": {
        "charge": -5,
        "formula": "C100H150N30O40S2",
        "uniprot": "P12345",
        "sbo": "SBO:0000247"
      }
    }
  }
}
```

### 5. Add Reactions

```json
{
  "reactions": [
    {
      "id": "R_My_Reaction",
      "name": "My network reaction",
      "description": "Description of what this reaction does",
      "equation": "1 ATP[mn] + 1 H2O[mn] --> 1 ADP[mn] + 1 Pi[mn]",
      "subsystem": "MyNetwork metabolism",
      "ec": "3.6.1.3",
      "gpr": "ENSG00000123456",
      "lower_bound": 0.0,
      "upper_bound": 1000.0,
      "sbo": "SBO:0000176"
    }
  ]
}
```

### 6. Run Pipeline

```bash
# Implement the pathway
python3 pathway_implementation.py -c config/config_mynetwork.json

# Validate (all 9 tests)
python3 validate.py 0 --config config/config_mynetwork.json --model model/endoA_250917_3_with_MyNetwork.json

# Visualize
python3 visualize.py 0 --config config/config_mynetwork.json --model model/endoA_250917_3_with_MyNetwork.json
```

The pipeline automatically:
- Creates new compartments
- Finds shared metabolites via annotations
- Creates pathway-specific metabolites
- Parses and adds reactions from equations
- Validates mass balance
- Tests flux capability
- Generates visualizations

---

## 🔧 Troubleshooting

### Issue: Config file not found

```bash
# Error: Configuration file not found: config/myconfig.json

# Solution: Ensure the file exists in the config folder
ls config/

# Or provide the full path
python3 pathway_implementation.py -c config/config_cytoskeleton.json
```

### Issue: Import errors (functions module not found)

```bash
# Error: ModuleNotFoundError: No module named 'functions'

# Solution: Ensure you're in the repository root
cd /path/to/implement_pathway

# Install dependencies
pip install -r requirements.txt
```

### Issue: Model file not found

```bash
# Verify model path in config file is correct
cat config/config_cytoskeleton.json | grep "base"

# Ensure the model file exists
ls -lh model/endoA_250917_3.json
````

## 🔧 Troubleshooting

### Issue: Config file not found

```bash
# Error: Configuration file not found: config/myconfig.json

# Solution: Ensure you're in the implement_pathway directory
cd implement_pathway

# List available config files
ls config/

# Use the correct path
python3 pathway_implementation.py -c config/config_cytoskeleton.json
```

### Issue: Import errors (functions module not found)

```bash
# Error: ModuleNotFoundError: No module named 'functions'

# Solution: Ensure the parent directory structure is correct
# The functions/ folder should be at the parent level:
# /
# ├── functions/
# └── implement_pathway/

# From implement_pathway directory, check parent:
ls ../functions/

# The scripts automatically add parent to path:
# sys.path.insert(0, str(Path(__file__).parent.parent))
```

### Issue: Model file not found

```bash
# Error: Model file not found: models/endoA_250917_3.json

# Verify model path in config file is correct
cat config/config_cytoskeleton.json | grep "base"

# Ensure the model file exists in parent models directory
ls ../models/endoA_250917_3.json

# Config files should use paths like "models/..." (not "../models/")
# These are resolved relative to project root by functions/config.py
```

### Issue: Invalid JSON in config file

```bash
# Error: Invalid JSON in config file

# Solution: Validate JSON syntax
python3 -m json.tool config/config_cytoskeleton.json
```

### Issue: Sequential pathway base model not found

```bash
# Error: Model file not found when running glycocalyx

# Solution: Run cytoskeleton first to create the base model
python3 pathway_implementation.py -c config/config_cytoskeleton.json

# Verify the output was created
ls ../models/endoA_250917_3_with_cytoskeleton.json

# Then run glycocalyx (which should reference the cytoskeleton output)
python3 pathway_implementation.py -c config/config_glycocalyx.json
```

### Issue: Metabolite not found

Check that metabolite names in equations match exactly those in `config.json`:
```python
import json
config = json.load(open('config.json'))
print("Available metabolites:")
print(config['metabolites']['targets'])
print(config['metabolites']['cytoskeleton_specific'].keys())
```

### Issue: Different results

```bash
# Verify Python version
python --version  # Should be 3.6+

# Verify COBRApy version
python -c "import cobra; print(cobra.__version__)"  # Should be >= 0.26.0
```

---

## 📝 License

This project is part of metabolic modeling research.

---

**Last Updated**: November 2025
