# Generate Human Metabolic Database

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![COBRApy](https://img.shields.io/badge/COBRApy-0.29+-green.svg)](https://opencobra.github.io/cobrapy/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](../LICENSE)

## Overview

This module provides a comprehensive pipeline for generating genome-scale metabolic models (GEMs) from KEGG pathway data. It automatically fetches reactions and metabolites from KEGG, resolves compound identities, calculates molecular formulas, performs mass balance validation, and exports SBML-compatible models for use with COBRApy and other metabolic modeling tools.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Architecture](#pipeline-architecture)
- [File Descriptions](#file-descriptions)
- [Configuration](#configuration)
- [Usage Examples](#usage-examples)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [API Reference](#api-reference)
- [Contributing](#contributing)

## Features

### Core Capabilities

- **KEGG Integration**: Automated fetching of reactions, compounds, and pathways from KEGG REST API
- **Compound Identification**: Multi-source metabolite identification using PubChem, KEGG, and custom databases
- **Formula Resolution**: Intelligent formula parsing including glycan composition calculation
- **Mass Balance Validation**: Automatic detection and reporting of mass-imbalanced reactions
- **SBML Export**: Export to SBML format compatible with COBRApy, Escher, and other tools
- **Gene Associations**: Integration of gene-protein-reaction (GPR) rules
- **Compartmentalization**: Support for multiple cellular compartments

### Recent Improvements (v2.0)

1. **Fixed Mass Balance Counting**: Corrected false positive reporting in mass balance validation
2. **KEGG URL Parsing Update**: Adapted to new KEGG stoichiometry URL format (`map` → `R`)
3. **Glycan Formula Support**: Added fallback mechanism for glycan formulas and composition-based calculation

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager
- Internet connection (for KEGG API access)

### Setup

1. Clone the repository:
```bash
git clone https://github.com/MarindeMasLab/THG_The-Human-GEM-protocol.git
cd THG_The-Human-GEM-protocol
```

2. Create a virtual environment:
```bash
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r generate_data-base/requirements.txt
```

### Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| cobra | ≥0.29.0 | Metabolic model manipulation and SBML I/O |
| pandas | ≥2.0.0 | Data manipulation and analysis |
| numpy | ≥1.24.0 | Numerical computations |
| scipy | ≥1.10.0 | Scientific computing |
| requests | ≥2.28.0 | HTTP requests to KEGG API |
| pubchempy | ≥1.0.4 | PubChem compound identification |
| tqdm | ≥4.65.0 | Progress bars |
| dill | ≥0.3.6 | Extended pickling support |
| xlsxwriter | ≥3.0.0 | Excel report generation |
| sympy | ≥1.12 | Symbolic mathematics for formula parsing |
| python-dotenv | ≥1.0.0 | Environment variable management |

## Quick Start

### Basic Usage

```python
from generate_db import run_pipeline

# Run the full pipeline with default settings
run_pipeline(
    pathway_file='../files/human_kegg_pathways.txt',
    output_dir='../models/',
    model_name='Human_GEM'
)
```

### Minimal Example

```python
from generate_db import (
    KEGGPWYS, 
    rxnsFromPathway, 
    getGeneAssoc,
    modelExchange
)

# Load pathways
kegg_pws = KEGGPWYS(rxnsFile='../files/pathway_reactions.tsv')

# Get reactions for a specific pathway
reactions = rxnsFromPathway('hsa00010')  # Glycolysis

# Build exchange reactions
exchanges = modelExchange(reactions)
```

## Pipeline Architecture

The database generation follows a multi-stage pipeline:

```
┌─────────────────────────────────────────────────────────────────┐
│                    STAGE 1: PATHWAY LOADING                     │
│  Load KEGG pathways → Extract reaction IDs → Cache results      │
└─────────────────────────────────────────────────────────────────┘
                                 ↓
┌─────────────────────────────────────────────────────────────────┐
│                  STAGE 2: REACTION PROCESSING                   │
│  Fetch reaction data → Parse equations → Extract compounds      │
└─────────────────────────────────────────────────────────────────┘
                                 ↓
┌─────────────────────────────────────────────────────────────────┐
│                 STAGE 3: COMPOUND IDENTIFICATION                │
│  KEGG lookup → PubChem search → Formula resolution → Glycans    │
└─────────────────────────────────────────────────────────────────┘
                                 ↓
┌─────────────────────────────────────────────────────────────────┐
│                   STAGE 4: MASS BALANCE CHECK                   │
│  Parse formulas → Count atoms → Validate balance → Report       │
└─────────────────────────────────────────────────────────────────┘
                                 ↓
┌─────────────────────────────────────────────────────────────────┐
│                    STAGE 5: MODEL ASSEMBLY                      │
│  Create COBRA model → Add reactions → Set bounds → Export SBML  │
└─────────────────────────────────────────────────────────────────┘
```

### Stage Details

#### Stage 1: Pathway Loading
- Reads pathway IDs from input file (e.g., `human_kegg_pathways.txt`)
- Queries KEGG REST API for pathway composition
- Extracts reaction IDs associated with each pathway
- Caches results to avoid redundant API calls

#### Stage 2: Reaction Processing
- Fetches detailed reaction data from KEGG
- Parses reaction equations to extract substrates and products
- Handles stoichiometry coefficients (including from HTML/image sources)
- Resolves EC numbers and pathway associations

#### Stage 3: Compound Identification
- Multi-source compound identification:
  - Primary: KEGG compound database
  - Secondary: PubChem by name/formula
  - Fallback: Custom compound database
- Formula resolution for glycans via COMPOSITION parsing
- Charge and compartment assignment

#### Stage 4: Mass Balance Validation
- Parses molecular formulas (including complex glycan formulas)
- Counts atoms on both sides of reactions
- Identifies and reports mass-imbalanced reactions
- Handles special cases (polymers, generic compounds)

#### Stage 5: Model Assembly
- Creates COBRApy Model object
- Adds metabolites with annotations (KEGG, PubChem, InChI)
- Adds reactions with GPR rules
- Sets reaction bounds (reversibility)
- Exports to SBML format

## File Descriptions

### Core Files

| File | Description |
|------|-------------|
| `generate_db.py` | Main pipeline script - orchestrates the entire database generation process |
| `../functions/function_bm_gdb.py` | KEGG API functions, formula parsing, and compound parameter extraction |
| `../functions/equations_bm_gdb.py` | Reaction equation parsing and stoichiometry extraction |
| `../functions/functions_mass_balance.py` | Mass balance validation and atom counting |
| `../functions/class_generate_database.py` | Core classes for database generation |

### Support Files

| File | Description |
|------|-------------|
| `../functions/function_metabolite_identification.py` | PubChem integration and metabolite lookup |
| `../functions/function_reac_identification.py` | Reaction identification utilities |
| `../functions/ensembl_client.py` | Ensembl API client for gene annotations |
| `../functions/error_tracker.py` | Error logging and tracking utilities |

### Input Files (in `../files/`)

| File | Description |
|------|-------------|
| `human_kegg_pathways.txt` | List of KEGG pathway IDs to process |
| `compartments_info.txt` | Cellular compartment definitions |
| `extra_compounds.txt` | Custom compound definitions |
| `extra_formula.txt` | Manual formula overrides |
| `special_compounds.txt` | Special handling rules for specific compounds |

## Configuration

### Environment Variables

Create a `.env` file in the project root:

```env
# KEGG API settings
KEGG_DELAY=0.5          # Delay between API calls (seconds)
KEGG_MAX_RETRIES=3      # Maximum retry attempts for failed requests

# Output settings
OUTPUT_FORMAT=sbml      # Output format (sbml, json, mat)
VERBOSE=true            # Enable verbose logging
```

### Configuration Options

```python
# In generate_db.py
CONFIG = {
    'kegg_delay': 0.5,           # API rate limiting
    'max_retries': 3,            # Failed request retries
    'cache_enabled': True,       # Enable result caching
    'mass_balance_strict': False, # Strict mass balance mode
    'include_orphan_rxns': True,  # Include orphan reactions
    'compartments': ['c', 'm', 'e', 'n', 'r', 'g', 'x', 'l'],
}
```

## Usage Examples

### Example 1: Full Pipeline Execution

```python
from generate_db import main

# Run complete pipeline
main(
    pathway_file='../files/human_kegg_pathways.txt',
    output_prefix='Human_GEM_v1',
    skip_cache=False,
    verbose=True
)
```

### Example 2: Single Pathway Processing

```python
from generate_db import rxnsFromPathway, processReactions

# Get reactions from glycolysis pathway
glycolysis_rxns = rxnsFromPathway('hsa00010')
print(f"Found {len(glycolysis_rxns)} reactions in glycolysis")

# Process reactions
processed = processReactions(glycolysis_rxns)
```

### Example 3: Compound Information Retrieval

```python
from functions.function_bm_gdb import getCompParamFromRestAPI

# Get compound parameters from KEGG
compound_data = getCompParamFromRestAPI('C00001')  # Water
print(f"Name: {compound_data['NAME']}")
print(f"Formula: {compound_data['FORMULA']}")
```

### Example 4: Mass Balance Check

```python
from functions.functions_mass_balance import check_mass_balance

# Check mass balance for a reaction
result = check_mass_balance(
    substrates={'C00001': 1, 'C00002': 1},  # H2O + ATP
    products={'C00008': 1, 'C00009': 1},    # ADP + Pi
    formulas={'C00001': 'H2O', 'C00002': 'C10H16N5O13P3', 
              'C00008': 'C10H15N5O10P2', 'C00009': 'H3O4P'}
)
print(f"Mass balanced: {result['balanced']}")
```

### Example 5: Glycan Formula Calculation

```python
from functions.function_bm_gdb import calculate_glycan_formula

# Calculate formula from composition
composition = "(Glc)3 (GlcNAc)2 (Man)9"
formula = calculate_glycan_formula(composition)
print(f"Glycan formula: {formula}")  # C80H132N2O61
```

## Output Files

### Generated Models

| File Pattern | Description |
|--------------|-------------|
| `Human_database_YYYYMMDD.xml` | Full SBML model |
| `Human_database_YYYYMMDD.json` | JSON format model |
| `Human_database_YYYYMMDD_report.xlsx` | Validation report |

### Report Contents

The Excel report (`*_report.xlsx`) contains:

1. **Summary Sheet**: Overall statistics
   - Total reactions, metabolites, genes
   - Mass balance pass/fail counts
   - Pathway coverage

2. **Reactions Sheet**: Detailed reaction information
   - Reaction ID, name, equation
   - EC number, pathway associations
   - Mass balance status

3. **Metabolites Sheet**: Compound information
   - KEGG ID, name, formula
   - Charge, compartment
   - External database cross-references

4. **Errors Sheet**: Issues encountered
   - Failed API calls
   - Unparseable formulas
   - Missing data warnings

## Troubleshooting

### Common Issues

#### 1. KEGG API Connection Errors

```
Error: Failed to fetch data from KEGG REST API
```

**Solution**: Check internet connection and increase retry count:
```python
CONFIG['max_retries'] = 5
CONFIG['kegg_delay'] = 1.0  # Increase delay
```

#### 2. Missing Formulas for Glycans

```
Warning: Empty formula for glycan G00123
```

**Solution**: The pipeline now automatically calculates formulas from COMPOSITION. If still missing, add to `../files/extra_formula.txt`:
```
G00123  C40H66N2O31
```

#### 3. Mass Balance Failures

```
Warning: Reaction R00001 is mass-imbalanced
```

**Debugging steps**:
1. Check if all compound formulas are resolved
2. Verify stoichiometry coefficients are correct
3. Look for generic compounds (R groups)

#### 4. Memory Issues with Large Models

```
MemoryError: Unable to allocate array
```

**Solution**: Process pathways in batches:
```python
for pathway_batch in chunks(pathways, 10):
    process_batch(pathway_batch)
    gc.collect()
```

### Debug Mode

Enable debug logging:

```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

### Cache Reset

Clear cached data if experiencing stale data issues:

```python
from generate_db import clear_cache
clear_cache()
```

## API Reference

### Main Functions

#### `run_pipeline(pathway_file, output_dir, model_name, **kwargs)`
Run the complete database generation pipeline.

**Parameters:**
- `pathway_file` (str): Path to file containing KEGG pathway IDs
- `output_dir` (str): Output directory for generated files
- `model_name` (str): Base name for output files
- `**kwargs`: Additional configuration options

**Returns:**
- `cobra.Model`: Generated metabolic model

#### `rxnsFromPathway(pathway_id)`
Extract reactions from a KEGG pathway.

**Parameters:**
- `pathway_id` (str): KEGG pathway ID (e.g., 'hsa00010')

**Returns:**
- `list`: List of reaction IDs

#### `getCompParamFromRestAPI(compound_id)`
Fetch compound parameters from KEGG REST API.

**Parameters:**
- `compound_id` (str): KEGG compound ID (e.g., 'C00001')

**Returns:**
- `dict`: Compound parameters including NAME, FORMULA, etc.

### Classes

#### `KEGGPWYS`
Class for managing KEGG pathway data.

```python
class KEGGPWYS:
    def __init__(self, rxnsFile=None):
        """Initialize with optional cached reactions file."""
        
    def get_reactions(self, pathway_id):
        """Get reactions for a pathway."""
        
    def get_all_compounds(self):
        """Get all unique compounds across pathways."""
```

#### `ModelBuilder`
Class for assembling COBRA models.

```python
class ModelBuilder:
    def __init__(self, model_id, model_name):
        """Initialize model builder."""
        
    def add_metabolite(self, met_id, **kwargs):
        """Add metabolite to model."""
        
    def add_reaction(self, rxn_id, **kwargs):
        """Add reaction to model."""
        
    def export_sbml(self, filepath):
        """Export model to SBML."""
```

## Contributing

### Development Setup

1. Fork the repository
2. Create a feature branch:
   ```bash
   git checkout -b feature/my-feature
   ```
3. Install development dependencies:
   ```bash
   pip install -r requirements-dev.txt
   ```
4. Run tests:
   ```bash
   pytest tests/
   ```

### Code Style

- Follow PEP 8 guidelines
- Use type hints for function signatures
- Document all public functions with docstrings
- Add tests for new functionality

### Submitting Changes

1. Ensure all tests pass
2. Update documentation if needed
3. Create a pull request with a clear description

## Citation

If you use this tool in your research, please cite:

```bibtex
@software{thg_human_gem_protocol,
  author = {Marinde Mas Lab},
  title = {THG: The Human GEM Protocol},
  year = {2024},
  url = {https://github.com/MarindeMasLab/THG_The-Human-GEM-protocol}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](../LICENSE) file for details.

## Acknowledgments

- KEGG database for metabolic pathway data
- The COBRApy team for the metabolic modeling framework
- Human Metabolic Atlas for reference annotations

---

**Last Updated**: December 2024  
**Version**: 2.0.0  
**Maintainer**: Marinde Mas Lab
