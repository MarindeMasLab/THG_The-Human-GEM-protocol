# Generate Human Metabolic Database

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![COBRApy](https://img.shields.io/badge/COBRApy-0.29+-green.svg)](https://opencobra.github.io/cobrapy/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](../LICENSE)

## Overview

This module provides a comprehensive pipeline for generating genome-scale metabolic models (GEMs) from KEGG pathway data, with optional extension via **Rhea** reactions. It automatically fetches reactions and metabolites from KEGG (and Rhea), resolves compound identities, calculates molecular formulas, performs mass balance validation, and exports SBML-compatible models for use with COBRApy and other metabolic modeling tools.

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
- [Key Functions and Classes](#key-functions-and-classes)
- [Contributing](#contributing)

## Features

### Core Capabilities

- **KEGG Integration**: Automated fetching of reactions, compounds, and pathways from KEGG REST API
- **Rhea Extension**: Optional integration of Rhea reactions with ChEBI metabolites, Reactome pathway grouping, and UniProt EC-to-gene mapping
- **Compound Identification**: Multi-source metabolite identification using PubChem, KEGG, ChEBI, and custom databases
- **Formula Resolution**: Intelligent formula parsing including glycan composition calculation
- **Mass Balance Validation**: Automatic detection and reporting of mass-imbalanced reactions
- **SBML Export**: Export to SBML format compatible with COBRApy, Escher, and other tools
- **Gene Associations**: Integration of gene-protein-reaction (GPR) rules
- **Compartmentalization**: Support for multiple cellular compartments

### Recent Improvements (v2.1)

1. **Rhea Extension Module**: New multi-database pipeline supporting:
   - Rhea reaction database integration via SPARQL queries
   - ChEBI metabolite resolution with formula and charge
   - Reactome pathway grouping with human-readable names
   - UniProt EC-to-gene mapping for GPR rules
   - Ensembl gene annotation integration
2. **Fixed Mass Balance Counting**: Corrected false positive reporting in mass balance validation
3. **KEGG URL Parsing Update**: Adapted to new KEGG stoichiometry URL format (`map` → `R`)
4. **Glycan Formula Support**: Added fallback mechanism for glycan formulas and composition-based calculation
5. **Improved Group IDs**: Clean SBML SId-compatible group identifiers

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
|---------|---------|----------|
| cobra | ≥0.29.0 | Metabolic model manipulation and SBML I/O |
| pandas | ≥2.0.0 | Data manipulation and analysis |
| numpy | ≥1.24.0 | Numerical computations |
| scipy | ≥1.10.0 | Scientific computing |
| requests | ≥2.28.0 | HTTP requests to KEGG/Rhea/UniProt APIs |
| pubchempy | ≥1.0.4 | PubChem compound identification |
| openpyxl | ≥3.0.0 | Excel file reading (compartment definitions) |
| tqdm | ≥4.65.0 | Progress bars |
| dill | ≥0.3.6 | Extended pickling support |
| xlsxwriter | ≥3.0.0 | Excel report generation |
| sympy | ≥1.12 | Symbolic mathematics for formula parsing |
| python-dotenv | ≥1.0.0 | Environment variable management |

## Quick Start

### Command Line Usage

The pipeline is executed as a script from the command line:

```bash
# Navigate to the generate_data-base directory
cd generate_data-base

# Run without Rhea extension (KEGG only)
python generate_db.py

# Run with Rhea extension (all human Rhea reactions)
python generate_db.py --rhea

# Run with specific Rhea IDs only
python generate_db.py --rhea 10040,10112,10116
```

### Environment Setup

Before running, set up BioCyc credentials in a `.env` file at the project root:

```env
# BioCyc credentials required for GPR (Gene-Protein-Reaction) data extraction
# You need to create an account at BioCyc from a licensed institution
BIOCYC_EMAIL=your.email@institution.edu
BIOCYC_PASSWORD=your_password
```

### Pathway Subset

To process only specific pathways, set the `PATHWAY_SUBSET` environment variable:

```bash
export PATHWAY_SUBSET="../files/human_kegg_pathways_subset.txt"
python generate_db.py
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
│              STAGE 5: GPR & sGPR RULE GENERATION                │
│  Fetch gene associations → Build GPR rules → Generate sGPR      │
└─────────────────────────────────────────────────────────────────┘
                                 ↓
┌─────────────────────────────────────────────────────────────────┐
│         STAGE 5b: RHEA EXTENSION (Optional)                     │
│  Query Rhea SPARQL → ChEBI metabolites → Reactome pathways      │
└─────────────────────────────────────────────────────────────────┘
                                 ↓
┌─────────────────────────────────────────────────────────────────┐
│          STAGE 6: ISOFORM-BASED COMPARTMENTALIZATION            │
│  Expand reactions by isoforms → Assign compartments → Localize  │
└─────────────────────────────────────────────────────────────────┘
                                 ↓
┌─────────────────────────────────────────────────────────────────┐
│                    STAGE 7: MODEL ASSEMBLY                      │
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

#### Stage 5: GPR & sGPR Rule Generation
- Fetches gene associations from KEGG for each reaction
- Builds Gene-Protein-Reaction (GPR) rules with AND/OR logic
- Generates simplified GPR (sGPR) rules for downstream analysis
- Maps EC numbers to genes via KEGG orthology
- Handles multi-enzyme complexes and isozyme alternatives
- Integrates with Ensembl for gene annotation and validation

#### Stage 5b: Rhea Extension (Optional)
- Queries Rhea database via SPARQL for additional reactions
- Resolves ChEBI metabolites with formulas and charges
- Maps reactions to Reactome pathways with human-readable names
- Retrieves EC numbers and maps to human genes via UniProt
- Integrates with existing KEGG-based reactions and pathways
- Configurable via `ENABLE_RHEA` set of Rhea master IDs

#### Stage 6: Isoform-Based Compartmentalization
- Expands reactions based on enzyme isoform localization
- Queries BioCyc for protein subcellular localization
- Creates compartment-specific copies of reactions
- Uses EndoA adjacency matrix to determine allowed transport pairs
- Transport reactions only created between adjacent compartments:
  - Cytosol ↔ Mitochondria (c ↔ m)
  - Mitochondria ↔ Inner Mitochondria (m ↔ i)
  - Cytosol ↔ ER, Golgi, Nucleus, etc.
- Handles isoform-specific metabolite pools
- Ensures metabolic connectivity across cellular compartments

#### Stage 7: Model Assembly
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
| `rhea_extension.py` | Rhea database integration - SPARQL queries, ChEBI metabolites, UniProt gene mapping |
| `transport_utils.py` | Transport reaction handling - compartment boundaries and EndoA adjacency matrix |
| `extend_model_with_rhea.py` | Standalone tool to extend existing SBML models with Rhea reactions |
| `make_model_from_pkl.py` | Utility to create COBRA model from pickle checkpoint files |
| `resume_db_gen.py` | Resume interrupted pipeline runs from checkpoint |

### Function Modules (in `../functions/`)

| File | Description |
|------|-------------|
| `class_generate_database.py` | Core classes: `reaction`, `compound`, `gene`, `gpr` |
| `function_bm_gdb.py` | KEGG API functions, formula parsing, compound parameter extraction |
| `equations_bm_gdb.py` | Reaction equation parsing and stoichiometry extraction |
| `functions_mass_balance.py` | Mass balance validation and atom counting |
| `pattern_generate_database.py` | Compartmentalization logic (`rxnSubcel`) |
| `ensembl_client.py` | Ensembl API client for gene annotations |
| `error_tracker.py` | Error logging and tracking utilities |

### Input Files (in `../files/`)

| File | Description |
|------|-------------|
| `human_kegg_pathways.txt` | List of KEGG pathway IDs to process |
| `ListOfCompartments_sept2024.xlsx` | Compartment definitions and EndoA adjacency matrix |
| `extra_compounds.txt` | Custom compound definitions |
| `extra_formula.txt` | Manual formula overrides |
| `special_compounds.txt` | Special handling rules for specific compounds |

## Configuration

### Environment Variables

Create a `.env` file in the project root:

```env
# BioCyc credentials required for GPR data extraction
# Requires an account from a licensed institution
BIOCYC_EMAIL=your.email@institution.edu
BIOCYC_PASSWORD=your_password
```

### Script Configuration (in `generate_db.py`)

Key configuration variables at the top of `generate_db.py`:

```python
# Debugging: limit number of reactions (None = no limit)
MAX_REACTIONS = None  # Set to integer for debugging

# Compartmentalization mode:
#   1 = RESTRICTED: Only use compartments from Excel lookup table
#   0 = UNRESTRICTED: Keep all compartment names as-is
IMPOSE_LOCATIONS = 1  # Default: restricted mode

# Rhea extension options:
#   ENABLE_RHEA = True        → Process ALL human Rhea reactions (~4700+)
#   ENABLE_RHEA = False       → Skip Rhea extension
#   ENABLE_RHEA = {'10040'}   → Process only specific Rhea IDs
ENABLE_RHEA = False  # Set to True for production
```

### Compartments

The pipeline uses 9 cellular compartments defined in `ListOfCompartments_sept2024.xlsx`:

| Abbreviation | Compartment |
|--------------|-------------|
| c | cytosol |
| m | mitochondria |
| i | inner mitochondria (intermembrane space) |
| n | nucleus |
| r | endoplasmic reticulum |
| g | golgi apparatus |
| l | lysosome |
| x | peroxisome |
| e | extracellular |

## Usage Examples

### Example 1: Full Pipeline (Command Line)

```bash
# Set up environment
cd generate_data-base
export BIOCYC_EMAIL="your.email@institution.edu"
export BIOCYC_PASSWORD="your_password"

# Run full pipeline with KEGG + Rhea
python generate_db.py --rhea
```

### Example 2: Process Specific Pathways

```bash
# Create a file with specific pathway IDs
echo -e "hsa00010\nhsa00020\nhsa00030" > ../files/my_pathways.txt

# Set environment variable and run
export PATHWAY_SUBSET="../files/my_pathways.txt"
python generate_db.py
```

### Example 3: Testing with Specific Rhea IDs

```bash
# Process only specific Rhea reactions for testing
python generate_db.py --rhea 10040,10112,10116,11436
```

### Example 4: Using the RheaExtender Programmatically

```python
from rhea_extension import RheaExtender

# Initialize extender with data structures from KEGG pipeline
extender = RheaExtender(
    MetList=MetList,
    MetIdent=MetIdent,
    MetEquiv=MetEquiv,
    RxnList=RxnList,
    RxnIdent=RxnIdent,
    GPRList=GPRList,
    GPRIdent=GPRIdent,
    GeneList=GeneList,
    GeneIdent=GeneIdent,
    PathNameRxn=PathNameRxn,
    impose_locations=1  # Use restricted compartmentalization
)

# Extend with specific Rhea IDs
extended_data = extender.extend_with_rhea(
    rhea_ids={'10040', '10112'},
    session=biocyc_session
)
```

## Output Files

### Generated Models

| File | Description |
|------|-------------|
| `models/Human_Database.xml` | Full SBML model (COBRApy compatible) |

### Error Reports

| File | Description |
|------|-------------|
| `logs/error_report.txt` | Detailed error report with all warnings and issues |
| `logs/error_report.json` | JSON format error data for programmatic analysis |
| `logs/generate_db.log` | Full pipeline execution log |

### Intermediate Files (in `files/`)

| File | Description |
|------|-------------|
| `pre_sbml_raw.pk` | Pickle checkpoint before compartmentalization |
| `pre_sbml_pos_comp.pk` | Pickle checkpoint after compartmentalization |
| `gene_location_cache.pkl` | Cached gene location data from BioCyc |

### Error Report Contents

The error report (`logs/error_report.txt` and `.json`) contains:

1. **Overall Statistics**:
   - Total errors by category
   - Success/failure rates
   - Processing summary

2. **Error Categories**:
   - API connection failures
   - Mass balance issues
   - Missing formula warnings
   - GPR retrieval errors
   - Compartmentalization issues

3. **Recommendations**:
   - Suggested fixes for common issues
   - Manual overrides needed

## Troubleshooting

### Common Issues

#### 1. KEGG API Connection Errors

```
Error: Failed to fetch data from KEGG REST API
```

**Solution**: Check internet connection. The pipeline includes automatic retries with delays.

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

**Solution**: Use the `MAX_REACTIONS` variable in `generate_db.py` to limit processing:
```python
MAX_REACTIONS = 100  # Process only first 100 reactions for testing
```

Or process a subset of pathways via environment variable:
```bash
export PATHWAY_SUBSET="../files/human_kegg_pathways_subset.txt"
```

### Debug Mode

Enable debug logging by modifying `generate_db.py`:

```python
logging.basicConfig(level=logging.DEBUG)  # Already default
# Change to logging.INFO for less verbose output
```

### Cache Reset

Clear cached data if experiencing stale data issues:

```bash
# Remove checkpoint and cache files
rm -f ../files/pre_sbml_raw.pk ../files/pre_sbml_pos_comp.pk
rm -f ../files/gene_location_cache.pkl
rm -f ../files/checkpoint_progress.pkl
```

## Key Functions and Classes

### `generate_db.py`

#### `cobra_reconstruction(ModName, ModID, MetList_CL, RxnList_CL, GeneList, PathNameRxn, LocVar, MetEquiv, MetList)`
Assemble the final COBRA model from processed data structures.

#### `process_reaction_gpr(rxn, gpr_list, gpr_ident, session, impose_locations)`
Process GPR (Gene-Protein-Reaction) associations for a reaction via BioCyc.

### `rhea_extension.py`

#### `RheaExtender`
Main class for extending KEGG-based models with Rhea reactions.

```python
extender = RheaExtender(
    MetList, MetIdent, MetEquiv, RxnList, RxnIdent,
    GPRList, GPRIdent, GeneList, GeneIdent, PathNameRxn,
    impose_locations=1
)
result = extender.extend_with_rhea(rhea_ids, session)
```

### `transport_utils.py`

#### `TransportBoundaryManager`
Manages compartment boundaries and transport reaction logic.

```python
manager = TransportBoundaryManager(excel_file)
allowed_pairs = manager.get_allowed_transport_pairs('mitochondria')
```

## Contributing

### Development Setup

1. Fork the repository
2. Create a feature branch:
   ```bash
   git checkout -b feature/my-feature
   ```
3. Install dependencies:
   ```bash
   pip install -r generate_data-base/requirements.txt
   ```

### Code Style

- Follow PEP 8 guidelines
- Use type hints for function signatures
- Document all public functions with docstrings

### Submitting Changes

1. Test your changes with a subset of pathways
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

**Last Updated**: December 2025  
**Version**: 2.1.0  
**Maintainer**: Marinde Mas Lab
