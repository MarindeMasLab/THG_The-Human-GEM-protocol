# Gap-Filling Pipeline for Metabolic Network Reconstruction

This module provides a comprehensive 3-phase gap-filling pipeline designed to connect isolated components in genome-scale metabolic models by adding **Putative Transport Reactions (PTRs)**. The pipeline minimizes the number of added reactions while maximizing network connectivity and reducing blocked reactions.

## Table of Contents

- [Overview](#overview)
- [Pipeline Architecture](#pipeline-architecture)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Detailed Usage](#detailed-usage)
  - [Phase 1: Candidate Generation](#phase-1-candidate-generation)
  - [Phase 2: Minimal Connector Selection](#phase-2-minimal-connector-selection)
  - [Phase 3: Blocked Reaction Optimization](#phase-3-blocked-reaction-optimization)
- [Phase 3 Strategies](#phase-3-strategies)
- [Output Files](#output-files)
- [CLI Reference](#cli-reference)
- [Examples](#examples)
- [File Structure](#file-structure)
- [Algorithm Details](#algorithm-details)

---

## Overview

Genome-scale metabolic models often contain **isolated components**—groups of reactions and metabolites that are not connected to the main metabolic network. These isolated components result from:

- Missing transport reactions between compartments
- Incomplete annotation of metabolic pathways
- Database gaps in reaction-metabolite associations

This pipeline identifies and adds **Putative Transport Reactions (PTRs)** to connect these isolated components while:

1. **Minimizing** the number of added reactions (parsimony)
2. **Maximizing** the reduction in blocked reactions (functionality)
3. **Prioritizing** high-confidence candidates (Type A > Type B > Type C)

---

## Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         GAP-FILLING PIPELINE                            │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  PHASE 1: Candidate Generation                                          │
│  ─────────────────────────────                                          │
│  • Identify dead-end metabolites (substrates/products)                  │
│  • Generate PTR candidates connecting same base metabolite              │
│  • Classify candidates: Type A, B, C                                    │
│  • Output: candidates_all.csv, deadends_summary.csv                     │
│                                                                         │
│                              ↓                                          │
│                                                                         │
│  PHASE 2: Minimal Connector Selection                                   │
│  ────────────────────────────────────                                   │
│  • Build component graph from candidates                                │
│  • Select minimal spanning tree (MST) to connect all components         │
│  • Prioritize Type A > Type B > Type C                                  │
│  • Output: selected_minimal_connectors.csv, phase2_model.json           │
│                                                                         │
│                              ↓                                          │
│                                                                         │
│  PHASE 3: Blocked Reaction Optimization                                 │
│  ──────────────────────────────────────                                 │
│  • Identify blocked reactions in each original component                │
│  • Select additional PTRs to unblock reactions                          │
│  • Multiple strategies: greedy, MILP, component-wise, sink_milp         │
│  • Output: phase3_selected_ptrs.csv, phase3_results.json                │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Installation

### Requirements

```bash
pip install -r requirements.txt
```

### Dependencies

- **COBRApy** (≥0.26.0): Constraint-based metabolic modeling
- **NetworkX** (≥2.6): Graph algorithms for component analysis
- **NumPy** (≥1.20): Numerical operations
- **SciPy** (≥1.7): MILP optimization (scipy.optimize.milp)
- **PuLP** (≥2.7): Alternative MILP solver interface (optional)

---

## Quick Start

### Run the full pipeline with recommended settings:

```bash
# Deadends objective (minimize dead-end metabolites)
python3 gapfill/gapfill.py run-all --pipeline deadends

# Blocked reactions objective with sink_milp strategy (RECOMMENDED)
python3 gapfill/gapfill.py run-all --pipeline blocked --strategy sink_milp \
    --parallel-components --workers-components 4 \
    --parallel-fba --workers-fba 2
```

### Run individual phases:

```bash
# Phase 1 only
python3 gapfill/gapfill.py phase1

# Phase 2 only
python3 gapfill/gapfill.py phase2 --mode minimal

# Phase 3 only
python3 gapfill/gapfill.py phase3 --mode blocked
```

---

## Detailed Usage

### Phase 1: Candidate Generation

Identifies dead-end metabolites and generates PTR candidates.

```bash
python3 gapfill/gapfill.py phase1 --model path/to/model.json --out output_dir/
```

**What it does:**
1. Loads the unconnected metabolic model
2. Identifies **dead-end metabolites**:
   - **Substrates**: Consumed but never produced
   - **Products**: Produced but never consumed
3. Groups metabolites by base ID (e.g., `MAM00001c`, `MAM00001m` → base `MAM00001`)
4. Generates **PTR candidates** connecting same base metabolites across compartments
5. Classifies candidates:
   - **Type A**: Both endpoints are dead-ends (highest confidence)
   - **Type B**: One endpoint is a dead-end
   - **Type C**: Neither endpoint is a dead-end (lowest confidence)

**Outputs:**
- `candidates_all.csv`: All PTR candidates with classification
- `deadends_summary.csv`: Dead-end metabolites by compartment
- `components_summary.csv`: Connected component sizes

---

### Phase 2: Minimal Connector Selection

Selects a minimal set of PTRs to connect all isolated components.

```bash
python3 gapfill/gapfill.py phase2 --mode minimal --candidates candidates_all.csv
```

**Modes:**
- `minimal`: MST-based minimal connector selection (default)
- `prioritized`: Priority-based selection (Type A → B → C)

**What it does:**
1. Builds a graph where nodes are components and edges are PTR candidates
2. Computes a **Minimum Spanning Tree (MST)** to connect all components
3. Selects one representative PTR per component pair (preferring Type A)
4. Adds selected PTRs to the model

**Outputs:**
- `selected_minimal_connectors.csv`: Selected PTRs for connectivity
- `THG-beta-batch_251106_phase2_minimal_connected.json`: Connected model

---

### Phase 3: Blocked Reaction Optimization

Selects additional PTRs to minimize blocked reactions.

```bash
python3 gapfill/gapfill.py phase3 --mode blocked --strategy sink_milp
```

**What it does:**
1. Identifies **blocked reactions** (reactions that cannot carry flux)
2. Tests each PTR candidate to see which blocked reactions it can unblock
3. Uses optimization (MILP/greedy) to select minimal PTRs with maximum coverage
4. Preserves **original component assignments** from Phase 1

---

## Phase 3 Strategies

| Strategy | Description | Speed | Quality | Recommended |
|----------|-------------|-------|---------|-------------|
| `exact` | Greedy per-candidate | Slow | Good | No |
| `hybrid` | Greedy with top-K prefiltering | Medium | Good | No |
| `hybrid_batch` | Batch greedy with checkpointing | Medium | Good | No |
| `milp` | Global MILP optimization | Very Slow | Optimal | No |
| `component_milp` | MILP per component | Fast | Good | Yes |
| `tiered_milp` | Tiered MILP (A→B→C) | Fast | Good | Yes |
| `sink_milp` | Temp sinks + hybrid parallelization | Fast | Best | **Yes** |

### Recommended: `sink_milp`

The `sink_milp` strategy is the most advanced and recommended approach:

```bash
python3 gapfill/gapfill.py run-all --pipeline blocked --strategy sink_milp \
    --components 4,5,6,7 \
    --parallel-components --workers-components 4 \
    --parallel-fba --workers-fba 2 \
    --solver-lp glpk
```

**Key features:**
- **Temporary sink/source approach**: Adds temporary exchange reactions for dead-end metabolites during FBA testing
- **Original component preservation**: Tracks which reactions belonged to which component before Phase 2 connected them
- **Hybrid parallelization**:
  - Small components (< 500 nodes): Processed in parallel
  - Large components: Processed sequentially with PTR reuse
- **Intra-component FBA parallelization**: Multiple FBA tests run in parallel within each component

---

## Output Files

### Phase 1 Outputs

| File | Description |
|------|-------------|
| `candidates_all.csv` | All PTR candidates with type classification |
| `deadends_summary.csv` | Dead-end metabolites by type and compartment |
| `components_summary.csv` | Component sizes and statistics |

### Phase 2 Outputs

| File | Description |
|------|-------------|
| `selected_minimal_connectors.csv` | Selected PTRs for connectivity |
| `*_phase2_minimal_connected.json` | Model with Phase 2 PTRs added |

### Phase 3 Outputs

| File | Description |
|------|-------------|
| `phase3_selected_ptrs.csv` | Final selected PTRs |
| `phase3_results.json` | Detailed results per component |

---

## CLI Reference

### Main CLI

```bash
python3 gapfill/gapfill.py <command> [options]
```

### Commands

| Command | Description |
|---------|-------------|
| `phase1` | Run Phase 1 candidate generation |
| `phase2` | Run Phase 2 connector selection |
| `phase3` | Run Phase 3 optimization |
| `run-all` | Run complete pipeline (Phase 1 → 2 → 3) |

### Global Options

| Option | Description | Default |
|--------|-------------|---------|
| `--model` | Input model JSON | Auto-detected |
| `--candidates` | Phase 1 candidates CSV | Auto-detected |
| `--out` | Output directory | `gapfill/files/` |

### Phase 3 Options (sink_milp)

| Option | Description | Default |
|--------|-------------|---------|
| `--strategy` | Optimization strategy | `exact` |
| `--components` | Comma-separated component IDs | All |
| `--parallel-components` | Enable parallel component processing | False |
| `--workers-components` | Workers for component parallelization | 4 |
| `--parallel-fba` | Enable parallel FBA testing | False |
| `--workers-fba` | Workers for FBA parallelization | 2 |
| `--small-threshold` | Size cutoff for parallel vs sequential | 500 |
| `--solver-lp` | LP solver (glpk, glpk_exact, scipy) | glpk |
| `--lambda` | MILP tradeoff weight | 0.01 |
| `--min-comp-size` | Minimum component size to process | 4 |
| `--max-components` | Maximum components to process | None |

---

## Examples

### Example 1: Full pipeline with default settings

```bash
python3 gapfill/gapfill.py run-all --pipeline blocked --strategy sink_milp
```

### Example 2: Process specific components with parallelization

```bash
python3 gapfill/gapfill.py run-all --pipeline blocked --strategy sink_milp \
    --components 4,5,6,7 \
    --parallel-components --workers-components 4 \
    --parallel-fba --workers-fba 2
```

### Example 3: Run Phase 3 standalone

```bash
python3 gapfill/phase3_sink_milp_original.py \
    --components 4,5,6,7 \
    --parallel-components --workers-components 2 \
    --parallel-fba --workers-fba 2 \
    --solver-lp glpk
```

### Example 4: Deadends objective

```bash
python3 gapfill/gapfill.py run-all --pipeline deadends --max 500
```

---

## File Structure

```
gapfill/
├── README.md                        # This file
├── requirements.txt                 # Python dependencies
├── gapfill.py                       # Main CLI orchestrator
│
├── phase1_connect_components.py     # Phase 1: Candidate generation
├── phase2_minimal_connector.py      # Phase 2: Minimal connector (MST)
├── phase2_prioritized_connector.py  # Phase 2: Prioritized connector
│
├── phase3_greedy_optimizer.py       # Phase 3: Greedy (deadends)
├── phase3_blocked_optimizer.py      # Phase 3: Greedy (blocked)
├── phase3_milp_optimizer.py         # Phase 3: Global MILP
├── phase3_component_milp.py         # Phase 3: Per-component MILP
├── phase3_tiered_milp.py            # Phase 3: Tiered MILP (A→B→C)
├── phase3_sink_milp_original.py     # Phase 3: Sink MILP (RECOMMENDED)
│
└── files/                           # Output directory
    ├── candidates_all.csv
    ├── components_summary.csv
    ├── deadends_summary.csv
    ├── selected_minimal_connectors.csv
    ├── phase3_selected_ptrs.csv
    └── phase3_results.json
```

---

## Algorithm Details

### PTR Candidate Classification

| Type | Definition | Confidence |
|------|------------|------------|
| **A** | Both metabolites are dead-ends | Highest |
| **B** | One metabolite is a dead-end | Medium |
| **C** | Neither metabolite is a dead-end | Lowest |

### Sink MILP Approach

The `sink_milp` strategy uses **temporary sink/source reactions** to test if a PTR can unblock reactions:

1. **Dead-end substrates** (consumed but not produced): Add temporary **source** reaction (→ met)
2. **Dead-end products** (produced but not consumed): Add temporary **sink** reaction (met →)

This allows FBA to find feasible flux distributions even when metabolites have nowhere to go, enabling accurate testing of PTR candidates.

### Hybrid Parallelization

The pipeline uses a two-level parallelization strategy:

1. **Inter-component parallelization**: Multiple small components processed simultaneously
2. **Intra-component parallelization**: Multiple FBA tests run in parallel within each component

Memory safety limits are enforced to prevent system crashes.

---

## Authors

- THG Team

## License

See the main repository LICENSE file.
