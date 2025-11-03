# Metabolite reac identification

The "metabolite_reac_identification" repository is a python module to enrich metabolite and reaction annotation of a reference human genome-scale metabolic model

## Overview

metabolite_reac_identification.py:
- Identifies metabolites by name and molecular formula (function_metabolite_identification.generate_met_annotation)
- Identifies reactions by unique combination of substrates and products (function_reac_identification.execute_jaccard)
- Correct/Enrich reference model with the information collected in the previous steps
- Writes the corrected/enriched new model (function_annotate_cobra_model.annotate_cobra_model)

## Running metabolite reac identification

### Standard Version (Single Pass)

To run the basic metabolite reac identification:

```bash
python3 metabolite_reac_identification.py
```

### Intelligent Multi-Pass Version (Recommended)

For better results with automatic retry of failures:

```bash
/path/to/venv/bin/python metabolite_reac_identification_with_retry.py
```

Or using the virtual environment:

```bash
/home/igor/Documents/THG_The-Human-GEM-protocol/.venv/bin/python metabolite_reac_identification_with_retry.py
```

**Run in background:**

```bash
nohup /home/igor/Documents/THG_The-Human-GEM-protocol/.venv/bin/python metabolite_reac_identification_with_retry.py > metabolite_run.log 2>&1 &
```

## Intelligent Multi-Pass System

The improved script (`metabolite_reac_identification_with_retry.py`) automatically retries failed metabolites intelligently:

### How It Works

1. **First Pass**: Attempts to annotate all metabolites from the model

2. **Smart Retry System**: 
   - Automatically retries ALL failures from the first pass
   - Continues retrying until one of the following conditions is met:
     - **Threshold reached**: Remaining failures < 10% of potential achievable annotations
     - **No improvement**: No new successful annotations in a round
     - **Max rounds reached**: Completes maximum number of retry rounds (default: 10)

3. **Intelligent Stopping Logic**:
   - Distinguishes between temporary API errors (worth retrying) and metabolites genuinely not in PubChem database
   - Maximizes annotation success rate while avoiding infinite retry loops

### Example

Given a model with 1,000 metabolites:
- After Round 1: 450 identified, 550 failed
- After Round 2: 50 more identified → 500 total identified, 500 still failed
- If 250 are genuinely not in PubChem, potential max = 750
- Threshold = 10% of 750 = 75
- If remaining failures ≤ 75 → Stop (these are likely not in PubChem)
- Otherwise → Continue retrying

### Configuration

Easily adjustable parameters at the top of `metabolite_reac_identification_with_retry.py`:

```python
MAX_RETRY_ROUNDS = 10    # Maximum number of retry attempts
STOP_THRESHOLD = 0.10    # Stop if remaining failures < 10% of achievable
```

### Benefits

✅ **Fully automated** - No manual intervention needed  
✅ **Maximizes success rate** - Keeps trying while there's improvement  
✅ **Efficient** - Stops when hitting diminishing returns  
✅ **Intelligent** - Distinguishes between API errors and missing data  
✅ **Complete pipeline** - After annotation, continues with reaction identification and model enrichment

### Monitoring Progress

Use the monitoring script to check progress at any time:

```bash
./monitor_progress.sh
```

This displays:
- Current progress (metabolites annotated / total)
- Success rate
- Recent successful annotations
- Error analysis with failure counts by type
- Recent log activity

### Output Files

- `reports/met_annotation.tsv` - Successfully annotated metabolites
- `reports/met_annotation_failures.tsv` - Failed metabolites with failure reasons
- `metabolite_run.log` - Detailed execution log

### Expected Runtime

- First pass: ~8 hours (for ~4,000 metabolites)
- Retry rounds: ~2-4 hours per round (depends on failure count)
- Total: ~10-15 hours for complete process

**Recommended**: Run overnight and check results in the morning

