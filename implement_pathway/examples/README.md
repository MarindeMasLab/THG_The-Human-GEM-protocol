# Examples Guide

This directory contains example scripts demonstrating how to use the pathway implementation pipeline.

## 📁 Directory Structure

```
examples/
├── README.md                    # This file
│
├── inputs/                      # Input files (START HERE)
│   ├── endoA_250917_3.json     # Base model (~16 MB)
│   └── config_example.json     # Example configuration
│
├── outputs/                     # Generated outputs
│   ├── data/                   # Metabolite databases
│   ├── reports/                # Implementation/validation reports
│   ├── figures/                # Visualizations (D3.js, Plotly, GraphViz)
│   └── *.json                  # Modified models
│
└── run_example.py              # Complete workflow (RECOMMENDED)
```

**Note**: All files in `outputs/` are git-ignored.

---

## 🚀 Quick Start

### Step 1: Run the Example

```bash
cd implement_pathway
python3 examples/run_example.py
```

This will:
1. Load base model from `examples/inputs/endoA_250917_3.json`
2. Load configuration from `examples/inputs/config_example.json`
3. Build metabolite ID database
4. Add Example_Pathway with 3 reactions:
   - ATP transport (extracellular → example compartment)
   - Protein synthesis (ATP → Protein + ADP + Pi)
   - Protein demand (Protein →)
5. Save outputs to `examples/outputs/`

**Expected outputs:**
```
examples/outputs/
├── data/
│   └── metabolite_id_database.json
├── reports/
│   └── 2025-11-10_example_implementation_report.txt
└── endoA_250917_3_with_Example_Pathway.json
```

### Step 2: Validate the Implementation

```bash
python3 validate.py 0 \
  --config examples/inputs/config_example.json \
  --model examples/outputs/endoA_250917_3_with_Example_Pathway.json
```

This generates:
- Mass balance report
- FBA verification
- Reaction flux analysis

Reports saved to `examples/outputs/reports/`.

### Step 3: Visualize the Network

```bash
python3 visualize.py 0 \
  --config examples/inputs/config_example.json \
  --model examples/outputs/endoA_250917_3_with_Example_Pathway.json
```

This generates three types of visualizations in `examples/outputs/figures/`:
- `Example_Pathway_network_d3.html` - Interactive D3.js network
- `Example_Pathway_network_plotly.html` - Interactive Plotly with FVA
- `Example_Pathway_network.png/svg/pdf` - Static GraphViz images

---

## 📖 Understanding the Example Configuration

The `config_example.json` defines a minimal pathway with 3 reactions:

### Reaction 1: ATP Transport
```json
{
  "reaction_abbreviation": "ATP_exp_transport",
  "reaction_name": "ATP transport to example compartment",
  "reaction_equation": "atp[e] <=> atp[exp]",
  "reversibility": true,
  "category": "Transport"
}
```

### Reaction 2: Protein Synthesis
```json
{
  "reaction_abbreviation": "PROT_exp_synthesis", 
  "reaction_name": "Generic protein synthesis in example compartment",
  "reaction_equation": "atp[exp] => prot_generic[exp] + adp[exp] + pi[exp]",
  "reversibility": false,
  "category": "Biosynthesis"
}
```

### Reaction 3: Protein Demand
```json
{
  "reaction_abbreviation": "DM_prot_generic_exp",
  "reaction_name": "Demand reaction for generic protein",
  "reaction_equation": "prot_generic[exp] =>",
  "reversibility": false,
  "category": "Demand"
}
```

This minimal example demonstrates:
- ✅ Transport reactions (reversible)
- ✅ Biosynthesis reactions (irreversible)
- ✅ Demand reactions (sink)
- ✅ Stoichiometry (1:1 transport, 1:3 synthesis)
- ✅ Compartmentalization

---

## 🔧 Customizing the Example

### Modify the Configuration

Edit `examples/inputs/config_example.json` to:
- Add more reactions
- Change metabolites
- Adjust stoichiometry
- Set flux bounds

```json
{
  "pathway_name": "My_Custom_Pathway",
  "compartment_abbreviation": "mycomp",
  "compartment_name": "My Custom Compartment",
  "reactions": [
    // Your custom reactions here
  ]
}
```

### Use a Different Model

Replace the model file in `inputs/`:
```bash
cd examples/inputs
cp /path/to/YOUR_MODEL.json ./
# Update config_example.json to reference YOUR_MODEL.json
```

### Change Output Locations

Edit paths in `run_example.py`:
```python
# Lines 23-26
SCRIPT_DIR = Path(__file__).parent
INPUTS_DIR = SCRIPT_DIR / "inputs"
OUTPUTS_DIR = SCRIPT_DIR / "outputs"
DATA_DIR = OUTPUTS_DIR / "data"
```

---

## 📊 Example Output Statistics

### Model Before Implementation
- **Metabolites**: 16,670
- **Reactions**: 19,300
- **Compartments**: 10 (e, c, m, n, x, r, l, v, g, p)
- **Size**: ~16 MB

### Model After Implementation
- **Metabolites**: 16,673 (+3)
- **Reactions**: 19,303 (+3)
- **Compartments**: 11 (+1: exp)
- **Size**: ~16 MB

### Added Components
- **New compartment**: `exp` (Example_Pathway)
- **New metabolites**:
  - `atp_exp` - ATP in example compartment
  - `prot_generic_exp` - Generic protein
  - `adp_exp` - ADP in example compartment
  - `pi_exp` - Inorganic phosphate
- **New reactions**:
  - `ATP_exp_transport` - Transport
  - `PROT_exp_synthesis` - Biosynthesis
  - `DM_prot_generic_exp` - Demand

---

## 🧪 Testing Different Scenarios

### Test 1: Flux Balance Analysis
```python
import cobra
model = cobra.io.load_json_model('examples/outputs/endoA_250917_3_with_Example_Pathway.json')

# Optimize
solution = model.optimize()
print(f"Objective value: {solution.objective_value}")

# Check protein synthesis flux
print(f"Protein synthesis flux: {solution.fluxes['PROT_exp_synthesis']}")
```

### Test 2: Flux Variability Analysis
```python
from cobra.flux_analysis import flux_variability_analysis

# FVA for example reactions
rxn_ids = ['ATP_exp_transport', 'PROT_exp_synthesis', 'DM_prot_generic_exp']
fva_result = flux_variability_analysis(model, rxn_ids)
print(fva_result)
```

### Test 3: Mass Balance Check
```bash
python3 validate.py 0 \
  --config examples/inputs/config_example.json \
  --model examples/outputs/endoA_250917_3_with_Example_Pathway.json
```

---

## 🆘 Troubleshooting

### Issue: "Database not found"
**Cause**: `metabolite_id_database.json` hasn't been generated yet

**Solution**: Run `run_example.py` first - it builds the database automatically

### Issue: "Model file not found"
**Cause**: Model file is missing from `inputs/`

**Solution**: Copy the model file:
```bash
cd examples/inputs
cp ../../../models/endoA_250917_3.json ./
```

### Issue: "Permission denied" in outputs/
**Cause**: Insufficient write permissions

**Solution**: 
```bash
chmod -R u+w examples/outputs/
```

### Issue: "Module not found" errors
**Cause**: Missing dependencies

**Solution**:
```bash
pip install -r requirements.txt
```

### Issue: Visualizations not generated
**Cause**: Missing visualization dependencies

**Solution**:
```bash
pip install plotly networkx graphviz
# For GraphViz system dependency (Ubuntu/Debian):
sudo apt-get install graphviz
```

---

## 🔗 Next Steps

After running the example successfully:

1. **Read Main Documentation**:
   - `../README.md` - Complete project documentation
   - `../INTEGRATION_GUIDE.md` - Integration instructions

2. **Try Real Pathways**:
   - Use `config/config_cytoskeleton.json` (46 reactions)
   - Use `config/config_glycocalyx.json` (glycocalyx pathway)
   - Create your own configuration

3. **Advanced Topics**:
   - Customize reaction categories for visualization
   - Add pathway-specific validation rules
   - Integrate with COBRApy analysis tools
   - Use the pipeline programmatically in your own code

---

## ✨ Benefits of the Examples Structure

1. **Isolated I/O**: All example inputs and outputs are contained in `examples/`
2. **No Interference**: Example runs don't affect main pipeline files
3. **Easy Cleanup**: Delete `outputs/` to reset without affecting main work
4. **Clear Learning Path**: Single recommended starting point (`run_example.py`)
5. **Minimal Example**: Only 3 reactions - easy to understand
6. **Complete Workflow**: Shows database building, implementation, and reporting
7. **Git-Friendly**: All outputs are automatically ignored by git

---

**Last Updated**: November 2025

