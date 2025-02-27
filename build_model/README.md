# Build model

The build model repository is a python module to enrich gene annotation, curate and expands a metabolic network model based on information currently available in multiple online databases (KEGG, CheBI, PubChem, LIPIDMAPS, Ensembl...)

build_model.py:
- Mass balance metabolic reactions (equations_build_model.mass_balance)
- Construct SGPRs and GPRs (gpr.auth_gpr.getGPR)
- Identifies the cellular location of the metabolic reactions (gpr.getLocation.getLocation)
- Additional functions necessary to run the script (function_build_model.py)
## Running build model

To run the model building:

```
python3 build_model.py
```
