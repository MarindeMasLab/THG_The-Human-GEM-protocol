# Generate databae

The generate database directory is a python module to generates a metabolic network by gathering currently available information from a large number of online databases.

generate_database.py:

## Running generate database

To run the database generation:

```
```
```markdown
# generate_data-base

This folder contains the scripts and small data files used to build the Human
metabolic database (SBML) from KEGG (and related sources). The main entry
point is `generate_db.py` which performs a multi-stage pipeline that:

- downloads pathway, reaction and compound data from KEGG REST endpoints;
- constructs reaction stoichiometry and attempts to mass-balance reactions;
- fetches and sanitizes GPR and subcellular localization (SGPR) data;
- compartmentalizes reactions and metabolites into cellular compartments;
- annotates metabolites and genes with external database identifiers;
- assembles a COBRA model and writes it as SBML (`models/Human_Database.xml`).

Contents
- `generate_db.py` — Main pipeline. Run this to regenerate the complete
	database from KEGG.
- `requirements.txt` — Python dependencies required by the scripts (keep this).
- `make_model_from_pkl.py`, `resume_db_gen.py` — helpers to rebuild a model from
	intermediate `*.pkl` files produced by `generate_db.py`.
- `gly2` — generated text file containing glycan formulas that include grouped
	repeating units (created by the pipeline; not required to start a clean run).

What was cleaned
- Runtime log files and compiled caches were removed from this folder to
	provide a clean state for fresh database generation runs. The main scripts
	and `requirements.txt` were preserved.

Prerequisites
- Python 3.10+ (this repository imports `cobra`, `dill`, `tqdm`, `requests`,
	etc.).
- Install the dependencies listed in this folder:

```bash
python3 -m pip install -r generate_data-base/requirements.txt
```

- Ensure outgoing network access to KEGG and Ensembl REST APIs (the script
	fetches data online). If you run behind a proxy set `HTTP_PROXY`/`HTTPS_PROXY`.

Running a clean build (from scratch)

1. Remove any checkpoint or intermediate pickle files to force a fresh run:

```bash
rm -f files/checkpoint_progress.pkl files/pre_sbml_*.pk
```

2. Run the generator from the repository root (so `functions/` is importable):

```bash
python3 generate_data-base/generate_db.py
```

Notes and tips
- The script saves periodic checkpoints to `files/checkpoint_progress.pkl` so
	that long runs can be resumed. Keep the checkpoint if you want to resume a
	partially completed run.
- For development and faster iteration set `MAX_REACTIONS` near the top of
	`generate_db.py` to a small integer (e.g. `100`) to limit the number of
	processed reactions.
- If you only have a finished pickle (`pre_sbml_pos_comp.pk` or
	`pre_sbml_raw.pk`) use `make_model_from_pkl.py` or `resume_db_gen.py` to
	reconstruct an SBML model without re-running the full KEGG/Ensembl fetches.
- `gly2` is an auxiliary file produced when glycan formulas include grouped
	repeating units; it will be recreated by the pipeline when needed and does
	not need to be present for a clean run.

Further actions I can do next
- remove `gly2` if you prefer no generated artifacts checked in (it will be
	recreated by the script when necessary);
- add a small `clean.sh` to remove checkpoints and generated files before a
	fresh run; or
- run a short dry-run (set `MAX_REACTIONS=50`) and report missing dependencies
	and external API errors.

```
