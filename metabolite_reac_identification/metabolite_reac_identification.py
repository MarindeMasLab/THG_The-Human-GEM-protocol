import cobra.io
import numpy as np
import re
import os
import sys

# Determine the current file's directory and the project root.
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.join(current_dir, "..")

# Add the project root to sys.path to access top-level folders like 'functions' and 'models'
if project_root not in sys.path:
    sys.path.append(project_root)

from functions.function_metabolite_identification import *
from functions.function_reac_identification import *
from functions.function_annotate_cobra_model import *

model = os.path.join(project_root, "models", "Human-GEM1_19.xml")

print(f"Number of reactions in model: {len(cobra.io.read_sbml_model(model).reactions)}")

database = os.path.join(project_root, "models", "Human_Database_old.xml")

print(
    f"Number of reactions in database: {len(cobra.io.read_sbml_model(database).reactions)}"
)

cobra_model = cobra.io.read_sbml_model(model)


# Metabolites

metabolites = gather_metabolites(cobra_model)
# Use 0.25s delay to ensure no more than 5 requests per second
# (accounting for up to 3 API calls per metabolite: get_cids, alternative name, from_cid)
generate_met_annotation(metabolites, delay_between_requests=0.25)
met_annotation = process_annotation()


# Reactions

reac = process_reac(model, "[A-Z]+[0-9]+[a-z]+[0-9]*", "MAM02040", "MAM02039")
reac = replace_met_id_by_met_kegg(reac, gather_kegg_metabolites(model))
reac_y = process_reac(database, "([A-Z][0-9]+_?[a-z]+[0-9]*)", "C00080", "C00001")
jaccard = execute_jaccard(reac_y, reac)
reac_annotation = process_jaccard(reac_y, reac, jaccard)


# SBML
logging.info("Annotating the model and writing SBML files...")

annotate_cobra_model(cobra_model, met_annotation, reac_annotation)
