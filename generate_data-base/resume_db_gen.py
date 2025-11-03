import os
import sys
import dill
import logging
import importlib.util
import traceback


# Setup logging
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
console_handler.setFormatter(formatter)
if not LOGGER.handlers:
    LOGGER.addHandler(console_handler)


# Determine the current file's directory and the project root.
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.join(current_dir, "..")

# Ensure project root is on sys.path so pickled objects referencing the
# local `functions` package can be imported during unpickling.
if project_root not in sys.path:
    sys.path.insert(0, project_root)

# Path to the pre-serialized pickle created by generate_db.py
pickle_path = os.path.join(project_root, "files", "pre_sbml_pos_comp.pk")
if not os.path.exists(pickle_path):
    LOGGER.error("Pickle file not found: %s", pickle_path)
    raise SystemExit(1)

# Load the pickle
with open(pickle_path, "rb") as f:
    # This file was serialized with dill.dump earlier in generate_db.py
    data = dill.load(f)

LOGGER.debug(
    "Loaded pre_sbml_pos_comp.pk keys: %s",
    list(data.keys()) if isinstance(data, dict) else type(data),
)

# Map the pickled contents to the names expected by the reconstruction routine.
# The original generator stores these keys: name,id,mets,mets_cl,met_equiv,reactions,reactions_cl,genes,pathways,loc
ModName = data.get("name") if isinstance(data, dict) else None
ModID = data.get("id") if isinstance(data, dict) else None
MetList = data.get("mets") if isinstance(data, dict) else {}
MetList_CL = data.get("mets_cl") if isinstance(data, dict) else {}
MetEquiv = data.get("met_equiv") if isinstance(data, dict) else {}
RxnList = data.get("reactions") if isinstance(data, dict) else {}
RxnList_CL = data.get("reactions_cl") if isinstance(data, dict) else {}
GeneList = data.get("genes") if isinstance(data, dict) else {}
PathNameRxn = data.get("pathways") if isinstance(data, dict) else {}
LocVar = data.get("loc") if isinstance(data, dict) else {}

# Ensure there is an output path. The original flow used Output=ModID or a models path.
if isinstance(ModID, str) and ModID.lower().endswith(".xml"):
    Output = ModID
else:
    os.makedirs(os.path.join(project_root, "models"), exist_ok=True)
    Output = os.path.join(project_root, "models", "Recovered_model_from_pickle.xml")

# Import helper functions from generate_db.py using importlib so we don't execute its __main__ block
generate_db_path = os.path.join(project_root, "generate_data-base", "generate_db.py")
if not os.path.exists(generate_db_path):
    # try the sibling dir (some checkouts use different naming)
    generate_db_path = os.path.join(project_root, "generate_db.py")

if not os.path.exists(generate_db_path):
    LOGGER.error(
        "Could not find generate_db.py to import helper functions. Expected at %s",
        generate_db_path,
    )
    raise SystemExit(1)

spec = importlib.util.spec_from_file_location("_generate_db_module", generate_db_path)
generate_db_mod = importlib.util.module_from_spec(spec)
try:
    spec.loader.exec_module(generate_db_mod)  # type: ignore
except Exception as e:
    LOGGER.error("Failed to import generate_db.py: %s", e)
    LOGGER.debug(traceback.format_exc())
    raise

# Grab helpers
sanitize_fn = getattr(generate_db_mod, "sanitize_loaded_reactions", None)
cobra_reconstruction_fn = getattr(generate_db_mod, "cobra_reconstruction", None)

# Ensure the imported module has a LOGGER object used by its functions
try:
    if not hasattr(generate_db_mod, "LOGGER"):
        generate_db_mod.LOGGER = LOGGER
except Exception:
    LOGGER.debug("Could not set LOGGER on generate_db_mod", exc_info=True)

if sanitize_fn is None or cobra_reconstruction_fn is None:
    LOGGER.error(
        "Required functions not found in generate_db.py: sanitize_loaded_reactions or cobra_reconstruction"
    )
    raise SystemExit(1)

# Sanitize reaction objects that may contain pickled lambdas/closures
try:
    sanitize_fn(RxnList, name="RxnList")
    sanitize_fn(RxnList_CL, name="RxnList_CL")
    LOGGER.info("Sanitization complete")
except Exception as e:
    LOGGER.warning("Failed to sanitize reactions: %s", e)
    LOGGER.debug(traceback.format_exc())

# Reconstruct cobra model and write SBML
try:
    # Import cobra locally to avoid importing it at module import time
    # (this keeps errors localized and makes it clearer which step fails).
    import cobra

    model = cobra_reconstruction_fn(
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
except Exception as e:
    LOGGER.error("cobra_reconstruction failed: %s", e)
    LOGGER.debug(traceback.format_exc())
    raise

try:
    cobra.io.write_sbml_model(model, Output)
    LOGGER.info("Wrote reconstructed SBML to %s", Output)
except Exception as e:
    LOGGER.error("Failed to write SBML model to %s: %s", Output, e)
    LOGGER.debug(traceback.format_exc())
    raise
