import os
import sys
import dill
import logging
import importlib.util
import traceback
from pprint import pformat


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


# --- Debug: print detailed information for a couple of metabolites from the pickle
def _safe_call(obj, attr):
    try:
        val = getattr(obj, attr)
        if callable(val):
            try:
                return val()
            except TypeError:
                # some callables require args; skip
                return f"<callable {attr} requires arguments>"
            except Exception as e:
                return f"<error calling {attr}: {e}>"
        else:
            return val
    except Exception as e:
        return f"<missing {attr}: {e}>"


def print_compound_info(label, compound):
    LOGGER.info("--- Compound info: %s ---", label)
    try:
        LOGGER.info("Type: %s", type(compound))
        LOGGER.info("Repr: %s", repr(compound))
        # If object stores attributes in __dict__, print them
        if hasattr(compound, "__dict__"):
            try:
                LOGGER.info("__dict__: %s", pformat(compound.__dict__))
            except Exception:
                LOGGER.info("Could not pretty-print __dict__", exc_info=True)

        # Try a set of commonly-used accessor methods safely
        probe_methods = [
            "ID1",
            "ID2",
            "Name",
            "Formula1",
            "Formula2",
            "Formula3",
            "Formula4",
            "PubChem",
            "CheBI",
            "GlyDB",
            "JCGGDB",
            "inchi",
            "inchikey",
            "LipidBank",
            "LIPIDMAPS",
            "Atom1",
            "charge",
            "Subcel",
        ]
        for m in probe_methods:
            LOGGER.info("%s(): %s", m, _safe_call(compound, m))

        # Show available attributes/methods (shortened)
        try:
            names = [n for n in dir(compound) if not n.startswith("__")]
            LOGGER.info(
                "dir(): %s", ", ".join(names[:50]) + ("..." if len(names) > 50 else "")
            )
        except Exception:
            LOGGER.debug("Could not list dir() for compound", exc_info=True)
    except Exception as e:
        LOGGER.error("Failed to print compound info for %s: %s", label, e)


# Pick a couple of metabolites from the raw and compartmentalized lists and print them
try:
    sample_raw = list(MetList.keys())[:2]
    sample_cl = list(MetList_CL.keys())[:2]
    if sample_raw:
        LOGGER.info("Printing up to 2 metabolites from raw MetList: %s", sample_raw)
        for k in sample_raw:
            print_compound_info(f"raw:{k}", MetList.get(k))
    else:
        LOGGER.info("No entries found in MetList to print")

    if sample_cl:
        LOGGER.info(
            "Printing up to 2 metabolites from compartmentalized MetList_CL: %s",
            sample_cl,
        )
        for k in sample_cl:
            print_compound_info(f"cl:{k}", MetList_CL.get(k))
    else:
        LOGGER.info("No entries found in MetList_CL to print")
except Exception as e:
    LOGGER.warning("Failed while printing sample metabolites: %s", e)


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
