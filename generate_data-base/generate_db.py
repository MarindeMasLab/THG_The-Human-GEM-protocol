#!/usr/bin/python
"""
TODO:
Currently lambda functions are pickled for the checkpoint files. It raises an error
for missing localc objects, therefore, a sanitization step is added after loading.
This approach just patches the problem, but it should be changed so
lambdas are not pickled at all.

Generate metabolic model database from KEGG pathways.

Logging Levels:
- DEBUG: Show all detailed information about reactions, compounds, mass balance, etc.
- INFO: Show progress information only
- WARNING: Show warnings and errors
- ERROR: Show only errors

To change logging level, modify the level parameter in logging.basicConfig():
    logging.basicConfig(level=logging.DEBUG)  # Current setting - verbose
    logging.basicConfig(level=logging.INFO)   # Less verbose
    logging.basicConfig(level=logging.WARNING) # Minimal output

"""
import copy
import logging
import os
import pickle
import re
import urllib
from functools import reduce
from typing import Dict, List

import cobra
from cobra.io import read_sbml_model
import dill
import sys
import pdb
from dotenv import load_dotenv

# Determine the current file's directory and the project root.
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.join(current_dir, "..")

# Add the project root to sys.path to access top-level folders like 'functions' and 'models'
if project_root not in sys.path:
    sys.path.append(project_root)

# reimports for type hints
from functions.class_generate_database import *
from functions.class_generate_database import compound as CompoundType
from functions.class_generate_database import gene as GeneType
from functions.class_generate_database import reaction as ReactionType
from functions.pattern_generate_database import *
from functions.function_bm_gdb import *
from functions.equations_bm_gdb import *
from functions.function_bm_gdb import batch_fetch_kegg_entries
from types import MethodType


# Module-level helper methods to avoid fragile lambdas/closures when
# attaching callable accessors to reaction objects. These read the
# raw lists stored on the reaction (subs/prods) so they are safe to
# rebind after deserialization and avoid referencing local names like S2.
def _rxn_substrate(self):
    """Return substrate list stored on reaction instance.

    Returns the attribute `subs` if present, otherwise an empty list.
    """
    return getattr(self, "subs", [])


def _rxn_product(self):
    """Return product list stored on reaction instance.

    Returns the attribute `prods` if present, otherwise an empty list.
    """
    return getattr(self, "prods", [])


def sanitize_loaded_reactions(rxn_dict, name="reactions"):
    """Sanitize reaction objects loaded from disk.

    Ensures each reaction has `subs`/`prods` attributes and that callable
    accessors (`Substrate`, `Product`, `SetSubstrate`, `SetProduct`) are
    bound to stable module-level methods instead of fragile lambdas.
    """
    if not isinstance(rxn_dict, dict):
        return
    for key, rxn in list(rxn_dict.items()):
        try:
            # Try calling existing accessors to get data
            got = False
            try:
                if hasattr(rxn, "Substrate") and callable(rxn.Substrate):
                    subs = rxn.Substrate()
                    got = True
                else:
                    subs = getattr(rxn, "subs", None)
            except NameError as e:
                # Known bad closure (e.g. lambda referencing S2) — fallback
                LOGGER.warning(
                    "Reaction %s: Substrate accessor raised NameError: %s",
                    key,
                    e,
                )
                subs = getattr(rxn, "subs", None)

            try:
                if hasattr(rxn, "Product") and callable(rxn.Product):
                    prods = rxn.Product()
                    got = True
                else:
                    prods = getattr(rxn, "prods", None)
            except NameError as e:
                LOGGER.warning(
                    "Reaction %s: Product accessor raised NameError: %s",
                    key,
                    e,
                )
                prods = getattr(rxn, "prods", None)

            # Ensure lists exist
            subs = subs if subs is not None else []
            prods = prods if prods is not None else []

            # Store raw lists and bind stable methods
            try:
                rxn.subs = subs
                rxn.prods = prods
            except Exception:
                # Some reaction objects may not allow attribute setting; skip
                LOGGER.debug(
                    "Could not set subs/prods on reaction %s", key, exc_info=True
                )

            try:
                rxn.Substrate = MethodType(_rxn_substrate, rxn)
                rxn.Product = MethodType(_rxn_product, rxn)
                rxn.SetSubstrate = MethodType(_rxn_substrate, rxn)
                rxn.SetProduct = MethodType(_rxn_product, rxn)
            except Exception:
                LOGGER.debug(
                    "Could not bind methods on reaction %s", key, exc_info=True
                )

        except Exception:
            LOGGER.warning(f"Failed to sanitize loaded reaction {key}", exc_info=True)


def cobra_reconstruction(
    model_name: str,
    model_id: str,
    metabolite_list: Dict[List, CompoundType],
    reaction_list: Dict[List, ReactionType],
    gene_list: Dict[List, GeneType],
    pathways: Dict[str, str],
    location_dict: Dict[str, str],
    metabolite_equivalent: Dict[str, str],
    metabolite_list_general: Dict[List, CompoundType],
) -> cobra.Model:
    """Reconstruction of gathered information given by using cobrapy.

    Parameters
    ----------
    model_name: str
    model_name: id
    metabolite_list: Dict[List, compound]
    reaction_list: Dict[List, reaction]
    gene_list: Dict[List, gene]
    pathways: Dict[str, str]
        Map from Pathway to Reaction identifiers (PathNameRxn).
    location_dict: Dict[str, str]
        generated from the Comparment_Cl. This stores human readable
        compartment names to comparment identifiers. The extracted identifiers
        in the reactions and metabolites contain the id and the comparment
        name, so this mapping is necessary to achieve a proper
        SBML-compatible identifier (no spaces).
    """
    model = cobra.Model(model_id or model_name, model_name or model_id)
    location_dict = {k.lower(): v for k, v in location_dict.items()}
    LOGGER.debug(
        "Starting cobra_reconstruction: metabolites=%s reactions=%s genes=%s pathways=%s loc=%s",
        (
            len(metabolite_list)
            if hasattr(metabolite_list, "__len__")
            else type(metabolite_list)
        ),
        (
            len(reaction_list)
            if hasattr(reaction_list, "__len__")
            else type(reaction_list)
        ),
        len(gene_list) if hasattr(gene_list, "__len__") else type(gene_list),
        len(pathways) if hasattr(pathways, "__len__") else type(pathways),
        (
            len(location_dict)
            if hasattr(location_dict, "__len__")
            else type(location_dict)
        ),
    )
    try:
        LOGGER.debug("Sample metabolite keys: %s", list(metabolite_list.keys())[:10])
    except Exception:
        LOGGER.debug("Could not list metabolite_list keys", exc_info=True)
    try:
        LOGGER.debug("Sample reaction keys: %s", list(reaction_list.keys())[:10])
    except Exception:
        LOGGER.debug("Could not list reaction_list keys", exc_info=True)
    # metabolites
    compounds = [
        cobra.Metabolite(
            compound.ID2() + "_" + location_dict.get(compound.Subcel.lower()),
            compound.Formula1(),
            compound.Name(),
            float(compound.charge()),
            location_dict.get(compound.Subcel.lower()),
        )
        for iden, compound in metabolite_list.items()
    ]
    model.add_metabolites(compounds)
    # store mapping (met identifier -> met.id in model) for reaction section
    met_mapping = {}
    # metabolite annotation
    for iden, compound in metabolite_list.items():
        model_met = model.metabolites.get_by_id(
            compound.ID2() + "_" + location_dict.get(compound.Subcel.lower())
        )
        met_mapping[iden] = model_met.id
        annotation = {
            "pubchem.compound": compound.PubChem(),
            "chebi.compound": compound.CheBI(),
            "glycomedb": compound.GlyDB(),
            "jcggdb": compound.JCGGDB(),
            "inchi": compound.inchi(),
            "inchikey": compound.inchikey(),
            "lipidbank": compound.LipidBank(),
            "lipidmaps": compound.LIPIDMAPS(),
        }
        alt_formulas = [compound.Formula2(), compound.Formula3(), compound.Formula4()]
        alt_formula = [
            form for form in alt_formulas if form != model_met.formula and form
        ]
        if alt_formula:
            model_met.annotation["glycan_formula"] = alt_formula[0]

        # only add non-empty annotation
        model_met.annotation = {k: v for k, v in annotation.items() if v}

    for x in model.metabolites:  # replace glycan formula by a sbml suitable format
        if x.id[0] == "G":
            try:
                print(x.id)
                compartment = [
                    y[0] for y in location_dict.items() if y[1] in x.id.split("_")[1]
                ][0]
                xth_metabolite_id = x.id.split("_")[0] + "_" + compartment
                x.formula = metabolite_list[xth_metabolite_id].Formula4()
            except Exception as e:
                import traceback

                LOGGER.error(f"Error updating glycan formula for {x.id}: {e}")
                LOGGER.error(traceback.format_exc())
                continue

    for (
        x
    ) in (
        model.metabolites
    ):  # eliminate potential discrepancies between metabolite id and reaction compounds ids
        if x.id.split("_")[0] in metabolite_equivalent.keys():
            x.id = (
                metabolite_list_general[metabolite_equivalent[x.id.split("_")[0]]].ID1()
                + "_"
                + x.compartment.lower()
            )

    def normalize_id(reac_id: str):
        iden, comp_desc = reac_id.split("_")
        comp = location_dict[comp_desc.lower()]
        return f"{iden}_{comp}"

    from functions.gpr.ast_gpr import sanitize_gpr

    # reactions
    reactions = [
        cobra.Reaction(
            normalize_id(iden), reac.Name(), "", 0 if reac.Termodyn() else -1000, 1000
        )
        for iden, reac in reaction_list.items()
    ]
    model.add_reactions(reactions)
    for iden, rxn in reaction_list.items():
        print(iden)
        reac = model.reactions.get_by_id(normalize_id(iden))
        rxn_id = rxn.ID
        if "_" in rxn_id:
            kegg_id = rxn_id.split("_")[0]
            comp_id = "_" + location_dict[rxn_id.split("_")[-1].lower()]
        else:
            kegg_id = rxn_id
            comp_id = ""
        try:
            # this may fail after serialization because these are overwritten
            # at runtime via a capturing lambda
            # substrates might come with positive coefficients
            substrates = [
                (-abs(convert_to_float(subs[0])), subs[1], subs[2])
                for subs in rxn.Substrate()
            ]
            products = [
                (abs(convert_to_float(prod[0])), prod[1], prod[2])
                for prod in rxn.Product()
            ]
            reac_compounds = products + substrates
        except Exception as e:
            # these metabolites do not have an specified comparment!
            # substrates might come with positive coefficients
            import traceback

            LOGGER.warning(
                f"Error getting reaction compounds from Product/Substrate methods, using subs/prods: {e}"
            )
            LOGGER.warning(traceback.format_exc())
            substrates = [
                (-abs(convert_to_float(subs[0])), subs[1], subs[2]) for subs in rxn.subs
            ]
            products = [
                (abs(convert_to_float(prod[0])), prod[1], prod[2]) for prod in rxn.prods
            ]
            reac_compounds = products + substrates
        metabolites = {
            (
                met_mapping[met[2]] if met[2] in met_mapping else f"{met[2]}{comp_id}"
            ): float(met[0])
            for met in reac_compounds
        }
        # there may be some metabolites that were not passed in metabolite list
        new_mets = []
        for met_id in metabolites:
            if met_id not in model.metabolites:
                # first check if we have them in a different comparment
                met_root, met_comp = met_id.split("_")[0], met_id.split("_")[1]
                maybe_mets = [m for m in model.metabolites if m.id.startswith(met_root)]
                if maybe_mets:
                    model_met = maybe_mets[0]
                    LOGGER.warning(
                        f"Reactant '{met_root}' was added in a different comparment '{met_comp}'."
                    )
                    new_met = cobra.Metabolite(
                        met_id,
                        model_met.formula,
                        model_met.name,
                        model_met.charge,
                        # might be a new compartment!
                        met_comp,
                    )
                    new_met.annotation = model_met.annotation
                else:
                    LOGGER.warning(
                        f"Reactant {met_id} was not found in any compartment. Creating new one!"
                    )
                    new_met = cobra.Metabolite(
                        met_id,
                        compartment=met_comp,
                    )
                new_mets.append(new_met)
        if new_mets:
            model.add_metabolites(new_mets)

        reac.add_metabolites(metabolites)
        ec = rxn.EC()
        reac.annotation = {"kegg.reaction": kegg_id, "ec-code": ec[0] if ec else ""}
        sgpr, gpr = rxn.GPR[0], rxn.GPR[1].replace("[", "").replace("]", "")
        if sgpr and sgpr != "[]":
            if "or" in gpr and "and" not in gpr:
                # GPRs scrapped from Kegg are added with ORs and the genes
                # may be repeated so they have to be deduplicated
                # TODO(carrascomj): should come from getGPR / getLocation
                gpr = " or ".join({gene for gene in gpr.split(" or ") if gene})
            # Sanitize GPR but guard against malformed strings from KEGG
            try:
                reac.gene_reaction_rule = sanitize_gpr(gpr)
            except Exception as e:
                import traceback

                LOGGER.warning(
                    f"Failed to sanitize GPR for reaction {rxn_id if 'rxn_id' in locals() else iden}: {gpr!r}: {e}"
                )
                LOGGER.warning(traceback.format_exc())
                # Fallback: leave gene reaction rule empty so processing continues
                reac.gene_reaction_rule = ""
            reac.annotation["sGPR"] = sgpr
        reac.id = kegg_id + comp_id
    # add a group per pathway
    model.add_groups([cobra.core.Group(group, group) for group in pathways])
    for group, members in pathways.items():
        # the members are the reactions in each pathway
        # TODO(carrascomj): reaction ids coming from paths are not in compartments
        model.groups.get_by_id(group).add_members(
            reduce(
                lambda x, y: x + y,
                [model.reactions.query(member) for member in members.split()],
                [],
            )
        )
    # gene annotation (genes were added with the GPRs)
    pat_enstp = re.compile("ENS[TP][0-9]+")
    for iden, gene in gene_list.items():
        print(iden)

        if not iden in model.genes:
            LOGGER.warning(f"Gene '{iden}' was not found. Creating new one!")
            model_gene = cobra.Gene(
                iden,
                gene.Name().replace("[", "").replace("]", "").replace("-", ""),
            )

        else:
            model_gene = model.genes.get_by_id(iden)
            model_gene.name = (
                gene.Name().replace("[", "").replace("]", "").replace("-", "")
            )

        # there may be other ensembl genes, we have to query them
        ensembl = gene.Ensg()
        ensembl_genes = []
        if ensembl:
            try:
                retrieved_ids = str(
                    urllib.request.urlopen(
                        "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=" + ensembl
                    ).read()
                )  # only add non-empty annotation
                ensembl_genes = [
                    str(x) for x in list(set(pat_enstp.findall(retrieved_ids)))
                ] + [ensembl]
            except Exception as e:
                import traceback

                LOGGER.warning(f"Error fetching ensembl genes for {ensembl}: {e}")
                LOGGER.warning(traceback.format_exc())
                ensembl_genes = []

        model_gene.annotation = {
            k: v
            for k, v in {
                "ensembl": ensembl_genes,
                "ncbigene": gene.Entrez(),
                "uniprot": gene.Uniprot(),
                "hgcn.symbol": model_gene.name,
            }.items()
            if v
        }
    return model


if __name__ == "__main__":

    # Load environment variables from .env file
    env_file = os.path.join(project_root, ".env")
    if os.path.exists(env_file):
        load_dotenv(env_file, override=True)
        print(f"Loaded environment variables from {env_file}")

    # Setup logging
    LOGGER = logging.getLogger(__name__)

    # Create formatters and handlers
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    console_handler.setFormatter(formatter)

    # File handler - write to the same log file
    log_file = os.path.join(project_root, "logs", "generate_db.log")
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    file_handler = logging.FileHandler(log_file, mode="a")  # Append mode
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    # Configure root logger
    logging.basicConfig(
        level=logging.DEBUG,  # Change to INFO, WARNING, or ERROR to reduce output
        handlers=[console_handler, file_handler],
    )

    session = setup_biocyc_session()

    #### Initial Parameters
    # Determine the current file's directory and the project root.
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.join(current_dir, "..")

    # Add the project root to sys.path to access top-level folders like 'functions' and 'models'
    if project_root not in sys.path:
        sys.path.append(project_root)

    ListOfPaths = os.path.join(project_root, "files", "human_kegg_pathways.txt")
    ModelCompounds = os.path.join(project_root, "files", "extra_compounds.txt")
    ExtraFormula = os.path.join(project_root, "files", "extra_formula.txt")
    ModelReactions = ""
    ModelGenes = ""
    EnsblDB = os.path.join(
        project_root, "files", "ensembl"
    )  # From Ensembl database: ensembl gene ID vs Entrez vs Name.
    time = 20  # Time to download url: Parameter defined in function getHtml
    Path = (
        open(ListOfPaths, "r").read().split("\n")
    )  # From analysis using metaboanalyst.
    Compound = open(ModelCompounds, "r").read().split("\n")
    EF = [_f for _f in open(ExtraFormula, "r").read().split("\n") if _f]
    Output = os.path.join(project_root, "models", "Human_Database.xml")  # Output model
    ModID = Output
    ModName = Output
    variablesFile = os.path.join(
        project_root, "files", "model_variables.pkl"
    )  # File where the working environment is saved
    specialCompounds = os.path.join(
        project_root, "files", "special_compounds.txt"
    )  # File where we save the IDs of the compounds with a (group)n in their formula
    open(specialCompounds, "w").close()  # Erase or create the file

    # Dictionary with extra compounds that can be added to mass balance the metabolic reactions
    extra_compound = {
        "H": "C00080",
        "H2O": "C00001",
        "Fe": "C00023",
        "Na": "C01330",
        "Ca": "C00076",
        "K": "C00238",
        "F": "C00023",
        "R": "C00000",
        "X": "C0000X",
    }

    #### Initial List and dictionaries
    PathList = {}
    RxnList = {}
    MetList = {}
    GPRList = {}
    MetEquiv = {}
    PathNameRxn = {}
    RxnEquiv = {}
    PathIdent = []
    RxnIdent = []
    MetIdent = []
    GPRIdent = []
    RxnIDList = []
    MetIDList = []

    #### List and dictionaries for the subcelular location annotation
    CSL_ID = {
        "extracellular": "e",
        "peroxisome": "x",
        "mitochondria": "m",
        "cytosol": "c",
        "lysosome": "l",
        "endoplasmic reticulum": "r",
        "golgi apparatus": "g",
        "nucleus": "n",
        "inner mitochondria": "i",
    }  # to keep the consistency between the DB and the initial compartments in Human1
    CSL_ID = compartment_file_to_dict()
    CSL_ID = dict((k.lower(), v.lower()) for k, v in CSL_ID.items())
    listOfID = []
    LocVar = {}
    RxnList_CL = {}
    MetList_CL = {}
    RxnIdent_CL = []
    MetIdent_CL = []
    Compartment_CL = []
    RxnList_Subcel = []

    ######### Checkpoint file for resuming progress ###########
    checkpoint_file = os.path.join(project_root, "files", "checkpoint_progress.pkl")

    # Try to load checkpoint if it exists
    start_pathway_index = 0
    if os.path.exists(checkpoint_file):
        try:
            with open(checkpoint_file, "rb") as f:
                checkpoint_data = dill.load(f)
                start_pathway_index = (
                    checkpoint_data.get("last_completed_pathway", 0) + 1
                )
                PathList = checkpoint_data.get("PathList", PathList)
                RxnList = checkpoint_data.get("RxnList", RxnList)
                MetList = checkpoint_data.get("MetList", MetList)
                GPRList = checkpoint_data.get("GPRList", GPRList)
                MetEquiv = checkpoint_data.get("MetEquiv", MetEquiv)
                PathNameRxn = checkpoint_data.get("PathNameRxn", PathNameRxn)
                RxnEquiv = checkpoint_data.get("RxnEquiv", RxnEquiv)
                PathIdent = checkpoint_data.get("PathIdent", PathIdent)
                RxnIdent = checkpoint_data.get("RxnIdent", RxnIdent)
                MetIdent = checkpoint_data.get("MetIdent", MetIdent)
                GPRIdent = checkpoint_data.get("GPRIdent", GPRIdent)
                RxnList_CL = checkpoint_data.get("RxnList_CL", RxnList_CL)
                MetList_CL = checkpoint_data.get("MetList_CL", MetList_CL)
                RxnIdent_CL = checkpoint_data.get("RxnIdent_CL", RxnIdent_CL)
                MetIdent_CL = checkpoint_data.get("MetIdent_CL", MetIdent_CL)
                Compartment_CL = checkpoint_data.get("Compartment_CL", Compartment_CL)
                LOGGER.info(
                    f"Resuming from pathway index {start_pathway_index} (pathway {start_pathway_index + 1}/{len(Path) - 1})"
                )
                print(
                    f"Resuming from pathway index {start_pathway_index} (pathway {start_pathway_index + 1}/{len(Path) - 1})"
                )
                # Sanitize any loaded reaction objects to remove fragile lambdas
                try:
                    sanitize_loaded_reactions(RxnList, name="RxnList")
                    sanitize_loaded_reactions(RxnList_CL, name="RxnList_CL")
                except Exception:
                    LOGGER.debug(
                        "Failed to sanitize reactions loaded from checkpoint",
                        exc_info=True,
                    )
        except Exception as e:
            LOGGER.warning(
                f"Could not load checkpoint file: {e}. Starting from beginning."
            )
            print(f"Could not load checkpoint file: {e}. Starting from beginning.")
            start_pathway_index = 0

    ######### Pathways ###########
    i = start_pathway_index
    while i < len(Path) - 1:
        # print(Path[i])
        #### Build network
        PathID = Path[i].split("\t")[0]
        PathName = Path[i].split("\t")[1]
        PathURL = "https://rest.kegg.jp/get/" + PathID + "/kgml"
        PathReferer = "https://www.kegg.jp/kegg-bin/show_pathway?" + PathID
        PathList[PathID] = pathway(PathURL, time, PathID, PathReferer, PathName)
        if PathList[PathID].Compounds():
            print(
                PathName
                + ": defined in human"
                + "("
                + str(i + 1)
                + "/"
                + str(len(Path) - 1)
                + ")"
            )
            PathNameRxn[PathName] = ""
            j = 0
            while j < len(PathList[PathID].Reactions()):
                try:
                    RxnID = PathList[PathID].Reactions()[j][0][0]
                    if not RxnID in RxnIdent and not RxnID in RxnEquiv:

                        ######### Define New Reaction ###########
                        LOGGER.debug("=" * 80)
                        LOGGER.debug(
                            f"Processing Reaction {RxnID} from pathway {PathName}"
                        )
                        LOGGER.debug("=" * 80)
                        RxnIdent.append(
                            RxnID
                        )  # Optimized: use append instead of concatenation
                        RxnURL = PathList[PathID].Reactions()[j][1]
                        RxnTermDyn = PathList[PathID].Reactions()[j][0][1]
                        RxnList[RxnID] = reaction(
                            RxnURL, time, RxnID, PathName, RxnTermDyn
                        )

                        # Debug: Show initial reaction data from KEGG
                        LOGGER.debug(f"Reaction Name: {RxnList[RxnID].Name()}")
                        LOGGER.debug(f"Reaction EC: {RxnList[RxnID].EC()}")
                        LOGGER.debug(f"Thermodynamic: {RxnTermDyn}")
                        LOGGER.debug(f"Initial Substrates (raw from KEGG):")
                        for sub in RxnList[RxnID].Substrate():
                            LOGGER.debug(
                                f"  - Coeff: {sub[0]}, Link: {sub[1][:50]}..., ID: {sub[2]}"
                            )
                        LOGGER.debug(f"Initial Products (raw from KEGG):")
                        for prod in RxnList[RxnID].Product():
                            LOGGER.debug(
                                f"  - Coeff: {prod[0]}, Link: {prod[1][:50]}..., ID: {prod[2]}"
                            )

                        ######### Check if all the compounds in the jth reaction are in the compound list ###########
                        RxnCmp = [x[2] for x in RxnList[RxnID].Substrate()] + [
                            x[2] for x in RxnList[RxnID].Product()
                        ]

                        LOGGER.debug(f"Compounds in reaction {RxnID}: {RxnCmp}")
                        LOGGER.debug(f"Total unique compounds: {len(set(RxnCmp))}")

                        # Collect new compound IDs for concurrent fetching
                        new_compounds_to_fetch = []
                        c = 0
                        while c < len(RxnCmp):
                            CompID = RxnCmp[c]
                            if not CompID in MetIdent and not CompID in MetEquiv:
                                MetIdent.append(
                                    CompID
                                )  # Optimized: use append instead of concatenation
                                new_compounds_to_fetch.append(CompID)
                            c = c + 1

                        LOGGER.debug(
                            f"New compounds to fetch: {new_compounds_to_fetch}"
                        )

                        # Fetch compounds using optimized batch + concurrent requests
                        if new_compounds_to_fetch:

                            # Separate glycans (G-prefix) from regular compounds (C-prefix)
                            glycans_to_fetch = [
                                cid
                                for cid in new_compounds_to_fetch
                                if cid.startswith("G")
                            ]
                            compounds_to_fetch = [
                                cid
                                for cid in new_compounds_to_fetch
                                if cid.startswith("C")
                            ]

                            if glycans_to_fetch:
                                LOGGER.debug(
                                    f"Fetching {len(glycans_to_fetch)} glycans: {glycans_to_fetch}"
                                )
                            if compounds_to_fetch:
                                LOGGER.debug(
                                    f"Fetching {len(compounds_to_fetch)} compounds: {compounds_to_fetch}"
                                )

                            batch_data = {}

                            # Fetch regular compounds using REST API
                            if compounds_to_fetch:
                                compound_batch_data = batch_fetch_kegg_entries(
                                    compounds_to_fetch,
                                    database="compound",
                                    batch_size=10,  # KEGG API limit per request
                                    max_workers=5,  # Number of concurrent batch requests
                                )
                                batch_data.update(compound_batch_data)

                            # Fetch glycans using REST API
                            if glycans_to_fetch:
                                glycan_batch_data = batch_fetch_kegg_entries(
                                    glycans_to_fetch,
                                    database="glycan",
                                    batch_size=10,  # KEGG API limit per request
                                    max_workers=5,  # Number of concurrent batch requests
                                )
                                batch_data.update(glycan_batch_data)

                            LOGGER.debug(
                                f"Batch data retrieved: {len(batch_data)} entries"
                            )

                            for CompID in new_compounds_to_fetch:
                                try:
                                    LOGGER.debug(f"Processing compound {CompID}...")
                                    if CompID in batch_data and batch_data[CompID]:
                                        # Use pre-fetched data
                                        LOGGER.debug(f"Using batch data for {CompID}")
                                        MetList[CompID] = compound.from_batch_data(
                                            CompID,
                                            batch_data[CompID],
                                            time,
                                            EF,
                                            specialCompounds,
                                        )
                                    else:
                                        # Fallback to individual fetch if concurrent fetch failed
                                        LOGGER.warning(
                                            f"Batch fetch failed for {CompID}, using individual fetch"
                                        )
                                        CompURL = "https://rest.kegg.jp/get/" + CompID
                                        MetList[CompID] = compound(
                                            CompURL, CompID, time, EF, specialCompounds
                                        )

                                    # Debug: Show compound details
                                    LOGGER.debug(f"Compound {CompID} details:")
                                    LOGGER.debug(f"  - ID1: {MetList[CompID].ID1()}")
                                    LOGGER.debug(f"  - ID2: {MetList[CompID].ID2()}")
                                    LOGGER.debug(f"  - Name: {MetList[CompID].Name()}")
                                    LOGGER.debug(
                                        f"  - Formula1: {MetList[CompID].Formula1()}"
                                    )
                                    LOGGER.debug(
                                        f"  - Formula2: {MetList[CompID].Formula2()}"
                                    )
                                    if (
                                        MetList[CompID].ID1()
                                        and MetList[CompID].ID1()[0] == "G"
                                    ):
                                        LOGGER.debug(
                                            f"  - [GLYCAN] Formula4 (reformulated): {MetList[CompID].Formula4()}"
                                        )
                                    LOGGER.debug(
                                        f"  - Atom composition: {MetList[CompID].Atom1()}"
                                    )

                                    # Handle ID equivalences
                                    if MetList[CompID].ID1() != MetList[CompID].ID2():
                                        LOGGER.debug(
                                            f"ID equivalence found: {CompID} -> {MetList[CompID].ID1()}"
                                        )
                                        MetIdent[len(MetIdent) - 1] = MetList[
                                            CompID
                                        ].ID1()
                                        MetEquiv[CompID] = MetList[CompID].ID1()
                                        # Direct assignment instead of deepcopy when possible
                                        MetList[MetList[CompID].ID1()] = MetList[CompID]
                                        del MetList[CompID]
                                except Exception as e:
                                    LOGGER.error(
                                        f"Error processing compound {CompID}: {e}"
                                    )
                                    import traceback

                                    LOGGER.error(traceback.format_exc())
                                    raise  # Re-raise to trigger outer exception handler

                        ######### Define Substrates, Products and New Compounds ###########
                        # Evaluate the relation between substrates and products #
                        LOGGER.debug(f"Calling getRxncons to evaluate reaction {RxnID}")
                        Rxn = getRxncons(
                            RxnList[RxnID],
                            time,
                            MetEquiv,
                            MetList,
                            MetIdent,
                            EF,
                            specialCompounds,
                        )
                        RxnList[RxnID] = Rxn  # Removed unnecessary deepcopy

                        # Check reaction ID
                        if RxnID != RxnList[RxnID].ID:
                            LOGGER.debug(
                                f"Reaction ID changed: {RxnID} -> {RxnList[RxnID].ID} (Glycan -> Compound equivalent)"
                            )
                            RxnEquiv[RxnID] = RxnList[
                                RxnID
                            ].ID  # Glycan Reaction : Compound Reaction
                            RxnIdent[len(RxnIdent) - 1] = RxnList[RxnID].ID
                            RxnList[RxnList[RxnID].ID] = RxnList[
                                RxnID
                            ]  # Removed unnecessary deepcopy
                            tmpID = RxnList[RxnList[RxnID].ID].ID
                            del RxnList[RxnID]  # remove the old reaction ID
                            RxnID = tmpID

                        ######### Mass Balance the reaction #########
                        ithRxn = RxnList[RxnID]
                        eq, mb_test = RxnParam2Eq(ithRxn, MetList, MetEquiv)
                        LOGGER.debug(f"Mass balance test result for {RxnID}: {mb_test}")
                        LOGGER.debug(f"Reaction equation: {eq}")
                        LibIni = WrapRxnSubsProdParam(ithRxn, MetList, MetEquiv)
                        if mb_test != 0:
                            IthRxnMB = mass_balance(eq, RxnID)
                            LOGGER.debug(f"Mass balance result: {IthRxnMB}")
                            if IthRxnMB[
                                4
                            ]:  # If new compounds have to be added to mass balance the reactions, then check if they need to be added to the network as compounds
                                LOGGER.debug(
                                    f"Adding extra compounds for mass balance: {IthRxnMB[4]}"
                                )
                                for x in IthRxnMB[4]:
                                    if (
                                        not extra_compound[x[0]] in MetIdent
                                        and not extra_compound[x[0]] in MetEquiv
                                    ):
                                        MetIdent.append(
                                            extra_compound[x[0]]
                                        )  # Optimized: use append
                                        if not x[0] in "R" and not x[0] in "X":
                                            extra_url = (
                                                "https://www.genome.jp/entry/"
                                                + extra_compound[x[0]]
                                            )
                                            MetList[extra_compound[x[0]]] = compound(
                                                extra_url,
                                                extra_compound[x[0]],
                                                time,
                                                EF,
                                                specialCompounds,
                                            )
                                        else:
                                            MetList[extra_compound[x[0]]] = (
                                                add_extra_compound(
                                                    x[0],
                                                    extra_compound,
                                                    time,
                                                    EF,
                                                    specialCompounds,
                                                )
                                            )
                        else:  # if the reaction cannot be mass balanced all the stoichimetric coef are assumed to be like in the original reaction
                            LOGGER.warning(
                                f"Reaction {RxnID} cannot be mass balanced - using original stoichiometry"
                            )
                            IthRxnMB = (
                                [float(x[0]) for x in ithRxn.Substrate()],
                                [float(x[0]) for x in ithRxn.Product()],
                                0,
                                0,
                                0,
                                0,
                                [
                                    MetEquiv[x[2]] if x[2] in MetEquiv else x[2]
                                    for x in ithRxn.Substrate()
                                ],
                                [
                                    MetEquiv[x[2]] if x[2] in MetEquiv else x[2]
                                    for x in ithRxn.Product()
                                ],
                                "",
                                "",
                                0,
                            )
                        LibEnd = UnwrapRxnSubsProdParam(IthRxnMB, LibIni, IthRxnMB)

                        # Add metabolites and stc coeff to reaction
                        S = list()
                        for x in LibEnd[0]:
                            S.append(
                                [
                                    str(LibEnd[0][x][0]),
                                    (
                                        "http://www.genome.jp/dbget-bin/www_bget?cpd:"
                                        + LibEnd[0][x][1]
                                    ),
                                    LibEnd[0][x][1],
                                ]
                            )
                        P = list()
                        for x in LibEnd[1]:
                            P.append(
                                [
                                    str(LibEnd[1][x][0]),
                                    (
                                        "http://www.genome.jp/dbget-bin/www_bget?cpd:"
                                        + LibEnd[1][x][1]
                                    ),
                                    LibEnd[1][x][1],
                                ]
                            )

                        # Debug: Show final reaction structure
                        LOGGER.debug(f"Final Reaction Structure for {RxnID}:")
                        LOGGER.debug(f"  Substrates:")
                        for sub in S:
                            met_id = sub[2]
                            met_name = (
                                MetList[met_id].Name()
                                if met_id in MetList
                                else "Unknown"
                            )
                            LOGGER.debug(f"    {sub[0]} {met_id} ({met_name})")
                        LOGGER.debug(f"  Products:")
                        for prod in P:
                            met_id = prod[2]
                            met_name = (
                                MetList[met_id].Name()
                                if met_id in MetList
                                else "Unknown"
                            )
                            LOGGER.debug(f"    {prod[0]} {met_id} ({met_name})")

                        S2 = copy.deepcopy(S)
                        P2 = copy.deepcopy(P)
                        # Log the substrate/product lists being assigned for this reaction
                        LOGGER.debug(
                            "Assigning substrates/products for reaction %s: substrates=%s products=%s",
                            RxnID,
                            S2,
                            P2,
                        )
                        # Bind module-level methods to the reaction instance so
                        # callable accessors are stable and don't close over
                        # local names (avoids NameError after pickling).
                        RxnList[RxnID].subs = S2
                        RxnList[RxnID].prods = P2
                        RxnList[RxnID].Substrate = MethodType(
                            _rxn_substrate, RxnList[RxnID]
                        )
                        RxnList[RxnID].Product = MethodType(
                            _rxn_product, RxnList[RxnID]
                        )
                        RxnList[RxnID].SetSubstrate = MethodType(
                            _rxn_substrate, RxnList[RxnID]
                        )
                        RxnList[RxnID].SetProduct = MethodType(
                            _rxn_product, RxnList[RxnID]
                        )

                        RxnList[RxnID].subs = S2
                        RxnList[RxnID].prods = P2
                        if not RxnID in PathNameRxn.get(PathName):
                            PathNameRxn[PathName] += RxnID + " "

                        ######### Define New GPR ###########
                        tmpGPR = ()
                        tmpSC = ()
                        for x in RxnList[RxnID].EC():
                            if x not in GPRIdent:
                                GPRIdent.append(x)  # Optimized: use append
                                GPRList[x] = gpr(x, session)
                            try:
                                gpr_result = GPRList[x].GprSubcell()
                                # Check if we got a valid tuple result (not empty string)
                                if (
                                    gpr_result
                                    and isinstance(gpr_result, tuple)
                                    and len(gpr_result) >= 4
                                ):
                                    tmpGPR = tmpGPR + gpr_result[0:2]
                                    tmpSC = tmpSC + gpr_result[2:4]
                            except Exception as e:
                                import traceback

                                LOGGER.warning(f"Could not process GPR for EC {x}: {e}")
                                LOGGER.warning(traceback.format_exc())
                                continue
                        # Reorganize S-GPRs and GPRs based on their specific location
                        tmpSC2 = [dict(), dict()]
                        reactio_compartment_list = list(
                            set(
                                [
                                    x.strip()
                                    for x in str([list(x.keys()) for x in tmpSC])
                                    .replace("[", "")
                                    .replace("]", "")
                                    .replace("'", "")
                                    .split(",")
                                ]
                            )
                        )
                        for x in reactio_compartment_list:
                            xth_tmp_gpr = [y for y in tmpSC if x in y.keys()]
                            tmp_xth_sgpr = ""
                            tmp_xth_gpr = ""
                            for y in range(int(len(xth_tmp_gpr) / 2)):
                                if not re.findall("^\[\]$", xth_tmp_gpr[y + y][x]):
                                    tmp_xth_sgpr += xth_tmp_gpr[y + y][x]
                                if not re.findall("^\[\]$", xth_tmp_gpr[y + y + 1][x]):
                                    tmp_xth_gpr += xth_tmp_gpr[y + y + 1][x]
                            tmp_xth_sgpr = (
                                str(
                                    set(
                                        tmp_xth_sgpr.replace("][", "] or [").split(
                                            " or "
                                        )
                                    )
                                )
                                .replace("'", "")
                                .replace("{", "")
                                .replace("}", "")
                                .replace(",", " or")
                            )
                            tmp_xth_gpr = (
                                str(
                                    set(
                                        tmp_xth_gpr.replace("][", "] or [").split(
                                            " or "
                                        )
                                    )
                                )
                                .replace("'", "")
                                .replace("{", "")
                                .replace("}", "")
                                .replace(",", " or")
                            )
                            tmpSC2[0][x] = tmp_xth_sgpr
                            tmpSC2[1][x] = tmp_xth_gpr

                        RxnList[RxnID].GPR = tmpGPR
                        RxnList[RxnID].Subcel = tmpSC2

                        ######### Expand the annotations based on the cellular location ###########
                        Compartment_CL, rxn_cl, comp_cl = rxnSubcel(
                            RxnList[RxnID],
                            RxnList_CL,
                            MetList_CL,
                            RxnIdent_CL,
                            MetIdent_CL,
                            Compartment_CL,
                            RxnList,
                            MetList,
                            MetEquiv,
                        )

                        RxnIdent_CL.extend(
                            rxn_cl
                        )  # Optimized: use extend instead of concatenation
                        MetIdent_CL.extend(
                            comp_cl
                        )  # Optimized: use extend instead of concatenation

                        print(
                            "Reaction("
                            + str(j + 1)
                            + "/"
                            + str(len(PathList[PathID].Reactions()))
                            + ")_Pathway"
                            + "("
                            + str(i + 1)
                            + "/"
                            + str(len(Path) - 1)
                            + ")"
                        )
                except Exception as e:
                    LOGGER.error(
                        f"Failed to process reaction {RxnID if 'RxnID' in locals() else 'unknown'}: {e}"
                    )
                    import traceback

                    LOGGER.error(traceback.format_exc())
                    continue
                j = j + 1
        else:
            print(
                PathName
                + ": not defined in human"
                + "("
                + str(i + 1)
                + "/"
                + str(len(Path) - 1)
                + ")"
            )

        # Save checkpoint after each pathway
        try:
            checkpoint_data = {
                "last_completed_pathway": i,
                "PathList": PathList,
                "RxnList": RxnList,
                "MetList": MetList,
                "GPRList": GPRList,
                "MetEquiv": MetEquiv,
                "PathNameRxn": PathNameRxn,
                "RxnEquiv": RxnEquiv,
                "PathIdent": PathIdent,
                "RxnIdent": RxnIdent,
                "MetIdent": MetIdent,
                "GPRIdent": GPRIdent,
                "RxnList_CL": RxnList_CL,
                "MetList_CL": MetList_CL,
                "RxnIdent_CL": RxnIdent_CL,
                "MetIdent_CL": MetIdent_CL,
                "Compartment_CL": Compartment_CL,
            }
            with open(checkpoint_file, "wb") as f:
                dill.dump(checkpoint_data, f)
            LOGGER.debug(f"Checkpoint saved after pathway {i + 1}/{len(Path) - 1}")
        except Exception as e:
            LOGGER.warning(f"Failed to save checkpoint: {e}")

        i = i + 1

    ######### Genes ###########
    GeneList = {}
    GeneIdent = []
    g = 0
    while g < len(GPRList):
        if GPRIdent[g] in GPRList.keys() and GPRList[GPRIdent[g]].GprSubcell():
            gene_matches = re.findall(
                "([A-Za-z0-9\-]+)",
                GPRList[GPRIdent[g]]
                .GprSubcell()[1]
                .replace("and", "")
                .replace("or", ""),
            )
            z = 0
            while z < len(gene_matches):
                if not gene_matches[z] in GeneIdent:
                    GeneIdent.append(gene_matches[z])  # Optimized: use append
                    GeneList[gene_matches[z]] = gene(gene_matches[z], EnsblDB)
                z = z + 1
        g = g + 1

    with open(os.path.join(project_root, "files", "pre_sbml_raw.pk"), "wb") as f:
        dill.dump(
            {
                "name": ModName,
                "id": ModID,
                "mets": MetList,
                "mets_cl": MetList_CL,
                "met_equiv": MetEquiv,
                "reactions": RxnList,
                "reactions_cl": RxnList_CL,
                "genes": GeneList,
                "pathways": PathNameRxn,
                "loc": LocVar,
            },
            f,
        )

    Compartment_CL = sorted(Compartment_CL)

    listOfID = list(CSL_ID.values())  # abbr. id
    LipidMasterlistOfID = []

    for CSL in Compartment_CL:
        CSL2 = re.sub(r"[^A-Za-z0-9 ]+", "", CSL)
        ModMaster = list(set(LipidMasterlistOfID + listOfID))
        if CSL_ID.get(CSL):
            ID = CSL_ID.get(CSL)
        elif len(re.sub(" $", "", re.sub("^ ", "", CSL)).split(" ")) > 1:
            ID = (
                (CSL2.split(" ")[0][0] + CSL2.split(" ")[1][0]).lower().replace(" ", "")
            )
        else:
            if len(CSL.split(" ")) > 1:
                ID = CSL2[0:3].lower().replace(" ", "")
            else:
                ID = CSL2[0:2].lower().replace(" ", "")
        if not CSL_ID.get(CSL) and ID in ModMaster:
            r = re.compile(ID)
            ID = ID + str(len(list(filter(r.match, ModMaster))) + 1)
        LocVar[CSL] = ""
        LocVar[CSL] += ID
        listOfID.append(ID)
        LipidMasterlistOfID.append(ID)

    with open(os.path.join(project_root, "files", "pre_sbml_pos_comp.pk"), "wb") as f:
        dill.dump(
            {
                "name": ModName,
                "id": ModID,
                "mets": MetList,
                "mets_cl": MetList_CL,
                "met_equiv": MetEquiv,
                "reactions": RxnList,
                "reactions_cl": RxnList_CL,
                "genes": GeneList,
                "pathways": PathNameRxn,
                "loc": LocVar,
            },
            f,
        )

    with open(os.path.join(project_root, "files", "pre_sbml_pos_comp.pk"), "rb") as f:
        # This file was serialized with dill.dump earlier in this script.
        # Use dill.load to correctly deserialize objects (and avoid
        # Python 2 -> 3 module name issues such as '__builtin__').
        data = dill.load(f)
    try:
        LOGGER.debug("Loaded pre_sbml_pos_comp.pk keys: %s", list(data.keys()))
        if isinstance(data, dict):
            for k in ("mets", "reactions", "genes", "pathways", "loc"):
                if k in data:
                    v = data[k]
                    try:
                        LOGGER.debug("%s: type=%s, len=%s", k, type(v), len(v))
                    except Exception:
                        LOGGER.debug("%s: type=%s", k, type(v))
        # Sanitize reactions that may have been loaded from older pickles
        try:
            sanitize_loaded_reactions(RxnList, name="RxnList")
            sanitize_loaded_reactions(RxnList_CL, name="RxnList_CL")
        except Exception:
            LOGGER.debug(
                "Failed to sanitize reactions after loading pre_sbml_pos_comp.pk",
                exc_info=True,
            )
    except Exception:
        LOGGER.debug("Could not introspect loaded pre_sbml_pos_comp.pk", exc_info=True)
    # with open('files/pre_sbml_pos_comp.pk', 'rb') as f:
    #  data = f.read()

    model = cobra_reconstruction(
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
    cobra.io.write_sbml_model(model, Output)

    # Remove checkpoint file after successful completion
    if os.path.exists(checkpoint_file):
        try:
            os.remove(checkpoint_file)
            LOGGER.info("Checkpoint file removed after successful completion")
            print("Database generation completed successfully!")
        except Exception as e:
            LOGGER.warning(f"Could not remove checkpoint file: {e}")
