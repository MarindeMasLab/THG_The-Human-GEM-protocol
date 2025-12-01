#!/usr/bin/python
"""
Generate metabolic model database from KEGG pathways.

Note: Lambda functions have been completely removed from this code to avoid
pickle serialization issues. Reaction methods use stable module-level functions
bound via MethodType, and list concatenation uses operator.add instead of lambdas.

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
import operator
import os
import pickle
import re
import requests
import urllib
from functools import reduce
from typing import Dict, List

import cobra
from cobra.io import read_sbml_model
import dill
import sys
import pdb
from dotenv import load_dotenv
from tqdm import tqdm

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
from functions.error_tracker import ErrorTracker
from types import MethodType
from functions.ensembl_client import fetch_ensembl_annotations

# Debugging flag: limit number of reactions to process (None = no limit)
# Set this to an integer to process at most that many reactions and then
# stop early. Useful when debugging to avoid long runs.
MAX_REACTIONS = None  # e.g. set to 100 for debugging


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
    
    This function provides backward compatibility for checkpoint files created
    before the lambda removal fix. It ensures each reaction has `subs`/`prods` 
    attributes and that callable accessors (`Substrate`, `Product`, etc.) are
    bound to stable module-level methods.
    
    Note: New checkpoints created with this fixed code won't need sanitization,
    but this function remains for loading legacy checkpoint files.
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
            # Ensure lists exist - convert None or empty string to empty list
            # Empty string can occur when reaction initialization fails (e.g., invalid URL)
            if subs is None or subs == "":
                # Check if there's a subs attribute we should use instead
                subs = getattr(rxn, "subs", [])
            if prods is None or prods == "":
                # Check if there's a prods attribute we should use instead
                prods = getattr(rxn, "prods", [])

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
    
    # First, collect all metabolites actually used in reactions
    LOGGER.info("Scanning reactions to find metabolites actually used")
    metabolites_in_reactions = set()
    for rxn_id, rxn in reaction_list.items():
        try:
            # Get reaction compartment from reaction ID (e.g., "R00703_cytosol" -> "cytosol")
            rxn_compartment = rxn_id.split("_", 1)[1] if "_" in rxn_id else ""
            
            # Get metabolites from substrates and products
            if hasattr(rxn, 'Substrate') and hasattr(rxn, 'Product'):
                substrates = [subs[2] for subs in rxn.Substrate()]
                products = [prod[2] for prod in rxn.Product()]
            else:
                # Fallback to subs/prods attributes
                substrates = [subs[2] for subs in rxn.subs] if hasattr(rxn, 'subs') else []
                products = [prod[2] for prod in rxn.prods] if hasattr(rxn, 'prods') else []
            
            # Build full metabolite IDs with compartment (to match metabolite_list keys)
            for met_id in substrates + products:
                full_met_id = f"{met_id}_{rxn_compartment}"
                metabolites_in_reactions.add(full_met_id)
        except Exception as e:
            LOGGER.warning(f"Could not extract metabolites from reaction {rxn_id}: {e}")
    
    LOGGER.info(f"Found {len(metabolites_in_reactions)} metabolites used in {len(reaction_list)} reactions")
    
    LOGGER.info(f"Found {len(metabolites_in_reactions)} unique metabolite IDs referenced in reactions")
    
    # metabolites - only add those that are actually used in reactions
    LOGGER.info(f"Creating metabolites from {len(metabolite_list)} total metabolite entries")
    # Create metabolites using explicit keyword arguments to avoid positional
    # argument ordering mistakes (name/formula/charge were previously swapped).
    # Use a dict to deduplicate by ID2+compartment to avoid duplicate metabolites
    # when equivalent metabolites (e.g., D00003 and C00007) have the same ID2
    # Only include metabolites that are referenced in reactions (prevent orphans)
    unique_metabolites = {}
    filtered_out = 0
    for iden, compound in metabolite_list.items():
        # Only add metabolites that are actually used in reactions
        if iden not in metabolites_in_reactions:
            filtered_out += 1
            continue
            
        met_id = compound.ID2 + "_" + location_dict.get(compound.Subcel.lower())
        if met_id not in unique_metabolites:
            unique_metabolites[met_id] = cobra.Metabolite(
                id=met_id,
                name=compound.Name,
                formula=compound.Formula1,
                charge=(
                    float(compound.charge)
                    if compound.charge is not None
                    and str(compound.charge) not in ["", "None"]
                    else None
                ),
                compartment=location_dict.get(compound.Subcel.lower()),
            )
    compounds = list(unique_metabolites.values())
    LOGGER.info(f"Created {len(compounds)} unique metabolites (filtered out {filtered_out} not used in reactions)")
    
    model.add_metabolites(compounds)
    LOGGER.info(f"Successfully added {len(model.metabolites)} metabolites to model")
    
    # Set compartment names (COBRA only sets IDs by default)
    # location_dict maps: compartment_name -> compartment_id
    # model.compartments should be: {compartment_id: compartment_name}
    # Create reverse mapping: compartment_id -> compartment_name
    compartment_id_to_name = {v: k for k, v in location_dict.items()}
    
    # IMPORTANT: model.compartments is a read-only property!
    # We need to set the private _compartments attribute instead
    model._compartments = compartment_id_to_name
    
    LOGGER.info(f"Set names for {len(model.compartments)} compartments")
    
    # store mapping (met identifier -> met.id in model) for reaction section
    met_mapping = {}
    # metabolite annotation - only annotate metabolites that are in the model
    LOGGER.info(f"Annotating metabolites (only those used in reactions)")
    annotated_count = 0
    skipped_count = 0
    for iden, compound in tqdm(
        metabolite_list.items(), desc="Annotating metabolites", unit="met"
    ):
        # Skip metabolites not in reactions (already filtered out above)
        if iden not in metabolites_in_reactions:
            skipped_count += 1
            continue
        
        met_id = compound.ID2 + "_" + location_dict.get(compound.Subcel.lower())
        try:
            model_met = model.metabolites.get_by_id(met_id)
        except KeyError:
            LOGGER.warning(f"Metabolite {met_id} not found in model (iden={iden})")
            continue
            
        met_mapping[iden] = model_met.id
        annotation = {
            "pubchem.compound": compound.PubChem,
            "chebi.compound": compound.CheBI,
            "glycomedb": compound.GlyDB,
            "jcggdb": compound.JCGGDB,
            "inchi": compound.inchi,
            "inchikey": compound.inchikey,
            "lipidbank": compound.LipidBank,
            "lipidmaps": compound.LIPIDMAPS,
        }
        alt_formulas = [compound.Formula2, compound.Formula3, compound.Formula4]
        alt_formula = [
            form for form in alt_formulas if form != model_met.formula and form
        ]
        if alt_formula:
            model_met.annotation["glycan_formula"] = alt_formula[0]

        # only add non-empty annotation
        model_met.annotation = {k: v for k, v in annotation.items() if v}
        annotated_count += 1
    
    LOGGER.info(f"Annotated {annotated_count} metabolites, skipped {skipped_count} not used in reactions")

    LOGGER.info("Processing glycan formulas")
    for x in tqdm(
        model.metabolites, desc="Processing glycans", unit="met"
    ):  # replace glycan formula by a sbml suitable format
        if x.id[0] == "G":
            try:
                LOGGER.debug(f"Processing glycan: {x.id}")
                compartment = [
                    y[0] for y in location_dict.items() if y[1] in x.id.split("_")[1]
                ][0]
                xth_metabolite_id = x.id.split("_")[0] + "_" + compartment
                x.formula = metabolite_list[xth_metabolite_id].Formula4
            except Exception as e:
                import traceback
                LOGGER.error(f"Error updating glycan formula for {x.id}: {e}")
                LOGGER.error(traceback.format_exc())
                continue
    LOGGER.info("Normalizing metabolite IDs using equivalency mapping")
    for x in tqdm(
        model.metabolites, desc="Normalizing metabolite IDs", unit="met"
    ):  # eliminate potential discrepancies between metabolite id and reaction compounds ids
        if x.id.split("_")[0] in metabolite_equivalent.keys():
            equiv_id = metabolite_equivalent[x.id.split("_")[0]]
            if equiv_id not in metabolite_list_general:
                LOGGER.warning(f"Metabolite ID '{equiv_id}' from metabolite_equivalent not found in metabolite_list_general for '{x.id}'")
                continue
            new_id = (
                metabolite_list_general[equiv_id].ID1
                + "_"
                + x.compartment.lower()
            )
            # Skip renaming if target ID already exists (canonical ID already in model)
            if new_id in model.metabolites:
                LOGGER.debug(f"Skipping normalization of '{x.id}' to '{new_id}' - target already exists")
                continue
            x.id = new_id            
    def normalize_id(reac_id: str):
        iden, comp_desc = reac_id.split("_")
        comp = location_dict[comp_desc.lower()]
        return f"{iden}_{comp}"

    from functions.gpr.ast_gpr import sanitize_gpr

    # reactions
    LOGGER.info(f"Adding {len(reaction_list)} reactions to the model")
    reactions = [
        cobra.Reaction(
            normalize_id(iden), reac.Name(), "", 0 if reac.Termodyn() else -1000, 1000
        )
        for iden, reac in reaction_list.items()
    ]
    model.add_reactions(reactions)
    LOGGER.info(f"Successfully added {len(model.reactions)} reactions")

    LOGGER.info("Processing reaction metabolites and annotations")
    for iden, rxn in tqdm(
        reaction_list.items(), desc="Processing reactions", unit="rxn"
    ):
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
                    # LOGGER.warning(
                    #     f"Reactant '{met_root}' was added in a different comparment '{met_comp}'."
                    # )
                    # Use keyword args to ensure fields are assigned correctly
                    new_met = cobra.Metabolite(
                        id=met_id,
                        name=model_met.name,
                        formula=model_met.formula,
                        charge=model_met.charge,
                        # might be a new compartment!
                        compartment=met_comp,
                    )
                    new_met.annotation = model_met.annotation
                else:
                    # LOGGER.warning(
                    #     f"Reactant {met_id} was not found in any compartment. Creating new one!"
                    # )
                    # Create a minimal metabolite with explicit compartment
                    new_met = cobra.Metabolite(
                        id=met_id,
                        compartment=met_comp,
                    )
                new_mets.append(new_met)
        if new_mets:
            model.add_metabolites(new_mets)

        reac.add_metabolites(metabolites)
        ec = rxn.EC()
        reac.annotation = {"kegg.reaction": kegg_id, "ec-code": ec[0] if ec else ""}
        
        # Handle reactions with empty or missing GPR data
        if rxn.GPR and len(rxn.GPR) >= 2:
            sgpr, gpr = rxn.GPR[0], rxn.GPR[1].replace("[", "").replace("]", "")
        else:
            sgpr, gpr = "", ""
            LOGGER.debug(f"Reaction {kegg_id} has no GPR data")
        
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
    LOGGER.info(f"Adding {len(pathways)} pathway groups")
    model.add_groups([cobra.core.Group(group, group) for group in pathways])
    LOGGER.info("Assigning reactions to pathway groups")
    for group, members in tqdm(
        pathways.items(), desc="Assigning pathways", unit="pathway"
    ):
        # the members are the reactions in each pathway
        # TODO(carrascomj): reaction ids coming from paths are not in compartments
        model.groups.get_by_id(group).add_members(
            reduce(
                operator.add,
                [model.reactions.query(member) for member in members.split()],
                [],
            )
        )
    # gene annotation (genes were added with the GPRs)
    LOGGER.info(f"Processing {len(gene_list)} genes")
    pat_enstp = re.compile("ENS[TP][0-9]+")
    # Batch-fetch Ensembl annotations using gene symbols/names instead of calling Ensg()
    # The gene names in gene_list are the gene symbols (e.g., "BRCA1", "TP53")
    LOGGER.info("Batch-fetching Ensembl annotations for gene symbols")
    gene_symbols = list(gene_list.keys())
    LOGGER.info(f"Attempting to fetch annotations for {len(gene_symbols)} gene symbols")

    # Try to fetch using gene symbols - the API can look up by symbol
    try:
        ensembl_annotations = fetch_ensembl_annotations(gene_symbols, max_workers=10)
        LOGGER.info(
            "Fetched Ensembl annotations for %d genes (out of %d requested)",
            len(ensembl_annotations),
            len(gene_symbols),
        )
    except Exception as e:
        LOGGER.warning("Batch Ensembl annotation fetch failed: %s", e)
        ensembl_annotations = {}

    LOGGER.info(f"Annotating {len(gene_list)} genes")

    for iden, gene in tqdm(gene_list.items(), desc="Annotating genes", unit="gene"):
        # Avoid printing every gene to stdout (very slow for large models).
        # Use debug logging so the output can be enabled when needed.
        # LOGGER.debug("Processing gene: %s", iden)

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

        # Try to use batched Ensembl annotations (gene symbol as key)
        ensembl_genes = []
        entrez_ids = []
        uniprots = []

        # Check if we have this gene symbol in our batched results
        if iden in ensembl_annotations:
            ann = ensembl_annotations[iden]
            if ann.get("ensembl"):
                ensembl_genes = [ann.get("ensembl")]
            if ann.get("entrez"):
                entrez_ids = ann.get("entrez")
            if ann.get("uniprot"):
                uniprots = ann.get("uniprot")

        # Fall back to local file lookups (Entrez/Uniprot from db file)
        # These are fast - just reading from local files
        if not entrez_ids:
            try:
                e = gene.Entrez()  # Fast - reads from local file
                if e:
                    entrez_ids = [e]
            except Exception:
                pass

        if not uniprots:
            try:
                u = gene.Uniprot()  # Fast - reads from local file
                if u:
                    uniprots = [u]
            except Exception:
                pass

        model_gene.annotation = {
            k: v
            for k, v in {
                "ensembl": ensembl_genes if ensembl_genes else None,
                "ncbigene": entrez_ids if entrez_ids else None,
                "uniprot": uniprots if uniprots else None,
                "hgcn.symbol": model_gene.name,
            }.items()
            if v
        }
    return model


def reaction_has_human_genes(rxn, time=20):
    """Check if reaction has any human genes via KEGG EC lookup.
    
    Args:
        rxn: Reaction object with EC() method
        time: Timeout for HTTP requests
        
    Returns:
        bool: True if at least one EC number has human genes, False otherwise
    """
    import requests
    
    ec_numbers = rxn.EC()
    if not ec_numbers:
        LOGGER.warning(f"Reaction {rxn.ID} has no EC numbers - including by default")
        return True  # Include reactions without EC (may have manual annotations)
    
    # Deduplicate EC numbers
    unique_ecs = set(ec_numbers)
    LOGGER.debug(f"Checking {len(unique_ecs)} unique EC numbers for human genes: {unique_ecs}")
    
    for ec in unique_ecs:
        try:
            url = f"https://rest.kegg.jp/link/hsa/ec:{ec}"
            response = requests.get(url, timeout=time)
            
            if response.status_code == 200 and response.text.strip():
                # Found at least one human gene for this EC
                gene_count = len(response.text.strip().split('\n'))
                LOGGER.debug(f"EC {ec} has {gene_count} human genes")
                return True
            else:
                LOGGER.debug(f"EC {ec} has no human genes")
        except Exception as e:
            LOGGER.warning(f"Error checking EC {ec} for human genes: {e}")
            # On error, include reaction to be safe
            return True
    
    LOGGER.info(f"Reaction {rxn.ID} excluded: no human genes found for any EC number")
    return False


def process_reaction_gpr(rxn, gpr_list, gpr_ident, session):
    """Process GPR data for a reaction and return GPR/Subcel data.
    
    Returns:
        tuple: (tmpGPR, tmpSC2, error_type) where:
            - tmpSC2 is [dict, dict] with compartment info
            - error_type is None (success), 'timeout', 'connection', or 'no_data'
    """
    tmpGPR = ()
    tmpSC = ()
    error_type = None
    has_timeout_error = False
    has_connection_error = False
    
    # Fetch GPR for each unique EC number (deduplicate to avoid redundant processing)
    for ec in set(rxn.EC()):
        if ec not in gpr_ident:
            gpr_ident.append(ec)
            try:
                gpr_list[ec] = gpr(ec, session)
            except (requests.exceptions.Timeout, TimeoutError) as e:
                LOGGER.warning(f"Timeout querying GPR for EC {ec}: {e}")
                has_timeout_error = True
                continue
            except (requests.exceptions.ConnectionError, requests.exceptions.RequestException) as e:
                LOGGER.warning(f"Connection error querying GPR for EC {ec}: {e}")
                has_connection_error = True
                continue
            except Exception as e:
                LOGGER.warning(f"Unknown error querying GPR for EC {ec}: {e}")
                continue
        
        try:
            gpr_result = gpr_list[ec].GprSubcell()
            # Check if we got a valid tuple result (not empty string)
            if (
                gpr_result
                and isinstance(gpr_result, tuple)
                and len(gpr_result) >= 4
            ):
                tmpGPR = tmpGPR + gpr_result[0:2]
                tmpSC = tmpSC + gpr_result[2:4]
        except (requests.exceptions.Timeout, TimeoutError) as e:
            LOGGER.warning(f"Timeout processing GPR for EC {ec}: {e}")
            has_timeout_error = True
            continue
        except (requests.exceptions.ConnectionError, requests.exceptions.RequestException) as e:
            LOGGER.warning(f"Connection error processing GPR for EC {ec}: {e}")
            has_connection_error = True
            continue
        except Exception as e:
            LOGGER.warning(f"Could not process GPR for EC {ec}: {e}")
            continue
    
    # Determine error type
    if has_timeout_error:
        error_type = 'timeout'
    elif has_connection_error:
        error_type = 'connection'
    
    # Reorganize S-GPRs and GPRs based on their specific location
    tmpSC2 = [dict(), dict()]
    
    # Debug: Log tmpSC structure
    LOGGER.debug(f"tmpSC has {len(tmpSC)} dicts")
    for idx, sc_dict in enumerate(tmpSC):
        if sc_dict:
            LOGGER.debug(f"  tmpSC[{idx}]: {list(sc_dict.keys())}")
    
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
            if not re.findall(r"^\[\]$", xth_tmp_gpr[y + y][x]):
                tmp_xth_sgpr += xth_tmp_gpr[y + y][x]
            if not re.findall(r"^\[\]$", xth_tmp_gpr[y + y + 1][x]):
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
    
    # Debug: log compartment information
    if tmpSC2[0]:
        LOGGER.debug(f"Compartments detected: {list(tmpSC2[0].keys())}")
    else:
        LOGGER.debug("No compartments detected in GPR data")
        # If no compartments and no timeout/connection errors, it's just no data available
        if error_type is None:
            error_type = 'no_data'
    
    return tmpGPR, tmpSC2, error_type


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
    # If PATHWAY_SUBSET is set, use that file instead (allows testing subset runs)
    subset_file = os.environ.get("PATHWAY_SUBSET")
    if subset_file:
        ListOfPaths = subset_file
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

    # Initialize error tracking system
    error_tracker = ErrorTracker()
    LOGGER.info("Initialized comprehensive error tracking system")

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
    gene_cache_file = os.path.join(project_root, "files", "gene_location_cache.pkl")

    # Try to load gene location cache
    from functions.class_generate_database import gene as GeneClass
    GeneClass.load_cache(gene_cache_file)

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
                # Restore the set of reactions with failed location detection
                failed_location_reactions = checkpoint_data.get("failed_location_reactions", set())
                if failed_location_reactions:
                    LOGGER.info(
                        f"Restored {len(failed_location_reactions)} reactions with previously failed location detection: "
                        f"{', '.join(sorted(failed_location_reactions))}"
                    )
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

    # Internal counters/flags used when MAX_REACTIONS is set
    processed_reactions = 0
    stop_processing = False

    # Track failed reactions for retry
    failed_reactions = []  # List of tuples: (RxnID, RxnURL, PathName, RxnTermDyn, error_msg)
    
    # Track reactions with failed location detection (should not be saved to checkpoint)
    failed_location_reactions = set()  # Set of RxnIDs where location detection failed
    
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

                        # Check if reaction has human genes before processing
                        if not reaction_has_human_genes(RxnList[RxnID], time):
                            LOGGER.info(f"Skipping {RxnID}: no human genes found")
                            # Remove from RxnIdent and RxnList
                            RxnIdent.remove(RxnID)
                            del RxnList[RxnID]
                            j = j + 1
                            continue

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
                                    LOGGER.debug(f"  - ID1: {MetList[CompID].ID1}")
                                    LOGGER.debug(f"  - ID2: {MetList[CompID].ID2}")
                                    LOGGER.debug(f"  - Name: {MetList[CompID].Name}")
                                    LOGGER.debug(
                                        f"  - Formula1: {MetList[CompID].Formula1}"
                                    )
                                    LOGGER.debug(
                                        f"  - Formula2: {MetList[CompID].Formula2}"
                                    )
                                    if (
                                        MetList[CompID].ID1
                                        and MetList[CompID].ID1[0] == "G"
                                    ):
                                        LOGGER.debug(
                                            f"  - [GLYCAN] Formula4 (reformulated): {MetList[CompID].Formula4}"
                                        )
                                    LOGGER.debug(
                                        f"  - Atom composition: {MetList[CompID].Atom1}"
                                    )

                                    # Handle ID equivalences
                                    if MetList[CompID].ID1 != MetList[CompID].ID2:
                                        LOGGER.debug(
                                            f"ID equivalence found: {CompID} -> {MetList[CompID].ID1}"
                                        )
                                        MetIdent[len(MetIdent) - 1] = MetList[
                                            CompID
                                        ].ID1
                                        MetEquiv[CompID] = MetList[CompID].ID1
                                        # Direct assignment instead of deepcopy when possible
                                        MetList[MetList[CompID].ID1] = MetList[CompID]
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
                        IthRxnMB = None
                        if mb_test != 0:
                            try:
                                IthRxnMB = mass_balance(eq, RxnID)
                                LOGGER.debug(f"Mass balance result: {IthRxnMB}")
                            except Exception as e:
                                LOGGER.warning(f"Mass balance failed for reaction {RxnID}: {e}")
                                error_tracker.add_error(
                                    'mass_balance',
                                    f"Mass balance failed: {str(e)}",
                                    context={'reaction_id': RxnID, 'equation': eq},
                                    level='warning'
                                )
                                IthRxnMB = None
                        
                        # Check if mass balance was successful
                        # IthRxnMB[10] is TestOfBalance: 0=balanced initially, 1=balanced after adjustment, 2=failed
                        # IthRxnMB[4] is NewSpecies[0]: list of new compounds to add (can be empty for already-balanced reactions)
                        mass_balance_succeeded = IthRxnMB and IthRxnMB[10] != 2
                        
                        if mass_balance_succeeded:
                            # If new compounds have to be added to mass balance the reactions,
                            # check if they need to be added to the network as compounds
                            if IthRxnMB[4]:
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
                            error_tracker.add_success('mass_balance')
                        else:  # if the reaction cannot be mass balanced all the stoichimetric coef are assumed to be like in the original reaction
                            LOGGER.warning(
                                f"Reaction {RxnID} cannot be mass balanced - using original stoichiometry"
                            )
                            error_tracker.add_error(
                                'mass_balance',
                                f"Cannot mass balance reaction - using original stoichiometry",
                                context={'reaction_id': RxnID, 'equation': eq},
                                level='warning'
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
                                MetList[met_id].Name if met_id in MetList else "Unknown"
                            )
                            LOGGER.debug(f"    {sub[0]} {met_id} ({met_name})")
                        LOGGER.debug(f"  Products:")
                        for prod in P:
                            met_id = prod[2]
                            met_name = (
                                MetList[met_id].Name if met_id in MetList else "Unknown"
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
                        # Note: subs and prods were already set above (lines 1047-1048)
                        if not RxnID in PathNameRxn.get(PathName):
                            PathNameRxn[PathName] += RxnID + " "

                        ######### Define New GPR ###########
                        tmpGPR, tmpSC2, error_type = process_reaction_gpr(
                            RxnList[RxnID], GPRList, GPRIdent, session
                        )
                        
                        RxnList[RxnID].GPR = tmpGPR
                        RxnList[RxnID].Subcel = tmpSC2

                        # Determine if location detection failed and how to handle it
                        # Following the 4 scenarios:
                        # 1. No human genes → already handled (skipped earlier)
                        # 2. Compartment found but not in Excel (restricted mode) → place in cytosol
                        # 3. Human reaction but no EC/location in databases → place in cytosol
                        # 4. Timeout/connection errors → add to retry queue, skip rxnSubcel
                        
                        should_skip_reaction = False
                        location_detection_failed = False
                        
                        if error_type in ('timeout', 'connection'):
                            # Case 4: Timeout or connection error → add to retry queue, skip processing
                            failed_location_reactions.add(RxnID)
                            should_skip_reaction = True
                            LOGGER.warning(
                                f"Reaction {RxnID} added to retry queue due to {error_type} error. "
                                "Will NOT be compartmentalized until location data is available."
                            )
                        elif not tmpSC2 or (isinstance(tmpSC2, (list, tuple)) and len(tmpSC2) >= 1 and 
                                          (not tmpSC2[0] or all(k == '' or not k for k in tmpSC2[0].keys()))):
                            # Case 2 & 3: No location data but valid human genes → default to cytosol
                            location_detection_failed = True
                            LOGGER.info(
                                f"No location data for reaction {RxnID} (error_type: {error_type}), will default to cytosol"
                            )

                        # Only process through rxnSubcel if NOT in retry queue
                        if not should_skip_reaction:
                            logging.debug(
                                f"len metlist_CL before rxnSubcel: {len(MetList_CL)}"
                            )

                            # TODO: rxnSubcel should also return MetList_CL - now it
                            # updates it in place
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
                            logging.debug(
                                f"len metlist_CL after rxnSubcel: {len(MetList_CL)}"
                            )

                            RxnIdent_CL.extend(
                                rxn_cl
                            )  # Optimized: use extend instead of concatenation
                            MetIdent_CL.extend(
                                comp_cl
                            )  # Optimized: use extend instead of concatenation
                        else:
                            # Reaction in retry queue - skip compartmentalization
                            LOGGER.info(
                                f"Skipping compartmentalization for reaction {RxnID} (in retry queue)"
                            )

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
                        # Increment processed reactions and check debug limit
                        try:
                            processed_reactions += 1
                        except NameError:
                            # If variables not present for some reason, initialize
                            processed_reactions = 1

                        if (
                            MAX_REACTIONS is not None
                            and processed_reactions >= MAX_REACTIONS
                        ):
                            LOGGER.info(
                                "Reached MAX_REACTIONS=%s. Stopping early for debugging.",
                                MAX_REACTIONS,
                            )
                            stop_processing = True
                            break
                except Exception as e:
                    # Extract RxnID safely
                    rxn_id = 'unknown'
                    rxn_url = None
                    rxn_termdyn = None
                    try:
                        rxn_id = PathList[PathID].Reactions()[j][0][0]
                        rxn_url = PathList[PathID].Reactions()[j][1]
                        rxn_termdyn = PathList[PathID].Reactions()[j][0][1]
                    except:
                        pass
                    
                    error_msg = str(e)
                    LOGGER.error(
                        f"Failed to process reaction {rxn_id} from pathway {PathName}: {error_msg}"
                    )
                    import traceback
                    LOGGER.error(traceback.format_exc())
                    
                    # Add to failed reactions list for retry
                    failed_reactions.append({
                        'RxnID': rxn_id,
                        'RxnURL': rxn_url,
                        'PathName': PathName,
                        'RxnTermDyn': rxn_termdyn,
                        'error': error_msg,
                        'attempt': 1
                    })
                    LOGGER.info(f"Added {rxn_id} to retry queue (total failed: {len(failed_reactions)})")
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
            # Filter out reactions with failed location detection (timeout/connection errors)
            # These are in the retry queue and will be reprocessed
            RxnList_CL_filtered = {
                rxn_key: rxn_obj 
                for rxn_key, rxn_obj in RxnList_CL.items()
                if not any(failed_rxn_id in rxn_key for failed_rxn_id in failed_location_reactions)
            }
            
            # Also filter metabolites - only keep those used in filtered reactions
            if failed_location_reactions:
                metabolites_in_filtered_rxns = set()
                for rxn_key, rxn_obj in RxnList_CL_filtered.items():
                    if '_' in rxn_key:
                        comp = '_'.join(rxn_key.split('_')[1:])
                        # Extract substrates
                        if hasattr(rxn_obj, 'subs') and rxn_obj.subs:
                            for sub in rxn_obj.subs:
                                if len(sub) >= 3:
                                    met_id = f"{sub[2]}_{comp}"
                                    metabolites_in_filtered_rxns.add(met_id)
                        # Extract products
                        if hasattr(rxn_obj, 'prods') and rxn_obj.prods:
                            for prod in rxn_obj.prods:
                                if len(prod) >= 3:
                                    met_id = f"{prod[2]}_{comp}"
                                    metabolites_in_filtered_rxns.add(met_id)
                
                MetList_CL_filtered = {
                    met_key: met_obj
                    for met_key, met_obj in MetList_CL.items()
                    if met_key in metabolites_in_filtered_rxns
                }
                
                LOGGER.info(
                    f"Excluding {len(failed_location_reactions)} reactions (retry queue) and "
                    f"{len(MetList_CL) - len(MetList_CL_filtered)} associated metabolites from checkpoint: "
                    f"{', '.join(sorted(failed_location_reactions))}"
                )
            else:
                MetList_CL_filtered = MetList_CL
            
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
                "RxnList_CL": RxnList_CL_filtered,  # Filtered reactions
                "MetList_CL": MetList_CL_filtered,   # Filtered metabolites (no orphans!)
                "RxnIdent_CL": RxnIdent_CL,
                "MetIdent_CL": MetIdent_CL,
                "Compartment_CL": Compartment_CL,
                "failed_location_reactions": failed_location_reactions,  # Track for retry
            }
            with open(checkpoint_file, "wb") as f:
                dill.dump(checkpoint_data, f)
            
            # Save gene location cache
            GeneClass.save_cache(gene_cache_file)
            
            LOGGER.debug(f"Checkpoint saved after pathway {i + 1}/{len(Path) - 1}")
        except Exception as e:
            LOGGER.warning(f"Failed to save checkpoint: {e}")

        # If a debugging stop was requested (MAX_REACTIONS reached), break out
        # after saving the checkpoint so partial work is preserved.
        if stop_processing:
            LOGGER.info(
                "Stopping processing after reaching MAX_REACTIONS=%s", MAX_REACTIONS
            )
            break

        i = i + 1

    # Retry mechanism for failed reactions
    if failed_reactions:
        LOGGER.info("="*80)
        LOGGER.info("RETRY PHASE: %d reactions failed during first pass", len(failed_reactions))
        LOGGER.info("="*80)
        
        max_retry_attempts = 3
        permanently_failed = []
        retry_round = 1
        
        while failed_reactions and retry_round <= max_retry_attempts:
            LOGGER.info("Retry round %d: Attempting to process %d failed reactions", 
                       retry_round, len(failed_reactions))
            
            # Copy list for iteration, we'll modify the original
            reactions_to_retry = failed_reactions[:]
            failed_reactions.clear()
            
            for failed_rxn in reactions_to_retry:
                rxn_id = failed_rxn['RxnID']
                rxn_url = failed_rxn['RxnURL']
                path_name = failed_rxn['PathName']
                rxn_termdyn = failed_rxn['RxnTermDyn']
                
                LOGGER.info("Retrying reaction %s (attempt %d)", rxn_id, retry_round + 1)
                
                try:
                    # Create reaction object with correct parameters
                    RxnList[rxn_id] = reaction(rxn_url, time, rxn_id, path_name, rxn_termdyn)
                    
                    # Process reaction (similar to main loop logic)
                    # Check compounds and add if missing
                    RxnCmp = [x[2] for x in RxnList[rxn_id].Substrate()] + [
                        x[2] for x in RxnList[rxn_id].Product()
                    ]
                    
                    for cmp_id in RxnCmp:
                        if cmp_id not in MetList and cmp_id not in extra_compound:
                            cmp_url = "https://www.genome.jp/entry/" + cmp_id
                            MetList[cmp_id] = compound(cmp_url, cmp_id, time, EF, specialCompounds)
                    
                    # Bind stable methods (avoiding lambda pickle issues)
                    S2 = copy.deepcopy([x for x in RxnList[rxn_id].Substrate()])
                    P2 = copy.deepcopy([x for x in RxnList[rxn_id].Product()])
                    RxnList[rxn_id].subs = S2
                    RxnList[rxn_id].prods = P2
                    RxnList[rxn_id].Substrate = MethodType(_rxn_substrate, RxnList[rxn_id])
                    RxnList[rxn_id].Product = MethodType(_rxn_product, RxnList[rxn_id])
                    RxnList[rxn_id].SetSubstrate = MethodType(_rxn_substrate, RxnList[rxn_id])
                    RxnList[rxn_id].SetProduct = MethodType(_rxn_product, RxnList[rxn_id])
                    
                    if path_name not in PathNameRxn:
                        PathNameRxn[path_name] = ""
                    if rxn_id not in PathNameRxn.get(path_name):
                        PathNameRxn[path_name] += rxn_id + " "
                    
                    # Process GPR - this is critical for rxnSubcel to work!
                    tmpGPR, tmpSC2, error_type = process_reaction_gpr(
                        RxnList[rxn_id], GPRList, GPRIdent, session
                    )
                    
                    # Assign GPR data to reaction (required for rxnSubcel)
                    RxnList[rxn_id].GPR = tmpGPR
                    RxnList[rxn_id].Subcel = tmpSC2
                    
                    # Check if still has timeout/connection errors
                    if error_type in ('timeout', 'connection'):
                        LOGGER.warning(
                            f"Retry for reaction {rxn_id} still has {error_type} error, keeping in retry queue"
                        )
                        continue  # Skip this reaction, will retry again next time
                    
                    # Successfully retrieved location data - can now compartmentalize
                    LOGGER.info(f"Retry successful for reaction {rxn_id}, compartmentalizing now")
                    
                    # Expand reaction to compartments
                    Compartment_CL, rxn_cl, comp_cl = rxnSubcel(
                        RxnList[rxn_id],
                        RxnList_CL,
                        MetList_CL,
                        RxnIdent_CL,
                        MetIdent_CL,
                        Compartment_CL,
                        RxnList,
                        MetList,
                        MetEquiv,
                    )
                    
                    RxnIdent_CL.extend(rxn_cl)
                    MetIdent_CL.extend(comp_cl)
                    
                    LOGGER.info("Successfully processed %s on retry", rxn_id)
                    
                except Exception as e:
                    # Still failed, track for next retry or permanent failure
                    error_msg = str(e)
                    LOGGER.warning("Reaction %s failed again on retry round %d: %s", 
                                 rxn_id, retry_round + 1, error_msg)
                    
                    failed_rxn['attempt'] = retry_round + 1
                    failed_rxn['error'] = error_msg
                    failed_reactions.append(failed_rxn)
            
            retry_round += 1
        
        # Any remaining failures are permanent
        if failed_reactions:
            permanently_failed = failed_reactions[:]
            LOGGER.error("="*80)
            LOGGER.error("PERMANENTLY FAILED REACTIONS: %d", len(permanently_failed))
            LOGGER.error("="*80)
            for failed_rxn in permanently_failed:
                LOGGER.error("Reaction %s from pathway %s failed after %d attempts. Last error: %s",
                           failed_rxn['RxnID'], failed_rxn['PathName'], 
                           failed_rxn['attempt'], failed_rxn['error'])
        else:
            LOGGER.info("="*80)
            LOGGER.info("All failed reactions successfully recovered during retry phase!")
            LOGGER.info("="*80)

    # Retry mechanism for reactions with location detection timeouts/connection errors
    if failed_location_reactions:
        LOGGER.info("="*80)
        LOGGER.info("LOCATION RETRY PHASE: %d reactions had timeout/connection errors during gene location lookup", 
                   len(failed_location_reactions))
        LOGGER.info("Gene cache now contains ~%d genes - retrying with fuller cache for instant lookups", 
                   GeneClass.get_cache_stats()['cache_size'])
        LOGGER.info("="*80)
        
        retry_success = 0
        retry_failed = 0
        still_failed_location_reactions = set()
        
        # Convert set to list for iteration
        reactions_to_retry = list(failed_location_reactions)
        
        for rxn_id in reactions_to_retry:
            LOGGER.info("Retrying gene location lookup for reaction %s", rxn_id)
            
            try:
                # Find the reaction object in RxnList
                if rxn_id not in RxnList:
                    LOGGER.warning("Reaction %s not found in RxnList, cannot retry", rxn_id)
                    still_failed_location_reactions.add(rxn_id)
                    retry_failed += 1
                    continue
                
                rxn_obj = RxnList[rxn_id]
                
                # Re-attempt gene location lookup with fuller cache
                tmpGPR, tmpSC2, error_type = process_reaction_gpr(
                    rxn_obj, GPRList, GPRIdent, session
                )
                
                # Update reaction with new data
                rxn_obj.GPR = tmpGPR
                rxn_obj.Subcel = tmpSC2
                
                # Check if still has timeout/connection errors
                if error_type in ('timeout', 'connection'):
                    LOGGER.warning(
                        f"Retry for reaction {rxn_id} still has {error_type} error"
                    )
                    still_failed_location_reactions.add(rxn_id)
                    retry_failed += 1
                    continue
                
                # Successfully retrieved location data - can now compartmentalize
                LOGGER.info(f"Gene location lookup successful for reaction {rxn_id}, compartmentalizing now")
                
                # Expand reaction to compartments using rxnSubcel
                Compartment_CL, rxn_cl, comp_cl = rxnSubcel(
                    rxn_obj,
                    RxnList_CL,
                    MetList_CL,
                    RxnIdent_CL,
                    MetIdent_CL,
                    Compartment_CL,
                    RxnList,
                    MetList,
                    MetEquiv,
                )
                
                RxnIdent_CL.extend(rxn_cl)
                MetIdent_CL.extend(comp_cl)
                
                LOGGER.info("Successfully compartmentalized reaction %s on retry", rxn_id)
                retry_success += 1
                
            except Exception as e:
                error_msg = str(e)
                LOGGER.error("Reaction %s failed during location retry: %s", rxn_id, error_msg)
                import traceback
                LOGGER.error(traceback.format_exc())
                still_failed_location_reactions.add(rxn_id)
                retry_failed += 1
        
        # Update the failed_location_reactions set with only those that still failed
        failed_location_reactions = still_failed_location_reactions
        
        # Log results
        LOGGER.info("="*80)
        LOGGER.info("LOCATION RETRY RESULTS:")
        LOGGER.info(f"  Successfully recovered: {retry_success} reactions")
        LOGGER.info(f"  Still failed: {retry_failed} reactions")
        if failed_location_reactions:
            LOGGER.info(f"  Permanently failed reactions: {', '.join(sorted(failed_location_reactions))}")
        LOGGER.info("="*80)
        
        if retry_success > 0:
            print(f"\nLocation retry recovered {retry_success} reactions with fuller gene cache!")
        if failed_location_reactions:
            print(f"\nWarning: {len(failed_location_reactions)} reactions still failed after retry:")
            print(f"  {', '.join(sorted(failed_location_reactions))}")

    ######### Genes ###########
    GeneList = {}
    GeneIdent = []
    g = 0
    while g < len(GPRList):
        if GPRIdent[g] in GPRList.keys() and GPRList[GPRIdent[g]].GprSubcell():
            gene_matches = re.findall(
                r"([A-Za-z0-9\-]+)",
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
    
    # Log gene location cache statistics
    cache_stats = GeneClass.get_cache_stats()
    LOGGER.info("="*80)
    LOGGER.info("Gene Location Cache Statistics:")
    LOGGER.info(f"  Total cache lookups: {cache_stats['hits'] + cache_stats['misses']}")
    LOGGER.info(f"  Cache hits: {cache_stats['hits']}")
    LOGGER.info(f"  Cache misses: {cache_stats['misses']}")
    LOGGER.info(f"  Database queries: {cache_stats['queries']}")
    LOGGER.info(f"  Cache hit rate: {cache_stats['hit_rate']}")
    LOGGER.info(f"  Unique genes cached: {cache_stats['cache_size']}")
    LOGGER.info("="*80)
    print("\nGene Location Cache Statistics:")
    print(f"  Cache hit rate: {cache_stats['hit_rate']}")
    print(f"  Unique genes cached: {cache_stats['cache_size']}")
    print(f"  Database queries saved: {cache_stats['hits']}")
    
    # Generate and save comprehensive error report
    LOGGER.info("="*80)
    LOGGER.info("Generating comprehensive error report...")
    print("\nGenerating error report...")
    
    # Save error report to file
    error_report_file = os.path.join(project_root, "logs", "error_report.txt")
    error_tracker.save_report(error_report_file, verbose=True)
    
    # Export JSON version for programmatic analysis
    error_json_file = os.path.join(project_root, "logs", "error_report.json")
    error_tracker.export_json(error_json_file)
    
    # Print summary to console
    print("\n" + "="*80)
    print("ERROR REPORT SUMMARY")
    print("="*80)
    summary = error_tracker.generate_summary(verbose=False)
    # Print just the summary section (not full details)
    summary_lines = summary.split('\n')
    in_summary = False
    for line in summary_lines:
        if 'OVERALL STATISTICS' in line:
            in_summary = True
        if in_summary:
            print(line)
        if 'RECOMMENDATIONS' in line and in_summary:
            # Print recommendations section too
            for remaining_line in summary_lines[summary_lines.index(line):]:
                print(remaining_line)
            break
    
    print(f"\nFull error report saved to: {error_report_file}")
    print(f"JSON error data saved to: {error_json_file}")
    LOGGER.info(f"Error reports generated: {error_report_file}, {error_json_file}")
    LOGGER.info("="*80)

    # # Remove checkpoint file after successful completion
    # if os.path.exists(checkpoint_file):
    #     try:
    #         os.remove(checkpoint_file)
    #         LOGGER.info("Checkpoint file removed after successful completion")
    #         print("Database generation completed successfully!")
    #     except Exception as e:
    #         LOGGER.warning(f"Could not remove checkpoint file: {e}")
