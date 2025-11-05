#!/usr/bin/env python3
"""
Batch-prefetch build_model script.

This script demonstrates an approach to prefetch required remote data
(KEGG reaction entries, GPR lookups per EC, Ensembl annotations) in
larger batches, cache them locally, and then process reactions using
the cached information instead of calling the remote functions for
each reaction.

It is intended as a replacement/companion for the existing
`build_model.py` processing loop and focuses on the batching + cache
strategy. It re-uses repository helper functions where available.

Usage: run from repository root. Adjust INPUT_MODEL and other paths.
"""
import os
import re
import time
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
import pickle

import requests
import cobra
import pandas as pd
import copy
from collections import ChainMap
from cobra import Reaction, Metabolite
from tqdm import tqdm

# Repository helpers
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.join(current_dir, "..")
if project_root not in os.sys.path:
    os.sys.path.append(project_root)

from functions.gpr.gpr_def import getGPR, setup_biocyc_session
from functions.gpr.get_location_def import getLocationnew as getLocation
from functions.ensembl_client import fetch_ensembl_annotations
from functions.function_bm_gdb import (
    meltGeneList,
    create_compartments_dict_bm,
    create_comp_abbreviations_dict_bm,
    compartment_file_to_dict_bm,
    update_comp_names_bm,
)

LOGGER = logging.getLogger("build_model_batch")
logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s %(levelname)s %(name)s %(message)s"
)

# Config (tweak as needed)
INPUT_MODEL = os.path.join(project_root, "models", "THG-beta1.1.1_251031.xml")
OUTPUT_MODEL = os.path.join(project_root, "models", "THG-beta-batch.xml")
ENSEMBL_CACHE_FILE = os.path.join(project_root, "files", "ensembl_cache_batch.pkl")
KEGG_CACHE_FILE = os.path.join(project_root, "files", "kegg_reaction_cache_batch.pkl")
GETGPR_CACHE_FILE = os.path.join(project_root, "files", "getgpr_cache_batch.pkl")

KEGG_BATCH_SIZE = 10
ENSEMBL_BATCH_SIZE = 10
MAX_WORKERS_GPR = 1


def load_cache(path):
    """Load a pickle cache safely.

    - On unpickling error: rename the corrupt file to <path>.corrupt.<ts> and return empty dict.
    - If the loaded object is a dict, normalize keys to stripped strings to avoid
      mismatches between e.g. bytes/ints and expected string keys like 'R00010'.
    """
    if os.path.exists(path):
        try:
            with open(path, "rb") as f:
                obj = pickle.load(f)
            # If it's a mapping, coerce keys to normalized str keys
            if isinstance(obj, dict):
                normalized = {}
                for k, v in obj.items():
                    try:
                        nk = str(k).strip()
                    except Exception:
                        nk = repr(k)
                    normalized[nk] = v
                LOGGER.debug("Loaded cache %s (entries=%d)", path, len(normalized))
                return normalized
            else:
                LOGGER.debug("Loaded cache %s of type %s", path, type(obj))
                return obj
        except Exception:
            LOGGER.exception(
                "Failed to load cache %s - renaming corrupt file and continuing with empty cache",
                path,
            )
            try:
                corrupt = path + ".corrupt.%d" % int(time.time())
                os.replace(path, corrupt)
                LOGGER.warning("Renamed corrupt cache to %s", corrupt)
            except Exception:
                LOGGER.exception("Failed to rename corrupt cache %s", path)
    return {}


def save_cache(path, obj):
    """Save cache atomically. Writes to a temp file then replaces the target.

    This reduces risk of producing a partially-written pickle if the process is
    interrupted. Any exception is logged but not raised to avoid crashing the
    batch run.
    """
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        tmp = path + ".tmp"
        with open(tmp, "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
        # Atomic replace
        os.replace(tmp, path)
        LOGGER.debug(
            "Saved cache %s (entries=%d)",
            path,
            len(obj) if isinstance(obj, dict) else 1,
        )
    except Exception:
        LOGGER.exception("Failed to save cache %s", path)


def batch_fetch_kegg_entries(
    kegg_ids, batch_size=50, max_workers=3, requests_per_second=3
):
    """Fetch KEGG entries for a list of KEGG reaction ids using rest.kegg.jp/get

    Uses parallel requests with rate limiting to speed up fetching while respecting
    API limits.

    Args:
        kegg_ids: List of KEGG reaction IDs to fetch
        batch_size: Number of IDs to fetch in a single request (KEGG API supports multiple)
        max_workers: Maximum number of parallel requests
        requests_per_second: Maximum requests per second (rate limit)

    Returns:
        dict {id: text} mapping KEGG IDs to their entry text
    """
    base = "https://rest.kegg.jp/get/"
    results = {}
    ids = list(kegg_ids)

    # Create batches
    batches = []
    for i in range(0, len(ids), batch_size):
        batch = ids[i : i + batch_size]
        batches.append(batch)

    # Rate limiting: minimum time between requests
    min_delay = 1.0 / requests_per_second
    last_request_time = [0.0]  # Use list to allow modification in nested function
    request_lock = __import__("threading").Lock()

    def fetch_batch(batch):
        """Fetch a single batch with rate limiting."""
        query = "+".join(batch)
        url = base + query

        # Rate limiting
        with request_lock:
            elapsed = time.time() - last_request_time[0]
            if elapsed < min_delay:
                time.sleep(min_delay - elapsed)
            last_request_time[0] = time.time()

        try:
            LOGGER.debug("Fetching KEGG batch: %s", batch[:5])  # Show first 5 IDs
            r = requests.get(url, timeout=30)
            r.raise_for_status()
            text = r.text

            # rest.kegg returns concatenated entries; split by 'ENTRY' lines
            batch_results = {}
            entries = re.split(r"\n(?=ENTRY\s+R)", text)

            for entry in entries:
                if not entry.strip():
                    continue
                m = re.search(r"ENTRY\s+(R[0-9]+)", entry)
                if m:
                    rid = m.group(1)
                    batch_results[rid] = entry

            fetched_count = len(batch_results)
            LOGGER.debug("Fetched %d/%d entries in batch", fetched_count, len(batch))

            # Warn if we didn't get all expected entries
            if fetched_count < len(batch):
                missing = [b for b in batch if b not in batch_results]
                LOGGER.warning(
                    "Only fetched %d/%d KEGG entries in batch. Missing: %s",
                    fetched_count,
                    len(batch),
                    missing,
                )

            return batch_results

        except Exception as e:
            LOGGER.exception("Failed fetching KEGG batch: %s", e)
            return {}

    # Fetch batches in parallel with rate limiting
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(fetch_batch, batch): batch for batch in batches}

        for future in tqdm(
            as_completed(futures), total=len(futures), desc="KEGG batches"
        ):
            batch_results = future.result()
            results.update(batch_results)

    return results


def extract_kegg_reaction_id(reac):
    # Prefer annotation 'kegg.reaction' if present
    ann = reac.annotation if hasattr(reac, "annotation") else {}
    try:
        k = ann.get("kegg.reaction")
        if k:
            # annotation may be list-like
            return str(k).split()[0]
    except Exception:
        pass
    # fallback: try to find Rxxxx in reaction id or name
    m = re.search(r"R[0-9]{5}", reac.id)
    if m:
        return m.group(0)
    m = re.search(r"R[0-9]{5}", getattr(reac, "name", "") or "")
    if m:
        return m.group(0)
    return None


def main():
    LOGGER.info("Loading model %s", INPUT_MODEL)
    model = cobra.io.read_sbml_model(INPUT_MODEL)

    # prepare session for getGPR
    session = setup_biocyc_session()

    # Load caches
    ensembl_cache = load_cache(ENSEMBL_CACHE_FILE)
    kegg_cache = load_cache(KEGG_CACHE_FILE)
    getgpr_cache = load_cache(GETGPR_CACHE_FILE)

    # Collect unique ECs and KEGG reaction ids
    ec_set = set()
    kegg_ids = set()
    for r in model.reactions:
        ec = r.annotation.get("ec-code") if hasattr(r, "annotation") else None
        if isinstance(ec, (list, tuple)):
            for e in ec:
                if e:
                    ec_set.add(e)
        elif ec:
            ec_set.add(ec)
        kid = extract_kegg_reaction_id(r)
        if kid:
            kegg_ids.add(kid)

    LOGGER.info("Unique ECs: %d  KEGG reactions: %d", len(ec_set), len(kegg_ids))

    # Batch fetch KEGG entries and update cache
    missing_kegg = [k for k in kegg_ids if k not in kegg_cache]
    if missing_kegg:
        LOGGER.info("Fetching %d missing KEGG entries in batches", len(missing_kegg))
        fetched = batch_fetch_kegg_entries(missing_kegg, batch_size=KEGG_BATCH_SIZE)
        LOGGER.info(
            "Successfully fetched %d/%d KEGG entries", len(fetched), len(missing_kegg)
        )

        # Mark entries that don't exist in KEGG with None to avoid refetching
        for kid in missing_kegg:
            if kid not in fetched:
                LOGGER.warning("KEGG ID %s not found in KEGG database", kid)
                kegg_cache[kid] = None  # Mark as checked but not found
            else:
                kegg_cache[kid] = fetched[kid]

        save_cache(KEGG_CACHE_FILE, kegg_cache)
        LOGGER.info("Updated cache now has %d entries", len(kegg_cache))
    else:
        LOGGER.info("All KEGG entries present in cache")

    # Save all not-found KEGG reactions to CSV (from both cache and newly fetched)
    not_found_kegg = [kid for kid in kegg_ids if kegg_cache.get(kid) is None]
    if not_found_kegg:
        csv_path = os.path.join(project_root, "files", "kegg_reactions_not_found.csv")
        try:
            df = pd.DataFrame({"KEGG_ID": sorted(not_found_kegg)})
            df.to_csv(csv_path, index=False)
            LOGGER.info(
                "Saved %d not-found KEGG reactions to %s",
                len(not_found_kegg),
                csv_path,
            )
        except Exception:
            LOGGER.exception("Failed to save not-found KEGG reactions to CSV")
    else:
        LOGGER.info("All requested KEGG reactions were found in database")

    # For ECs, call getGPR once and cache per EC
    missing_ec = [e for e in ec_set if e and e not in getgpr_cache]
    LOGGER.info("getGPR: %d ECs missing from cache", len(missing_ec))

    def call_getGPR(ec_value, session_obj=None, verbose=True):
        """Wrapper around getGPR that supports a verbose flag.

        If verbose is False, temporarily silence the gpr_def logger for quieter output.
        Returns (ec, result) where result is the value returned by getGPR or None on failure.
        """
        logger_name = "functions.gpr.gpr_def"
        gpr_logger = logging.getLogger(logger_name)
        # prev_disabled = gpr_logger.disabled
        # prev_level = gpr_logger.level
        try:
            if not verbose:
                gpr_logger.disabled = True
            # Inform user (thread-safe with tqdm) that this EC has started processing
            try:
                tqdm.write(f"getGPR START: {ec_value}")
            except Exception:
                LOGGER.debug("getGPR START: %s", ec_value)
            # ensure we have a session per call if none provided
            if session_obj is None:
                session_local = setup_biocyc_session()
            else:
                session_local = session_obj
            res = getGPR(ec_value, session_local)
            # Short, safe summary for the result to avoid cluttering output
            try:
                summary = "None"
                if res is None:
                    summary = "FAILED"
                elif isinstance(res, (list, tuple)) and len(res) > 1:
                    names = res[1]
                    if isinstance(names, (list, tuple)):
                        summary = f"genes={len(names)} sample={names[:5]}"
                    else:
                        summary = f"genes=1 sample={names}"
                else:
                    summary = str(type(res))
                try:
                    tqdm.write(f"getGPR DONE : {ec_value} -> {summary}")
                except Exception:
                    LOGGER.info("getGPR DONE: %s -> %s", ec_value, summary)
            except Exception:
                # Fallback minimal message
                try:
                    tqdm.write(f"getGPR DONE : {ec_value}")
                except Exception:
                    LOGGER.info("getGPR DONE: %s", ec_value)
            return ec_value, res
        except Exception:
            LOGGER.exception("getGPR failed for %s", ec_value)
            return ec_value, None
        finally:
            # restore logger state
            try:
                gpr_logger.disabled = prev_disabled
                gpr_logger.setLevel(prev_level)
            except Exception:
                pass

    # Run getGPR in parallel to speed up many EC lookups; silence individual calls here
    if missing_ec:

        with ThreadPoolExecutor(max_workers=MAX_WORKERS_GPR) as ex:
            futures = {
                ex.submit(call_getGPR, ec, session, True): ec for ec in missing_ec
            }
            # Iterate over completed futures and print a concise per-future status
            for fut in tqdm(
                as_completed(futures), total=len(futures), desc="getGPR ECs"
            ):
                try:
                    # Add timeout to prevent hanging indefinitely
                    # If a single EC takes more than 120 seconds, skip it
                    ec_ret, res = fut.result(timeout=120)
                    if res is not None:
                        getgpr_cache[ec_ret] = res
                    # Print a short completion line (thread-safe)
                    try:
                        if res is None:
                            tqdm.write(f"getGPR result: {ec_ret} -> FAILED/None")
                        else:
                            names = None
                            try:
                                names = res[1]
                            except Exception:
                                pass
                            if isinstance(names, (list, tuple)):
                                tqdm.write(
                                    f"getGPR result: {ec_ret} -> genes={len(names)} sample={names[:5]}"
                                )
                            else:
                                tqdm.write(f"getGPR result: {ec_ret} -> {type(res)}")
                    except Exception:
                        LOGGER.debug("Completed getGPR for %s", ec_ret)
                except TimeoutError:
                    # Get the EC number from the futures dict
                    ec_num = futures.get(fut, "unknown")
                    LOGGER.error("getGPR timed out after 120s for EC: %s", ec_num)
                    try:
                        tqdm.write(f"getGPR TIMEOUT: {ec_num}")
                    except Exception:
                        pass
                except Exception:
                    LOGGER.exception("getGPR future failed")
    save_cache(GETGPR_CACHE_FILE, getgpr_cache)

    # Gather all unique gene symbols from getgpr results
    gene_symbols = set()
    for ec, res in getgpr_cache.items():
        try:
            # res[1] holds gene symbols in original build_model usage
            genes = res[1]
            for g in genes:
                if g:
                    gene_symbols.add(g)
        except Exception:
            continue

    LOGGER.info("Unique gene symbols to annotate: %d", len(gene_symbols))

    # Build a mapping gene_symbol -> list of biocyc ids from getgpr_cache to
    # speed up lookups later when creating genelist2. Normalize symbols to strings.
    gene_to_biocyc = {}
    for ec, res in getgpr_cache.items():
        try:
            names = res[1]
            biocyc_ids = res[2]
            for idx, name in enumerate(names):
                if not name:
                    continue
                try:
                    bid = biocyc_ids[idx] if idx < len(biocyc_ids) else None
                except Exception:
                    bid = None
                gene_to_biocyc.setdefault(str(name), []).append(bid)
        except Exception:
            continue

    # Batch fetch Ensembl annotations for genes not present in cache
    to_fetch = [g for g in gene_symbols if g and g not in ensembl_cache]
    LOGGER.info("Fetching %d Ensembl annotations in batches", len(to_fetch))
    if to_fetch:
        # fetch_ensembl_annotations supports batch fetching; run in smaller batches
        try:
            batches = [
                to_fetch[i : i + ENSEMBL_BATCH_SIZE]
                for i in range(0, len(to_fetch), ENSEMBL_BATCH_SIZE)
            ]
            for batch in tqdm(batches, desc="Ensembl batches"):
                try:
                    fetched = fetch_ensembl_annotations(
                        batch, batch_size=len(batch), max_workers=10
                    )
                    if fetched:
                        ensembl_cache.update(fetched)
                        save_cache(ENSEMBL_CACHE_FILE, ensembl_cache)
                except Exception:
                    LOGGER.exception("fetch_ensembl_annotations sub-batch failed")
        except Exception:
            LOGGER.exception("fetch_ensembl_annotations batch failed")

    # Build a mapping from unique GPR strings to getLocation outputs to avoid repeated calls
    gpr_to_location = {}

    # We'll iterate over getgpr_cache values and for each distinct GPR string call getLocation once
    unique_gpr_values = set()
    ec_to_gprvals = {}
    for ec, res in getgpr_cache.items():
        try:
            # res[-1] in previous script was a dict of sGPR->gpr string values
            sGPRdict = res[-1] if len(res) >= 4 else {}
            if isinstance(sGPRdict, dict):
                vals = list(sGPRdict.values())
            else:
                vals = [sGPRdict]
            ec_to_gprvals[ec] = vals
            unique_gpr_values.update([v for v in vals if v])
        except Exception:
            continue

    LOGGER.info("Unique GPR strings to run location on: %d", len(unique_gpr_values))

    # For each unique gpr string, determine the genelist1/genelist2 and call getLocation once
    for gpr in tqdm(list(unique_gpr_values), desc="GPR -> location"):
        try:
            # derive genelist1 and genelist2: find which known gene symbols appear in the gpr
            # Use whole-word matching to avoid accidental substring matches (e.g. "ACCC" in other tokens)
            # Use a proper word-boundary and case-insensitive matching so
            # gene symbols like 'CRAT' are found inside GPR strings such as
            # '[CRAT*1]'. The previous code used a double-escaped "\\b"
            # which searches for a literal backslash + 'b' and never matched.
            genelist1 = [
                g
                for g in gene_symbols
                if g
                and re.search(
                    r"\b" + re.escape(str(g)) + r"\b", gpr, flags=re.IGNORECASE
                )
            ]
            # genelist2: map gene symbols to biocyc ids using the precomputed mapping
            genelist2 = []
            for g in genelist1:
                bids = gene_to_biocyc.get(str(g)) or []
                # extend with all known biocyc ids for the gene (if any)
                for b in bids:
                    if b:
                        genelist2.append(b)

            if not genelist1:
                LOGGER.debug(
                    "No gene_symbol matches found for GPR: %r (sample length=%d)",
                    gpr[:200] if isinstance(gpr, str) else gpr,
                    len(gpr) if isinstance(gpr, str) else 0,
                )

            # Call getLocation once for this gpr
            LOGGER.debug("Calling getLocation for gpr (len genes=%d)", len(genelist1))
            loc = getLocation(
                gpr,
                genelist1,
                genelist2,
                1,
                os.path.join(project_root, "files", "dict_compartments.pkl"),
                session,
                ensembl_cache,
            )
            gpr_to_location[gpr] = loc
        except Exception:
            LOGGER.exception("getLocation failed for gpr %r", gpr)

    # Now process reactions using the cached data and expand into compartment-specific reactions
    LOGGER.info("Processing reactions using cached lookups and expanding compartments")

    # Prepare compartment mappings (reuse logic from original build_model)
    excel_file = os.path.join(project_root, "files", "ListOfCompartments_sept2024.xlsx")
    compartments_sheet_name = "Def-Compartments"
    column_index = 2
    abbreviations_sheet_name = "Endo1a_abb"
    location_pkl_file = os.path.join(project_root, "files", "dict_compartments.pkl")
    comp_abb_file = os.path.join(
        project_root, "files", "compartments_abbreviations.pkl"
    )

    # Ensure compartment pickle exists (function will create it if missing)
    try:
        create_compartments_dict_bm(
            excel_file, compartments_sheet_name, column_index, location_pkl_file
        )
    except Exception:
        LOGGER.debug(
            "create_compartments_dict_bm failed or file already present", exc_info=True
        )
    try:
        create_comp_abbreviations_dict_bm(
            excel_file, abbreviations_sheet_name, comp_abb_file
        )
    except Exception:
        LOGGER.debug(
            "create_comp_abbreviations_dict_bm failed or file already present",
            exc_info=True,
        )

    # Read abbreviations mapping (abbrev -> full name)
    try:
        comp_dict = pd.read_pickle(comp_abb_file)
    except Exception:
        comp_dict = {}

    # create mapping full_name -> abbrev and abbrev -> full name
    try:
        c, comp_dict = compartment_file_to_dict_bm(
            excel_file, compartments_sheet_name, comp_dict, comp_abb_file
        )
    except Exception:
        # Fallback: try the simpler compartment_file_to_dict if available
        try:
            from functions.function_bm_gdb import compartment_file_to_dict

            c = compartment_file_to_dict()
        except Exception:
            c = {}

    model2 = model.copy()
    variables = defaultdict(list)

    # Precompute helper structures from original model
    RxnPath = defaultdict(list)
    for p in model.groups:
        for r in p.members:
            RxnPath[r.id].append(p.id)

    # Rxn: mapping to detect existing MAR ids per reaction equation root
    Rxn = defaultdict(list)
    for x in tqdm(model.reactions, desc="processing reactions"):
        Rxn[re.sub("[a-z]+", "", x.reaction)].extend(
            [[x.id, re.findall("[a-z]+", x.reaction)[0]]]
        )

    # Starting MAR ID
    try:
        RxnID = (
            int(re.findall("[0-9]+", sorted([x.id for x in model.reactions])[-1])[0])
            + 1
        )
    except Exception:
        RxnID = 1

    ListOfMetFrom = {x.id: x.formula for x in model.metabolites}

    # counters used in original script
    n = 1
    eList = []
    RxnList = []
    listA, listB, ComptoID2 = list(), list(), {}

    # Process each reaction: build listOfgeneList5 from cached getgpr/getLocation results
    for x in model.reactions:
        x2 = model2.reactions.get_by_id(x.id)
        annotation2 = x2.annotation if hasattr(x2, "annotation") else {}
        bounds = x2.bounds

        # collect ECs
        ec_value = annotation2.get("ec-code") if annotation2 else None
        ecs = []
        if isinstance(ec_value, (list, tuple)):
            ecs = [e for e in ec_value if e]
        elif ec_value:
            ecs = [ec_value]

        listOfgeneList5 = []
        for ec in ecs:
            if ec in getgpr_cache:
                res = getgpr_cache[ec]
                # res[-1] is rxn_gpr_dict mapping reaction ids -> gpr strings
                sgpr_dict = res[-1] if len(res) >= 6 else {}
                if isinstance(sgpr_dict, dict) and sgpr_dict:
                    # iterate over unique gpr strings for that EC
                    seen_vals = set()
                    for key, gpr_expr in sgpr_dict.items():
                        if not gpr_expr or gpr_expr in seen_vals:
                            continue
                        seen_vals.add(gpr_expr)
                        # fetch precomputed location result
                        loc = gpr_to_location.get(gpr_expr)
                        if not loc:
                            # attempt to reconstruct genelist1/genelist2 heuristically
                            genelist1 = [
                                g for g in ensembl_cache.keys() if g in gpr_expr
                            ]
                            genelist2 = []
                            for e2, r2 in getgpr_cache.items():
                                try:
                                    names = r2[1]
                                    biocyc_ids = r2[2]
                                    for idx, name in enumerate(names):
                                        if name in genelist1 and idx < len(biocyc_ids):
                                            genelist2.append(biocyc_ids[idx])
                                except Exception:
                                    continue
                            try:
                                loc = getLocation(
                                    gpr_expr,
                                    genelist1,
                                    genelist2,
                                    1,
                                    location_pkl_file,
                                    session,
                                    ensembl_cache,
                                )
                            except Exception:
                                loc = ({}, {}, {}, {}, 0)
                            gpr_to_location[gpr_expr] = loc
                        listOfgeneList5.append(loc)

        # if we have any location-derived gene lists, melt them into geneList5 and apply expansion
        if listOfgeneList5:
            geneList5 = meltGeneList(listOfgeneList5)
            variables[x.id] = geneList5

            # species and species3 are as in original script
            species = [s.id for s in x2.reactants] + [s.id for s in x2.products]
            species3 = [
                (
                    {s.split(" ")[0]: -1}
                    if len(s.split(" ")) == 1
                    else {s.split(" ")[1]: -abs(float(s.split(" ")[0]))}
                )
                for s in re.split(" --> | <=> ", re.sub("[a-z]+", "", x2.reaction))[
                    0
                ].split(" + ")
            ] + [
                (
                    {s.split(" ")[0]: 1}
                    if len(s.split(" ")) == 1
                    else {s.split(" ")[1]: float(s.split(" ")[0])}
                )
                for s in re.split(" --> | <=> ", re.sub("[a-z]+", "", x2.reaction))[
                    1
                ].split(" + ")
            ]
            species3 = dict(ChainMap(*species3))

            for CSL in [*geneList5[0]]:
                CSL2 = CSL
                if CSL not in c:
                    CSL = "cytosol"

                # compute ID abbreviation
                if CSL in c:
                    ID = c[CSL]
                elif len(re.sub(r" $", "", re.sub(r"^ ", "", CSL)).split(" ")) > 1:
                    ID = (
                        (CSL.strip().split(" ")[0][0] + CSL.strip().split(" ")[1][0])
                        .lower()
                        .replace(" ", "")
                    )
                else:
                    if len(CSL.split(" ")) > 1:
                        ID = CSL[0:3].lower().replace(" ", "")
                    else:
                        ID = CSL[0:2].lower().replace(" ", "")
                if not c.get(CSL) and ID in listA + listB:
                    r = re.compile(ID)
                    ID = ID + str(len(list(filter(r.match, listA + listB))) + 1)
                ComptoID2[CSL] = ID
                listA.append(ID)
                listB.append(ID)

                # check duplicate with existing reaction variants
                if ID not in [
                    rxn[1] for rxn in Rxn[re.sub("[a-z]+[0-9]*", "", x2.reaction)]
                ]:
                    reaction2 = Reaction("MAR" + str(RxnID + n))
                    reaction2.lower_bound = bounds[0]
                    reaction2.upper_bound = bounds[1]

                    for species2 in [
                        [re.sub(r"[a-z][0-9]*", ComptoID2[CSL], s), s] for s in species
                    ]:
                        try:
                            m = model2.metabolites.get_by_id(species2[0])
                        except KeyError:
                            src = model2.metabolites.get_by_id(species2[1])
                            m = Metabolite(
                                species2[0],
                                charge=src.charge,
                                formula=src.formula,
                                name=src.name,
                                compartment=ID,
                            )
                            m.annotation = src.annotation

                        if m.id not in model2.metabolites:
                            model2.metabolites.add(m)

                        if not model2.compartments.get(ID):
                            model2.compartments[ID] = CSL2

                        # add stoichiometry
                        key = re.sub(r"[a-z]+[0-9]*", "", species2[0])
                        stoich = species3.get(key, 1)
                        reaction2.add_metabolites({m: stoich})

                    # set gene_reaction_rule
                    try:
                        gr_rule = geneList5[2].get(CSL2, "")
                        if "or" in gr_rule and "and" in gr_rule:
                            reaction2.gene_reaction_rule = " and ".join(
                                ["(" + v + ")" for v in gr_rule.split(" and ")]
                            )
                        else:
                            reaction2.gene_reaction_rule = gr_rule
                    except Exception:
                        reaction2.gene_reaction_rule = ""

                    reaction2.annotation = annotation2.copy() if annotation2 else {}

                    # set sGPR annotation if available
                    try:
                        comp = list(reaction2.compartments)
                        compartment_name = (
                            comp_dict.get(comp[0], comp[0]) if comp else ""
                        )
                        reaction2.annotation["sGPR"] = geneList5[0].get(
                            compartment_name, ""
                        )
                    except Exception:
                        pass

                    model2.add_reactions([reaction2])

                    # annotate genes with ensembl and hgnc.symbol when available
                    try:
                        for gene in re.split(
                            r" or | and ", reaction2.gene_reaction_rule
                        ):
                            gene = gene.replace(")", "").replace("(", "").strip()
                            if gene:
                                ensemble = ensembl_cache.get(gene)
                                if ensemble and isinstance(ensemble, dict):
                                    model2.genes.get_by_id(gene).annotation[
                                        "ensembl"
                                    ] = ensemble.get("ensembl")
                                else:
                                    model2.genes.get_by_id(gene).annotation[
                                        "ensembl"
                                    ] = ensemble
                                try:
                                    model2.genes.get_by_id(gene).annotation[
                                        "hgnc.symbol"
                                    ] = geneList5[3].get(gene, gene)
                                except Exception:
                                    pass
                    except Exception:
                        LOGGER.debug(
                            "Failed adding gene annotations for reaction %s",
                            reaction2.id,
                            exc_info=True,
                        )

                    # add to original pathway group if possible
                    try:
                        if RxnPath.get(x.id):
                            model2.groups.get_by_id(RxnPath[x.id][0]).members.add(
                                reaction2
                            )
                    except Exception:
                        pass

                    RxnList.append(reaction2.reaction)
                    n += 1
                else:
                    # existing reaction: try to set sGPR on existing reaction
                    try:
                        model2.reactions.get_by_id(
                            x2.id
                        ).gene_reaction_rule = geneList5[2].get(
                            CSL2, model2.reactions.get_by_id(x2.id).gene_reaction_rule
                        )
                        comp = list(model2.reactions.get_by_id(x2.id).compartments)
                        compartment_name = (
                            comp_dict.get(comp[0], comp[0]) if comp else ""
                        )
                        model2.reactions.get_by_id(x2.id).annotation["sGPR"] = (
                            geneList5[0].get(compartment_name, "")
                        )
                    except Exception:
                        pass

        else:
            # no variables for this reaction
            continue

    # Persist caches and variables
    save_cache(ENSEMBL_CACHE_FILE, ensembl_cache)
    save_cache(KEGG_CACHE_FILE, kegg_cache)
    save_cache(GETGPR_CACHE_FILE, getgpr_cache)
    save_cache(os.path.join(project_root, "files", "variables_batch.pkl"), variables)

    # Write modified model
    LOGGER.info("Writing output model to %s", OUTPUT_MODEL)
    cobra.io.write_sbml_model(model2, OUTPUT_MODEL)
    LOGGER.info("Done. Wrote %s", OUTPUT_MODEL)


if __name__ == "__main__":
    main()
