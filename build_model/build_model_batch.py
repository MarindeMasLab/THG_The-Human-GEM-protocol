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
from datetime import datetime
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
import pickle

import requests
import cobra
import pandas as pd
import copy
import xml.etree.ElementTree as ET
from collections import ChainMap
from cobra import Reaction, Metabolite
from tqdm import tqdm
import urllib.parse
import sympy

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
from functions.equations_bm_gdb import *

LOGGER = logging.getLogger("build_model_batch")
logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s %(levelname)s %(name)s %(message)s"
)

# Config (tweak as needed)
INPUT_MODEL = os.path.join(project_root, "models", "THG-beta1.1.1_251031.xml")
OUTPUT_MODEL = os.path.join(project_root, "models", "THG-beta-batch.xml")

# create a date tag like "_yymmdd"
_date_tag = "_" + datetime.now().strftime("%y%m%d")

OUTPUT_MODEL_TMP = os.path.join(
    project_root, "models", f"THG-beta-batch-tmp{_date_tag}.xml"
)
OUTPUT_MODEL_FINAL = os.path.join(
    project_root, "models", f"THG-beta-batch{_date_tag}.xml"
)
OUTPUT_ERRORS = os.path.join(project_root, "files", f"Rxn2Fix{_date_tag}.txt")
ENSEMBL_CACHE_FILE = os.path.join(
    project_root, "files", "caches", "ensembl_cache_batch.pkl"
)
KEGG_CACHE_FILE = os.path.join(
    project_root, "files", "caches", "kegg_reaction_cache_batch.pkl"
)
GETGPR_CACHE_FILE = os.path.join(
    project_root, "files", "caches", "getgpr_cache_batch.pkl"
)

KEGG_BATCH_SIZE = 10
ENSEMBL_BATCH_SIZE = 10
MAX_WORKERS_GPR = 10


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
            # LOGGER.debug("Fetching KEGG batch: %s", batch[:5])  # Show first 5 IDs
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
            # LOGGER.debug("Fetched %d/%d entries in batch", fetched_count, len(batch))

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
    # # fallback: try to find Rxxxx in reaction id or name
    # m = re.search(r"R[0-9]{5}", reac.id)
    # if m:
    #     return m.group(0)
    # m = re.search(r"R[0-9]{5}", getattr(reac, "name", "") or "")
    # if m:
    #     return m.group(0)
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
    # Include annotation fields from the model
    not_found_kegg = [kid for kid in kegg_ids if kegg_cache.get(kid) is None]
    if not_found_kegg:
        csv_path = os.path.join(project_root, "files", "kegg_reactions_not_found.csv")
        try:
            # Build a mapping from KEGG ID to reaction objects
            kegg_to_reactions = defaultdict(list)
            for r in model.reactions:
                kid = extract_kegg_reaction_id(r)
                if kid and kid in not_found_kegg:
                    kegg_to_reactions[kid].append(r)

            # Collect data for CSV
            csv_data = []
            for kid in sorted(not_found_kegg):
                reactions = kegg_to_reactions.get(kid, [])
                if reactions:
                    # Use the first reaction if multiple have the same KEGG ID
                    r = reactions[0]
                    row = {
                        "KEGG_ID": kid,
                        "Reaction_ID": r.id,
                        "Reaction_Name": getattr(r, "name", ""),
                        "Reaction_Equation": r.reaction,
                    }
                    # Add all annotation fields
                    if hasattr(r, "annotation"):
                        for key, value in r.annotation.items():
                            # Convert lists/tuples to comma-separated strings
                            if isinstance(value, (list, tuple)):
                                row[f"annotation_{key}"] = ", ".join(
                                    str(v) for v in value
                                )
                            else:
                                row[f"annotation_{key}"] = str(value)
                    csv_data.append(row)
                else:
                    # KEGG ID without matching reaction (shouldn't happen, but handle it)
                    csv_data.append({"KEGG_ID": kid})

            df = pd.DataFrame(csv_data)
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

    def call_getGPR(ec_value, session_obj=None, verbose=False):
        """Wrapper around getGPR that supports a verbose flag.

        If verbose is False, temporarily silence the gpr_def logger for quieter output.
        Returns (ec, result) where result is the value returned by getGPR or None on failure.
        """
        logger_name = "functions.gpr.gpr_def"
        gpr_logger = logging.getLogger(logger_name)
        prev_disabled = getattr(gpr_logger, "disabled", False)
        prev_level = getattr(gpr_logger, "level", logging.NOTSET)
        try:
            if not verbose:
                gpr_logger.disabled = True
            # ensure we have a session per call if none provided
            if session_obj is None:
                session_local = setup_biocyc_session()
            else:
                session_local = session_obj
            res = getGPR(ec_value, session_local)
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
                ex.submit(call_getGPR, ec, session, verbose=False): ec
                for ec in missing_ec
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
                except TimeoutError:
                    # Get the EC number from the futures dict
                    ec_num = futures.get(fut, "unknown")
                    LOGGER.error("getGPR timed out after 120s for EC: %s", ec_num)
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
    # Precompute genelist1/genelist2 for each GPR so we can prefetch the
    # external pages (BioCyc, UniProt, KEGG) in batches and cache them.
    gpr_to_genelists = {}
    all_biocyc_ids = set()
    all_gene_symbols = set()

    for gpr in unique_gpr_values:
        try:
            genelist1 = [
                g
                for g in gene_symbols
                if g
                and re.search(
                    r"\b" + re.escape(str(g)) + r"\b", gpr, flags=re.IGNORECASE
                )
            ]
            genelist2 = []
            for g in genelist1:
                for b in gene_to_biocyc.get(str(g)) or []:
                    if b:
                        genelist2.append(b)

            gpr_to_genelists[gpr] = (genelist1, genelist2)
            all_biocyc_ids.update([b.upper() for b in genelist2 if b])
            all_gene_symbols.update([g.upper() for g in genelist1 if g])
        except Exception:
            LOGGER.exception("Failed preparing genelist for gpr %r", gpr)

    LOGGER.info(
        "Prefetching pages for %d Biocyc IDs and %d gene symbols",
        len(all_biocyc_ids),
        len(all_gene_symbols),
    )

    # Use BioVelo webservice to resolve gene symbols to BioCyc IDs in batches
    BIOVELO_BATCH_SIZE = 50
    BIOVELO_CACHE_FILE = os.path.join(
        project_root, "files", "caches", "biovelo_cache.pkl"
    )
    BIOVELO_LOC_CACHE_FILE = os.path.join(
        project_root, "files", "caches", "biovelo_location_cache.pkl"
    )

    try:
        biovelo_cache = load_cache(BIOVELO_CACHE_FILE)
    except Exception:
        biovelo_cache = {}

    # Determine which symbols still need resolution (uppercase normalization)
    symbols_to_resolve = [s for s in all_gene_symbols if s not in biovelo_cache]
    if symbols_to_resolve:
        LOGGER.info(
            "Resolving %d gene symbols via BioVelo in batches", len(symbols_to_resolve)
        )
        batches = [
            symbols_to_resolve[i : i + BIOVELO_BATCH_SIZE]
            for i in range(0, len(symbols_to_resolve), BIOVELO_BATCH_SIZE)
        ]

        def _fetch_biovelo(batch):
            # Build BioVelo query: find HUMAN genes whose name matches any of the symbols in the batch
            # Query: [x:x<-HUMAN^^genes, x^name in ("G1","G2",...)]
            names = ",".join(f'"{n}"' for n in batch)
            query = f"[x:x<-HUMAN^^genes, x^name in ({names})]"
            q = urllib.parse.quote(query, safe="")
            url = f"https://websvc.biocyc.org/xmlquery?query={q}&detail=low"
            try:
                r = requests.get(url, timeout=30)
                r.raise_for_status()
                text = r.text
                return url, batch, text
            except Exception:
                LOGGER.debug("BioVelo fetch failed for batch %s", batch, exc_info=True)
                return url, batch, ""

        # Run batch queries in parallel
        with ThreadPoolExecutor(max_workers=5) as ex:
            futures = {ex.submit(_fetch_biovelo, b): b for b in batches}
            for fut in tqdm(
                as_completed(futures), total=len(futures), desc="BioVelo batches"
            ):
                url, batch, text = fut.result()
                # store response in cache keyed by URL and by individual symbols
                biovelo_cache[url] = text
                # Parse XML to extract product/protein frameids robustly (NP_/XP_/ENSG_/HS_/AT_/etc.)
                ids = set()
                if text:
                    try:
                        root = ET.fromstring(text)
                        # Look for Gene elements and their product -> Protein frameids
                        for g in root.findall(".//Gene"):
                            # look for product elements containing Protein children
                            for prod in g.findall(".//product"):
                                p = prod.find("Protein")
                                if p is not None and p.get("frameid"):
                                    ids.add(p.get("frameid"))
                                else:
                                    # fallback: any child element with frameid
                                    for child in prod:
                                        if child is not None and child.get("frameid"):
                                            ids.add(child.get("frameid"))

                        # As additional fallback, search for common frameid patterns in the raw text
                        if not ids:
                            # capture ENSG, NP_, XP_, HS, AT, G66, and similar frameids
                            ids.update(
                                re.findall(
                                    r"(ENSG[0-9]+|NP_[0-9]+\.?[0-9]*|XP_[0-9]+\.?[0-9]*|HS[0-9]+|AT[0-9Gg\-]+|G66-[0-9]+|RNA[0-9]+)",
                                    text,
                                )
                            )
                    except Exception:
                        LOGGER.debug(
                            "Failed to parse BioVelo XML for batch %s",
                            batch,
                            exc_info=True,
                        )

                # Ensure we store a list (possibly empty) for each requested symbol
                for s in batch:
                    biovelo_cache[s] = list(ids)

        # save biovelo cache
        try:
            save_cache(BIOVELO_CACHE_FILE, biovelo_cache)
        except Exception:
            LOGGER.exception("Failed to save BioVelo cache")

    # After resolving gene symbols to Biocyc IDs, fetch subcellular location
    # information using the BioCyc getxml endpoint. This retrieves protein objects
    # with their location data (CCO/CCI frameids) which we'll resolve to names later.
    try:
        biovelo_loc_cache = load_cache(BIOVELO_LOC_CACHE_FILE)
    except Exception:
        biovelo_loc_cache = {}

    # Build list of Biocyc IDs we may want to query (already normalized to UPPER)
    ids_to_query = [
        b for b in sorted(all_biocyc_ids) if b and b not in biovelo_loc_cache
    ]
    if ids_to_query:
        LOGGER.info(
            "Fetching protein locations via getxml for %d Biocyc IDs", len(ids_to_query)
        )

        def _fetch_protein_location(bid):
            """Fetch protein object using getxml and extract location frameids.

            Returns (bid, list_of_location_dicts) where each dict has:
                - frameid: the CCO/CCI compartment frameid
                - orgid: the organism ID
            """
            try:
                # Try several candidate getxml URL forms. Some Biocyc IDs in the
                # caches are in the older 'G-12345' format while the webservice
                # uses 'HSxxxxx' for human proteins. Build a list of candidates
                # that includes HS variants when appropriate.
                candidates = []
                bstr = str(bid)
                # If bid already contains an org prefix, try it directly first
                if ":" in bstr:
                    candidates.append(
                        f"https://websvc.biocyc.org/getxml?{bstr}&detail=full"
                    )

                # If it's a legacy G-12345 or G12345 id, produce HS variants
                m = re.match(r"^G-?(\d+)$", bstr, flags=re.IGNORECASE)
                if m:
                    digits = m.group(1)
                    hs_plain = f"HS{digits}"
                    hs_padded = f"HS{digits.zfill(5)}"
                    candidates.append(
                        f"https://websvc.biocyc.org/getxml?HUMAN:{hs_plain}&detail=full"
                    )
                    if hs_padded != hs_plain:
                        candidates.append(
                            f"https://websvc.biocyc.org/getxml?HUMAN:{hs_padded}&detail=full"
                        )

                # Default attempts: HUMAN:bid then plain bid
                candidates.append(
                    f"https://websvc.biocyc.org/getxml?HUMAN:{bstr}&detail=full"
                )
                candidates.append(
                    f"https://websvc.biocyc.org/getxml?{bstr}&detail=full"
                )

                locations = []
                for url in candidates:
                    try:
                        r = session.get(url, timeout=30)
                        if r.status_code == 404:
                            LOGGER.debug("getxml 404 for %s", url)
                            continue
                        r.raise_for_status()
                        xml_text = r.text
                        # Parse XML and extract <location>/<cco> entries
                        try:
                            root = ET.fromstring(xml_text)
                        except Exception:
                            LOGGER.debug(
                                "Failed to parse XML from %s", url, exc_info=True
                            )
                            continue

                        for loc_elem in root.findall(".//location"):
                            for cco in loc_elem.findall(".//cco"):
                                frameid = cco.get("frameid")
                                orgid = cco.get("orgid")
                                if frameid:
                                    locations.append(
                                        {"frameid": frameid, "orgid": orgid}
                                    )

                        if locations:
                            return bid, locations

                        # If no locations, try to follow product->Protein frameids inside this object
                        prod_frameids = []
                        for prod in root.findall(".//product"):
                            p = prod.find("Protein")
                            if p is not None and p.get("frameid"):
                                prod_frameids.append(p.get("frameid"))
                            else:
                                for child in prod:
                                    if child is not None and child.get("frameid"):
                                        prod_frameids.append(child.get("frameid"))

                        for pf in prod_frameids:
                            pf = str(pf)
                            # try common forms for product frameid as well
                            pf_candidates = [
                                f"https://websvc.biocyc.org/getxml?HUMAN:{pf}&detail=full",
                                f"https://websvc.biocyc.org/getxml?{pf}&detail=full",
                            ]
                            for purl in pf_candidates:
                                try:
                                    r2 = session.get(purl, timeout=30)
                                    if r2.status_code == 404:
                                        continue
                                    r2.raise_for_status()
                                    try:
                                        root2 = ET.fromstring(r2.text)
                                    except Exception:
                                        continue
                                    for loc_elem in root2.findall(".//location"):
                                        for cco in loc_elem.findall(".//cco"):
                                            frameid = cco.get("frameid")
                                            orgid = cco.get("orgid")
                                            if frameid:
                                                locations.append(
                                                    {"frameid": frameid, "orgid": orgid}
                                                )
                                    if locations:
                                        return bid, locations
                                except Exception:
                                    LOGGER.debug(
                                        "getxml product fetch failed for %s via %s",
                                        (pf, purl),
                                        exc_info=True,
                                    )

                        # nothing found for this URL; try next candidate
                    except Exception:
                        LOGGER.debug(
                            "getxml location fetch failed for %s", (url,), exc_info=True
                        )

                return bid, []
            except Exception:
                LOGGER.debug("getxml location fetch failed for %s", bid, exc_info=True)
                return bid, []

        # Run in parallel with modest concurrency
        with ThreadPoolExecutor(max_workers=6) as ex:
            futures = {ex.submit(_fetch_protein_location, b): b for b in ids_to_query}
            for fut in tqdm(
                as_completed(futures), total=len(futures), desc="Protein locations"
            ):
                bid, locations = fut.result()
                biovelo_loc_cache[bid] = locations

        # persist location cache
        try:
            save_cache(BIOVELO_LOC_CACHE_FILE, biovelo_loc_cache)
        except Exception:
            LOGGER.exception("Failed to save BioVelo location cache")

    # Collect all unique compartment frameids (CCO/CCI) from the location data
    # so we can batch-resolve them to human-readable names
    compartment_frameids = set()
    for bid, locations in biovelo_loc_cache.items():
        if isinstance(locations, list):
            for loc in locations:
                if isinstance(loc, dict) and loc.get("frameid"):
                    compartment_frameids.add(
                        (loc.get("orgid", "HUMAN"), loc["frameid"])
                    )

    LOGGER.info(
        "Found %d unique compartment frameids to resolve", len(compartment_frameids)
    )

    # Cache for compartment frameid -> common name mapping
    COMPARTMENT_NAME_CACHE_FILE = os.path.join(
        project_root, "files", "compartment_name_cache.pkl"
    )
    try:
        compartment_name_cache = load_cache(COMPARTMENT_NAME_CACHE_FILE)
    except Exception:
        compartment_name_cache = {}

    # Resolve compartment frameids to common names in batch
    frameids_to_resolve = [
        (org, fid)
        for org, fid in compartment_frameids
        if f"{org}:{fid}" not in compartment_name_cache
    ]

    if frameids_to_resolve:
        LOGGER.info(
            "Resolving %d compartment names via getxml", len(frameids_to_resolve)
        )

        def _resolve_compartment_name(org_fid_tuple):
            """Resolve a compartment frameid to its common name using getxml."""
            org, frameid = org_fid_tuple
            key = f"{org}:{frameid}"
            try:
                url = f"https://websvc.biocyc.org/getxml?{org}:{frameid}&detail=low"
                response = session.get(url, timeout=10)
                response.raise_for_status()
                xml_text = response.text

                # Parse and extract common-name
                root = ET.fromstring(xml_text)
                common_name_elem = root.find(".//common-name")
                if common_name_elem is not None and common_name_elem.text:
                    return key, common_name_elem.text

                # Fallback: try to find it in string elements
                strings = [el.text for el in root.findall(".//string") if el.text]
                if strings:
                    return key, strings[0]

                return key, frameid  # Return frameid if we can't resolve it
            except Exception:
                LOGGER.debug(
                    "Failed to resolve compartment name for %s", key, exc_info=True
                )
                return key, frameid

        # Resolve in parallel
        with ThreadPoolExecutor(max_workers=10) as ex:
            futures = {
                ex.submit(_resolve_compartment_name, ofid): ofid
                for ofid in frameids_to_resolve
            }
            for fut in tqdm(
                as_completed(futures), total=len(futures), desc="Compartment names"
            ):
                key, common_name = fut.result()
                compartment_name_cache[key] = common_name

        # Save compartment name cache
        try:
            save_cache(COMPARTMENT_NAME_CACHE_FILE, compartment_name_cache)
        except Exception:
            LOGGER.exception("Failed to save compartment name cache")

    # Build location mappings from cached BioCyc data instead of calling getLocation
    # This constructs the same data structure that getLocation returns:
    # (RuleLoc, RuleLoc2, RuleLoc3, RuleLoc4) where:
    # - RuleLoc: {compartment: SGPR with *N notation}
    # - RuleLoc2: {compartment: GPR without *N notation}
    # - RuleLoc3: {compartment: GPR with Ensembl IDs}
    # - RuleLoc4: {Ensembl ID: [gene names]}

    LOGGER.info(
        "Building location mappings from cached BioCyc data for %d GPRs",
        len(gpr_to_genelists),
    )

    for gpr, (genelist1, genelist2) in tqdm(
        gpr_to_genelists.items(), desc="GPR -> location"
    ):
        try:
            if not genelist1:
                LOGGER.debug(
                    "No gene_symbol matches found for GPR: %r (sample length=%d)",
                    gpr[:200] if isinstance(gpr, str) else gpr,
                    len(gpr) if isinstance(gpr, str) else 0,
                )
                # Return empty structures like getLocation does on error
                gpr_to_location[gpr] = ({}, {}, {}, {}, 0)
                continue

            # Build gene -> compartments mapping from cached BioCyc location data
            gene_compartments = defaultdict(set)
            for gene_name in genelist1:
                # Find corresponding BioCyc IDs for this gene
                biocyc_ids = gene_to_biocyc.get(str(gene_name), [])
                for bid in biocyc_ids:
                    if bid and bid.upper() in biovelo_loc_cache:
                        locations = biovelo_loc_cache[bid.upper()]
                        if isinstance(locations, list):
                            for loc_info in locations:
                                if isinstance(loc_info, dict) and loc_info.get(
                                    "frameid"
                                ):
                                    org = loc_info.get("orgid", "HUMAN")
                                    frameid = loc_info["frameid"]
                                    key = f"{org}:{frameid}"
                                    # Resolve frameid to common name
                                    compartment_name = compartment_name_cache.get(
                                        key, "cytosol"
                                    )
                                    gene_compartments[gene_name].add(compartment_name)

            # If no locations found for any gene, default to cytosol
            if not gene_compartments:
                logging.warning(
                    "No compartment data found for GPR %r; defaulting to cytosol", gpr
                )
                for gene_name in genelist1:
                    gene_compartments[gene_name].add("cytosol")

            # Build the 4-tuple return structure similar to getLocation
            RuleLoc = {}  # {compartment: SGPR with *N}
            RuleLoc2 = {}  # {compartment: GPR without *N}
            RuleLoc3 = {}  # {compartment: GPR with Ensembl IDs}
            RuleLoc4 = defaultdict(list)  # {Ensembl ID: [gene names]}

            # Group genes by compartment
            compartment_genes = defaultdict(list)
            for gene_name in genelist1:
                compartments = gene_compartments.get(gene_name, {"cytosol"})
                for comp in compartments:
                    compartment_genes[comp].append(gene_name)

            # Build the location-specific GPRs
            for compartment, genes in compartment_genes.items():
                if not genes:
                    continue

                # Build GPR string for this compartment
                # Extract the GPR pattern from original gpr for these genes
                gpr_clean = re.sub(r"\*[0-9]+", "", gpr)

                # Simple approach: if all genes in genelist1 are in this compartment,
                # use the full GPR; otherwise build a partial GPR
                if set(genes) == set(genelist1):
                    # All genes are in this compartment
                    RuleLoc[compartment] = gpr  # with *N notation
                    RuleLoc2[compartment] = gpr_clean  # without *N notation
                else:
                    # Partial GPR - just the genes in this compartment
                    partial_gpr = " or ".join(f"[{g}]" for g in genes)
                    RuleLoc[compartment] = partial_gpr
                    RuleLoc2[compartment] = partial_gpr

                # Build Ensembl ID version
                ensembl_parts = []
                for gene_name in genes:
                    ens_data = ensembl_cache.get(gene_name)
                    if ens_data:
                        if isinstance(ens_data, dict):
                            ens_id = ens_data.get("ensembl", gene_name)
                        else:
                            ens_id = ens_data
                        ensembl_parts.append(str(ens_id))
                        RuleLoc4[str(ens_id)].append(gene_name)
                    else:
                        ensembl_parts.append(gene_name)
                        RuleLoc4[gene_name].append(gene_name)

                RuleLoc3[compartment] = " or ".join(ensembl_parts)

            # Return 5-tuple like getLocation: (RuleLoc, RuleLoc2, RuleLoc3, RuleLoc4, status)
            # status=0 means successful location determination
            gpr_to_location[gpr] = (RuleLoc, RuleLoc2, RuleLoc3, dict(RuleLoc4), 0)

        except Exception:
            LOGGER.exception("Failed building location for gpr %r", gpr)
            gpr_to_location[gpr] = ({}, {}, {}, {}, 1)  # status=1 indicates failure

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
    # Also get the BioCyc-to-Endo1a mapping (col0_to_col2)
    try:
        c, comp_dict, biocyc_to_endo1a = compartment_file_to_dict_bm(
            excel_file, compartments_sheet_name, comp_dict, comp_abb_file
        )
    except Exception:
        # Fallback: try the simpler compartment_file_to_dict if available
        try:
            from functions.function_bm_gdb import compartment_file_to_dict

            c = compartment_file_to_dict()
            biocyc_to_endo1a = {}
        except Exception:
            c = {}
            biocyc_to_endo1a = {}

    model2 = model.copy()
    # Create a case-insensitive view of the compartments mapping so lookups
    # succeed whether the compartment name from BioCyc is capitalized or not.
    # c maps full compartment name -> abbreviation. Build c_lower for lookups.
    try:
        c_lower = {k.lower(): v for k, v in c.items()}
    except Exception:
        c_lower = {}

    # Create a lowercase version of the BioCyc-to-Endo1a mapping for lookups
    try:
        biocyc_lower = {k.lower(): v.lower() for k, v in biocyc_to_endo1a.items()}
    except Exception:
        biocyc_lower = {}

    # Pre-compute compartment mappings once for all unique compartments we might encounter
    # This avoids repeated mapping logic for every reaction
    LOGGER.info("Pre-computing compartment mappings...")
    compartment_mapping_cache = {}

    def map_compartment_once(CSL_original):
        """Map a BioCyc compartment name to a standardized Endo1-a compartment.

        This function implements the 3-tier mapping strategy:
        1. Direct BioCyc-to-Endo1a mapping from Excel
        2. Direct match in Endo1-a compartments
        3. Fuzzy matching for common patterns
        4. Default to cytosol
        """
        if CSL_original in compartment_mapping_cache:
            return compartment_mapping_cache[CSL_original]

        CSL = CSL_original
        lookup = CSL.lower() if isinstance(CSL, str) else CSL
        mapping_method = "unknown"

        # Check if this is a BioCyc compartment name that needs mapping
        if lookup in biocyc_lower:
            # Map BioCyc name to Endo1-a name
            mapped_endo1a = biocyc_lower[lookup]
            # Now find the canonical form of this Endo1-a compartment
            if mapped_endo1a in c_lower:
                orig = next((k for k in c.keys() if k.lower() == mapped_endo1a), None)
                if orig:
                    CSL = orig
                    mapping_method = f"BioCyc->Endo1a ({CSL_original} -> {CSL})"
            else:
                # Endo1-a compartment not in c_lower, use it anyway
                CSL = mapped_endo1a
                mapping_method = f"BioCyc->Endo1a-unchecked ({CSL_original} -> {CSL})"
        elif lookup in c_lower:
            # Direct match in Excel Endo1-a compartments
            orig = next((k for k in c.keys() if k.lower() == lookup), None)
            if orig:
                CSL = orig
                mapping_method = "direct"
        else:
            # No direct mapping found - try fuzzy matching for common BioCyc compartment patterns
            matched = False
            if isinstance(CSL, str):
                csl_lower = CSL.lower()
                # Check for vesicle-related compartments
                if any(
                    keyword in csl_lower
                    for keyword in ["vesicle", "endosome", "granule", "vacuole"]
                ):
                    if "vesicle" in c_lower:
                        CSL = next(k for k in c.keys() if k.lower() == "vesicle")
                        matched = True
                        mapping_method = f"fuzzy-vesicle ({CSL_original} -> {CSL})"
                # Check for Golgi-related compartments
                elif "golgi" in csl_lower:
                    if "golgi apparatus" in c_lower:
                        CSL = next(
                            k for k in c.keys() if k.lower() == "golgi apparatus"
                        )
                        matched = True
                        mapping_method = f"fuzzy-golgi ({CSL_original} -> {CSL})"
                # Check for ER-related compartments
                elif (
                    "endoplasmic reticulum" in csl_lower
                    or "endoplasmic-reticulum" in csl_lower
                ):
                    if "endoplasmic reticulum" in c_lower:
                        CSL = next(
                            k for k in c.keys() if k.lower() == "endoplasmic reticulum"
                        )
                        matched = True
                        mapping_method = f"fuzzy-ER ({CSL_original} -> {CSL})"
                # Check for mitochondria-related compartments
                elif "mitochondri" in csl_lower:
                    # Try to match to inner mitochondria for membrane/matrix compartments
                    if any(
                        keyword in csl_lower
                        for keyword in ["inner", "matrix", "lumen", "cristae"]
                    ):
                        if "inner mitochondria" in c_lower:
                            CSL = next(
                                k for k in c.keys() if k.lower() == "inner mitochondria"
                            )
                            matched = True
                            mapping_method = (
                                f"fuzzy-mito-inner ({CSL_original} -> {CSL})"
                            )
                    elif "mitochondria" in c_lower:
                        CSL = next(k for k in c.keys() if k.lower() == "mitochondria")
                        matched = True
                        mapping_method = f"fuzzy-mito ({CSL_original} -> {CSL})"
                # Check for lysosome-related compartments
                elif "lyso" in csl_lower:
                    if "lysosome" in c_lower:
                        CSL = next(k for k in c.keys() if k.lower() == "lysosome")
                        matched = True
                        mapping_method = f"fuzzy-lyso ({CSL_original} -> {CSL})"
                # Check for peroxisome-related compartments
                elif "peroxisom" in csl_lower:
                    if "peroxisome" in c_lower:
                        CSL = next(k for k in c.keys() if k.lower() == "peroxisome")
                        matched = True
                        mapping_method = f"fuzzy-perox ({CSL_original} -> {CSL})"
                # Check for nucleus-related compartments
                elif "nucle" in csl_lower and "nucleolus" not in csl_lower:
                    if "nucleus" in c_lower:
                        CSL = next(k for k in c.keys() if k.lower() == "nucleus")
                        matched = True
                        mapping_method = f"fuzzy-nucleus ({CSL_original} -> {CSL})"
                # Check for plasma/cell membrane
                elif "plasma membrane" in csl_lower or "cell membrane" in csl_lower:
                    if "cell membrane" in c_lower:
                        CSL = next(k for k in c.keys() if k.lower() == "cell membrane")
                        matched = True
                        mapping_method = f"fuzzy-membrane ({CSL_original} -> {CSL})"
                # Check for extracellular
                elif "extracellular" in csl_lower:
                    if "extracellular" in c_lower:
                        CSL = next(k for k in c.keys() if k.lower() == "extracellular")
                        matched = True
                        mapping_method = f"fuzzy-extra ({CSL_original} -> {CSL})"
                # Check for cytoskeleton
                elif "cytoskeleton" in csl_lower or "cytoskeletal" in csl_lower:
                    if "cytoskeleton" in c_lower:
                        CSL = next(k for k in c.keys() if k.lower() == "cytoskeleton")
                        matched = True
                        mapping_method = f"fuzzy-cytosk ({CSL_original} -> {CSL})"

            if not matched:
                LOGGER.warning(
                    "Compartment '%s' not found in mapping, defaulting to cytosol", CSL
                )
                CSL = "cytosol"
                mapping_method = f"unmapped->cytosol ({CSL_original})"

        compartment_mapping_cache[CSL_original] = (CSL, mapping_method)
        return CSL, mapping_method

    # Collect all unique compartments from gpr_to_location to pre-map them
    all_compartments_to_map = set()
    for loc_tuple in gpr_to_location.values():
        if loc_tuple and len(loc_tuple) >= 1:
            RuleLoc = loc_tuple[0]
            if isinstance(RuleLoc, dict):
                all_compartments_to_map.update(RuleLoc.keys())

    LOGGER.info("Found %d unique compartments to map", len(all_compartments_to_map))

    # Pre-map all compartments
    mapped_summary = defaultdict(int)
    for comp in all_compartments_to_map:
        mapped_comp, method = map_compartment_once(comp)
        mapped_summary[method] += 1

    LOGGER.info("Compartment mapping summary:")
    for method, count in sorted(mapped_summary.items(), key=lambda x: -x[1]):
        LOGGER.info("  %s: %d compartments", method, count)

    variables = defaultdict(list)

    # Track any exchange reactions that have EC annotations so we can warn and inspect
    exchange_ec_warnings = []

    def warn_if_ec_on_exchange(rxn_obj, annotation_obj, context=""):
        """Log a warning if an exchange reaction has an EC annotation.

        Returns True if a warning was logged.
        """
        try:
            if not annotation_obj:
                return False
            ec = annotation_obj.get("ec-code")
            if not ec:
                return False
            # Heuristic: exchange if no reactants or no products, or typical EX id
            is_exchange_local = (
                len(rxn_obj.reactants) == 0
                or len(rxn_obj.products) == 0
                or str(rxn_obj.id).upper().startswith("EX_")
                or str(rxn_obj.id).upper().startswith("EX")
            )
            if is_exchange_local:
                LOGGER.warning(
                    "EC annotation on exchange reaction %s (%s): %r",
                    rxn_obj.id,
                    context,
                    ec,
                )
                exchange_ec_warnings.append(
                    {"reaction": rxn_obj.id, "context": context, "ec": ec}
                )
                return True
        except Exception:
            LOGGER.exception("Failed checking EC on exchange for %s", rxn_obj.id)
        return False

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

    # Track compartment usage for diagnostic purposes
    biocyc_compartments_seen = defaultdict(int)
    mapped_compartments = defaultdict(int)

    # Process each reaction: build listOfgeneList5 from cached getgpr/getLocation results
    for x in model.reactions:
        x2 = model2.reactions.get_by_id(x.id)
        annotation2 = x2.annotation if hasattr(x2, "annotation") else {}
        bounds = x2.bounds

        # collect ECs
        ec_value = annotation2.get("ec-code") if annotation2 else None
        # Warn early if the source reaction is an exchange and has EC annotation
        try:
            if annotation2 and ec_value:
                warn_if_ec_on_exchange(x2, annotation2, context="original")
        except Exception:
            LOGGER.debug(
                "Failed to check original reaction EC on exchange", exc_info=True
            )
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

            # Use explicit metabolite stoichiometries from the Reaction object
            # instead of fragile string parsing. This avoids incorrect coefficients
            # (e.g., exchange reactions being parsed incorrectly as coefficient 2.0).
            species = [s.id for s in x2.reactants] + [s.id for s in x2.products]
            species_coeffs = {m.id: coeff for m, coeff in x2.metabolites.items()}

            for CSL in [*geneList5[0]]:
                CSL2 = CSL
                biocyc_compartments_seen[CSL2] += 1

                # Use pre-computed compartment mapping
                CSL, mapping_method = map_compartment_once(CSL2)
                mapped_compartments[mapping_method] += 1

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

                        # add stoichiometry using original metabolite id mapping
                        # default to 1 if unexpected
                        stoich = species_coeffs.get(species2[1], 1)
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

                    # Preserve annotations but warn if this looks like an exchange reaction
                    annotation_copy = annotation2.copy() if annotation2 else {}
                    warn_if_ec_on_exchange(x2, annotation_copy, context="original->new")
                    reaction2.annotation = annotation_copy

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
                                # Only set ensembl annotation if we have valid data (not None)
                                if ensemble and isinstance(ensemble, dict):
                                    ensembl_id = ensemble.get("ensembl")
                                    if ensembl_id:
                                        model2.genes.get_by_id(gene).annotation[
                                            "ensembl"
                                        ] = ensembl_id
                                elif ensemble:
                                    # ensemble is a non-dict value (string, etc.)
                                    model2.genes.get_by_id(gene).annotation[
                                        "ensembl"
                                    ] = ensemble
                                # If ensemble is None, don't set the annotation

                                try:
                                    hgnc_symbol = geneList5[3].get(gene, gene)
                                    if hgnc_symbol:
                                        model2.genes.get_by_id(gene).annotation[
                                            "hgnc.symbol"
                                        ] = hgnc_symbol
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
                        existing_rxn = model2.reactions.get_by_id(x2.id)
                        # If this is an exchange/boundary reaction, make sure we don't keep ec-code
                        is_exchange_existing = (
                            len(existing_rxn.reactants) == 0
                            or len(existing_rxn.products) == 0
                            or len(species) == 1
                            or str(existing_rxn.id).lower().startswith("ex")
                        )
                        # Warn if existing exchange reaction has EC annotation but do not alter it
                        if is_exchange_existing and existing_rxn.annotation:
                            warn_if_ec_on_exchange(
                                existing_rxn,
                                existing_rxn.annotation,
                                context="existing",
                            )

                        existing_rxn.gene_reaction_rule = geneList5[2].get(
                            CSL2, existing_rxn.gene_reaction_rule
                        )
                        comp = list(model2.reactions.get_by_id(x2.id).compartments)
                        compartment_name = (
                            comp_dict.get(comp[0], comp[0]) if comp else ""
                        )
                        existing_rxn.annotation["sGPR"] = geneList5[0].get(
                            compartment_name, ""
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

    # Log compartment mapping summary
    LOGGER.info("=" * 60)
    LOGGER.info("Compartment Mapping Summary")
    LOGGER.info("=" * 60)
    LOGGER.info("BioCyc compartments seen (%d unique):", len(biocyc_compartments_seen))
    for comp, count in sorted(biocyc_compartments_seen.items(), key=lambda x: -x[1])[
        :20
    ]:
        LOGGER.info("  %s: %d times", comp, count)
    if len(biocyc_compartments_seen) > 20:
        LOGGER.info("  ... and %d more", len(biocyc_compartments_seen) - 20)

    LOGGER.info("\nCompartment mappings applied (%d unique):", len(mapped_compartments))
    for mapping, count in sorted(mapped_compartments.items(), key=lambda x: -x[1])[:30]:
        LOGGER.info("  %s: %d times", mapping, count)
    if len(mapped_compartments) > 30:
        LOGGER.info("  ... and %d more", len(mapped_compartments) - 30)
    LOGGER.info("=" * 60)

    # Before saving the model, put a name to the newly added compartments
    LOGGER.info(
        "Updating compartment names. Before: %d compartments", len(model2.compartments)
    )
    LOGGER.info("Compartments: %s", sorted(model2.compartments.keys()))
    model2 = update_comp_names_bm(model2, comp_abb_file)
    LOGGER.info("After: %d compartments", len(model2.compartments))
    LOGGER.info("Compartments: %s", sorted(model2.compartments.keys()))

    # Evaluate and correct the model
    model4 = copy.deepcopy(model2)

    ## Check for duplicate reactions and eliminate
    reactions_and_metabolites = [
        (x.id, str([(y.id, x.metabolites[y]) for y in x.metabolites]))
        for x in [x for x in model4.reactions]
    ]

    pattern = [str(x[1]) for x in reactions_and_metabolites]
    pattern = [x for x in pattern if pattern.count(x) > 1]
    pattern = list(set(pattern))

    threshold = len(model.reactions)
    duplicate_reactions = []
    removed_reactions = []
    for count, x in enumerate(pattern):
        # print(count, "/", len(pattern))
        index = [
            i for i, y in enumerate([n[1] for n in reactions_and_metabolites]) if y == x
        ]

        if len(index) > 1:  # duplicate reactions
            duplicate_reactions.append(index)
            reactions_in_original_model = [x for x in index if x < threshold]

            if (
                len(reactions_in_original_model) < 1
            ):  # if no reaction in the original model, then keep the reaction with the lowest index
                index = index[1 : len(index)]
            else:  # otherwise, keep the reactions in the original model and delete the newly added reactions
                index = [x for x in index if x > threshold]

            for y in index:
                for g in model4.groups:  # remove removed reactions from groups
                    if reactions_and_metabolites[y][0] in [m.id for m in g.members]:
                        g.members.remove(
                            model4.reactions.get_by_id(reactions_and_metabolites[y][0])
                        )
                model4.reactions.remove(
                    model4.reactions.get_by_id(reactions_and_metabolites[y][0])
                )
                removed_reactions.append(reactions_and_metabolites[y][0])

    ## Check for unconnected metabolites
    metabolites_in_model = [x.id for x in model4.metabolites]
    metabolites_in_reactions = [
        [y.id for y in x.metabolites.keys()] for x in model4.reactions
    ]

    metabolites_in_reactions = []
    for x in model4.reactions:
        for y in x.metabolites:
            metabolites_in_reactions.append(y.id)
    metabolites_in_reactions = list(set(metabolites_in_reactions))

    isolated_metabolites = [
        x for x in metabolites_in_model if not x in metabolites_in_reactions
    ]
    for x in isolated_metabolites:
        model4.metabolites.remove(x)

    ## Check for unbalanced reactions
    list_reaction_balance = []
    for x in tqdm(model4.reactions, desc="Checking reaction balance"):
        # print(x.id)
        eq = x.reaction
        eq_test = re.sub("<=>", "->", re.sub("-->", "->", eq)).strip()
        eq = re.sub(
            "(^| )[0-9\.]+", "", re.sub("<=>", "->", re.sub("-->", "->", eq))
        ).strip()
        species = [s.id for s in x.reactants] + [s.id for s in x.products]
        try:
            all_met_have_formula_test = min(
                [1 if ListOfMetFrom[x] else 0 for x in species]
            )
        except:
            all_met_have_formula_test = 0
        if all_met_have_formula_test == 1:
            FromList = list()
            for y in species:
                eq_test = re.sub(y, ListOfMetFrom[y], eq_test)
                eq = re.sub(y, ListOfMetFrom[y], eq)

            if 1 < len(species) < 25:
                try:
                    mass_balance_test = test_reaction_balance(eq_test)
                    if not mass_balance_test[0]:
                        MB = mass_balance(eq, "R")

                        # check if the reaction is in the original model
                        if x.id in [r.id for r in model.reactions]:
                            original_model_reaction = model.reactions.get_by_id(
                                x.id
                            ).reaction
                        else:
                            original_model_reaction = ""

                        list_reaction_balance.append(
                            [x.id, mass_balance_test, MB, original_model_reaction]
                        )
                except Exception as e:
                    # Log the error but continue processing other reactions
                    LOGGER.warning(
                        "Failed to check mass balance for reaction %s: %s", x.id, str(e)
                    )
                    # Optionally add to error list with exception info
                    list_reaction_balance.append(
                        [x.id, (False, f"Error: {str(e)}"), None, ""]
                    )

    # Save unbalanced reactions to error file
    if list_reaction_balance:
        LOGGER.info(
            "Found %d unbalanced reactions, writing to %s",
            len(list_reaction_balance),
            OUTPUT_ERRORS,
        )
        try:
            with open(OUTPUT_ERRORS, "w") as f:
                f.write("Reaction_ID\tBalance_Test\tMass_Balance\tOriginal_Reaction\n")
                for rxn_data in list_reaction_balance:
                    rxn_id, balance_test, mb, orig = rxn_data
                    f.write(f"{rxn_id}\t{balance_test}\t{mb}\t{orig}\n")
        except Exception:
            LOGGER.exception("Failed to write error file %s", OUTPUT_ERRORS)
    else:
        LOGGER.info("All reactions are balanced!")

    # Write modified model
    LOGGER.info("Writing output model to %s", OUTPUT_MODEL_FINAL)
    cobra.io.write_sbml_model(model4, OUTPUT_MODEL_FINAL)
    LOGGER.info("Done. Wrote %s", OUTPUT_MODEL_FINAL)

    # Save any EC-on-exchange warnings to CSV for review
    if exchange_ec_warnings:
        try:
            warn_path = os.path.join(
                project_root, "files", f"exchange_ec_warnings{_date_tag}.csv"
            )
            df_warn = pd.DataFrame(exchange_ec_warnings)
            df_warn.to_csv(warn_path, index=False)
            LOGGER.info(
                "Saved %d EC-on-exchange warnings to %s",
                len(exchange_ec_warnings),
                warn_path,
            )
        except Exception:
            LOGGER.exception("Failed to save EC-on-exchange warnings CSV")


if __name__ == "__main__":
    main()
