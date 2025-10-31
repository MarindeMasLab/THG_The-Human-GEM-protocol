import difflib
import logging
import re
import time
import requests

import cobra
import pandas as pd
import pubchempy as pcp
import numpy as np
import pickle
from collections import defaultdict
from typing import List, Dict, Optional, Tuple
import os
import sys

# Determine the current file's directory and the project root.
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.join(current_dir, "..")

# Add the project root to sys.path to access top-level folders like 'functions' and 'models'
if project_root not in sys.path:
    sys.path.append(project_root)

from functions.function_reac_identification import *

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)

LOGGER = logging.getLogger(__name__)
H_PATTERN = re.compile(r"H[0-9]+")
CHARGE_PATTERN = re.compile(r"[-\+][0-9]+")

# Rate limiting: Track last request time to ensure max 5 requests per second
_last_request_time = 0
_min_request_interval = 0.2  # Minimum 0.2 seconds between API calls (max 5 req/sec)


def _rate_limit():
    """Ensure we don't exceed 5 requests per second to PubChem API."""
    global _last_request_time
    current_time = time.time()
    time_since_last = current_time - _last_request_time

    if time_since_last < _min_request_interval:
        sleep_time = _min_request_interval - time_since_last
        time.sleep(sleep_time)

    _last_request_time = time.time()


def setup_proxy():
    """
    Setup proxy for PubChem requests to avoid IP rate limiting.

    Options:
    1. Set HTTP_PROXY and HTTPS_PROXY environment variables
    2. Use a proxy service (requires subscription)
    3. Disable VPN and use home internet
    4. Set DISABLE_PROXY=1 to disable proxy usage

    To use a proxy, set environment variables before running:
        export HTTP_PROXY="http://proxy-server:port"
        export HTTPS_PROXY="http://proxy-server:port"

    To disable proxy:
        export DISABLE_PROXY=1

    Or uncomment and configure the proxy dict below.
    """
    # Check if proxy is disabled
    if os.environ.get("DISABLE_PROXY") == "1":
        print("Proxy disabled by DISABLE_PROXY environment variable")
        LOGGER.info("Proxy disabled by environment variable")
        return None

    proxy = None

    # Option 1: Read from environment variables
    http_proxy = os.environ.get("HTTP_PROXY") or os.environ.get("http_proxy")
    https_proxy = os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy")

    if http_proxy or https_proxy:
        proxy = {}
        if http_proxy:
            proxy["http"] = http_proxy
        if https_proxy:
            proxy["https"] = https_proxy
        print(f"Using proxy configuration: {proxy}")
        LOGGER.info(f"Using proxy: {proxy}")

    # Option 2: Manually configure proxy (DISABLED by default - proxy is broken)
    # Uncomment and update the proxy below if you have a working proxy
    # if not proxy:
    #     proxy = {
    #         "http": "http://your-proxy-server:port",
    #         "https": "http://your-proxy-server:port",
    #     }
    #     print(f"Using hardcoded proxy: {proxy}")
    #     LOGGER.info(f"Using hardcoded proxy: {proxy}")

    if not proxy:
        print("No proxy configured - using direct connection")
        LOGGER.info("No proxy configured - using direct connection")

    return proxy


# Setup proxy at module load
PROXY_CONFIG = setup_proxy()


def global_met_annotation_file():
    global annotation_file
    annotation_file = os.path.join(
        project_root, "metabolite_reac_identification", "reports", "met_annotation.tsv"
    )
    return annotation_file


def gather_metabolites(model: cobra.Model) -> List[Tuple[str, str, str, str]]:
    """Gather metabolites info from a `model`, filtering out from `identifiers`."""
    met_list = []
    for met in model.metabolites:
        met_tuple = (met.name, met.formula, met.annotation, met.id[:-1])
        if met_tuple not in met_list:
            met_list.append(met_tuple)
    return met_list


def formula_similarity(formula_left: str, formula_right: str) -> float:
    """Compute similiratiy based on the stdlib SequenceMatcher."""
    pruned_left = H_PATTERN.sub("H", CHARGE_PATTERN.sub("", atom(formula_left)))
    pruned_right = H_PATTERN.sub("H", CHARGE_PATTERN.sub("", atom(formula_right)))
    return difflib.SequenceMatcher(
        lambda x: x == " ", pruned_left, pruned_right
    ).ratio()


def identify_metabolite(
    name: str,
    formula: str,
    iden: str,
    threshold: float = 0.82,
    max_retries: int = 3,
    use_proxy: bool = False,
) -> Optional[str]:
    """Try to match metabolite info to a pubchem compound to get the annotation.

    Parameters
    ----------
    name: str
    formula: str
    iden: str
    threshold: float, default=0.82
        minimum similarity ratio between the molecular formulas. The default (0.82)
        was decided by performing sensibility tests of the metabolites in Human 1.
    max_retries: int, default=3
        maximum number of retry attempts when PubChem server is busy
    use_proxy: bool, default=False
        whether to use proxy configuration if available

    Returns
    -------
    met_result: Optional[str]
        Tab-separated str. If None, the metabolite was not identified.

    """
    # Add rate limiting to avoid PubChem API throttling
    time.sleep(0.2)  # 200ms delay between requests
    
    CompoundID = None
    met_result = None

    # Configure requests session with proxy if available
    if use_proxy and PROXY_CONFIG:
        import urllib.request

        if "http" in PROXY_CONFIG:
            os.environ["HTTP_PROXY"] = PROXY_CONFIG["http"]
        if "https" in PROXY_CONFIG:
            os.environ["HTTPS_PROXY"] = PROXY_CONFIG["https"]

    # Retry logic with exponential backoff for PubChem API calls
    retry_count = 0
    while retry_count < max_retries:
        try:
            CompoundID = pcp.get_cids(name.strip(), "name")
            break
        except Exception as e:
            if "503" in str(e) or "ServerBusy" in str(e):
                retry_count += 1
                if retry_count < max_retries:
                    wait_time = 2 ** retry_count  # Exponential backoff: 2, 4, 8 seconds
                    LOGGER.warning(f"PubChem busy, waiting {wait_time}s before retry {retry_count}/{max_retries}")
                    time.sleep(wait_time)
                else:
                    LOGGER.warning(f"Failed after {max_retries} retries for {name}")
                    return None
            else:
                LOGGER.warning(f"Error getting CID for {name}: {e}")
                return None

    # Retry with alternative name format if first attempt failed
    if not CompoundID:
        LOGGER.warning("get_cids did not work")
        time.sleep(0.2)
        try:
            CompoundID = pcp.get_cids(re.sub(r"[\(\)]", "", name), "name")
        except Exception as e:
            LOGGER.warning(f"get_cids extra error: {e}")
            return None
    
    # Fetch compound details if we have a CID
    if CompoundID:
        try:
            time.sleep(0.2)
            Compound = pcp.Compound.from_cid(CompoundID)
        except Exception as e:
            LOGGER.warning(f"Error getting compound from CID: {e}")
            return None
            
        molecular_formula = Compound.molecular_formula

        if formula_similarity(formula, molecular_formula) >= threshold:
            try:
                time.sleep(0.2)
                Compound = pcp.Compound.from_cid(CompoundID)
            except Exception as e:
                LOGGER.warning(f"Error getting detailed compound info: {e}")
                return None
            CompoundID = Compound.cid
            synonyms = Compound.synonyms
            Metabolite = sorted(
                [
                    [
                        difflib.SequenceMatcher(
                            None, name.lower(), synonym.lower()
                        ).ratio(),
                        synonym,
                    ]
                    for synonym in synonyms
                ]
            )[-1][1]
            Kegg = "".join(list(filter(re.compile("[CG][0-9]{5}$").match, synonyms)))
            LIPIDMAPSID = "".join(
                list(filter(re.compile("L[A-Z]{3}[0-9]+$").match, synonyms))
            )
            CHEBI = "".join(list(filter(re.compile("CHEBI:[0-9]+$").match, synonyms)))
            inchi = Compound.inchi
            inchikey = Compound.inchikey
            met_result = (
                Metabolite
                + "\t"
                + ""
                + "\t"
                + molecular_formula
                + "\t"
                + LIPIDMAPSID
                + "\t"
                + Kegg
                + "\t"
                + CHEBI
                + "\t"
                + str(CompoundID)
                + "\t"
                + ""
                + "\t"
                + inchikey
                + "\t"
                + inchi
                + "\t"
                + ""
                + "\t"
                + iden
                + "\n"
            )
    else:
        LOGGER.warning("get_cids extra did not work")
    return met_result


def generate_met_annotation(
    met_list: List[Tuple[str, str, str, str]],
    out: str = global_met_annotation_file(),
    delay_between_requests: float = 1.0,
    checkpoint_interval: int = 10,
    resume: bool = True,
) -> Tuple[List, List]:
    """Identify and write metabolite annotations to file (Tab-separated).

    Parameters
    ----------
    met_list: list[tuple[str, str, str, str]]
        list of metabolites to annotate
    out: str
        output file path
    delay_between_requests: float, default=1.0
        delay in seconds between PubChem API requests to avoid rate limiting
    checkpoint_interval: int, default=10
        save progress after every N metabolites
    resume: bool, default=True
        if True, resume from checkpoint file if it exists

    Returns
    -------
    unnanotated: list[tuple[str, str, str, str]]
        list of unnanotated metabolites
    anotated: list[tuple[str, str, str, str]]
        list of anotated metabolites
    """
    unnanotated, annotated = [], []
    api_failures = []  # Track API/temporary failures separately
    total = len(met_list)
    
    # Failure tracking file
    failure_file = out.replace('.tsv', '_failures.tsv')
    
    print(f"\n{'='*70}")
    print(f"Starting metabolite annotation for {total} metabolites")
    print(f"{'='*70}\n")

    with open(out, "w") as f, open(failure_file, "w") as fail_f:
        # Write header for failure file
        fail_f.write("name\tformula\tidentifier\tfailure_reason\n")
        
        for idx, (name, formula, annotation, iden) in enumerate(met_list, 1):
            if idx % 100 == 0:
                print(f"[{idx}/{total}] {(idx/total)*100:.1f}% | Success: {len(annotated)}, Failed: {len(unnanotated)}, API errors: {len(api_failures)}")
            
            try:
                # Track if this is an API failure
                api_error = False
                failure_reason = "unknown"
                
                result = identify_metabolite(name, formula, iden)
                met = [name, formula, annotation, iden]
                
                if result is None:
                    # Try to determine failure reason from recent logs
                    # Since we can't easily capture it, we'll mark for retry
                    unnanotated.append(met)
                    failure_reason = "not_found_or_api_error"
                    fail_f.write(f"{name}\t{formula}\t{iden}\t{failure_reason}\n")
                else:
                    f.write(result)
                    f.flush()  # Ensure data is written immediately
                    annotated.append(met)
                    
            except Exception as e:
                error_str = str(e)
                met = [name, formula, annotation, iden]
                unnanotated.append(met)
                
                # Categorize error
                if "503" in error_str or "ServerBusy" in error_str or "HTTP Error" in error_str:
                    failure_reason = "api_error"
                    api_failures.append(met)
                elif "timeout" in error_str.lower():
                    failure_reason = "timeout"
                    api_failures.append(met)
                else:
                    failure_reason = f"exception: {error_str[:50]}"
                
                fail_f.write(f"{name}\t{formula}\t{iden}\t{failure_reason}\n")
                continue
    
    print(f"\n{'='*70}")
    print(f"ANNOTATION COMPLETE")
    print(f"{'='*70}")
    print(f"Total metabolites:          {total}")
    print(f"Successfully annotated:     {len(annotated)} ({(len(annotated)/total)*100:.1f}%)")
    print(f"Failed (likely not in DB):  {len(unnanotated)-len(api_failures)} ({((len(unnanotated)-len(api_failures))/total)*100:.1f}%)")
    print(f"Failed (API/temp errors):   {len(api_failures)} ({(len(api_failures)/total)*100:.1f}%)")
    print(f"\nResults saved to:      {out}")
    print(f"Failures saved to:     {failure_file}")
    print(f"{'='*70}\n")
    
    return annotated, unnanotated


def remove_null_value(d):
    return {
        a: c
        for a, b in d.items()
        if (c := (b if not isinstance(b, dict) else remove_null_value(b)))
    }


def process_annotation(annotation_file: str = global_met_annotation_file()) -> Dict:
    """Generate metabolite annotation file.

    Parameters
    ----------
    annotation_file: str
        Tab-separated file, generated with `generate_met_annotation`
    """
    print(f"Processing annotation file: {annotation_file}")
    variableFile = pd.read_csv(annotation_file, sep="\t", header=None)
    variableFile[5] = variableFile[5].str.extract("(CHEBI:[0-9]+)", expand=True)
    variableFile = variableFile.rename(
        columns={
            3: "lipidmaps",
            4: "kegg.compound",
            5: "chebi",
            6: "pubchem.compound",
            8: "inchikey",
            9: "inchi",
        }
    )
    variableFile = variableFile.fillna(str())

    annotation = defaultdict(dict)

    for _, row in variableFile.iterrows():
        annotation[row[11]] = row.drop([0, 1, 2, 7, 10, 11]).to_dict()
    annotation = remove_null_value(annotation)
    return annotation


def atom(Formula):
    try:
        if "C" in Formula:
            C = re.findall(r"(C[0-9]+)", Formula)
            if not C:
                C = ["C1"]
        else:
            C = []
        if "H" in Formula:
            H = re.findall(r"(H\+*[0-9]+)", Formula)
            if not H:
                H = ["H1"]
        else:
            H = []
        if "O" in Formula:
            O = re.findall(r"(O[0-9]+)", Formula)
            if not O:
                O = ["O1"]
        else:
            O = []
        if "N" in Formula:
            N = re.findall(r"(N[0-9]+)", Formula)
            if not N:
                N = ["N1"]
        else:
            N = []
        if "P" in Formula:
            P = re.findall(r"(P[0-9]+)", Formula)
            if not P:
                P = ["P1"]
        else:
            P = []
        if "S" in Formula:
            S = re.findall(r"(Se?[0-9]+)", Formula)
            if not S:
                S = ["S1"]
        else:
            S = []
        if "I" in Formula:
            I = re.findall(r"(I[0-9]+)", Formula)
            if not I:
                I = ["I1"]
        else:
            I = []
        if "F" in Formula:
            F = re.findall(r"(Fe*[0-9]+)", Formula)
            if not F:
                F = ["F1"]
        else:
            F = []
        if "R" in Formula:
            R = re.findall(r"(R[0-9]+)", Formula)
            if not R:
                R = ["R1"]
        else:
            R = []
        return "".join(str(i) for i in C + H + O + N + P + S + I + F + R)
    except Exception as e:
        print(e)
