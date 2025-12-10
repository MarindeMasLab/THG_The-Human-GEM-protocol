from cobra.io import read_sbml_model, write_sbml_model
import cobra
import os
import urllib.request, urllib.error, urllib.parse
import re
import urllib.request, urllib.parse, urllib.error
import requests
import copy
import time
import traceback
import itertools
import pubchempy as pcp
import string
import pickle
from collections import defaultdict
from itertools import zip_longest
from cobra import Model, Reaction, Metabolite
from collections import ChainMap
import dill
import pdb
from typing import List, Optional, Tuple
import sys
import logging

LOGGER = logging.getLogger(__name__)

# import the functions
from functions.gpr.gpr_def import getGPR, setup_biocyc_session

from functions.equations_bm_gdb import *


# retrieve website with function from getgpr


def get_html(request_url: str, session: Optional[requests.Session] = None) -> str:
    """Fetch an html by perfoming a GET HTTPS request, maybe with session."""
    if session is not None:
        return session.get(request_url, timeout=30).text
    else:
        return requests.get(request_url, timeout=30).text


def create_dict(gene, e):
    n = ""
    global my_dict
    if "my_dict" not in globals():
        my_dict = {}
    elif gene and e:
        my_dict[gene] = e
    elif gene and not e:
        n = my_dict[gene]
    return my_dict, n


# create a function that takes a string and a dictionary and replaces the keys in the string with the corresponding values in the dictionary
def multiple_replace_new(dictionary, text):
    # Remove any leading or trailing whitespace from the text
    text = text.strip()
    # put everything in lowercase
    text = text.lower()
    # Search for the text in the dictionary
    for key in dictionary:
        # If the key is found in the text, replace it with the corresponding value
        if key == text:
            text = text.replace(key, dictionary[key])
        if key != text:
            text = text
    return text


def normalize_synonyms_only(compartment_raw):
    """Normalize only synonym compartment names without mapping to restricted list.
    
    This is used in unrestricted mode to merge synonyms (e.g., cytoplasm→cytosol)
    while keeping distinct compartments separate (e.g., "endoplasmic reticulum" 
    vs "endoplasmic reticulum membrane").
    
    Args:
        compartment_raw: Raw compartment string from BioCyc/database
        
    Returns:
        Normalized compartment name (lowercase) with synonyms merged
    """
    if not compartment_raw or not isinstance(compartment_raw, str):
        return "cytosol"
    
    comp_lower = compartment_raw.strip().lower()
    
    # Direct synonym mappings (only merge true synonyms, not related compartments)
    SYNONYM_MAP = {
        # Cytoplasm/Cytosol are the same
        "cytoplasm": "cytosol",
        "intracellular": "cytosol",
        
        # Nuclear variations (but keep specific sub-compartments separate)
        "nuclear": "nucleus",
        "nuclei": "nucleus",
        
        # Mitochondria/Mitochondrion (singular/plural)
        "mitochondrion": "mitochondria",
        
        # ER variations (but keep "endoplasmic reticulum membrane" separate)
        "er": "endoplasmic reticulum",
        "sarcoplasmic reticulum": "endoplasmic reticulum",
        
        # Golgi variations
        "golgi complex": "golgi apparatus",
        "golgi body": "golgi apparatus",
        
        # Plasma membrane variants (true synonyms)
        "plasma membrane": "cell membrane",
        
        # Extracellular variants
        "extracellular space": "extracellular",
        "extracellular region": "extracellular",
        
        # Peroxisome plural/singular
        "peroxisomes": "peroxisome",
        
        # Lysosome plural/singular
        "lysosomes": "lysosome",
        "lysosomal": "lysosome",
        
        # Cytoskeleton variants
        "cytoskeletal": "cytoskeleton",
        
        # Glycocalyx variants
        "glycocalix": "glycocalix",
        "cell coat": "glycocalix",
        
        # Vesicle variations (keep specific types separate unless clearly synonyms)
        "vesicles": "vesicle",
        
        # Mitochondrial sub-compartments (specific terms to normalize)
        "mitochondrial outer membrane": "mitochondrial membrane",
        
        # Membrane generic (map to cell membrane only if no other qualifier)
        # Note: Do NOT map "mitochondrial membrane", "ER membrane" etc here
    }
    
    # Check for exact match in synonym map
    if comp_lower in SYNONYM_MAP:
        return SYNONYM_MAP[comp_lower]
    
    # Check for patterns that need normalization but aren't exact matches
    # Only normalize if it's JUST "membrane" with no other qualifier
    if comp_lower == "membrane":
        return "cell membrane"
    
    # Return as-is (keeps "endoplasmic reticulum membrane", "nucleoplasm", etc. separate)
    return comp_lower


def normalize_compartment_name(compartment_raw, standard_compartments_dict, biocyc_to_standard_mapping=None):
    """Normalize a BioCyc/raw compartment name to standard compartment using fuzzy matching.
    
    This implements a robust 3-tier strategy:
    1. Direct match in BioCyc-to-Standard mapping (column A → column C from Excel)
    2. Fuzzy match against BioCyc names, then map to standard
    3. Fuzzy match against standard compartments
    4. Default to cytosol if no match
    
    Args:
        compartment_raw: Raw compartment string from BioCyc/database
        standard_compartments_dict: Dictionary of valid compartments (from bb.pickle)
        biocyc_to_standard_mapping: Optional dict mapping BioCyc names to standard names (column 0 → column 2)
        
    Returns:
        Normalized compartment name (lowercase)
    """
    if not compartment_raw or not isinstance(compartment_raw, str):
        return "cytosol"
    
    # Build case-insensitive lookup of standard compartments
    std_lower = {k.lower(): k.lower() for k in standard_compartments_dict.values()}
    
    comp_lower = compartment_raw.strip().lower()
    
    # Tier 1: Direct match in BioCyc-to-Standard mapping (Excel column A → C)
    if biocyc_to_standard_mapping and comp_lower in biocyc_to_standard_mapping:
        mapped = biocyc_to_standard_mapping[comp_lower]
        if mapped in std_lower:
            return std_lower[mapped]
        return mapped  # Return even if not in std_lower (trust the Excel mapping)
    
    # Tier 2: Fuzzy match against BioCyc source names, then use mapping
    # BUT: Exclude generic terms that are too short or could match incorrectly
    # Special case: Handle mitochondrial compartments BEFORE generic "membrane" matching
    if "mitochondri" in comp_lower or "respiratory chain" in comp_lower:
        LOGGER.debug(f"Tier 2.5 mitochondria: Prioritizing mitochondrial check for '{compartment_raw}'")
        # This will be handled properly in Tier 4, but we need to skip Tier 2 fuzzy matching
        # that might incorrectly match "membrane" → "cell membrane"
        pass  # Skip to Tier 3/4
    elif biocyc_to_standard_mapping:
        # Try fuzzy matching against BioCyc names (keys in the mapping)
        # Require longer matches to avoid generic "membrane" matching
        for biocyc_name, standard_name in biocyc_to_standard_mapping.items():
            # Check if the raw compartment contains the BioCyc name pattern
            # Use minimum length of 15 to avoid generic terms like "membrane"
            if len(biocyc_name) > 14 and biocyc_name in comp_lower:
                if standard_name in std_lower:
                    LOGGER.debug(f"Tier 2: Fuzzy BioCyc match: '{compartment_raw}' → '{biocyc_name}' → '{standard_name}'")
                    return std_lower[standard_name]
                return standard_name
    
    # Tier 3: Direct match in standard compartments
    if comp_lower in std_lower:
        return std_lower[comp_lower]
    
    # Tier 4: Fuzzy matching for common BioCyc patterns (fallback)
    # Get standard compartments values as lowercase list for case-insensitive matching
    std_values_lower = [v.lower() for v in standard_compartments_dict.values()] if standard_compartments_dict else []
    
    # Check for mitochondria-related compartments (most specific first)
    # This MUST be checked before general membrane check to avoid "mitochondrial membrane" → "cell membrane"
    if "mitochondri" in comp_lower or "respiratory chain" in comp_lower:
        LOGGER.info(f"Tier 4 mitochondria: Checking '{compartment_raw}'")
        # Try to match to inner mitochondria for membrane/matrix/inner/respiratory chain compartments
        inner_keywords = ["inner", "matrix", "lumen", "cristae", "intermembrane", "respiratory chain"]
        matching_inner = any(keyword in comp_lower for keyword in inner_keywords)
        if matching_inner:
            has_inner_mito = "inner mitochondria" in std_values_lower
            LOGGER.info(f"Tier 4 mitochondria inner: '{compartment_raw}' has inner keywords={matching_inner}, 'inner mitochondria' in std_values? {has_inner_mito}, std_values_lower={std_values_lower}")
            if has_inner_mito:
                LOGGER.info(f"Tier 4: '{compartment_raw}' → 'inner mitochondria' via mitochondrial keywords")
                return "inner mitochondria"
        # General mitochondria (fallback for other mitochondrial compartments including just "mitochondrial membrane")
        has_mito = any("mitochondria" in v for v in std_values_lower)
        LOGGER.info(f"Tier 4 mitochondria general: '{compartment_raw}' has 'mitochondria' in std_values? {has_mito}")
        if has_mito:
            LOGGER.info(f"Tier 4: '{compartment_raw}' → 'mitochondria' via mitochondrial keywords")
            return "mitochondria"
    
    # Check for nucleus-related compartments (including nucleolus)
    if "nucle" in comp_lower:
        if std_values_lower and any("nucleus" in v for v in std_values_lower):
            return "nucleus"
    
    # Check for ER-related compartments
    if "endoplasmic" in comp_lower or "endoplasmic-reticulum" in comp_lower:
        if any("endoplasmic reticulum" in v for v in std_values_lower):
            return "endoplasmic reticulum"
    
    # Check for Golgi-related compartments
    if "golgi" in comp_lower:
        if any("golgi apparatus" in v for v in std_values_lower):
            return "golgi apparatus"
    
    # Check for lysosome-related compartments
    if "lyso" in comp_lower:
        if any("lysosome" in v for v in std_values_lower):
            return "lysosome"
    
    # Check for peroxisome-related compartments
    if "peroxisom" in comp_lower:
        if any("peroxisome" in v for v in std_values_lower):
            return "peroxisome"
    
    # Check for plasma/cell membrane (but not mitochondrial/ER/Golgi membranes)
    if ("plasma membrane" in comp_lower or "cell membrane" in comp_lower or comp_lower == "membrane") and \
       "mitochondri" not in comp_lower and "endoplasmic" not in comp_lower and "golgi" not in comp_lower:
        if any("cell membrane" in v for v in std_values_lower):
            return "cell membrane"
    
    # Check for extracellular
    if "extracellular" in comp_lower:
        if any("extracellular" in v for v in std_values_lower):
            return "extracellular"
    
    # Check for vesicle-related compartments
    if any(keyword in comp_lower for keyword in ["vesicle", "endosome", "granule", "vacuole"]):
        if any("vesicle" in v for v in std_values_lower):
            return "vesicle"
    
    # Check for cytoskeleton
    if "cytoskeleton" in comp_lower or "cytoskeletal" in comp_lower:
        if any("cytoskeleton" in v for v in std_values_lower):
            return "cytoskeleton"
    
    # Check for cytosol/cytoplasm
    if "cytosol" in comp_lower or "cytoplasm" in comp_lower:
        if any("cytosol" in v for v in std_values_lower):
            return "cytosol"
    
    # Tier 5: Default to cytosol for unrecognized compartments
    LOGGER.debug(f"Compartment '{compartment_raw}' not recognized, defaulting to cytosol")
    return "cytosol"


def getLocationnew(
    gpr,
    genelist1,
    genelist2,
    impose_locations,
    location_dict_file,
    session: Optional[requests.Session] = None,
    ensembl_cache: Optional[dict] = None,
    excel_file: Optional[str] = None,
    compartments_sheet_name: Optional[str] = None,
    comp_dict: Optional[dict] = None,
    comp_abb_file: Optional[str] = None,
):
    """Finds subcellular location.
    Using a SGPR (or GPR), it identifies the corresponding cellular locations and adapts the location-specific SGPRs accordingly.
    i.e.:
        - input GPR: a*1 or b*1 (tipically output from getGPR function)
        - identified cellular location: cytosol, mitochondria
        - identified isoforms in the citosol: a
        - identified isoforms in the mitochondria: b
        - output (cellular location: SGPR): {cytosol: a*1} {mitochondria: b*1}

    Inputs
    ----------
    gpr: (s)gpr string , tipically output from getGPR function
    genelist1: list of gene names (i.e. ['ADH1B', 'ADH6', 'ADH4'])
    genelist2: list of gene biocyc IDs (i.e. ['HS08983', 'HS10600', 'HS06569'])
    impose_locations: boolean. 0: free search. 1: limit the cellular locations to a list of locations and all the locations that do not fit within this list are considered cytosol.
    location_dict_file: file name (".pickle" extension) with the list of compartments

    Output
    -------
    RuleLoc: {cellular location: SGPR}, i.e: {'Cytosol': '[CES2*1]'}
    RuleLoc2: {cellular location: GPR}, i.e: {'Cytosol': '[CES2]'}
    RuleLoc3: {cellular location: GPR with Ensembl IDs}, i.e: {'Cytosol': 'ENSG00000172831'}
    RuleLoc4:  {genes Ensemble ID: genes name}, i.e {'ENSG00000172831': ['CES2']})
    """
    if session is None:
        session = requests
    
    # Import gene class for location caching
    from functions.class_generate_database import gene as GeneClass
    
    try:
        gprgpr = ""
        OtherLocations = ["Other locations"]
        gpr2 = gpr
        gpr = re.sub(r"\*[0-9]+", "", gpr)
        urls0 = [
            x.replace("[", "")
            .replace("(", "")
            .replace("]", "")
            .replace(")", "")
            .replace(",", "")
            .replace(" ", "")
            for x in gpr2.split("or")
        ]
        hh = 1
        biocyc_mapping = None  # Initialize BioCyc-to-Standard mapping
        if impose_locations == 1:
            Var = open(location_dict_file, "rb")  # Endo1b_variables.pkl
            dictt = pickle.load(Var)
            hh = 0  # dict and not all
            # print(location_dict_file + " is used")
            
            # Try to load BioCyc-to-Endo1a mapping from parameters or hardcoded path
            try:
                from functions.function_bm_gdb import compartment_file_to_dict_bm
                import os
                import pandas as pd
                
                # Use parameters if provided, otherwise fall back to hardcoded paths
                if excel_file is None:
                    current_dir = os.path.dirname(os.path.abspath(__file__))
                    project_root = os.path.join(current_dir, "..", "..")
                    excel_file = os.path.join(project_root, "files", "ListOfCompartments_sept2024.xlsx")
                    comp_abb_file_path = os.path.join(project_root, "files", "compartments_info.txt")
                    compartments_sheet_name = "Def-Compartments"
                    comp_dict_local = {} if comp_dict is None else comp_dict
                else:
                    comp_abb_file_path = comp_abb_file if comp_abb_file else ""
                    comp_dict_local = {} if comp_dict is None else comp_dict
                    if not compartments_sheet_name:
                        compartments_sheet_name = "Def-Compartments"
                
                if os.path.exists(excel_file):
                    # Load BioCyc-to-Endo1a mapping
                    _, _, biocyc_mapping = compartment_file_to_dict_bm(
                        excel_file, compartments_sheet_name, comp_dict_local, comp_abb_file_path
                    )
                    LOGGER.info(f"Loaded BioCyc-to-Standard mapping with {len(biocyc_mapping)} entries from {excel_file}")
                else:
                    LOGGER.warning(f"Excel file not found: {excel_file}")
            except Exception as e:
                LOGGER.warning(f"Could not load BioCyc-to-Standard mapping: {e}")
                biocyc_mapping = None
        LocationList = []
        Locations = []
        uniprot = ""
        biocyc = ""
        n = 0

        for n in range(len(urls0)):  # isoforms
            a = urls0[n].split("and")
            m = 0
            d = ""
            for m in range(len(a)):

                uniprot = ""
                biocyc = ""

                p = 1
                b = (
                    "http://www.genome.jp/dbget-bin/www_bget?sp:"
                    + re.sub(r"\*[0-9]+", "", a[m])
                    + "_HUMAN"
                )  # location in genome net human
                # bb = str(getHtml(b, session))  # .decode('utf-8')
                try:
                    bb = str(urllib.request.urlopen(b, timeout=30).read())
                except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
                    bb = ""
                dd = re.findall("GO:[0-9]+.+?C:(.+?);", bb)

                if not dd:
                    ddd = re.search("(SUBCELLULAR LOCATION:.*)", "")
                    if ddd:
                        dd = [ddd.group(1)]
                if not dd or dd:
                    if genelist1:
                        if isinstance(genelist1, str):
                            genelist11 = re.findall(r"\[(.+?)\]", genelist1)
                            genelist11 = [
                                gene.replace("(", "")
                                .replace(")", "")
                                .replace("[", "")
                                .replace("]", "")
                                for gene in genelist11
                            ]
                        if not isinstance(genelist1, str):
                            genelist11 = genelist1
                        # Use a try-except to handle genes not in the list
                        try:
                            index = [x.upper() for x in genelist11].index(
                                re.sub(r"\*[0-9]+", "", a[m].upper())
                            )
                        except ValueError:
                            # Gene not found in list, skip this gene
                            gene_name_clean = re.sub(r"\*[0-9]+", "", a[m].upper())
                            LOGGER.warning(
                                f"Gene {gene_name_clean} not found in genelist11, skipping"
                            )
                            continue

                        # Get gene symbol and BioCyc ID for caching
                        gene_symbol = re.sub(r"\*[0-9]+", "", a[m].upper())
                        biocyc_id = genelist2[index].upper() if index < len(genelist2) else None
                        
                        # Check cache first
                        cached_locations = GeneClass.get_locations(gene_symbol, biocyc_id, session)
                        if cached_locations is not None:
                            # Use cached locations
                            LOGGER.debug(f"Using cached locations for {gene_symbol}: {cached_locations}")
                            SubUnLoc = cached_locations
                            Locations = Locations + cached_locations
                            # Skip database queries since we have cached data
                            d = sorted(set(SubUnLoc))
                            m = len(a)  # Stop processing this complex
                            break
                        
                        # Cache miss - proceed with database queries
                        SubUnLoc = []

                        # try to retrieve website with function from getgpr (authenticated BioCyc)

                        # for humancyc
                        humancyc = (
                            "https://biocyc.org/gene?orgid=HUMAN&id="
                            + genelist2[index].upper()
                        )
                        # for metacyc
                        metacyc = (
                            "http://biocyc.org/gene?orgid=META&id="
                            + genelist2[index].upper()
                        )
                        btry = get_html(humancyc, session)
                        bb_human = str(btry)
                        btry = get_html(metacyc, session)
                        bb_meta = str(btry)

                        # Extract the section between "Locations" and "Reactions"
                        location_section_human = re.search(
                            r"(?<=\nLocations)(.*?)(?=>\nReactions)",
                            bb_human,
                            re.DOTALL,
                        )
                        location_section_meta = re.search(
                            r"(?<=\nLocations)(.*?)(?=>\nReactions)", bb_meta, re.DOTALL
                        )

                        if location_section_human:
                            compartments = re.findall(
                                r"[\w\s-]+(?=\s*<a href)",
                                location_section_human.group(0),
                            )
                            # remove any \n in the compartments
                            compartments = [
                                compartment.replace("\n", "")
                                for compartment in compartments
                            ]
                            # remove any leading or trailing whitespace
                            compartments = [
                                compartment.strip() for compartment in compartments
                            ]
                            dd.extend(compartments)

                        elif location_section_meta:
                            compartments = re.findall(
                                r"[\w\s-]+(?=\s*<a href)",
                                location_section_human.group(0),
                            )
                            # remove any \n in the compartments
                            compartments = [
                                compartment.replace("\n", "")
                                for compartment in compartments
                            ]
                            # remove any leading or trailing whitespace
                            compartments = [
                                compartment.strip() for compartment in compartments
                            ]
                            dd.extend(compartments)

                        dd = list(set(dd))  # remove duplicates

                        # Early exit: if BioCyc found location data, skip UniProt
                        if dd:
                            gene_name_clean = re.sub(r"\*[0-9]+", "", a[m])
                            LOGGER.info(f"Gene {gene_name_clean} location found in BioCyc: {dd}")
                            # Process BioCyc locations immediately for caching
                            SubUnLoc = []
                            
                            if hh == 0:
                                # With compartment restrictions: normalize to standard compartments
                                std_comps = dictt if 'dictt' in locals() else {}
                                biocyc_map = biocyc_mapping if 'biocyc_mapping' in locals() else None
                                for loc in dd:
                                    if loc and loc not in OtherLocations and not "GO" in loc and not "PubMed" in loc:
                                        # Normalize using fuzzy matching (improved with BioCyc-to-Standard mapping)
                                        normalized_loc = normalize_compartment_name(loc, std_comps, biocyc_map)
                                        SubUnLoc.append(normalized_loc)
                            else:
                                # Without compartment restrictions: keep raw BioCyc compartment names
                                # but normalize ALL synonyms (cytoplasm→cytosol, nuclear→nucleus, etc.)
                                for loc in dd:
                                    if loc and loc not in OtherLocations and not "GO" in loc and not "PubMed" in loc:
                                        # Normalize only synonyms, keep distinct compartments separate
                                        normalized_loc = normalize_synonyms_only(loc)
                                        SubUnLoc.append(normalized_loc)
                            
                            # Cache the BioCyc locations (normalized or raw depending on hh)
                            # Deduplicate the list (e.g., after cytoplasm→cytosol merge)
                            SubUnLoc = sorted(set(SubUnLoc))
                            if SubUnLoc and gene_symbol:
                                GeneClass.cache_location(gene_symbol, SubUnLoc, biocyc_id)
                                LOGGER.debug(f"Cached locations for {gene_symbol}: {SubUnLoc}")
                            
                            # Add to global Locations list
                            Locations = Locations + SubUnLoc
                            d = sorted(set(SubUnLoc))
                            m = len(a)  # Stop processing this complex
                            break  # Exit gene loop since we have locations

                        # Only query UniProt if no location data found in BioCyc
                        if not dd:
                            UniProtKB = ""
                            GeneID = re.sub(r"\*[0-9-]+", "", a[m])
                            https = (
                                "https://www.uniprot.org/uniprot/?query="
                                + GeneID
                                + "&sort=score"
                            )
                            try:
                                url = str(urllib.request.urlopen(https, timeout=30).read())
                            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
                                url = ""
                            if re.search(
                                'uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>'
                                + GeneID
                                + "_HUMAN",
                                url,
                            ):
                                UniProtKB = re.search(
                                    'uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>'
                                    + GeneID
                                    + "_HUMAN",
                                    url,
                                ).group(1)
                            elif re.search(
                                '<a href="\/uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>[A-Z0-9-]+_HUMAN.+?<div class="gene-names"><span class="shortName">(<strong>[a-zA-Z0-9-]+<\/strong>[a-zA-Z0-9-, ]*)<\/span><\/div><\/td><td>?',
                                url,
                                re.IGNORECASE,
                            ):
                                search = re.search(
                                    '<a href="\/uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>[A-Z0-9-]+_HUMAN.+?<div class="gene-names"><span class="shortName">(<strong>[a-zA-Z0-9-]+<\/strong>[a-zA-Z0-9-, ]*)<\/span><\/div><\/td><td>?',
                                    url,
                                    re.IGNORECASE,
                                )
                                if re.findall(GeneID, search.group(2), re.IGNORECASE):
                                    UniProtKB = search.group(1)
                            elif re.search(
                                '<a href="\/uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>[A-Z0-9-]+_MOUSE.+?<div class="gene-names"><span class="shortName">(<strong>[a-zA-Z0-9-]+<\/strong>[a-zA-Z0-9-, ]*)<\/span><\/div><\/td><td>?',
                                url,
                                re.IGNORECASE,
                            ):
                                search = re.search(
                                    '<a href="\/uniprot\/([A-Z0-9-]+)">[A-Z0-9]+<\/a><\/td><td>[A-Z0-9-]+_MOUSE.+?<div class="gene-names"><span class="shortName">(<strong>[a-zA-Z0-9-]+<\/strong>[a-zA-Z0-9-, ]*)<\/span><\/div><\/td><td>?',
                                    url,
                                    re.IGNORECASE,
                                )
                                if re.findall(GeneID, search.group(2), re.IGNORECASE)[
                                    0
                                ]:
                                    UniProtKB = search.group(1)
                            cc = UniProtKB
                            b = (
                                "https://www.uniprot.org/uniprot/"
                                + cc
                                + "#subcellular_location"
                            )
                            try:
                                bb = str(urllib.request.urlopen(b, timeout=30).read())
                            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
                                bb = ""
                            ddd = re.findall(
                                'class="[a-zA-Z_ ]+"><h6>([a-zA-Z ]+)</h6>', bb
                            )
                            if ddd == ["Other locations"]:
                                ddd = re.findall(
                                    'locations*/SL-[0-9]+">([a-zA-Z ]+) </a>', bb
                                )
                            dd.extend(ddd)
                            uniprot = "1"

                dd = [i.strip() for i in dd]
                dd = [d for d in dd if d]

                if dd:
                    ee = [
                        re.sub("}[.,;]", "}@", x, flags=re.DOTALL)
                        .replace("  ", "")
                        .replace("@ ", "")
                        .replace(" {", "{")
                        .replace("SUBCELLULAR LOCATION: ", "")
                        .split("}")[0]
                        for x in dd
                    ]
                    if uniprot or biocyc:
                        ee = re.findall(r"([A-Za-z ]+)", ee[0])
                        uniprot = ""
                    SubUnLoc = []
                    x = 0

                    for x in range(len(dd)):
                        if (
                            dd[x]
                            and not dd[x] in OtherLocations
                            and not "GO" in dd[x]
                            and not "PubMed" in dd[x]
                        ):
                            ee = dd
                            ff = ee[x].split("{")[0]
                            if biocyc:
                                ff = re.findall(r"([A-Za-z0-9 ,;\-\(\)]+)", ff)[0]
                                if re.findall(r",", ff):
                                    ee = ee + ff.split(",")[1:]
                                    ff = ff.split(",")[0]
                                biocyc = ""
                            if re.findall(
                                ":", ff
                            ):  # if there is ":" in the subcellular location, the real location is after :
                                a = ff.split(":")
                                ff = a[1]
                            if hh == 0:

                                ff = ff.lower()
                                ff = multiple_replace_new(
                                    dictt, ff.split("{")[0]
                                )  # replace the keys in the string with the corresponding values in the dictionary

                            if re.findall("Dendriti", ff, re.IGNORECASE):
                                ff = "Dendrite"
                            if (
                                re.match("membrane", ff, re.IGNORECASE)
                                or re.findall("plasma membrane", ff, re.IGNORECASE)
                                or re.findall(
                                    "integral component of membrane", ff, re.IGNORECASE
                                )
                            ):
                                ff = "cell membrane"
                            if re.findall("eroxisom", ff):
                                ff = "Peroxisome"
                            if re.findall("itochondri", ff) or re.findall(
                                "mitochondri", ff
                            ):
                                ff = "Mitochondria"
                            if re.findall("lysosom", ff) or re.findall("Lysosom", ff):
                                ff = "Lysosome"
                            if re.findall("olgi", ff):
                                ff = "Golgi apparatus"
                            if re.findall("xtracel", ff):
                                ff = "Extracellular"
                            if re.findall("ndoplasm", ff) or re.findall("ndosom", ff):
                                ff = "Endoplasmic reticulum"
                            if (
                                re.findall("ytosol", ff)
                                or re.findall("ytoplasm", ff)
                                or re.findall("ntracel", ff)
                            ):
                                ff = "Cytosol"
                            if (
                                re.findall("ucle", ff)
                                or re.findall("enter", ff)
                                or re.findall("entro", ff)
                                or re.findall("entri", ff)
                                or re.findall("pindle", ff)
                                or re.findall("RNA", ff)
                                or re.findall("DNA", ff)
                                or re.findall("SMN complex", ff)
                                or re.findall("RISC complex", ff)
                                or re.findall("axon", ff)
                                or re.findall("Axon", ff)
                            ):
                                ff = "Nucleus"
                            if hh == 1 and not ff.lower():  # important
                                ff = "Cytosol"
                            if hh == 0 and not ff.lower() in list(
                                set(dictt.values())
                            ):  # important
                                ff = "Cytosol"
                            if re.findall("\\\\n", ff) or re.findall(
                                "\.", ff
                            ):  # important
                                ff = "Cytosol"
                            if not re.findall("Note=", ff):
                                ff = ff.lower()
                                # Apply synonym normalization in unrestricted mode too
                                if hh == 1:
                                    ff = normalize_synonyms_only(ff)
                                SubUnLoc = SubUnLoc + [ff]
                                Locations = Locations + [ff]

                            # transform everything in lowercase

                    d = sorted(set(SubUnLoc))

                    # Cache the successfully queried locations
                    if SubUnLoc and gene_symbol:
                        GeneClass.cache_location(gene_symbol, SubUnLoc, biocyc_id)
                        LOGGER.debug(f"Cached locations for {gene_symbol}: {SubUnLoc}")

                    m = len(
                        a
                    )  # once it is defined a cellular location the process stops because is assumed that all the subunit of the same complex are in the same place
                if not dd:
                    m = m + 1
                    if not gpr in gprgpr:
                        gprgpr += gpr + "\n"
                    d = [
                        "cytosol"
                    ]  # by default, if there is not anotated location, the reaction is located into the cytosol
                    p = 0
                    Locations = Locations + [d][0]

            LocationList = LocationList + [d]
            n = n + 1
        Locations = list(set(Locations))
        l = 0
        RuleLoc = {}
        RuleLoc2 = {}
        RuleLoc3 = {}
        RuleLoc4 = {}
        while l < len(Locations):
            L = Locations[l]
            k = 0
            LocGPR = "["
            while k < len(LocationList):
                g = 0
                while g < len(LocationList[k]):
                    if L == LocationList[k][g]:
                        LocGPR = LocGPR + urls0[k] + "] or ["
                    g = g + 1
                k = k + 1
            LocGPR = LocGPR[:-5].replace("and", " and ").replace("  ", " ")
            RuleLoc[Locations[l]] = LocGPR
            RuleLoc2[Locations[l]] = re.sub("\*[0-9]+", "", LocGPR)
            LocGPR2 = re.sub("\*[0-9]+", "", LocGPR)
            LocGPR2 = re.sub("\[", "", LocGPR2)
            LocGPR2 = re.sub("\]", "", LocGPR2)
            m = LocGPR2.split(" ")[1::2]
            LocGPR2 = LocGPR2.split(" ")[::2]
            GPR6 = create_dict(str(), str())[0]
            for i in range(len(LocGPR2)):
                iiii = ""
                if LocGPR2[i] in GPR6:
                    iiii = create_dict(LocGPR2[i], str())[1]
                else:
                    if LocGPR2[i]:
                        iiii = LocGPR2[i]
                        iiii2 = ""
                        iiii2 = ""
                        try:
                            iiii2 = re.search(
                                "gene=([A-Z0-9]+)",
                                str(
                                    urllib.request.urlopen(
                                        "https://www.genome.jp/dbget-bin/www_bget?hsa+"
                                        + LocGPR2[i],
                                        timeout=30
                                    ).read()
                                ),
                            )  # .group(1)
                        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
                            iiii2 = None
                        if not iiii2:
                            try:
                                iiii2 = re.search(
                                    "(ENSG[0-9]+)",
                                    str(
                                        urllib.request.urlopen(
                                            "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g="
                                            + LocGPR2[i],
                                            timeout=30
                                        ).read()
                                    ),
                                )
                            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
                                iiii2 = None
                        if iiii2:
                            iiii = iiii2.group(1)
                        # If an external ensembl cache was provided, prefer it to avoid
                        # expensive web requests. The cache is expected to map gene
                        # identifiers (symbols or other ids) to Ensembl IDs or to a
                        # dictionary containing an 'ensembl' key.
                        if ensembl_cache is not None:
                            try:
                                if LocGPR2[i] in ensembl_cache:
                                    val = ensembl_cache[LocGPR2[i]]
                                    # fetch_ensembl_annotations returns either a dict of
                                    # annotation fields or a single Ensembl id string in
                                    # some contexts; handle both
                                    if isinstance(val, dict):
                                        maybe_ens = val.get("ensembl")
                                        if maybe_ens:
                                            iiii = maybe_ens
                                    elif isinstance(val, str):
                                        iiii = val
                            except Exception:
                                # fall back to existing lookup behaviour on any error
                                pass
                RuleLoc4[iiii] = list()
                RuleLoc4[iiii].append(LocGPR2[i])
                create_dict(LocGPR2[i], iiii)
                LocGPR2[i] = iiii
            LocGPR2 = " ".join(
                [m + " " + str(n) for m, n in zip_longest(LocGPR2, m, fillvalue="")]
            )[:-1]
            RuleLoc3[Locations[l]] = LocGPR2
            l = l + 1

        if not RuleLoc4 == {"": [""]}:
            LOGGER.debug(f"getLocationnew returning RuleLoc with {len(RuleLoc)} compartments: {list(RuleLoc.keys())}")
            return (
                RuleLoc,
                RuleLoc2,
                RuleLoc3,
                RuleLoc4,
            )  # , p # p indicates that the location couldn't be determined and citosol has been put instead

    except Exception as e:
        # return ""
        print("exception occurred")
        print(e)
        # On error return empty structures with the same shape as the happy-path return
        # so callers can safely iterate/extend the result instead of receiving None.
        # Returning empty dicts preserves the expected (RuleLoc, RuleLoc2, RuleLoc3, RuleLoc4)
        return ({}, {}, {}, {})
