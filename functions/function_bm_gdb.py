# -*- coding: utf-8 -*-
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
import pandas as pd
import pdb
import os

from functions.equations_bm_gdb import *
from functions.class_generate_database import *
from functions.pattern_generate_database import *

from functions.gpr.auth_gpr import getGPR, setup_biocyc_session

# 👇 import the addtional function
from functions.gpr.get_location_def import getLocationnew as getLocation


# Optimized batch fetching functions for KEGG API
def parse_kegg_flat_file_to_html(flat_file_text, entry_id):
    """
    Convert KEGG REST API flat file format to HTML format expected by existing parsers.
    This is much more efficient than fetching HTML pages individually.

    Parameters:
    -----------
    flat_file_text : str
        Raw KEGG flat file format text from REST API
    entry_id : str
        KEGG entry ID (e.g., 'C00001')

    Returns:
    --------
    str : HTML formatted text compatible with existing HTML parsers
    """
    lines = flat_file_text.strip().split("\n")

    # Extract fields from flat file
    fields = {}
    current_field = None
    current_value = []

    for line in lines:
        # New field starts if line doesn't begin with whitespace
        if line and not line[0].isspace():
            # Save previous field
            if current_field:
                fields[current_field] = "\n".join(current_value)

            # Parse new field
            field_name = line[:12].strip()
            field_value = line[12:].strip()

            if field_name:
                current_field = field_name
                current_value = [field_value] if field_value else []
            else:
                current_field = None
                current_value = []
        else:
            # Continuation of current field
            if current_field:
                current_value.append(line.strip())

    # Don't forget last field
    if current_field:
        fields[current_field] = "\n".join(current_value)

    # Build HTML mimicking genome.jp structure
    html_parts = []

    # Title (required for entry ID detection)
    entry_type = "COMPOUND" if entry_id.startswith("C") else "GLYCAN"
    html_parts.append(f"<title>KEGG {entry_type}: {entry_id}</title>")

    # Name field (critical for parser)
    if "NAME" in fields:
        name_value = fields["NAME"].replace("\n", "; ").rstrip(";")
        html_parts.append(
            f"""Name</span></th>
<td class="td21 defd"><div class="cel"><div class="cel">{name_value}<br>
</div></div></td></tr>"""
        )

    # Formula field (embedded in special div structure)
    if "FORMULA" in fields:
        formula = fields["FORMULA"]
        html_parts.append(
            f'<div class="cel"><div class="cel">Formula [{formula}]</div></div>'
        )

    # Composition field (for glycans)
    if "COMPOSITION" in fields:
        composition = fields["COMPOSITION"]
        html_parts.append(
            f'<div class="cel"><div class="cel">Composition: {composition}<br>'
        )
        html_parts.append(f'<input type="hidden">{composition}<br>')

    # Remark field (for "Same as" links)
    if "REMARK" in fields:
        remark = fields["REMARK"]
        # Check for "Same as" pattern
        same_as_match = re.search(r"Same as:\s*([CDG][0-9]+)", remark)
        if same_as_match:
            same_as_id = same_as_match.group(1)
            html_parts.append(f'Remark\nSame as:\n<a href="#">{same_as_id}</a>')

    # DBLINKS field (for external database references)
    if "DBLINKS" in fields:
        dblinks_text = fields["DBLINKS"]

        # PubChem
        pubchem_matches = re.findall(r"PubChem:\s*(\d+)", dblinks_text)
        for pc_id in pubchem_matches:
            html_parts.append(f"PubChem: sid={pc_id}")

        # ChEBI
        chebi_matches = re.findall(r"ChEBI:\s*(\d+)", dblinks_text)
        for chebi_id in chebi_matches:
            html_parts.append(f'chebiId=CHEBI:{chebi_id}"')

        # LipidBank
        lipidbank_matches = re.findall(r"LipidBank:\s*([A-Z0-9]+)", dblinks_text)
        for lipid_id in lipidbank_matches:
            html_parts.append(f'LipidBank\n<a href="#" id="{lipid_id}"')

        # LIPID MAPS
        lipidmaps_matches = re.findall(r"LIPID MAPS:\s*([A-Z]+[0-9]+)", dblinks_text)
        for lm_id in lipidmaps_matches:
            html_parts.append(f"LMID={lm_id}")

        # GlycomeDB
        glycomedb_matches = re.findall(r"GlycomeDB:\s*([A-Z0-9]+)", dblinks_text)
        for glyc_id in glycomedb_matches:
            html_parts.append(f'GlycomeDB\n<a href="#?glycomeId={glyc_id}"')

        # JCGGDB (part of GlycomeDB entries)
        jcgg_matches = re.findall(r"(JCGG-[A-Z0-9]+)", dblinks_text)
        for jcgg_id in jcgg_matches:
            html_parts.append(f'GlycomeDB\n<a href="#">{jcgg_id}<')

    return "\n".join(html_parts)


def batch_fetch_kegg_entries(
    entry_ids, database="compound", batch_size=10, max_workers=5
):
    """
    Fetch multiple KEGG entries using REST API batch requests with concurrent execution.

    Combines two optimization strategies:
    1. Batch requests: Fetch up to 10 entries per HTTP request (KEGG API limit)
    2. Concurrent batches: Fetch multiple batches in parallel

    Parameters:
    -----------
    entry_ids : list
        List of KEGG entry IDs to fetch
    database : str
        Database type ('compound', 'reaction', 'glycan', etc.)
    batch_size : int
        Number of entries to fetch per batch request (max 10 for KEGG API)
    max_workers : int
        Maximum number of concurrent batch requests

    Returns:
    --------
    dict : Mapping of entry_id -> parsed page content (pseudo-HTML format)
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    results = {}

    # Split entry_ids into batches
    batches = []
    for i in range(0, len(entry_ids), batch_size):
        batches.append(entry_ids[i : i + batch_size])

    def fetch_batch(batch):
        """Fetch a single batch of entries"""
        try:
            # Build batch request with proper database prefix
            if database == "compound":
                batch_str = "+".join([f"cpd:{eid}" for eid in batch])
            elif database == "glycan":
                batch_str = "+".join([f"gl:{eid}" for eid in batch])
            elif database == "reaction":
                batch_str = "+".join([f"rn:{eid}" for eid in batch])
            else:
                batch_str = "+".join(batch)

            # Use KEGG REST API batch get
            url = f"https://rest.kegg.jp/get/{batch_str}"
            response = urllib.request.urlopen(url).read()

            # Handle both bytes and str
            content = (
                response.decode("utf-8") if isinstance(response, bytes) else response
            )

            # Parse the batched response into individual entries
            # KEGG separates entries with "///"
            entries = content.split("///")

            # Map entries back to IDs - return flat file format directly
            batch_results = {}
            for j, entry_id in enumerate(batch):
                if j < len(entries) and entries[j].strip():
                    # Return the flat file format directly
                    # The compound class will detect this and use the appropriate parser
                    batch_results[entry_id] = entries[j].strip()
                else:
                    batch_results[entry_id] = None

            return batch_results

        except Exception as e:
            print(f"Warning: Batch fetch failed for {batch}: {e}")
            import traceback

            print(traceback.format_exc())
            return {eid: None for eid in batch}

    # Fetch batches concurrently
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_batch = {
            executor.submit(fetch_batch, batch): batch for batch in batches
        }

        for future in as_completed(future_to_batch):
            batch_results = future.result()
            results.update(batch_results)

    return results


def compartment_file_to_dict_bm(
    excel_file, compartments_sheet_name, comp_dict, pickle_file
):
    # ComEquiv is the comp_dict but reversed (keys and values are swapped) -> compartment name: abbreviation
    ComEquiv = {v: k for k, v in comp_dict.items()}
    CompList = list(ComEquiv.values())
    Compartment_CL = list(
        set(
            sorted(
                [
                    x
                    for x in pd.read_excel(
                        excel_file,
                        sheet_name=compartments_sheet_name,
                        header=None,
                        skiprows=lambda x: x in [0, 1],
                    )[0]
                ]
            )
        )
    )
    Compartment_CL = [
        x.strip() for x in Compartment_CL
    ]  # remove leading and trailing whitespaces

    # Put everything in lower case
    CompList = [x.lower() for x in CompList]
    Compartment_CL = [x.lower() for x in Compartment_CL]

    for CSL in Compartment_CL:
        CompList = list(set(CompList))
        if CSL in ComEquiv:
            ID = ComEquiv[CSL]
        elif len(re.sub(" $", "", re.sub("^ ", "", CSL)).split(" ")) > 1:
            ID = (CSL.split(" ")[0][0] + CSL.split(" ")[1][0]).lower().replace(" ", "")
        else:
            if len(CSL.split(" ")) > 1:
                ID = CSL[0:3].lower().replace(" ", "")
            else:
                ID = CSL[0:2].lower().replace(" ", "")
        if not CSL in ComEquiv and ID in CompList:
            r = re.compile(ID)
            ID = ID + str(len(list(filter(r.match, CompList))) + 1)
        ComEquiv[CSL] = ""
        ComEquiv[CSL] += ID
        CompList.append(ID)

    # Original mapping: full compartment name -> abbreviation.
    name_to_abbr = {k: v for k, v in sorted(ComEquiv.items(), key=lambda item: item[1])}
    # New mapping: abbreviation -> full compartment name.
    abbr_to_name = {v: k for k, v in name_to_abbr.items()}

    # Save the dictionary as a pickle file
    pd.to_pickle(abbr_to_name, pickle_file)

    # Return both mappings.
    return name_to_abbr, abbr_to_name


# This function creates a dictionary with the abbreviations and the complete names of the compartments´
# Should be used only once to create the dictionary (in BM) and then the dictionary should be saved as a pickle file
def create_comp_abbreviations_dict_bm(excel_file, sheet_name, pickle_file):
    df = pd.read_excel(excel_file, sheet_name=sheet_name, header=None, index_col=None)
    # In this excel the first column is the complete compartment names and the second column is the abbreviations
    compartment_dict = dict(zip(df[1], df[0]))  # Abbreviation: Complete name
    pd.to_pickle(compartment_dict, pickle_file)


# Pickle dict with the abbreviations will be used in BM to update the compartment names and generate the sGPR annotations
def update_comp_names_bm(model, pickle_file):
    comp_dict = pd.read_pickle(pickle_file)
    model.compartments = {key: comp_dict.get(key, "") for key in model.compartments}
    return model


def create_compartments_dict_bm(excel_file, sheet_name, column_index, pickle_file_name):
    """
    Creates a dictionary mapping compartment names to their corresponding values in a model
    from an Excel file and saves it as a pickle file.

    Parameters:
    ----------
    excel_file : str
        The path to the Excel file containing compartment data.
    sheet_name : str
        The name of the sheet within the Excel file to read. In the Excel file, the first column should contain the compartment names.
    column_index : int
        The zero-based index of the column in the Excel file where mapping to your desired model compartments is located.
        For example, if the mapping of the compartment of the model are in the third column, pass `2`.

    Returns:
    -------
    dict
        A dictionary mapping compartment names (keys) to their corresponding names in the model (values).

    Notes:
    -----
    - The first column of the sheet is assumed to contain the compartment names.
    - The dictionary keys and values are converted to lowercase and stripped of whitespace.
    - Saves the dictionary as a pickle file in the current working directory.
    """
    # Read the Excel file and specified sheet
    compartments_df = pd.read_excel(excel_file, sheet_name=sheet_name)

    # Skip the first row and extract the required columns
    compartments_df = compartments_df.iloc[1:, [0, column_index]]

    # Create the dictionary
    dict_compartments = {}
    for i in range(len(compartments_df)):
        key = compartments_df.iloc[i, 0].lower().strip()
        value = compartments_df.iloc[i, 1].lower().strip()
        dict_compartments[key] = value

    # Save the dictionary to a pickle file
    with open(pickle_file_name, "wb") as f:
        pickle.dump(dict_compartments, f)

    # Return the dictionary for further use
    return dict_compartments


def compartment_file_to_dict():
    # Determine the current file's directory and the project root.
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.join(current_dir, "..")
    excel_file_path = os.path.join(project_root, "files", "ListOfCompartments.xlsx")

    ComEquiv = {
        "Glycosylphosphatidylinositol-N-acetylglucosaminyltransferase (GPI-GnT) complex": "gc",
        "P-body": "pb",
        "Extracellular": "e",
        "Peroxisome": "x",
        "Mitochondria": "m",
        "Cytosol": "c",
        "Lysosome": "l",
        "Endoplasmic reticulum": "r",
        "Golgi apparatus": "g",
        "Nucleus": "n",
        "Inner mitochondria": "i",
    }  # to keep the consistency between the DB and the initial compartments in H1
    CompList = list(ComEquiv.values())
    Compartment_CL = list(
        set(
            sorted(
                [
                    x[0].upper() + x[1:]
                    for x in pd.read_excel(
                        excel_file_path,
                        sheet_name=2,
                        header=None,
                        skiprows=lambda x: x in [0, 1],
                    )[0]
                ]
            )
        )
    )
    Compartment_CL = [x.strip() for x in Compartment_CL]

    for CSL in Compartment_CL:
        CSL = CSL[0].upper() + CSL[1:]
        CompList = list(set(CompList))
        if CSL in ComEquiv:
            ID = ComEquiv[CSL]
        elif len(re.sub(" $", "", re.sub("^ ", "", CSL)).split(" ")) > 1:
            ID = (CSL.split(" ")[0][0] + CSL.split(" ")[1][0]).lower().replace(" ", "")
        else:
            if len(CSL.split(" ")) > 1:
                ID = CSL[0:3].lower().replace(" ", "")
            else:
                ID = CSL[0:2].lower().replace(" ", "")
        if not CSL in ComEquiv and ID in CompList:
            r = re.compile(ID)
            ID = ID + str(len(list(filter(r.match, CompList))) + 1)
        ComEquiv[CSL] = ""
        ComEquiv[CSL] += ID
        CompList.append(ID)
    return {k: v for k, v in sorted(ComEquiv.items(), key=lambda item: item[1])}


def recondict():
    # Determine the current file's directory and the project root.
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.join(current_dir, "..")
    pathTo = os.path.join(project_root, "models", "Human-GEM_2022-06-21.xml")
    with open(pathTo) as m:
        m = m.read()
        defaultdictt = defaultdict(list)
        for species in re.findall(" +<species metaid.+?<.species>", m, re.DOTALL):
            try:
                iii = re.findall(
                    'fbc:charge="(.+?)" fbc:chemicalFormula="(.+?)">',
                    species,
                    re.DOTALL,
                )[0]
                i, ii = iii[0], iii[1]
                for res in re.findall(
                    '[a-z]+/[a-zA-Z0-9:/.]+[:/]([a-zA-Z0-9]+?)"/>', species
                ):
                    if not res in defaultdictt:
                        defaultdictt[res] = [ii, i]
            except Exception:
                continue
        return defaultdictt


def DefEnsblDB(Ensbl):
    with open(Ensbl) as DB:
        DB = DB.read()
        for line in DB.split("\n"):
            try:
                GeneAss = re.search("(.+?)_HUMAN", line.split()[0]).group(1)
                Ensemble = re.search("(ENSG\d+)", line.split()[1]).group(1)
                if not GeneAss in GeneEnsbl:
                    GeneEnsbl[GeneAss] = ""
                if GeneAss in GeneEnsbl:
                    GeneEnsbl[GeneAss] = Ensemble
            except Exception:
                continue
    return GeneEnsbl


def getGPR22(ec, time):
    try:
        NCBI_ID = "9606"
        page1 = "https://biocyc.org/META/NEW-IMAGE?type=EC-NUMBER&object=EC-" + ec
        page = str(urllib.request.urlopen(page1).read())
        page_cp = str(urllib.request.urlopen(page1).read())
        page_cp = str(urllib.request.urlopen(page1).read())
        page = str(page)
        url = page
        ss = []
        url = ""  # delete
        if re.search("class.+?Homo sapiens", url):
            if re.search(
                "<b>Species:<\/b> <i>Homo sapiens<\/i><br>.+?<b>Genes*:<\/b> *([A-Za-z/0-9-]+ *[A-Za-z/0-9-]*)",
                url,
            ):
                urls0 = re.search(
                    "<b>Species:<\/b> <i>Homo sapiens<\/i><br>.+?<b>Genes*:<\/b> *([A-Za-z/0-9-]+ *[A-Za-z/0-9-]*)",
                    url,
                ).group()
                if not NCBI_ID in urls0:
                    if re.findall(
                        "<b>Species:</b> <i>Homo sapiens</i><br> <b>Genes:</b> ([A-Za-z/0-9-]+, [A-Za-z/0-9-]+,* *[A-Za-z/0-9-]*,* *[A-Za-z/0-9-]*)<br>",
                        url,
                    ):
                        h = "h"
                    elif re.findall(
                        "<b>Species:<\/b> <i>Homo sapiens<\/i><br>.+?<b>Genes*:<\/b> *([A-Za-z/0-9-]+) *[A-Za-z/0-9-]*<br>",
                        url,
                    ):
                        urls0 = re.findall(
                            "<b>Species:<\/b> <i>Homo sapiens<\/i><br>.+?<b>Genes*:<\/b> *([A-Za-z/0-9-_]+) *[A-Za-z/0-9-_ ]*[A-Za-z/0-9-_ ]*<br>",
                            url,
                        )
                        for i in urls0:
                            if re.match("HS[0-9]{5}", i):
                                ss.extend([(i, i)])
                            if not re.search("HS[0-9]{5}", i) and not re.search(
                                "G[0-9]*-[0-9]{5}", i
                            ):
                                if re.search((i) + "*    +(HS[0-9]+)", url):
                                    r = re.search((i) + "*    +(HS[0-9]+)", url).group(
                                        0
                                    )
                                elif re.search((i) + "*    +(G[0-9]*-[0-9]+)", url):
                                    r = re.search(
                                        (i) + "*    +(G[0-9]*-[0-9]+)", url
                                    ).group(0)
                                elif re.search((i) + "*    +(HP_RS[0-9]*)", url):
                                    r = re.search(
                                        (i) + "*    +(HP_RS[0-9]*)", url
                                    ).group(0)
                                r = re.findall(
                                    "([A-Za-z/0-9-_]+) +([A-Za-z/0-9-_]+)", r
                                )
                                ss.extend(r)
        urls0 = ss

        if urls0:
            urls1 = [
                i for i in reversed(sorted([x[0:][0] for x in urls0], key=len))
            ]  # Sort genes by name lengh
            Ls = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
            urls2 = []
            urls22 = []
            for x in urls1:
                for y in urls0:
                    if re.findall(r"\('" + str(x) + "',", str(y)):
                        urls2 = urls2 + [y[1]]
                        urls22 = urls22 + [(y[0], y[1], Ls.pop(0))]
            a = page.replace(r"\n", r"").replace(
                r"Enzymes and Genes:", r"\nEnzymes and Genes:"
            )
            a2 = re.findall(r"Enzymes and Genes:.*[\S\s]+", a)
            b = (
                a2[0]
                .replace("<br> <a href=", "\n<br> <a href=")
                .replace("</a>)", "</a>)\n")
            )
            c = re.findall(
                r"<br> <a href=.*[\S\s].*>Homo sapiens</a>", b
            )  # evaluate the enzymes active in human
            if not c:
                c = re.findall(
                    r"<br> <a href=.*[\S\s].*>Mus musculus</a>", b
                )  # evaluate the enzymes active in mouse
            dict = {
                "<SUB>": "*",
                "<@SUB>(": " and (",
                "<@SUB>[": " and [",
                "][": "] and [",
                ")(": " and ",
                "'": "",
                ",": "",
                ";": "",
                "<BR>": "",
            }
            gpr2 = []
            i = 0
            while i < len(c):
                d = ""
                c2 = c[i].replace(" ", "\n")
                c3 = c2.replace('"\n', '" \n')
                isourl = re.findall(
                    r"(http://biocyc.org/META/NEW-IMAGE\?type=ENZYME.*)\" ",
                    c3.replace("/META", "http://biocyc.org/META"),
                )  # if it is a complex
                if not isourl:
                    isourl = re.findall(
                        r"<br>\n<a\nhref=\"(http://biocyc.org/gene\?orgid.*)\" ",
                        c3.replace("/gene?orgid", "http://biocyc.org/gene?orgid"),
                    )  # if it is not a complex
                isopage = getHtml(isourl[0], time)
                isopage = str(isopage)
                d1 = isopage.replace("\n", " ").replace("</a>", "\n</a>")
                d2 = ""
                if d2:
                    continue
                else:
                    isourl2 = re.findall(
                        b"/gene-tab.*META&tab=SUMMARY", isopage.encode("utf-8")
                    )
                    if isourl2:
                        isourl2 = ["http://biocyc.org" + isourl2[0].decode("utf-8")]
                        isopage2 = getHtml(isourl2[0], time)
                        isopage2 = str(isopage2)
                        d2 = re.findall(
                            r"Subunit Composition.+?\[(.+?)</",
                            re.sub(
                                "]*<SUB>",
                                "*",
                                isopage2.replace("<tr><td", "\n").replace(
                                    "</td></tr>", "\n"
                                ),
                            ),
                        )
                if not d2:
                    d1 = page.replace(r"\n", r" ").replace(r"</a>", r"</a>\n")
                    d2 = re.findall(
                        '\/META\/NEW-IMAGE\?type=REACTION&object=(.+?)"',
                        d1.replace(" <br>", ""),
                    )
                    d2 = d2[0]
                synonim = (
                    str(
                        re.findall(
                            "<p class=ecoparagraph>  Synonyms:  (.+?) <\/p>",
                            str(d1).replace(" #</p>", "\n"),
                        )
                    )
                    .replace("'", "")
                    .replace("]", "")
                    .replace("[", "")
                    .replace("  ", "")
                )
                if synonim:
                    synonim = synonim.split(",")
                    synonim = [
                        _f
                        for _f in [
                            re.sub(" $", "", re.sub("^ ", "", x)) for x in synonim
                        ]
                        if _f
                    ]
                if not d:  # defines d
                    d = [
                        _f
                        for _f in re.findall(
                            "[\S]+",
                            str(d2)
                            .replace("</SUB>", "<@SUB>")
                            .replace("/", " ")
                            .replace("(", "[")
                            .replace(")", "]"),
                        )
                        if _f
                    ]
                j = 0
                e = [""] * len(d)
                while j < len(d):
                    z = 0
                    while z < len(urls1):
                        # defines IsVar
                        IsVar = ""
                        if not re.findall("/NEW-IMAGE\?TYPE=REACTION", d[j].upper()):
                            if re.findall(urls1[z].upper(), d[j].upper()):
                                IsVar = [d[j].upper()]
                            if not IsVar and synonim:
                                dict3 = {}
                                for x in synonim:
                                    dict3[x.upper()] = urls1[z].upper()
                                IsVar = [
                                    _f
                                    for _f in [
                                        re.findall(
                                            x.upper()
                                            .replace("\\", "\\\\")
                                            .replace("(", "\(")
                                            .replace(")", "\)")
                                            .replace('"', '"')
                                            .replace("'", "'"),
                                            d[j].upper(),
                                        )
                                        for x in synonim
                                    ]
                                    if _f
                                ]
                                if IsVar:
                                    IsVar = [
                                        multiple_replace(dict3, x) for x in IsVar[0]
                                    ]
                                del dict3
                            if not IsVar:
                                IsVar = re.findall(urls1[z].upper(), d[j].upper())
                        if IsVar:
                            h = (
                                multiple_replace(dict, IsVar[0])
                                .replace("<@SUB>", "")
                                .replace("] and [", "]*1 and [")
                                .replace("] or [", "]*1 or [")
                                .replace(")", "")
                                .replace("(", "")
                            )
                            h = re.sub("]$", "]*1", h)
                            if not re.findall(r"]", h) or not re.findall(r"\[", h):
                                h = (
                                    "[" + str(h) + "]*1"
                                )  # to adapt the case of single gene reaction association
                                h = "[" + str(h) + "]*1"
                                isourl2 = re.findall(
                                    b"/gene-tab.*META&tab=SUMMARY",
                                    isopage.encode("utf-8"),
                                )
                                if isourl2:
                                    isourl2 = [
                                        "http://biocyc.org" + isourl2[0].decode("utf-8")
                                    ]
                                    isopage2 = getHtml(isourl2[0], time)
                                    isopage2 = str(isopage2)
                                    if re.findall(
                                        r"Subunit Composition.+?\[(.+?)</",
                                        re.sub(
                                            "]<SUB>",
                                            "*",
                                            isopage2.replace("<tr><td", "\n").replace(
                                                "</td></tr>", "\n"
                                            ),
                                        ),
                                    ):
                                        z = len(
                                            urls1
                                        )  # if we are analyzing a complex, the gpr is found, stop the iteration
                                        s = h
                                    if not re.findall(
                                        r"Subunit Composition.+?\[(.+?)</",
                                        re.sub(
                                            "]<SUB>",
                                            "*",
                                            isopage2.replace("<tr><td", "\n").replace(
                                                "</td></tr>", "\n"
                                            ),
                                        ),
                                    ):
                                        s = ""
                                        z += 1
                                else:
                                    s = ""
                                    z += 1
                            else:
                                NestLevel = re.findall(r"(?=(\]\*[0-9]+\]\*))", h)
                                h0 = ParseNestedParen(str(h), len(NestLevel))[0]
                                Mult = re.findall(r"([\]\*0-9]+$)", h)
                                if Mult:
                                    Mult = eval(
                                        str(Mult[0].split("*")[1:])
                                        .replace("'", "")
                                        .replace("]", "")
                                        .replace("[", "")
                                        .replace(", ", "*")
                                    )
                                else:
                                    Mult = 1
                                h = h0
                                for x in urls1:
                                    h = h.replace(x, "@")
                                h = [_f for _f in h.split("@") if _f]
                                dict4 = {
                                    "\\": "",
                                    '"': "",
                                    "'": "",
                                    "?": "",
                                    ")": "",
                                    "(": "",
                                    "]": "",
                                    "[": "",
                                    "*": "",
                                    "$": "",
                                    ".": "",
                                }
                                hh = []
                                for x in h:
                                    for w in urls1:
                                        for y in urls1:
                                            x1 = multiple_replace(dict4, x)
                                            m1 = str(y) + x1 + str(w)
                                            m2 = str(y) + x1
                                            if re.findall(
                                                m1, multiple_replace(dict4, h0)
                                            ) or re.findall(
                                                m2 + "$", multiple_replace(dict4, h0)
                                            ):
                                                t = re.findall(r"([0-9]+)", str(x))
                                                if t:
                                                    t = str(
                                                        eval(
                                                            str(t)
                                                            .replace("'", "")
                                                            .replace('"', "")
                                                            .replace("]", "")
                                                            .replace("[", "")
                                                            .replace(", ", "*")
                                                        )
                                                    )
                                                else:
                                                    t = "1"
                                                hh = hh + [str(y) + "*" + t]
                                hh = sorted(set([_f for _f in hh if _f]))
                                if not hh:
                                    hh = [str(h0) + "*1"]
                                if not hh:
                                    hh = [
                                        "[" + str(multiple_replace(dict4, h0)) + "]*1"
                                    ]  # no hh means that the string only contains the gene name
                                h22 = str(hh).replace("'", "").replace(", ", " and ")
                                UniqGene = sorted(
                                    set([x.split("*")[0] for x in hh])
                                )  # if some of the subunits of the complex encoded by the same gene it is necessary to change"S and S" by 2*A
                                GeneCount = [
                                    (x.split("*")[0], x.split("*")[1])
                                    for x in h22.replace("[", "")
                                    .replace("]", "")
                                    .split(" and ")
                                ]
                                UniqCount = []
                                for x in UniqGene:
                                    count = "0"
                                    for y in GeneCount:
                                        x = x.replace("[", "").replace("]", "")  #######
                                        if re.findall(x, y[0]):
                                            count = count + "+" + y[1]
                                    UniqCount = UniqCount + [[x, str(eval(count))]]
                                s = (
                                    "["
                                    + str(
                                        [str(x).replace("', '", "*") for x in UniqCount]
                                    )
                                    .replace("'", "")
                                    .replace('"', "")
                                    .replace(", ", " and ")
                                    + "]*"
                                    + str(Mult)
                                )
                                z = len(urls1)
                        else:
                            s = ""
                            z += 1
                    e[j] = s
                    j += 1
                e = sorted(set([_f for _f in e if _f]))
                if e:
                    gpr = (
                        str(e)
                        .replace("'[", "")
                        .replace("]'", "")
                        .replace("'", "")
                        .replace(", ", " and ")
                    )
                else:
                    gpr = str(urls0[i][0]) + "*1"
                gpr2 = gpr2 + [
                    str(gpr).replace("'[", "").replace("]'", "").replace("'", "")
                ]
                i += 1
            urls3 = [
                str(sorted(set(gpr2)))
                .replace("'[", "")
                .replace("]'", "")
                .replace("'", "")
                .replace(", ", " or ")
            ]
            testIs = [
                _f for _f in [x.replace("[", "").replace("]", "") for x in urls3] if _f
            ]
            if not urls3 or not testIs:
                urls3 = [str(x) + "*1" for x in urls1]
            g = []
            r = 0
            while r < len(urls3):
                GsprGenes = sorted(
                    set(
                        re.findall(
                            r"[A-Za-z0-9/-]+",
                            re.sub(
                                r"\*[0-9]+",
                                "",
                                str(urls3[r])
                                .replace(" and ", " ")
                                .replace(" or ", " "),
                            ),
                        )
                    )
                )
                GsprGenes = [
                    x.upper().replace("[", "").replace("]", "") for x in GsprGenes
                ]
                intersection = int(
                    len(GsprGenes)
                    - len(set(GsprGenes).intersection([x.upper() for x in urls1]))
                )
                if intersection < 1:
                    g = g + [
                        str(urls3[r].upper())
                        .replace(" AND ", " and ")
                        .replace(" OR ", " or ")
                    ]
                r += 1
            urls3 = (
                str(g)
                .replace("'', ", "")
                .replace("', ''", ")")
                .replace("', '", ") or (")
                .replace("'", "(")
                .replace("](]", "])]")
                .replace('"', "")
                .replace(" *", "*")
                .replace("[ ", "[")
                .replace(" ]", "]")
                .replace("(]", ")]")
                .replace("[)", "[(")
            )
            w = []
            for x in str(urls3).split(" or "):
                dict4 = {}
                for y in x.split(" and "):
                    y2 = re.sub("(\(|\)|\]|\[)", "", y)
                    a = re.findall(r"([A-Za-z0-9\-/]+)\*([0-9\*]+)", y2)
                    if not a[0][0] in dict4:
                        dict4[a[0][0]] = eval(a[0][1])
                    else:
                        dict4[a[0][0]] = eval(str(dict4[a[0][0]]) + "+" + a[0][1])
                z = [(x + "*" + str(dict4[x])) for x in dict4]
                w = w + [z]
            urls3 = (
                str(
                    sorted(
                        set([str(x).replace(", ", " and ").replace("'", "") for x in w])
                    )
                )
                .replace("'", "")
                .replace(", ", " or ")
                .replace("[", "([")
                .replace("]", "])")
            )
            urls4 = re.sub(r"\*[0-9]+", "", urls3)
        else:
            GPRURL22 = "http://www.genome.jp/dbget-bin/www_bget?ec:" + ec
            GPRPage2 = getHtml(GPRURL22, time)
            GPRPage2 = str(GPRPage2)
            urls0 = re.search("HSA.+?<table", GPRPage2)
            if urls0:
                urls0 = urls0.group()
                urls0 = re.findall("\((.+?)\)", urls0)
            if not urls0:
                urls0 = re.findall('mmu:[0-9]+">[0-9]+<\/a>\((.+?)\)', GPRPage2)
            if urls0:
                urls1 = [x[0:] for x in urls0]
                urls2 = urls1
                for i in range(len(urls2)):
                    oo = re.findall(urls2[i] + " +([HPRS_G]+[0-9]+)", page_cp)
                    if oo:
                        urls2[i] = oo[0]
                        # print('hel',urls2)
                urls3 = (
                    "[(["
                    + str(urls1)
                    .replace("[", "")
                    .replace("]", "")
                    .replace(", ", "*1]) or ([")
                    .replace("'", "")
                    + "*1])]"
                )
                urls4 = re.sub("\*[0-9]+", "", str(urls3))
            else:
                urls1 = ""
                urls2 = ""
                urls3 = ""
                urls4 = ""
        return urls0, urls1, urls2, urls3, urls4
    except Exception as e:
        raise
        return ""
    except Exception as error:
        print(error)


def getGPR_old(page, ec, time):
    try:
        NCBI_ID = "9606"
        if not page:
            page = "https://biocyc.org/META/NEW-IMAGE?type=EC-NUMBER&object=EC-" + ec
            page = str(urllib.request.urlopen(page).read())
        page = str(page)
        page_cp = page
        url = page
        ss = []
        if re.search("class.+?Homo sapiens", url):
            if re.search(
                "<b>Species:<\/b> <i>Homo sapiens<\/i><br>.+?<b>Genes*:<\/b> *([A-Za-z/0-9-]+ *[A-Za-z/0-9-]*)",
                url,
            ):
                urls0 = re.search(
                    "<b>Species:<\/b> <i>Homo sapiens<\/i><br>.+?<b>Genes*:<\/b> *([A-Za-z/0-9-]+ *[A-Za-z/0-9-]*)",
                    url,
                ).group()
                if not NCBI_ID in urls0:
                    if re.findall(
                        "<b>Species:</b> <i>Homo sapiens</i><br> <b>Genes:</b> ([A-Za-z/0-9-]+, [A-Za-z/0-9-]+,* *[A-Za-z/0-9-]*,* *[A-Za-z/0-9-]*)<br>",
                        url,
                    ):
                        h = "h"
                    elif re.findall(
                        "<b>Species:<\/b> <i>Homo sapiens<\/i><br>.+?<b>Genes*:<\/b> *([A-Za-z/0-9-]+) *[A-Za-z/0-9-]*<br>",
                        url,
                    ):
                        urls0 = re.findall(
                            "<b>Species:<\/b> <i>Homo sapiens<\/i><br>.+?<b>Genes*:<\/b> *([A-Za-z/0-9-_]+) *[A-Za-z/0-9-_ ]*[A-Za-z/0-9-_ ]*<br>",
                            url,
                        )
                        for i in urls0:
                            if re.match("HS[0-9]{5}", i):
                                ss.extend([(i, i)])
                            if not re.search("HS[0-9]{5}", i) and not re.search(
                                "G[0-9]*-[0-9]{5}", i
                            ):
                                if re.search((i) + "*    +(HS[0-9]+)", url):
                                    r = re.search((i) + "*    +(HS[0-9]+)", url).group(
                                        0
                                    )
                                elif re.search((i) + "*    +(G[0-9]*-[0-9]+)", url):
                                    r = re.search(
                                        (i) + "*    +(G[0-9]*-[0-9]+)", url
                                    ).group(0)
                                elif re.search((i) + "*    +(HP_RS[0-9]*)", url):
                                    r = re.search(
                                        (i) + "*    +(HP_RS[0-9]*)", url
                                    ).group(0)
                                r = re.findall(
                                    "([A-Za-z/0-9-_]+) +([A-Za-z/0-9-_]+)", r
                                )
                                ss.extend(r)
        urls0 = ss
        if urls0:
            urls1 = [
                i for i in reversed(sorted([x[0:][0] for x in urls0], key=len))
            ]  # Sort genes by name lengh
            Ls = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
            urls2 = []
            urls22 = []
            for x in urls1:
                for y in urls0:
                    if re.findall(r"\('" + str(x) + "',", str(y)):
                        urls2 = urls2 + [y[1]]
                        urls22 = urls22 + [(y[0], y[1], Ls.pop(0))]
            a = page.replace(r"\n", r"").replace(
                r"Enzymes and Genes:", r"\nEnzymes and Genes:"
            )
            a2 = re.findall(r"Enzymes and Genes:.*[\S\s]+", a)
            b = (
                a2[0]
                .replace("<br> <a href=", "\n<br> <a href=")
                .replace("</a>)", "</a>)\n")
            )
            c = re.findall(
                r"<br> <a href=.*[\S\s].*>Homo sapiens</a>", b
            )  # evaluate the enzymes active in human
            if not c:
                c = re.findall(
                    r"<br> <a href=.*[\S\s].*>Mus musculus</a>", b
                )  # evaluate the enzymes active in mouse
            dict = {
                "<SUB>": "*",
                "<@SUB>(": " and (",
                "<@SUB>[": " and [",
                "][": "] and [",
                ")(": " and ",
                "'": "",
                ",": "",
                ";": "",
                "<BR>": "",
            }
            gpr2 = []
            i = 0
            while i < len(c):
                d = ""
                c2 = c[i].replace(" ", "\n")
                c3 = c2.replace('"\n', '" \n')
                isourl = re.findall(
                    r"(http://biocyc.org/META/NEW-IMAGE\?type=ENZYME.*)\" ",
                    c3.replace("/META", "http://biocyc.org/META"),
                )  # if it is a complex
                if not isourl:
                    isourl = re.findall(
                        r"<br>\n<a\nhref=\"(http://biocyc.org/gene\?orgid.*)\" ",
                        c3.replace("/gene?orgid", "http://biocyc.org/gene?orgid"),
                    )  # if it is not a complex
                isopage = getHtml(isourl[0], time)
                isopage = str(isopage)
                d1 = isopage.replace("\n", " ").replace("</a>", "\n</a>")
                d2 = ""
                if d2:
                    continue
                else:
                    isourl2 = re.findall(
                        b"/gene-tab.*META&tab=SUMMARY", isopage.encode("utf-8")
                    )
                    if isourl2:
                        isourl2 = ["http://biocyc.org" + isourl2[0].decode("utf-8")]
                        isopage2 = getHtml(isourl2[0], time)
                        isopage2 = str(isopage2)
                        d2 = re.findall(
                            r"Subunit Composition.+?\[(.+?)</",
                            re.sub(
                                "]*<SUB>",
                                "*",
                                isopage2.replace("<tr><td", "\n").replace(
                                    "</td></tr>", "\n"
                                ),
                            ),
                        )
                if not d2:
                    d1 = page.replace(r"\n", r" ").replace(r"</a>", r"</a>\n")
                    d2 = re.findall(
                        '\/META\/NEW-IMAGE\?type=REACTION&object=(.+?)"',
                        d1.replace(" <br>", ""),
                    )
                    d2 = d2[0]
                synonim = (
                    str(
                        re.findall(
                            "<p class=ecoparagraph>  Synonyms:  (.+?) <\/p>",
                            str(d1).replace(" #</p>", "\n"),
                        )
                    )
                    .replace("'", "")
                    .replace("]", "")
                    .replace("[", "")
                    .replace("  ", "")
                )
                if synonim:
                    synonim = synonim.split(",")
                    synonim = [
                        _f
                        for _f in [
                            re.sub(" $", "", re.sub("^ ", "", x)) for x in synonim
                        ]
                        if _f
                    ]
                if not d:  # defines d
                    d = [
                        _f
                        for _f in re.findall(
                            "[\S]+",
                            str(d2)
                            .replace("</SUB>", "<@SUB>")
                            .replace("/", " ")
                            .replace("(", "[")
                            .replace(")", "]"),
                        )
                        if _f
                    ]
                j = 0
                e = [""] * len(d)
                while j < len(d):
                    z = 0
                    while z < len(urls1):
                        # defines IsVar
                        IsVar = ""
                        if not re.findall("/NEW-IMAGE\?TYPE=REACTION", d[j].upper()):
                            if re.findall(urls1[z].upper(), d[j].upper()):
                                IsVar = [d[j].upper()]
                            if not IsVar and synonim:
                                dict3 = {}
                                for x in synonim:
                                    dict3[x.upper()] = urls1[z].upper()
                                IsVar = [
                                    _f
                                    for _f in [
                                        re.findall(
                                            x.upper()
                                            .replace("\\", "\\\\")
                                            .replace("(", "\(")
                                            .replace(")", "\)")
                                            .replace('"', '"')
                                            .replace("'", "'"),
                                            d[j].upper(),
                                        )
                                        for x in synonim
                                    ]
                                    if _f
                                ]
                                if IsVar:
                                    IsVar = [
                                        multiple_replace(dict3, x) for x in IsVar[0]
                                    ]
                                del dict3
                            if not IsVar:
                                IsVar = re.findall(urls1[z].upper(), d[j].upper())
                        if IsVar:
                            h = (
                                multiple_replace(dict, IsVar[0])
                                .replace("<@SUB>", "")
                                .replace("] and [", "]*1 and [")
                                .replace("] or [", "]*1 or [")
                                .replace(")", "")
                                .replace("(", "")
                            )
                            h = re.sub("]$", "]*1", h)
                            if not re.findall(r"]", h) or not re.findall(r"\[", h):
                                h = (
                                    "[" + str(h) + "]*1"
                                )  # to adapt the case of single gene reaction association
                                h = "[" + str(h) + "]*1"
                                isourl2 = re.findall(
                                    b"/gene-tab.*META&tab=SUMMARY",
                                    isopage.encode("utf-8"),
                                )
                                if isourl2:
                                    isourl2 = [
                                        "http://biocyc.org" + isourl2[0].decode("utf-8")
                                    ]
                                    isopage2 = getHtml(isourl2[0], time)
                                    isopage2 = str(isopage2)
                                    if re.findall(
                                        r"Subunit Composition.+?\[(.+?)</",
                                        re.sub(
                                            "]<SUB>",
                                            "*",
                                            isopage2.replace("<tr><td", "\n").replace(
                                                "</td></tr>", "\n"
                                            ),
                                        ),
                                    ):
                                        z = len(
                                            urls1
                                        )  # if we are analyzing a complex, the gpr is found, stop the iteration
                                        s = h
                                    if not re.findall(
                                        r"Subunit Composition.+?\[(.+?)</",
                                        re.sub(
                                            "]<SUB>",
                                            "*",
                                            isopage2.replace("<tr><td", "\n").replace(
                                                "</td></tr>", "\n"
                                            ),
                                        ),
                                    ):
                                        s = ""
                                        z += 1
                                else:
                                    s = ""
                                    z += 1
                            else:
                                NestLevel = re.findall(r"(?=(\]\*[0-9]+\]\*))", h)
                                h0 = ParseNestedParen(str(h), len(NestLevel))[0]
                                Mult = re.findall(r"([\]\*0-9]+$)", h)
                                if Mult:
                                    Mult = eval(
                                        str(Mult[0].split("*")[1:])
                                        .replace("'", "")
                                        .replace("]", "")
                                        .replace("[", "")
                                        .replace(", ", "*")
                                    )
                                else:
                                    Mult = 1
                                h = h0
                                for x in urls1:
                                    h = h.replace(x, "@")
                                h = [_f for _f in h.split("@") if _f]
                                dict4 = {
                                    "\\": "",
                                    '"': "",
                                    "'": "",
                                    "?": "",
                                    ")": "",
                                    "(": "",
                                    "]": "",
                                    "[": "",
                                    "*": "",
                                    "$": "",
                                    ".": "",
                                }
                                hh = []
                                for x in h:
                                    for w in urls1:
                                        for y in urls1:
                                            x1 = multiple_replace(dict4, x)
                                            m1 = str(y) + x1 + str(w)
                                            m2 = str(y) + x1
                                            if re.findall(
                                                m1, multiple_replace(dict4, h0)
                                            ) or re.findall(
                                                m2 + "$", multiple_replace(dict4, h0)
                                            ):
                                                t = re.findall(r"([0-9]+)", str(x))
                                                if t:
                                                    t = str(
                                                        eval(
                                                            str(t)
                                                            .replace("'", "")
                                                            .replace('"', "")
                                                            .replace("]", "")
                                                            .replace("[", "")
                                                            .replace(", ", "*")
                                                        )
                                                    )
                                                else:
                                                    t = "1"
                                                hh = hh + [str(y) + "*" + t]
                                hh = sorted(set([_f for _f in hh if _f]))
                                if not hh:
                                    hh = [str(h0) + "*1"]
                                if not hh:
                                    hh = [
                                        "[" + str(multiple_replace(dict4, h0)) + "]*1"
                                    ]  # no hh means that the string only contains the gene name
                                h22 = str(hh).replace("'", "").replace(", ", " and ")
                                UniqGene = sorted(
                                    set([x.split("*")[0] for x in hh])
                                )  # if some of the subunits of the complex encoded by the same gene it is necessary to change"S and S" by 2*A
                                GeneCount = [
                                    (x.split("*")[0], x.split("*")[1])
                                    for x in h22.replace("[", "")
                                    .replace("]", "")
                                    .split(" and ")
                                ]
                                UniqCount = []
                                for x in UniqGene:
                                    count = "0"
                                    for y in GeneCount:
                                        x = x.replace("[", "").replace("]", "")  #######
                                        if re.findall(x, y[0]):
                                            count = count + "+" + y[1]
                                    UniqCount = UniqCount + [[x, str(eval(count))]]
                                s = (
                                    "["
                                    + str(
                                        [str(x).replace("', '", "*") for x in UniqCount]
                                    )
                                    .replace("'", "")
                                    .replace('"', "")
                                    .replace(", ", " and ")
                                    + "]*"
                                    + str(Mult)
                                )
                                z = len(urls1)
                        else:
                            s = ""
                            z += 1
                    e[j] = s
                    j += 1
                e = sorted(set([_f for _f in e if _f]))
                if e:
                    gpr = (
                        str(e)
                        .replace("'[", "")
                        .replace("]'", "")
                        .replace("'", "")
                        .replace(", ", " and ")
                    )
                else:
                    gpr = str(urls0[i][0]) + "*1"
                gpr2 = gpr2 + [
                    str(gpr).replace("'[", "").replace("]'", "").replace("'", "")
                ]
                i += 1
            urls3 = [
                str(sorted(set(gpr2)))
                .replace("'[", "")
                .replace("]'", "")
                .replace("'", "")
                .replace(", ", " or ")
            ]
            testIs = [
                _f for _f in [x.replace("[", "").replace("]", "") for x in urls3] if _f
            ]
            if not urls3 or not testIs:
                urls3 = [str(x) + "*1" for x in urls1]
            g = []
            r = 0
            while r < len(urls3):
                GsprGenes = sorted(
                    set(
                        re.findall(
                            r"[A-Za-z0-9/-]+",
                            re.sub(
                                r"\*[0-9]+",
                                "",
                                str(urls3[r])
                                .replace(" and ", " ")
                                .replace(" or ", " "),
                            ),
                        )
                    )
                )
                GsprGenes = [
                    x.upper().replace("[", "").replace("]", "") for x in GsprGenes
                ]
                intersection = int(
                    len(GsprGenes)
                    - len(set(GsprGenes).intersection([x.upper() for x in urls1]))
                )
                if intersection < 1:
                    g = g + [
                        str(urls3[r].upper())
                        .replace(" AND ", " and ")
                        .replace(" OR ", " or ")
                    ]
                r += 1
            urls3 = (
                str(g)
                .replace("'', ", "")
                .replace("', ''", ")")
                .replace("', '", ") or (")
                .replace("'", "(")
                .replace("](]", "])]")
                .replace('"', "")
                .replace(" *", "*")
                .replace("[ ", "[")
                .replace(" ]", "]")
                .replace("(]", ")]")
                .replace("[)", "[(")
            )
            w = []
            for x in str(urls3).split(" or "):
                dict4 = {}
                for y in x.split(" and "):
                    y2 = re.sub("(\(|\)|\]|\[)", "", y)
                    a = re.findall(r"([A-Za-z0-9\-/]+)\*([0-9\*]+)", y2)
                    if not a[0][0] in dict4:
                        dict4[a[0][0]] = eval(a[0][1])
                    else:
                        dict4[a[0][0]] = eval(str(dict4[a[0][0]]) + "+" + a[0][1])
                z = [(x + "*" + str(dict4[x])) for x in dict4]
                w = w + [z]
            urls3 = (
                str(
                    sorted(
                        set([str(x).replace(", ", " and ").replace("'", "") for x in w])
                    )
                )
                .replace("'", "")
                .replace(", ", " or ")
                .replace("[", "([")
                .replace("]", "])")
            )
            urls4 = re.sub(r"\*[0-9]+", "", urls3)
        else:
            GPRURL22 = "http://www.genome.jp/dbget-bin/www_bget?ec:" + ec
            GPRPage2 = getHtml(GPRURL22, time)
            GPRPage2 = str(GPRPage2)
            # print(GPRURL22)
            # urls0 = sorted(set(str(re.findall(r'(\([A-Za-z0-9]+\))', str(re.findall(r'hsa:............................................',GPRPage2.decode('utf-8'))))).replace("(","").replace(")","").replace("'","").replace(" ","").replace("[","").replace("]","").split(",")))
            # if not urls0[0]:
            # 	urls0 = sorted(set(str(re.findall('(\([A-Za-z0-9]+\))', str(re.findall('mmu:............................................',GPRPage2.decode('utf-8'))))).replace("(","").replace(")","").replace("'","").replace(" ","").replace("[","").replace("]","").split(",")))
            # urls0 = re.findall('hsa:[0-9]+">[0-9]+<\/a>\((.+?)\)',GPRPage2)
            urls0 = re.search("HSA.+?<table", GPRPage2)
            if urls0:
                urls0 = urls0.group()
                urls0 = re.findall("\((.+?)\)", urls0)
            if not urls0:
                urls0 = re.findall('mmu:[0-9]+">[0-9]+<\/a>\((.+?)\)', GPRPage2)
            if urls0:
                urls1 = [x[0:] for x in urls0]
                urls2 = urls1
                for i in range(len(urls2)):
                    oo = re.findall(urls2[i] + " +([HPRS_G]+[0-9]+)", page_cp)
                    if oo:
                        urls2[i] = oo[0]
                        # print('hel',urls2)
                urls3 = (
                    "[(["
                    + str(urls1)
                    .replace("[", "")
                    .replace("]", "")
                    .replace(", ", "*1]) or ([")
                    .replace("'", "")
                    + "*1])]"
                )
                urls4 = re.sub("\*[0-9]+", "", str(urls3))
                # print(urls4)
            else:
                urls1 = ""
                urls2 = ""
                urls3 = ""
                urls4 = ""
        return urls0, urls1, urls2, urls3, urls4
    except Exception as e:
        raise
        return ""
    except Exception as error:
        print(error)


####################### Warm up #######################
""""Download a HTML code"""


def getHtml(url, timeout, referer=False, file_data=[], additional_data={}):
    try:
        if additional_data:
            url += "?" + urllib.parse.urlencode(additional_data)
        if file_data:
            with open(file_data[1], "rb") as f:
                response = requests.post(url, files={file_data[0]: f})
                return response.tgext
        else:
            req = urllib.request.Request(url)
        if referer:
            req.add_header("Referer", referer)
        return urllib.request.urlopen(req).read()
    except Exception as e:
        time.sleep(timeout)
        return ""
    return ""


""""Download a HTMLS code"""


def getHtmlS(url, timeout):
    try:
        hdr = {
            "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.11 (KHTML, like Gecko) Chrome/23.0.1271.64 Safari/537.11",
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
            #       'Accept-Charset': 'ISO-8859-1,utf-8;q=0.7,*;q=0.3',
            "Accept-Encoding": "none",
            "Accept-Language": "en-US,en;q=0.8",
            "Connection": "keep-alive",
        }
        req = urllib.request.Request(url, headers=hdr)
        return urllib.request.urlopen(req).read()
    except Exception as e:
        time.sleep(timeout)
        print('Exception "' + str(e) + '" in getHtml with URL "' + url + '"')
        return ""
    return ""


""""Multiple Replacement"""


def multiple_replace(dict, text):
    # Create a regular expression  from the dictionary keys
    regex = re.compile("(%s)" % "|".join(map(re.escape, list(dict.keys()))))
    # For each match, look-up corresponding value in dictionary
    return regex.sub(lambda mo: dict[mo.string[mo.start() : mo.end()]], text)


"""Generate strings contained in nested (), indexing i = level"""


def ParseNestedParen(string, level):
    if len(re.findall("\[", string)) == len(re.findall("\]", string)):
        LeftRightIndex = [
            x
            for x in zip(
                [Left.start() + 1 for Left in re.finditer("\[", string)],
                reversed([Right.start() for Right in re.finditer("\]", string)]),
            )
        ]
    elif len(re.findall("\[", string)) > len(re.findall("\]", string)):
        return ParseNestedParen(string + "]", level)
    elif len(re.findall("\[", string)) < len(re.findall("\]", string)):
        return ParseNestedParen("[" + string, level)
    else:
        return "fail"
    return [string[LeftRightIndex[level][0] : LeftRightIndex[level][1]]]


####################### Path Ident. #######################
""""Path: Extract the links from a HTML page"""


def getLinkPath(page):
    try:
        urls0 = list(set(re.findall(r'name="rn:(R[0-9]+)" type="([a-z]+)', page)))
        urls1 = [
            x[0].replace("R", "http://www.kegg.jp/dbget-bin/www_bget?rn:R")
            for x in urls0
        ]
        urls21 = list(
            set(
                re.findall(
                    r'name="cpd:(C[0-9]+)" type="compound"\n.*link="([a-zA-Z0-9\:\.\+\-\_\?\/]+)',
                    page,
                )
            )
        )
        urls22 = list(
            set(
                re.findall(
                    r'name="gl:(G[0-9]+)" type="compound"\n.*link="([a-zA-Z0-9\:\.\+\-\_\?\/]+)',
                    page,
                )
            )
        )
        urls2 = urls21 + urls22
        urls3 = list(zip(urls0, urls1))
        return urls2, urls3
    except Exception as e:
        print(e)
        return ""


####################### Reaction Ident. #######################
""""Reaction: Extract the links from a HTML page"""


def getReacParam(page, time):
    try:
        # page = urllib.request.urlopen(page).read()
        urls0 = page.replace(b'href="', b"http://www.genome.jp")
        # urls0=str(urls0)
        # print(page)
        ID = re.findall(r"KEGG REACTION: ([A-Z0-9]+)", urls0.decode("utf-8"))
        # print(ID)
        # Associated compounds and linksg
        # 		urls11 = re.findall(r'(http://www.genome.jp[^\'" >]+).>(C\w+)<', urls0)  #Associated compounds and links
        # 		urls12 = re.findall(r'(http://www.genome.jp[^\'" >]+).>(G\w+)<', urls0)  #Associated compounds and links
        # 		urls1 = urls11+urls12
        a0 = (
            urls0.replace(b"\n", b"")
            .replace(b"<tr><th", b"\n<tr><th")
            .replace(b"</td></tr>", b"</td></tr>\n")
        )
        b0 = re.findall(r"<nobr>Equation</nobr>.*", a0.decode("utf-8"))
        b0 = re.findall(r"Equation.*", a0.decode("utf-8"))
        # b0 = re.findall(r'Equation.+?reaction', a0.decode('utf-8'),re.DOTALL)
        urls11 = [
            ("http://www.genome.jp/dbget-bin/www_bget?cpd:" + str(x), str(x))
            for x in sorted(set(re.findall(r"C[0-9]+", b0[0])))
        ]
        # 		urls11 = [('http://www.genome.jp/dbget-bin/www_bget?cpd:'+str(x),str(x)) for x in sorted(set(re.findall(r'G[0-9]+',page)))]
        urls12 = [
            ("http://www.genome.jp/dbget-bin/www_bget?gl:" + str(x), str(x))
            for x in sorted(set(re.findall(r"G[0-9]+", b0[0])))
        ]
        urls1 = urls11 + urls12
        # print(urls1)
        # Associated reagents and links
        # 		a1 = re.findall(r'(hidden">.*<a http://www.genome.jp[^\'" >]+.*</a><br>)', urls0)[0]
        # 		b1 = re.findall(r'(.*)&lt',a1)[0]
        # 		c1 = b1.replace(" ","").replace("hidden\">","</a>").replace("</a><","</a>1<").replace(">+<",">+1<").replace(">+",">").replace(">n<",">1<")
        # 		urls21 = re.findall(r'</a>([0-9]+)<a(http://www.genome.jp[^\'" >]+)">(C[0-9]+)',c1)
        # 		urls22 = re.findall(r'</a>([0-9]+)<a(http://www.genome.jp[^\'" >]+)">(G[0-9]+)',c1)
        a1 = re.findall(r"(.*)&lt", b0[0])[0]
        b1 = (
            a1.replace(" ", "")
            .replace('hidden">', "</a>")
            .replace("</a><", "</a>1<")
            .replace(">+<", ">+1<")
            .replace(">+", ">")
            .replace(">n<", ">1<")
        )
        c11 = (
            re.sub("<ahttp://www.genome\.jp/dbget-bin/www_bget\?cpd:C[0-9]+", "", b1)
            .replace(">C", ">1C")
            .replace(">", "\n")
        )
        d11 = ""
        if c11:
            d11 = re.findall(r"([0-9]+)(C[0-9]+)", c11)
        urls21 = [
            (
                str(x[0]),
                "http://www.genome.jp/dbget-bin/www_bget?cpd:" + str(x[1]),
                str(x[1]),
            )
            for x in d11
        ]
        c12 = (
            re.sub("<ahttp://www.genome\.jp/dbget-bin/www_bget\?gl:G[0-9]+", "", b1)
            .replace(">G", ">1G")
            .replace(">", "\n")
        )
        d12 = ""
        if c12:
            d12 = re.findall(r"([0-9]+)(G[0-9]+)", c12)
        urls22 = [
            (
                str(x[0]),
                "http://www.genome.jp/dbget-bin/www_bget?gl:" + str(x[1]),
                str(x[1]),
            )
            for x in d12
        ]
        urls2 = urls21 + urls22
        # Associated products and links
        # 		a2 = re.findall(r'(hidden">.*<a http://www.genome.jp[^\'" >]+.*</a><br>)', urls0)[0]
        # 		b2 = re.findall(r'&lt;(.*)',a2)[0]
        # 		c2 = b2.replace(" ","").replace("=>","</a>").replace("</a><","</a>1<").replace(">+<",">+1<").replace(">+",">").replace(">n<",">1<")
        # 		urls31 = re.findall(r'</a>([0-9]+)<a(http://www.genome.jp[^\'" >]+)">(C[0-9]+)',c2)
        # 		urls32 = re.findall(r'</a>([0-9]+)<a(http://www.genome.jp[^\'" >]+)">(G[0-9]+)',c2)
        a2 = re.findall(r"&lt;(.*)", b0[0])[0]
        b2 = (
            a2.replace(" ", "")
            .replace('hidden">', "</a>")
            .replace("</a><", "</a>1<")
            .replace(">+<", ">+1<")
            .replace(">+", ">")
            .replace(">n<", ">1<")
        )
        c21 = (
            re.sub("<ahttp://www.genome\.jp/dbget-bin/www_bget\?cpd:C[0-9]+", "", b2)
            .replace(">C", ">1C")
            .replace(">", "\n")
        )
        d21 = ""
        if c21:
            d21 = re.findall(r"([0-9]+)(C[0-9]+)", c21)
        urls31 = [
            (
                str(x[0]),
                "http://www.genome.jp/dbget-bin/www_bget?cpd:" + str(x[1]),
                str(x[1]),
            )
            for x in d21
        ]
        c22 = (
            re.sub("<ahttp://www.genome\.jp/dbget-bin/www_bget\?gl:G[0-9]+", "", b2)
            .replace(">G", ">1G")
            .replace(">", "\n")
        )
        d22 = ""
        if c22:
            d22 = re.findall(r"([0-9]+)(G[0-9]+)", c22)
        urls32 = [
            (
                str(x[0]),
                "http://www.genome.jp/dbget-bin/www_bget?gl:" + str(x[1]),
                str(x[1]),
            )
            for x in d22
        ]
        urls3 = urls31 + urls32
        # Name
        a3 = re.findall(
            '<nobr>Name.*\n.*hidden">([^\n]+)<br>\n</div></div></td></tr>',
            urls0.decode("utf-8"),
        )
        if a3:
            b3 = a3[0].replace("<br>\n", "").replace(" ", "_")
            # c3 = re.findall('[A-Za-z0-9\_\-\+:\/\\\]+', b3)
            # urls4 = [x.replace("_"," ").replace(";","") for x in c3]
            urls4 = [b3]
        else:
            a31 = re.findall(
                '(<nobr>Definition[\S\s]+:hidden">[A-Za-z0-9+=;,-<>& ()\[\]]+)<br>',
                urls0.decode("utf-8"),
            )
            if a31:
                a3 = re.findall(
                    r":hidden\">.*:hidden\">(.*)", a31[0].split("<br>\n</div>")[0]
                )
            # 			a3 = re.findall(r':hidden\">(.*)',re.findall('(<nobr>Definition[\S\s]+:hidden">[A-Za-z0-9+=;,-<>& ()\[\]]+)<br>', urls0)[0].split("<br>\n</div>")[0])
            if a3:
                # 				b = a[0].replace("<br>\n" , "").replace(" " , "_").replace(";" , " ")
                # 				c = re.findall('[A-Za-z0-9\_\-\+:\/\\\]+', b)
                # 				urls4 = [x.replace("_"," ") for x in c]
                urls4 = [a3[0].replace(" ", "_")]
            else:
                urls4 = ""
        # EC
        # 		urls5 = re.findall(r"Enzyme[\S\s]+?.*(http[a-zA-Z0-9./:?_-]+)\">([0-9.]+)", urls0)
        ##		if not urls5:
        ##			urls5 = re.findall(r"hidden\">([0-9.-]+)<br>", urls0)
        # 		if not urls5: # If there is not E.C. number, then check in all the databases
        # 			a4 = getHtml('http://www.genome.jp/dbget-bin/get_linkdb?-t+alldb+rn:'+ID,time)
        # 			b4 = re.findall(r"hsa:.*EC:([0-9.]+)[\)\]]", a4)
        # 			if not b4:
        # 				b4 = re.findall(r"hsa:.*EC:([0-9.]+)", a4)
        # 			if b4:
        # 				c4 = 'http://www.genome.jp/dbget-bin/www_bget?ec:'+b4[0]
        # 				urls5 = [[c4,b4[0]]]
        # 			if not urls5: # If there is not E.C. number, then check in the ortology
        # 				a4 = re.findall(r"Orthology[\S\s]+?.*(http[a-zA-Z0-9./:?_-]+)\"", urls0)
        # 				b4 = getHtml(a4[0],time)
        # 				c4 = re.findall(r"hsa:([0-9]+)", b4) # check the gene
        # 				d4 = getHtml("http://www.genome.jp/dbget-bin/www_bget?hsa:"+c4[0],time)
        # 				e4 = re.findall(r"EC:([0-9.]+)[\)\]]", d4)
        # 				f4 = 'http://www.genome.jp/dbget-bin/www_bget?ec:'+d4[0]
        # 				urls5 = [[f4,e4[0]]]
        # 				if not urls5: # If there is not E.C. number, then check in uniprot
        # 					a4 = re.findall(r"(http://www.uniprot.org/uniprot/[A-Za-z0-9]+)\"", d4)
        # 					b4 = getHtml(a4[0],time) #Uniprot
        # 					c4 = re.findall(r"EC:.*EC/([0-9.]+)\"", f4)
        # 					d4 = 'http://www.genome.jp/dbget-bin/www_bget?ec:'+c4[0]
        # 					urls5 = [[d4,c4[0]]]
        a = urls0.replace(b"\n", b"").replace(
            b"Enzymes and Genes:", b"\nEnzymes and Genes:"
        )
        b = a.replace(b"<nobr>Enzyme</nobr>", b"\n<nobr>Enzyme</nobr>").replace(
            b"</div></td></tr>", b"</div></td></tr>\n"
        )
        c = re.findall(r"Enzyme.+?entry.([0-9.]+)", b.decode("utf8"), re.DOTALL)
        alt = re.findall(r"hidden\">([0-9.-]+)<br>", urls0.decode("utf8"))
        urls5 = ""
        if c:
            urls5 = re.findall(r"(http[a-zA-Z0-9./:?_-]+)\">([0-9.]+)", c[0])
        # 		urls5 = re.findall(r"Enzyme[\S\s]+?.*(http[a-zA-Z0-9./:?_-]+)\">([0-9.]+)", urls0)
        if not urls5:  # If there is not E.C. number, then check in all the databases
            a4 = getHtml(
                "http://www.genome.jp/dbget-bin/get_linkdb?-t+alldb+rn:" + str(ID), time
            )
            # 			b4 = re.findall(r"hsa:.*EC:([0-9.]+)[\)\]]", a4)
            # Handle both str and bytes
            a4_str = a4.decode("utf-8") if isinstance(a4, bytes) else a4
            b4 = sorted(
                set(
                    str(
                        [
                            _f
                            for _f in [
                                re.findall(r"[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+", x[0])
                                for x in [
                                    _f
                                    for _f in [
                                        re.findall(r"EC:(.*)", x)
                                        for x in [
                                            x for x in re.findall(r"hsa:.*", a4_str)
                                        ]
                                    ]
                                    if _f
                                ]
                            ]
                            if _f
                        ]
                    )
                    .replace("[", "")
                    .replace("]", "")
                    .replace("'", "")
                    .split(", ")
                )
            )
            if not b4:
                b4 = re.findall(r"hsa:.*EC:([0-9.]+)", a4)
                if not alt:
                    alt = re.findall(r"hsa:.*EC:([0-9.-]+)", a4)
            if b4[0]:
                # 				c4 = 'http://www.genome.jp/dbget-bin/www_bget?ec:'+b4[0]
                # 				urls5 = [[c4,b4[0]]]
                urls5 = []
                for y in b4:
                    urls5 = urls5 + [
                        ["http://www.genome.jp/dbget-bin/www_bget?ec:" + str(y), str(y)]
                    ]
            del a4
            del b4
        if not urls5:  # If there is not E.C. number, then check in the ortology I
            a4 = re.findall(
                r"Orthology.*",
                urls0.decode("utf-8")
                .replace("\n", "")
                .replace("<nobr>Orthology</nobr>", "\n<nobr>Orthology</nobr>")
                .replace("</table></td></tr>", "</table></td></tr>\n"),
            )
            if a4:
                b4 = re.findall(">([0-9\.]+)</a>", a4[0])
                if b4:
                    urls5 = []
                    for y in b4:
                        urls5 = urls5 + [
                            [
                                "http://www.genome.jp/dbget-bin/www_bget?ec:" + str(y),
                                str(y),
                            ]
                        ]
        if not urls5:  # If there is not E.C. number, then check in the ortology II
            a4 = re.findall(
                r"Orthology[\S\s]+?.*(http[a-zA-Z0-9./:?_-]+)\"", urls0.decode("utf-8")
            )
            if a4:
                b4 = getHtml(a4[0], time)
                c4 = re.findall(b"hsa:([0-9]+)", b4)  # check the gene
                del b4
                if c4:
                    d4 = getHtml(
                        "http://www.genome.jp/dbget-bin/www_bget?hsa:" + str(c4[0]),
                        str(time),
                    )
                    e4 = re.findall(r"EC:([0-9.]+)[\)\]]", d4.decode("utf-8"))
                    if not alt:
                        alt = re.findall(r"EC:([0-9.-]+)[\)\]]", d4.decode("utf-8"))
                    if e4:
                        f4 = "http://www.genome.jp/dbget-bin/www_bget?ec:" + d4[0]
                        urls5 = [[f4, e4[0]]]
                        del f4
                    del e4
            if not urls5:
                d4 = ""
            del a4
        if not urls5:  # If there is not E.C. number, then check in uniprot
            if d4:
                a4 = re.findall(r"(http://www.uniprot.org/uniprot/[A-Za-z0-9]+)\"", d4)
                if a4:
                    b4 = getHtml(a4[0], time)  # Uniprot
                    c4 = re.findall(r"EC:.*EC/([0-9.]+)\"", b4)
                    if c4:
                        d4 = "http://www.genome.jp/dbget-bin/www_bget?ec:" + c4[0]
                        urls5 = [[d4, c4[0]]]
        if (
            not urls5
        ):  # If there is not E.C. number, then check in pair reaction in KEGG
            if re.findall(b"<nobr>RPair</nobr>", urls0):
                a4 = (
                    urls0.replace("\n", "")
                    .replace("Enzymes and Genes:", "\n<nobr>RPair</nobr>")
                    .replace("</table></td></tr>", "</table></td></tr>\n")
                )
                b4 = re.findall(r"<nobr>RPair</nobr>.*", a4)[0].replace(
                    "http", "\nhttp"
                )
                c4 = re.findall(r'(http.*)">RP', b4)
                if c4:
                    urls5 = []
                    for x in c4:
                        d4 = getHtml(x, time)
                        e4 = d4.replace("\n", "").replace(
                            "Enzymes and Genes:", "\nEnzymes and Genes:"
                        )
                        f4 = e4.replace(
                            "<nobr>Enzyme</nobr>", "\n<nobr>Enzyme</nobr>"
                        ).replace("</div></td></tr>", "</div></td></tr>\n")
                        g4 = re.findall(b"(Enzyme.*\<\/tr>)", f4)
                        if g4:
                            g4 = [
                                _f
                                for _f in [
                                    re.findall(b'(http.*)">(.*)', x)
                                    for x in g4[0]
                                    .replace('href="', "\nhttp://www.genome.jp")
                                    .split("</a>")
                                ]
                                if _f
                            ]
                            for y in g4:
                                if y[0][0] not in [x[0] for x in urls5]:
                                    urls5 = urls5 + [[y[0][0], y[0][1]]]
        if not urls5:
            if alt:
                urls5 = [
                    ["http://www.genome.jp/dbget-bin/www_bget?ec:" + alt[0], alt[0]]
                ]
            else:
                urls5 = [["", ""]]
        # Equivalent reaction in case the original reaction includes Glycans
        if urls32 or urls22:
            a5 = re.findall(
                'Remark[\S\s]+Same as:[\S\s]+">(R[0-9]+)</a>', urls0.decode("utf-8")
            )
            if a5:
                urls6 = a5[0]
            if not a5:
                urls6 = "1"
        else:
            urls6 = "0"
        # Are glycans involved in the reaction?
        if urls12:
            urls7 = 1
        else:
            urls7 = 0
        # Are metabolites involved in the reaction?
        if urls11:
            urls8 = 1
        else:
            urls8 = 0
        # print(urls5)
        # 		session = setup_biocyc_session()
        # 		print('ec', urls5[0][1])
        # 		urls9 = getGPR(urls5[0][1],session) # new
        # 		print('gpr', urls9)
        ##		urls9 = getGPR_old(str(),urls5[0][1],20) # old
        ##		urls9 = getLocation(urls9[3], urls9[4], urls9[2],20)[2]  # old
        # 		urls9 = getLocation(urls9[3], urls9[4], urls9[2],1, 'pkl/Human_variables.pkl',session)[2]  # old
        # 		print('location', urls9)
        print(urls1, urls2, urls3, urls4, urls5, urls6, urls7, urls8)
        return [urls1, urls2, urls3, urls4, urls5, urls6, urls7, urls8]
    # 		return [urls1 , urls2 , urls3 , urls4 , urls5 , urls6 , urls7 , urls8, urls9]
    except Exception as e:
        print(traceback.format_exc())
        # print('Exception "'+str(e)+'" in getReacParam') !!!!!
        return ""


""""Reaction: Evaluate the consistency between the species of the reaction"""


def getRxncons(rxn, time, MetEquiv, MetList, MetIdent, EF, specialCompounds):
    from functions.class_generate_database import compound, reaction

    RxnCmp = [x[2] for x in rxn.Substrate()] + [x[2] for x in rxn.Product()]
    # Cmptest = [x for x in RxnCmp if "C" in x]
    Glytest = [x for x in RxnCmp if "G" in x]
    if not Glytest:  # if there is only Compounds don't do any change
        reaction2 = rxn
    else:  # if there are Glycans
        # 1st check if there is an equivalent reaction with compounds
        if rxn.Equivalent() != "1":
            RxnID2 = rxn.Equivalent()
            RxnURL2 = "http://www.kegg.jp/dbget-bin/www_bget?rn:" + rxn.Equivalent()
            PathName2 = rxn.Pathway()
            RxnTermDyn2 = rxn.Termodyn()
            reaction2 = reaction(RxnURL2, time, RxnID2, PathName2, RxnTermDyn2)
            reaction2.Equivalent = lambda: str(rxn.ID)
            # use the most informative EC number set between the equivalent reactions
            test_ec_ref = 0  # 0 assumes that both EC numbers are the same
            for x in reaction2.EC():
                for y in rxn.EC():
                    test_ec = (len(re.findall("\.", x)) - len(re.findall("\-", x))) - (
                        len(re.findall("\.", y)) - len(re.findall("\-", y))
                    )  # test_ec has valuese between -3 and 3. 1.test_ec = 0: quality of EC is equal in both reactions, test_ec < 0 quality of EC is higher in x, test_ec > 0: quality of EC is higher in y
                    if test_ec < test_ec_ref:
                        test_ec_ref = test_ec
            if test_ec_ref == 0:
                test_ec_ref = len(reaction2.EC()) - len(rxn.EC())
            if test_ec_ref < 0:
                reaction2.EC = copy.deepcopy(rxn.EC)
            ######### Check if all the compounds in the jth reaction are in the compound list ###########
            RxnCmp2 = [x[2] for x in reaction2.Substrate()] + [
                x[2] for x in reaction2.Product()
            ]
            c = 0
            while c < len(RxnCmp2):
                CompID = RxnCmp2[c]
                if not CompID in MetIdent and not CompID in MetEquiv:
                    MetIdent = MetIdent + [CompID]
                    if CompID[0] == "C":
                        CURL = "http://www.kegg.jp/dbget-bin/www_bget?cpd:" + CompID
                    if CompID[0] == "G":
                        CURL = "http://www.kegg.jp/dbget-bin/www_bget?gl:" + CompID
                    MetList[CompID] = compound(CURL, CompID, time, EF, specialCompounds)
                    if MetList[CompID].ID1 != MetList[CompID].ID2:
                        MetIdent[len(MetIdent) - 1] = MetList[CompID].ID1
                        MetEquiv[CompID] = MetList[CompID].ID1
                        MetList[MetList[CompID].ID1] = copy.deepcopy(
                            MetList[CompID]
                        )  # Change the reference in the dictionary to account for the 1th ID
                        del MetList[CompID]
                c = c + 1
        # 2st check if the glycans have associated compounds
        else:
            # Substrates
            Sini = copy.deepcopy(rxn.Substrate())
            GlyS = [x[2] for x in rxn.Substrate()]
            s = 0
            count = 0
            while s < len(GlyS):
                if (
                    GlyS[s] in MetEquiv
                ):  # If the Glycan have an alternative compund, then, replace it in the reaction
                    NewS = MetEquiv[GlyS[s]]
                    count = count + 1
                    Sini[s][1] = Sini[s][1].replace("gl:" + GlyS[s], "cpd:" + NewS)
                    Sini[s][2] = Sini[s][2].replace(GlyS[s], NewS)
                s = s + 1
            # Products
            Pini = copy.deepcopy(rxn.Product())
            GlyP = [x[2] for x in rxn.Product()]
            p = 0
            while p < len(GlyP):
                if (
                    GlyP[p] in MetEquiv
                ):  # If the Glycan have an alternative compund, then, replace it in the reaction
                    NewP = MetEquiv[GlyP[p]]
                    count = count + 1
                    Pini[p][1] = Pini[p][1].replace("gl:" + GlyP[p], "cpd:" + NewP)
                    Pini[p][2] = Pini[p][2].replace(GlyP[p], NewP)
                p = p + 1
            # Analize the results of substrates and products
            if (
                len(GlyS) + len(GlyP)
            ) - count == 0:  # all the glycans have an associated compound
                ######### Check if all the compounds in the jth reaction are in the compound list ###########
                RxnCmp = [x[2] for x in Pini] + [x[2] for x in Sini]
                c = 0
                while c < len(RxnCmp):
                    CompID = RxnCmp[c]
                    if not CompID in MetIdent and not CompID in MetEquiv:
                        MetIdent = MetIdent + [CompID]
                        if CompID[0] == "C":
                            CURL = "http://www.kegg.jp/dbget-bin/www_bget?cpd:" + CompID
                        if CompID[0] == "G":
                            CURL = "http://www.kegg.jp/dbget-bin/www_bget?gl:" + CompID
                        MetList[CompID] = compound(
                            CURL, CompID, time, EF, specialCompounds
                        )
                        if MetList[CompID].ID1 != MetList[CompID].ID2:
                            MetIdent[len(MetIdent) - 1] = MetList[CompID].ID1
                            MetEquiv[CompID] = MetList[CompID].ID1
                            MetList[MetList[CompID].ID1] = copy.deepcopy(
                                MetList[CompID]
                            )  # Change the reference in the dictionary to account for the 1th ID
                            del MetList[CompID]
                    c = c + 1
                ######### Replace the old glycans by the new compounds ###########
                reaction2 = copy.deepcopy(rxn)
                reaction2.SetProduct(Pini)
                reaction2.SetSubstrate(Sini)
            else:  # not all the glycans have an associated compound
                if (
                    len(RxnCmp) - len(Glytest) == 0
                ):  # if the reaction only has glycans don't do anything because it will be transformed latelly to be mass balanced
                    reaction2 = rxn
                    reaction2.MBTest = (
                        lambda: "1"
                    )  # if all the component of the reaction are glycans then, it is necessary to transform their composition in order to be mass balanced
                else:  # assume that the stoichometry is 1 for all the elements of the reaction
                    reaction2 = rxn
                    reaction2.MBTest = (
                        lambda: "0"
                    )  # modify this parameter to avoid mass balance because it is assumed  a 1 to 1 stoichometry (default value = R -then the reaction is balanced)
    return reaction2


### Called directly by meltGeneList
### Collapse list of genelists
### For model builing
def meltGene(geneList3):
    geneList4 = dict()
    for i in [x for x in geneList3]:
        for CSL, gene in i.items():
            # Normalize gene to string so concatenation below cannot fail when gene is list/tuple
            if isinstance(gene, (list, tuple, set)):
                # join multiple gene names with ' or ' (choose or because these represent alternatives)
                gene_str = " or ".join([str(x) for x in gene])
            else:
                gene_str = str(gene)

            if not CSL in geneList4:
                geneList4[CSL] = gene_str
            elif CSL in geneList4 and gene_str != geneList4[CSL]:
                geneList4[CSL] += " and " + gene_str
    return geneList4


### Collapse list of genelists
### For reactions with more than one EC specified
### For model builing
def meltGeneList(listOfgenelists):
    geneList = list()
    n = 0
    for i in range(4):
        geneList2 = []
        for i2 in range(len(listOfgenelists)):
            geneList3 = listOfgenelists[i2][i]
            geneList2.append(geneList3)
        geneList.append(meltGene(geneList2))
    # geneList.append(sorted(list(set([listOfgenelists[x][4] for x in range(len(listOfgenelists))])))[-1])
    return geneList


def getLocation_old(gpr, genelist1, genelist2, time):
    # 	if genelist2: print(genelist2)
    # 	print(gpr,genelist1,genelist2,time)
    # print(genelist2)
    gprgpr = ""
    ppList = list()
    OtherLocations = ["Other locations"]
    try:
        # print(genelist1)
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
        # print(0, urls0)
        # 		dict = {"Cell" : "Cytosol", "Membrane" : "Cytosol", "Lipid-anchor" : "Cytosol", "Cytosol membrane" : "Cytosol", "Multi-pass membrane protein" : "Cytosol", "Single-pass membrane protein" : "Cytosol","terminal bouton":"Cytosol","catalytic complex":"Cytosol", "endocytic vesicle":"Cytosol", "presynapse":"Cytosol","neuron projection":"Cytosol","synapse":"Cytosol","glutamatergic synapse":"Cytosol","lamellipodium":"Cytosol","postsynapse":"Cytosol","transport vesicle":"Cytosol","AMPA glutamate receptor complex":"Cytosol","caveola":"Cytosol","excitatory synapse":"Cytosol", " Mitochondrion inner membrane" : "Mitochondrion", "Intermembrane side" : "Mitochondrion", "Endoplasmic reticulum membrane" : "Endoplasmic"} # expand this dict
        Path = "pkl/Human_variables.pkl"  # Endo1a_variables.pkl #Human_variables.pkl
        Var = open(Path, "rb")  # Endo1b_variables.pkl
        hh = 1  # dict and not all
        dictt = pickle.load(Var)
        if hh == 0:
            print(Path + " is used")
            # print(dictt.keys())
        LocationList = []
        Locations = []
        uniprot = ""
        biocyc = ""
        n = 0
        # if urls0: print(urls0)
        for n in range(len(urls0)):  # isoforms
            # print(n)
            # print(urls0)
            # print(urls0[n])
            a = urls0[n].split("and")
            m = 0
            d = ""
            # while m < len(a): # subunits
            for m in range(len(a)):
                # print str(n)+'_'+str(m)+'_00'
                p = 1
                b = (
                    "http://www.genome.jp/dbget-bin/www_bget?sp:"
                    + re.sub(r"\*[0-9]+", "", a[m])
                    + "_HUMAN"
                )  # location in genome net human
                # print(b)
                bb = str(getHtml(b, time))  # .decode('utf-8')
                # cc = bb.replace("\nCC","").replace("-!-","\n")
                # dd = re.findall(r"GO:[0-9]+.*;.*C:(.*);",cc)
                dd = re.findall("GO:[0-9]+.+?C:(.+?);", bb)
                # print str(n)+'_'+str(m)+'_A'
                # print(dd,1000)
                # if dd:
                # 	print(b)
                # 	print(dd)
                if not dd:
                    ddd = re.search("(SUBCELLULAR LOCATION:.*)", "")
                    if ddd:
                        dd = [ddd.group(1)]
                        # print(b)
                        # print(dd)
                        # print str(n)+'_'+str(m)+'_A2'
                if not dd or dd:
                    # print(77777)
                    # b = 'http://www.genome.jp/dbget-bin/www_bget?sp:'+re.sub(r"\*[0-9]+" , "", a[m])+'_MOUSE' #location in genome net mouse
                    # print(b)
                    # bb =  str(getHtml(b,time))
                    # cc = bb.replace("\nCC","").replace("-!-","\n")
                    # dd = re.findall(r"GO:[0-9]+.*;.*C:(.*);",cc) # GANC
                    # if dd:
                    # 	print(b)
                    # 	print('mus')
                    # print str(n)+'_'+str(m)+'_B'
                    # if not dd:
                    # 	ddd = re.search('(SUBCELLULAR LOCATION:.*)',cc)
                    # 	if ddd:
                    # 			dd = [ddd.group(1)]
                    # 		ddd = re.search('(SUBCELLULAR LOCATION:.*)',cc)
                    # 		print(b)
                    # print(dd)
                    # print(ddd)
                    # print str(n)+'_'+str(m)+'_B2'
                    # dd =[]
                    if genelist1:
                        # print(a[m].upper())
                        # if not type(genelist1) == '''<class 'str'>''': print(type(genelist1))
                        if isinstance(genelist1, str):
                            genelist11 = re.findall("\[(.+?)\]", genelist1)
                            genelist11 = [
                                gene.replace("(", "")
                                .replace(")", "")
                                .replace("[", "")
                                .replace("]", "")
                                for gene in genelist11
                            ]
                        if not isinstance(genelist1, str):
                            genelist11 = genelist1
                        index = [x.upper() for x in genelist11].index(
                            re.sub(r"\*[0-9]+", "", a[m].upper())
                        )
                        # print(index)
                        # print(genelist2)

                        # print(genelist11[index].upper())
                        b = (
                            "http://biocyc.org/gene?orgid=META&id="
                            + genelist2[index].upper()
                        )  # location in BioCyc human #########################
                        # print(b)
                        # print(genelist11[index].upper(), re.sub(r"\*[0-9]+" , "", a[m].upper()))
                        # if genelist11[index].upper() != re.sub(r"\*[0-9]+" , "", a[m].upper()):    # !!!????
                        # if b:
                        bb = str(getHtml(b, time))
                        if bb:
                            # print(88)
                            # print(d)
                            if re.search("Location", bb, flags=re.DOTALL):
                                # print(b)
                                ddd1 = re.findall("Locations?.+?Reactions?", bb)[0]
                                ddd = ddd = [
                                    re.sub(r"\\n", "", i)
                                    for i in re.findall(
                                        "([\\\\n+]?[a-z ]+[a-zA-Z() ]+) <a", ddd1
                                    )
                                ]
                                if not ddd:
                                    ddd = [
                                        re.sub(r"\\n", "", i)
                                        for i in re.findall(
                                            "(\\\\[a-zA-Z, ()]+)</", ddd1
                                        )
                                    ]
                                    # print(ddd)
                                    ddd = [
                                        re.sub("^ ", "", i) for i in ddd[0].split(",")
                                    ]
                                # print(ddd)
                                dd.extend(ddd)
                            # 							bb = bytes(bb)
                            # 							cc = cc.encode('utf-8')
                            # cc = str(bb).replace(r"\n",r"").replace(r"<td align=RIGHT valign=TOP class=\"label\">",r"\n<td align=RIGHT valign=TOP class=\"label\"")
                            # 							dd = re.search('<td align=RIGHT valign=TOP class=\"label\"[\S\s]+Location(.*)',cc)
                            # dd = re.findall(r'<td align=RIGHT valign=TOP class=\"label\"[\S\s]+Location(.*)',cc)
                            # if dd:
                            # print(b)
                            # 	#print(dd)
                            # dd=''
                            biocyc = "1"
                            # print str(n)+'_'+str(m)+'_D'
                            if bb:  # location in Uniprot
                                biocyc = ""
                                # cc = str(bb).replace(r"\n",r"").replace(r"http://www.uniprot.org",r"\nhttp://www.uniprot.org").replace(r"nbsp",r"nbsp\n")
                                # cc2 = re.findall(r'(http://www.uniprot.org.*)\">',cc)
                                # print(cc2)
                                # if cc2:
                                # 	cc3 =  getHtml(cc2[0],time)
                                # 	cc4 = cc3.replace("\n","").replace("</span>Subcellular location","\n</span>Subcellular location").replace("&#xd;","\n")
                                # 								#	dd = re.search('</span>Subcellular location(.*)',cc4)
                                # 	dd = re.findall('</span>Subcellular location(.*)',cc4)
                                # 		if dd:
                                # 			print(cc2[0])
                                # 			#print(dd)
                                if re.search("uniprot/([A-Z0-9]+)", bb):
                                    cc = re.search("uniprot/([A-Z0-9]+)", bb).group(1)
                                    b = (
                                        "https://rest.uniprot.org/uniprotkb/"
                                        + cc
                                        + ".txt"
                                    )
                                    # print(b)
                                    bb = str(getHtml(b, time))
                                    ddd = re.findall("GO:[0-9]+.+?C:(.+?);", bb)
                                    # print(ddd)
                                    ddd = [x for x in ddd if not "GO" in x]
                                    dd.extend(ddd)
                                    # if ddd:
                                    # 	print(b)
                                    # 	print(ddd)
                                    uniprot = "1"
                                    # print str(n)+'_'+str(m)+'_E'
                        if not dd:
                            UniProtKB = ""
                            GeneID = re.sub(r"\*[0-9-]+", "", a[m])
                            # print(GeneID)
                            # print(GeneID)
                            https = (
                                "https://www.uniprot.org/uniprot/?query="
                                + GeneID
                                + "&sort=score"
                            )
                            url = str(urllib.request.urlopen(https).read())
                            # print(https)
                            if re.search(
                                'uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>'
                                + GeneID
                                + "_HUMAN",
                                url,
                            ):
                                # print(0)
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
                                # print(1)
                                search = re.search(
                                    '<a href="\/uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>[A-Z0-9-]+_HUMAN.+?<div class="gene-names"><span class="shortName">(<strong>[a-zA-Z0-9-]+<\/strong>[a-zA-Z0-9-, ]*)<\/span><\/div><\/td><td>?',
                                    url,
                                    re.IGNORECASE,
                                )
                                # print(search.group())
                                # print(search.group(2))
                                if re.findall(GeneID, search.group(2), re.IGNORECASE):
                                    # print(2)
                                    UniProtKB = search.group(1)
                            elif re.search(
                                '<a href="\/uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>[A-Z0-9-]+_MOUSE.+?<div class="gene-names"><span class="shortName">(<strong>[a-zA-Z0-9]+<\/strong>[a-zA-Z0-9-, ]*)<\/span><\/div><\/td><td>?',
                                url,
                                re.IGNORECASE,
                            ):
                                # print(3)
                                search = re.search(
                                    '<a href="\/uniprot\/([A-Z0-9-]+)">[A-Z0-9]+<\/a><\/td><td>[A-Z0-9-]+_MOUSE.+?<div class="gene-names"><span class="shortName">(<strong>[a-zA-Z0-9-]+<\/strong>[a-zA-Z0-9-, ]*)<\/span><\/div><\/td><td>?',
                                    url,
                                    re.IGNORECASE,
                                )
                                if re.findall(GeneID, search.group(2), re.IGNORECASE)[
                                    0
                                ]:
                                    # print(4)
                                    UniProtKB = search.group(1)
                            cc = UniProtKB
                            # print(cc)
                            b = (
                                "https://www.uniprot.org/uniprot/"
                                + cc
                                + "#subcellular_location"
                            )
                            # print(b)
                            bb = str(getHtml(b, time))
                            ddd = re.findall(
                                'class="[a-zA-Z_ ]+"><h6>([a-zA-Z ]+)</h6>', bb
                            )
                            # print(ddd)
                            # if ddd:
                            # 	print(b)
                            if ddd == ["Other locations"]:
                                ddd = re.findall(
                                    'locations*/SL-[0-9]+">([a-zA-Z ]+) </a>', bb
                                )
                                # print(ddd)
                            if GeneID == "Uox":  # uniprot
                                ddd = ["Peroxisome", "Mitochondrion"]  # mouse
                            if GeneID == "NME1-NME2":  # uniprot
                                ddd = ["cytosol"]
                            if GeneID == "CKMT1A" or GeneID == "CKMT1B":  # uniprot
                                ddd = ["mitochondrion"]
                            if GeneID == "TRMT11":  # uniprot
                                ddd = ["cytosol"]
                            dd.extend(ddd)
                            uniprot = "1"
                dd = [i.strip() for i in dd]
                dd = [d for d in dd if d]
                # print(dd)
                if dd:
                    # print(dd)
                    # for i in dd:
                    # 	print(i)
                    # print('he')
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
                        ee = re.findall(r"([A-Za-z ]+)", ee[0])  # ?????!!!!!
                        uniprot = ""
                    SubUnLoc = []
                    x = 0
                    for x in range(len(dd)):
                        if (
                            dd[x]
                            and not dd[x] in OtherLocations
                            and not "GO" in dd[x]
                            and not "PubMed" in dd[x]
                        ):  # ['3.5.1.6']
                            ff = ee[x].split("{")[0]
                            if biocyc:
                                ff = re.findall(r"([A-Za-z0-9 ,;\-\(\)]+)", ff)[0]
                                # print(ff)
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
                                # print(ff)
                                ff = multiple_replace(
                                    dictt, ff.split("{")[0]
                                )  ### !!!!!!??????
                            # print(ff)
                            # print(dictt)
                            # print(ff, dictt.get(ff))
                            if re.findall("Dendriti", ff, re.IGNORECASE):
                                ff = "Dendrite"
                            if (
                                re.match("membrane", ff, re.IGNORECASE)
                                or re.findall("plasma membrane", ff, re.IGNORECASE)
                                or re.findall(
                                    "integral component of membrane", ff, re.IGNORECASE
                                )
                            ):
                                ff = "Plasma membrane"
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
                            # print(ff,list(set(dictt.values())))
                            if hh == 1 and not ff.lower() in list(
                                set(dictt.keys())
                            ):  # important
                                print(ff)
                                # print(sorted(set(dd)))
                                ff = "Cytosol"
                            # if hh == 0 and not ff in list(set(dictt.values())): #do not know if this is nessecary
                            # 	print(1240, ff)
                            # 	#print(sorted(set(dd)))
                            # 	ff = 'Cytosol'
                            if re.findall("\\\\n", ff) or re.findall(
                                "\.", ff
                            ):  # important
                                # print(ff)
                                # print(sorted(set(dd)))
                                ff = "Cytosol"
                                # print(ff)
                                # print()
                            # if hh == 1 and not ff in list(set(dictt.keys())):
                            # 	print(ff)
                            # 	print(sorted(set(dd)))
                    d = sorted(set(SubUnLoc))
                    # print(d)
                    LocList = [
                        "Extracellular",
                        "Peroxisome",
                        "Mitochondria",
                        "Cytosol",
                        "Lysosome",
                        "Endoplasmic reticulum",
                        "Golgi apparatus",
                        "Nucleus",
                        "Inner mitochondria",
                    ]
                    # print([x for x in d if x not in LocList])
                    m = len(
                        a
                    )  # once it is defined a cellular location the process stops because is assumed that all the subunit of the same complex are in the same place
                    # if not dd:
                    # 	print(re.sub(r"\*[0-9]+" , "", a[m]))
                    m = m + 1
                    # if not d and gpr or genelist1 or genelist2:
                    if not gpr in gprgpr:
                        # print()
                        # print('No location was found')
                        # print(gpr, genelist1, genelist2)
                        # print(genelist1,)
                        # print(genelist2)
                        # print()
                        gprgpr += gpr + "\n"
                    d = [
                        "Cytosol"
                    ]  # by default, if there is not anotated location, the reaction is located into the cytoso
                    # d = ['']
                    p = 0
                    # ppList.extend(p)
                    Locations = Locations + [d][0]
                # print 'not dd'
            LocationList = LocationList + [d]
            # print(LocationList)
            # print str(n)+'_'+str(m)+'_000'
            n = n + 1
        Locations = list(set(Locations))
        l = 0
        RuleLoc = {}
        RuleLoc2 = {}
        RuleLoc3 = {}
        RuleLoc4 = {}
        # print()
        # if not RuleLoc4 == {'': ['']}:
        # 	print()
        # 	print(RuleLoc , RuleLoc2,RuleLoc3,RuleLoc4)
        # 	print()
        return (
            RuleLoc,
            RuleLoc2,
            RuleLoc3,
            RuleLoc4,
            p,
        )  # p only with programa_3_2  just a single number and not a list because of line 1256 ca 					m = len(a) # once it is defined a cellular location the process stops because is assumed that all the subunit of the same complex are in the same place
        # print(0, RuleLoc , RuleLoc2)
    except Exception as e:
        print(traceback.format_exc())
        # RuleLoc = {}
        # RuleLoc2 = {}
        # RuleLoc['Cytosol'] = ''
        # RuleLoc2['Cytosol'] = ''
        # print('Exception "'+str(e)+'" in getLocation with following arguments:')
        # print('gpr: '+str(gpr))
        # print('genelist1: '+str(genelist1))
        # print('genelist2: '+str(genelist2))
        # print(RuleLoc , RuleLoc2)
        return RuleLoc, RuleLoc2, RuleLoc3, RuleLoc4


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


def getFormula(page, time, EF, specialCompounds, RxnID):
    try:
        # This function now expects a flat file string instead of HTML
        formula = ""

        # Extract formula from FORMULA field
        formula_match = re.search(r"FORMULA\s+([A-Za-z0-9()]+)", page)
        if formula_match:
            formula = formula_match.group(1)

        # Fallback to COMPOSITION for glycans
        if not formula:
            composition_match = re.search(
                r"COMPOSITION\s+([\s\S]+?)(?=\n[A-Z]|$)", page
            )
            if composition_match:
                formula = composition_match.group(1).strip()

        # If there is a (group)n in the formula, we save the compound ID
        if formula and ")n" in formula:
            entry_match = re.search(r"ENTRY\s+([CG][0-9]+)", page)
            if entry_match:
                ident = entry_match.group(1)
                with open(specialCompounds, "a") as f:
                    f.write(ident + "\n")

        # The rest of the logic for parsing and normalizing the formula string
        if formula:
            CoreForm = re.sub(r"\(", "", re.sub(r"\)n[0-9\-]*", "", formula))
            AtomList = sorted(set(re.findall(r"[A-Z][a-z]?", CoreForm)))
            NewForm = []
            for atom in AtomList:
                # Count atoms with explicit numbers
                count = sum(int(n) for n in re.findall(atom + r"([0-9]+)", CoreForm))
                # Count atoms without explicit numbers
                count += len(re.findall(atom + r"(?![0-9])", CoreForm))
                if count > 0:
                    NewForm.append(f"{atom}{count}")

            return "".join(NewForm)

        # Fallback to external file if no formula found in entry
        entry_match = re.search(r"ENTRY\s+([CG][0-9]+)", page)
        if entry_match:
            ident = entry_match.group(1)
            for line in EF:
                if line.startswith(ident):
                    return line.split("\t")[2]

        return ""

    except Exception as e:
        print(traceback.format_exc())
        print(f'Exception "{e}" in getFormula')
        return ""


def parse_kegg_flat_file(flat_file_text):
    """
    Parse KEGG flat file format (from REST API) into a dictionary.

    Parameters:
    -----------
    flat_file_text : str
        Raw KEGG flat file format text from REST API

    Returns:
    --------
    dict : Dictionary with field names as keys and their content as values
    """
    lines = flat_file_text.strip().split("\n")
    fields = {}
    current_field = None
    current_value = []

    for line in lines:
        # New field starts if line doesn't begin with whitespace
        if line and not line[0].isspace():
            # Save previous field
            if current_field:
                fields[current_field] = "\n".join(current_value)

            # Parse new field
            field_name = line[:12].strip()
            field_value = line[12:].strip()

            if field_name:
                current_field = field_name
                current_value = [field_value] if field_value else []
            else:
                current_field = None
                current_value = []
        else:
            # Continuation of current field
            if current_field:
                current_value.append(line.strip())

    # Don't forget last field
    if current_field:
        fields[current_field] = "\n".join(current_value)

    return fields


def getCompParamFromRestAPI(flat_file_text, ident, time, EF, specialCompounds, RxnID):
    """
    Extract compound parameters from KEGG REST API flat file format.
    This is a new implementation that works directly with REST API data.

    Parameters:
    -----------
    flat_file_text : str
        Raw KEGG flat file format text from REST API
    ident : str
        Compound identifier (e.g., 'C00001')
    time : int
        Timeout for external API calls
    EF : list
        Extra formulas list
    specialCompounds : str
        Path to special compounds file
    RxnID : str
        Reaction ID (for context)

    Returns:
    --------
    tuple : Same format as getCompParam
    """
    try:
        defaultdict = recondict()

        # Parse flat file
        fields = parse_kegg_flat_file(flat_file_text)

        # Initialize return values
        urls01 = ident  # principal identifier
        urls02 = ident  # secondary identifier
        urls3 = ""  # Name
        urls21 = ""  # Formula1
        urls22 = ""  # Formula2
        urls4 = []  # PubChem SID
        urls4aa = ""  # PubChem CID
        urls5 = []  # ChEBI
        urls6 = []  # LIPIDMAPS
        urls7 = []  # LipidBank
        urls8 = []  # GlycomeDB
        urls9 = []  # JCGGDB
        urls11 = "0"  # charge
        urls21a = ""  # adjusted formula
        urls21aa = ""  # reformulated formula
        inchikey = ""
        inchi = ""
        urls1 = []  # Associated reactions

        # Extract NAME
        if "NAME" in fields:
            # Take first name, remove semicolon if present
            names = fields["NAME"].split(";")
            urls3 = names[0].strip()
        else:
            urls3 = ident  # Fallback to identifier

        # Extract FORMULA
        if "FORMULA" in fields:
            urls21 = fields["FORMULA"].strip()
            urls22 = urls21
        else:
            # Fallback to EF file or empty
            for line in EF:
                if line.startswith(ident):
                    urls21 = line.split("\t")[2] if len(line.split("\t")) > 2 else ""
                    urls22 = urls21
                    break

        # Handle COMPOSITION for glycans
        if "COMPOSITION" in fields:
            composition = fields["COMPOSITION"]
            # For glycans, combine name with composition
            if not urls3 or urls3 == ident:
                urls3 = f"{ident} ({composition})"
            else:
                urls3 = f"{urls3} ({composition})"

        # Check for REMARK "Same as" (glycan -> compound mapping)
        if "REMARK" in fields:
            same_as_match = re.search(r"Same as:\s*([CDG][0-9]+)", fields["REMARK"])
            if same_as_match:
                urls01 = same_as_match.group(1)
                urls02 = ident
                # Fetch the primary compound data
                try:
                    primary_url = f"https://rest.kegg.jp/get/{urls01}"
                    primary_response = urllib.request.urlopen(primary_url).read()
                    primary_text = (
                        primary_response.decode("utf-8")
                        if isinstance(primary_response, bytes)
                        else primary_response
                    )
                    primary_fields = parse_kegg_flat_file(primary_text)

                    if "NAME" in primary_fields:
                        urls3 = primary_fields["NAME"].split(";")[0].strip()
                    if "FORMULA" in primary_fields:
                        urls21 = primary_fields["FORMULA"].strip()
                except Exception as e:
                    print(f"Warning: Could not fetch primary compound {urls01}: {e}")

        # Extract DBLINKS
        if "DBLINKS" in fields:
            dblinks_text = fields["DBLINKS"]

            # PubChem
            pubchem_matches = re.findall(r"PubChem:\s*(\d+)", dblinks_text)
            urls4 = pubchem_matches

            # ChEBI
            chebi_matches = re.findall(r"ChEBI:\s*(\d+)", dblinks_text)
            urls5 = chebi_matches

            # LIPID MAPS
            lipidmaps_matches = re.findall(
                r"LIPID MAPS:\s*([A-Z]+[0-9]+)", dblinks_text
            )
            urls6 = lipidmaps_matches

            # LipidBank
            lipidbank_matches = re.findall(r"LipidBank:\s*([A-Z0-9]+)", dblinks_text)
            urls7 = lipidbank_matches

            # GlycomeDB
            glycomedb_matches = re.findall(r"GlycomeDB:\s*([A-Z0-9]+)", dblinks_text)
            urls8 = glycomedb_matches

            # JCGGDB
            jcgg_matches = re.findall(r"(JCGG-[A-Z0-9]+)", dblinks_text)
            urls9 = jcgg_matches

        # Try to get PubChem CID and InChI info
        try:
            if urls3 and urls3 != ident:
                urls4aa_list = pcp.get_cids(urls3, "name")
                if urls4aa_list:
                    urls4aa = str(urls4aa_list[0])
                    compound = pcp.Compound.from_cid(urls4aa_list[0])
                    inchikey = compound.inchikey
                    inchi = compound.inchi
        except Exception:
            pass

        # Handle charge from defaultdict
        urlss = (
            [urls01, urls02]
            + urls4
            + ([urls4aa] if urls4aa else [])
            + urls5
            + urls6
            + urls7
            + urls8
            + urls9
        )
        urls21a = urls21

        # helper: parse formula into element->count dict
        def parse_formula_to_dict(f):
            d = {}
            try:
                if not f:
                    return d
                # remove brackets and spaces, keep + and - if present
                fclean = re.sub(r"[()\[\]\s]", "", f)
                parts = re.findall(r"([A-Z][a-z]?)(\+?\-?\d*)", fclean)
                for el, num in parts:
                    if num is None or num == "":
                        n = 1
                    else:
                        try:
                            n = int(num)
                        except Exception:
                            n = 1
                    d[el] = d.get(el, 0) + n
            except Exception:
                return {}
            return d

        for j in urlss:
            try:
                if str(j) in defaultdict and defaultdict[str(j)]:
                    urls10a = defaultdict[str(j)]
                    urls11 = urls10a[1]
                    known_formula = urls10a[0]
                    # compare parsed formulas (ignore ordering)
                    fk = parse_formula_to_dict(known_formula)
                    fc = parse_formula_to_dict(urls21)
                    if fk and fc and fk == fc:
                        # matched composition -> prefer KEGG known formula
                        urls21a = known_formula
                        if urls21 == urls22:
                            urls22 = urls21a
                        else:
                            urls21 = urls21a
            except Exception:
                # be defensive: continue without failing
                continue

        # Handle 'R' group in formula
        bb = [x for x in [urls01, urls02] if recondict()[x]]
        if bb and "R" in urls21:
            if bb[0] in defaultdict:
                urls211 = re.findall("[A-Z]", urls21)
                urls10a = defaultdict[bb[0]]
                urls10a1 = re.findall("[A-Z]", urls10a[0])
                if (
                    all(v in urls211 for v in urls10a1)
                    and all(v in urls10a1 for v in urls211)
                    and re.findall("C([0-9]*)", urls10a[0])
                    > re.findall("C([0-9]*)", urls21)
                ):
                    if len(urls211) == len(urls10a1):
                        urls11 = urls10a[1]
                        # prefer the known KEGG formula (atom2 not available here)
                        urls21a = urls10a[0]
                        if urls21 == urls22:
                            urls21, urls22 = urls21a, urls21a
                        else:
                            urls21 = urls21a
                    else:
                        urls21 = urls10a[0]

        # Handle formulas with groups
        if urls21 and "(" in urls21:
            Path = "gly2"
            with open(Path, "a") as Formula_file:
                try:
                    urls21aa = ReformulationFormula(
                        [x for x in [urls21, urls22] if "(" in x][0]
                    )
                except Exception:
                    urls21aa = ""
                Formula_file.write(
                    ident + "\t" + urls21 + "\t" + urls22 + "\t" + urls21aa + "\n"
                )

        if urls21a == "":
            urls21a = urls21

        # Build urls1 (associated reactions) from REACTION field
        if "REACTION" in fields:
            reactions = fields["REACTION"].split()
            urls1 = [
                (f"http://www.genome.jp/entry/{rxn}", rxn)
                for rxn in reactions
                if rxn.startswith("R")
            ]
            urls1 = [list(zip([x[0]], [x[1]])) for x in urls1]

        return (
            list(zip([urls01], [urls02])),
            list(zip([urls21a], [urls22])),
            urls1,
            [urls3],
            urls4,
            urls5,
            urls6,
            urls7,
            urls8,
            urls9,
            urls11,
            urls21aa,
            inchikey,
            inchi,
            urls4aa,
        )

    except Exception as e:
        print(traceback.format_exc())
        print('Exception "' + str(e) + '" in getCompParamFromRestAPI')
        print(ident)
        return ""


""""Compound: Extract the links from a HTML page or REST API data"""


def getCompParam(page, ident, time, EF, specialCompounds, RxnID):
    try:
        # page can be either HTML (legacy) or REST API formatted data
        defaultdict = recondict()

        # Extract reaction links
        urls1 = [
            list(zip([x[0].replace('href="', "http://www.genome.jp")], [x[1]]))
            for x in re.findall(r'(href=\"[^\'" >]+).>(\w+)<', page)
        ]  # R removed

        # If the compound is a Glycan ...
        GlyTest = re.findall(r"KEGG GLYCAN: G[0-9]+", page)
        if GlyTest:
            urls01 = re.findall(
                'Remark[\S\s]+Same as:[\S\s]+">(C[0-9]+)<\/a>', page
            )  # principal identifier
            if urls01:  # If the compound is a Glycan and it has an alternative compound
                urls02 = ident  # secondary identifier
                urls22 = getFormula(page, time, EF, specialCompounds, None)
                urls8 = re.findall(
                    'GlycomeDB[\S\s]+?glycomeId=([A-Z0-9]+)"', page
                )  # GlycomeDB
                urls9 = re.findall(
                    "GlycomeDB[\S\s]+?>(JCGG-[A-Z0-9]+)<", page
                )  # JCGGDB
                urls01 = urls01[0]
                page2 = getHtml(
                    "http://www.genome.jp/dbget-bin/www_bget?cpd:" + urls01, time
                )  # page of the primary compoun
                urls3 = re.findall(
                    'Name<.span><.th>.n<td class="td21 defd"><div class="cel"><div class="cel">(.+?);?<br>',
                    str(page2),
                    re.DOTALL,
                )[0].replace(
                    "\\", ""
                )  # # #Name
                urls21 = getFormula(page2, time, EF, specialCompounds, None)
                page = page2
            else:  # If the compound is a Glycan and it has not an alternative compound
                urls01 = ident  # principal identifier
                urls02 = ident  # secondary identifier
                urls21 = getFormula(
                    page, time, EF, specialCompounds, None
                )  # principal formula¨
                urls22 = urls21  # secondary formula
                dd = re.findall(
                    "Composition.+?(\(.+?)<", page, re.DOTALL
                )  # name of glycan
                ddd = re.findall(
                    'Name<.span><.th>.n<td class="td21 defd"><div class="cel"><div class="cel">(.+?)[;<]',
                    page,
                    re.DOTALL,
                )  # name of glycan
                if ddd:
                    dd = ddd[0] + " (" + dd[0] + ")"
                else:
                    dd = dd[0] + " (" + urls01 + ")"
                urls3 = dd.replace("\\", "")
                urls8 = re.findall(
                    'GlycomeDB[\S\s]+?glycomeId=([A-Z0-9]+)"', page
                )  # GlycomeDB
                urls9 = re.findall(
                    "GlycomeDB[\S\s]+?>(JCGG-[A-Z0-9]+)<", page
                )  # JCGGDB
        else:  # If the compound is not a Glycan
            urls01 = ident  # principal identifier
            urls02 = ident  # secondary identifier
            urls21 = getFormula(page, time, EF, specialCompounds, None)
            urls22 = urls21  # secondary formula
            # Extract compound name with better error handling
            name_match = re.findall(
                'Name<.span><.th>.n<td class="td21 defd"><div class="cel"><div class="cel">(.+?);?<br>',
                page,
                re.DOTALL,
            )
            if name_match:
                urls3 = name_match[0].replace("\\", "")
            else:
                # Fallback: try to extract from title or use identifier
                title_match = re.search(r"<title>KEGG [A-Z]+: ([^<]+)</title>", page)
                if title_match:
                    urls3 = title_match.group(1).strip()
                else:
                    urls3 = ident  # Use identifier as last resort
            urls8 = re.findall(
                'GlycomeDB[\S\s]+?glycomeId=([A-Z0-9]+)"', page
            )  # GlycomeDB
            urls9 = re.findall("GlycomeDB[\S\s]+?>(JCGG-[A-Z0-9]+)<", page)  # JCGGDB
        urls3 = urls3.replace("&gt;", ">")
        page = str(page)
        # PubChem sid
        compound = ""
        inchikey = ""
        inchi = ""
        urls4 = re.findall(r"PubChem[\S\s]+?sid=([0-9]+)", page)
        try:
            urls4aa = pcp.get_cids(urls01, "name")
        except Exception:
            urls4aa = list()
        if urls4aa:
            compound = pcp.Compound.from_cid(urls4aa)
            inchikey = compound.inchikey
            inchi = compound.inchi
        # CheBI
        urls5 = re.findall("chebiId=CHEBI:([0-9]+)", page)
        if len(urls5) > 1:
            urls5 = [urls5[-1]]
        # LIPIDMAPS
        urls6 = re.findall("LMID=([A-Z]+[0-9]+)", page)
        if len(urls6) > 1:
            urls6 = [urls6[-1]]
        # LipidBank
        urls7 = re.findall('LipidBank[\S\s]+?id=([A-Z0-9]+)"', page)
        del GlyTest
        urlss = (
            [urls01]
            + [urls02]
            + urls4
            + urls4aa
            + urls5
            + urls6
            + urls7
            + urls8
            + urls9
        )
        if urls4aa:
            urls4aa = urls4aa[0]
        else:
            urls4aa = ""
        urls11 = "0"
        urls21a = ""
        urls21aa = ""

        # helper: parse formula into element->count dict (reused from above)
        def parse_formula_to_dict(f):
            d = {}
            try:
                if not f:
                    return d
                fclean = re.sub(r"[()\[\]\s]", "", f)
                parts = re.findall(r"([A-Z][a-z]?)(\+?\-?\d*)", fclean)
                for el, num in parts:
                    if num is None or num == "":
                        n = 1
                    else:
                        try:
                            n = int(num)
                        except Exception:
                            n = 1
                    d[el] = d.get(el, 0) + n
            except Exception:
                return {}
            return d

        for j in urlss:
            try:
                if str(j) in defaultdict and defaultdict[str(j)]:
                    urls10a = defaultdict[str(j)]
                    urls11 = urls10a[1]
                    known_formula = urls10a[0]
                    fk = parse_formula_to_dict(known_formula)
                    fc = parse_formula_to_dict(urls21)
                    if fk and fc and fk == fc:
                        urls21a = known_formula
                        if urls21 == urls22:
                            urls22 = urls21a
            except Exception:
                continue
        bb = [x for x in [urls01, urls02] if recondict()[x]]
        if bb and "R" in urls21:
            if bb[0] in defaultdict:
                urls211 = re.findall("[A-Z]", urls21)
                # 		if urls01 ==urls02 and urls01 in defaultdict:
                urls10a = defaultdict[bb[0]]
                urls10a1 = re.findall("[A-Z]", urls10a[0])
                if (
                    all(v in urls211 for v in urls10a1)
                    and all(v in urls10a1 for v in urls211)
                    and re.findall("C([0-9]*)", urls10a[0])
                    > re.findall("C([0-9]*)", urls21)
                ):
                    if len(urls211) == len(urls10a1):
                        urls11 = urls10a[1]
                        # prefer KEGG known formula (atom2 unavailable)
                        urls21a = urls10a[0]
                        if urls21 == urls22:
                            urls21, urls22 = urls21a, urls21a
                        else:
                            urls21 = urls21a
                    else:
                        urls21 = urls10a[0]
        # 		urls21a,urls22=g,g
        if urls21:
            if "(" in urls21:
                # if "(" in Formula: # in case glycan
                Path = "gly2"
                with open(Path, "a") as Formula_file:
                    try:
                        urls21aa = ReformulationFormula(
                            [x for x in [urls21, urls22] if "(" in x][0]
                        )
                    except Exception:
                        urls21aa = ""
                    Formula_file.write(
                        ident + "\t" + urls21 + "\t" + urls22 + "\t" + urls21aa + "\n"
                    )
        if urls21a == "":
            urls21a = urls21
        # urls22 = urls21a
        urls4aa = str(urls4aa)
        return (
            list(zip([urls01], [urls02])),
            list(zip([urls21a], [urls22])),
            urls1,
            [urls3],
            urls4,
            urls5,
            urls6,
            urls7,
            urls8,
            urls9,
            urls11,
            urls21aa,
            inchikey,
            inchi,
            urls4aa,
        )
        # (ID1,ID2), (Formula1, Formula2), Associated reactions, Name, PubChem, CheBI, LIPIDMAPS, LipidBank, GlycomeDB, JCGGDB,charge,inchikey,inchi
    except Exception as e:
        print(traceback.format_exc())
        print('Exception "' + str(e) + '" in getCompParam')
        print(ident)
        return ""
