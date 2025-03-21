from cobra.io import read_sbml_model,write_sbml_model
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

#import the functions
from functions.gpr.gpr_def import getGPR, setup_biocyc_session

from functions.function_bm_gdb import *
from functions.equations_bm_gdb import *


#retrieve website with function from getgpr

def get_html(request_url: str, session: Optional[requests.Session] = None) -> str:
    """Fetch an html by perfoming a GET HTTPS request, maybe with session."""
    if session is not None:
        return session.get(request_url).text
    else:
        return requests.get(request_url).text


def create_dict(gene,e):
    n=''
    global my_dict
    if 'my_dict' not in globals():
        my_dict = {}
    elif gene and e:
        my_dict[gene]=e 
    elif gene and not e:
        n= my_dict[gene]
    return my_dict,n

#create a function that takes a string and a dictionary and replaces the keys in the string with the corresponding values in the dictionary
def multiple_replace_new(dictionary, text):
    #Remove any leading or trailing whitespace from the text
    text = text.strip()
    #put everything in lowercase
    text = text.lower()
    #Search for the text in the dictionary
    for key in dictionary:
        #If the key is found in the text, replace it with the corresponding value
        if key == text:
            text = text.replace(key, dictionary[key])
        if key != text:
            text = text
    return text

def getLocationnew(
    gpr,
    genelist1,
    genelist2,
    impose_locations,
    location_dict_file,
    session: Optional[requests.Session] = None,
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
        if impose_locations == 1:
            Var = open(location_dict_file, "rb")  # Endo1b_variables.pkl
            dictt = pickle.load(Var)
            hh = 0  # dict and not all
            # print(location_dict_file + " is used")
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
                bb = str(urllib.request.urlopen(b).read())
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
                        index = [x.upper() for x in genelist11].index(
                            re.sub(r"\*[0-9]+", "", a[m].upper())
                        )
                        b = (
                            "http://biocyc.org/gene?orgid=META&id="
                            + genelist2[index].upper()
                        )  # location in BioCyc human #########################
                        # bb = str(getHtml(b, session))
                        try:
                            bb = str(urllib.request.urlopen(b).read())
                        except urllib.error.HTTPError:
                            bb = ''
                        if not bb:
                            b = (
                                "https://biocyc.org/gene?orgid=HUMAN&id="
                                + genelist2[index].upper()
                            )  # location in BioCyc human #########################
                            #bb = str(getHtml(b, session))
                            try:
                                bb = str(urllib.request.urlopen(b).read())
                            except urllib.error.HTTPError:
                                bb = ''
                            # bb = str(urllib.request.urlopen(b).read())

                        if bb:
                            if re.search("Location", bb, flags=re.DOTALL):
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
                                    ddd = [
                                        re.sub("^ ", "", i) for i in ddd[0].split(",")
                                    ]
                                dd.extend(ddd)

                            biocyc = "1"

                            if bb:  # location in Uniprot

                                biocyc = ""
                                if re.search("uniprot/([A-Z0-9]+)", bb):
                                    cc = re.search("uniprot/([A-Z0-9]+)", bb).group(1)
                                    b = (
                                        "https://rest.uniprot.org/uniprotkb/"
                                        + cc
                                        + ".txt"
                                    )
                                    # bb = str(getHtml(b, session))
                                    bb = str(urllib.request.urlopen(b).read())
                                    ddd = re.findall("GO:[0-9]+.+?C:(.+?);", bb)
                                    ddd = [x for x in ddd if "GO" not in x]
                                    dd.extend(ddd)
                                    uniprot = "1"

                        #try to retrieve website with function from getgpr

                        #for humancyc
                        humancyc="https://biocyc.org/gene?orgid=HUMAN&id="+ genelist2[index].upper()
                        #for metacyc
                        metacyc="http://biocyc.org/gene?orgid=META&id="+ genelist2[index].upper()
                        btry=get_html(humancyc, session)  
                        bb_human = str(btry)
                        btry=get_html(metacyc, session)
                        bb_meta = str(btry)

                        # Extract the section between "Locations" and "Reactions"
                        location_section_human = re.search(r"(?<=\nLocations)(.*?)(?=>\nReactions)", bb_human, re.DOTALL)
                        location_section_meta = re.search(r"(?<=\nLocations)(.*?)(?=>\nReactions)", bb_meta, re.DOTALL)

                        if location_section_human:
                            compartments = re.findall(r"[\w\s-]+(?=\s*<a href)", location_section_human.group(0))
                            #remove any \n in the compartments 
                            compartments = [compartment.replace("\n", "") for compartment in compartments]
                            #remove any leading or trailing whitespace
                            compartments = [compartment.strip() for compartment in compartments]
                            dd.extend(compartments)

                        elif location_section_meta:
                            compartments = re.findall(r"[\w\s-]+(?=\s*<a href)", location_section_human.group(0))
                            #remove any \n in the compartments 
                            compartments = [compartment.replace("\n", "") for compartment in compartments]
                            #remove any leading or trailing whitespace
                            compartments = [compartment.strip() for compartment in compartments]
                            dd.extend(compartments)
                                                
                        dd = list(set(dd)) #remove duplicates

                        if not dd:
                            UniProtKB = ""
                            GeneID = re.sub(r"\*[0-9-]+", "", a[m])
                            https = (
                                "https://www.uniprot.org/uniprot/?query="
                                + GeneID
                                + "&sort=score"
                            )
                            # url = str(
                            #     getHtml(https, session)
                            # )  # str(urllib.request.urlopen(https).read())
                            url = str(urllib.request.urlopen(https).read())
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
                            # bb = str(getHtml(b, session))
                            bb = str(urllib.request.urlopen(b).read())
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
                                ff = multiple_replace_new(dictt, ff.split("{")[0]) #replace the keys in the string with the corresponding values in the dictionary


                            if re.findall("Dendriti", ff, re.IGNORECASE):
                                ff = "Dendrite"
                            if (
                                re.match("membrane", ff, re.IGNORECASE)
                                or re.findall("plasma membrane", ff, re.IGNORECASE)
                                or re.findall("integral component of membrane", ff, re.IGNORECASE)
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
                                SubUnLoc = SubUnLoc + [ff]
                                Locations = Locations + [ff]
                            

                            #transform everything in lowercase
                            

                    d = sorted(set(SubUnLoc))

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
                        iiii2 = re.search(
                            "gene=([A-Z0-9]+)",
                            str(
                                urllib.request.urlopen(
                                    "https://www.genome.jp/dbget-bin/www_bget?hsa+"
                                    + LocGPR2[i]
                                ).read()
                            ),
                        )  # .group(1)
                        if not iiii2:
                            iiii2 = re.search(
                                "(ENSG[0-9]+)",
                                str(
                                    urllib.request.urlopen(
                                        "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g="
                                        + LocGPR2[i]
                                    ).read()
                                ),
                            )
                        if iiii2:
                            iiii = iiii2.group(1)
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


