# -*- coding: utf-8 -*-

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
import sys

# Determine the current file's directory and the project root.
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.join(current_dir, "..")

# Add the project root to sys.path to access top-level folders like 'functions' and 'models'
if project_root not in sys.path:
    sys.path.append(project_root)

from functions.function_bm_gdb import *
from functions.equations_bm_gdb import *

from functions.gpr.gpr_def import getGPR, setup_biocyc_session

# 👇 import the addtional function
from functions.gpr.get_location_def import getLocationnew as getLocation
from functions.ensembl_client import fetch_ensembl_annotations
from datetime import datetime

if __name__ == "__main__":

    USE_CHECKPOINT = False

    session = setup_biocyc_session()

    # Ensembl annotation cache (persisted across runs to speed up lookups)
    ensembl_cache_file = os.path.join(project_root, "files", "ensembl_cache.pkl")
    ensembl_cache = {}
    try:
        if os.path.exists(ensembl_cache_file):
            with open(ensembl_cache_file, "rb") as _f:
                ensembl_cache = pickle.load(_f)
                if not isinstance(ensembl_cache, dict):
                    ensembl_cache = dict(ensembl_cache)
            print(f"Loaded ensembl cache with {len(ensembl_cache)} entries")
    except Exception as _e:
        print("Warning: could not load ensembl cache:", _e)

    input_model_path = os.path.join(project_root, "models", "THG-beta1.1.1_251031.xml")
    # input_model_path = os.path.join(project_root, "models", "Human_Database_old.xml")
    model = read_sbml_model(input_model_path)

    # create a date tag like "_yymmdd"
    _date_tag = "_" + datetime.now().strftime("%y%m%d")

    variablesFile = os.path.join(project_root, "files", f"variablesFile{_date_tag}.pkl")
    Output = os.path.join(project_root, "models", f"model_bb{_date_tag}.xml")
    Output2 = os.path.join(project_root, "files", f"Rxn2Fix{_date_tag}.txt")
    Output3 = os.path.join(project_root, "models", f"THG-beta2-tmp{_date_tag}.xml")
    Output4 = os.path.join(project_root, "models", f"THG-beta2{_date_tag}.xml")
    Output5 = os.path.join(project_root, "models", f"THG-beta1{_date_tag}.xml")

    handle = open(Output2, "w")

    # Parameters about the compartments (change according to your model)
    excel_file = os.path.join(
        project_root, "files", "ListOfCompartments_sept2024.xlsx"
    )  # Excel file containing the compartments information
    compartments_sheet_name = "Def-Compartments"  # Sheet name containing the compartments information, where the first column corresponds to the compartment name and the rest of the columns correspond to the mapping of this compartment in different models

    # Column index 2 is for endoA
    column_index = (
        2  # Column with the compartment mapping corresponding to your desired model
    )
    abbreviations_sheet_name = (
        "Endo1a_abb"  # Sheet name containing the abbreviations of the compartments
    )

    # Create a dictionary mapping compartment names to their corresponding values in a model
    location_pkl_file = os.path.join(project_root, "files", "dict_compartments.pkl")
    create_compartments_dict_bm(
        excel_file, compartments_sheet_name, column_index, location_pkl_file
    )  # create the pkl file: dict_compartments.pkl

    # Create a dictionary with the abbreviations of the compartments as keys and the full names as values
    comp_abb_file = os.path.join(
        project_root, "files", "compartments_abbreviations.pkl"
    )
    create_comp_abbreviations_dict_bm(
        excel_file, abbreviations_sheet_name, comp_abb_file
    )  # create the pkl file: compartments_abbreviations.pkl
    comp_dict = pd.read_pickle(comp_abb_file)

    c, comp_dict = compartment_file_to_dict_bm(
        excel_file, compartments_sheet_name, comp_dict, comp_abb_file
    )  # c is the mapping with the full name as keys, abbreviations as values and comp_dict is the mapping with the abbreviations as keys and the full name as values
    variables = defaultdict(list)
    listA, listB, ComptoID2 = list(c.values()), [], {}

    # --- Resume support: if a previous run saved variables or the temporary model, load them
    processed_rxns = set()

    if USE_CHECKPOINT:
        # If the temporary output model exists, load it and skip any reactions already present
        try:
            if os.path.exists(Output3):
                print(
                    f"Found existing temporary model {Output3}, attempting to load to resume..."
                )
                model2 = read_sbml_model(Output3)
                processed_rxns = set([r.id for r in model2.reactions])
                print(
                    f"Resuming: {len(processed_rxns)} reactions already present in {Output3}; these will be skipped."
                )
            else:
                model2 = None
        except Exception as _e:
            print("Warning: could not load existing temporary model to resume:", _e)
            model2 = None

        # If model2 wasn't loaded from the temp file, fall back to copying the original model
        if model2 is None:
            model2 = model.copy()

        # If a variables file exists from a previous run, load it so we reuse cached GPR/location lookups
        try:
            if os.path.exists(variablesFile):
                print(
                    f"Found existing variables file {variablesFile}, attempting to load cached variables..."
                )
                with open(variablesFile, "rb") as savedVariables:
                    loaded = dill.load(savedVariables)
                # ensure defaultdict behavior
                if not isinstance(loaded, defaultdict):
                    loaded = defaultdict(list, loaded)
                variables = loaded
                print(f"Loaded {len(variables)} cached entries from {variablesFile}.")
        except Exception as _e:
            print("Warning: could not load variables file to resume:", _e)
    else:
        # Checkpoints disabled: start fresh
        print(
            "USE_CHECKPOINT is False -> not loading checkpoint files; starting from original model and empty cache."
        )
        model2 = model.copy()
        variables = defaultdict(list)

    n = 1
    eList = []
    RxnList = []

    RxnPath = defaultdict(list)
    for p in model.groups:
        for r in p.members:
            RxnPath[r.id].append(p.id)

    Rxn = defaultdict(list)

    for x in model.reactions:
        Rxn[re.sub("[a-z]+", "", x.reaction)].extend(
            [[x.id, re.findall("[a-z]+", x.reaction)[0]]]
        )

    RxnID = (
        int(re.findall("[0-9]+", sorted([x.id for x in model.reactions])[-1])[0]) + 1
    )

    ListOfMetFrom = {x.id: x.formula for x in model.metabolites}

    Atom_ID = {
        "H2O": "MAM02040",
        "H": "MAM02039",
        "Fe": "MAM01821",
        "X": "MAM01823",
        "R": "MAM03573",
        "Na": "MAM02519",
        "K": "MAM02200",
        "Ca": "MAM01413",
    }

    rr = [x.id for x in model.reactions]
    m2 = [x.id for x in model.metabolites]
    geneList = [x.id for x in model.genes]

    print(f"{len(model.reactions)} reactions initially")
    print(f"{len(model.metabolites)} metabolites initially")
    print(f"{len(model.genes)} genes initially")

    list_of_compartments = []
    model2 = model.copy()

    with model:
        for reac_count, x in enumerate(model.reactions[0:]):

            x2 = model2.reactions.get_by_id(x.id)

            print(f"Reaction({reac_count}/{len(rr)})")
            try:
                listOfgeneList5 = []
                if "ec-code" in x2.annotation:
                    EC = x.annotation["ec-code"]
                else:
                    EC = ""
                if type(EC) == str:
                    EC = [EC]
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
                mb = x2.check_mass_balance()
                bounds = x2.bounds
                gpr = x2.gpr
                annotation2 = x2.annotation
                locations = list(
                    set([model2.metabolites.get_by_id(s).compartment for s in species])
                )

                eq = x2.reaction
                eq_test = re.sub("<=>", "->", re.sub("-->", "->", eq)).strip()
                eq = re.sub(
                    "(^| )[0-9\.]+", "", re.sub("<=>", "->", re.sub("-->", "->", eq))
                ).strip()

                all_met_have_formula_test = min(
                    [1 if ListOfMetFrom[x] else 0 for x in species]
                )
                if all_met_have_formula_test == 1:
                    FromList = list()
                    for y in species:
                        eq = re.sub(y, ListOfMetFrom[y], eq)
                        eq_test = re.sub(y, ListOfMetFrom[y], eq_test)

                if (
                    len(species) > 1
                    and all_met_have_formula_test == 1
                    and not test_reaction_balance(eq_test)[0]
                ):  # first condition avoids exchange and sink reactions, second condition avoid to mass balance eq with some missing formula # Previous, worng: if len(FromList) > 1 and not '' in FromList: # MAR13082

                    MB = mass_balance(eq, "R")

                    NewSpecies = MB[6] + MB[7]  # HH
                    Stch = [-abs(x) for x in MB[0]] + MB[1]
                    if MB[4] and MB[-1] != 2:
                        [
                            x2.add_metabolites(
                                {
                                    model.metabolites.get_by_id(
                                        Atom_ID[MB[4][0]] + locations[0]
                                    ): MB[5][i]
                                }
                            )
                            for i in range(len(MB[4]))
                        ]
                        species = [s.id for s in x2.reactants] + [
                            s.id for s in x2.products
                        ]  #
                    if MB[-1] != 2:  # Original: MB[-2]
                        for y in x2.metabolites:
                            x2.metabolites[y] = Stch[NewSpecies.index(y.formula)]
                            z = sympy.symbols("z")
                            expr = sympy.Eq(
                                x2.metabolites[y] + z, Stch[NewSpecies.index(y.formula)]
                            )
                            sol = sympy.solve(expr)
                            x2.add_metabolites(
                                {
                                    model.metabolites.get_by_id(y.id): round(
                                        float(sol[0]), 1
                                    )
                                }
                            )

                        species3 = [
                            (
                                {s.split(" ")[0]: -1}
                                if len(s.split(" ")) == 1
                                else {s.split(" ")[1]: -abs(float(s.split(" ")[0]))}
                            )
                            for s in re.split(
                                " --> | <=> ", re.sub("[a-z]+", "", x2.reaction)
                            )[0].split(" + ")
                        ] + [
                            (
                                {s.split(" ")[0]: 1}
                                if len(s.split(" ")) == 1
                                else {s.split(" ")[1]: float(s.split(" ")[0])}
                            )
                            for s in re.split(
                                " --> | <=> ", re.sub("[a-z]+", "", x2.reaction)
                            )[1].split(" + ")
                        ]
                        species3 = dict(ChainMap(*species3))

                if (
                    EC[0] and len(locations) == 1
                ):  # first: only reactions with EC can be analyzed, second: transport reactions are not analyzed

                    for ec in EC:

                        print(ec)
                        if not ec in variables:
                            new_gpr = getGPR(ec, session)  # new function

                            if new_gpr[-1]:

                                # Keep only unique sGPRs: if a sGPR is repeated, keep only the first occurrence and remove the rest
                                # Create a copy of the dictionary to iterate over
                                updated_dict = new_gpr[-1].copy()

                                # Create a set to keep track of seen values
                                seen_values = set()

                                # Iterate over the copy of the dictionary
                                for key, value in new_gpr[-1].items():
                                    if value in seen_values:
                                        # If the value has been seen, remove the corresponding key from the original dictionary
                                        updated_dict.pop(key)
                                    else:
                                        # If the value hasn't been seen, add it to the set
                                        seen_values.add(value)

                                new_gpr = new_gpr[:-1] + (updated_dict,)

                                # separate the complete gpr into individual gprs
                                # the individual gprs are the values of the new_gpr[-1] dictionary

                                genelist1 = []
                                genelist2 = []

                                for key, value in new_gpr[
                                    -1
                                ].items():  # separate the complete gpr into individual gprs
                                    gpr = value
                                    print("gpr: ", gpr)

                                    # in new_gpr[1] and new_gpr[2] the gene names and biocyc ids are stored, but we need to keep only the ones that appear in the gpr
                                    for gene in new_gpr[1]:
                                        if gene in gpr:
                                            genelist1.append(gene)
                                            # find position of gene in new_gpr[1]
                                            index = new_gpr[1].index(gene)
                                            # add the corresponding biocyc id to genelist2
                                            genelist2.append(new_gpr[2][index])

                                    # Prefetch Ensembl annotations for the genes in this
                                    # GPR to avoid many per-gene web calls inside
                                    # getLocation/getLocationnew. We update a local
                                    # cache and persist it to disk so subsequent runs
                                    # are faster.
                                    try:
                                        if genelist1:
                                            unseen = [
                                                g
                                                for g in genelist1
                                                if g and g not in ensembl_cache
                                            ]
                                            if unseen:
                                                # fetch in batches using the shared client
                                                try:
                                                    fetched = fetch_ensembl_annotations(
                                                        unseen,
                                                        batch_size=50,
                                                        max_workers=10,
                                                    )
                                                    if fetched:
                                                        ensembl_cache.update(fetched)
                                                        try:
                                                            with open(
                                                                ensembl_cache_file, "wb"
                                                            ) as _f:
                                                                pickle.dump(
                                                                    ensembl_cache,
                                                                    _f,
                                                                    protocol=pickle.HIGHEST_PROTOCOL,
                                                                )
                                                        except Exception as _e:
                                                            print(
                                                                "Warning: could not persist ensembl cache:",
                                                                _e,
                                                            )
                                                except Exception as _e:
                                                    print(
                                                        "Warning: ensembl prefetch failed:",
                                                        _e,
                                                    )
                                    except Exception:
                                        # If anything goes wrong, proceed without cache
                                        pass

                                    # call the getLocation function, provide the cache so
                                    # it can avoid remote Ensembl lookups
                                    new_locations = getLocation(
                                        gpr,
                                        genelist1,
                                        genelist2,
                                        1,
                                        location_pkl_file,
                                        session,
                                        ensembl_cache,
                                    )
                                    print("new_locations: ", new_locations)
                                    genelist1 = []
                                    genelist2 = []
                                    variables[ec].extend(new_gpr)
                                    variables[ec].extend(new_locations)

                        else:
                            new_gpr = variables[ec][0:6]
                            print(new_gpr)

                            new_locations = variables[ec][6:]
                            print(new_locations)
                        if new_gpr[-1]:
                            listOfgeneList5.append(variables[ec][6:])
                    if listOfgeneList5 and not x2.reaction in RxnList:
                        print(listOfgeneList5)
                        geneList5 = meltGeneList(listOfgeneList5)
                        variables[x2.id] = geneList5
                        for CSL in [*variables[x2.id][0]]:

                            CSL2 = CSL

                            if not CSL in c:
                                print()
                                print(0, CSL)
                                print()
                                CSL = "cytosol"

                            listAB = list(set(listB + listA))
                            if CSL in c:
                                ID = c[CSL]
                            elif (
                                len(re.sub(" $", "", re.sub("^ ", "", CSL)).split(" "))
                                > 1
                            ):
                                ID = (
                                    (
                                        CSL.strip().split(" ")[0][0]
                                        + CSL.strip().split(" ")[1][0]
                                    )
                                    .lower()
                                    .replace(" ", "")
                                )
                            else:
                                if len(CSL.split(" ")) > 1:
                                    ID = CSL[0:3].lower().replace(" ", "")
                                else:
                                    ID = CSL[0:2].lower().replace(" ", "")
                            if not c.get(CSL) and ID in listAB:
                                r = re.compile(ID)
                                ID = ID + str(len(list(filter(r.match, listAB))) + 1)
                            ComptoID2[CSL] = ""
                            ComptoID2[CSL] += ID
                            listA.append(ID)
                            listB.append(ID)

                            if not ID in [
                                rxn[1]
                                for rxn in Rxn[re.sub("[a-z]+[0-9]*", "", x2.reaction)]
                            ]:
                                reaction2 = Reaction("MAR" + str(RxnID + n))
                                reaction2.lower_bound = bounds[0]
                                reaction2.upper_bound = bounds[1]
                                for species2 in [
                                    [re.sub("[a-z][0-9]*", ComptoID2[CSL], x), x]
                                    for x in species
                                ]:
                                    try:
                                        m = model2.metabolites.get_by_id(species2[0])
                                    except KeyError:
                                        m = Metabolite(
                                            species2[0],
                                            charge=model2.metabolites.get_by_id(
                                                species2[1]
                                            ).charge,
                                            formula=model2.metabolites.get_by_id(
                                                species2[1]
                                            ).formula,
                                            name=model2.metabolites.get_by_id(
                                                species2[1]
                                            ).name,
                                            compartment=ID,
                                        )
                                        m.annotation = model2.metabolites.get_by_id(
                                            species2[1]
                                        ).annotation

                                    # Identify if the metabolite exists and if not add it to the model
                                    if not m.id in model2.metabolites:
                                        model2.metabolites.add(
                                            m
                                        )  # add metabolite m to the model

                                    # Add compartment name if not exists
                                    if not model2.compartments[ID]:
                                        model2.compartments[ID] = CSL2

                                    # Add the metabolite "species2" to the reaction x
                                    reaction2.add_metabolites(
                                        {
                                            m: species3[
                                                re.sub("[a-z]+[0-9]*", "", species2[0])
                                            ]
                                        }
                                    )

                                if (
                                    "or" in variables[x2.id][2][CSL2]
                                    and "and" in variables[x2.id][2][CSL2]
                                ):
                                    reaction2.gene_reaction_rule = " and ".join(
                                        [
                                            "(" + v + ")"
                                            for v in variables[x2.id][2][CSL2].split(
                                                " and "
                                            )
                                        ]
                                    )
                                else:
                                    reaction2.gene_reaction_rule = variables[x2.id][2][
                                        CSL2
                                    ]

                                reaction2.annotation = annotation2
                                # add the sGPR to the annotation of the reaction

                                # extract the compartment from the reaction
                                comp = reaction2.compartments

                                # transform set to list
                                comp = list(comp)

                                # get the compartment name from the compartment id
                                compartment_name = comp_dict[comp[0]]

                                if "sGPR" in reaction2.annotation:
                                    print(
                                        "before adding sGPR: ",
                                        reaction2.annotation["sGPR"],
                                    )

                                try:
                                    reaction2.annotation["sGPR"] = variables[x2.id][0][
                                        compartment_name
                                    ]
                                except KeyError:
                                    reaction2.annotation["sGPR"] = ""
                                print(
                                    "after adding sGPR: ", reaction2.annotation["sGPR"]
                                )

                                model2.add_reactions([reaction2])

                                print(model2.reactions.get_by_id(reaction2.id))
                                print(len(model2.reactions))
                                for gene in re.split(
                                    " or | and ", reaction2.gene_reaction_rule
                                ):

                                    gene = gene.replace(")", "").replace("(", "")
                                    if not gene in variables:
                                        # Prefer cached Ensembl annotations when available
                                        ensemble = None
                                        try:
                                            if ensembl_cache and gene in ensembl_cache:
                                                val = ensembl_cache[gene]
                                                if isinstance(val, dict):
                                                    maybe = val.get("ensembl")
                                                    if maybe:
                                                        ensemble = [maybe]
                                                elif isinstance(val, str):
                                                    ensemble = [val]
                                        except Exception:
                                            ensemble = None

                                        if not ensemble:
                                            try:
                                                ensemble = sorted(
                                                    list(
                                                        set(
                                                            re.findall(
                                                                "ENS[A-Z][0-9]+",
                                                                str(
                                                                    urllib.request.urlopen(
                                                                        "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g="
                                                                        + gene
                                                                    ).read()
                                                                ),
                                                            )
                                                        )
                                                    )
                                                )
                                            except Exception:
                                                ensemble = []

                                        variables[gene] = ensemble
                                    else:
                                        ensemble = variables[gene]

                                    model2.genes.get_by_id(gene).annotation[
                                        "ensembl"
                                    ] = ensemble
                                    model2.genes.get_by_id(gene).annotation[
                                        "hgnc.symbol"
                                    ] = variables[x2.id][3][gene]

                                model2.groups.get_by_id(RxnPath[x2.id][0]).members.add(
                                    reaction2
                                )

                                RxnList.append(reaction2.reaction)
                                n = n + 1
                            else:
                                x2.gene_reaction_rule = variables[x2.id][2][CSL2]

                                # extract the compartment from the reaction
                                comp = x2.compartments

                                # transform set to list
                                comp = list(comp)

                                # get the compartment name
                                compartment_name = comp_dict[comp[0]]

                                if "sGPR" in x2.annotation:
                                    print("before adding sGPR: ", x2.annotation["sGPR"])

                                # upload the reaction in the model
                                try:
                                    model2.reactions.get_by_id(x2.id).annotation[
                                        "sGPR"
                                    ] = variables[x2.id][0][compartment_name]
                                    print(
                                        "after adding sGPR: ",
                                        variables[x2.id][0][compartment_name],
                                    )
                                except KeyError:
                                    if (
                                        "sGPR"
                                        in model2.reactions.get_by_id(x2.id).annotation
                                    ):
                                        # print("before adding sGPR: (key error) ", model2.reactions.get_by_id(x2.id).annotation['sGPR'])
                                        pass  # do nothing
                                    else:
                                        model2.reactions.get_by_id(x2.id).annotation[
                                            "sGPR"
                                        ] = ""

                                    print(
                                        "after adding sGPR: (key error) ",
                                        model2.reactions.get_by_id(x2.id).annotation[
                                            "sGPR"
                                        ],
                                    )

                                for gene in re.split(
                                    " or | and ", x2.gene_reaction_rule
                                ):
                                    gene = gene.replace(")", "").replace("(", "")
                                    if not gene in variables:
                                        # Try cache first
                                        ensemble = None
                                        try:
                                            if ensembl_cache and gene in ensembl_cache:
                                                val = ensembl_cache[gene]
                                                if isinstance(val, dict):
                                                    maybe = val.get("ensembl")
                                                    if maybe:
                                                        ensemble = [maybe]
                                                elif isinstance(val, str):
                                                    ensemble = [val]
                                        except Exception:
                                            ensemble = None

                                        if not ensemble:
                                            try:
                                                ensemble = sorted(
                                                    list(
                                                        set(
                                                            re.findall(
                                                                "ENS[A-Z][0-9]+",
                                                                str(
                                                                    urllib.request.urlopen(
                                                                        "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g="
                                                                        + gene
                                                                    ).read()
                                                                ),
                                                            )
                                                        )
                                                    )
                                                )
                                            except Exception:
                                                ensemble = []

                                        variables[gene] = ensemble
                                    else:
                                        ensemble = variables[gene]

                                    model2.genes.get_by_id(gene).annotation[
                                        "hgnc.symbol"
                                    ] = variables[x2.id][3][gene]
                                    model2.genes.get_by_id(gene).annotation[
                                        "ensembl"
                                    ] = ensemble

            except Exception as e:
                print(e)
                print(traceback.format_exc())
                print(x2.id)
                eList.append(x2.id)
                handle.write(x2.id + "\n")

            print("final", len(model.reactions))
            print("final2", len(model2.reactions))

        print()
        print(f"{len(model2.reactions)-len(rr)} new reactions")
        print([rxn.id for rxn in model2.reactions if not rxn.id in rr])

        print(f"{len(model2.metabolites)-len(m2)} new metabolites")
        print([m.id for m in model2.metabolites if not m.id in m2])

        print(f"{len(model2.genes)-len(geneList)} new genes")
        print([g.id for g in model2.genes if not g.id in geneList])

        print()
        print("An error was caused for the following reactions")

        print(eList)

        # Before saving the model, put a name to the newly added compartments

        print("Updating compartment names. Before:")
        print(model2.compartments)
        model2 = update_comp_names_bm(model2, comp_abb_file)
        print("After:")
        print(model2.compartments)

        print("Saving variable environment")
        with open(variablesFile, "wb") as savedVariables:
            for variable in [variables]:
                dill.dump(variable, savedVariables, protocol=pickle.HIGHEST_PROTOCOL)
        print("Variable environment saved in " + variablesFile)

        print("Saving model")

        write_sbml_model(model, Output)
        write_sbml_model(model2, Output3)

        print("The model has been successfully written to " + Output + "\n")

        handle.close()

    model3 = read_sbml_model(Output3)
    model3_sec_copy = copy.deepcopy(model3)

    ########################### Evaluate and Correct the model ############################
    model = read_sbml_model(input_model_path)

    model4 = copy.deepcopy(model3)

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
        print(count, "/", len(pattern))
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

    ## Check for umbalanced reactions

    list_reaction_balance = []
    for x in model4.reactions:
        print(x.id)
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

    # Writte the sbml model:

    write_sbml_model(model4, Output4)  # Human1.5

    ########################### Generate THG_beta1 from THG_beta2 ############################

    model5 = copy.deepcopy(model4)

    # Reactions

    list_of_reactions_in_thgb2_not_in_thgb1 = [
        x.id for x in model4.reactions if not x.id in [y.id for y in model.reactions]
    ]

    for x in list_of_reactions_in_thgb2_not_in_thgb1:
        for g in model5.groups:  # remove removed reactions from groups
            if x in [m.id for m in g.members]:
                g.members.remove(model5.reactions.get_by_id(x))
        model5.reactions.remove(x)

    # Metabolites

    list_of_metabolites_in_thgb2_not_in_thgb1 = [
        x.id
        for x in model4.metabolites
        if not x.id in [y.id for y in model.metabolites]
    ]

    for x in list_of_metabolites_in_thgb2_not_in_thgb1:
        model5.metabolites.remove(x)

    ## Check for umbalanced reactions

    list_reaction_balance = []
    for x in model5.reactions:
        print(x.id)
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
                mass_balance_test = test_reaction_balance(eq_test)
                if not mass_balance_test[0]:
                    MB = mass_balance(eq, "R")

                    if x.id in [r.id for r in model.reactions]:
                        original_model_reaction = model.reactions.get_by_id(
                            x.id
                        ).reaction
                    else:
                        original_model_reaction = ""

                    list_reaction_balance.append(
                        [x.id, mass_balance_test, MB, original_model_reaction]
                    )

    # Write the sbml model:

    write_sbml_model(model5, Output5)  # Human1.2

    # Test if sbml models can be read

    test1 = read_sbml_model(Output4)
    test2 = read_sbml_model(Output5)
