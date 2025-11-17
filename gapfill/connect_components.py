from cobra.io import read_sbml_model, load_json_model
from cobra.util import create_stoichiometric_matrix
from cobra.manipulation.delete import prune_unused_reactions, prune_unused_metabolites
import networkx as nx
from tqdm import tqdm

import logging


def get_unconnected_components(model):

    model, _ = prune_unused_metabolites(model)

    logging.info(f"Number of metabolites: {len(model.metabolites)}")
    logging.info(f"Number of reactions: {len(model.reactions)}")

    # (S is not used but keep a small sanity check if needed)
    # S = create_stoichiometric_matrix(model)
    # Build a directed bipartite graph: metabolites <-> reactions
    G = nx.DiGraph()

    # Add metabolite nodes
    for met in model.metabolites:
        G.add_node(met.id, bipartite=0)

    # Add reaction nodes
    for rxn in model.reactions:
        G.add_node(rxn.id, bipartite=1)

    # Add edges: metabolite -> reaction (if consumed), reaction -> metabolite (if produced)
    for rxn in model.reactions:
        for met, coeff in rxn.metabolites.items():
            if coeff < 0:
                G.add_edge(met.id, rxn.id)  # consumed
            elif coeff > 0:
                G.add_edge(rxn.id, met.id)  # produced

    # Precompute compartment mappings for metabolites and reactions
    met_compartments = {met.id: met.compartment for met in model.metabolites}
    # Reaction compartments aren't stored directly on Reaction objects in a
    # consistent way; infer them from the compartments of participating metabolites
    rxn_compartments = {
        rxn.id: set(m.compartment for m in rxn.metabolites) for rxn in model.reactions
    }

    # Get all weakly connected components
    components = list(nx.weakly_connected_components(G))
    components = sorted(components, key=len, reverse=True)  # Order by size descending
    logging.info(f"Number of weakly connected components: {len(components)}")

    # Create a set of reaction IDs for fast lookup
    reaction_ids = set(rxn_compartments.keys())

    # iterate components with a progress bar
    for idx, comp in enumerate(tqdm(components, desc="components")):
        # Count reactions in this component by comparing node IDs
        reaction_count = sum(1 for node in comp if node in reaction_ids)

        # Efficiently collect compartments present in this component
        met_ids = comp & set(met_compartments.keys())
        rxn_ids = comp & reaction_ids
        compartments = set()
        compartments.update(met_compartments[mid] for mid in met_ids)
        for rid in rxn_ids:
            compartments.update(rxn_compartments[rid])

        logging.info(
            f"Component {idx+1}: {len(comp)} nodes, {reaction_count} reactions, compartments: {sorted(compartments)}"
        )

    # Get the largest weakly connected component (if any)
    if len(components) == 0:
        logging.info("No components found.")
        return components

    largest_component = components[0]  # Already sorted
    logging.info(f"Largest weakly connected component size: {len(largest_component)}")
    logging.info(f"Total nodes in network: {G.number_of_nodes()}")

    if len(largest_component) == G.number_of_nodes():
        logging.info("The network is fully connected.")
    else:
        logging.info("The network is not fully connected.")

    return components


def connect_components(model, allowed_connections=None, conservative=False):
    """
    Connect disconnected components in a metabolic model by adding transport reactions.

    Parameters:
    - model: COBRA model to connect
    - allowed_connections: List of allowed compartment pairs for transport (optional)
    - conservative: If True, only create transport reactions for metabolites that already
      participate in transport reactions. If False (default), create transport reactions
      for all metabolites that appear in multiple compartments across disconnected components.
    """

    # If you want to use the manual option, uncomment and edit the following line:
    allowed_compartment_pairs = [
        ("a", "ca"),
        ("a", "cb"),
        ("a", "cj"),
        ("a", "ci"),
        ("a", "ck"),
        ("a", "c"),
        ("a", "r"),
        ("a", "e"),
        ("a", "v"),
        ("ca", "a"),
        ("ca", "ci"),
        ("ca", "ck"),
        ("ca", "c"),
        ("ca", "e"),
        ("ca", "v"),
        ("cb", "a"),
        ("cb", "ck"),
        ("cb", "e"),
        ("cb", "v"),
        ("cj", "a"),
        ("cj", "ck"),
        ("cj", "e"),
        ("cj", "v"),
        ("ci", "a"),
        ("ci", "ca"),
        ("ci", "ck"),
        ("ci", "c"),
        ("ci", "e"),
        ("ck", "a"),
        ("ck", "ca"),
        ("ck", "cb"),
        ("ck", "cj"),
        ("ck", "ci"),
        ("c", "a"),
        ("c", "ca"),
        ("c", "cb"),
        ("c", "cj"),
        ("c", "ci"),
        ("c", "ck"),
        ("c", "r"),
        ("c", "g"),
        ("c", "l"),
        ("c", "m"),
        ("c", "n"),
        ("c", "x"),
        ("c", "v"),
        ("r", "c"),
        ("e", "a"),
        ("e", "ca"),
        ("e", "cb"),
        ("e", "cj"),
        ("e", "ci"),
        ("e", "l"),
        ("y", "a"),
        ("y", "ca"),
        ("y", "e"),
        ("g", "c"),
        ("i", "g"),
        ("l", "c"),
        ("m", "c"),
        ("m", "n"),
        ("n", "c"),
        ("x", "c"),
        ("v", "a"),
        ("v", "ca"),
        ("v", "cb"),
        ("v", "cj"),
        ("v", "r"),
    ]

    components = get_unconnected_components(model)
    # Optimized matching: group by prefix after stripping lowercase letters (near O(n))
    import re
    import json
    import csv
    from collections import defaultdict
    import os

    # Build list of metabolite ids
    met_ids = [met.id for met in model.metabolites]

    # Compute canonical base IDs by removing lowercase letters, and the prefix (base without final char)
    # Using re.sub in Python is fast; this avoids the quadratic pairwise loop below
    base_ids = [re.sub(r"[a-z]+$", "", mid) for mid in met_ids]
    prefixes = [b for b in base_ids]

    # Group metabolites by prefix. Only items sharing the same prefix can match
    buckets = defaultdict(list)
    for mid, pref, base in zip(met_ids, prefixes, base_ids):
        buckets[pref].append((mid, base))

    # Map metabolite id -> component number (component numbering matches earlier prints: idx+1)
    met_to_comp = {}
    met_ids_set = set(met_ids)
    for idx, comp in enumerate(components):
        # Only consider metabolite nodes inside each component
        for node in comp:
            if node in met_ids_set:
                met_to_comp[node] = idx + 1

    # --- Filter metabolites: only keep those that participate in any transport reaction (between any compartments) ---
    if conservative:
        # Build set of base metabolite IDs involved in transport reactions (remove compartment suffix)
        met_with_transport_base = set()
        for rxn in model.reactions:
            mets = rxn.metabolites
            if len(mets) == 2:
                items = list(mets.items())
                (metA, coeffA), (metB, coeffB) = items[0], items[1]
                # Check for transport stoichiometry (-1, +1) or (+1, -1)
                if (coeffA, coeffB) == (-1, 1) or (coeffA, coeffB) == (1, -1):
                    # Remove compartment suffix using regex
                    metA_base = re.sub(r"[a-z]+$", "", metA.id)
                    metB_base = re.sub(r"[a-z]+$", "", metB.id)
                    met_with_transport_base.add(metA_base)
                    met_with_transport_base.add(metB_base)

        # Filter met_ids to only those whose base ID is in met_with_transport_base
        met_ids = [
            mid
            for mid in met_ids
            if re.sub(r"[a-z]+$", "", mid) in met_with_transport_base
        ]

        # Recompute base_ids and prefixes for filtered metabolites
        base_ids = [re.sub(r"[a-z]+$", "", mid) for mid in met_ids]
        prefixes = [b for b in base_ids]

        # Group metabolites by prefix again
        buckets = defaultdict(list)
        for mid, pref, base in zip(met_ids, prefixes, base_ids):
            buckets[pref].append((mid, base))
    else:
        # Non-conservative mode: use all metabolites
        pass

    # Build matches_list from buckets but keep only inter-component pairs
    matches_list = []
    for pref, items in buckets.items():
        if len(items) < 2:
            continue
        n = len(items)
        for i in range(n):
            for j in range(i + 1, n):
                met1, base1 = items[i]
                met2, base2 = items[j]
                # verify base equality (defensive check)
                if base1 != base2:
                    continue
                comp1 = met_to_comp.get(met1)
                comp2 = met_to_comp.get(met2)
                # ignore pairs where either metabolite isn't assigned to a component
                if comp1 is None or comp2 is None:
                    continue
                # only care about metabolites in different components
                if comp1 == comp2:
                    continue
                matches_list.append(
                    {
                        "met1": met1,
                        "met2": met2,
                        "base_id1": base1,
                        "base_id2": base2,
                        "component1": comp1,
                        "component2": comp2,
                    }
                )

    print(
        f"Found {len(matches_list)} inter-component matches (grouped by prefix after removing compartment suffix{', filtered for metabolites with any transport reaction by base ID' if conservative else ''})."
    )

    # --- Filter matches by compartment suffix (trailing lowercase letters in metabolite id) ---
    # Compartment for a metabolite is defined as the trailing lowercase letters in its id.
    # The variable `allowed_compartment_pairs` should be an iterable of 2-tuples like [('a','ca'), ('a','cb'), ...]
    # If `allowed_compartment_pairs` is not defined in the notebook, the cell will try to load a CSV at
    # 'data/compartment_pairs.csv' (two columns: comp1,comp2). If neither is provided, no filtering is applied
    # and a short instruction is printed.

    # helper to extract trailing lowercase suffix
    def comp_suffix(met_id):
        m = re.search(r"[a-z]+$", met_id)
        return m.group(0) if m else ""

    # build suffix map for all metabolites appearing in matches
    suffix_map = {}
    for m in met_ids:
        suffix_map[m] = comp_suffix(m)

    # obtain allowed pairs: prefer explicit argument, otherwise use the function-local default
    try:
        if allowed_connections:
            _allowed_pairs = set(tuple(x) for x in allowed_connections)
        else:
            _allowed_pairs = set(tuple(x) for x in allowed_compartment_pairs)
    except Exception as e:
        logging.warning(
            f"Could not parse allowed_connections; falling back to default allowed_compartment_pairs. Error: {e}"
        )
        try:
            _allowed_pairs = set(tuple(x) for x in allowed_compartment_pairs)
        except Exception:
            _allowed_pairs = set()

    # If no allowed pairs are available, skip compartment-pair filtering (keep matches_list as-is)
    if not _allowed_pairs:
        print(
            "No allowed compartment pairs provided; skipping compartment-pair filtering."
        )
        # keep matches_list as-is
    else:
        # normalize allowed pairs to include both orders (comp1,comp2) and (comp2,comp1)
        allowed_both = set()
        for c1, c2 in _allowed_pairs:
            allowed_both.add((c1, c2))
            allowed_both.add((c2, c1))

        filtered = []
        for m in matches_list:
            s1 = suffix_map.get(m["met1"], "")
            s2 = suffix_map.get(m["met2"], "")
            if (s1, s2) in allowed_both:
                # include suffix info for clarity
                m["suffix1"] = s1
                m["suffix2"] = s2
                filtered.append(m)
        matches_list = filtered
        print(f"After compartment-pair filtering: {len(matches_list)} matches remain.")

    # Create transport reactions for each match in matches_list and add them to the model
    # Reaction stoichiometry: met1 -> met2 (consumes met1 in its compartment, produces met2 in its compartment)
    # Reactions are named TRANS_<base>_<suf1>_<suf2>_<i> and made reversible by default.
    from cobra import Reaction
    import csv

    created_reactions = []
    failed = []

    for i, m in enumerate(matches_list):
        met1 = m["met1"]
        met2 = m["met2"]
        base = m.get("base_id1", m.get("base_id2", ""))
        s1 = m.get("suffix1", "")
        s2 = m.get("suffix2", "")

        # build a safe reaction id
        raw_id = f"TRANS_{base}_{s1}_{s2}_{i}"
        rid = re.sub(r"[^0-9A-Za-z_]", "_", raw_id)

        # skip if reaction already exists (check by id)
        try:
            _ = model.reactions.get_by_id(rid)
            # reaction id already present
            continue
        except KeyError:
            pass

        try:
            met_obj1 = model.metabolites.get_by_id(met1)
            met_obj2 = model.metabolites.get_by_id(met2)
        except KeyError as e:
            failed.append((rid, met1, met2, f"missing_metabolite: {e}"))
            continue

        # create reaction
        rxn = Reaction(rid)
        rxn.name = f"transport {met1} -> {met2}"
        # stoichiometry: met1 consumed, met2 produced
        rxn.add_metabolites({met_obj1: -1.0, met_obj2: 1.0})
        # make reversible (allow transport both ways)
        rxn.lower_bound = -1000.0
        rxn.upper_bound = 1000.0
        # annotate
        rxn.annotation = {
            "created_from_match": "true",
            "base_id": base,
            "suffix1": s1,
            "suffix2": s2,
            "component1": m.get("component1"),
            "component2": m.get("component2"),
        }

        try:
            model.add_reactions([rxn])
            created_reactions.append(
                {
                    "reaction_id": rid,
                    "met1": met1,
                    "met2": met2,
                    "base": base,
                    "suffix1": s1,
                    "suffix2": s2,
                    "component1": m.get("component1"),
                    "component2": m.get("component2"),
                }
            )
        except Exception as e:
            failed.append((rid, met1, met2, str(e)))

    print(
        f"Created {len(created_reactions)} transport reactions; {len(failed)} failures."
    )

    return model, created_reactions


if __name__ == "__main__":

    from cobra.io import read_sbml_model, write_sbml_model
    import csv

    model_name = "THG-beta-batch_251106"

    model = read_sbml_model(f"models/{model_name}.xml")

    model, created_reactions = connect_components(model)

    write_sbml_model(model, f"models/{model_name}_connected.xml")

    # Save created reactions to CSV
    if created_reactions:
        csv_filename = f"models/{model_name}_created_reactions.csv"
        with open(csv_filename, "w", newline="") as csvfile:
            fieldnames = created_reactions[0].keys()
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(created_reactions)
        print(f"Saved {len(created_reactions)} created reactions to {csv_filename}")
    else:
        print("No reactions were created.")
