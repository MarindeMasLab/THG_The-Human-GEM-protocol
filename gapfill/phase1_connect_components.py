"""Phase 1 candidate generator (renamed)

This is a copy of the previous `connect_components.py` with the same behavior
but now named `phase1_connect_components.py` to be called from `gapfill.py`.
"""

from cobra.io import load_json_model
from cobra.manipulation.delete import prune_unused_metabolites
import networkx as nx
from collections import defaultdict
import re
import os
import csv


def generate_phase1_outputs(model, allowed_connections=None, conservative=False, out_dir=None):
    """Generate Phase-1 candidate list and dead-end summaries.

    See original script for details.
    """
    base_dir = os.path.dirname(__file__)
    if out_dir is None:
        out_dir = os.path.join(base_dir, "files")
    os.makedirs(out_dir, exist_ok=True)

    model, _ = prune_unused_metabolites(model)

    G = nx.DiGraph()
    for met in model.metabolites:
        G.add_node(met.id, bipartite=0)
    for rxn in model.reactions:
        G.add_node(rxn.id, bipartite=1)
        for met, coeff in rxn.metabolites.items():
            if coeff < 0:
                G.add_edge(met.id, rxn.id)
            elif coeff > 0:
                G.add_edge(rxn.id, met.id)

    components = list(nx.weakly_connected_components(G))
    components = sorted(components, key=len, reverse=True)

    met_ids = [m.id for m in model.metabolites]
    met_ids_set = set(met_ids)
    met_to_comp = {}
    for idx, comp in enumerate(components):
        for node in comp:
            if node in met_ids_set:
                met_to_comp[node] = idx + 1

    met_produced = {m.id: False for m in model.metabolites}
    met_consumed = {m.id: False for m in model.metabolites}
    for rxn in model.reactions:
        for m, coeff in rxn.metabolites.items():
            if coeff > 0:
                met_produced[m.id] = True
            if coeff < 0:
                met_consumed[m.id] = True

    met_type = {}
    for m in model.metabolites:
        produced = met_produced.get(m.id, False)
        consumed = met_consumed.get(m.id, False)
        if produced and not consumed:
            met_type[m.id] = "product"
        elif consumed and not produced:
            met_type[m.id] = "substrate"
        else:
            met_type[m.id] = "none"

    base_ids = [re.sub(r"[a-z]+$", "", mid) for mid in met_ids]
    buckets = defaultdict(list)
    for mid, base in zip(met_ids, base_ids):
        buckets[base].append((mid, base))

    if conservative:
        transport_bases = set()
        for rxn in model.reactions:
            mets = list(rxn.metabolites.items())
            if len(mets) == 2:
                (mA, cA), (mB, cB) = mets
                if (cA, cB) == (-1, 1) or (cA, cB) == (1, -1):
                    transport_bases.add(re.sub(r"[a-z]+$", "", mA.id))
                    transport_bases.add(re.sub(r"[a-z]+$", "", mB.id))
        buckets = {b: items for b, items in buckets.items() if b in transport_bases}

    def comp_suffix(mid):
        m = re.search(r"[a-z]+$", mid)
        return m.group(0) if m else ""

    default_allowed = [
        ("a", "ca"), ("a", "cb"), ("a", "cj"), ("a", "ci"), ("a", "ck"),
        ("a", "c"), ("a", "r"), ("a", "e"), ("a", "v"), ("ca", "a"),
        ("ca", "ci"), ("ca", "ck"), ("ca", "c"), ("ca", "e"), ("ca", "v"),
        ("cb", "a"), ("cb", "ck"), ("cb", "e"), ("cb", "v"), ("cj", "a"),
        ("cj", "ck"), ("cj", "e"), ("cj", "v"), ("ci", "a"), ("ci", "ca"),
        ("ci", "ck"), ("ci", "c"), ("ci", "e"), ("ck", "a"), ("ck", "ca"),
        ("ck", "cb"), ("ck", "cj"), ("ck", "ci"), ("c", "a"), ("c", "ca"),
        ("c", "cb"), ("c", "cj"), ("c", "ci"), ("c", "ck"), ("c", "r"),
        ("c", "g"), ("c", "l"), ("c", "m"), ("c", "n"), ("c", "x"), ("c", "v"),
        ("r", "c"), ("e", "a"), ("e", "ca"), ("e", "cb"), ("e", "cj"), ("e", "ci"),
        ("e", "l"), ("y", "a"), ("y", "ca"), ("y", "e"), ("g", "c"), ("i", "g"),
        ("l", "c"), ("m", "c"), ("m", "n"), ("n", "c"), ("x", "c"), ("v", "a"),
        ("v", "ca"), ("v", "cb"), ("v", "cj"), ("v", "r"),
    ]

    if allowed_connections:
        allowed_pairs = set(tuple(x) for x in allowed_connections)
    else:
        allowed_pairs = set(default_allowed)
    allowed_both = set()
    for a, b in allowed_pairs:
        allowed_both.add((a, b))
        allowed_both.add((b, a))

    matches = []
    for base, items in buckets.items():
        if len(items) < 2:
            continue
        n = len(items)
        for i in range(n):
            for j in range(i + 1, n):
                met1, _ = items[i]
                met2, _ = items[j]
                comp1 = met_to_comp.get(met1)
                comp2 = met_to_comp.get(met2)
                if comp1 is None or comp2 is None or comp1 == comp2:
                    continue
                s1 = comp_suffix(met1)
                s2 = comp_suffix(met2)
                if allowed_both and (s1, s2) not in allowed_both:
                    continue
                t1 = met_type.get(met1, "none")
                t2 = met_type.get(met2, "none")
                if (t1 == "substrate" and t2 == "product") or (t1 == "product" and t2 == "substrate"):
                    cand_type = "A"
                elif t1 == t2 and t1 in ("substrate", "product"):
                    cand_type = "B"
                else:
                    cand_type = "C"

                matches.append({
                    "met1": met1,
                    "met2": met2,
                    "base": base,
                    "suffix1": s1,
                    "suffix2": s2,
                    "component1": comp1,
                    "component2": comp2,
                    "type": cand_type,
                    "met1_type": t1,
                    "met2_type": t2,
                })

    cand_path = os.path.join(out_dir, "candidates_all.csv")
    with open(cand_path, "w", newline="") as cf:
        fieldnames = [
            "met1",
            "met2",
            "base",
            "suffix1",
            "suffix2",
            "component1",
            "component2",
            "type",
            "met1_type",
            "met2_type",
        ]
        writer = csv.DictWriter(cf, fieldnames=fieldnames)
        writer.writeheader()
        for m in matches:
            writer.writerow({k: m.get(k, "") for k in fieldnames})

    dead_path = os.path.join(out_dir, "deadends_summary.csv")
    with open(dead_path, "w", newline="") as df:
        fieldnames = ["met_id", "type", "component", "compartment"]
        writer = csv.DictWriter(df, fieldnames=fieldnames)
        writer.writeheader()
        for met in model.metabolites:
            t = met_type.get(met.id, "none")
            if t in ("product", "substrate"):
                comp = met_to_comp.get(met.id, "")
                writer.writerow({"met_id": met.id, "type": t, "component": comp, "compartment": met.compartment})

    comp_path = os.path.join(out_dir, "components_summary.csv")
    with open(comp_path, "w", newline="") as pf:
        fieldnames = ["component_id", "num_nodes"]
        writer = csv.DictWriter(pf, fieldnames=fieldnames)
        writer.writeheader()
        for idx, comp in enumerate(components):
            writer.writerow({"component_id": idx + 1, "num_nodes": len(comp)})

    return cand_path, dead_path, comp_path


if __name__ == "__main__":
    default_model = os.path.join(os.path.dirname(__file__), "..", "models", "base", "THG-beta-batch_251106_unconnected.json")
    default_model = os.path.normpath(default_model)
    if not os.path.exists(default_model):
        print(f"Default model not found at {default_model}. Please provide a model path.")
    else:
        print(f"Loading model {default_model}")
        m = load_json_model(default_model)
        cand, dead, comp = generate_phase1_outputs(m)
        print(f"Wrote candidates: {cand}\ndeadends: {dead}\ncomponents: {comp}")
