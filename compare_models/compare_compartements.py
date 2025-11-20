"""
Compare reaction counts per compartment between two cobra models.

Produces a summary printed to stdout and CSV files under "reports/".
Two comparisons are produced:
 - Raw models
 - Models without blocked reactions (computed using cobra.flux_analysis.find_blocked_reactions)

Metrics per compartment: number of reactions in model A, number in model B, number in both, and number matching exactly by stoichiometry (and ID)
"""

import csv
import os
import sys
from collections import defaultdict
from typing import Dict, Set, Tuple

from cobra.io import read_sbml_model, load_json_model
from cobra.flux_analysis import find_blocked_reactions


def reactions_by_compartment(model) -> Dict[str, Set[str]]:
    """Return mapping compartment -> set of reaction IDs that include that compartment (any metabolite in that compartment).
    Reaction counts will be based on reaction.id membership.
    """
    comps = defaultdict(set)
    for rxn in model.reactions:
        rxn_comps = set([m.compartment for m in rxn.metabolites])
        if len(rxn_comps) == 0:
            rxn_comps = {""}
        for c in rxn_comps:
            comps[c].add(rxn.id)
    return comps


def normalize_stoichiometry(rxn):
    """Return a frozenset of (met_id, stoich) pairs for a reaction useful for stoichiometry equality checks."""
    return frozenset(((m.id, float(coeff)) for m, coeff in rxn.metabolites.items()))


def comp_compare(model_a, model_b) -> dict:
    """Compute comparison metrics per compartment between model_a and model_b.

    Returns a dict where keys are compartments and values are a dict with counts:
      - n_a: number of reactions in model_a for compartment
      - n_b: number of reactions in model_b for compartment
      - intersect_id: number of reactions with identical reaction IDs present in both models in that compartment
      - intersect_stoich: number of reactions with identical ID and identical stoichiometry
    """
    comps_a = reactions_by_compartment(model_a)
    comps_b = reactions_by_compartment(model_b)

    all_comps = set(comps_a.keys()).union(comps_b.keys())

    # Precompute reaction stoichiometries for fast lookup
    stoich_a = {r.id: normalize_stoichiometry(r) for r in model_a.reactions}
    stoich_b = {r.id: normalize_stoichiometry(r) for r in model_b.reactions}

    out = {}
    for c in sorted(all_comps):
        set_a = comps_a.get(c, set())
        set_b = comps_b.get(c, set())
        intersect_id = set_a.intersection(set_b)
        intersect_stoich = set()
        for rid in intersect_id:
            if rid in stoich_a and rid in stoich_b and stoich_a[rid] == stoich_b[rid]:
                intersect_stoich.add(rid)

        out[c] = {
            "n_a": len(set_a),
            "n_b": len(set_b),
            "intersect_id": len(intersect_id),
            "intersect_stoich": len(intersect_stoich),
            "rxns_a": set_a,
            "rxns_b": set_b,
            "match_ids": sorted(list(intersect_id)),
            "match_stoich_ids": sorted(list(intersect_stoich)),
        }
    # also provide some summary
    out["_summary"] = {
        "total_rxns_a": len(model_a.reactions),
        "total_rxns_b": len(model_b.reactions),
    }
    return out


def remove_blocked_reactions(model):
    """Return a copy of the input model with universally blocked reactions removed.
    Uses find_blocked_reactions(model, open_exchanges=True)
    """
    blocked = set(find_blocked_reactions(model, open_exchanges=True))
    # Make a copy of the model
    from cobra import Model

    new_model = model.copy()
    for rid in sorted(blocked):
        if rid in new_model.reactions:
            new_model.reactions.remove(new_model.reactions.get_by_id(rid))
    return new_model


def print_comparison(comp_dict, title="Comparison"):
    """Print a nicely formatted table for human-readable output.

    The table shows per-compartment counts and two percentage columns:
      - overlap_a_pct: fraction of model A reactions in that compartment that are present in model B
      - overlap_b_pct: fraction of model B reactions in that compartment that are present in model A
    """
    print("\n" + "=" * 120)
    print(f"  {title}")
    print("=" * 120)

    headers = [
        "Compartment",
        "Model A",
        "Model B",
        "Only A",
        "Only B",
        "Both (ID)",
        "Both (Stoich)",
        "A Overlap %",
        "B Overlap %",
    ]

    # Map headers to dict keys
    header_keys = [
        "compartment",
        "n_model_a",
        "n_model_b",
        "only_in_a",
        "only_in_b",
        "in_both_by_id",
        "in_both_by_stoich",
        "overlap_a_pct",
        "overlap_b_pct",
    ]

    # Alignment: left for compartment, right for all numbers/percentages
    alignments = ["left"] + ["right"] * 8

    # Compute column widths and build rows
    rows = []
    max_widths = {h: len(h) for h in headers}
    total_a = comp_dict.get("_summary", {}).get("total_rxns_a", 0)
    total_b = comp_dict.get("_summary", {}).get("total_rxns_b", 0)

    for c, d in comp_dict.items():
        if c == "_summary":
            continue
        only_a = d["n_a"] - d["intersect_id"]
        only_b = d["n_b"] - d["intersect_id"]
        overlap_a_pct = 100 * d["intersect_id"] / d["n_a"] if d["n_a"] else 0
        overlap_b_pct = 100 * d["intersect_id"] / d["n_b"] if d["n_b"] else 0

        row_data = {
            "compartment": c if c != "" else "UNASSIGNED",
            "n_model_a": str(d["n_a"]),
            "n_model_b": str(d["n_b"]),
            "only_in_a": str(only_a),
            "only_in_b": str(only_b),
            "in_both_by_id": str(d["intersect_id"]),
            "in_both_by_stoich": str(d["intersect_stoich"]),
            "overlap_a_pct": f"{overlap_a_pct:5.1f}%",
            "overlap_b_pct": f"{overlap_b_pct:5.1f}%",
        }
        rows.append(row_data)

        # Update widths
        for i, h in enumerate(headers):
            max_widths[h] = max(max_widths[h], len(row_data[header_keys[i]]))

    # Print header with proper alignment
    header_parts = []
    for i, h in enumerate(headers):
        if alignments[i] == "left":
            header_parts.append(h.ljust(max_widths[h]))
        else:
            header_parts.append(h.rjust(max_widths[h]))

    print(" | ".join(header_parts))
    print("-+-".join(["-" * max_widths[h] for h in headers]))

    # Print rows with proper alignment
    for row in rows:
        row_parts = []
        for i, h in enumerate(headers):
            val = row[header_keys[i]]
            if alignments[i] == "left":
                row_parts.append(val.ljust(max_widths[h]))
            else:
                row_parts.append(val.rjust(max_widths[h]))
        print(" | ".join(row_parts))

    # Print summary
    print("=" * 120)
    print(
        f"SUMMARY: Model A = {total_a:,} rxns  |  Model B = {total_b:,} rxns  |  Compartments = {len(rows)}"
    )
    print("=" * 120)


def save_csv(comp_dict, outpath):
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "compartment",
                "n_model_a",
                "n_model_b",
                "only_in_a",
                "only_in_b",
                "n_in_both_by_id",
                "n_in_both_by_stoich",
                "overlap_a_pct",
                "overlap_b_pct",
                "match_ids",
                "match_stoich_ids",
            ]
        )
        for c, d in comp_dict.items():
            if c == "_summary":
                continue
            only_a = d["n_a"] - d["intersect_id"]
            only_b = d["n_b"] - d["intersect_id"]
            overlap_a_pct = 100 * d["intersect_id"] / d["n_a"] if d["n_a"] else 0
            overlap_b_pct = 100 * d["intersect_id"] / d["n_b"] if d["n_b"] else 0
            w.writerow(
                [
                    c,
                    d["n_a"],
                    d["n_b"],
                    only_a,
                    only_b,
                    d["intersect_id"],
                    d["intersect_stoich"],
                    f"{overlap_a_pct:.1f}%",
                    f"{overlap_b_pct:.1f}%",
                    ";".join(d.get("match_ids", [])),
                    ";".join(d.get("match_stoich_ids", [])),
                ]
            )
        w.writerow([])
        w.writerow(
            [
                "SUMMARY",
                comp_dict["_summary"]["total_rxns_a"],
                comp_dict["_summary"]["total_rxns_b"],
            ]
        )


def main(model_name_1, model_name_2, outdir=None):
    path_a = os.path.abspath(f"{model_name_1}.json")
    path_b = os.path.abspath(f"{model_name_2}.json")
    outdir = (
        os.path.abspath(outdir)
        if outdir
        else os.path.join(os.path.dirname(__file__), "reports")
    )
    os.makedirs(outdir, exist_ok=True)

    print(f"Loading models: {path_a} and {path_b}")
    model_a = load_json_model(path_a)
    model_b = load_json_model(path_b)

    print("Performing raw comparison...")
    raw = comp_compare(model_a, model_b)
    print_comparison(raw, title="Raw models")
    save_csv(raw, os.path.join(outdir, "compartments_comparison_raw.csv"))

    print("Removing blocked reactions from both models and comparing again...")
    model_a_nb = remove_blocked_reactions(model_a)
    model_b_nb = remove_blocked_reactions(model_b)
    nb = comp_compare(model_a_nb, model_b_nb)
    print_comparison(nb, title="Models without blocked reactions")
    save_csv(nb, os.path.join(outdir, "compartments_comparison_no_blocked.csv"))

    print(f"Reports written to {outdir}")


if __name__ == "__main__":

    # model_name_1 = "models/THG-beta-expanded_251118"
    model_name_1 = "models/THG-beta-expanded_251118_transcriptomics"
    model_name_2 = "models/endoA_251119"
    output_dir = "compare_models/reports"

    main(model_name_1, model_name_2, output_dir)
