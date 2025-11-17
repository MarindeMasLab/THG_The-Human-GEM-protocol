import cobra
from cobra.io import read_sbml_model, write_sbml_model
import csv
import pandas as pd


def load_candidates(csv_path):
    """Load candidate reactions from simple text file (one reaction ID per line)."""
    with open(csv_path, "r") as f:
        return set(line.strip() for line in f if line.strip())


def load_minimal_candidates(csv_path):
    """Load minimal candidate reactions from CSV file.

    Expects format: candidate_reaction,active_in_any_solution
    Returns set of reaction IDs that are active (active_in_any_solution == '1')
    """
    candidates = set()
    with open(csv_path, "r") as f:
        reader = csv.reader(f)
        next(reader, None)  # Skip header
        for row in reader:
            if len(row) >= 2 and row[1] == "1":
                candidates.add(row[0])
    return candidates


# Compartment abbreviation dictionary (reuse from earlier)
dict_compartment_abbr = {
    "a": "Cell membrane",
    "c": "Cytosol",
    "r": "Endoplasmic reticulum",
    "e": "Extracellular",
    "y": "Glycocalix",
    "g": "Golgi apparatus",
    "i": "Inner mitochondria",
    "l": "Lysosome",
    "m": "Mitochondria",
    "n": "Nucleus",
    "x": "Peroxisome",
    "v": "Vesicle",
    "ca": "Cell membrane Apical",
    "cb": "Cell membrane Basal",
    "cj": "Cell membrane Junktion/Adhesion",
    "ci": "Cilium",
    "ck": "Cytoskeleton",
}


def compartment_counts(model, rxn_ids=None):
    counts = {}
    if rxn_ids is None:
        rxns = model.reactions
    else:
        rxns = [model.reactions.get_by_id(rid) for rid in rxn_ids]
    for rxn in rxns:
        temp_compartments = rxn.compartments
        if len(temp_compartments) == 1:
            compartment = list(temp_compartments)[0]
        else:
            compartment = "transport"
        counts[compartment] = counts.get(compartment, 0) + 1
    return counts


if __name__ == "__main__":
    # --- Define paths here ---
    model_name = "endoA_250911_clean"
    model_path = f"models/{model_name}.xml"
    candidates_path = f"scripts/created_transport_reactions_{model_name}.csv"
    minimal_candidates_path = f"minimal_transporters_pfba_{model_name}.csv"

    print(f"Loading model from {model_path}")
    model = read_sbml_model(model_path)

    # print("Determining blocked reactions")
    # blocked_rxns_start = cobra.flux_analysis.find_blocked_reactions(model, processes=60)
    # print(f"Blocked reactions at the start {len(blocked_rxns_start)}")

    print(f"Loading minimal candidate reactions from {minimal_candidates_path}")
    minimal_candidates = load_minimal_candidates(minimal_candidates_path)

    print(f"Removing non-minimal candidate reactions...")
    all_candidates = load_candidates(candidates_path)
    non_minimal = all_candidates - minimal_candidates

    print(f"Total candidate reactions: {len(all_candidates)}")
    print(f"Active (minimal) reactions: {len(minimal_candidates)}")
    print(f"Non-active reactions: {len(non_minimal)}")

    # Write minimal candidates to CSV
    minimal_out_path = f"minimal_candidates_{model_name}.csv"
    with open(minimal_out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["reaction_id"])
        for rid in sorted(minimal_candidates):
            writer.writerow([rid])
    print(f"Wrote minimal candidate reactions to {minimal_out_path}")

    model.remove_reactions(
        [
            model.reactions.get_by_id(rid)
            for rid in non_minimal
            if rid in model.reactions
        ]
    )

    write_sbml_model(model, "models/endoA_250915.xml")

    print(f"Running find_blocked_reactions...")
    blocked_rxns = cobra.flux_analysis.find_blocked_reactions(model, processes=60)
    print(
        f"Blocked reactions after removing non-minimal candidates: {len(blocked_rxns)}"
    )

    rxns = [rxn.id for rxn in model.reactions]
    comp_counts = compartment_counts(model, rxns)
    comp_counts_blocked = compartment_counts(model, blocked_rxns)

    rows = []
    for comp in sorted(comp_counts):
        row = {
            "Abbreviation": comp,
            "Full Name": dict_compartment_abbr.get(comp, comp),
            "Total in Model": comp_counts.get(comp, 0),
            "Blocked": comp_counts_blocked.get(comp, 0),
        }

        rows.append(row)
    df_compartments_huge = pd.DataFrame(rows)
    df_compartments_huge = df_compartments_huge.sort_values(
        by="Blocked", ascending=False
    ).reset_index(drop=True)
    print(df_compartments_huge)
    # for rxn in blocked_rxns:
    #     print(rxn)
