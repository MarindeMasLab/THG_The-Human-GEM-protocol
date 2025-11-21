from compaction import full_compaction
import logging


def remove_loops(model, no_blocked_reactions=False):
    """Remove loops from the model by performing full compaction."""
    _, loop_reactions = full_compaction(model, no_blocked_reactions=False)

    logging.info(f"Identified {len(loop_reactions)} infeasible loops to remove.")

    loop_individual_reactions = []

    for compacted_reaction in loop_reactions:
        # Remove all brackets "()"
        cleaned_reaction = compacted_reaction.id.replace("(", "").replace(")", "")
        # Split by # and @
        parts = cleaned_reaction.split("#")
        for part in parts:
            subparts = part.split("@")
            loop_individual_reactions.extend(subparts)
    logging.info(f"Compacted loops identified: {[reaction.id for reaction in loop_reactions]}")
    logging.info(f"Individual reactions to remove: {loop_individual_reactions}")
    logging.info(
        f"Removing {len(loop_individual_reactions)} reactions involved in infeasible loops."
    )

    model.remove_reactions(loop_individual_reactions)
    return model, loop_reactions, loop_individual_reactions


if __name__ == "__main__":
    from cobra.io import load_json_model, save_json_model
    import csv

    logging.basicConfig(level=logging.INFO)

    # model_name = "THG-beta-expanded_251118_transcriptomics"
    model_name = "iMM904"

    model = load_json_model(f"models/{model_name}.json")

    model.solver = "gurobi"
    model_no_loops, loops, loop_reactions = remove_loops(
        model, no_blocked_reactions=False
    )

    # Save comprehensive loop analysis
    with open(f"models/{model_name}_loop_analysis.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)

        # Write summary header
        writer.writerow(["Loop Removal Analysis Summary"])
        writer.writerow(["Total compacted loops identified", len(loops)])
        writer.writerow(["Total individual reactions to remove", len(loop_reactions)])
        writer.writerow([])

        # Write compacted loops section
        writer.writerow(["Compacted Loop ID", "Compacted Loop Full ID"])
        for loop in loops:
            writer.writerow([loop.id, loop.annotation.get("full_id", "")])

        writer.writerow([])

        # Write individual reactions section
        writer.writerow(["Individual Reactions to Remove"])
        for reac in loop_reactions:
            writer.writerow([reac])

    save_json_model(model_no_loops, f"models/{model_name}_no_loops.json")
