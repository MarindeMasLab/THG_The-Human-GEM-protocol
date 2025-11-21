import logging
import numpy as np


def are_reactions_proportional(reaction1, reaction2):
    """
    Check if two reactions have stoichiometric coefficients that are proportional
    (multiplied by some factor).

    Parameters
    ----------
    reaction1 : cobra.Reaction
    reaction2 : cobra.Reaction

    Returns
    -------
    tuple
        (is_proportional, factor, is_reversed)
        is_proportional: bool - whether the reactions are proportional
        factor: float - the proportionality factor (reaction2 = factor * reaction1)
        is_reversed: bool - whether the reactions are reversed versions
    """
    tol = 1e-9

    # Check if reactions have the same metabolites
    if set(reaction1.metabolites.keys()) != set(reaction2.metabolites.keys()):
        return False, 0, False

    # Check if reactions are identical
    if reaction1.metabolites == reaction2.metabolites:
        return True, 1.0, False

    # Check if reactions are reversed
    if reaction1.metabolites == {
        met: -coeff for met, coeff in reaction2.metabolites.items()
    }:
        return True, -1.0, True

    # Check if reactions are proportional
    met_list = list(reaction1.metabolites.keys())
    if not met_list:
        return False, 0, False

    # Calculate the ratio for the first metabolite
    coeff1 = reaction1.metabolites[met_list[0]]
    coeff2 = reaction2.metabolites[met_list[0]]

    if abs(coeff1) < tol:  # Avoid division by zero
        return False, 0, False

    ratio = coeff2 / coeff1
    is_reversed = ratio < 0

    # Check if all metabolites maintain the same ratio
    for met in met_list[1:]:
        coeff1 = reaction1.metabolites[met]
        coeff2 = reaction2.metabolites[met]

        if abs(coeff1) < tol:
            if abs(coeff2) < tol:
                continue  # Both are essentially zero
            else:
                return False, 0, False

        current_ratio = coeff2 / coeff1
        if abs(current_ratio - ratio) > tol:
            return False, 0, False

    return True, ratio, is_reversed


def combine_identical_reactions(model_in, non_comp_id=None):
    """
    Combines reactions in a cobrapy Model that have the same stoichiometry
    or proportional stoichiometry recursively.

    Parameters
    ----------
    model_in : cobra.Model
        The cobra Model to modify.
    non_comp_id : set or list, optional
        Set of reaction IDs to exclude from combination.

    Returns
    ----------
    cobra.Model
        The modified cobra Model.
    """
    if non_comp_id is None:
        non_comp_id = set()
    else:
        non_comp_id = set(non_comp_id)

    tol = 1e-9
    model = model_in.copy()
    num_reactions_start = len(model.reactions)
    removed_reactions = []

    combined_count = 0  # Counter for combined reactions in this iteration

    for i, base_reaction in tqdm(
        enumerate(model.reactions), desc="Combining identical reactions"
    ):
        if base_reaction in removed_reactions:
            continue
        # Skip if base reaction is a boundary reaction
        if base_reaction.boundary:
            continue
        # Skip if base reaction is in non_comp_id
        if base_reaction.id in non_comp_id:
            continue
        for x in range(i + 1, len(model.reactions)):
            if model.reactions[x] in removed_reactions:
                continue
            # Skip if comparing reaction is a boundary reaction
            if model.reactions[x].boundary:
                continue
            # Skip if comparing reaction is in non_comp_id
            if model.reactions[x].id in non_comp_id:
                continue

            # Check if reactions are proportional
            is_proportional, factor, is_reversed = are_reactions_proportional(
                base_reaction, model.reactions[x]
            )

            if is_proportional:

                base_gpr = base_reaction.gene_reaction_rule
                add_gpr = model.reactions[x].gene_reaction_rule
                new_gpr = ""
                if base_gpr == "" and add_gpr != "":
                    new_gpr = model.reactions[x].gene_reaction_rule
                elif base_gpr != "" and add_gpr != "":
                    new_gpr = f"{base_gpr} or " f"({add_gpr})"
                elif base_gpr != "" and add_gpr == "":
                    new_gpr = base_gpr
                base_reaction.gene_reaction_rule = new_gpr
                removed_reactions.append(model.reactions[x])

                # Combine the bounds, accounting for the scaling factor
                if is_reversed:
                    # For reversed reactions: base + (-factor*other) with flipped bounds
                    base_reaction.bounds = (
                        min(
                            max(
                                base_reaction.lower_bound
                                - abs(factor) * model.reactions[x].upper_bound,
                                -1000,
                            ),
                            1000,
                        ),
                        max(
                            min(
                                base_reaction.upper_bound
                                - abs(factor) * model.reactions[x].lower_bound,
                                1000,
                            ),
                            -1000,
                        ),
                    )
                else:
                    # For same-direction reactions: base + (factor*other)
                    base_reaction.bounds = (
                        max(
                            base_reaction.lower_bound
                            + abs(factor) * model.reactions[x].lower_bound,
                            -1000,
                        ),
                        min(
                            base_reaction.upper_bound
                            + abs(factor) * model.reactions[x].upper_bound,
                            1000,
                        ),
                    )
                logging.debug(
                    f"Combined reactions: {base_reaction.id} and {model.reactions[x].id}, factor={factor}, reversed={is_reversed}"
                )

                # Combine reaction names
                base_reaction.name = f"{base_reaction.name} + {model.reactions[x].name}"

                # Get the full IDs from annotations if available
                # The full_id annotation is created by the compaction process
                base_id = base_reaction._annotation.get("full_id", base_reaction.id)
                other_id = model.reactions[x]._annotation.get(
                    "full_id", model.reactions[x].id
                )

                # Check if reactions lack full_id when they should have one
                if (
                    len(base_reaction.id) == 100
                    and "full_id" not in base_reaction._annotation
                ):
                    logging.warning(
                        f"Reaction {base_reaction.id} has truncated ID but no full_id in annotations"
                    )
                if (
                    len(model.reactions[x].id) == 100
                    and "full_id" not in model.reactions[x]._annotation
                ):
                    logging.warning(
                        f"Reaction {model.reactions[x].id} has truncated ID but no full_id in annotations"
                    )

                # Create new ID using full IDs when available
                if (
                    "@" in base_id
                    or "@" in other_id
                    or "#" in base_id
                    or "#" in other_id
                ):
                    new_id = f"({base_id})@({other_id})"
                else:
                    new_id = f"{base_id}@{other_id}"

                if len(new_id) > 100:
                    base_reaction.id = new_id[:100]
                    base_reaction._annotation["full_id"] = new_id
                else:
                    base_reaction.id = new_id
                    # Keep full_id annotation if original reactions had long IDs
                    if len(base_id) == 100 or len(other_id) == 100:
                        base_reaction._annotation["full_id"] = new_id
                combined_count += 1

    # Remove the combined reactions
    model.remove_reactions(removed_reactions)

    # Recursively call the function if reactions were combined in this iteration
    if combined_count > 0:
        logging.info(
            f"Combining identical reactions: combined {combined_count} pairs, removed {len(removed_reactions)} reactions"
        )
        model, additional_removed = combine_identical_reactions(
            model, non_comp_id=non_comp_id
        )  # Recursive call
        removed_reactions.extend(additional_removed)
    else:
        logging.debug("No more identical reactions to combine.")

    return model, removed_reactions


import cobra
import cobra.manipulation
import logging
import pickle
from tqdm import tqdm
from decimal import Decimal, getcontext

# Set precision for Decimal calculations
getcontext().prec = 50


def reaction_comb(
    base_rxn: cobra.Reaction, add_rxn: cobra.Reaction, flip=False
) -> None:
    """
    Combines two reactions into a single reaction in place.

    Parameters
    ----------
    base_rxn : cobra.Reaction
        The base reaction.
    add_rxn : cobra.Reaction
        The reaction to be added to the base reaction.
    flip : bool, optional
        Flag indicating whether to flip the add_rxn direction. Default is False.
    """
    # Combine gene-reaction rules
    new_rule = (
        base_rxn.grr_temp
        if hasattr(base_rxn, "grr_temp")
        else base_rxn.gene_reaction_rule
    )
    if (
        "@" in base_rxn.id
        or "@" in add_rxn.id
        or "#" in base_rxn.id
        or "#" in add_rxn.id
    ):
        new_id = f"({base_rxn.id})@({add_rxn.id})"
    else:
        new_id = f"{base_rxn.id}@{add_rxn.id}"

    rule1 = add_rxn.gene_reaction_rule.strip()
    logging.debug(new_id)
    if rule1 and new_rule:
        new_rule = f"{new_rule} and ({rule1})"
    elif rule1:
        new_rule = rule1

    base_rxn.grr_temp = new_rule

    base_rxn.name = f"{base_rxn.name} + {add_rxn.name}"

    # Get the full IDs from annotations if available
    # The full_id annotation is created by the compaction process
    base_id = base_rxn._annotation.get("full_id", base_rxn.id)
    add_id = add_rxn._annotation.get("full_id", add_rxn.id)

    # Check if the reaction lacks full_id when it should have one
    if len(base_rxn.id) == 100 and "full_id" not in base_rxn._annotation:
        logging.warning(
            f"Reaction {base_rxn.id} has truncated ID but no full_id in annotations"
        )
    if len(add_rxn.id) == 100 and "full_id" not in add_rxn._annotation:
        logging.warning(
            f"Reaction {add_rxn.id} has truncated ID but no full_id in annotations"
        )

    # Create new ID using full IDs when available
    if "@" in base_id or "@" in add_id or "#" in base_id or "#" in add_id:
        new_id = f"({base_id})@({add_id})"
    else:
        new_id = f"{base_id}@{add_id}"

    if len(new_id) > 100:
        base_rxn.id = new_id[:100]
        base_rxn._annotation["full_id"] = new_id
    else:
        base_rxn.id = new_id
        # Keep full_id annotation if original reactions had long IDs
        if len(base_id) == 100 or len(add_id) == 100:
            base_rxn._annotation["full_id"] = new_id
        elif "full_id" in base_rxn._annotation:
            # Remove full_id only if neither original reaction had a truncated ID
            del base_rxn._annotation["full_id"]

    if flip == False:
        # Combine metabolites with higher precision
        base_rxn.add_metabolites(
            {met: float(Decimal(coeff)) for met, coeff in add_rxn.metabolites.items()},
            combine=True,
        )

        # Adjust bounds with higher precision
        base_rxn.bounds = (
            float(max(Decimal(base_rxn.lower_bound), Decimal(add_rxn.lower_bound))),
            float(min(Decimal(base_rxn.upper_bound), Decimal(add_rxn.upper_bound))),
        )
    else:
        # if flip True
        # Flip the reaction
        base_rxn.add_metabolites(
            {met: -coeff for met, coeff in add_rxn.metabolites.items()}, combine=True
        )
        base_rxn.bounds = (
            float(max(Decimal(base_rxn.lower_bound), -Decimal(add_rxn.upper_bound))),
            float(min(Decimal(base_rxn.upper_bound), -Decimal(add_rxn.lower_bound))),
        )

    # check if the reaction has no metabolites
    if not base_rxn.metabolites:
        logging.warning(f"Reaction {base_rxn.id} has no metabolites")


def simple_compact(model: cobra.Model, non_comp_id=None) -> cobra.Model:
    """
    Compacts the metabolic model by combining reactions connected via metabolites that have exactly two reactions.

    Parameters
    ----------
    model : cobra.Model
        The COBRApy metabolic model to compact.
    non_comp_id : list, optional
        List of reaction IDs to exclude from compaction.

    Returns
    -------
    cobra.Model
        The compacted metabolic model.
    """
    if non_comp_id is None:
        non_comp_id = set(reaction.id for reaction in model.boundary)
        objective_reaction = next(
            (
                rxn
                for rxn in model.reactions
                if rxn.flux_expression == model.objective.expression
            ),
            None,
        )
        if objective_reaction:
            logging.info(f"Objective reaction: {objective_reaction.id}")
            non_comp_id.add(objective_reaction.id)
        else:
            logging.warning("No objective reaction found")

    loop_iteration = 0
    while True:
        loop_iteration += 1
        logging.info(f"Loop iteration: {loop_iteration}")
        compaction_flag = False

        # Prune unused metabolites
        model, _ = cobra.manipulation.prune_unused_metabolites(model)

        # Optimize model once per loop
        optimization = model.optimize()
        if optimization.objective_value == 0 or optimization.status == "infeasible":
            logging.error("Model is infeasible after pruning.")
            raise ValueError("Model is infeasible")

        metabolites = list(model.metabolites)
        if logging.getLogger().level == logging.INFO:
            metabolites = tqdm(metabolites, desc="Processing metabolites")
        for metabolite in metabolites:
            # Check for metabolites with only one reaction (dead-end)
            if len(metabolite.reactions) == 1:
                reaction = list(metabolite.reactions)[0]
                # Only remove if it's not a boundary reaction
                if reaction not in model.boundary and reaction.id not in non_comp_id:
                    logging.debug(f"Removing dead-end reaction: {reaction.id}")
                    model.remove_reactions([reaction])
                    compaction_flag = True
                continue

            if len(metabolite.reactions) != 2:
                continue

            reactions = list(metabolite.reactions)
            reactions_ids = {rxn.id for rxn in reactions}

            if reactions_ids & non_comp_id:
                continue

            reaction_0, reaction_1 = reactions
            metab_coeff_0 = Decimal(reaction_0.get_coefficient(metabolite))
            metab_coeff_1 = Decimal(reaction_1.get_coefficient(metabolite))

            # Determine reaction directions and combine accordingly
            if (metab_coeff_0 < 0 and metab_coeff_1 > 0) or (
                metab_coeff_0 > 0 and metab_coeff_1 < 0
            ):
                if abs(metab_coeff_0) == abs(metab_coeff_1):
                    base_rxn, add_rxn = (
                        (reaction_0, reaction_1)
                        if metab_coeff_0 > 0
                        else (reaction_1, reaction_0)
                    )
                    reaction_comb(base_rxn, add_rxn)
                    model.remove_reactions([add_rxn])
                    compaction_flag = True
                else:
                    # Normalize coefficients with higher precision
                    for rxn in [reaction_0, reaction_1]:
                        coeff = Decimal(rxn.get_coefficient(metabolite))
                        norm_factor = abs(coeff)
                        rxn.add_metabolites(
                            {
                                met: float(Decimal(c) / norm_factor)
                                for met, c in rxn.metabolites.items()
                            },
                            combine=False,
                        )
                        rxn.bounds = (
                            float(Decimal(rxn.lower_bound) * norm_factor),
                            float(Decimal(rxn.upper_bound) * norm_factor),
                        )

                    base_rxn, add_rxn = (
                        (reaction_0, reaction_1)
                        if metab_coeff_0 > 0
                        else (reaction_1, reaction_0)
                    )
                    reaction_comb(base_rxn, add_rxn)
                    model.remove_reactions([add_rxn])
                    compaction_flag = True
            elif (metab_coeff_0 < 0 and metab_coeff_1 < 0) or (
                metab_coeff_0 > 0 and metab_coeff_1 > 0
            ):
                # Flip reaction_0 coefficients with higher precision
                for met, coeff in reaction_0.metabolites.items():
                    reaction_0.add_metabolites(
                        {met: float(Decimal(coeff) / abs(metab_coeff_0))}, combine=False
                    )
                reaction_0.bounds = (
                    float(Decimal(reaction_0.lower_bound) * abs(metab_coeff_0)),
                    float(Decimal(reaction_0.upper_bound) * abs(metab_coeff_0)),
                )

                # Normalize reaction_1 coefficients with higher precision
                for met, coeff in reaction_1.metabolites.items():
                    reaction_1.add_metabolites(
                        {met: float(Decimal(coeff) / abs(metab_coeff_1))}, combine=False
                    )
                reaction_1.bounds = (
                    float(Decimal(reaction_1.lower_bound) * abs(metab_coeff_1)),
                    float(Decimal(reaction_1.upper_bound) * abs(metab_coeff_1)),
                )

                reaction_comb(reaction_0, reaction_1, flip=True)
                model.remove_reactions([reaction_1])
                compaction_flag = True
            else:
                logging.warning(
                    f"Metabolite coefficients are not in the expected format: "
                    f"{metab_coeff_0}, {metab_coeff_1} in reactions {reaction_0.id}, {reaction_1.id}"
                )

        if not compaction_flag:
            logging.info(
                f"Iteration {loop_iteration}: No further compaction possible. Final model has {len(model.reactions)} reactions and {len(model.metabolites)} metabolites."
            )
            for reaction in model.reactions:
                # Restore gene_reaction_rule for compacted reactions
                # It is done this way to improve performance
                if hasattr(reaction, "grr_temp"):
                    reaction.gene_reaction_rule = reaction.grr_temp
                    del reaction.grr_temp
            break
        else:
            logging.info(
                f"Iteration {loop_iteration}: Compacted model to {len(model.reactions)} reactions and {len(model.metabolites)} metabolites."
            )
    model, _ = cobra.manipulation.prune_unused_metabolites(model)
    model, loop_reactions = cobra.manipulation.prune_unused_reactions(model)
    if loop_reactions:
        logging.info(f"Removed {len(loop_reactions)} loop reactions")
        logging.debug(
            f"Loop reactions removed: {[loop_reaction.id for loop_reaction in loop_reactions]}"
        )
    model.repair()
    return model, loop_reactions


def full_compaction(
    model: cobra.Model, non_comp_id=None, no_blocked_reactions=False
) -> cobra.Model:
    """
    Performs full compaction of the metabolic model by first combining identical reactions
    and then applying simple compaction.

    Parameters
    ----------
    model : cobra.Model
        The COBRApy metabolic model to compact.
    non_comp_id : list, optional
        List of reaction IDs to exclude from compaction.
    no_blocked_reactions : bool, optional
        If True, it is assumed there are no blocked reactions in the model and 
        they are not removed before compaction.

    Returns
    -------
    cobra.Model
        The fully compacted metabolic model.
    """
    temp_model = model.copy()

    if no_blocked_reactions == False:
        logging.info("Removing blocked reactions before compaction")
        from cobra.flux_analysis import find_blocked_reactions

        blocked_reactions = find_blocked_reactions(temp_model)
        temp_model.remove_reactions(blocked_reactions)
        logging.info(f"Removed {len(blocked_reactions)} blocked reactions")

    initial_rxns = len(temp_model.reactions)
    initial_mets = len(temp_model.metabolites)
    logging.info(
        f"Starting full compaction: {initial_rxns} reactions, {initial_mets} metabolites"
    )

    len_rxns = initial_rxns
    loop_reactions = []
    compaction_round = 0

    while True:
        compaction_round += 1
        logging.info(f"\n=== Full compaction round {compaction_round} ===")

        temp_model, loop_reactions_0 = simple_compact(
            temp_model, non_comp_id=non_comp_id
        )
        loop_reactions.extend(loop_reactions_0)
        temp_model, _ = combine_identical_reactions(temp_model, non_comp_id=non_comp_id)
        new_len_rxns = len(temp_model.reactions)

        if new_len_rxns == len_rxns:
            logging.info(f"Full compaction complete after {compaction_round} rounds")
            break

        logging.info(
            f"Round {compaction_round} reduced reactions from {len_rxns} to {new_len_rxns}"
        )
        len_rxns = new_len_rxns

    final_rxns = len(temp_model.reactions)
    final_mets = len(temp_model.metabolites)

    logging.info("\n" + "=" * 60)
    logging.info("FULL COMPACTION SUMMARY")
    logging.info("=" * 60)
    logging.info(
        f"Reactions:   {initial_rxns} → {final_rxns} ({initial_rxns - final_rxns} removed, {100*(initial_rxns-final_rxns)/initial_rxns:.1f}% reduction)"
    )
    logging.info(
        f"Metabolites: {initial_mets} → {final_mets} ({initial_mets - final_mets} removed, {100*(initial_mets-final_mets)/initial_mets:.1f}% reduction)"
    )
    logging.info(f"Loop reactions removed: {len(loop_reactions)}")
    if loop_reactions:
        logging.debug(f"Loop reaction IDs: {[rxn.id for rxn in loop_reactions]}")
    logging.info("=" * 60)

    return temp_model, loop_reactions


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    from cobra.io import load_json_model, save_json_model

    model = load_json_model("pipeline/models/endoC_final_5.json")
