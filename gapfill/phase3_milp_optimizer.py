"""Phase-3 MILP-based blocked-reaction coverage optimizer.

This approach is based on the insight that adding a reversible PTR can only
expand the feasible flux space, never shrink it. Therefore:
1. We only need to check blocked reactions (not all reactions).
2. A blocked reaction becomes unblocked if there exists a feasible steady-state
   solution with that reaction carrying flux.
3. We formulate a MILP to maximize unblocked reactions using minimum PTRs.

Algorithm:
1. Identify blocked reactions B in the starting model.
2. For each candidate PTR p, determine which blocked reactions it can unblock
   (precompute coverage matrix).
3. Solve a MILP: maximize coverage with minimum PTRs.
4. Add selected PTRs to the model.
5. Iterate: recompute blocked set, repeat until no improvement.

Outputs:
 - `gapfill/files/selected_phase3_milp_connectors.csv`
 - `models/base/THG-beta-batch_251106_phase3_milp_connected.json`
"""

import csv
import os
import json
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from cobra.io import load_json_model, save_json_model
from cobra import Reaction, Model


def load_candidates(path):
    """Load candidate connectors from Phase-1 CSV."""
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            try:
                r["component1"] = (
                    int(r.get("component1"))
                    if r.get("component1") not in (None, "")
                    else None
                )
                r["component2"] = (
                    int(r.get("component2"))
                    if r.get("component2") not in (None, "")
                    else None
                )
            except Exception:
                pass
            rows.append(r)
    return rows


def add_transport(model, sel, idx, prefix="MILP"):
    """Add a reversible transport reaction to the model."""
    met1 = sel["met1"]
    met2 = sel["met2"]
    base = sel.get("base", "")
    s1 = sel.get("suffix1", "")
    s2 = sel.get("suffix2", "")
    raw_id = f"{prefix}_TRANS_{base}_{s1}_{s2}_{idx}"
    rid = "".join(ch if (ch.isalnum() or ch == "_") else "_" for ch in raw_id)
    try:
        m1 = model.metabolites.get_by_id(met1)
        m2 = model.metabolites.get_by_id(met2)
    except KeyError:
        return None
    rxn = Reaction(rid)
    rxn.name = f"MILP transport {met1} -> {met2}"
    rxn.add_metabolites({m1: -1.0, m2: 1.0})
    rxn.lower_bound = -1000.0
    rxn.upper_bound = 1000.0
    rxn.annotation = {"phase3_milp_selected": "true"}
    model.add_reactions([rxn])
    return rid


def find_blocked_reactions(model, reaction_list=None, zero_cutoff=None):
    """
    Find blocked reactions by testing if each can carry flux.
    A reaction is blocked if it cannot carry any flux in either direction.

    This is faster than cobra.flux_analysis.find_blocked_reactions for
    a subset of reactions because we only test the specified reactions.
    """
    if zero_cutoff is None:
        zero_cutoff = model.tolerance * 10
    if reaction_list is None:
        reaction_list = [r.id for r in model.reactions]

    blocked = set()

    for rid in reaction_list:
        try:
            rxn = model.reactions.get_by_id(rid)
        except KeyError:
            continue

        # Skip exchange/demand/sink reactions (typically not "blocked" in the usual sense)
        # but include them if they're in the list

        can_carry_flux = False

        # Test forward direction
        if rxn.upper_bound > zero_cutoff:
            with model:
                rxn.lower_bound = zero_cutoff  # force positive flux
                try:
                    sol = model.optimize()
                    if sol.status == "optimal":
                        can_carry_flux = True
                except Exception:
                    pass

        # Test reverse direction if forward failed
        if not can_carry_flux and rxn.lower_bound < -zero_cutoff:
            with model:
                rxn.upper_bound = -zero_cutoff  # force negative flux
                try:
                    sol = model.optimize()
                    if sol.status == "optimal":
                        can_carry_flux = True
                except Exception:
                    pass

        if not can_carry_flux:
            blocked.add(rid)

    return blocked


def find_blocked_reactions_fast(model, zero_cutoff=None):
    """
    Find all blocked reactions using FVA-like approach but more efficiently.
    Uses cobra's built-in find_blocked_reactions if available.
    """
    # Use model's tolerance if not specified (avoids cutoff < tolerance error)
    if zero_cutoff is None:
        zero_cutoff = model.tolerance * 10  # slightly above tolerance

    try:
        from cobra.flux_analysis import find_blocked_reactions as cobra_find_blocked

        blocked = cobra_find_blocked(model, zero_cutoff=zero_cutoff)
        return set(blocked)
    except Exception:
        # Fallback to our implementation
        return find_blocked_reactions(
            model, [r.id for r in model.reactions], zero_cutoff
        )


def test_candidate_coverage(model_json_path, candidate, blocked_rxns, zero_cutoff=None):
    """
    Test which blocked reactions a candidate PTR can unblock.
    This function is designed to be called in parallel.

    Returns: (candidate_key, set of unblocked reaction IDs)
    """
    model = load_json_model(model_json_path)
    if zero_cutoff is None:
        zero_cutoff = model.tolerance * 10

    # Add the candidate PTR
    met1 = candidate["met1"]
    met2 = candidate["met2"]
    base = candidate.get("base", "")

    try:
        m1 = model.metabolites.get_by_id(met1)
        m2 = model.metabolites.get_by_id(met2)
    except KeyError:
        return ((met1, met2), set())

    rxn = Reaction("_TEST_PTR_")
    rxn.add_metabolites({m1: -1.0, m2: 1.0})
    rxn.lower_bound = -1000.0
    rxn.upper_bound = 1000.0
    model.add_reactions([rxn])

    # Test which blocked reactions can now carry flux
    unblocked = set()
    for rid in blocked_rxns:
        try:
            r = model.reactions.get_by_id(rid)
        except KeyError:
            continue

        can_carry = False

        # Test forward
        if r.upper_bound > zero_cutoff:
            with model:
                r.lower_bound = zero_cutoff
                try:
                    sol = model.optimize()
                    if sol.status == "optimal":
                        can_carry = True
                except Exception:
                    pass

        # Test reverse
        if not can_carry and r.lower_bound < -zero_cutoff:
            with model:
                r.upper_bound = -zero_cutoff
                try:
                    sol = model.optimize()
                    if sol.status == "optimal":
                        can_carry = True
                except Exception:
                    pass

        if can_carry:
            unblocked.add(rid)

    return ((met1, met2), unblocked)


def compute_coverage_matrix(
    model, candidates, blocked_rxns, model_json_path=None, n_workers=None, verbose=True
):
    """
    Compute the coverage matrix: for each candidate, which blocked reactions can it unblock.

    Returns: dict mapping (met1, met2) -> set of unblocked reaction IDs
    """
    if model_json_path is None:
        # Save model to temp file for parallel workers
        import tempfile

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            model_json_path = f.name
        save_json_model(model, model_json_path)
        cleanup_temp = True
    else:
        cleanup_temp = False

    coverage = {}
    total = len(candidates)

    if verbose:
        print(
            f"Computing coverage matrix for {total} candidates over {len(blocked_rxns)} blocked reactions..."
        )

    # Use sequential processing for now (parallel can be added later)
    # Parallel processing with ProcessPoolExecutor can cause issues with some solvers
    for i, cand in enumerate(candidates):
        if verbose and (i + 1) % 50 == 0:
            print(f"  Progress: {i+1}/{total} candidates...")

        key, unblocked = test_candidate_coverage(model_json_path, cand, blocked_rxns)
        coverage[key] = unblocked

    if cleanup_temp:
        try:
            os.remove(model_json_path)
        except Exception:
            pass

    if verbose:
        total_coverage = sum(len(v) for v in coverage.values())
        non_empty = sum(1 for v in coverage.values() if v)
        print(
            f"Coverage computed: {non_empty}/{total} candidates can unblock at least one reaction"
        )
        print(f"Total (candidate, blocked_rxn) pairs covered: {total_coverage}")

    return coverage


def solve_coverage_milp(
    coverage,
    blocked_rxns,
    candidates,
    max_ptrs=None,
    min_unblocked=None,
    tradeoff_lambda=0.01,
    verbose=True,
):
    """
    Solve the set-cover MILP to select PTRs that maximize unblocked reactions
    while minimizing the number of PTRs.

    Formulation:
    - y_p ∈ {0,1}: whether PTR p is selected
    - z_b ∈ {0,1}: whether blocked reaction b is unblocked
    - Constraints: z_b ≤ Σ_{p covers b} y_p
    - Objective: maximize Σ z_b - λ Σ y_p

    Returns: list of selected candidate dicts
    """
    try:
        import gurobipy as gp
        from gurobipy import GRB

        solver = "gurobi"
    except ImportError:
        try:
            from scipy.optimize import milp, LinearConstraint, Bounds
            import numpy as np

            solver = "scipy"
        except ImportError:
            try:
                import pulp

                solver = "pulp"
            except ImportError:
                raise ImportError(
                    "No MILP solver available. Install gurobipy, scipy, or pulp."
                )

    # Build candidate index
    cand_idx = {}
    for i, c in enumerate(candidates):
        key = (c["met1"], c["met2"])
        cand_idx[key] = i

    # Build blocked reaction index
    blocked_list = list(blocked_rxns)
    blocked_idx = {b: i for i, b in enumerate(blocked_list)}

    n_cands = len(candidates)
    n_blocked = len(blocked_list)

    if verbose:
        print(
            f"Setting up MILP: {n_cands} candidate vars, {n_blocked} blocked-rxn vars"
        )

    if solver == "gurobi":
        return _solve_gurobi(
            coverage,
            blocked_list,
            blocked_idx,
            candidates,
            cand_idx,
            max_ptrs,
            min_unblocked,
            tradeoff_lambda,
            verbose,
        )
    elif solver == "pulp":
        return _solve_pulp(
            coverage,
            blocked_list,
            blocked_idx,
            candidates,
            cand_idx,
            max_ptrs,
            min_unblocked,
            tradeoff_lambda,
            verbose,
        )
    else:
        return _solve_scipy(
            coverage,
            blocked_list,
            blocked_idx,
            candidates,
            cand_idx,
            max_ptrs,
            min_unblocked,
            tradeoff_lambda,
            verbose,
        )


def _solve_gurobi(
    coverage,
    blocked_list,
    blocked_idx,
    candidates,
    cand_idx,
    max_ptrs,
    min_unblocked,
    tradeoff_lambda,
    verbose,
):
    """Solve using Gurobi."""
    import gurobipy as gp
    from gurobipy import GRB

    n_cands = len(candidates)
    n_blocked = len(blocked_list)

    m = gp.Model("coverage")
    if not verbose:
        m.Params.OutputFlag = 0

    # Variables
    y = m.addVars(n_cands, vtype=GRB.BINARY, name="y")  # PTR selection
    z = m.addVars(n_blocked, vtype=GRB.BINARY, name="z")  # blocked rxn unblocked

    # Constraints: z_b <= sum of y_p for p that covers b
    for b_id, b_i in blocked_idx.items():
        covering_cands = []
        for key, unblocked_set in coverage.items():
            if b_id in unblocked_set and key in cand_idx:
                covering_cands.append(cand_idx[key])

        if covering_cands:
            m.addConstr(
                z[b_i] <= gp.quicksum(y[p] for p in covering_cands), name=f"cov_{b_i}"
            )
        else:
            m.addConstr(z[b_i] == 0, name=f"nocov_{b_i}")

    # Optional constraints
    if max_ptrs is not None:
        m.addConstr(
            gp.quicksum(y[p] for p in range(n_cands)) <= max_ptrs, name="max_ptrs"
        )

    if min_unblocked is not None:
        m.addConstr(
            gp.quicksum(z[b] for b in range(n_blocked)) >= min_unblocked,
            name="min_unblocked",
        )

    # Objective: maximize unblocked - lambda * PTRs
    m.setObjective(
        gp.quicksum(z[b] for b in range(n_blocked))
        - tradeoff_lambda * gp.quicksum(y[p] for p in range(n_cands)),
        GRB.MAXIMIZE,
    )

    m.optimize()

    if m.Status != GRB.OPTIMAL:
        if verbose:
            print(f"MILP did not find optimal solution. Status: {m.Status}")
        return []

    # Extract selected candidates
    selected = []
    for i, c in enumerate(candidates):
        if y[i].X > 0.5:
            selected.append(c)

    if verbose:
        n_unblocked = sum(1 for b in range(n_blocked) if z[b].X > 0.5)
        print(
            f"MILP solution: {len(selected)} PTRs selected, {n_unblocked} reactions unblocked"
        )

    return selected


def _solve_pulp(
    coverage,
    blocked_list,
    blocked_idx,
    candidates,
    cand_idx,
    max_ptrs,
    min_unblocked,
    tradeoff_lambda,
    verbose,
):
    """Solve using PuLP."""
    import pulp

    n_cands = len(candidates)
    n_blocked = len(blocked_list)

    prob = pulp.LpProblem("coverage", pulp.LpMaximize)

    # Variables
    y = [pulp.LpVariable(f"y_{i}", cat="Binary") for i in range(n_cands)]
    z = [pulp.LpVariable(f"z_{i}", cat="Binary") for i in range(n_blocked)]

    # Constraints
    for b_id, b_i in blocked_idx.items():
        covering_cands = []
        for key, unblocked_set in coverage.items():
            if b_id in unblocked_set and key in cand_idx:
                covering_cands.append(cand_idx[key])

        if covering_cands:
            prob += z[b_i] <= pulp.lpSum(y[p] for p in covering_cands)
        else:
            prob += z[b_i] == 0

    if max_ptrs is not None:
        prob += pulp.lpSum(y) <= max_ptrs

    if min_unblocked is not None:
        prob += pulp.lpSum(z) >= min_unblocked

    # Objective
    prob += pulp.lpSum(z) - tradeoff_lambda * pulp.lpSum(y)

    prob.solve(pulp.PULP_CBC_CMD(msg=0 if not verbose else 1))

    if prob.status != pulp.LpStatusOptimal:
        if verbose:
            print(
                f"MILP did not find optimal solution. Status: {pulp.LpStatus[prob.status]}"
            )
        return []

    selected = []
    for i, c in enumerate(candidates):
        if y[i].value() > 0.5:
            selected.append(c)

    if verbose:
        n_unblocked = sum(1 for b in range(n_blocked) if z[b].value() > 0.5)
        print(
            f"MILP solution: {len(selected)} PTRs selected, {n_unblocked} reactions unblocked"
        )

    return selected


def _solve_scipy(
    coverage,
    blocked_list,
    blocked_idx,
    candidates,
    cand_idx,
    max_ptrs,
    min_unblocked,
    tradeoff_lambda,
    verbose,
):
    """Solve using scipy.optimize.milp (requires scipy >= 1.9)."""
    import numpy as np
    from scipy.optimize import milp, LinearConstraint, Bounds

    n_cands = len(candidates)
    n_blocked = len(blocked_list)
    n_vars = n_cands + n_blocked  # y variables then z variables

    # Objective: minimize -sum(z) + lambda*sum(y)  (we minimize, so negate maximize)
    c = np.zeros(n_vars)
    c[:n_cands] = tradeoff_lambda  # coefficient for y (PTR penalty)
    c[n_cands:] = -1.0  # coefficient for z (reward for unblocking)

    # Constraints: z_b - sum(y_p for p covering b) <= 0
    A_rows = []
    for b_id, b_i in blocked_idx.items():
        row = np.zeros(n_vars)
        row[n_cands + b_i] = 1.0  # z_b

        for key, unblocked_set in coverage.items():
            if b_id in unblocked_set and key in cand_idx:
                row[cand_idx[key]] = -1.0  # -y_p

        A_rows.append(row)

    A_ub = np.array(A_rows) if A_rows else np.zeros((0, n_vars))
    b_ub = np.zeros(len(A_rows))

    # Additional constraints
    A_eq_rows = []
    b_eq = []

    if max_ptrs is not None:
        row = np.zeros(n_vars)
        row[:n_cands] = 1.0
        A_rows.append(row)
        b_ub = np.append(b_ub, max_ptrs)
        A_ub = np.vstack([A_ub, row]) if A_ub.size > 0 else row.reshape(1, -1)

    # Bounds: all variables binary [0, 1]
    bounds = Bounds(lb=np.zeros(n_vars), ub=np.ones(n_vars))

    # Integrality: all variables are integers
    integrality = np.ones(n_vars, dtype=int)

    constraints = LinearConstraint(A_ub, -np.inf, b_ub) if A_ub.size > 0 else None

    result = milp(c, constraints=constraints, bounds=bounds, integrality=integrality)

    if not result.success:
        if verbose:
            print(f"MILP did not find optimal solution: {result.message}")
        return []

    selected = []
    for i, cand in enumerate(candidates):
        if result.x[i] > 0.5:
            selected.append(cand)

    if verbose:
        n_unblocked = sum(1 for b in range(n_blocked) if result.x[n_cands + b] > 0.5)
        print(
            f"MILP solution: {len(selected)} PTRs selected, {n_unblocked} reactions unblocked"
        )

    return selected


def run_phase3(
    candidates_csv,
    starting_model_json=None,
    out_dir=None,
    phase3_model_out=None,
    max_ptrs=None,
    min_unblocked=None,
    tradeoff_lambda=0.01,
    max_iterations=10,
    verbose=True,
):
    """
    Run the MILP-based Phase-3 optimizer with iterative refinement.

    Args:
        candidates_csv: Path to Phase-1 candidates CSV
        starting_model_json: Path to starting model (Phase-2 output)
        out_dir: Output directory for results
        phase3_model_out: Path for output model JSON
        max_ptrs: Maximum number of PTRs to add (None = no limit)
        min_unblocked: Minimum number of reactions to unblock (None = no limit)
        tradeoff_lambda: Weight for PTR penalty in objective (higher = fewer PTRs)
        max_iterations: Maximum number of iterative refinement rounds
        verbose: Print progress messages

    Returns:
        (output_csv_path, output_model_path, num_selected)
    """
    if out_dir is None:
        out_dir = os.path.join(os.path.dirname(__file__), "files")
    os.makedirs(out_dir, exist_ok=True)

    candidates = load_candidates(candidates_csv)
    if verbose:
        print(f"Loaded {len(candidates)} candidates from {candidates_csv}")

    if starting_model_json is None:
        starting_model_json = os.path.normpath(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "models",
                "base",
                "THG-beta-batch_251106_phase2_minimal_connected.json",
            )
        )

    model = load_json_model(starting_model_json)
    if verbose:
        print(
            f"Loaded model with {len(model.reactions)} reactions, {len(model.metabolites)} metabolites"
        )

    # Find initial blocked reactions
    if verbose:
        print("Finding initially blocked reactions...")
    blocked_rxns = find_blocked_reactions_fast(model)
    if verbose:
        print(f"Found {len(blocked_rxns)} blocked reactions")

    all_selected = []
    remaining_candidates = candidates.copy()
    iteration = 0

    while iteration < max_iterations and blocked_rxns and remaining_candidates:
        iteration += 1
        if verbose:
            print(f"\n=== Iteration {iteration} ===")
            print(
                f"Blocked reactions: {len(blocked_rxns)}, Remaining candidates: {len(remaining_candidates)}"
            )

        # Save model to temp file for coverage computation
        import tempfile

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            temp_model_path = f.name
        save_json_model(model, temp_model_path)

        # Compute coverage matrix
        coverage = compute_coverage_matrix(
            model,
            remaining_candidates,
            blocked_rxns,
            model_json_path=temp_model_path,
            verbose=verbose,
        )

        # Clean up temp file
        try:
            os.remove(temp_model_path)
        except Exception:
            pass

        # Filter out candidates with no coverage
        effective_candidates = [
            c
            for c in remaining_candidates
            if coverage.get((c["met1"], c["met2"]), set())
        ]

        if not effective_candidates:
            if verbose:
                print(
                    "No candidates can unblock any remaining blocked reactions. Stopping."
                )
            break

        if verbose:
            print(f"Effective candidates (with coverage): {len(effective_candidates)}")

        # Solve MILP
        selected = solve_coverage_milp(
            coverage,
            blocked_rxns,
            effective_candidates,
            max_ptrs=max_ptrs,
            min_unblocked=min_unblocked,
            tradeoff_lambda=tradeoff_lambda,
            verbose=verbose,
        )

        if not selected:
            if verbose:
                print("MILP returned no selections. Stopping.")
            break

        # Add selected PTRs to model
        for i, sel in enumerate(selected):
            rid = add_transport(model, sel, len(all_selected) + i, prefix="MILP")
            if rid:
                sel["reaction_id"] = rid
                all_selected.append(sel)

        if verbose:
            print(f"Added {len(selected)} PTRs to model")

        # Remove selected from remaining candidates
        selected_keys = {(s["met1"], s["met2"]) for s in selected}
        remaining_candidates = [
            c
            for c in remaining_candidates
            if (c["met1"], c["met2"]) not in selected_keys
        ]

        # Recompute blocked reactions
        new_blocked = find_blocked_reactions_fast(model)
        unblocked_count = len(blocked_rxns) - len(new_blocked)

        if verbose:
            print(f"Reactions unblocked this iteration: {unblocked_count}")
            print(f"Remaining blocked: {len(new_blocked)}")

        if len(new_blocked) >= len(blocked_rxns):
            if verbose:
                print("No improvement in blocked reactions. Stopping.")
            break

        blocked_rxns = new_blocked

    # Write output CSV
    out_csv = os.path.join(out_dir, "selected_phase3_milp_connectors.csv")
    with open(out_csv, "w", newline="") as f:
        fieldnames = ["reaction_id", "met1", "met2", "type", "base"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for sel in all_selected:
            writer.writerow(sel)

    # Save model
    if phase3_model_out is None:
        phase3_model_out = os.path.normpath(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "models",
                "base",
                "THG-beta-batch_251106_phase3_milp_connected.json",
            )
        )
    save_json_model(model, phase3_model_out)

    if verbose:
        print(f"\n=== Final Summary ===")
        print(f"Total PTRs added: {len(all_selected)}")
        print(f"Final blocked reactions: {len(blocked_rxns)}")
        print(f"Output CSV: {out_csv}")
        print(f"Output model: {phase3_model_out}")

    return out_csv, phase3_model_out, len(all_selected)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="MILP-based Phase-3 blocked reaction optimizer"
    )
    parser.add_argument("--candidates", default=None, help="Phase-1 candidates CSV")
    parser.add_argument("--model", default=None, help="Starting model JSON")
    parser.add_argument("--out", default=None, help="Output directory")
    parser.add_argument("--max-ptrs", type=int, default=None, help="Max PTRs to add")
    parser.add_argument(
        "--lambda",
        dest="tradeoff_lambda",
        type=float,
        default=0.01,
        help="Tradeoff weight (higher = fewer PTRs)",
    )
    parser.add_argument("--max-iter", type=int, default=10, help="Max iterations")
    args = parser.parse_args()

    cand_csv = args.candidates or os.path.join(
        os.path.dirname(__file__), "files", "candidates_all.csv"
    )
    if not os.path.exists(cand_csv):
        print("Candidates CSV not found:", cand_csv)
    else:
        out_csv, model_out, n = run_phase3(
            cand_csv,
            starting_model_json=args.model,
            out_dir=args.out,
            max_ptrs=args.max_ptrs,
            tradeoff_lambda=args.tradeoff_lambda,
            max_iterations=args.max_iter,
        )
        print(
            f"Phase3 MILP selected {n} connectors, wrote {out_csv} and model {model_out}"
        )
