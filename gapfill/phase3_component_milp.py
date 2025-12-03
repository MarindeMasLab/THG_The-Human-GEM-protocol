"""Phase-3 Component-wise MILP optimizer.

This approach solves the gap-filling problem by treating each isolated component
separately, dramatically reducing the computational burden.

Key insight: The main component (ID=1) is already well-connected. We only need
to maximize connectivity of isolated components (ID >= 2).

Algorithm:
1. Identify reactions belonging to each component (via their metabolites)
2. For each isolated component (in decreasing size order):
   a. Find blocked reactions IN that component
   b. Find candidate PTRs that connect TO that component
   c. Compute coverage matrix (small: ~few hundred candidates × ~few hundred blocked)
   d. Solve MILP to maximize unblocked with minimum PTRs
   e. Add selected PTRs to model
   f. Remove used candidates from pool
3. Output final model with all selected PTRs

This is much faster than the global MILP because:
- Each sub-problem is much smaller
- We can parallelize across components if needed
- Candidates are filtered to only those relevant to each component
"""

import csv
import os
import json
import time
import pickle
import networkx as nx
from collections import defaultdict
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


def load_components_summary(path):
    """Load component summary CSV."""
    components = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            cid = int(r["component_id"])
            size = int(r["num_nodes"])
            components[cid] = size
    return components


def build_metabolite_to_component_map(model, components_file=None, bipartite_pkl=None):
    """
    Build a mapping from metabolite IDs to their component IDs.

    If bipartite_pkl is provided, load the pre-computed component assignments.
    Otherwise, rebuild from scratch.
    """
    met_to_comp = {}

    # Try to load from pickle if available
    if bipartite_pkl and os.path.exists(bipartite_pkl):
        with open(bipartite_pkl, "rb") as f:
            data = pickle.load(f)
        # The pickle should contain node-to-component mapping
        if isinstance(data, dict) and "node_to_component" in data:
            met_to_comp = data["node_to_component"]
            return {k: v for k, v in met_to_comp.items() if k.startswith("MAM")}

    # Rebuild bipartite graph and compute components
    G = nx.Graph()
    for rxn in model.reactions:
        rid = rxn.id
        G.add_node(rid, bipartite=0)  # reaction node
        for met in rxn.metabolites:
            mid = met.id
            G.add_node(mid, bipartite=1)  # metabolite node
            G.add_edge(rid, mid)

    # Find connected components
    components = list(nx.connected_components(G))
    # Sort by size descending (largest = main component = ID 1)
    components.sort(key=len, reverse=True)

    for comp_id, comp_nodes in enumerate(components, start=1):
        for node in comp_nodes:
            if node.startswith("MAM"):  # metabolite
                met_to_comp[node] = comp_id

    return met_to_comp


def get_reactions_in_component(model, component_id, met_to_comp):
    """Get all reactions that have metabolites in the given component."""
    rxns_in_comp = set()
    for rxn in model.reactions:
        for met in rxn.metabolites:
            if met_to_comp.get(met.id) == component_id:
                rxns_in_comp.add(rxn.id)
                break
    return rxns_in_comp


def get_candidates_for_component(candidates, component_id, used_candidates=None):
    """
    Get candidates that connect TO the given component.
    A candidate connects to component X if component1 == X or component2 == X.

    Excludes candidates that only connect within the main component (1-to-1).
    """
    if used_candidates is None:
        used_candidates = set()

    relevant = []
    for c in candidates:
        key = (c["met1"], c["met2"])
        if key in used_candidates:
            continue

        c1 = c.get("component1")
        c2 = c.get("component2")

        # Skip if both endpoints are in main component (1)
        if c1 == 1 and c2 == 1:
            continue

        # Include if either endpoint is in target component
        if c1 == component_id or c2 == component_id:
            relevant.append(c)

    return relevant


def find_blocked_reactions_in_set(model, reaction_ids, zero_cutoff=None):
    """
    Find which reactions in the given set are blocked.
    A reaction is blocked if it cannot carry flux in either direction.
    """
    if zero_cutoff is None:
        zero_cutoff = model.tolerance * 10

    blocked = set()

    for rid in reaction_ids:
        try:
            rxn = model.reactions.get_by_id(rid)
        except KeyError:
            continue

        can_carry_flux = False

        # Test forward direction
        if rxn.upper_bound > zero_cutoff:
            with model:
                rxn.lower_bound = zero_cutoff
                try:
                    sol = model.optimize()
                    if sol.status == "optimal":
                        can_carry_flux = True
                except Exception:
                    pass

        # Test reverse direction if forward failed
        if not can_carry_flux and rxn.lower_bound < -zero_cutoff:
            with model:
                rxn.upper_bound = -zero_cutoff
                try:
                    sol = model.optimize()
                    if sol.status == "optimal":
                        can_carry_flux = True
                except Exception:
                    pass

        if not can_carry_flux:
            blocked.add(rid)

    return blocked


def add_transport(model, sel, idx, prefix="CMILP"):
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
    rxn.name = f"Component MILP transport {met1} -> {met2}"
    rxn.add_metabolites({m1: -1.0, m2: 1.0})
    rxn.lower_bound = -1000.0
    rxn.upper_bound = 1000.0
    rxn.annotation = {"phase3_component_milp": "true"}
    model.add_reactions([rxn])
    return rid


def test_candidate_unblocks(model, candidate, blocked_rxns, zero_cutoff=None):
    """
    Test which blocked reactions a candidate PTR can unblock.
    Returns set of unblocked reaction IDs.
    """
    if zero_cutoff is None:
        zero_cutoff = model.tolerance * 10

    met1 = candidate["met1"]
    met2 = candidate["met2"]

    try:
        m1 = model.metabolites.get_by_id(met1)
        m2 = model.metabolites.get_by_id(met2)
    except KeyError:
        return set()

    # Add temporary PTR
    with model:
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

    return unblocked


def compute_coverage_for_component(model, candidates, blocked_rxns, verbose=True):
    """
    Compute coverage matrix for a single component's candidates and blocked reactions.

    Returns: dict mapping (met1, met2) -> set of unblocked reaction IDs
    """
    coverage = {}
    total = len(candidates)
    n_blocked = len(blocked_rxns)

    if verbose:
        print(
            f"    Computing coverage: {total} candidates × {n_blocked} blocked reactions"
        )

    start_time = time.time()

    for i, cand in enumerate(candidates):
        if verbose and (i + 1) % 20 == 0:
            elapsed = time.time() - start_time
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            eta = (total - i - 1) / rate if rate > 0 else 0
            print(
                f"      Progress: {i+1}/{total} ({rate:.1f} cand/sec, ETA: {eta:.0f}s)"
            )

        key = (cand["met1"], cand["met2"])
        unblocked = test_candidate_unblocks(model, cand, blocked_rxns)
        coverage[key] = unblocked

    if verbose:
        elapsed = time.time() - start_time
        non_empty = sum(1 for v in coverage.values() if v)
        total_pairs = sum(len(v) for v in coverage.values())
        print(
            f"    Coverage done in {elapsed:.1f}s: {non_empty}/{total} candidates useful, {total_pairs} total unblocks"
        )

    return coverage


def solve_component_milp(
    coverage, blocked_rxns, candidates, tradeoff_lambda=0.01, verbose=True
):
    """
    Solve MILP for a single component.

    Objective: maximize unblocked reactions - lambda * number of PTRs
    """
    # Try different solvers in order of preference
    try:
        import gurobipy

        return _solve_gurobi(
            coverage, blocked_rxns, candidates, tradeoff_lambda, verbose
        )
    except ImportError:
        pass

    try:
        import pulp

        return _solve_pulp(coverage, blocked_rxns, candidates, tradeoff_lambda, verbose)
    except ImportError:
        pass

    try:
        from scipy.optimize import milp

        return _solve_scipy(
            coverage, blocked_rxns, candidates, tradeoff_lambda, verbose
        )
    except ImportError:
        pass

    raise ImportError("No MILP solver available. Install gurobipy, pulp, or scipy>=1.9")


def _solve_gurobi(coverage, blocked_rxns, candidates, tradeoff_lambda, verbose):
    """Solve using Gurobi (fastest commercial solver)."""
    import gurobipy as gp
    from gurobipy import GRB

    blocked_list = list(blocked_rxns)
    blocked_idx = {b: i for i, b in enumerate(blocked_list)}
    cand_idx = {(c["met1"], c["met2"]): i for i, c in enumerate(candidates)}

    n_cands = len(candidates)
    n_blocked = len(blocked_list)

    m = gp.Model("component_coverage")
    m.Params.OutputFlag = 1 if verbose else 0

    # Variables
    y = m.addVars(n_cands, vtype=GRB.BINARY, name="y")  # PTR selection
    z = m.addVars(n_blocked, vtype=GRB.BINARY, name="z")  # blocked rxn unblocked

    # Constraints: z_b <= sum of y_p for p that covers b
    for b_id, b_i in blocked_idx.items():
        covering = [
            cand_idx[key]
            for key, unblocked in coverage.items()
            if b_id in unblocked and key in cand_idx
        ]
        if covering:
            m.addConstr(z[b_i] <= gp.quicksum(y[p] for p in covering))
        else:
            m.addConstr(z[b_i] == 0)

    # Objective: maximize unblocked - lambda * PTRs
    m.setObjective(
        gp.quicksum(z[b] for b in range(n_blocked))
        - tradeoff_lambda * gp.quicksum(y[p] for p in range(n_cands)),
        GRB.MAXIMIZE,
    )

    m.optimize()

    if m.Status != GRB.OPTIMAL:
        return []

    selected = [c for i, c in enumerate(candidates) if y[i].X > 0.5]

    if verbose:
        n_unblocked = sum(1 for b in range(n_blocked) if z[b].X > 0.5)
        print(f"    MILP: {len(selected)} PTRs → {n_unblocked} unblocked")

    return selected


def _solve_pulp(coverage, blocked_rxns, candidates, tradeoff_lambda, verbose):
    """Solve using PuLP (open source)."""
    import pulp

    blocked_list = list(blocked_rxns)
    blocked_idx = {b: i for i, b in enumerate(blocked_list)}
    cand_idx = {(c["met1"], c["met2"]): i for i, c in enumerate(candidates)}

    n_cands = len(candidates)
    n_blocked = len(blocked_list)

    prob = pulp.LpProblem("component_coverage", pulp.LpMaximize)

    y = [pulp.LpVariable(f"y_{i}", cat="Binary") for i in range(n_cands)]
    z = [pulp.LpVariable(f"z_{i}", cat="Binary") for i in range(n_blocked)]

    for b_id, b_i in blocked_idx.items():
        covering = [
            cand_idx[key]
            for key, unblocked in coverage.items()
            if b_id in unblocked and key in cand_idx
        ]
        if covering:
            prob += z[b_i] <= pulp.lpSum(y[p] for p in covering)
        else:
            prob += z[b_i] == 0

    prob += pulp.lpSum(z) - tradeoff_lambda * pulp.lpSum(y)

    prob.solve(pulp.PULP_CBC_CMD(msg=1 if verbose else 0))

    if prob.status != pulp.LpStatusOptimal:
        return []

    selected = [c for i, c in enumerate(candidates) if y[i].value() > 0.5]

    if verbose:
        n_unblocked = sum(1 for b in range(n_blocked) if z[b].value() > 0.5)
        print(f"    MILP: {len(selected)} PTRs → {n_unblocked} unblocked")

    return selected


def _solve_scipy(coverage, blocked_rxns, candidates, tradeoff_lambda, verbose):
    """Solve using scipy.optimize.milp."""
    import numpy as np
    from scipy.optimize import milp, LinearConstraint, Bounds

    blocked_list = list(blocked_rxns)
    blocked_idx = {b: i for i, b in enumerate(blocked_list)}
    cand_idx = {(c["met1"], c["met2"]): i for i, c in enumerate(candidates)}

    n_cands = len(candidates)
    n_blocked = len(blocked_list)
    n_vars = n_cands + n_blocked

    # Objective: minimize -sum(z) + lambda*sum(y)
    c = np.zeros(n_vars)
    c[:n_cands] = tradeoff_lambda
    c[n_cands:] = -1.0

    # Constraints: z_b - sum(y_p covering b) <= 0
    A_rows = []
    for b_id, b_i in blocked_idx.items():
        row = np.zeros(n_vars)
        row[n_cands + b_i] = 1.0
        for key, unblocked in coverage.items():
            if b_id in unblocked and key in cand_idx:
                row[cand_idx[key]] = -1.0
        A_rows.append(row)

    A_ub = np.array(A_rows) if A_rows else np.zeros((1, n_vars))
    b_ub = np.zeros(len(A_rows)) if A_rows else np.zeros(1)

    bounds = Bounds(lb=np.zeros(n_vars), ub=np.ones(n_vars))
    integrality = np.ones(n_vars, dtype=int)

    constraints = LinearConstraint(A_ub, -np.inf, b_ub)

    result = milp(c, constraints=constraints, bounds=bounds, integrality=integrality)

    if not result.success:
        return []

    selected = [c for i, c in enumerate(candidates) if result.x[i] > 0.5]

    if verbose:
        n_unblocked = sum(1 for b in range(n_blocked) if result.x[n_cands + b] > 0.5)
        print(f"    MILP: {len(selected)} PTRs → {n_unblocked} unblocked")

    return selected


def run_phase3_component_wise(
    candidates_csv,
    starting_model_json=None,
    components_csv=None,
    out_dir=None,
    phase3_model_out=None,
    tradeoff_lambda=0.01,
    min_component_size=4,
    max_components=None,
    verbose=True,
):
    """
    Run Phase-3 optimization component by component.

    Args:
        candidates_csv: Path to Phase-1 candidates CSV
        starting_model_json: Path to starting model (Phase-2 output)
        components_csv: Path to components summary CSV
        out_dir: Output directory
        phase3_model_out: Output model path
        tradeoff_lambda: Weight for PTR penalty (higher = fewer PTRs)
        min_component_size: Only process components with at least this many nodes
        max_components: Max number of components to process (None = all)
        verbose: Print progress

    Returns:
        (output_csv_path, output_model_path, total_selected)
    """
    # Setup paths
    base_dir = os.path.dirname(__file__)
    if out_dir is None:
        out_dir = os.path.join(base_dir, "files")
    os.makedirs(out_dir, exist_ok=True)

    if starting_model_json is None:
        starting_model_json = os.path.normpath(
            os.path.join(
                base_dir,
                "..",
                "models",
                "base",
                "THG-beta-batch_251106_phase2_minimal_connected.json",
            )
        )

    if components_csv is None:
        components_csv = os.path.join(out_dir, "components_summary.csv")

    # Load data
    candidates = load_candidates(candidates_csv)
    components = load_components_summary(components_csv)
    model = load_json_model(starting_model_json)

    if verbose:
        print(f"Loaded {len(candidates)} candidates")
        print(f"Loaded {len(components)} components")
        print(
            f"Model: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites"
        )

    # Build metabolite-to-component mapping
    if verbose:
        print("Building metabolite-to-component mapping...")
    met_to_comp = build_metabolite_to_component_map(model)
    if verbose:
        print(f"  Mapped {len(met_to_comp)} metabolites to components")

    # Get isolated components (exclude main component = ID 1)
    isolated_comps = [
        (cid, size)
        for cid, size in components.items()
        if cid > 1 and size >= min_component_size
    ]
    isolated_comps.sort(key=lambda x: x[1], reverse=True)  # largest first

    if max_components:
        isolated_comps = isolated_comps[:max_components]

    if verbose:
        print(
            f"\nProcessing {len(isolated_comps)} isolated components (size >= {min_component_size})"
        )
        print(f"Main component (ID=1) has {components.get(1, 0)} nodes - EXCLUDED")

    # Track used candidates and all selections
    used_candidates = set()
    all_selected = []
    component_results = []

    total_start = time.time()

    for comp_idx, (comp_id, comp_size) in enumerate(isolated_comps):
        if verbose:
            print(f"\n{'='*60}")
            print(
                f"Component {comp_id} ({comp_idx+1}/{len(isolated_comps)}): {comp_size} nodes"
            )
            print(f"{'='*60}")

        # Get reactions in this component
        rxns_in_comp = get_reactions_in_component(model, comp_id, met_to_comp)
        if verbose:
            print(f"  Reactions in component: {len(rxns_in_comp)}")

        if not rxns_in_comp:
            if verbose:
                print(f"  No reactions found, skipping")
            continue

        # Find blocked reactions in this component
        blocked_in_comp = find_blocked_reactions_in_set(model, rxns_in_comp)
        if verbose:
            print(f"  Blocked reactions: {len(blocked_in_comp)}")

        if not blocked_in_comp:
            if verbose:
                print(f"  No blocked reactions, skipping")
            continue

        # Get candidates for this component
        comp_candidates = get_candidates_for_component(
            candidates, comp_id, used_candidates
        )
        if verbose:
            print(f"  Available candidates: {len(comp_candidates)}")

        if not comp_candidates:
            if verbose:
                print(f"  No candidates available, skipping")
            continue

        # Compute coverage matrix (the expensive part, but now much smaller!)
        coverage = compute_coverage_for_component(
            model, comp_candidates, blocked_in_comp, verbose
        )

        # Filter to candidates with any coverage
        effective_candidates = [
            c for c in comp_candidates if coverage.get((c["met1"], c["met2"]), set())
        ]

        if not effective_candidates:
            if verbose:
                print(f"  No candidates can unblock any reactions, skipping")
            continue

        if verbose:
            print(f"  Effective candidates: {len(effective_candidates)}")

        # Solve MILP for this component
        selected = solve_component_milp(
            coverage, blocked_in_comp, effective_candidates, tradeoff_lambda, verbose
        )

        if not selected:
            if verbose:
                print(f"  MILP returned no selections")
            continue

        # Add selected PTRs to model
        for i, sel in enumerate(selected):
            rid = add_transport(model, sel, len(all_selected) + i)
            if rid:
                sel["reaction_id"] = rid
                sel["target_component"] = comp_id
                all_selected.append(sel)
                used_candidates.add((sel["met1"], sel["met2"]))

        # Track results
        component_results.append(
            {
                "component_id": comp_id,
                "component_size": comp_size,
                "reactions_in_comp": len(rxns_in_comp),
                "blocked_before": len(blocked_in_comp),
                "candidates_available": len(comp_candidates),
                "ptrs_selected": len(selected),
            }
        )

        if verbose:
            # Check improvement
            new_blocked = find_blocked_reactions_in_set(model, rxns_in_comp)
            print(
                f"  After adding PTRs: {len(new_blocked)} still blocked (was {len(blocked_in_comp)})"
            )

    total_elapsed = time.time() - total_start

    # Write output CSV
    out_csv = os.path.join(out_dir, "selected_phase3_component_milp.csv")
    with open(out_csv, "w", newline="") as f:
        fieldnames = ["reaction_id", "met1", "met2", "type", "base", "target_component"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for sel in all_selected:
            writer.writerow(sel)

    # Write component results summary
    results_csv = os.path.join(out_dir, "phase3_component_results.csv")
    with open(results_csv, "w", newline="") as f:
        fieldnames = [
            "component_id",
            "component_size",
            "reactions_in_comp",
            "blocked_before",
            "candidates_available",
            "ptrs_selected",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in component_results:
            writer.writerow(r)

    # Save model
    if phase3_model_out is None:
        phase3_model_out = os.path.normpath(
            os.path.join(
                base_dir,
                "..",
                "models",
                "base",
                "THG-beta-batch_251106_phase3_component_milp.json",
            )
        )
    save_json_model(model, phase3_model_out)

    if verbose:
        print(f"\n{'='*60}")
        print(f"FINAL SUMMARY")
        print(f"{'='*60}")
        print(f"Total time: {total_elapsed:.1f}s")
        print(f"Components processed: {len(component_results)}")
        print(f"Total PTRs selected: {len(all_selected)}")
        print(f"Output CSV: {out_csv}")
        print(f"Component results: {results_csv}")
        print(f"Output model: {phase3_model_out}")

    return out_csv, phase3_model_out, len(all_selected)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Component-wise MILP Phase-3 optimizer"
    )
    parser.add_argument("--candidates", default=None, help="Phase-1 candidates CSV")
    parser.add_argument("--model", default=None, help="Starting model JSON")
    parser.add_argument("--components", default=None, help="Components summary CSV")
    parser.add_argument("--out", default=None, help="Output directory")
    parser.add_argument(
        "--lambda",
        dest="tradeoff_lambda",
        type=float,
        default=0.01,
        help="Tradeoff weight (higher = fewer PTRs)",
    )
    parser.add_argument(
        "--min-size", type=int, default=4, help="Minimum component size to process"
    )
    parser.add_argument(
        "--max-components",
        type=int,
        default=None,
        help="Max number of components to process",
    )
    args = parser.parse_args()

    cand_csv = args.candidates or os.path.join(
        os.path.dirname(__file__), "files", "candidates_all.csv"
    )

    if not os.path.exists(cand_csv):
        print("Candidates CSV not found:", cand_csv)
    else:
        run_phase3_component_wise(
            cand_csv,
            starting_model_json=args.model,
            components_csv=args.components,
            out_dir=args.out,
            tradeoff_lambda=args.tradeoff_lambda,
            min_component_size=args.min_size,
            max_components=args.max_components,
        )
