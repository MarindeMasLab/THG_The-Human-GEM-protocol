"""Phase-3 Tiered Component-wise MILP optimizer.

This approach solves the gap-filling problem using a lexicographic/tiered strategy:
1. First, use only Type A candidates (both endpoints are dead-ends)
2. Then, use Type B candidates (one endpoint is dead-end) with remaining blocked reactions
3. Finally, use Type C candidates (neither endpoint is dead-end) for any remaining

This avoids arbitrary weight choices while respecting the biological intuition that
Type A candidates are most likely to represent missing transporters.

For each tier, we solve component-by-component to keep the problem tractable.
"""

import csv
import os
import time
import pickle
import networkx as nx
from collections import defaultdict
from cobra.io import load_json_model, save_json_model
from cobra import Reaction, Model


def load_candidates(path):
    """Load candidate connectors from Phase-1 CSV."""
    rows = []
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        for r in reader:
            try:
                r['component1'] = int(r.get('component1')) if r.get('component1') not in (None, '') else None
                r['component2'] = int(r.get('component2')) if r.get('component2') not in (None, '') else None
            except Exception:
                pass
            rows.append(r)
    return rows


def load_components_summary(path):
    """Load component summary CSV."""
    components = {}
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        for r in reader:
            cid = int(r['component_id'])
            size = int(r['num_nodes'])
            components[cid] = size
    return components


def build_metabolite_to_component_map(model):
    """Build mapping from metabolite IDs to their component IDs."""
    G = nx.Graph()
    for rxn in model.reactions:
        rid = rxn.id
        G.add_node(rid, bipartite=0)
        for met in rxn.metabolites:
            mid = met.id
            G.add_node(mid, bipartite=1)
            G.add_edge(rid, mid)
    
    components = list(nx.connected_components(G))
    components.sort(key=len, reverse=True)  # largest = main component = ID 1
    
    met_to_comp = {}
    for comp_id, comp_nodes in enumerate(components, start=1):
        for node in comp_nodes:
            if node.startswith('MAM'):
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


def filter_candidates_by_type(candidates, candidate_type):
    """Filter candidates by type (A, B, or C)."""
    return [c for c in candidates if c.get('type', '').upper() == candidate_type.upper()]


def get_candidates_for_component(candidates, component_id, used_candidates=None):
    """Get candidates that connect TO the given component."""
    if used_candidates is None:
        used_candidates = set()
    
    relevant = []
    for c in candidates:
        key = (c['met1'], c['met2'])
        if key in used_candidates:
            continue
        
        c1 = c.get('component1')
        c2 = c.get('component2')
        
        # Skip if both endpoints are in main component (1)
        if c1 == 1 and c2 == 1:
            continue
        
        # Include if either endpoint is in target component
        if c1 == component_id or c2 == component_id:
            relevant.append(c)
    
    return relevant


def find_blocked_reactions_in_set(model, reaction_ids, zero_cutoff=None):
    """Find which reactions in the given set are blocked."""
    if zero_cutoff is None:
        zero_cutoff = model.tolerance * 10
    
    blocked = set()
    
    for rid in reaction_ids:
        try:
            rxn = model.reactions.get_by_id(rid)
        except KeyError:
            continue
        
        can_carry_flux = False
        
        if rxn.upper_bound > zero_cutoff:
            with model:
                rxn.lower_bound = zero_cutoff
                try:
                    sol = model.optimize()
                    if sol.status == 'optimal':
                        can_carry_flux = True
                except Exception:
                    pass
        
        if not can_carry_flux and rxn.lower_bound < -zero_cutoff:
            with model:
                rxn.upper_bound = -zero_cutoff
                try:
                    sol = model.optimize()
                    if sol.status == 'optimal':
                        can_carry_flux = True
                except Exception:
                    pass
        
        if not can_carry_flux:
            blocked.add(rid)
    
    return blocked


def add_transport(model, sel, idx, prefix='TIER'):
    """Add a reversible transport reaction to the model."""
    met1 = sel['met1']
    met2 = sel['met2']
    base = sel.get('base', '')
    s1 = sel.get('suffix1', '')
    s2 = sel.get('suffix2', '')
    ctype = sel.get('type', 'X')
    raw_id = f"{prefix}{ctype}_TRANS_{base}_{s1}_{s2}_{idx}"
    rid = ''.join(ch if (ch.isalnum() or ch == '_') else '_' for ch in raw_id)
    
    try:
        m1 = model.metabolites.get_by_id(met1)
        m2 = model.metabolites.get_by_id(met2)
    except KeyError:
        return None
    
    rxn = Reaction(rid)
    rxn.name = f"Tiered MILP transport {met1} -> {met2} (Type {ctype})"
    rxn.add_metabolites({m1: -1.0, m2: 1.0})
    rxn.lower_bound = -1000.0
    rxn.upper_bound = 1000.0
    rxn.annotation = {'phase3_tiered_milp': 'true', 'candidate_type': ctype}
    model.add_reactions([rxn])
    return rid


def test_candidate_unblocks(model, candidate, blocked_rxns, zero_cutoff=None):
    """Test which blocked reactions a candidate PTR can unblock."""
    if zero_cutoff is None:
        zero_cutoff = model.tolerance * 10
    
    met1 = candidate['met1']
    met2 = candidate['met2']
    
    try:
        m1 = model.metabolites.get_by_id(met1)
        m2 = model.metabolites.get_by_id(met2)
    except KeyError:
        return set()
    
    with model:
        rxn = Reaction('_TEST_PTR_')
        rxn.add_metabolites({m1: -1.0, m2: 1.0})
        rxn.lower_bound = -1000.0
        rxn.upper_bound = 1000.0
        model.add_reactions([rxn])
        
        unblocked = set()
        for rid in blocked_rxns:
            try:
                r = model.reactions.get_by_id(rid)
            except KeyError:
                continue
            
            can_carry = False
            
            if r.upper_bound > zero_cutoff:
                with model:
                    r.lower_bound = zero_cutoff
                    try:
                        sol = model.optimize()
                        if sol.status == 'optimal':
                            can_carry = True
                    except Exception:
                        pass
            
            if not can_carry and r.lower_bound < -zero_cutoff:
                with model:
                    r.upper_bound = -zero_cutoff
                    try:
                        sol = model.optimize()
                        if sol.status == 'optimal':
                            can_carry = True
                    except Exception:
                        pass
            
            if can_carry:
                unblocked.add(rid)
    
    return unblocked


def compute_coverage_for_component(model, candidates, blocked_rxns, verbose=True):
    """Compute coverage matrix for candidates and blocked reactions."""
    coverage = {}
    total = len(candidates)
    n_blocked = len(blocked_rxns)
    
    if verbose:
        print(f"      Coverage: {total} candidates × {n_blocked} blocked")
    
    start_time = time.time()
    
    for i, cand in enumerate(candidates):
        if verbose and (i + 1) % 50 == 0:
            elapsed = time.time() - start_time
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            eta = (total - i - 1) / rate if rate > 0 else 0
            print(f"        {i+1}/{total} ({rate:.1f}/sec, ETA: {eta:.0f}s)")
        
        key = (cand['met1'], cand['met2'])
        unblocked = test_candidate_unblocks(model, cand, blocked_rxns)
        coverage[key] = unblocked
    
    if verbose:
        elapsed = time.time() - start_time
        non_empty = sum(1 for v in coverage.values() if v)
        print(f"      Done in {elapsed:.1f}s: {non_empty}/{total} useful")
    
    return coverage


def solve_milp(coverage, blocked_rxns, candidates, tradeoff_lambda=0.01, verbose=True):
    """Solve MILP for coverage optimization."""
    # Try solvers in order
    try:
        import gurobipy
        return _solve_gurobi(coverage, blocked_rxns, candidates, tradeoff_lambda, verbose)
    except ImportError:
        pass
    
    try:
        import pulp
        return _solve_pulp(coverage, blocked_rxns, candidates, tradeoff_lambda, verbose)
    except ImportError:
        pass
    
    try:
        from scipy.optimize import milp
        return _solve_scipy(coverage, blocked_rxns, candidates, tradeoff_lambda, verbose)
    except ImportError:
        pass
    
    raise ImportError("No MILP solver available")


def _solve_gurobi(coverage, blocked_rxns, candidates, tradeoff_lambda, verbose):
    import gurobipy as gp
    from gurobipy import GRB
    
    blocked_list = list(blocked_rxns)
    blocked_idx = {b: i for i, b in enumerate(blocked_list)}
    cand_idx = {(c['met1'], c['met2']): i for i, c in enumerate(candidates)}
    
    n_cands = len(candidates)
    n_blocked = len(blocked_list)
    
    m = gp.Model("coverage")
    m.Params.OutputFlag = 1 if verbose else 0
    
    y = m.addVars(n_cands, vtype=GRB.BINARY, name="y")
    z = m.addVars(n_blocked, vtype=GRB.BINARY, name="z")
    
    for b_id, b_i in blocked_idx.items():
        covering = [cand_idx[key] for key, unblocked in coverage.items() 
                   if b_id in unblocked and key in cand_idx]
        if covering:
            m.addConstr(z[b_i] <= gp.quicksum(y[p] for p in covering))
        else:
            m.addConstr(z[b_i] == 0)
    
    m.setObjective(
        gp.quicksum(z[b] for b in range(n_blocked)) - tradeoff_lambda * gp.quicksum(y[p] for p in range(n_cands)),
        GRB.MAXIMIZE
    )
    
    m.optimize()
    
    if m.Status != GRB.OPTIMAL:
        return []
    
    selected = [c for i, c in enumerate(candidates) if y[i].X > 0.5]
    
    if verbose:
        n_unblocked = sum(1 for b in range(n_blocked) if z[b].X > 0.5)
        print(f"      MILP: {len(selected)} PTRs → {n_unblocked} unblocked")
    
    return selected


def _solve_pulp(coverage, blocked_rxns, candidates, tradeoff_lambda, verbose):
    import pulp
    
    blocked_list = list(blocked_rxns)
    blocked_idx = {b: i for i, b in enumerate(blocked_list)}
    cand_idx = {(c['met1'], c['met2']): i for i, c in enumerate(candidates)}
    
    n_cands = len(candidates)
    n_blocked = len(blocked_list)
    
    prob = pulp.LpProblem("coverage", pulp.LpMaximize)
    
    y = [pulp.LpVariable(f"y_{i}", cat='Binary') for i in range(n_cands)]
    z = [pulp.LpVariable(f"z_{i}", cat='Binary') for i in range(n_blocked)]
    
    for b_id, b_i in blocked_idx.items():
        covering = [cand_idx[key] for key, unblocked in coverage.items() 
                   if b_id in unblocked and key in cand_idx]
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
        print(f"      MILP: {len(selected)} PTRs → {n_unblocked} unblocked")
    
    return selected


def _solve_scipy(coverage, blocked_rxns, candidates, tradeoff_lambda, verbose):
    import numpy as np
    from scipy.optimize import milp, LinearConstraint, Bounds
    
    blocked_list = list(blocked_rxns)
    blocked_idx = {b: i for i, b in enumerate(blocked_list)}
    cand_idx = {(c['met1'], c['met2']): i for i, c in enumerate(candidates)}
    
    n_cands = len(candidates)
    n_blocked = len(blocked_list)
    n_vars = n_cands + n_blocked
    
    c = np.zeros(n_vars)
    c[:n_cands] = tradeoff_lambda
    c[n_cands:] = -1.0
    
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
        print(f"      MILP: {len(selected)} PTRs → {n_unblocked} unblocked")
    
    return selected


def process_tier_for_component(model, tier_candidates, blocked_rxns, 
                                tradeoff_lambda, verbose=True):
    """Process a single tier (A, B, or C) for a single component."""
    if not tier_candidates or not blocked_rxns:
        return [], blocked_rxns
    
    # Compute coverage
    coverage = compute_coverage_for_component(model, tier_candidates, blocked_rxns, verbose)
    
    # Filter to candidates with any coverage
    effective = [c for c in tier_candidates if coverage.get((c['met1'], c['met2']), set())]
    
    if not effective:
        if verbose:
            print(f"      No effective candidates")
        return [], blocked_rxns
    
    if verbose:
        print(f"      Effective: {len(effective)} candidates")
    
    # Solve MILP
    selected = solve_milp(coverage, blocked_rxns, effective, tradeoff_lambda, verbose)
    
    return selected, coverage


def run_phase3_tiered(
    candidates_csv,
    starting_model_json=None,
    components_csv=None,
    out_dir=None,
    phase3_model_out=None,
    tradeoff_lambda=0.01,
    min_component_size=4,
    max_components=None,
    verbose=True
):
    """
    Run Phase-3 with tiered optimization: Type A → Type B → Type C.
    
    For each isolated component:
      1. Find blocked reactions
      2. Try to unblock with Type A candidates (MILP)
      3. Try remaining blocked with Type B candidates (MILP)  
      4. Try remaining blocked with Type C candidates (MILP)
    """
    # Setup paths
    base_dir = os.path.dirname(__file__)
    if out_dir is None:
        out_dir = os.path.join(base_dir, 'files')
    os.makedirs(out_dir, exist_ok=True)
    
    if starting_model_json is None:
        starting_model_json = os.path.normpath(os.path.join(
            base_dir, '..', 'models', 'base', 
            'THG-beta-batch_251106_phase2_minimal_connected.json'))
    
    if components_csv is None:
        components_csv = os.path.join(out_dir, 'components_summary.csv')
    
    # Load data
    all_candidates = load_candidates(candidates_csv)
    components = load_components_summary(components_csv)
    model = load_json_model(starting_model_json)
    
    # Split candidates by type
    type_a = filter_candidates_by_type(all_candidates, 'A')
    type_b = filter_candidates_by_type(all_candidates, 'B')
    type_c = filter_candidates_by_type(all_candidates, 'C')
    
    if verbose:
        print(f"Loaded {len(all_candidates)} candidates:")
        print(f"  Type A (both dead-end): {len(type_a)}")
        print(f"  Type B (one dead-end):  {len(type_b)}")
        print(f"  Type C (neither):       {len(type_c)}")
        print(f"Model: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")
    
    # Build metabolite-to-component mapping
    if verbose:
        print("Building metabolite-to-component mapping...")
    met_to_comp = build_metabolite_to_component_map(model)
    
    # Get isolated components
    isolated_comps = [(cid, size) for cid, size in components.items() 
                      if cid > 1 and size >= min_component_size]
    isolated_comps.sort(key=lambda x: x[1], reverse=True)
    
    if max_components:
        isolated_comps = isolated_comps[:max_components]
    
    if verbose:
        print(f"\nProcessing {len(isolated_comps)} components (size >= {min_component_size})")
        print(f"Main component (ID=1, {components.get(1, 0)} nodes) EXCLUDED\n")
    
    # Track state
    used_candidates = set()
    all_selected = []
    tier_stats = {'A': 0, 'B': 0, 'C': 0}
    component_results = []
    
    total_start = time.time()
    
    for comp_idx, (comp_id, comp_size) in enumerate(isolated_comps):
        if verbose:
            print(f"\n{'='*60}")
            print(f"Component {comp_id} ({comp_idx+1}/{len(isolated_comps)}): {comp_size} nodes")
            print(f"{'='*60}")
        
        # Get reactions in this component
        rxns_in_comp = get_reactions_in_component(model, comp_id, met_to_comp)
        if verbose:
            print(f"  Reactions in component: {len(rxns_in_comp)}")
        
        if not rxns_in_comp:
            continue
        
        # Find blocked reactions
        blocked = find_blocked_reactions_in_set(model, rxns_in_comp)
        initial_blocked = len(blocked)
        
        if verbose:
            print(f"  Initially blocked: {len(blocked)}")
        
        if not blocked:
            continue
        
        comp_selected = []
        
        # Process each tier
        for tier_name, tier_candidates in [('A', type_a), ('B', type_b), ('C', type_c)]:
            if not blocked:
                break
            
            # Get candidates for this component and tier
            comp_tier_cands = get_candidates_for_component(tier_candidates, comp_id, used_candidates)
            
            if verbose:
                print(f"\n  --- Tier {tier_name}: {len(comp_tier_cands)} candidates available ---")
            
            if not comp_tier_cands:
                if verbose:
                    print(f"    No candidates available")
                continue
            
            # Process this tier
            selected, coverage = process_tier_for_component(
                model, comp_tier_cands, blocked, tradeoff_lambda, verbose)
            
            if not selected:
                continue
            
            # Add selected PTRs to model
            for sel in selected:
                rid = add_transport(model, sel, len(all_selected))
                if rid:
                    sel['reaction_id'] = rid
                    sel['target_component'] = comp_id
                    sel['tier'] = tier_name
                    all_selected.append(sel)
                    comp_selected.append(sel)
                    used_candidates.add((sel['met1'], sel['met2']))
                    tier_stats[tier_name] += 1
            
            # Recompute blocked
            new_blocked = find_blocked_reactions_in_set(model, rxns_in_comp)
            unblocked_count = len(blocked) - len(new_blocked)
            
            if verbose:
                print(f"    Added {len(selected)} PTRs, unblocked {unblocked_count} reactions")
                print(f"    Remaining blocked: {len(new_blocked)}")
            
            blocked = new_blocked
        
        # Record component results
        final_blocked = len(find_blocked_reactions_in_set(model, rxns_in_comp))
        component_results.append({
            'component_id': comp_id,
            'component_size': comp_size,
            'reactions_in_comp': len(rxns_in_comp),
            'blocked_before': initial_blocked,
            'blocked_after': final_blocked,
            'ptrs_selected': len(comp_selected),
            'type_a_used': sum(1 for s in comp_selected if s.get('tier') == 'A'),
            'type_b_used': sum(1 for s in comp_selected if s.get('tier') == 'B'),
            'type_c_used': sum(1 for s in comp_selected if s.get('tier') == 'C'),
        })
    
    total_elapsed = time.time() - total_start
    
    # Write outputs
    out_csv = os.path.join(out_dir, 'selected_phase3_tiered_milp.csv')
    with open(out_csv, 'w', newline='') as f:
        fieldnames = ['reaction_id', 'met1', 'met2', 'type', 'tier', 'base', 'target_component']
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        for sel in all_selected:
            writer.writerow(sel)
    
    results_csv = os.path.join(out_dir, 'phase3_tiered_results.csv')
    with open(results_csv, 'w', newline='') as f:
        fieldnames = ['component_id', 'component_size', 'reactions_in_comp', 
                     'blocked_before', 'blocked_after', 'ptrs_selected',
                     'type_a_used', 'type_b_used', 'type_c_used']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in component_results:
            writer.writerow(r)
    
    if phase3_model_out is None:
        phase3_model_out = os.path.normpath(os.path.join(
            base_dir, '..', 'models', 'base',
            'THG-beta-batch_251106_phase3_tiered_milp.json'))
    save_json_model(model, phase3_model_out)
    
    if verbose:
        print(f"\n{'='*60}")
        print(f"FINAL SUMMARY")
        print(f"{'='*60}")
        print(f"Total time: {total_elapsed:.1f}s ({total_elapsed/60:.1f} min)")
        print(f"Components processed: {len(component_results)}")
        print(f"Total PTRs selected: {len(all_selected)}")
        print(f"  Type A: {tier_stats['A']}")
        print(f"  Type B: {tier_stats['B']}")
        print(f"  Type C: {tier_stats['C']}")
        print(f"Output CSV: {out_csv}")
        print(f"Component results: {results_csv}")
        print(f"Output model: {phase3_model_out}")
    
    return out_csv, phase3_model_out, len(all_selected)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Tiered MILP Phase-3 optimizer (A → B → C)')
    parser.add_argument('--candidates', default=None, help='Phase-1 candidates CSV')
    parser.add_argument('--model', default=None, help='Starting model JSON')
    parser.add_argument('--components', default=None, help='Components summary CSV')
    parser.add_argument('--out', default=None, help='Output directory')
    parser.add_argument('--lambda', dest='tradeoff_lambda', type=float, default=0.01,
                       help='Tradeoff weight (higher = fewer PTRs)')
    parser.add_argument('--min-size', type=int, default=4,
                       help='Minimum component size to process')
    parser.add_argument('--max-components', type=int, default=None,
                       help='Max number of components to process')
    args = parser.parse_args()
    
    cand_csv = args.candidates or os.path.join(os.path.dirname(__file__), 'files', 'candidates_all.csv')
    
    if not os.path.exists(cand_csv):
        print('Candidates CSV not found:', cand_csv)
    else:
        run_phase3_tiered(
            cand_csv,
            starting_model_json=args.model,
            components_csv=args.components,
            out_dir=args.out,
            tradeoff_lambda=args.tradeoff_lambda,
            min_component_size=args.min_size,
            max_components=args.max_components
        )
