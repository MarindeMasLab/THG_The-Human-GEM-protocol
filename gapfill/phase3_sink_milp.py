"""Phase-3 Component-wise MILP optimizer with temporary sink/source approach.

This approach adds temporary sink/source reactions for dead-end metabolites when
testing if a candidate PTR can enable flux through a blocked reaction.

Key insight: A blocked reaction may be blocked because metabolites have nowhere
to go, not because of missing transport. By adding temporary exchanges, we test
whether the PTR enables the pathway connectivity, assuming exchanges exist.

Temporary reaction rules:
- Dead-end substrate (consumed but not produced): add source (→ met)
- Dead-end product (produced but not consumed): add sink (met →)
- Neither: add reversible exchange (↔ met)
"""

import csv
import os
import time
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
    components.sort(key=len, reverse=True)
    
    met_to_comp = {}
    rxn_to_comp = {}
    for comp_id, comp_nodes in enumerate(components, start=1):
        for node in comp_nodes:
            if node.startswith('MAM'):
                met_to_comp[node] = comp_id
            else:
                rxn_to_comp[node] = comp_id
    
    return met_to_comp, rxn_to_comp, components


def get_deadend_info(model):
    """
    Identify dead-end metabolites and their types.
    
    Returns dict: met_id -> {'type': 'substrate'|'product'|'both'|'none', 
                             'producing_rxns': [...], 'consuming_rxns': [...]}
    """
    deadend_info = {}
    
    for met in model.metabolites:
        producing = []  # reactions that produce this metabolite
        consuming = []  # reactions that consume this metabolite
        
        for rxn in met.reactions:
            coeff = rxn.metabolites[met]
            lb, ub = rxn.lower_bound, rxn.upper_bound
            
            # Check if reaction can produce this metabolite
            can_produce = (coeff > 0 and ub > 0) or (coeff < 0 and lb < 0)
            # Check if reaction can consume this metabolite
            can_consume = (coeff < 0 and ub > 0) or (coeff > 0 and lb < 0)
            
            if can_produce:
                producing.append(rxn.id)
            if can_consume:
                consuming.append(rxn.id)
        
        # Determine dead-end type
        if not producing and consuming:
            de_type = 'substrate'  # consumed but not produced
        elif producing and not consuming:
            de_type = 'product'  # produced but not consumed
        elif not producing and not consuming:
            de_type = 'isolated'  # completely isolated
        else:
            de_type = 'none'  # has both production and consumption
        
        if de_type != 'none':
            deadend_info[met.id] = {
                'type': de_type,
                'producing_rxns': producing,
                'consuming_rxns': consuming
            }
    
    return deadend_info


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


def add_temp_sink_source(model, met_id, deadend_info):
    """
    Add temporary sink/source reaction for a dead-end metabolite.
    
    Returns the reaction ID added, or None if not a dead-end.
    """
    if met_id not in deadend_info:
        return None
    
    info = deadend_info[met_id]
    de_type = info['type']
    
    try:
        met = model.metabolites.get_by_id(met_id)
    except KeyError:
        return None
    
    rxn_id = f'_TEMP_EXCH_{met_id}_'
    rxn = Reaction(rxn_id)
    rxn.name = f'Temporary exchange for {met_id}'
    
    if de_type == 'substrate':
        # Consumed but not produced -> add source (→ met)
        rxn.add_metabolites({met: 1.0})
        rxn.lower_bound = 0
        rxn.upper_bound = 1000.0
    elif de_type == 'product':
        # Produced but not consumed -> add sink (met →)
        rxn.add_metabolites({met: -1.0})
        rxn.lower_bound = 0
        rxn.upper_bound = 1000.0
    elif de_type == 'isolated':
        # Completely isolated -> add reversible exchange
        rxn.add_metabolites({met: -1.0})
        rxn.lower_bound = -1000.0
        rxn.upper_bound = 1000.0
    else:
        return None
    
    model.add_reactions([rxn])
    return rxn_id


def test_candidate_with_temp_sinks(model, candidate, blocked_rxn_id, deadend_info, zero_cutoff=None):
    """
    Test if PTR + temporary sinks can enable flux through a blocked reaction.
    
    1. Add the candidate PTR
    2. Add temporary sink/source for dead-end metabolites in the blocked reaction
    3. Test if blocked reaction can now carry flux
    
    Returns True if the reaction can carry flux with the PTR + temp sinks.
    """
    if zero_cutoff is None:
        zero_cutoff = model.tolerance * 10
    
    met1 = candidate['met1']
    met2 = candidate['met2']
    
    try:
        m1 = model.metabolites.get_by_id(met1)
        m2 = model.metabolites.get_by_id(met2)
        blocked_rxn = model.reactions.get_by_id(blocked_rxn_id)
    except KeyError:
        return False
    
    with model:
        # Add candidate PTR
        ptr = Reaction('_TEST_PTR_')
        ptr.add_metabolites({m1: -1.0, m2: 1.0})
        ptr.lower_bound = -1000.0
        ptr.upper_bound = 1000.0
        model.add_reactions([ptr])
        
        # Add temporary sinks/sources for dead-end metabolites in the blocked reaction
        temp_rxns = []
        for met in blocked_rxn.metabolites:
            temp_id = add_temp_sink_source(model, met.id, deadend_info)
            if temp_id:
                temp_rxns.append(temp_id)
        
        # Also add temp sinks for the PTR metabolites if they're dead-ends
        for mid in [met1, met2]:
            if mid not in [m.id for m in blocked_rxn.metabolites]:
                temp_id = add_temp_sink_source(model, mid, deadend_info)
                if temp_id:
                    temp_rxns.append(temp_id)
        
        # Test if blocked reaction can carry flux
        can_carry = False
        
        if blocked_rxn.upper_bound > zero_cutoff:
            with model:
                blocked_rxn.lower_bound = zero_cutoff
                try:
                    sol = model.optimize()
                    if sol.status == 'optimal':
                        can_carry = True
                except Exception:
                    pass
        
        if not can_carry and blocked_rxn.lower_bound < -zero_cutoff:
            with model:
                blocked_rxn.upper_bound = -zero_cutoff
                try:
                    sol = model.optimize()
                    if sol.status == 'optimal':
                        can_carry = True
                except Exception:
                    pass
        
        return can_carry


def compute_coverage_with_temp_sinks(model, candidates, blocked_rxns, deadend_info, verbose=True):
    """
    Compute coverage matrix using temporary sink/source approach.
    
    Returns: dict mapping (met1, met2) -> set of unblocked reaction IDs
    """
    coverage = {}
    total = len(candidates)
    n_blocked = len(blocked_rxns)
    
    if verbose:
        print(f"      Coverage (with temp sinks): {total} candidates × {n_blocked} blocked")
    
    start_time = time.time()
    
    for i, cand in enumerate(candidates):
        if verbose and (i + 1) % 10 == 0:
            elapsed = time.time() - start_time
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            eta = (total - i - 1) / rate if rate > 0 else 0
            print(f"        {i+1}/{total} ({rate:.2f}/sec, ETA: {eta:.0f}s)")
        
        key = (cand['met1'], cand['met2'])
        unblocked = set()
        
        for rid in blocked_rxns:
            if test_candidate_with_temp_sinks(model, cand, rid, deadend_info):
                unblocked.add(rid)
        
        coverage[key] = unblocked
    
    if verbose:
        elapsed = time.time() - start_time
        non_empty = sum(1 for v in coverage.values() if v)
        total_pairs = sum(len(v) for v in coverage.values())
        print(f"      Done in {elapsed:.1f}s: {non_empty}/{total} candidates useful, {total_pairs} total unblocks")
    
    return coverage


def add_transport(model, sel, idx, prefix='SINK'):
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
    rxn.name = f"Temp-sink MILP transport {met1} -> {met2} (Type {ctype})"
    rxn.add_metabolites({m1: -1.0, m2: 1.0})
    rxn.lower_bound = -1000.0
    rxn.upper_bound = 1000.0
    rxn.annotation = {'phase3_sink_milp': 'true', 'candidate_type': ctype}
    model.add_reactions([rxn])
    return rid


def solve_milp(coverage, blocked_rxns, candidates, tradeoff_lambda=0.01, verbose=True):
    """Solve MILP for coverage optimization."""
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
    
    raise ImportError("No MILP solver available (need pulp or scipy)")


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
    
    selected = [c for i, c in enumerate(candidates) if y[i].value() and y[i].value() > 0.5]
    
    if verbose:
        n_unblocked = sum(1 for b in range(n_blocked) if z[b].value() and z[b].value() > 0.5)
        print(f"      MILP: {len(selected)} PTRs → {n_unblocked} reactions can be unblocked")
    
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
        print(f"      MILP: {len(selected)} PTRs → {n_unblocked} reactions can be unblocked")
    
    return selected


def get_candidates_for_component(candidates, component_id, rxn_to_comp, used_candidates=None):
    """Get candidates relevant to a component (either endpoint connects to it)."""
    if used_candidates is None:
        used_candidates = set()
    
    relevant = []
    for c in candidates:
        key = (c['met1'], c['met2'])
        if key in used_candidates:
            continue
        
        c1 = c.get('component1')
        c2 = c.get('component2')
        
        # Include if either endpoint connects to target component
        if c1 == component_id or c2 == component_id:
            relevant.append(c)
    
    return relevant


def run_test_on_component(
    candidates_csv,
    starting_model_json,
    target_component_index=2,  # 1-indexed, 1=main, 2=largest isolated, etc.
    out_dir=None,
    tradeoff_lambda=0.01,
    verbose=True
):
    """
    Run a test on a specific component using temp-sink approach.
    
    Args:
        target_component_index: 1=main component, 2=largest isolated, 3=second largest, etc.
    """
    base_dir = os.path.dirname(__file__)
    if out_dir is None:
        out_dir = os.path.join(base_dir, 'files')
    os.makedirs(out_dir, exist_ok=True)
    
    # Load data
    candidates = load_candidates(candidates_csv)
    model = load_json_model(starting_model_json)
    
    if verbose:
        print(f"Loaded {len(candidates)} candidates")
        print(f"Model: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")
    
    # Build component mapping
    if verbose:
        print("Building component mapping...")
    met_to_comp, rxn_to_comp, components = build_metabolite_to_component_map(model)
    
    if verbose:
        print(f"Found {len(components)} components:")
        for i, comp in enumerate(components[:10], 1):
            rxns = sum(1 for n in comp if not n.startswith('MAM'))
            mets = sum(1 for n in comp if n.startswith('MAM'))
            print(f"  Component {i}: {len(comp)} nodes ({rxns} rxns, {mets} mets)")
    
    # Get target component
    if target_component_index < 1 or target_component_index > len(components):
        print(f"Invalid component index {target_component_index}")
        return None
    
    target_comp = components[target_component_index - 1]
    rxns_in_comp = [n for n in target_comp if not n.startswith('MAM')]
    mets_in_comp = [n for n in target_comp if n.startswith('MAM')]
    
    if verbose:
        print(f"\n{'='*60}")
        print(f"Testing Component {target_component_index}: {len(target_comp)} nodes")
        print(f"  Reactions: {len(rxns_in_comp)}")
        print(f"  Metabolites: {len(mets_in_comp)}")
        print(f"{'='*60}")
    
    # Get dead-end info
    if verbose:
        print("\nAnalyzing dead-end metabolites...")
    deadend_info = get_deadend_info(model)
    
    deadends_in_comp = {m: deadend_info[m] for m in mets_in_comp if m in deadend_info}
    if verbose:
        by_type = defaultdict(int)
        for info in deadends_in_comp.values():
            by_type[info['type']] += 1
        print(f"  Dead-ends in component: {len(deadends_in_comp)}")
        for t, c in sorted(by_type.items()):
            print(f"    {t}: {c}")
    
    # Find blocked reactions in component
    if verbose:
        print("\nFinding blocked reactions...")
    blocked = find_blocked_reactions_in_set(model, rxns_in_comp)
    
    if verbose:
        print(f"  Blocked: {len(blocked)} / {len(rxns_in_comp)}")
        if blocked:
            print(f"  Blocked reaction IDs: {list(blocked)[:10]}{'...' if len(blocked) > 10 else ''}")
    
    if not blocked:
        print("No blocked reactions in this component!")
        return None
    
    # Get candidates for this component
    comp_candidates = get_candidates_for_component(candidates, target_component_index, rxn_to_comp)
    
    if verbose:
        print(f"\nCandidates for this component: {len(comp_candidates)}")
        by_type = defaultdict(int)
        for c in comp_candidates:
            by_type[c.get('type', '?')] += 1
        for t, cnt in sorted(by_type.items()):
            print(f"  Type {t}: {cnt}")
    
    if not comp_candidates:
        print("No candidates available for this component!")
        return None
    
    # Compute coverage with temp sinks
    if verbose:
        print("\nComputing coverage matrix (with temp sinks)...")
    
    coverage = compute_coverage_with_temp_sinks(model, comp_candidates, blocked, deadend_info, verbose)
    
    # Filter to effective candidates
    effective = [c for c in comp_candidates if coverage.get((c['met1'], c['met2']), set())]
    
    if verbose:
        print(f"\nEffective candidates (can unblock ≥1 reaction): {len(effective)}")
    
    if not effective:
        print("No candidates can unblock any reactions even with temp sinks!")
        return {'blocked': len(blocked), 'effective_candidates': 0, 'selected': 0}
    
    # Solve MILP
    if verbose:
        print("\nSolving MILP...")
    
    selected = solve_milp(coverage, blocked, effective, tradeoff_lambda, verbose)
    
    # Report results
    result = {
        'component_index': target_component_index,
        'component_size': len(target_comp),
        'reactions_in_comp': len(rxns_in_comp),
        'blocked_before': len(blocked),
        'candidates_available': len(comp_candidates),
        'effective_candidates': len(effective),
        'selected': len(selected),
        'selected_details': selected
    }
    
    if verbose:
        print(f"\n{'='*60}")
        print("RESULTS")
        print(f"{'='*60}")
        print(f"Component {target_component_index}: {len(target_comp)} nodes, {len(rxns_in_comp)} reactions")
        print(f"Blocked reactions: {len(blocked)}")
        print(f"Effective candidates: {len(effective)}")
        print(f"Selected PTRs: {len(selected)}")
        if selected:
            print("Selected PTR details:")
            for s in selected[:10]:
                print(f"  {s['met1']} <-> {s['met2']} (Type {s.get('type', '?')})")
    
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Test temp-sink MILP approach on a component')
    parser.add_argument('--candidates', default=None, help='Phase-1 candidates CSV')
    parser.add_argument('--model', default=None, help='Starting model JSON')
    parser.add_argument('--component', type=int, default=3, 
                       help='Component index (1=main, 2=largest isolated, etc.)')
    parser.add_argument('--lambda', dest='tradeoff_lambda', type=float, default=0.01)
    args = parser.parse_args()
    
    base_dir = os.path.dirname(__file__)
    cand_csv = args.candidates or os.path.join(base_dir, 'files', 'candidates_all.csv')
    model_json = args.model or os.path.normpath(os.path.join(
        base_dir, '..', 'models', 'base', 'THG-beta-batch_251106_phase2_minimal_connected.json'))
    
    if not os.path.exists(cand_csv):
        print('Candidates CSV not found:', cand_csv)
    elif not os.path.exists(model_json):
        print('Model not found:', model_json)
    else:
        run_test_on_component(
            cand_csv,
            model_json,
            target_component_index=args.component,
            tradeoff_lambda=args.tradeoff_lambda
        )
