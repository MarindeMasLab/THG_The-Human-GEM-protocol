"""Phase-3 Sink MILP - Uses ORIGINAL component assignments from Phase 1.

This version preserves the original component structure from the unconnected model,
so we can track which reactions belonged to each isolated component and test if
additional PTRs can unblock them (even though they're now topologically connected).
"""

import csv
import json
import os
import time
import hashlib
import pickle
import networkx as nx
from collections import defaultdict
from cobra.io import load_json_model, save_json_model
from cobra import Reaction, Configuration


# Available LP solvers through optlang (used by COBRApy)
# Note: gurobi and cplex require commercial licenses
VALID_LP_SOLVERS = ["glpk", "glpk_exact", "scipy", "gurobi", "cplex"]


def set_lp_solver(solver_name):
    """Set the LP solver used by COBRApy for FBA.

    Args:
        solver_name: One of 'glpk', 'glpk_exact', 'scipy', 'gurobi', 'cplex'
                     - glpk: Default, good balance of speed and reliability
                     - glpk_exact: Uses exact arithmetic, slower but more precise
                     - scipy: Pure Python fallback, no external dependencies
                     - gurobi: Commercial solver, very fast (requires license)
                     - cplex: Commercial solver, very fast (requires license)

    Returns:
        True if solver was set successfully, False otherwise
    """
    if solver_name not in VALID_LP_SOLVERS:
        print(
            f"Warning: Unknown solver '{solver_name}'. Valid options: {VALID_LP_SOLVERS}"
        )
        return False

    try:
        Configuration().solver = solver_name
        print(f"LP solver set to: {solver_name}")
        return True
    except Exception as e:
        print(f"Warning: Could not set solver '{solver_name}': {e}")
        print(f"Make sure the solver is installed and properly configured.")
        if solver_name == "gurobi":
            print("  For Gurobi: pip install gurobipy and ensure license is active")
        elif solver_name == "cplex":
            print("  For CPLEX: pip install cplex and ensure license is active")
        return False


def compute_model_hash(model_json_path):
    """Compute hash of model file for cache validation."""
    with open(model_json_path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()[:12]


def get_coverage_cache_path(cache_dir, model_hash, component_id):
    """Get path to coverage cache file for a specific component."""
    os.makedirs(cache_dir, exist_ok=True)
    return os.path.join(cache_dir, f"coverage_{model_hash}_comp{component_id}.pkl")


def load_coverage_from_cache(cache_path, candidate_keys, verbose=True):
    """Load coverage matrix from cache if valid.

    Returns: coverage dict or None if cache invalid/missing
    """
    if not os.path.exists(cache_path):
        return None

    try:
        with open(cache_path, "rb") as f:
            cached = pickle.load(f)

        # Validate cache
        if cached.get("candidate_keys") != candidate_keys:
            if verbose:
                print("    Cache invalid: candidate set changed")
            return None

        if verbose:
            print(f"    Loaded coverage from cache ({os.path.basename(cache_path)})")
            age_mins = (time.time() - cached.get("timestamp", 0)) / 60
            print(f"    Cache age: {age_mins:.1f} minutes")

        return cached["coverage"]
    except Exception as e:
        if verbose:
            print(f"    Cache load failed: {e}")
        return None


def save_coverage_to_cache(cache_path, coverage, candidate_keys, verbose=True):
    """Save coverage matrix to cache."""
    try:
        cache_data = {
            "coverage": coverage,
            "candidate_keys": candidate_keys,
            "timestamp": time.time(),
            "num_candidates": len(candidate_keys),
            "num_entries": sum(len(v) for v in coverage.values()),
        }
        with open(cache_path, "wb") as f:
            pickle.dump(cache_data, f)

        if verbose:
            size_mb = os.path.getsize(cache_path) / (1024 * 1024)
            print(f"    Saved coverage to cache ({size_mb:.2f} MB)")
    except Exception as e:
        if verbose:
            print(f"    Cache save failed: {e}")


def add_transport(model, sel, idx, prefix="SINK"):
    """Add a PTR (pseudo-transport-reaction) to the model.

    Args:
        model: COBRApy model
        sel: PTR dict with 'met1', 'met2', 'base', 'suffix1', 'suffix2', etc.
        idx: Index for unique reaction ID
        prefix: Prefix for reaction ID (default: "SINK")

    Returns:
        reaction ID if successful, None if metabolites not found
    """
    met1 = sel["met1"]
    met2 = sel["met2"]
    base = sel.get("base", "")
    s1 = sel.get("suffix1", "")
    s2 = sel.get("suffix2", "")

    # Generate unique reaction ID
    raw_id = f"{prefix}_TRANS_{base}_{s1}_{s2}_{idx}"
    rid = "".join(ch if (ch.isalnum() or ch == "_") else "_" for ch in raw_id)

    # Get metabolites
    try:
        m1 = model.metabolites.get_by_id(met1)
        m2 = model.metabolites.get_by_id(met2)
    except KeyError:
        return None

    # Create transport reaction
    rxn = Reaction(rid)
    rxn.name = f"Phase3 sink MILP transport {met1} <-> {met2}"
    rxn.add_metabolites({m1: -1.0, m2: 1.0})
    rxn.lower_bound = -1000.0
    rxn.upper_bound = 1000.0
    rxn.annotation = {
        "phase3_selected": "true",
        "strategy": "sink_milp",
        "met1": met1,
        "met2": met2,
    }

    model.add_reactions([rxn])
    return rid


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
    """Load component summary CSV from Phase 1."""
    components = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            cid = int(r["component_id"])
            size = int(r["num_nodes"])
            components[cid] = size
    return components


def build_original_component_map(unconnected_model):
    """
    Build component mapping from the ORIGINAL unconnected model.
    This preserves which reactions/metabolites belonged to which component
    before Phase 2 connected them.
    """
    G = nx.Graph()
    for rxn in unconnected_model.reactions:
        rid = rxn.id
        G.add_node(rid, bipartite=0)
        for met in rxn.metabolites:
            mid = met.id
            G.add_node(mid, bipartite=1)
            G.add_edge(rid, mid)

    components = list(nx.connected_components(G))
    components.sort(key=len, reverse=True)  # largest = main = ID 1

    met_to_comp = {}
    rxn_to_comp = {}
    comp_sizes = {}

    for comp_id, comp_nodes in enumerate(components, start=1):
        comp_sizes[comp_id] = len(comp_nodes)
        for node in comp_nodes:
            if node.startswith("MAM"):
                met_to_comp[node] = comp_id
            else:
                rxn_to_comp[node] = comp_id

    return met_to_comp, rxn_to_comp, comp_sizes


def get_deadend_info(model):
    """
    Identify dead-end metabolites and their types.
    """
    deadend_info = {}

    for met in model.metabolites:
        producing = []
        consuming = []

        for rxn in met.reactions:
            coeff = rxn.metabolites[met]
            lb, ub = rxn.lower_bound, rxn.upper_bound

            can_produce = (coeff > 0 and ub > 0) or (coeff < 0 and lb < 0)
            can_consume = (coeff < 0 and ub > 0) or (coeff > 0 and lb < 0)

            if can_produce:
                producing.append(rxn.id)
            if can_consume:
                consuming.append(rxn.id)

        if not producing and consuming:
            de_type = "substrate"
        elif producing and not consuming:
            de_type = "product"
        elif not producing and not consuming:
            de_type = "isolated"
        else:
            de_type = "none"

        if de_type != "none":
            deadend_info[met.id] = {
                "type": de_type,
                "producing_rxns": producing,
                "consuming_rxns": consuming,
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
                    if sol.status == "optimal":
                        can_carry_flux = True
                except Exception:
                    pass

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


def add_temp_sink_source(model, met_id, deadend_info):
    """Add temporary sink/source for a dead-end metabolite.

    Only adds exchanges for true dead-ends:
    - substrate: consumed but never produced → add SOURCE
    - product: produced but never consumed → add SINK

    Does NOT add exchanges for 'isolated' metabolites - those have both
    producers and consumers but are blocked due to topological isolation,
    which is what the PTR candidates are meant to solve.
    """
    if met_id not in deadend_info:
        return None

    info = deadend_info[met_id]
    de_type = info["type"]

    # Only handle true dead-ends (substrate/product), not isolated
    if de_type not in ("substrate", "product"):
        return None

    try:
        met = model.metabolites.get_by_id(met_id)
    except KeyError:
        return None

    rxn_id = f"_TEMP_EXCH_{met_id}_"

    # Check if already added
    if rxn_id in model.reactions:
        return rxn_id

    rxn = Reaction(rxn_id)
    rxn.name = f"Temporary exchange for {met_id}"

    if de_type == "substrate":
        # Dead-end substrate: needs source (→ met)
        rxn.add_metabolites({met: 1.0})
        rxn.lower_bound = 0
        rxn.upper_bound = 1000.0
    elif de_type == "product":
        # Dead-end product: needs sink (met →)
        rxn.add_metabolites({met: -1.0})
        rxn.lower_bound = 0
        rxn.upper_bound = 1000.0

    model.add_reactions([rxn])
    return rxn_id


def test_candidate_with_temp_sinks(
    model, candidate, blocked_rxn_id, deadend_info, zero_cutoff=None
):
    """
    Test if PTR + temporary sinks can enable flux through a blocked reaction.

    NOTE: This function assumes temp sinks are already added to the model.
    """
    if zero_cutoff is None:
        zero_cutoff = model.tolerance * 10

    met1 = candidate["met1"]
    met2 = candidate["met2"]

    try:
        m1 = model.metabolites.get_by_id(met1)
        m2 = model.metabolites.get_by_id(met2)
        blocked_rxn = model.reactions.get_by_id(blocked_rxn_id)
    except KeyError:
        return False

    with model:
        # Add candidate PTR
        ptr = Reaction("_TEST_PTR_")
        ptr.add_metabolites({m1: -1.0, m2: 1.0})
        ptr.lower_bound = -1000.0
        ptr.upper_bound = 1000.0
        model.add_reactions([ptr])

        # Test if blocked reaction can carry flux
        can_carry = False

        if blocked_rxn.upper_bound > zero_cutoff:
            with model:
                blocked_rxn.lower_bound = zero_cutoff
                try:
                    sol = model.optimize()
                    if sol.status == "optimal":
                        can_carry = True
                except Exception:
                    pass

        if not can_carry and blocked_rxn.lower_bound < -zero_cutoff:
            with model:
                blocked_rxn.upper_bound = -zero_cutoff
                try:
                    sol = model.optimize()
                    if sol.status == "optimal":
                        can_carry = True
                except Exception:
                    pass

        return can_carry


def test_candidate_batch(model, candidate, blocked_rxn_ids, zero_cutoff=None):
    """
    Test which blocked reactions can carry flux when a PTR is added.

    Uses simple FBA (not FVA): we only need to know IF a reaction can carry flux,
    not its full range. This is 2x faster than FVA.

    For each reaction: try max(rxn), if > 0 it's unblocked.
    If not, try min(rxn), if < 0 it's unblocked.
    Worst case: 2 FBAs per reaction, often just 1.

    NOTE: Assumes temp sinks are already added to the model.

    Returns: set of reaction IDs that become unblocked
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

    unblocked = set()

    with model:
        # Add candidate PTR once
        ptr = Reaction("_TEST_PTR_")
        ptr.add_metabolites({m1: -1.0, m2: 1.0})
        ptr.lower_bound = -1000.0
        ptr.upper_bound = 1000.0
        model.add_reactions([ptr])

        # Test each blocked reaction with simple FBA
        for rid in blocked_rxn_ids:
            try:
                rxn = model.reactions.get_by_id(rid)
            except KeyError:
                continue

            can_carry = False

            # Try forward direction first (max)
            if rxn.upper_bound > zero_cutoff:
                model.objective = rxn.id
                model.objective_direction = "max"
                try:
                    val = model.slim_optimize()
                    if val is not None and val > zero_cutoff:
                        can_carry = True
                except Exception:
                    pass

            # If forward didn't work, try reverse direction (min)
            if not can_carry and rxn.lower_bound < -zero_cutoff:
                model.objective = rxn.id
                model.objective_direction = "min"
                try:
                    val = model.slim_optimize()
                    if val is not None and val < -zero_cutoff:
                        can_carry = True
                except Exception:
                    pass

            if can_carry:
                unblocked.add(rid)

    return unblocked


# Global variable for parallel worker (set by initializer)
_worker_model = None


def _init_worker(model_json_path, temp_sink_mets, deadend_types, solver_lp="glpk"):
    """Initialize worker process with its own model copy.

    Args:
        model_json_path: Path to model JSON file
        temp_sink_mets: List of metabolite IDs that need temp sinks
        deadend_types: List of dead-end types ('substrate' or 'product')
        solver_lp: LP solver to use for FBA
    """
    global _worker_model

    # Set LP solver for this worker process
    from cobra import Configuration

    try:
        Configuration().solver = solver_lp
    except Exception:
        pass  # Fall back to default

    _worker_model = load_json_model(model_json_path)

    # Add temp sinks to worker's model
    for met_id, de_type in zip(temp_sink_mets, deadend_types):
        if met_id not in _worker_model.metabolites:
            continue
        met = _worker_model.metabolites.get_by_id(met_id)
        rxn_id = f"_TEMP_EXCH_{met_id}_"
        if rxn_id in _worker_model.reactions:
            continue
        rxn = Reaction(rxn_id)
        if de_type == "substrate":
            rxn.add_metabolites({met: 1.0})
            rxn.bounds = (0, 1000.0)
        elif de_type == "product":
            rxn.add_metabolites({met: -1.0})
            rxn.bounds = (0, 1000.0)
        _worker_model.add_reactions([rxn])


def _test_single_reaction(args):
    """Worker function to test if a single reaction can carry flux with a PTR."""
    rid, met1, met2, zero_cutoff = args
    global _worker_model

    if _worker_model is None:
        return (rid, False)

    try:
        m1 = _worker_model.metabolites.get_by_id(met1)
        m2 = _worker_model.metabolites.get_by_id(met2)
        rxn = _worker_model.reactions.get_by_id(rid)
    except KeyError:
        return (rid, False)

    with _worker_model:
        # Add candidate PTR
        ptr = Reaction("_TEST_PTR_")
        ptr.add_metabolites({m1: -1.0, m2: 1.0})
        ptr.lower_bound = -1000.0
        ptr.upper_bound = 1000.0
        _worker_model.add_reactions([ptr])

        can_carry = False

        # Try forward direction
        if rxn.upper_bound > zero_cutoff:
            _worker_model.objective = rxn.id
            _worker_model.objective_direction = "max"
            try:
                val = _worker_model.slim_optimize()
                if val is not None and val > zero_cutoff:
                    can_carry = True
            except Exception:
                pass

        # Try reverse direction
        if not can_carry and rxn.lower_bound < -zero_cutoff:
            _worker_model.objective = rxn.id
            _worker_model.objective_direction = "min"
            try:
                val = _worker_model.slim_optimize()
                if val is not None and val < -zero_cutoff:
                    can_carry = True
            except Exception:
                pass

    return (rid, can_carry)


def test_candidate_batch_parallel(
    model_json_path,
    candidate,
    blocked_rxn_ids,
    temp_sink_mets,
    deadend_types,
    n_workers=None,
    zero_cutoff=1e-8,
    pool=None,
):
    """
    Parallel version of test_candidate_batch.

    Tests all blocked reactions in parallel using multiprocessing.
    Each worker has its own model copy to avoid conflicts.

    Args:
        model_json_path: Path to the model JSON file (workers load their own copy)
        candidate: PTR candidate dict with 'met1' and 'met2'
        blocked_rxn_ids: List of reaction IDs to test
        temp_sink_mets: List of metabolite IDs that have temp sinks
        deadend_types: List of types ('substrate' or 'product') for each temp sink met
        n_workers: Number of parallel workers (default: 2, max 4 to limit memory)
        zero_cutoff: Threshold for considering flux as non-zero
        pool: Optional pre-created Pool to reuse (avoids reinitializing workers)

    Returns: set of reaction IDs that become unblocked
    """
    from multiprocessing import Pool, cpu_count

    # MEMORY SAFETY: Limit workers to avoid RAM exhaustion
    # Each worker loads ~500MB model, so limit to 4 max
    if n_workers is None:
        n_workers = min(2, cpu_count())  # Default to 2 workers
    n_workers = min(n_workers, 4)  # Hard cap at 4

    if n_workers <= 1 or len(blocked_rxn_ids) < 4:
        # Fall back to sequential for small problems
        model = load_json_model(model_json_path)
        # Add temp sinks
        for met_id, de_type in zip(temp_sink_mets, deadend_types):
            if met_id not in model.metabolites:
                continue
            met = model.metabolites.get_by_id(met_id)
            rxn_id = f"_TEMP_EXCH_{met_id}_"
            if rxn_id in model.reactions:
                continue
            rxn = Reaction(rxn_id)
            if de_type == "substrate":
                rxn.add_metabolites({met: 1.0})
                rxn.bounds = (0, 1000.0)
            elif de_type == "product":
                rxn.add_metabolites({met: -1.0})
                rxn.bounds = (0, 1000.0)
            model.add_reactions([rxn])
        return test_candidate_batch(model, candidate, blocked_rxn_ids, zero_cutoff)

    met1 = candidate["met1"]
    met2 = candidate["met2"]

    # Prepare arguments for workers
    args_list = [(rid, met1, met2, zero_cutoff) for rid in blocked_rxn_ids]

    # Run in parallel - use provided pool or create a new one
    unblocked = set()

    if pool is not None:
        # Reuse existing pool
        results = pool.map(_test_single_reaction, args_list)
        for rid, can_carry in results:
            if can_carry:
                unblocked.add(rid)
    else:
        # Create new pool (more expensive, but self-contained)
        with Pool(
            processes=n_workers,
            initializer=_init_worker,
            initargs=(model_json_path, temp_sink_mets, deadend_types),
        ) as pool:
            results = pool.map(_test_single_reaction, args_list)
            for rid, can_carry in results:
                if can_carry:
                    unblocked.add(rid)

    return unblocked


def prepare_model_with_temp_sinks(model, deadend_info):
    """
    Add ALL temporary sinks to the model ONCE.
    Returns the modified model (in-place modification).
    """
    added = []
    for mid in deadend_info:
        rxn_id = add_temp_sink_source(model, mid, deadend_info)
        if rxn_id:
            added.append(rxn_id)
    return added


def build_component_graph(model, component_rxn_ids):
    """
    Build a graph of the component for reachability analysis.

    Returns a NetworkX graph where:
    - Nodes are metabolites and reactions within the component
    - Edges connect reactions to their metabolites
    """
    G = nx.Graph()

    for rid in component_rxn_ids:
        try:
            rxn = model.reactions.get_by_id(rid)
        except KeyError:
            continue

        G.add_node(rid, node_type="reaction")
        for met in rxn.metabolites:
            G.add_node(met.id, node_type="metabolite")
            G.add_edge(rid, met.id)

    return G


def get_reachable_reactions(graph, metabolite_id, component_rxn_ids):
    """
    Find which reactions in the component can be reached from a metabolite.

    Uses BFS from the metabolite to find all connected reactions.
    """
    if metabolite_id not in graph:
        return set()

    # BFS to find all reachable nodes
    reachable = set(nx.node_connected_component(graph, metabolite_id))

    # Filter to only reactions in our component
    return reachable & component_rxn_ids


def precompute_candidate_reachability(graph, candidates, blocked_rxns):
    """
    Pre-compute which blocked reactions each candidate can potentially affect.

    A candidate PTR (m1 <-> m2) can only affect a blocked reaction R if:
    - m1 or m2 is graph-reachable to R within the component

    This allows us to skip FBA tests for (candidate, reaction) pairs
    that have no graph path, reducing FBAs by 50-90%.

    Returns: dict mapping (met1, met2) -> set of potentially affected reaction IDs
    """
    blocked_set = set(blocked_rxns)
    reachability = {}

    for cand in candidates:
        met1, met2 = cand["met1"], cand["met2"]
        key = (met1, met2)

        # Find reactions reachable from either endpoint
        reachable_from_m1 = get_reachable_reactions(graph, met1, blocked_set)
        reachable_from_m2 = get_reachable_reactions(graph, met2, blocked_set)

        reachability[key] = reachable_from_m1 | reachable_from_m2

    return reachability


def compute_coverage_with_temp_sinks(
    model,
    candidates,
    blocked_rxns,
    deadend_info,
    component_rxn_ids=None,
    verbose=True,
    early_termination=True,
    stagnation_limit=20,
    parallel=False,
    n_workers=None,
    model_json_path=None,
    solver_lp="glpk",
):
    """
    Compute coverage matrix using temporary sink/source approach.

    OPTIMIZATIONS:
    1. Add all temp sinks ONCE (not per candidate)
    2. Graph pre-filtering: skip (candidate, reaction) pairs with no graph path
    3. Early termination: stop when all blocked reactions are covered
    4. Stagnation detection: stop if no new coverage for N candidates in a row
    5. Parallel FBA testing (optional): test multiple reactions in parallel

    Args:
        model: cobra model with temp sinks NOT yet added
        candidates: list of candidate PTRs (should be sorted by priority: A, B, C)
        blocked_rxns: set of blocked reaction IDs
        deadend_info: dict from get_deadend_info()
        component_rxn_ids: set of reaction IDs in the original component (for graph)
        verbose: print progress
        early_termination: stop when 100% coverage reached
        stagnation_limit: stop if no new coverage for this many consecutive candidates
        parallel: use parallel FBA testing
        n_workers: number of parallel workers (default: CPU count)
        model_json_path: path to model JSON (required for parallel mode)
        solver_lp: LP solver to use for FBA (glpk, glpk_exact, scipy)
    """
    coverage = {}
    total = len(candidates)
    n_blocked = len(blocked_rxns)

    if verbose:
        print(
            f"      Coverage (with temp sinks): {total} candidates × {n_blocked} blocked"
        )

    # Prepare for parallel mode if enabled
    pool = None  # Will be created once and reused
    temp_sink_mets = []
    deadend_types = []

    if parallel:
        if model_json_path is None:
            if verbose:
                print(
                    "      Warning: parallel mode requires model_json_path, falling back to sequential"
                )
            parallel = False
        else:
            # Prepare temp sink info for workers
            for mid, info in deadend_info.items():
                if info["type"] in ("substrate", "product"):
                    temp_sink_mets.append(mid)
                    deadend_types.append(info["type"])

            # MEMORY SAFETY: Limit workers to avoid RAM exhaustion
            from multiprocessing import Pool, cpu_count

            if n_workers is None:
                n_workers = 2  # Conservative default
            
            # allow user to brick their computer if they want
            # n_workers = min(n_workers, 4)  # Hard cap at 4 workers (~2GB RAM)

            if verbose:
                print(f"      Parallel mode: {n_workers} workers (memory-safe limit)")

            # Create ONE pool and reuse for all candidates
            # Pass solver_lp to workers so they use the same LP solver
            pool = Pool(
                processes=n_workers,
                initializer=_init_worker,
                initargs=(model_json_path, temp_sink_mets, deadend_types, solver_lp),
            )

    # Add ALL temp sinks to model ONCE (outside the loop) - for sequential mode
    if not parallel:
        if verbose:
            print(f"      Adding {len(deadend_info)} temp sinks to model...")
        temp_rxns = prepare_model_with_temp_sinks(model, deadend_info)
        if verbose:
            print(f"      Added {len(temp_rxns)} temp exchange reactions")

    # OPTIMIZATION 1: Build component graph for reachability pre-filtering
    if component_rxn_ids:
        if verbose:
            print(f"      Building component graph for pre-filtering...")
        comp_graph = build_component_graph(model, component_rxn_ids)
        reachability = precompute_candidate_reachability(
            comp_graph, candidates, blocked_rxns
        )

        # Stats on pre-filtering
        total_pairs = total * n_blocked
        reachable_pairs = sum(len(r) for r in reachability.values())
        if verbose:
            print(
                f"      Pre-filter: {reachable_pairs}/{total_pairs} pairs reachable ({100*reachable_pairs/total_pairs:.1f}%)"
            )
    else:
        reachability = None

    start_time = time.time()
    blocked_list = list(blocked_rxns)

    # OPTIMIZATION 3: Track cumulative coverage for early termination
    cumulative_coverage = set()
    skipped_by_prefilter = 0
    skipped_by_early_term = 0
    stagnation_counter = 0  # Track consecutive candidates with no new coverage

    for i, cand in enumerate(candidates):
        key = (cand["met1"], cand["met2"])

        # OPTIMIZATION 3a: Early termination - 100% coverage
        if early_termination and cumulative_coverage == blocked_rxns:
            skipped_by_early_term = total - i
            if verbose:
                print(
                    f"      Early termination: 100% coverage reached after {i} candidates"
                )
                print(f"      Skipping remaining {skipped_by_early_term} candidates")
            # Fill remaining coverage entries as empty (not tested)
            for j in range(i, total):
                remaining_key = (candidates[j]["met1"], candidates[j]["met2"])
                coverage[remaining_key] = set()
            break

        # OPTIMIZATION 3b: Stagnation detection - no new coverage for N candidates
        if stagnation_limit and stagnation_counter >= stagnation_limit:
            skipped_by_early_term = total - i
            if verbose:
                pct = 100 * len(cumulative_coverage) / n_blocked if n_blocked > 0 else 0
                print(
                    f"      Stagnation: no new coverage for {stagnation_limit} candidates"
                )
                print(
                    f"      Current coverage: {len(cumulative_coverage)}/{n_blocked} ({pct:.1f}%)"
                )
                print(f"      Skipping remaining {skipped_by_early_term} candidates")
            for j in range(i, total):
                remaining_key = (candidates[j]["met1"], candidates[j]["met2"])
                coverage[remaining_key] = set()
            break

        if verbose and (i + 1) % 10 == 0:
            elapsed = time.time() - start_time
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            eta = (total - i - 1) / rate if rate > 0 else 0
            pct_covered = (
                100 * len(cumulative_coverage) / n_blocked if n_blocked > 0 else 0
            )
            print(
                f"        {i+1}/{total} ({rate:.2f}/sec, ETA: {eta:.0f}s, coverage: {pct_covered:.1f}%)"
            )

        # OPTIMIZATION 1: Only test reachable reactions
        if reachability:
            rxns_to_test = list(reachability.get(key, set()))
            skipped_by_prefilter += n_blocked - len(rxns_to_test)
        else:
            rxns_to_test = blocked_list

        if not rxns_to_test:
            coverage[key] = set()
            continue

        # Use batch testing - adds PTR once, tests filtered reactions
        if parallel and pool is not None:
            # Use reusable pool for parallel FBA testing
            met1 = cand["met1"]
            met2 = cand["met2"]
            args_list = [
                (rid, met1, met2, model.tolerance * 10) for rid in rxns_to_test
            ]
            results = pool.map(_test_single_reaction, args_list)
            unblocked = set(rid for rid, can_carry in results if can_carry)
        else:
            unblocked = test_candidate_batch(model, cand, rxns_to_test)

        coverage[key] = unblocked

        # Track coverage growth for stagnation detection
        prev_coverage_size = len(cumulative_coverage)
        cumulative_coverage |= unblocked
        if len(cumulative_coverage) > prev_coverage_size:
            stagnation_counter = 0  # Reset on progress
        else:
            stagnation_counter += 1

    if verbose:
        elapsed = time.time() - start_time
        non_empty = sum(1 for v in coverage.values() if v)
        total_pairs = sum(len(v) for v in coverage.values())
        print(
            f"      Done in {elapsed:.1f}s: {non_empty}/{total} candidates useful, {total_pairs} total unblocks"
        )
        if skipped_by_prefilter > 0:
            print(f"      Pre-filter saved: {skipped_by_prefilter} FBA pairs skipped")
        if skipped_by_early_term > 0:
            print(
                f"      Early termination saved: {skipped_by_early_term} candidates not tested"
            )
        pct_covered = 100 * len(cumulative_coverage) / n_blocked if n_blocked > 0 else 0
        print(
            f"      Final coverage: {len(cumulative_coverage)}/{n_blocked} blocked reactions ({pct_covered:.1f}%)"
        )

    # Clean up parallel pool
    if pool is not None:
        pool.close()
        pool.join()

    return coverage


def solve_greedy(coverage, blocked_rxns, candidates, verbose=True):
    """
    Greedy set-cover solver: iteratively pick the candidate that covers the most
    uncovered blocked reactions. Much faster than MILP for large problems.
    """
    uncovered = set(blocked_rxns)
    selected = []

    while uncovered:
        # Find candidate that covers most uncovered reactions
        best_cand = None
        best_coverage = set()

        for c in candidates:
            key = (c["met1"], c["met2"])
            cov = coverage.get(key, set())
            new_cov = cov & uncovered
            if len(new_cov) > len(best_coverage):
                best_coverage = new_cov
                best_cand = c

        if not best_cand or not best_coverage:
            break  # No more progress possible

        selected.append(best_cand)
        uncovered -= best_coverage

        if verbose:
            pct = 100 * (len(blocked_rxns) - len(uncovered)) / len(blocked_rxns)
            print(
                f"      Greedy: selected {best_cand['met1']}<->{best_cand['met2']} "
                f"(+{len(best_coverage)} rxns, total coverage: {pct:.1f}%)"
            )

    if verbose:
        total_covered = len(blocked_rxns) - len(uncovered)
        print(
            f"      Greedy: {len(selected)} PTRs → {total_covered}/{len(blocked_rxns)} reactions covered"
        )

    return selected


def solve_milp(coverage, blocked_rxns, candidates, tradeoff_lambda=0.01, verbose=True):
    """Solve MILP for coverage optimization."""
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

    raise ImportError("No MILP solver available")


def _solve_pulp(coverage, blocked_rxns, candidates, tradeoff_lambda, verbose):
    import pulp

    blocked_list = list(blocked_rxns)
    blocked_idx = {b: i for i, b in enumerate(blocked_list)}
    cand_idx = {(c["met1"], c["met2"]): i for i, c in enumerate(candidates)}

    n_cands = len(candidates)
    n_blocked = len(blocked_list)

    prob = pulp.LpProblem("coverage", pulp.LpMaximize)

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

    selected = [
        c for i, c in enumerate(candidates) if y[i].value() and y[i].value() > 0.5
    ]

    if verbose:
        n_unblocked = sum(
            1 for b in range(n_blocked) if z[b].value() and z[b].value() > 0.5
        )
        print(
            f"      MILP: {len(selected)} PTRs → {n_unblocked} reactions can be unblocked"
        )

    return selected


def _solve_scipy(coverage, blocked_rxns, candidates, tradeoff_lambda, verbose):
    import numpy as np
    from scipy.optimize import milp, LinearConstraint, Bounds

    blocked_list = list(blocked_rxns)
    blocked_idx = {b: i for i, b in enumerate(blocked_list)}
    cand_idx = {(c["met1"], c["met2"]): i for i, c in enumerate(candidates)}

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
        print(
            f"      MILP: {len(selected)} PTRs → {n_unblocked} reactions can be unblocked"
        )

    return selected


def get_candidates_for_original_component(
    candidates, component_id, used_candidates=None
):
    """
    Get candidates that connect TO the given ORIGINAL component.
    Uses component1/component2 fields from Phase-1 candidates CSV.
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


def run_test_on_original_component(
    candidates_csv,
    phase2_model_json,
    unconnected_model_json,
    target_component_id=2,  # Original component ID from Phase 1
    out_dir=None,
    phase3_model_out=None,
    tradeoff_lambda=0.01,
    max_candidates=None,  # Limit candidates for testing
    solver="milp",  # 'milp', 'greedy', or 'dynamic'
    solver_lp="glpk",  # LP solver: 'glpk', 'glpk_exact', 'scipy'
    sample_blocked=None,  # Sample N blocked reactions (for large components)
    parallel=False,  # Use parallel FBA testing
    n_workers=None,  # Number of parallel workers
    use_cache=True,  # Use coverage cache
    cache_dir=None,  # Cache directory
    verbose=True,
):
    """
    Run test on an ORIGINAL component (using Phase-1 component assignments).

    Args:
        target_component_id: Original component ID from Phase 1 (1=main, 2=largest isolated, etc.)
        max_candidates: Limit number of candidates for quick testing
        solver: 'milp' for optimal solution, 'greedy' for fast approximation, 'dynamic' for auto
        solver_lp: LP solver for FBA - 'glpk' (default), 'glpk_exact' (more precise), 'scipy' (fallback)
        sample_blocked: Sample N blocked reactions (speeds up large components)
        parallel: Use parallel FBA testing (speeds up intra-component processing)
        n_workers: Number of parallel workers (default: CPU count)
    """
    # Set LP solver for FBA before loading any models
    if solver_lp:
        set_lp_solver(solver_lp)

    base_dir = os.path.dirname(__file__)
    if out_dir is None:
        out_dir = os.path.join(base_dir, "files")
    os.makedirs(out_dir, exist_ok=True)

    # Load data
    candidates = load_candidates(candidates_csv)
    phase2_model = load_json_model(phase2_model_json)
    unconnected_model = load_json_model(unconnected_model_json)

    if verbose:
        print(f"Loaded {len(candidates)} candidates")
        print(f"Phase-2 model: {len(phase2_model.reactions)} reactions")
        print(f"Unconnected model: {len(unconnected_model.reactions)} reactions")

    # Build ORIGINAL component mapping from unconnected model
    if verbose:
        print("\nBuilding ORIGINAL component mapping (from unconnected model)...")
    met_to_comp, rxn_to_comp, comp_sizes = build_original_component_map(
        unconnected_model
    )

    if verbose:
        print(f"Found {len(comp_sizes)} original components:")
        sorted_comps = sorted(comp_sizes.items(), key=lambda x: x[1], reverse=True)
        for cid, size in sorted_comps[:10]:
            rxns_in_comp = sum(1 for r, c in rxn_to_comp.items() if c == cid)
            print(f"  Original Component {cid}: {size} nodes ({rxns_in_comp} rxns)")

    # Get reactions originally in target component
    rxns_in_orig_comp = [
        rid for rid, cid in rxn_to_comp.items() if cid == target_component_id
    ]
    mets_in_orig_comp = [
        mid for mid, cid in met_to_comp.items() if cid == target_component_id
    ]

    if verbose:
        print(f"\n{'='*60}")
        print(f"Testing ORIGINAL Component {target_component_id}")
        print(f"  Original size: {comp_sizes.get(target_component_id, 0)} nodes")
        print(f"  Reactions: {len(rxns_in_orig_comp)}")
        print(f"  Metabolites: {len(mets_in_orig_comp)}")
        print(f"{'='*60}")

    if not rxns_in_orig_comp:
        print("No reactions in this component!")
        return None

    # Get dead-end info from Phase-2 model (current state)
    if verbose:
        print("\nAnalyzing dead-end metabolites in Phase-2 model...")
    deadend_info = get_deadend_info(phase2_model)

    deadends_in_orig = {
        m: deadend_info[m] for m in mets_in_orig_comp if m in deadend_info
    }
    if verbose:
        by_type = defaultdict(int)
        for info in deadends_in_orig.values():
            by_type[info["type"]] += 1
        print(f"  Dead-ends (original component mets): {len(deadends_in_orig)}")
        for t, c in sorted(by_type.items()):
            print(f"    {t}: {c}")

    # Find blocked reactions (in Phase-2 model, but originally from target component)
    # Only check reactions that exist in Phase-2 model
    rxns_to_check = [r for r in rxns_in_orig_comp if r in phase2_model.reactions]

    if verbose:
        print(
            f"\nFinding blocked reactions among {len(rxns_to_check)} original component reactions..."
        )
    blocked = find_blocked_reactions_in_set(phase2_model, rxns_to_check)

    if verbose:
        print(f"  Blocked: {len(blocked)} / {len(rxns_to_check)}")
        if blocked and len(blocked) <= 20:
            print(f"  Blocked IDs: {list(blocked)}")

    if not blocked:
        print("No blocked reactions in this original component!")
        return {"blocked": 0, "effective_candidates": 0, "selected": 0}

    # Sample blocked reactions if requested (speeds up large components)
    original_blocked = blocked
    if sample_blocked and len(blocked) > sample_blocked:
        import random

        blocked_list = list(blocked)
        random.shuffle(blocked_list)
        blocked = set(blocked_list[:sample_blocked])
        if verbose:
            print(
                f"  Sampled {len(blocked)} blocked reactions for faster coverage computation"
            )

    # Get candidates for this ORIGINAL component
    comp_candidates = get_candidates_for_original_component(
        candidates, target_component_id
    )

    if verbose:
        print(
            f"\nCandidates for original component {target_component_id}: {len(comp_candidates)}"
        )
        by_type = defaultdict(int)
        for c in comp_candidates:
            by_type[c.get("type", "?")] += 1
        for t, cnt in sorted(by_type.items()):
            print(f"  Type {t}: {cnt}")

    if not comp_candidates:
        print("No candidates available!")
        return {"blocked": len(blocked), "effective_candidates": 0, "selected": 0}

    # Limit candidates for testing if requested
    if max_candidates and len(comp_candidates) > max_candidates:
        # Prioritize Type A, then B, then C
        type_a = [c for c in comp_candidates if c.get("type") == "A"]
        type_b = [c for c in comp_candidates if c.get("type") == "B"]
        type_c = [c for c in comp_candidates if c.get("type") == "C"]

        comp_candidates = (type_a + type_b + type_c)[:max_candidates]
        if verbose:
            print(f"  Limited to {len(comp_candidates)} candidates for testing")

    # Sort candidates by priority (A > B > C) for early termination optimization
    comp_candidates = sorted(
        comp_candidates,
        key=lambda c: {"A": 0, "B": 1, "C": 2}.get(c.get("type", "C"), 3),
    )

    # Calculate adaptive stagnation limit based on problem characteristics
    num_type_a = sum(1 for c in comp_candidates if c.get("type") == "A")
    patience_limit = max(20, len(blocked) // 3)  # Scales with blocked reactions
    stagnation_limit = max(num_type_a, patience_limit)  # Always test all Type A

    if verbose:
        print(
            f"  Stagnation limit: {stagnation_limit} (Type A: {num_type_a}, patience: {patience_limit})"
        )

    # Compute coverage with temp sinks (with caching)
    if verbose:
        print("\nComputing coverage matrix (with temp sinks)...")

    coverage = None
    candidate_keys = {(c["met1"], c["met2"]) for c in comp_candidates}

    # Try to load from cache
    if use_cache:
        if cache_dir is None:
            cache_dir = os.path.join(out_dir, "cache")

        model_hash = compute_model_hash(phase2_model_json)
        cache_path = get_coverage_cache_path(cache_dir, model_hash, target_component_id)
        coverage = load_coverage_from_cache(cache_path, candidate_keys, verbose)

    # Compute if not cached
    if coverage is None:
        # Pass component_rxn_ids for graph pre-filtering optimization
        component_rxn_ids = set(rxns_to_check)
        coverage = compute_coverage_with_temp_sinks(
            phase2_model,
            comp_candidates,
            blocked,
            deadend_info,
            component_rxn_ids=component_rxn_ids,
            verbose=verbose,
            early_termination=True,
            stagnation_limit=stagnation_limit,
            parallel=parallel,
            n_workers=n_workers,
            model_json_path=phase2_model_json,
            solver_lp=solver_lp,
        )

        # Save to cache
        if use_cache:
            save_coverage_to_cache(cache_path, coverage, candidate_keys, verbose)

    # Filter to effective candidates
    effective = [
        c for c in comp_candidates if coverage.get((c["met1"], c["met2"]), set())
    ]

    if verbose:
        print(f"\nEffective candidates (can unblock ≥1 reaction): {len(effective)}")

    if not effective:
        print("No candidates can unblock any reactions even with temp sinks!")
        return {"blocked": len(blocked), "effective_candidates": 0, "selected": 0}

    # Solve for PTR selection
    actual_solver = solver

    if solver == "dynamic":
        # Automatically choose solver based on problem characteristics
        n_eff = len(effective)
        n_blocked = len(blocked)
        problem_size = n_eff * n_blocked

        # Check if coverage is simple (few candidates cover everything)
        total_coverage = sum(len(v) for v in coverage.values())
        avg_coverage = total_coverage / n_eff if n_eff > 0 else 0

        if n_eff <= 10:
            # Very few candidates - MILP is fast
            actual_solver = "milp"
            reason = f"few effective candidates ({n_eff})"
        elif problem_size > 10000:
            # Large problem - use greedy
            actual_solver = "greedy"
            reason = f"large problem size ({n_eff}×{n_blocked}={problem_size})"
        elif avg_coverage > n_blocked * 0.5:
            # High overlap - MILP can find better solution
            actual_solver = "milp"
            reason = f"high coverage overlap (avg {avg_coverage:.1f} per candidate)"
        else:
            # Default to greedy for safety
            actual_solver = "greedy"
            reason = "default for medium-sized problems"

        if verbose:
            print(f"\nDynamic solver selection: {actual_solver} ({reason})")

    if actual_solver == "greedy":
        if verbose:
            print("\nUsing greedy solver...")
        selected = solve_greedy(coverage, blocked, effective, verbose)
    else:
        if verbose:
            print("\nSolving MILP...")
        selected = solve_milp(coverage, blocked, effective, tradeoff_lambda, verbose)

    # Report results
    result = {
        "original_component_id": target_component_id,
        "original_component_size": comp_sizes.get(target_component_id, 0),
        "reactions_in_orig_comp": len(rxns_in_orig_comp),
        "blocked_in_phase2": len(blocked),
        "candidates_available": len(comp_candidates),
        "effective_candidates": len(effective),
        "selected": len(selected),
        "selected_details": selected,
    }

    if verbose:
        print(f"\n{'='*60}")
        print("RESULTS")
        print(f"{'='*60}")
        print(
            f"Original Component {target_component_id}: {comp_sizes.get(target_component_id, 0)} nodes"
        )
        print(f"Reactions from original component: {len(rxns_in_orig_comp)}")
        print(f"Still blocked in Phase-2 model: {len(blocked)}")
        print(f"Effective candidates: {len(effective)}")
        print(f"Selected PTRs: {len(selected)}")
        if selected:
            print("Selected PTR details:")
            for s in selected[:10]:
                print(f"  {s['met1']} <-> {s['met2']} (Type {s.get('type', '?')})")

    return result


def process_small_component(args):
    """Worker function for processing a single small component."""
    cid, size, cand_csv, p2_json, unc_json, params = args
    try:
        result = run_test_on_original_component(
            cand_csv,
            p2_json,
            unc_json,
            target_component_id=cid,
            tradeoff_lambda=params["lambda"],
            solver=params["solver"],
            solver_lp=params["solver_lp"],
            sample_blocked=params["sample_blocked"],
            parallel=False,  # No nested parallelization
            n_workers=None,
            use_cache=params.get("use_cache", True),
            cache_dir=params.get("cache_dir"),
            verbose=False,  # Reduce output in parallel mode
        )
        return (cid, size, result)
    except Exception as e:
        return (cid, size, {"error": str(e)})


def run_phase3_all_components(
    candidates_csv,
    phase2_model_json,
    unconnected_model_json,
    out_dir,
    phase3_model_out=None,
    tradeoff_lambda=0.01,
    solver="dynamic",
    solver_lp="glpk",
    parallel_fba=False,
    n_workers_fba=2,
    parallel_components=False,
    n_workers_components=4,
    small_component_threshold=500,
    min_component_size=4,
    component_ids=None,
    sample_blocked=None,
    use_cache=True,
    cache_dir=None,
    verbose=True,
):
    """
    Run Phase-3 on ALL original components using hybrid parallelization.

    HYBRID APPROACH:
    1. Small components (< threshold): Process in parallel (independent)
    2. Large components: Process sequentially, reusing PTRs from all previous

    Args:
        parallel_fba: Enable intra-component FBA parallelization
        n_workers_fba: Workers for FBA testing within each component
        parallel_components: Enable inter-component parallelization for small components
        n_workers_components: Workers for component processing
        small_component_threshold: Component size threshold for parallel vs sequential
        min_component_size: Skip components smaller than this
        component_ids: List of specific component IDs to process (None = all)

    Returns:
        dict with all_selected PTRs, per-component results, and summary stats
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed
    import json

    base_dir = os.path.dirname(__file__)
    if out_dir is None:
        out_dir = os.path.join(base_dir, "files")
    os.makedirs(out_dir, exist_ok=True)

    # Set LP solver
    if solver_lp:
        set_lp_solver(solver_lp)

    # Load models once to get component info
    if verbose:
        print("Loading models to analyze components...")
    unconnected_model = load_json_model(unconnected_model_json)

    # Build component mapping
    met_to_comp, rxn_to_comp, comp_sizes = build_original_component_map(
        unconnected_model
    )

    # Sort components by size (descending), apply filters
    components_to_process = []
    for cid, size in sorted(comp_sizes.items(), key=lambda x: x[1], reverse=True):
        # If specific component IDs provided, only process those
        if component_ids is not None:
            if cid not in component_ids:
                continue
        else:
            # Default filtering: skip main component (1) and tiny ones
            if cid == 1:
                continue
            if size < min_component_size:
                continue
        components_to_process.append((cid, size))

    if verbose:
        print(f"\nComponents to process: {len(components_to_process)}")
        if component_ids:
            print(f"  Filtered to: {component_ids}")
        else:
            print(
                f"  Skipped: component 1 (main), components < {min_component_size} nodes"
            )

        small = [c for c in components_to_process if c[1] < small_component_threshold]
        large = [c for c in components_to_process if c[1] >= small_component_threshold]
        print(f"  Small (< {small_component_threshold} nodes): {len(small)}")
        print(f"  Large (≥ {small_component_threshold} nodes): {len(large)}")

    # Separate small and large components
    small_components = [
        (cid, size)
        for cid, size in components_to_process
        if size < small_component_threshold
    ]
    large_components = [
        (cid, size)
        for cid, size in components_to_process
        if size >= small_component_threshold
    ]

    # Track all selected PTRs and results
    all_selected = []
    used_candidates = set()  # (met1, met2) pairs already selected
    component_results = {}

    start_time = time.time()

    # ==========================================================================
    # PHASE A: Process small components in parallel (independent)
    # ==========================================================================
    if small_components and parallel_components and n_workers_components > 1:
        if verbose:
            print(f"\n{'='*60}")
            print(
                f"PHASE A: Processing {len(small_components)} small components in parallel"
            )
            print(f"         Workers: {n_workers_components}")
            print(f"{'='*60}")

        # For parallel processing, we can't share state, so each component is independent
        # We'll merge results afterward

        # Prepare arguments for all small components
        params = {
            "lambda": tradeoff_lambda,
            "solver": solver,
            "solver_lp": solver_lp,
            "sample_blocked": sample_blocked,
            "use_cache": use_cache,
            "cache_dir": cache_dir,
        }

        small_args = [
            (
                cid,
                size,
                candidates_csv,
                phase2_model_json,
                unconnected_model_json,
                params,
            )
            for cid, size in small_components
        ]

        # Process in parallel
        with ProcessPoolExecutor(max_workers=n_workers_components) as executor:
            futures = {
                executor.submit(process_small_component, args): args[0]
                for args in small_args
            }

            completed = 0
            for future in as_completed(futures):
                cid = futures[future]
                completed += 1
                try:
                    result_cid, result_size, result = future.result()
                    # Use 'selected_details' (list) not 'selected' (count)
                    sel = result.get("selected_details", [])
                    if not isinstance(sel, (list, set)):
                        sel = []
                    blocked = result.get("blocked", 0)
                    if not isinstance(blocked, int):
                        blocked = 0
                    component_results[result_cid] = result
                    # Collect selected PTRs
                    for ptr in sel:
                        key = (ptr["met1"], ptr["met2"])
                        if key not in used_candidates:
                            used_candidates.add(key)
                            all_selected.append(ptr)
                    if verbose:
                        n_sel = (
                            result.get("selected", 0)
                            if isinstance(result.get("selected"), int)
                            else len(sel)
                        )
                        print(
                            f"  [{completed}/{len(small_components)}] Component {result_cid} "
                            f"({result_size} nodes): {n_sel} PTRs, {blocked} blocked"
                        )
                except Exception as e:
                    if verbose:
                        print(
                            f"  [{completed}/{len(small_components)}] Component {cid}: ERROR - {e}"
                        )
                    component_results[cid] = {"error": str(e)}

        if verbose:
            print(
                f"\nPhase A complete: {len(all_selected)} PTRs selected from small components"
            )

    elif small_components:
        # Sequential processing of small components (if parallel disabled)
        if verbose:
            print(f"\n{'='*60}")
            print(
                f"PHASE A: Processing {len(small_components)} small components sequentially"
            )
            print(f"{'='*60}")

        for i, (cid, size) in enumerate(small_components):
            if verbose:
                print(
                    f"\n[{i+1}/{len(small_components)}] Component {cid} ({size} nodes)"
                )

            result = run_test_on_original_component(
                candidates_csv,
                phase2_model_json,
                unconnected_model_json,
                target_component_id=cid,
                tradeoff_lambda=tradeoff_lambda,
                solver=solver,
                solver_lp=solver_lp,
                sample_blocked=sample_blocked,
                parallel=parallel_fba,
                n_workers=n_workers_fba,
                use_cache=use_cache,
                cache_dir=cache_dir,
                verbose=verbose,
            )
            component_results[cid] = result

            # Collect selected PTRs
            if "selected_details" in result and result["selected_details"]:
                for ptr in result["selected_details"]:
                    key = (ptr["met1"], ptr["met2"])
                    if key not in used_candidates:
                        used_candidates.add(key)
                        all_selected.append(ptr)

    # ==========================================================================
    # PHASE B: Process large components sequentially with PTR reuse
    # ==========================================================================
    if large_components:
        if verbose:
            print(f"\n{'='*60}")
            print(
                f"PHASE B: Processing {len(large_components)} large components sequentially"
            )
            print(f"         (with PTR reuse from previous components)")
            print(f"{'='*60}")

        # Reload candidates to filter out already-used ones
        all_candidates = load_candidates(candidates_csv)

        for i, (cid, size) in enumerate(large_components):
            if verbose:
                print(
                    f"\n[{i+1}/{len(large_components)}] Component {cid} ({size} nodes)"
                )
                print(f"  Previously selected PTRs: {len(all_selected)}")
                print(f"  Used candidates to skip: {len(used_candidates)}")

            # Get candidates for this component, excluding already-used ones
            comp_candidates = get_candidates_for_original_component(
                all_candidates, cid, used_candidates=used_candidates
            )

            if verbose:
                print(f"  Available candidates: {len(comp_candidates)}")

            if not comp_candidates:
                if verbose:
                    print(f"  No remaining candidates for this component")
                component_results[cid] = {
                    "blocked": 0,
                    "effective_candidates": 0,
                    "selected": [],
                    "skipped": "no_candidates",
                }
                continue

            # Run analysis on this component
            result = run_test_on_original_component(
                candidates_csv,
                phase2_model_json,
                unconnected_model_json,
                target_component_id=cid,
                tradeoff_lambda=tradeoff_lambda,
                solver=solver,
                solver_lp=solver_lp,
                sample_blocked=sample_blocked,
                parallel=parallel_fba,
                n_workers=n_workers_fba,
                use_cache=use_cache,
                cache_dir=cache_dir,
                verbose=verbose,
            )
            # Use 'selected_details' (list) not 'selected' (count)
            sel = result.get("selected_details", [])
            if not isinstance(sel, (list, set)):
                sel = []
            blocked = result.get("blocked", 0)
            if not isinstance(blocked, int):
                blocked = 0
            component_results[cid] = result

            # Collect newly selected PTRs
            new_ptrs = 0
            for ptr in sel:
                key = (ptr["met1"], ptr["met2"])
                if key not in used_candidates:
                    used_candidates.add(key)
                    all_selected.append(ptr)
                    new_ptrs += 1

            if verbose:
                n_sel = (
                    result.get("selected", 0)
                    if isinstance(result.get("selected"), int)
                    else len(sel)
                )
                print(f"  New PTRs selected: {new_ptrs} (total for component: {n_sel})")
                print(f"  Total PTRs so far: {len(all_selected)}")

    # ==========================================================================
    # SUMMARY
    # ==========================================================================
    elapsed = time.time() - start_time

    if verbose:
        print(f"\n{'='*60}")
        print("PHASE 3 COMPLETE - ALL COMPONENTS")
        print(f"{'='*60}")
        print(f"Total components processed: {len(component_results)}")
        print(f"Total PTRs selected: {len(all_selected)}")
        print(f"Total time: {elapsed:.1f}s")

        # Summary by component
        print(f"\nPer-component summary:")
        for cid in sorted(component_results.keys()):
            r = component_results[cid]
            if "error" in r:
                print(f"  Component {cid}: ERROR - {r['error']}")
            else:
                n_sel = (
                    r.get("selected", 0) if isinstance(r.get("selected"), int) else 0
                )
                blocked = (
                    r.get("blocked", 0) if isinstance(r.get("blocked"), int) else 0
                )
                print(f"  Component {cid}: {n_sel} PTRs, {blocked} blocked reactions")

    # Save results
    output = {
        "all_selected": all_selected,
        "component_results": {str(k): v for k, v in component_results.items()},
        "summary": {
            "total_components": len(component_results),
            "total_ptrs": len(all_selected),
            "elapsed_seconds": elapsed,
            "small_components": len(small_components),
            "large_components": len(large_components),
        },
    }

    # Save selected PTRs to CSV
    ptrs_csv = os.path.join(out_dir, "phase3_selected_ptrs.csv")
    with open(ptrs_csv, "w", newline="") as f:
        if all_selected:
            fieldnames = list(all_selected[0].keys())
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_selected)
    if verbose:
        print(f"\nSaved selected PTRs to: {ptrs_csv}")

    # Save full results as JSON
    results_json = os.path.join(out_dir, "phase3_results.json")

    # Convert any sets to lists for JSON serialization
    def convert_for_json(obj):
        if isinstance(obj, set):
            return list(obj)
        elif isinstance(obj, dict):
            return {k: convert_for_json(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_for_json(v) for v in obj]
        return obj

    with open(results_json, "w") as f:
        json.dump(convert_for_json(output), f, indent=2)
    if verbose:
        print(f"Saved full results to: {results_json}")

    # ==========================================================================
    # BUILD AND SAVE MODEL
    # ==========================================================================
    if verbose:
        print(f"\n{'='*60}")
        print("BUILDING FINAL MODEL")
        print(f"{'='*60}")

    # Load phase2 model (already connected)
    model = load_json_model(phase2_model_json)
    if verbose:
        print(f"Loaded phase2 model: {len(model.reactions)} reactions")

    # Add all selected PTRs to the model
    added_count = 0
    for idx, ptr in enumerate(all_selected):
        rid = add_transport(model, ptr, idx, prefix="SINK")
        if rid:
            ptr["reaction_id"] = rid
            added_count += 1
        else:
            if verbose:
                print(
                    f"  Warning: Could not add PTR {ptr.get('met1')} <-> {ptr.get('met2')}"
                )

    if verbose:
        print(f"Added {added_count} PTRs to model")
        print(f"Final model: {len(model.reactions)} reactions")

    # Save the model
    if phase3_model_out is None:
        phase3_model_out = os.path.join(out_dir, "phase3_sink_milp_model.json")

    save_json_model(model, phase3_model_out)
    if verbose:
        print(f"\nSaved final model to: {phase3_model_out}")

    # Add model path to output
    output["model_path"] = phase3_model_out

    return output


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Phase 3: Minimize blocked reactions with PTRs"
    )

    # Input files
    parser.add_argument("--candidates", default=None, help="Phase-1 candidates CSV")
    parser.add_argument(
        "--phase2-model", default=None, help="Phase-2 connected model JSON"
    )
    parser.add_argument(
        "--unconnected-model", default=None, help="Original unconnected model JSON"
    )

    # Mode selection
    parser.add_argument(
        "--all-components",
        action="store_true",
        help="Process ALL components (hybrid parallelization)",
    )
    parser.add_argument(
        "--components",
        type=str,
        default=None,
        help='Comma-separated component IDs to process (e.g., "5,6,7")',
    )
    parser.add_argument(
        "--component",
        type=int,
        default=None,
        help="Single component ID to process (1=main, 2=largest isolated, etc.)",
    )

    # Algorithm parameters
    parser.add_argument(
        "--max-candidates",
        type=int,
        default=None,
        help="Limit candidates for quick testing (single component mode)",
    )
    parser.add_argument(
        "--lambda",
        dest="tradeoff_lambda",
        type=float,
        default=0.01,
        help="Tradeoff: coverage vs PTR count",
    )
    parser.add_argument(
        "--solver",
        choices=["milp", "greedy", "dynamic"],
        default="dynamic",
        help="Solver: milp (optimal), greedy (fast), dynamic (auto)",
    )
    parser.add_argument(
        "--solver-lp",
        choices=["glpk", "glpk_exact", "scipy", "gurobi", "cplex"],
        default="glpk",
        help="LP solver for FBA: glpk (default), glpk_exact, scipy",
    )
    parser.add_argument(
        "--sample-blocked",
        type=int,
        default=None,
        help="Sample N blocked reactions (for large components)",
    )

    # Intra-component parallelization (FBA testing)
    parser.add_argument(
        "--parallel-fba",
        action="store_true",
        help="Parallel FBA testing within each component",
    )
    parser.add_argument(
        "--workers-fba",
        type=int,
        default=2,
        help="Workers for FBA testing (default: 2)",
    )

    # Inter-component parallelization (hybrid approach)
    parser.add_argument(
        "--parallel-components",
        action="store_true",
        help="Parallel processing of small components",
    )
    parser.add_argument(
        "--workers-components",
        type=int,
        default=4,
        help="Workers for component processing (default: 4)",
    )
    parser.add_argument(
        "--small-threshold",
        type=int,
        default=500,
        help="Component size threshold for parallel vs sequential (default: 500)",
    )
    parser.add_argument(
        "--min-size",
        type=int,
        default=4,
        help="Skip components smaller than this (default: 4)",
    )

    # Legacy single-component options (for backward compatibility)
    parser.add_argument(
        "--parallel", action="store_true", help="(Legacy) Same as --parallel-fba"
    )
    parser.add_argument(
        "--workers", type=int, default=None, help="(Legacy) Same as --workers-fba"
    )

    args = parser.parse_args()

    # Handle legacy arguments
    if args.parallel:
        args.parallel_fba = True
    if args.workers:
        args.workers_fba = args.workers

    base_dir = os.path.dirname(__file__)
    cand_csv = args.candidates or os.path.join(base_dir, "files", "candidates_all.csv")
    phase2_json = args.phase2_model or os.path.normpath(
        os.path.join(
            base_dir,
            "..",
            "models",
            "base",
            "THG-beta-batch_251106_phase2_minimal_connected.json",
        )
    )
    unconnected_json = args.unconnected_model or os.path.normpath(
        os.path.join(
            base_dir, "..", "models", "base", "THG-beta-batch_251106_unconnected.json"
        )
    )

    for f, name in [
        (cand_csv, "Candidates"),
        (phase2_json, "Phase2 model"),
        (unconnected_json, "Unconnected model"),
    ]:
        if not os.path.exists(f):
            print(f"{name} not found: {f}")
            exit(1)

    # Parse component IDs if provided
    component_ids = None
    if args.components:
        component_ids = [int(c.strip()) for c in args.components.split(",")]

    if args.all_components or args.components:
        # Process multiple components with hybrid parallelization
        run_phase3_all_components(
            cand_csv,
            phase2_json,
            unconnected_json,
            tradeoff_lambda=args.tradeoff_lambda,
            solver=args.solver,
            solver_lp=args.solver_lp,
            sample_blocked=args.sample_blocked,
            parallel_fba=args.parallel_fba,
            n_workers_fba=args.workers_fba,
            parallel_components=args.parallel_components,
            n_workers_components=args.workers_components,
            small_component_threshold=args.small_threshold,
            min_component_size=args.min_size,
            component_ids=component_ids,
            verbose=True,
        )
    else:
        # Single component mode (default: component 6 for testing)
        component = args.component if args.component is not None else 6
        run_test_on_original_component(
            cand_csv,
            phase2_json,
            unconnected_json,
            target_component_id=component,
            max_candidates=args.max_candidates,
            tradeoff_lambda=args.tradeoff_lambda,
            solver=args.solver,
            solver_lp=args.solver_lp,
            sample_blocked=args.sample_blocked,
            parallel=args.parallel_fba,
            n_workers=args.workers_fba,
        )
