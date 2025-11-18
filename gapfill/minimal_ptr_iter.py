"""
PFBA-based approach to identify minimal transporters from a candidate list.
This script performs, for every reaction (except candidate transporters), two PFBA runs
(maximize and minimize that reaction) in batches and parallel, collects the fluxes
of candidate transport reactions across all solutions, and reports which candidate
reactions are never active within a specified tolerance.

Usage (example at bottom): adjust `model_name` and `candidate_csv` as needed.
"""

import os
import csv
from typing import Set, List, Tuple, Any
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm
import time

import cobra
from cobra.flux_analysis import pfba
from optlang.symbolics import Zero


NUM_CORES = 18
DEBUG = False
# Number of sequential iterations to run the whole PFBA scanning procedure.
# Each iteration will run find_minimal_transporters_pfba using the set of
# candidate transporters returned as "required" from the previous iteration.
# Set to 1 for the original (single-run) behaviour.
NUM_ITERATIONS = 2
# If True, stop early when the set of required transporters doesn't change
# between iterations.
STOP_ON_CONVERGENCE = True


def fix_bounds_pfba_pair(
    cobra_model,
    candidate_indices,
    obj_rxn_ix,
    lb,
    ub,
    tol=None,
    debug_timing=False,
):
    """Perform both max and min parsimonious FBA optimizations efficiently.

    Sets up the parsimonious model once and reuses it for both max and min
    optimizations to save computation time.

    Parameters
    ----------
    cobra_model : cobra.Model
        The COBRA model to be used.
    candidate_indices : list
        List of indices for candidate reactions to include in parsimonious objective.
    obj_rxn_ix : int
        Index of the objective reaction.
    lb : array_like
        Lower bounds for all reactions in the model.
    ub : array_like
        Upper bounds for all reactions in the model.
    tol : float, optional
        Tolerance for bound adjustments. If None, model's tolerance is used.
        Defaults to None.
    show_progress : bool, optional
        Whether to print detailed timing information. Defaults to False.

    Returns
    -------
    tuple
        Tuple of (max_fluxes, min_fluxes) as numpy arrays.

    """
    if tol is None:
        tol = cobra_model.tolerance

    with cobra_model as m:
        # First optimization: maximize
        if debug_timing:
            fba_start = time.time()
            print(f"  FBA: Starting max optimization of objective")

        m.objective = m.reactions[obj_rxn_ix]
        opt_max = m.optimize("max")
        solution_max = opt_max.fluxes.to_numpy()

        if debug_timing:
            fba_time = time.time() - fba_start
            print(f"  FBA: max optimization completed in {fba_time:.4f} seconds")

        # Second optimization: minimize
        if debug_timing:
            fba_start = time.time()
            print(f"  FBA: Starting min optimization of objective")

        opt_min = m.optimize("min")
        solution_min = opt_min.fluxes.to_numpy()

        if debug_timing:
            fba_time = time.time() - fba_start
            print(f"  FBA: min optimization completed in {fba_time:.4f} seconds")

        # Setup parsimonious model once
        if debug_timing:
            setup_start = time.time()
            print(f"  pFBA: Setting up parsimonious optimization (shared)")

        m.objective = m.solver.interface.Objective(Zero, direction="min", sloppy=False)

        # Pre-calculate all variables and bounds for both solutions
        max_objective_vars = []
        min_objective_vars = []

        # Process candidate reactions for both solutions
        for ix in candidate_indices:
            try:
                # For max solution
                flux_val_max = solution_max[ix]
                if abs(flux_val_max) >= tol:
                    if flux_val_max >= 0:
                        max_objective_vars.append(m.reactions[ix].forward_variable)
                    else:
                        max_objective_vars.append(m.reactions[ix].reverse_variable)

                # For min solution
                flux_val_min = solution_min[ix]
                if abs(flux_val_min) >= tol:
                    if flux_val_min >= 0:
                        min_objective_vars.append(m.reactions[ix].forward_variable)
                    else:
                        min_objective_vars.append(m.reactions[ix].reverse_variable)

            except Exception as e:
                if debug_timing:
                    print(f"Warning: Issue with reaction {ix}: {e}")

        if debug_timing:
            setup_time = time.time() - setup_start
            print(f"  pFBA: Setup completed in {setup_time:.4f} seconds")

        # Run parsimonious optimization for max solution
        if debug_timing:
            parsimonious_start = time.time()
            print(f"  pFBA: Running parsimonious optimization (max)")

        # Set bounds for max solution
        obj_flux_max = solution_max[obj_rxn_ix]
        obj_lower_bound = float(
            min(max(obj_flux_max - tol, lb[obj_rxn_ix]), ub[obj_rxn_ix])
        )
        obj_upper_bound = float(
            max(min(obj_flux_max + tol, ub[obj_rxn_ix]), lb[obj_rxn_ix])
        )
        m.reactions[obj_rxn_ix].bounds = (obj_lower_bound, obj_upper_bound)

        for ix in candidate_indices:
            flux_val = solution_max[ix]
            if abs(flux_val) < tol:
                m.reactions[ix].bounds = (0, 0)
            elif flux_val >= 0:
                m.reactions[ix].bounds = (max(0, lb[ix]), min(flux_val + tol, ub[ix]))
            else:
                m.reactions[ix].bounds = (max(flux_val - tol, lb[ix]), min(0, ub[ix]))

        if max_objective_vars:
            m.objective.set_linear_coefficients({v: 1.0 for v in max_objective_vars})

        opt_max_pfba = m.optimize()
        fluxes_max = opt_max_pfba.fluxes.to_numpy()

        if debug_timing:
            parsimonious_time = time.time() - parsimonious_start
            print(
                f"  pFBA: Max parsimonious optimization completed in {parsimonious_time:.4f} seconds ({len(max_objective_vars)} variables)"
            )

        # Run parsimonious optimization for min solution
        if debug_timing:
            parsimonious_start = time.time()
            print(f"  pFBA: Running parsimonious optimization (min)")

        # Reset objective coefficients
        if max_objective_vars:
            m.objective.set_linear_coefficients({v: 0.0 for v in max_objective_vars})

        # Set bounds for min solution
        obj_flux_min = solution_min[obj_rxn_ix]
        obj_lower_bound = float(
            min(max(obj_flux_min - tol, lb[obj_rxn_ix]), ub[obj_rxn_ix])
        )
        obj_upper_bound = float(
            max(min(obj_flux_min + tol, ub[obj_rxn_ix]), lb[obj_rxn_ix])
        )
        m.reactions[obj_rxn_ix].bounds = (obj_lower_bound, obj_upper_bound)

        for ix in candidate_indices:
            flux_val = solution_min[ix]
            if abs(flux_val) < tol:
                m.reactions[ix].bounds = (0, 0)
            elif flux_val >= 0:
                m.reactions[ix].bounds = (max(0, lb[ix]), min(flux_val + tol, ub[ix]))
            else:
                m.reactions[ix].bounds = (max(flux_val - tol, lb[ix]), min(0, ub[ix]))

        if min_objective_vars:
            m.objective.set_linear_coefficients({v: 1.0 for v in min_objective_vars})

        opt_min_pfba = m.optimize()
        fluxes_min = opt_min_pfba.fluxes.to_numpy()

        if debug_timing:
            parsimonious_time = time.time() - parsimonious_start
            print(
                f"  pFBA: Min parsimonious optimization completed in {parsimonious_time:.4f} seconds ({len(min_objective_vars)} variables)"
            )

    return fluxes_max, fluxes_min


def find_minimal_transporters_pfba(
    cobra_model: Any,
    candidate_transporters: Set[str],
    clean: bool = True,
    n_jobs: int = 1,
    batch_size: int = None,
    tol: float = 1e-9,
    model_path: str = None,
    save_fluxes_path: str = "all_fluxes.npy",
) -> Tuple[Set[str], np.ndarray, List[str]]:
    """
    PFBA-based scanning approach.

    For every reaction in the model excluding candidate_transporters, perform two
    PFBA optimizations:
      - maximize the reaction flux
      - minimize the reaction flux

    Collect the fluxes of the candidate transporters from each PFBA solution into
    a matrix with shape (num_solutions, num_candidate_transporters).

    Args:
        cobra_model: a loaded cobrapy model
        candidate_transporters: set of reaction IDs considered as candidate transporters
        clean: if True, consider all reactions as targets; otherwise exclude blocked reactions
        n_jobs: parallel workers
        batch_size: number of target reactions per batch (default: ceil(num_targets / n_jobs))
        tol: tolerance below which a flux is considered inactive
        model_path: path to the model file for parallel processing
        save_fluxes_path: optional path to save flux matrix as .npy file

    Returns:
        (required_transporters, flux_matrix, candidate_list)
          required_transporters: set of candidate reaction IDs that are active in at least one solution
          flux_matrix: numpy array (num_solutions x num_candidates) with candidate fluxes
          candidate_list: ordered list of candidate reaction IDs corresponding to columns
    """
    if cobra is None or pfba is None:
        raise RuntimeError(
            "`cobra` and its `pfba` function are required to run this script."
        )

    reaction_ids = [r.id for r in cobra_model.reactions]

    # Determine target reactions to scan
    if clean:
        target_ids = [rid for rid in reaction_ids if rid not in candidate_transporters]
    else:
        blocked = set(cobra.flux_analysis.find_blocked_reactions(cobra_model))
        target_ids = [
            rid
            for rid in reaction_ids
            if rid not in candidate_transporters and rid not in blocked
        ]

    print(f"Total reactions: {len(reaction_ids)}; targets to scan: {len(target_ids)}")

    # Create an index map for candidate reactions to extract fluxes
    rxn_index = {rid: i for i, rid in enumerate(reaction_ids)}

    # Filter candidate list to only reactions present in the model to avoid KeyError
    missing_candidates = [rid for rid in candidate_transporters if rid not in rxn_index]
    if missing_candidates:
        print(
            f"Warning: the following candidate reactions are not present in the model and will be skipped: {missing_candidates}"
        )

    candidate_list = sorted([rid for rid in candidate_transporters if rid in rxn_index])
    num_candidates = len(candidate_list)
    if num_candidates == 0:
        print(
            "Warning: no candidate_transporters are present in the model. Returning empty results."
        )
        return set(), np.zeros((0, 0)), candidate_list

    candidate_indices = [rxn_index[rid] for rid in candidate_list]

    # Prepare batches
    if n_jobs < 1:
        n_jobs = 1
    if batch_size is None:
        import math

        batch_size = math.ceil(len(target_ids) / n_jobs)
    num_batches = int(np.ceil(len(target_ids) / batch_size))
    batches = [
        target_ids[i * batch_size : (i + 1) * batch_size] for i in range(num_batches)
    ]

    print(f"Running PFBA in {num_batches} batches with up to {n_jobs} parallel jobs...")

    # Use a process-based backend and have each worker load the model from disk to
    # avoid pickling heavy cobra.Model objects. The worker function is defined at
    # module level below as `_solve_batch`.

    # Run in parallel using loky (processes). Each worker will load the model from `model_path`.
    if n_jobs < 1:
        n_jobs = 1
    if model_path is None and n_jobs > 1:
        print(
            "Warning: running with multiple processes but no `model_path` provided; falling back to copying/ pickling the model which may be slow or fail."
        )

    parallel = Parallel(n_jobs=n_jobs, backend="loky")
    job_results = parallel(
        delayed(_solve_batch)(
            batch,
            candidate_list,
            num_candidates,
            model_path,
            tol,
            i == 0,  # Always show progress for first worker
            DEBUG,  # Pass DEBUG flag for timing
        )
        for i, batch in enumerate(tqdm(batches))
    )

    # Stack all results
    if len(job_results) == 0:
        all_fluxes = np.zeros((0, num_candidates))
    else:
        all_fluxes = np.vstack(job_results)

    # Identify candidate reactions that are active in any solution
    active_mask = np.any(np.abs(all_fluxes) > tol, axis=0)
    required_transporters = {candidate_list[i] for i, a in enumerate(active_mask) if a}

    inactive_candidates = [
        candidate_list[i] for i, a in enumerate(active_mask) if not a
    ]

    print(f"Total PFBA solutions: {all_fluxes.shape[0]}")
    print(
        f"Candidates active in at least one solution: {len(required_transporters)}/{num_candidates}"
    )

    # Save flux matrix if requested
    if save_fluxes_path is not None:
        np.save(save_fluxes_path, all_fluxes)
        print(f"Saved flux matrix to {save_fluxes_path}")

    return required_transporters, all_fluxes, candidate_list


def _solve_batch(
    batch: List[str],
    candidate_list: List[str],
    num_candidates: int,
    model_path: str,
    tol: float,
    show_progress: bool = False,
    debug_timing: bool = False,
) -> np.ndarray:
    """Worker function executed in a separate process.

    Loads the model from `model_path` if provided, otherwise raises an error.
    Runs custom pFBA maximize/minimize for each reaction ID in `batch` using
    fix_bounds_pfba and returns a (2*len(batch), num_candidates) array with
    candidate reaction fluxes.
    """
    if model_path is None:
        raise RuntimeError(
            "model_path must be provided when running workers in separate processes"
        )

    # Load the model inside the worker process to avoid pickling
    from cobra.io import read_sbml_model, load_json_model

    if debug_timing:
        start_time = time.time()
        print(f"Worker: Loading model from {model_path}")

    # local_model = read_sbml_model(model_path)
    local_model = load_json_model(model_path)

    if debug_timing:
        load_time = time.time() - start_time
        print(f"Worker: Model loaded in {load_time:.2f} seconds")

    local_reaction_ids = [r.id for r in local_model.reactions]
    # sanity: ensure all candidate_list ids exist in local model
    missing = [rid for rid in candidate_list if rid not in local_reaction_ids]
    if missing:
        # Return zeros for missing candidates but warn
        print(f"Warning: missing candidate reactions in local model: {missing}")

    # Get reaction indices and bounds
    rxn_index = {rid: i for i, rid in enumerate(local_reaction_ids)}
    candidate_indices = [rxn_index[rid] for rid in candidate_list if rid in rxn_index]

    # Get bounds for all reactions
    lb = np.array([r.lower_bound for r in local_model.reactions])
    ub = np.array([r.upper_bound for r in local_model.reactions])

    results = []
    for tid in tqdm(batch, desc="Processing reactions", disable=not show_progress):
        if debug_timing:
            reaction_start = time.time()

        # Maximize target
        try:
            target_rxn = local_model.reactions.get_by_id(tid)
            obj_rxn_ix = rxn_index[tid]
        except KeyError:
            # reaction not found in local model - append zeros for both max/min
            results.append(np.zeros(num_candidates))
            results.append(np.zeros(num_candidates))
            continue

        # Use the new paired optimization function for efficiency
        fluxes_max, fluxes_min = fix_bounds_pfba_pair(
            local_model,
            candidate_indices,
            obj_rxn_ix,
            lb,
            ub,
            tol,
            debug_timing,
        )

        candidate_fluxes_max = np.array(
            [
                fluxes_max[rxn_index[rid]] if rid in rxn_index else 0.0
                for rid in candidate_list
            ],
            dtype=float,
        )
        results.append(candidate_fluxes_max)

        candidate_fluxes_min = np.array(
            [
                fluxes_min[rxn_index[rid]] if rid in rxn_index else 0.0
                for rid in candidate_list
            ],
            dtype=float,
        )
        results.append(candidate_fluxes_min)

        if debug_timing:
            reaction_time = time.time() - reaction_start
            print(f"Worker: Processed {tid} in {reaction_time:.3f} seconds")

    return np.vstack(results) if results else np.zeros((0, num_candidates))


# --- Example usage when executed as a script ---
# This script analyzes the transport reactions created by connect_components.py
# to identify which ones are truly required for connecting disconnected components.
if __name__ == "__main__":
    if cobra is None:
        print("This script requires cobra. Please install cobrapy to run it.")
        raise SystemExit(1)

    from cobra.io import read_sbml_model, load_json_model

    model_name = "THG-beta-batch_251106"

    # Load the model directly from models/{model_name}.json (preferred), or
    # from models/{model_name}.xml and convert to JSON if only SBML is available.
    model_json = f"models/{model_name}.json"
    model_xml = f"models/{model_name}.xml"

    if os.path.exists(model_json):
        print(f"Loading model JSON from {model_json}")
        model = load_json_model(model_json)
        model_path = model_json
    elif os.path.exists(model_xml):
        print(
            f"Loading SBML model from {model_xml} and converting to JSON at {model_json}"
        )
        model = read_sbml_model(model_xml)
        cobra.io.save_json_model(model, model_json)
        model_path = model_json
    else:
        print(f"Model file {model_json} or {model_xml} not found. Exiting.")
        raise SystemExit(1)

    # Load candidate reactions from connect_components.py output (required)
    candidate_csv = f"models/{model_name}_created_reactions.csv"
    if not os.path.exists(candidate_csv):
        print(f"Candidate file {candidate_csv} not found. Exiting.")
        raise SystemExit(1)

    # Load initial candidate reactions
    current_candidates = set()
    with open(candidate_csv, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Get reaction ID from the appropriate column
            if "reaction_id" in row:
                rid = row["reaction_id"].strip()
            else:
                # Fallback for old format (just reaction IDs)
                rid = row[0].strip() if row else ""
            if rid:
                current_candidates.add(rid)

    if len(current_candidates) == 0:
        print("No candidate reactions found in candidate CSV. Exiting.")
        raise SystemExit(1)

    previous_candidates = None

    # Run multiple iterations sequentially, updating the candidate set to the
    # required transporters found in the previous iteration. This allows
    # progressive pruning until convergence or until NUM_ITERATIONS is reached.
    for iteration in range(1, max(1, NUM_ITERATIONS) + 1):
        print(f"\n=== Iteration {iteration}/{NUM_ITERATIONS} ===")
        print(f"Candidates this iteration: {len(current_candidates)}")

        required_transporters, flux_matrix, candidate_list = (
            find_minimal_transporters_pfba(
                cobra_model=model,
                candidate_transporters=current_candidates,
                clean=False,
                n_jobs=NUM_CORES,
                batch_size=None,
                tol=1e-9,
                model_path=model_path,
                save_fluxes_path=f"flux_matrix_{model_name}_iter{iteration}.npy",
            )
        )

        print("--- Results (iteration {}) ---".format(iteration))
        print(f"Required transporters: {len(required_transporters)}")
        for rid in sorted(required_transporters):
            print(rid)

        out_csv = f"minimal_transporters_pfba_{model_name}_iter{iteration}.csv"
        with open(out_csv, "w", newline="") as cf:
            writer = csv.writer(cf)
            writer.writerow(["candidate_reaction", "active_in_any_solution"])
            for i, rid in enumerate(candidate_list):
                writer.writerow([rid, int(rid in required_transporters)])

        print(f"Wrote results to {out_csv}")

        # Prepare for next iteration
        previous_candidates = set(current_candidates)
        current_candidates = set(required_transporters)

        # If no candidates remain, stop early
        if not current_candidates:
            print("No required transporters remain after this iteration; stopping.")
            break

        # Optional convergence check
        if STOP_ON_CONVERGENCE and previous_candidates == current_candidates:
            print(
                "Candidate set unchanged from previous iteration; converged. Stopping."
            )
            break

    print("\nAll iterations complete. Final required transporters:")
    for rid in sorted(current_candidates):
        print(rid)
