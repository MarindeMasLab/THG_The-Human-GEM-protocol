"""Main gapfill CLI to run Phase1/Phase2/Phase3 with different objectives.

Usage examples:
  # Individual phases
    python3 gapfill/gapfill.py phase1 --model <unconnected_model.json>
    python3 gapfill/gapfill.py phase2 --mode minimal
    python3 gapfill/gapfill.py phase3 --mode deadends
    python3 gapfill/gapfill.py phase3 --mode blocked --max 200

  # Full pipeline - deadends objective
    python3 gapfill/gapfill.py run-all --pipeline deadends --model <unconnected_model.json>

    # Full pipeline - blocked reactions objective (choose Phase 3 mode):
  
    # Option 1: component_milp - Single MILP per component (all types together, optimal)
    python3 gapfill/gapfill.py run-all --pipeline blocked --phase3-mode component_milp --model <unconnected_model.json>
  
    # Option 2: tiered_milp - Tiered MILP per component (Type A → B → C, prioritizes A)
    python3 gapfill/gapfill.py run-all --pipeline blocked --phase3-mode tiered_milp --model <unconnected_model.json>
  
    # Option 3: sink_milp - Temp sinks + original components + hybrid parallelization (RECOMMENDED)
    python3 gapfill/gapfill.py run-all --pipeline blocked --phase3-mode sink_milp --components 4,5,6,7 --model <unconnected_model.json>

    # Option 4: sink_milp_original - Legacy sink MILP implementation (with model export)
    python3 gapfill/gapfill.py run-all --pipeline blocked --phase3-mode sink_milp_original --components 4,5,6,7 --model <unconnected_model.json> \
                        --phase3-model gapfill/files/phase3_sink_test.json
  
  # Limit to specific components (faster testing):
    python3 gapfill/gapfill.py run-all --pipeline blocked --strategy component_milp --min-comp-size 200 --max-components 3 --model <unconnected_model.json>

Phase3 strategies for blocked pipeline:
  - exact: Greedy per-candidate (slow, not recommended)
  - hybrid: Greedy with top-K prefiltering
  - hybrid_batch: Batch greedy with checkpointing
  - milp: Global MILP (very slow for large models)
  - component_milp: MILP per component, all types
  - tiered_milp: MILP per component, A→B→C tiers
  - sink_milp: Temp sinks + original components + hybrid parallelization (RECOMMENDED)
"""
import argparse
import os
import sys
import importlib
from cobra import Configuration


def main():
    parser = argparse.ArgumentParser(description='Gapfill pipeline orchestrator')
    sub = parser.add_subparsers(dest='cmd')

    p1 = sub.add_parser('phase1', help='Run phase1 candidate generation')
    p1.add_argument('--model', required=True, help='Input model JSON for phase1')
    p1.add_argument('--out', help='Output directory for phase1 files')

    p2 = sub.add_parser('phase2', help='Run phase2 connector selection')
    p2.add_argument('--mode', choices=['minimal','prioritized'], default='minimal')
    p2.add_argument('--candidates', help='Phase1 candidates CSV')
    p2.add_argument('--model', help='Unconnected model JSON')
    p2.add_argument('--out', help='Output dir for phase2 files')

    p3 = sub.add_parser('phase3', help='Run phase3 optimizer')
    p3.add_argument('--mode', choices=['deadends','blocked'], default='deadends')
    p3.add_argument('--candidates', help='Phase1 candidates CSV')
    p3.add_argument('--model', help='Starting model JSON (phase2 output)')
    p3.add_argument('--out', help='Output dir for phase3 files')
    p3.add_argument('--max', type=int, default=500, help='Max additions for phase3')
    p3.add_argument('--solver-lp', choices=['glpk','glpk_exact','scipy','gurobi','cplex','hybrid','highs'], default=None, help='LP solver for Phase-3 run (deadends/blocked)')

    pall = sub.add_parser('run-all', help='Run full pipeline: phase1 -> phase2 -> phase3')
    pall.add_argument('--pipeline', choices=['deadends','blocked'], default='deadends', help='Objective pipeline')
    pall.add_argument('--candidates', help='Phase1 candidates CSV (skip phase1 if exists)')
    pall.add_argument('--model', required=True, help='Unconnected model JSON (phase1 input)')
    pall.add_argument('--out', help='Output dir for phase files')
    pall.add_argument('--max', type=int, default=500, help='Max additions for phase3')
    pall.add_argument('--force-phase1', action='store_true', help='Force re-run of phase1 even if candidates CSV exists')
    pall.add_argument('--phase2-mode', choices=['minimal','prioritized'], default='minimal', help='Phase2 connector strategy')
    pall.add_argument('--phase3-mode', choices=['exact','hybrid','hybrid_batch','milp','component_milp','tiered_milp','sink_milp','sink_milp_original'], default=None, 
                     help='Phase3 strategy: exact=greedy, milp=global MILP, component_milp=per-component MILP, tiered_milp=A→B→C, sink_milp=temp sinks+hybrid parallel, sink_milp_original=legacy sink MILP with model export')
    # Deprecated alias for backward compatibility
    pall.add_argument('--strategy', choices=['exact','hybrid','hybrid_batch','milp','component_milp','tiered_milp','sink_milp','sink_milp_original'], default='exact', help=argparse.SUPPRESS)
    pall.add_argument('--topk', type=int, default=500, help='Top-K candidates to prefilter for hybrid strategies')
    pall.add_argument('--batch-size', type=int, default=50, help='Batch size for hybrid_batch strategy')
    pall.add_argument('--lambda', dest='tradeoff_lambda', type=float, default=0.01, help='MILP tradeoff weight (higher = fewer PTRs)')
    pall.add_argument('--max-iter', type=int, default=10, help='Max iterations for MILP strategy')
    pall.add_argument('--min-comp-size', type=int, default=4, help='Minimum component size for component_milp strategy')
    pall.add_argument('--max-components', type=int, default=None, help='Max number of components to process (component_milp)')
    # sink_milp specific options
    pall.add_argument('--components', type=str, default=None, help='Comma-separated component IDs for sink_milp (e.g., "4,5,6,7")')
    pall.add_argument('--parallel-components', action='store_true', help='Enable parallel processing of small components (sink_milp)')
    pall.add_argument('--workers-components', type=int, default=4, help='Workers for parallel component processing (sink_milp)')
    pall.add_argument('--parallel-fba', action='store_true', help='Enable parallel FBA testing within components (sink_milp)')
    pall.add_argument('--workers-fba', type=int, default=2, help='Workers for parallel FBA testing (sink_milp)')
    pall.add_argument('--small-threshold', type=int, default=500, help='Component size threshold for parallel vs sequential (sink_milp)')
    pall.add_argument('--phase3-model', dest='phase3_model_out', default=None, help='Output path for Phase-3 model (only for sink_milp_original)')
    pall.add_argument('--solver-lp', choices=['glpk','glpk_exact','scipy','gurobi','cplex','hybrid','highs'], default='glpk', help='LP solver for FBA (sink_milp). gurobi/cplex require licenses. hybrid uses HiGHS via optlang. highs uses direct HiGHS (~4x faster than glpk!).')

    args = parser.parse_args()

    if args.cmd == 'phase1':
        # ensure local package directory is on path so imports work when running script directly
        script_dir = os.path.dirname(__file__)
        if script_dir not in sys.path:
            sys.path.insert(0, script_dir)
        phase1_connect_components = importlib.import_module('phase1_connect_components')
        model = os.path.normpath(args.model)
        if not os.path.exists(model):
            print('Phase1 model not found:', model); sys.exit(1)
        from cobra.io import load_json_model as _load
        mc = _load(model)
        out = args.out
        cand, dead, comp = phase1_connect_components.generate_phase1_outputs(mc, out_dir=out)
        print('Phase1 outputs:', cand, dead, comp)
        return

    if args.cmd == 'phase2':
        cand = args.candidates or os.path.join(os.path.dirname(__file__), 'files', 'candidates_all.csv')
        cand = os.path.normpath(cand)
        if not os.path.exists(cand):
            print('Candidates CSV not found:', cand); sys.exit(1)
        if args.mode == 'minimal':
            sys.path.insert(0, os.path.dirname(__file__))
            phase2_minimal = importlib.import_module('phase2_minimal_connector')
            run_phase2 = phase2_minimal.run_phase2
            out_csv, model_out, n = run_phase2(cand, unconnected_model_json=args.model, out_dir=args.out)
            print('Phase2 minimal wrote', out_csv, 'model', model_out, 'n=', n)
        else:
            sys.path.insert(0, os.path.dirname(__file__))
            phase2_prior = importlib.import_module('phase2_prioritized_connector')
            run_phase2 = phase2_prior.run_phase2
            out_csv, model_out, n = run_phase2(cand, unconnected_model_json=args.model, out_dir=args.out)
            print('Phase2 prioritized wrote', out_csv, 'model', model_out, 'n=', n)
        return

    if args.cmd == 'phase3':
        # Set LP solver globally if provided
        if args.solver_lp:
            try:
                Configuration().solver = args.solver_lp
                print(f"LP solver set to: {args.solver_lp}")
            except Exception as e:
                print(f"Warning: could not set solver '{args.solver_lp}': {e}")

        cand = args.candidates or os.path.join(os.path.dirname(__file__), 'files', 'candidates_all.csv')
        cand = os.path.normpath(cand)
        if not os.path.exists(cand):
            print('Candidates CSV not found:', cand); sys.exit(1)
        if args.mode == 'deadends':
            sys.path.insert(0, os.path.dirname(__file__))
            ph3_dead = importlib.import_module('phase3_greedy_optimizer')
            run_phase3 = ph3_dead.run_phase3
            out_csv, model_out, n = run_phase3(cand, phase2_model_json=args.model, out_dir=args.out, max_additions=args.max)
            print('Phase3 deadends wrote', out_csv, 'model', model_out, 'n=', n)
        else:
            sys.path.insert(0, os.path.dirname(__file__))
            ph3_block = importlib.import_module('phase3_blocked_optimizer')
            run_phase3 = ph3_block.run_phase3
            out_csv, model_out, n = run_phase3(cand, starting_model_json=args.model, out_dir=args.out, max_additions=args.max)
            print('Phase3 blocked wrote', out_csv, 'model', model_out, 'n=', n)
        return

    if args.cmd == 'run-all':
        # Ensure local modules importable
        script_dir = os.path.dirname(__file__)
        if script_dir not in sys.path:
            sys.path.insert(0, script_dir)

        # Set LP solver globally if provided (applies to Phase-3 steps)
        if args.solver_lp:
            try:
                Configuration().solver = args.solver_lp
                print(f"LP solver set to: {args.solver_lp}")
            except Exception as e:
                print(f"Warning: could not set solver '{args.solver_lp}': {e}")

        cand = args.candidates or os.path.join(os.path.dirname(__file__), 'files', 'candidates_all.csv')
        cand = os.path.normpath(cand)
        model = os.path.normpath(args.model)
        out_dir = args.out

        # Phase1: generate candidates if missing or forced
        if args.force_phase1 or not os.path.exists(cand):
            print('Running Phase-1 candidate generation...')
            phase1 = importlib.import_module('phase1_connect_components')
            from cobra.io import load_json_model as _load
            mc = _load(model)
            cand_path, dead, comps = phase1.generate_phase1_outputs(mc, out_dir=out_dir)
            cand = cand_path
            print('Phase-1 done ->', cand)
        else:
            print('Found existing candidates CSV:', cand)

        # Phase2: connector selection (minimal or prioritized)
        if args.phase2_mode == 'minimal':
            print('Running Phase-2 minimal connector selection...')
            phase2_min = importlib.import_module('phase2_minimal_connector')
            out_csv2, phase2_model, n2 = phase2_min.run_phase2(cand, unconnected_model_json=model, out_dir=out_dir)
            print(f'Phase-2 minimal: selected {n2} connectors -> model {phase2_model}')
        else:
            print('Running Phase-2 prioritized connector selection...')
            phase2_prior = importlib.import_module('phase2_prioritized_connector')
            out_csv2, phase2_model, n2 = phase2_prior.run_phase2(cand, unconnected_model_json=model, out_dir=out_dir)
            print(f'Phase-2 prioritized: selected {n2} connectors -> model {phase2_model}')

        # Phase3 mode resolution (phase3-mode preferred, strategy kept for backward compatibility)
        phase3_mode = args.phase3_mode or args.strategy or 'exact'

        # Phase3: objective-specific
        if args.pipeline == 'deadends':
            print('Running Phase-3 (deadends minimizer)...')
            ph3_dead = importlib.import_module('phase3_greedy_optimizer')
            out_csv3, phase3_model, n3 = ph3_dead.run_phase3(cand, phase2_model_json=phase2_model, out_dir=out_dir, max_additions=args.max)
            print(f'Phase-3 deadends: added {n3} connectors -> model {phase3_model}')
        else:
            # blocked pipeline - choose strategy
            if phase3_mode == 'tiered_milp':
                print('Running Phase-3 (Tiered MILP: Type A → B → C per component)...')
                ph3_tiered = importlib.import_module('phase3_tiered_milp')
                components_csv = os.path.join(out_dir or os.path.join(os.path.dirname(__file__), 'files'), 'components_summary.csv')
                out_csv3, phase3_model, n3 = ph3_tiered.run_phase3_tiered(
                    cand, starting_model_json=phase2_model,
                    components_csv=components_csv, out_dir=out_dir,
                    tradeoff_lambda=args.tradeoff_lambda,
                    min_component_size=args.min_comp_size,
                    max_components=args.max_components, verbose=True)
                print(f'Phase-3 tiered MILP: added {n3} connectors -> model {phase3_model}')
            elif phase3_mode == 'component_milp':
                print('Running Phase-3 (Component-wise MILP optimizer - FASTEST)...')
                ph3_cmilp = importlib.import_module('phase3_component_milp')
                components_csv = os.path.join(out_dir or os.path.join(os.path.dirname(__file__), 'files'), 'components_summary.csv')
                out_csv3, phase3_model, n3 = ph3_cmilp.run_phase3_component_wise(
                    cand, starting_model_json=phase2_model, 
                    components_csv=components_csv, out_dir=out_dir,
                    tradeoff_lambda=args.tradeoff_lambda,
                    min_component_size=args.min_comp_size,
                    max_components=args.max_components, verbose=True)
                print(f'Phase-3 component MILP: added {n3} connectors -> model {phase3_model}')
            elif phase3_mode == 'milp':
                print('Running Phase-3 (Global MILP blocked-reaction optimizer)...')
                ph3_milp = importlib.import_module('phase3_milp_optimizer')
                out_csv3, phase3_model, n3 = ph3_milp.run_phase3(
                    cand, starting_model_json=phase2_model, out_dir=out_dir,
                    max_ptrs=args.max, tradeoff_lambda=args.tradeoff_lambda,
                    max_iterations=args.max_iter, verbose=True)
                print(f'Phase-3 MILP: added {n3} connectors -> model {phase3_model}')
            elif phase3_mode == 'sink_milp':
                print('Running Phase-3 (Sink MILP: temp sinks + original components + hybrid parallelization)...')
                ph3_sink = importlib.import_module('phase3_sink_milp_original')
                # Parse component IDs if provided
                component_ids = None
                if args.components:
                    component_ids = [int(c.strip()) for c in args.components.split(',')]
                # Run the sink MILP strategy
                result = ph3_sink.run_phase3_all_components(
                    cand, 
                    phase2_model_json=phase2_model,
                    unconnected_model_json=model,
                    out_dir=out_dir or os.path.join(os.path.dirname(__file__), 'files'),
                    tradeoff_lambda=args.tradeoff_lambda,
                    solver='dynamic',
                    solver_lp=args.solver_lp,
                    parallel_fba=args.parallel_fba,
                    n_workers_fba=args.workers_fba,
                    parallel_components=args.parallel_components,
                    n_workers_components=args.workers_components,
                    small_component_threshold=args.small_threshold,
                    min_component_size=args.min_comp_size,
                    component_ids=component_ids,
                    verbose=True
                )
                n3 = result.get('summary', {}).get('total_ptrs', len(result.get('all_selected', [])))
                out_csv3 = os.path.join(out_dir or os.path.join(os.path.dirname(__file__), 'files'), 'phase3_selected_ptrs.csv')
                phase3_model = None  # sink_milp doesn't create a new model yet
                print(f'Phase-3 sink MILP: selected {n3} PTRs -> {out_csv3}')
            elif phase3_mode == 'sink_milp_original':
                print('Running Phase-3 (Sink MILP original implementation with model export)...')
                ph3_sink_orig = importlib.import_module('phase3_sink_milp_original')
                component_ids = None
                if args.components:
                    component_ids = [int(c.strip()) for c in args.components.split(',')]
                result = ph3_sink_orig.run_phase3_all_components(
                    cand,
                    phase2_model_json=phase2_model,
                    unconnected_model_json=model,
                    out_dir=out_dir or os.path.join(os.path.dirname(__file__), 'files'),
                    tradeoff_lambda=args.tradeoff_lambda,
                    solver='dynamic',
                    solver_lp=args.solver_lp,
                    parallel_fba=args.parallel_fba,
                    n_workers_fba=args.workers_fba,
                    parallel_components=args.parallel_components,
                    n_workers_components=args.workers_components,
                    small_component_threshold=args.small_threshold,
                    min_component_size=args.min_comp_size,
                    component_ids=component_ids,
                    verbose=True,
                    phase3_model_out=args.phase3_model_out
                )
                n3 = result.get('summary', {}).get('total_ptrs', len(result.get('all_selected', [])))
                out_csv3 = os.path.join(out_dir or os.path.join(os.path.dirname(__file__), 'files'), 'phase3_selected_ptrs.csv')
                phase3_model = result.get('phase3_model', args.phase3_model_out)
                print(f'Phase-3 sink MILP original: selected {n3} PTRs -> {out_csv3}, model {phase3_model or "<not written>"}')
            else:
                print('Running Phase-3 (blocked-reaction minimizer)...')
                ph3_block = importlib.import_module('phase3_blocked_optimizer')
                out_csv3, phase3_model, n3 = ph3_block.run_phase3(cand, starting_model_json=phase2_model, out_dir=out_dir, max_additions=args.max, strategy=phase3_mode, topk=args.topk, batch_size=args.batch_size)
                print(f'Phase-3 blocked: added {n3} connectors -> model {phase3_model}')

        print('Run-all pipeline complete.')
        return

    parser.print_help()


if __name__ == '__main__':
    main()
