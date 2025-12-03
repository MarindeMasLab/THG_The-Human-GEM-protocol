"""Phase-3 blocked-reaction optimizer

Greedily selects candidate connectors that minimize the number of blocked
reactions. At each step, evaluates remaining candidates by tentatively adding
them and computing `find_blocked_reactions`; picks the candidate with largest
decrease in blocked reactions. Stops when no candidate gives improvement or
when `max_additions` is reached.

Outputs:
 - `gapfill/files/selected_phase3_blocked_connectors.csv`
 - `models/base/THG-beta-batch_251106_phase3_blocked_connected.json`
"""
import csv
import os
import json
import time
from copy import deepcopy
from cobra.io import load_json_model, save_json_model
from cobra import Reaction
from cobra.flux_analysis import find_blocked_reactions


def load_candidates(path):
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


def add_transport(model, sel, idx, prefix='PH3B'):
    met1 = sel['met1']; met2 = sel['met2']; base = sel.get('base','')
    s1 = sel.get('suffix1',''); s2 = sel.get('suffix2','')
    raw_id = f"{prefix}_TRANS_{base}_{s1}_{s2}_{idx}"
    rid = ''.join(ch if (ch.isalnum() or ch == '_') else '_' for ch in raw_id)
    try:
        m1 = model.metabolites.get_by_id(met1)
        m2 = model.metabolites.get_by_id(met2)
    except KeyError:
        return None
    rxn = Reaction(rid)
    rxn.name = f"phase3 blocked transport {met1} -> {met2}"
    rxn.add_metabolites({m1: -1.0, m2: 1.0})
    rxn.lower_bound = -1000.0
    rxn.upper_bound = 1000.0
    rxn.annotation = {'phase3_blocked_selected': 'true'}
    model.add_reactions([rxn])
    return rid


def run_phase3(candidates_csv,
               starting_model_json=None,
               out_dir=None,
               phase3_model_out=None,
               max_additions=500,
               strategy='exact',
               topk=None,
               batch_size=None,
               resume=True):
    if out_dir is None:
        out_dir = os.path.join(os.path.dirname(__file__), 'files')
    os.makedirs(out_dir, exist_ok=True)

    candidates = load_candidates(candidates_csv)

    # checkpoint file paths
    checkpoint_model = os.path.join(out_dir, 'phase3_blocked_checkpoint_model.json')
    checkpoint_remaining = os.path.join(out_dir, 'phase3_blocked_remaining_candidates.csv')
    checkpoint_selected = os.path.join(out_dir, 'phase3_blocked_selected_so_far.csv')
    checkpoint_meta = os.path.join(out_dir, 'phase3_blocked_checkpoint_meta.json')

    # helper: compute simple dead-end benefit score from a model
    def compute_deadend_sets(model):
        produced = set()
        consumed = set()
        for rxn in model.reactions:
            for met, sto in rxn.metabolites.items():
                if sto > 0:
                    produced.add(met.id)
                if sto < 0:
                    consumed.add(met.id)
        produced_only = set([m.id for m in model.metabolites if m.id in produced and m.id not in consumed])
        consumed_only = set([m.id for m in model.metabolites if m.id in consumed and m.id not in produced])
        return produced_only, consumed_only

    def candidate_deadend_benefit(c, prod_only, cons_only):
        m1 = c['met1']; m2 = c['met2']
        b = 0
        if (m1 in prod_only and m2 in cons_only) or (m2 in prod_only and m1 in cons_only):
            return 2
        if m1 in prod_only or m1 in cons_only:
            b += 1
        if m2 in prod_only or m2 in cons_only:
            b += 1
        return b

    if starting_model_json is None:
        starting_model_json = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'models', 'base', 'THG-beta-batch_251106_phase2_minimal_connected.json'))

    model = None
    # resume from checkpoint if requested and available
    if resume and os.path.exists(checkpoint_model) and os.path.exists(checkpoint_remaining) and os.path.exists(checkpoint_selected):
        try:
            model = load_json_model(checkpoint_model)
            # load remaining candidates from checkpoint
            remaining = load_candidates(checkpoint_remaining)
            # load selected so far
            sel = []
            with open(checkpoint_selected, newline='') as f:
                r = csv.DictReader(f)
                for row in r:
                    try:
                        row['benefit_blocked'] = int(row.get('benefit_blocked', 0))
                    except Exception:
                        pass
                    sel.append(row)
            selected = sel
            idx = len(selected)
            # recompute base blocked on checkpoint model
            try:
                base_blocked = len(find_blocked_reactions(model))
            except Exception:
                base_blocked = 0
            # try to read meta
            meta = {}
            if os.path.exists(checkpoint_meta):
                try:
                    with open(checkpoint_meta) as mf:
                        meta = json.load(mf)
                except Exception:
                    meta = {}
            batch_idx = meta.get('batch_idx', 0)
            print(f"Resuming Phase-3 from checkpoint: {checkpoint_model}, selected={idx}, remaining={len(remaining)}")
        except Exception as e:
            print('Failed to load checkpoint, starting fresh:', e)
            model = load_json_model(starting_model_json)
            remaining = [c for c in candidates]
            selected = []
            idx = 0
            base_blocked = len(find_blocked_reactions(model))
            batch_idx = 0
    else:
        model = load_json_model(starting_model_json)

    # initial blocked reactions count
    if model is None:
        model = load_json_model(starting_model_json)
    if 'selected' not in locals():
        selected = []
    if 'idx' not in locals():
        idx = 0
    if 'batch_idx' not in locals():
        batch_idx = 0
    if 'remaining' not in locals():
        remaining = [c for c in candidates]
    if 'base_blocked' not in locals():
        try:
            base_blocked = len(find_blocked_reactions(model))
        except Exception:
            base_blocked = 0

    if strategy in ('hybrid', 'hybrid_batch') and topk:
        # Score candidates and keep topk
        def type_score(c):
            return {'A': 0, 'B': 1, 'C': 2}.get(c.get('type', 'C'), 2)

        scored = []
        for c in remaining:
            b = candidate_deadend_benefit(c, *compute_deadend_sets(model))
            # deterministic tiebreakers: base, met1, met2
            scored.append((type_score(c), -b, c.get('base',''), c.get('met1',''), c.get('met2',''), c))
        scored.sort()
        remaining = [t[5] for t in scored[:min(topk, len(scored))]]

    # helper to persist a checkpoint
    def save_checkpoint():
        try:
            save_json_model(model, checkpoint_model)
        except Exception as e:
            print('Warning: failed to save checkpoint model:', e)
        try:
            with open(checkpoint_remaining, 'w', newline='') as f:
                if remaining:
                    fieldnames = list(remaining[0].keys())
                else:
                    fieldnames = ['met1','met2','base','type','component1','component2']
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                for r in remaining:
                    writer.writerow(r)
        except Exception as e:
            print('Warning: failed to save checkpoint remaining:', e)
        try:
            with open(checkpoint_selected, 'w', newline='') as f:
                fieldnames = ['reaction_id','met1','met2','type','benefit_blocked']
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                for r in selected:
                    writer.writerow(r)
        except Exception as e:
            print('Warning: failed to save checkpoint selected:', e)
        try:
            meta = {'timestamp': time.time(), 'base_blocked': base_blocked, 'idx': idx, 'batch_idx': batch_idx, 'strategy': strategy, 'topk': topk, 'batch_size': batch_size}
            with open(checkpoint_meta, 'w') as mf:
                json.dump(meta, mf)
        except Exception as e:
            print('Warning: failed to save checkpoint meta:', e)

    # If using hybrid_batch, chunk into batches
    if strategy == 'hybrid_batch' and batch_size and batch_size > 0:
        # process batches in order; allow iterative improvements inside batches
        batches = [remaining[i:i+batch_size] for i in range(0, len(remaining), batch_size)]
        # resume support: start from saved batch_idx
        while idx < max_additions and batch_idx < len(batches):
            batch = batches[batch_idx]
            # test batch as a whole
            mcopy = model.copy()
            for c in batch:
                add_transport(mcopy, c, idx, prefix='PH3B_TMPB')
            try:
                br = find_blocked_reactions(mcopy)
                blocked_count = len(br)
            except Exception:
                blocked_count = base_blocked
            delta_batch = base_blocked - blocked_count
            if delta_batch <= 0:
                # skip this batch
                batch_idx += 1
                # save checkpoint after skipping a batch so progress isn't lost
                try:
                    save_checkpoint()
                except Exception:
                    pass
                continue
            # batch shows promise: refine inside batch
            inner_remaining = [c for c in batch]
            while idx < max_additions and inner_remaining:
                best = None
                best_delta = 0
                for c in inner_remaining:
                    mcopy = model.copy()
                    rid = add_transport(mcopy, c, idx, prefix='PH3B_TMP')
                    if rid is None:
                        continue
                    try:
                        br = find_blocked_reactions(mcopy)
                        bcount = len(br)
                    except Exception:
                        continue
                    delta = base_blocked - bcount
                    if delta > best_delta:
                        best_delta = delta
                        best = c
                if best is None or best_delta <= 0:
                    break
                # apply best to real model
                rid = add_transport(model, best, idx, prefix='PH3B')
                if rid is None:
                    inner_remaining = [c for c in inner_remaining if c is not best]
                    continue
                selected.append({'reaction_id': rid, 'met1': best['met1'], 'met2': best['met2'], 'type': best.get('type',''), 'benefit_blocked': best_delta})
                idx += 1
                try:
                    base_blocked = len(find_blocked_reactions(model))
                except Exception:
                    pass
                # remove best from both inner_remaining and batch
                inner_remaining = [c for c in inner_remaining if not (c['met1']==best['met1'] and c['met2']==best['met2'])]
                batch = [c for c in batch if not (c['met1']==best['met1'] and c['met2']==best['met2'])]
                # save checkpoint after each successful addition inside batch
                try:
                    save_checkpoint()
                except Exception:
                    pass
            # move to next batch
            batch_idx += 1
            # save checkpoint after finishing a batch
            try:
                save_checkpoint()
            except Exception:
                pass
    else:
        # fallback: exact greedy evaluation over remaining candidates (original behavior)
        while idx < max_additions and remaining:
            best = None
            best_delta = 0
            for i, c in enumerate(remaining):
                mcopy = model.copy()
                rid = add_transport(mcopy, c, idx, prefix='PH3B_TMP')
                if rid is None:
                    continue
                try:
                    br = find_blocked_reactions(mcopy)
                    blocked_count = len(br)
                except Exception:
                    continue
                delta = base_blocked - blocked_count
                if delta > best_delta:
                    best_delta = delta
                    best = c
            if best is None or best_delta <= 0:
                break
            rid = add_transport(model, best, idx, prefix='PH3B')
            if rid is None:
                remaining = [c for c in remaining if c is not best]
                continue
            selected.append({'reaction_id': rid, 'met1': best['met1'], 'met2': best['met2'], 'type': best.get('type',''), 'benefit_blocked': best_delta})
            idx += 1
            try:
                base_blocked = len(find_blocked_reactions(model))
            except Exception:
                pass
            remaining = [c for c in remaining if not (c['met1']==best['met1'] and c['met2']==best['met2'])]
            # save checkpoint after each addition in exact/hybrid strategy
            try:
                save_checkpoint()
            except Exception:
                pass

    out_csv = os.path.join(out_dir, 'selected_phase3_blocked_connectors.csv')
    with open(out_csv, 'w', newline='') as f:
        fieldnames = ['reaction_id','met1','met2','type','benefit_blocked']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in selected:
            writer.writerow(r)

    if phase3_model_out is None:
        phase3_model_out = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'models', 'base', 'THG-beta-batch_251106_phase3_blocked_connected.json'))
    save_json_model(model, phase3_model_out)

    # remove checkpoint files on successful completion
    try:
        if os.path.exists(checkpoint_model):
            os.remove(checkpoint_model)
        if os.path.exists(checkpoint_remaining):
            os.remove(checkpoint_remaining)
        if os.path.exists(checkpoint_selected):
            os.remove(checkpoint_selected)
        if os.path.exists(checkpoint_meta):
            os.remove(checkpoint_meta)
    except Exception:
        pass

    return out_csv, phase3_model_out, len(selected)


if __name__ == '__main__':
    cand_csv = os.path.join(os.path.dirname(__file__), 'files', 'candidates_all.csv')
    if not os.path.exists(cand_csv):
        print('Candidates CSV not found:', cand_csv)
    else:
        out_csv, model_out, n = run_phase3(cand_csv)
        print(f'Phase3 blocked-selected {n} connectors, wrote {out_csv} and model {model_out}')
