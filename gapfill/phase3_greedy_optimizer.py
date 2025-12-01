"""Phase-3 greedy optimizer

Greedily selects Phase-1 candidate transport reactions to maximally reduce
dead-end metabolites. Uses the Phase-2 model as the starting point and adds
reactions until no candidate provides positive marginal benefit.

Outputs:
 - `gapfill/files/selected_phase3_connectors.csv`
 - `models/base/THG-beta-batch_251106_phase3_connected.json`

Run: `python3 gapfill/phase3_greedy_optimizer.py`
"""
import csv
import os
from collections import defaultdict
from cobra.io import load_json_model, save_json_model
from cobra import Reaction


def load_candidates(path):
    rows = []
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        for r in reader:
            # component fields may be int-like
            if 'component1' in r:
                try:
                    r['component1'] = int(r['component1'])
                    r['component2'] = int(r['component2'])
                except Exception:
                    pass
            rows.append(r)
    return rows


def compute_deadend_sets(model):
    produced = set()
    consumed = set()
    for rxn in model.reactions:
        for met, sto in rxn.metabolites.items():
            if sto > 0:
                produced.add(met.id)
            elif sto < 0:
                consumed.add(met.id)
    produced_only = set([m.id for m in model.metabolites if m.id in produced and m.id not in consumed])
    consumed_only = set([m.id for m in model.metabolites if m.id in consumed and m.id not in produced])
    return produced_only, consumed_only


def benefit_of_candidate(candidate, prod_only, cons_only):
    m1 = candidate['met1']; m2 = candidate['met2']
    b = 0
    # if linking produced-only to consumed-only in either direction resolves both
    if (m1 in prod_only and m2 in cons_only) or (m2 in prod_only and m1 in cons_only):
        b = 2
    else:
        if m1 in prod_only or m1 in cons_only:
            b += 1
        if m2 in prod_only or m2 in cons_only:
            b += 1 if (m2 not in (m1,)) else 0
    return b


def add_transport_to_model(model, sel, idx):
    met1 = sel['met1']; met2 = sel['met2']; base = sel.get('base','')
    s1 = sel.get('suffix1',''); s2 = sel.get('suffix2','')
    raw_id = f"PH3_TRANS_{base}_{s1}_{s2}_{idx}"
    rid = ''.join(ch if (ch.isalnum() or ch == '_') else '_' for ch in raw_id)
    try:
        m1 = model.metabolites.get_by_id(met1)
        m2 = model.metabolites.get_by_id(met2)
    except KeyError:
        return None
    rxn = Reaction(rid)
    rxn.name = f"phase3 transport {met1} -> {met2}"
    rxn.add_metabolites({m1: -1.0, m2: 1.0})
    rxn.lower_bound = -1000.0
    rxn.upper_bound = 1000.0
    rxn.annotation = {'phase3_selected': 'true'}
    model.add_reactions([rxn])
    return rid


def run_phase3(candidates_csv,
               phase2_model_json=None,
               out_dir=None,
               phase3_model_out=None,
               max_additions=2000):
    if out_dir is None:
        out_dir = os.path.join(os.path.dirname(__file__), 'files')
    os.makedirs(out_dir, exist_ok=True)

    candidates = load_candidates(candidates_csv)

    if phase2_model_json is None:
        phase2_model_json = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'models', 'base', 'THG-beta-batch_251106_phase2_connected.json'))
    model = load_json_model(phase2_model_json)

    remaining = [c for c in candidates]
    selected = []

    idx = 0
    # Greedy loop
    while idx < max_additions and remaining:
        prod_only, cons_only = compute_deadend_sets(model)
        # compute benefit per candidate
        best = None
        best_b = 0
        for c in remaining:
            b = benefit_of_candidate(c, prod_only, cons_only)
            if b > best_b:
                best_b = b
                best = c
        if best is None or best_b <= 0:
            break
        # add best
        rid = add_transport_to_model(model, best, idx)
        if rid is None:
            # remove candidate if metabolite missing
            remaining = [c for c in remaining if c is not best]
            continue
        selrec = dict(reaction_id=rid,
                      met1=best['met1'], met2=best['met2'], base=best.get('base',''),
                      suffix1=best.get('suffix1',''), suffix2=best.get('suffix2',''),
                      component1=best.get('component1',''), component2=best.get('component2',''),
                      type=best.get('type',''), benefit=best_b)
        selected.append(selrec)
        idx += 1
        # remove candidates that now share the same exact met pair
        remaining = [c for c in remaining if not (c['met1']==best['met1'] and c['met2']==best['met2'])]

    # write selected CSV
    out_csv = os.path.join(out_dir, 'selected_phase3_connectors.csv')
    with open(out_csv, 'w', newline='') as f:
        fieldnames = ['reaction_id','met1','met2','base','suffix1','suffix2','component1','component2','type','benefit']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in selected:
            writer.writerow(r)

    if phase3_model_out is None:
        phase3_model_out = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'models', 'base', 'THG-beta-batch_251106_phase3_connected.json'))
    save_json_model(model, phase3_model_out)

    return out_csv, phase3_model_out, len(selected)


if __name__ == '__main__':
    cand_csv = os.path.join(os.path.dirname(__file__), 'files', 'candidates_all.csv')
    if not os.path.exists(cand_csv):
        print('Candidates CSV not found:', cand_csv)
    else:
        out_csv, model_out, n = run_phase3(cand_csv)
        print(f'Phase3 selected {n} connectors, wrote {out_csv} and model {model_out}')
