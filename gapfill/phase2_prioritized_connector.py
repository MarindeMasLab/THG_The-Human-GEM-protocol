"""Phase-2 prioritized connector selection

Selects a small set of transport reactions from Phase-1 candidates using the
priority: A (substrate<->product) -> B (same-type deadends) -> fallback MST
on component graph. Writes `selected_connectors.csv` to `gapfill/files/` and
adds the reactions to the unconnected model, saving a Phase-2 model JSON.

Run from repository root: `python3 gapfill/phase2_prioritized_connector.py`
"""
import csv
import os
from collections import defaultdict
from cobra.io import load_json_model, save_json_model
from cobra import Reaction


class UnionFind:
    def __init__(self, elements):
        self.parent = {e: e for e in elements}

    def find(self, x):
        p = self.parent.get(x, x)
        if p != x:
            p = self.find(p)
            self.parent[x] = p
        return p

    def union(self, a, b):
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return False
        self.parent[rb] = ra
        return True


def load_candidates(path):
    rows = []
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        for r in reader:
            # ensure components are ints
            r['component1'] = int(r['component1'])
            r['component2'] = int(r['component2'])
            rows.append(r)
    return rows


def group_by_component_pair(candidates):
    groups = defaultdict(list)
    for c in candidates:
        a = c['component1']; b = c['component2']
        key = (min(a,b), max(a,b))
        groups[key].append(c)
    return groups


def pick_representative(cands):
    # Simple deterministic pick: prefer Type A then B then C, else lexicographic
    def score(c):
        t = c.get('type','C')
        pri = {'A':0,'B':1,'C':2}.get(t,2)
        return (pri, c.get('base',''), c.get('met1',''), c.get('met2',''))
    return sorted(cands, key=score)[0]


def run_phase2(candidates_csv,
               unconnected_model_json=None,
               out_dir=None,
               phase2_model_out=None):
    if out_dir is None:
        out_dir = os.path.join(os.path.dirname(__file__), 'files')
    os.makedirs(out_dir, exist_ok=True)

    candidates = load_candidates(candidates_csv)

    # list components present
    comps = set()
    for c in candidates:
        comps.add(c['component1']); comps.add(c['component2'])
    comps = sorted(comps)

    uf = UnionFind(comps)

    groups = group_by_component_pair(candidates)

    selected = []

    # Phase 2a: Type A (substrate<->product)
    for pair, lits in sorted(groups.items()):
        a, b = pair
        # if already connected, skip
        if uf.find(a) == uf.find(b):
            continue
        a_cands = [c for c in lits if c.get('type') == 'A']
        if a_cands:
            chosen = pick_representative(a_cands)
            selected.append(chosen)
            uf.union(a, b)

    # Phase 2b: Type B (same-type deadends)
    for pair, lits in sorted(groups.items()):
        a, b = pair
        if uf.find(a) == uf.find(b):
            continue
        b_cands = [c for c in lits if c.get('type') == 'B']
        if b_cands:
            chosen = pick_representative(b_cands)
            selected.append(chosen)
            uf.union(a, b)

    # Phase 2c: fallback MST-like greedy on remaining pairs
    for pair, lits in sorted(groups.items()):
        a, b = pair
        if uf.find(a) == uf.find(b):
            continue
        chosen = pick_representative(lits)
        selected.append(chosen)
        uf.union(a, b)

    # Now add selected reactions to the model
    if unconnected_model_json is None:
        unconnected_model_json = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'models', 'base', 'THG-beta-batch_251106_unconnected.json'))
    model = load_json_model(unconnected_model_json)

    created = []
    idx = 0
    for sel in selected:
        met1 = sel['met1']; met2 = sel['met2']; base = sel.get('base','')
        s1 = sel.get('suffix1',''); s2 = sel.get('suffix2','')
        raw_id = f"TRANS_{base}_{s1}_{s2}_{idx}"
        rid = ''.join(ch if (ch.isalnum() or ch == '_') else '_' for ch in raw_id)
        idx += 1
        try:
            m1 = model.metabolites.get_by_id(met1)
            m2 = model.metabolites.get_by_id(met2)
        except KeyError:
            # skip if metabolite missing
            continue
        rxn = Reaction(rid)
        rxn.name = f"transport {met1} -> {met2}"
        rxn.add_metabolites({m1: -1.0, m2: 1.0})
        rxn.lower_bound = -1000.0
        rxn.upper_bound = 1000.0
        rxn.annotation = {'phase2_selected': 'true', 'base': base, 'suffix1': s1, 'suffix2': s2, 'component1': sel['component1'], 'component2': sel['component2']}
        model.add_reactions([rxn])
        created.append({
            'reaction_id': rid,
            'met1': met1,
            'met2': met2,
            'base': base,
            'suffix1': s1,
            'suffix2': s2,
            'component1': sel['component1'],
            'component2': sel['component2'],
            'type': sel.get('type','')
        })

    # write selected_connectors.csv
    out_csv = os.path.join(out_dir, 'selected_connectors.csv')
    with open(out_csv, 'w', newline='') as f:
        fieldnames = ['reaction_id','met1','met2','base','suffix1','suffix2','component1','component2','type']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in created:
            writer.writerow(r)

    # save phase2 model
    if phase2_model_out is None:
        phase2_model_out = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'models', 'base', 'THG-beta-batch_251106_phase2_connected.json'))
    save_json_model(model, phase2_model_out)

    return out_csv, phase2_model_out, len(created)


if __name__ == '__main__':
    cand_csv = os.path.join(os.path.dirname(__file__), 'files', 'candidates_all.csv')
    if not os.path.exists(cand_csv):
        print('Candidates CSV not found:', cand_csv)
    else:
        out_csv, model_out, n = run_phase2(cand_csv)
        print(f'Phase2 selected {n} connectors, wrote {out_csv} and model {model_out}')
