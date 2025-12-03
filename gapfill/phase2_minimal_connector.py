"""Phase-2 minimal connector: connect all components with minimum connectors.

Selects a minimal set of candidate transport reactions to connect all components
present in `candidates_all.csv`. Preference order when choosing a representative
for a component-pair: Type A > Type B > Type C.

Outputs:
 - `gapfill/files/selected_minimal_connectors.csv`
 - model saved to `models/base/THG-beta-batch_251106_phase2_minimal_connected.json`
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
            r['component1'] = int(r['component1'])
            r['component2'] = int(r['component2'])
            rows.append(r)
    return rows


def representative(cands):
    # prefer A then B then C, deterministic tiebreaker
    order = {'A': 0, 'B': 1, 'C': 2}
    return sorted(cands, key=lambda c: (order.get(c.get('type','C'), 2), c.get('base',''), c.get('met1',''), c.get('met2','')))[0]


def run_phase2(candidates_csv,
               unconnected_model_json=None,
               out_dir=None,
               phase2_model_out=None):
    if out_dir is None:
        out_dir = os.path.join(os.path.dirname(__file__), 'files')
    os.makedirs(out_dir, exist_ok=True)

    cands = load_candidates(candidates_csv)

    comps = set()
    pairs = defaultdict(list)
    for c in cands:
        a = c['component1']; b = c['component2']
        comps.add(a); comps.add(b)
        key = (min(a,b), max(a,b))
        pairs[key].append(c)

    comps = sorted(comps)
    uf = UnionFind(comps)

    # Build list of edges (component pairs) with chosen representative
    edges = []
    for pair, lst in pairs.items():
        rep = representative(lst)
        edges.append((pair, rep))

    # Sort edges by candidate priority (A first) to encourage biologically plausible connectors
    def edge_score(item):
        pair, rep = item
        t = rep.get('type','C')
        order = {'A':0,'B':1,'C':2}
        return (order.get(t,2), rep.get('base',''))

    edges = sorted(edges, key=edge_score)

    selected = []
    # Kruskal-like: pick edges until all components connected
    for (a,b), rep in edges:
        if uf.find(a) != uf.find(b):
            selected.append(rep)
            uf.union(a,b)

    # add to model
    if unconnected_model_json is None:
        unconnected_model_json = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'models', 'base', 'THG-beta-batch_251106_unconnected.json'))
    model = load_json_model(unconnected_model_json)

    created = []
    idx = 0
    for sel in selected:
        met1 = sel['met1']; met2 = sel['met2']; base = sel.get('base','')
        s1 = sel.get('suffix1',''); s2 = sel.get('suffix2','')
        raw_id = f"TRANS_MIN_{base}_{s1}_{s2}_{idx}"
        rid = ''.join(ch if (ch.isalnum() or ch == '_') else '_' for ch in raw_id)
        idx += 1
        try:
            m1 = model.metabolites.get_by_id(met1)
            m2 = model.metabolites.get_by_id(met2)
        except KeyError:
            continue
        rxn = Reaction(rid)
        rxn.name = f"transport {met1} -> {met2}"
        rxn.add_metabolites({m1: -1.0, m2: 1.0})
        rxn.lower_bound = -1000.0
        rxn.upper_bound = 1000.0
        rxn.annotation = {'phase2_minimal': 'true', 'base': base, 'suffix1': s1, 'suffix2': s2}
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

    out_csv = os.path.join(out_dir, 'selected_minimal_connectors.csv')
    with open(out_csv, 'w', newline='') as f:
        fieldnames = ['reaction_id','met1','met2','base','suffix1','suffix2','component1','component2','type']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in created:
            writer.writerow(r)

    if phase2_model_out is None:
        phase2_model_out = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'models', 'base', 'THG-beta-batch_251106_phase2_minimal_connected.json'))
    save_json_model(model, phase2_model_out)

    return out_csv, phase2_model_out, len(created)


if __name__ == '__main__':
    cand_csv = os.path.join(os.path.dirname(__file__), 'files', 'candidates_all.csv')
    if not os.path.exists(cand_csv):
        print('Candidates CSV not found:', cand_csv)
    else:
        out_csv, model_out, n = run_phase2(cand_csv)
        print(f'Phase2 minimal selected {n} connectors, wrote {out_csv} and model {model_out}')
