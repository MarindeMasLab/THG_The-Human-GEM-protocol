#!/usr/bin/env python3
"""List KEGG reaction IDs per pathway.

Writes a TSV `files/pathway_reactions.tsv` with columns:
pathway_id\tpathway_name\treaction_ids(space-separated)

By default runs only the first N pathways (dry-run). Use --all to run all.
"""
import os
import sys
import argparse

PROJECT_ROOT = os.path.dirname(os.path.dirname(__file__))
# Ensure project root is on sys.path so `functions` package can be imported
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from functions.class_generate_database import pathway
PATH_FILE = os.path.join(PROJECT_ROOT, 'files', 'human_kegg_pathways.txt')
OUT_FILE = os.path.join(PROJECT_ROOT, 'files', 'pathway_reactions.tsv')

parser = argparse.ArgumentParser()
parser.add_argument('--all', action='store_true', help='Process all pathways')
parser.add_argument('--limit', type=int, default=5, help='Dry-run limit (ignored with --all)')
parser.add_argument('--timeout', type=int, default=20, help='Timeout for KEGG requests (seconds)')
args = parser.parse_args()

if not os.path.exists(PATH_FILE):
    print('Pathways file not found:', PATH_FILE)
    sys.exit(1)

lines = [l for l in open(PATH_FILE, 'r').read().splitlines() if l.strip()]

# Prepare output
os.makedirs(os.path.join(PROJECT_ROOT, 'files'), exist_ok=True)
with open(OUT_FILE + '.tmp', 'w') as out:
    out.write('#pathway_id\tpathway_name\treaction_ids\n')
    count = 0
    for i, line in enumerate(lines):
        parts = line.split('\t')
        if len(parts) < 2:
            continue
        pid, pname = parts[0].strip(), parts[1].strip()
        # Respect dry-run limit
        if not args.all and count >= args.limit:
            break
        url = f'https://rest.kegg.jp/get/{pid}/kgml'
        referer = f'https://www.kegg.jp/kegg-bin/show_pathway?{pid}'
        try:
            p = pathway(url, args.timeout, pid, referer, pname)
            rxns = p.Reactions() or []
            # If no reactions found in this pathway, follow any linked map entries
            # (e.g. <entry ... name="path:hsa00512" ...>) one level deep and
            # merge their reactions. This helps when KGML uses pathway-maps
            # to delegate content instead of listing reactions directly.
            if not rxns:
                import re

                maps = re.findall(r'name="path:(hsa[0-9]{5})"', p.pagina or '')
                extra_rxns = []
                seen_map = set()
                for m in maps:
                    if m == pid or m in seen_map:
                        continue
                    seen_map.add(m)
                    try:
                        sub_url = f'https://rest.kegg.jp/get/{m}/kgml'
                        subp = pathway(sub_url, args.timeout, m, referer, m)
                        srx = subp.Reactions() or []
                        extra_rxns.extend(srx)
                    except Exception:
                        # ignore failures fetching linked maps
                        pass
                if extra_rxns:
                    # merge into rxns (preserve original shape: list of tuples)
                    rxns = extra_rxns
            # rxns entries are tuples: ((rn_id, type), url)
            rxn_ids = [r[0][0] for r in rxns if r and r[0] and r[0][0]]
            out.write(f"{pid}\t{pname}\t{' '.join(rxn_ids)}\n")
            print(f'[{i+1}/{len(lines)}] {pid} {pname}: {len(rxn_ids)} reactions')
        except Exception as e:
            print(f'Failed to process {pid} {pname}:', e)
            out.write(f"{pid}\t{pname}\tERROR\n")
        count += 1

# atomic replace
os.replace(OUT_FILE + '.tmp', OUT_FILE)
print('Wrote', OUT_FILE)
