#!/usr/bin/env python3
"""
Print compartments resolved from BioCyc caches.

Reads:
 - files/biovelo_location_cache.pkl (mapping Biocyc ID -> list of location dicts)
 - files/compartment_name_cache.pkl (mapping 'ORG:FRAMEID' -> common name)

Outputs a sorted list of unique compartment names and a mapping sample.
"""
import os
import pickle
from collections import defaultdict

project_root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
biovelo_loc_file = os.path.join(project_root, "files", "caches", "biovelo_location_cache.pkl")
compartment_name_file = os.path.join(
    project_root, "files", "compartment_name_cache.pkl"
)


def load_pickle(path):
    if not os.path.exists(path):
        return None
    try:
        with open(path, "rb") as f:
            return pickle.load(f)
    except Exception as e:
        print(f"Failed to load {path}: {e}")
        return None


biovelo_loc = load_pickle(biovelo_loc_file) or {}
compartment_names = load_pickle(compartment_name_file) or {}

# Collect frameids found in biovelo_loc
frameids = set()
for bid, locations in biovelo_loc.items():
    if isinstance(locations, list):
        for loc in locations:
            if isinstance(loc, dict):
                org = loc.get("orgid", "HUMAN") or "HUMAN"
                fid = loc.get("frameid")
                if fid:
                    frameids.add((org, fid))

# Map frameids to resolved names
resolved = {}
for org, fid in sorted(frameids):
    key = f"{org}:{fid}"
    name = compartment_names.get(key, None)
    resolved[key] = name or fid

# Unique compartment common names
unique_names = set(resolved.values())

print("BioCyc compartments summary")
print("=================================")
print(f"Biovelo location entries: {len(biovelo_loc)}")
print(f"Unique frameids found: {len(frameids)}")
print(f"Unique resolved compartment names: {len(unique_names)}")
print()
print("Resolved compartment mapping (org:frameid -> name):")
for k in sorted(resolved.keys()):
    print(f"{k} -> {resolved[k]}")

print()
print("Unique compartment names (sorted):")
for n in sorted(unique_names):
    print(n)

# Also print a small sample of locations for a few Biocyc IDs
print()
print("Sample Biocyc ID -> locations")
count = 0
for bid, locations in sorted(biovelo_loc.items()):
    if not locations:
        continue
    print(f"{bid}:")
    for loc in locations:
        org = loc.get("orgid", "HUMAN")
        fid = loc.get("frameid")
        key = f"{org}:{fid}"
        print(f"  - frameid={fid} ({key}) -> {compartment_names.get(key, fid)}")
    count += 1
    if count >= 10:
        break

if not frameids:
    print(
        "No frameids found in biovelo_location_cache.pkl - maybe the cache is missing or empty."
    )
