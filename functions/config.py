"""Compatibility wrapper for configuration helpers.

New code should import from :mod:`thg_protocol.config`. This module keeps the
legacy ``functions.config`` import path working during the package migration.
"""

import sys
from pathlib import Path

src_path = Path(__file__).resolve().parents[1] / "src"
if src_path.exists() and str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

from thg_protocol.config import (  # noqa: E402,F401
    get_compartments,
    get_model_paths,
    get_project_root,
    load_config,
    resolve_all_compartments,
    resolve_compartment_abbreviation,
)

__all__ = [
    "get_compartments",
    "get_model_paths",
    "get_project_root",
    "load_config",
    "resolve_all_compartments",
    "resolve_compartment_abbreviation",
]
