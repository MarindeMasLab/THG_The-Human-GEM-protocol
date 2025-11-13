"""
Functions package for pathway implementation pipeline

This package provides generic pathway builder functions for adding
any pathway/compartment to metabolic models.
"""

# Configuration management
from .config import (
    load_config,
    get_model_paths,
    resolve_all_compartments
)

# Generic pathway builder
from .pathway_builder import (
    add_compartment,
    create_compartment_metabolites,
    create_compartment_reactions,
    check_pathway_exists,
    find_metabolite_robust,
    find_metabolite_by_annotation,
    get_next_metabolite_id,
    get_next_reaction_id
)
