#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extend an existing KEGG-based model with Rhea reactions, or create a new Rhea-only model.

This script can either:
1. Load an existing SBML model (generated from KEGG) and extend it with Rhea reactions
2. Create a brand new model containing only Rhea reactions (standalone mode)

Usage:
    # Extend existing model with all human Rhea reactions
    python extend_model_with_rhea.py <input_model.xml> [output_model.xml] [--rhea-ids ID1,ID2,...]
    
    # Create standalone Rhea-only model (no KEGG model required)
    python extend_model_with_rhea.py --standalone output_model.xml [--rhea-ids ID1,ID2,...]

Example:
    # Extend existing KEGG model with all human Rhea reactions
    python extend_model_with_rhea.py ../models/Human_database_20251201_fixed.xml
    
    # Create new model with specific Rhea IDs only
    python extend_model_with_rhea.py --standalone ../models/Rhea_model.xml --rhea-ids 10040,10112,10116
    
    # Create new model with ALL human Rhea reactions
    python extend_model_with_rhea.py --standalone ../models/Human_Rhea_model.xml

Author: Copilot / Igor
Date: December 2025
"""

import os
import sys
import re
import logging
import argparse
from datetime import datetime
from typing import Dict, List, Set, Any, Optional
from collections import defaultdict
from tqdm import tqdm

# Add project root to path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.join(current_dir, "..")
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import cobra
from cobra.io import read_sbml_model, write_sbml_model

# Import existing classes and functions
from functions.class_generate_database import reaction, compound, gene, gpr
from functions.gpr.auth_gpr import setup_biocyc_session
from rhea_extension import RheaExtender

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
LOGGER = logging.getLogger(__name__)


# =============================================================================
# CREATE EMPTY DATA STRUCTURES FOR STANDALONE MODE
# =============================================================================

def create_empty_data_structures() -> dict:
    """
    Create empty data structures for standalone Rhea model creation.
    
    This allows running the Rhea extension without any existing KEGG model.
    
    Returns:
        Dictionary containing all empty data structures needed by RheaExtender
    """
    LOGGER.info("Creating empty data structures for standalone Rhea model")
    
    return {
        'MetList': {},
        'MetIdent': [],
        'MetEquiv': {},
        'RxnList': {},
        'RxnIdent': [],
        'GPRList': {},
        'GPRIdent': [],
        'GeneList': {},
        'GeneIdent': [],
        'RxnList_CL': {},
        'MetList_CL': {},
        'RxnIdent_CL': [],
        'MetIdent_CL': [],
        'Compartment_CL': [],  # Will be populated by Rhea extension
        'PathNameRxn': {},
    }


def create_empty_cobra_model(model_id: str = "Rhea_Human_Model") -> cobra.Model:
    """
    Create an empty COBRA model with standard compartments.
    
    Args:
        model_id: ID for the new model
        
    Returns:
        Empty COBRA model with standard compartments configured
    """
    LOGGER.info(f"Creating new empty COBRA model: {model_id}")
    
    model = cobra.Model(model_id)
    
    # Set standard human compartments
    model.compartments = {
        'c': 'cytosol',
        'm': 'mitochondria',
        'n': 'nucleus',
        'r': 'endoplasmic reticulum',
        'g': 'golgi apparatus',
        'l': 'lysosome',
        'x': 'peroxisome',
        'e': 'extracellular',
    }
    
    return model


# =============================================================================
# EXTRACT DATA STRUCTURES FROM SBML MODEL
# =============================================================================

def extract_metabolites_from_model(model: cobra.Model) -> tuple:
    """
    Extract metabolite information from an existing COBRA model.
    
    Returns:
        tuple: (MetList, MetIdent, MetEquiv) where:
            - MetList: Dict[str, compound-like] mapping base IDs to compound info
            - MetIdent: List of all metabolite base IDs
            - MetEquiv: Dict for equivalent metabolite mappings
    """
    MetList = {}
    MetIdent = []
    MetEquiv = {}  # Will be populated if needed
    
    LOGGER.info(f"Extracting metabolites from model ({len(model.metabolites)} total)")
    
    for met in tqdm(model.metabolites, desc="Extracting metabolites"):
        # Extract base ID (without compartment suffix)
        base_id = met.id.rsplit('_', 1)[0] if '_' in met.id else met.id
        
        # Skip if we already have this base metabolite
        if base_id in MetIdent:
            continue
            
        MetIdent.append(base_id)
        
        # Create a simple compound-like object with the info we need
        MetList[base_id] = MetaboliteProxy(
            id=base_id,
            name=met.name or "",
            formula=met.formula or "",
            charge=met.charge if met.charge is not None else 0,
            annotation=met.annotation or {},
        )
    
    LOGGER.info(f"  Extracted {len(MetIdent)} unique metabolites")
    return MetList, MetIdent, MetEquiv


def extract_reactions_from_model(model: cobra.Model) -> tuple:
    """
    Extract reaction information from an existing COBRA model.
    
    Returns:
        tuple: (RxnList, RxnIdent) where:
            - RxnList: Dict[str, reaction-like] mapping base IDs to reaction info
            - RxnIdent: List of all reaction base IDs
    """
    RxnList = {}
    RxnIdent = []
    
    LOGGER.info(f"Extracting reactions from model ({len(model.reactions)} total)")
    
    for rxn in tqdm(model.reactions, desc="Extracting reactions"):
        # Extract base ID (without compartment suffix)
        base_id = rxn.id.rsplit('_', 1)[0] if '_' in rxn.id else rxn.id
        
        # Skip exchange/demand/sink reactions
        if rxn.id.startswith(('EX_', 'DM_', 'SK_')):
            continue
        
        # Skip if we already have this base reaction
        if base_id in RxnIdent:
            continue
            
        RxnIdent.append(base_id)
        
        # Extract EC numbers from annotation
        ec_numbers = []
        if rxn.annotation:
            ec_list = rxn.annotation.get('ec-code', [])
            if isinstance(ec_list, str):
                ec_list = [ec_list]
            ec_numbers = [ec for ec in ec_list if ec]
        
        # Extract KEGG reaction IDs from annotation
        kegg_ids = []
        if rxn.annotation:
            kegg_list = rxn.annotation.get('kegg.reaction', [])
            if isinstance(kegg_list, str):
                kegg_list = [kegg_list]
            kegg_ids = [k for k in kegg_list if k]
        
        # Create a simple reaction-like object
        RxnList[base_id] = ReactionProxy(
            id=base_id,
            name=rxn.name or "",
            ec_numbers=ec_numbers,
            kegg_ids=kegg_ids,
            annotation=rxn.annotation or {},
            gpr=str(rxn.gene_reaction_rule) if rxn.gene_reaction_rule else "",
        )
    
    LOGGER.info(f"  Extracted {len(RxnIdent)} unique reactions")
    return RxnList, RxnIdent


def extract_gpr_from_model(model: cobra.Model) -> tuple:
    """
    Extract GPR information from an existing COBRA model.
    
    Returns:
        tuple: (GPRList, GPRIdent, GeneList, GeneIdent)
    """
    GPRList = {}
    GPRIdent = []
    GeneList = {}
    GeneIdent = []
    
    LOGGER.info(f"Extracting GPR/gene info from model ({len(model.genes)} genes)")
    
    # Extract genes
    for gene_obj in model.genes:
        if gene_obj.id and gene_obj.id not in GeneIdent:
            GeneIdent.append(gene_obj.id)
            GeneList[gene_obj.id] = GeneProxy(
                id=gene_obj.id,
                name=gene_obj.name or "",
                annotation=gene_obj.annotation or {},
            )
    
    # Extract EC numbers and map to GPRs
    for rxn in model.reactions:
        if rxn.annotation:
            ec_list = rxn.annotation.get('ec-code', [])
            if isinstance(ec_list, str):
                ec_list = [ec_list]
            
            for ec in ec_list:
                if ec and ec not in GPRIdent:
                    GPRIdent.append(ec)
                    GPRList[ec] = GPRProxy(
                        ec=ec,
                        gpr_rule=str(rxn.gene_reaction_rule) if rxn.gene_reaction_rule else "",
                    )
    
    LOGGER.info(f"  Extracted {len(GeneIdent)} genes, {len(GPRIdent)} EC-GPR mappings")
    return GPRList, GPRIdent, GeneList, GeneIdent


def extract_compartmentalized_data(model: cobra.Model) -> tuple:
    """
    Extract compartmentalized reaction/metabolite information.
    
    Returns:
        tuple: (RxnList_CL, MetList_CL, RxnIdent_CL, MetIdent_CL, Compartment_CL)
    """
    RxnList_CL = {}
    MetList_CL = {}
    RxnIdent_CL = []
    MetIdent_CL = []
    Compartment_CL = list(model.compartments.keys()) if model.compartments else []
    
    LOGGER.info("Extracting compartmentalized data")
    
    # For now, we store reaction IDs that exist per compartment
    for rxn in model.reactions:
        if rxn.id not in RxnIdent_CL:
            RxnIdent_CL.append(rxn.id)
            # Store a simple reference
            RxnList_CL[rxn.id] = ReactionProxy(
                id=rxn.id,
                name=rxn.name or "",
                ec_numbers=[],
                kegg_ids=[],
                annotation=rxn.annotation or {},
                gpr=str(rxn.gene_reaction_rule) if rxn.gene_reaction_rule else "",
            )
    
    for met in model.metabolites:
        if met.id not in MetIdent_CL:
            MetIdent_CL.append(met.id)
            MetList_CL[met.id] = MetaboliteProxy(
                id=met.id,
                name=met.name or "",
                formula=met.formula or "",
                charge=met.charge if met.charge is not None else 0,
                annotation=met.annotation or {},
            )
    
    LOGGER.info(f"  Found {len(RxnIdent_CL)} compartmentalized reactions")
    LOGGER.info(f"  Found {len(MetIdent_CL)} compartmentalized metabolites")
    LOGGER.info(f"  Compartments: {Compartment_CL}")
    
    return RxnList_CL, MetList_CL, RxnIdent_CL, MetIdent_CL, Compartment_CL


def extract_pathways_from_model(model: cobra.Model) -> Dict[str, str]:
    """
    Extract pathway groupings from an existing model.
    
    Returns:
        Dict mapping reaction IDs to pathway names
    """
    PathNameRxn = {}
    
    LOGGER.info("Extracting pathway groups")
    
    for group in model.groups:
        pathway_name = group.name or group.id
        for member in group.members:
            if hasattr(member, 'id'):
                # Extract base ID
                base_id = member.id.rsplit('_', 1)[0] if '_' in member.id else member.id
                if base_id not in PathNameRxn:
                    PathNameRxn[base_id] = pathway_name
                else:
                    # Append to existing pathway (space-separated)
                    if pathway_name not in PathNameRxn[base_id]:
                        PathNameRxn[base_id] = f"{PathNameRxn[base_id]} {pathway_name}"
    
    LOGGER.info(f"  Found {len(PathNameRxn)} reaction-pathway mappings from {len(model.groups)} groups")
    return PathNameRxn


# =============================================================================
# PROXY CLASSES
# =============================================================================

class MetaboliteProxy:
    """Simple proxy class to represent a metabolite for Rhea extension compatibility."""
    
    def __init__(self, id: str, name: str, formula: str, charge: int, annotation: dict):
        self._id = id
        self._name = name
        self._formula = formula
        self._charge = charge
        self._annotation = annotation
    
    def ID(self) -> str:
        return self._id
    
    def Name(self) -> str:
        return self._name
    
    def Formula(self) -> str:
        return self._formula
    
    def Charge(self) -> int:
        return self._charge
    
    def Annotation(self) -> dict:
        return self._annotation


class ReactionProxy:
    """Simple proxy class to represent a reaction for Rhea extension compatibility."""
    
    def __init__(self, id: str, name: str, ec_numbers: list, kegg_ids: list, 
                 annotation: dict, gpr: str):
        self._id = id
        self._name = name
        self._ec_numbers = ec_numbers
        self._kegg_ids = kegg_ids
        self._annotation = annotation
        self._gpr = gpr
    
    def ID(self) -> str:
        return self._id
    
    def Name(self) -> str:
        return self._name
    
    def EC(self) -> list:
        return self._ec_numbers
    
    def KEGG_IDs(self) -> list:
        return self._kegg_ids
    
    def Annotation(self) -> dict:
        return self._annotation
    
    def GPR(self) -> str:
        return self._gpr


class GeneProxy:
    """Simple proxy class to represent a gene for Rhea extension compatibility."""
    
    def __init__(self, id: str, name: str, annotation: dict):
        self._id = id
        self._name = name
        self._annotation = annotation
    
    def ID(self) -> str:
        return self._id
    
    def Name(self) -> str:
        return self._name
    
    def Annotation(self) -> dict:
        return self._annotation


class GPRProxy:
    """Simple proxy class to represent a GPR mapping for Rhea extension compatibility."""
    
    def __init__(self, ec: str, gpr_rule: str):
        self._ec = ec
        self._gpr_rule = gpr_rule
    
    def EC(self) -> str:
        return self._ec
    
    def GprSubcell(self):
        """Return GPR data in expected format, or empty if no data."""
        if self._gpr_rule:
            # Return tuple format expected by pipeline: (sGPR, GPR, sgpr_dict, gpr_dict)
            return (self._gpr_rule, self._gpr_rule, {}, {})
        return None


# =============================================================================
# MAIN EXTENSION FUNCTION
# =============================================================================

def extend_model_with_rhea(
    input_model_path: Optional[str] = None,
    output_model_path: Optional[str] = None,
    rhea_ids: Optional[Set[str]] = None,
    impose_locations: int = 1,
    standalone: bool = False,
    model_id: str = "Rhea_Human_Model",
) -> cobra.Model:
    """
    Load an existing SBML model and extend it with Rhea reactions,
    or create a new Rhea-only model in standalone mode.
    
    Args:
        input_model_path: Path to existing SBML model (not needed if standalone=True)
        output_model_path: Path for output (default: adds _rhea_extended suffix or Rhea_model.xml)
        rhea_ids: Optional set of specific Rhea IDs to add (None = all human Rhea)
        impose_locations: Compartmentalization mode (1=restricted, 0=unrestricted)
        standalone: If True, create a new model without loading existing KEGG model
        model_id: ID for the new model (used in standalone mode)
        
    Returns:
        Extended or new COBRA model
    """
    LOGGER.info("=" * 80)
    if standalone:
        LOGGER.info("RHEA MODEL CREATION - Standalone mode (no KEGG model)")
    else:
        LOGGER.info("RHEA EXTENSION - Extending existing KEGG model")
    LOGGER.info("=" * 80)
    
    if standalone:
        # Standalone mode: create empty structures and new model
        LOGGER.info(f"\n[Step 1] Creating new empty model (standalone mode)")
        model = create_empty_cobra_model(model_id)
        data = create_empty_data_structures()
        
        MetList = data['MetList']
        MetIdent = data['MetIdent']
        MetEquiv = data['MetEquiv']
        RxnList = data['RxnList']
        RxnIdent = data['RxnIdent']
        GPRList = data['GPRList']
        GPRIdent = data['GPRIdent']
        GeneList = data['GeneList']
        GeneIdent = data['GeneIdent']
        RxnList_CL = data['RxnList_CL']
        MetList_CL = data['MetList_CL']
        RxnIdent_CL = data['RxnIdent_CL']
        MetIdent_CL = data['MetIdent_CL']
        Compartment_CL = data['Compartment_CL']
        PathNameRxn = data['PathNameRxn']
        
        LOGGER.info(f"  Model ID: {model.id}")
        LOGGER.info(f"  Compartments: {list(model.compartments.keys())}")
    else:
        # Extension mode: load existing model
        if not input_model_path:
            raise ValueError("input_model_path is required when standalone=False")
            
        LOGGER.info(f"\n[Step 1] Loading existing model: {input_model_path}")
        model = read_sbml_model(input_model_path)
        LOGGER.info(f"  Model: {model.id}")
        LOGGER.info(f"  Reactions: {len(model.reactions)}")
        LOGGER.info(f"  Metabolites: {len(model.metabolites)}")
        LOGGER.info(f"  Genes: {len(model.genes)}")
        LOGGER.info(f"  Groups: {len(model.groups)}")
        
        # Step 2: Extract data structures
        LOGGER.info("\n[Step 2] Extracting data structures from model...")
        
        MetList, MetIdent, MetEquiv = extract_metabolites_from_model(model)
        RxnList, RxnIdent = extract_reactions_from_model(model)
        GPRList, GPRIdent, GeneList, GeneIdent = extract_gpr_from_model(model)
        RxnList_CL, MetList_CL, RxnIdent_CL, MetIdent_CL, Compartment_CL = extract_compartmentalized_data(model)
        PathNameRxn = extract_pathways_from_model(model)
    
    # Step 3: Setup BioCyc session
    step_num = 2 if standalone else 3
    LOGGER.info(f"\n[Step {step_num}] Setting up BioCyc session...")
    session = setup_biocyc_session()
    
    # Step 4: Initialize Rhea extender
    step_num += 1
    LOGGER.info(f"\n[Step {step_num}] Initializing Rhea extender...")
    rhea_extender = RheaExtender(cache_dir=os.path.join(current_dir, '.rhea_cache'))
    
    # Step 5: Run Rhea extension
    step_num += 1
    LOGGER.info(f"\n[Step {step_num}] Running Rhea extension...")
    
    # Load special compounds file
    specialCompounds = os.path.join(project_root, "files", "special_compounds.txt")
    
    # Determine Rhea IDs to process
    limit_rhea_ids = rhea_ids if rhea_ids else None
    if limit_rhea_ids:
        LOGGER.info(f"  Limiting to {len(limit_rhea_ids)} specific Rhea IDs")
    else:
        LOGGER.info("  Processing ALL human Rhea reactions" + ("" if standalone else " not in KEGG"))
    
    rhea_result = rhea_extender.run(
        RxnList=RxnList,
        MetList=MetList,
        RxnIdent=RxnIdent,
        MetIdent=MetIdent,
        MetEquiv=MetEquiv,
        GPRList=GPRList,
        GPRIdent=GPRIdent,
        GeneList=GeneList,
        GeneIdent=GeneIdent,
        RxnList_CL=RxnList_CL,
        MetList_CL=MetList_CL,
        RxnIdent_CL=RxnIdent_CL,
        MetIdent_CL=MetIdent_CL,
        Compartment_CL=Compartment_CL,
        PathNameRxn=PathNameRxn,
        session=session,
        time=20,
        EF=[],
        specialCompounds=specialCompounds,
        limit_rhea_ids=limit_rhea_ids,
        impose_locations=impose_locations,
    )
    
    reactions_added = rhea_result.get('reactions_added', 0)
    LOGGER.info(f"\n  Rhea extension added {reactions_added} reactions")
    
    # Step 6: Add Rhea reactions directly to the existing model
    LOGGER.info("\n[Step 6] Adding Rhea reactions to existing model...")
    
    # Build location dict from model compartments
    location_dict = {}
    for comp_id, comp_name in model.compartments.items():
        location_dict[comp_name.lower()] = comp_id
    # Also add direct mappings
    for comp_id in model.compartments:
        location_dict[comp_id] = comp_id
    
    # Add new metabolites from Rhea (those in MetList_CL but not in model)
    new_mets_added = 0
    for met_id, met_obj in MetList_CL.items():
        if met_id not in [m.id for m in model.metabolites]:
            try:
                # Extract compartment from ID
                if '_' in met_id:
                    base_id, comp_name = met_id.rsplit('_', 1)
                    comp_id = location_dict.get(comp_name.lower(), comp_name)
                else:
                    base_id = met_id
                    comp_id = 'c'
                
                # Check if this is a Rhea/ChEBI metabolite
                if hasattr(met_obj, 'ID') and callable(met_obj.ID):
                    formula = met_obj.Formula() if hasattr(met_obj, 'Formula') and callable(met_obj.Formula) else ""
                    charge = met_obj.Charge() if hasattr(met_obj, 'Charge') and callable(met_obj.Charge) else 0
                    name = met_obj.Name() if hasattr(met_obj, 'Name') and callable(met_obj.Name) else base_id
                    annotation = met_obj.Annotation() if hasattr(met_obj, 'Annotation') and callable(met_obj.Annotation) else {}
                else:
                    # MetaboliteProxy
                    formula = met_obj._formula if hasattr(met_obj, '_formula') else ""
                    charge = met_obj._charge if hasattr(met_obj, '_charge') else 0
                    name = met_obj._name if hasattr(met_obj, '_name') else base_id
                    annotation = met_obj._annotation if hasattr(met_obj, '_annotation') else {}
                
                # Create cobra metabolite
                cobra_met = cobra.Metabolite(
                    id=met_id,
                    name=name,
                    formula=formula,
                    charge=charge,
                    compartment=comp_id,
                )
                if annotation:
                    cobra_met.annotation = annotation
                
                model.add_metabolites([cobra_met])
                new_mets_added += 1
            except Exception as e:
                LOGGER.warning(f"Could not add metabolite {met_id}: {e}")
    
    LOGGER.info(f"  Added {new_mets_added} new metabolites")
    
    # Add new reactions from Rhea (those in RxnList_CL but not in model)
    new_rxns_added = 0
    new_pathways = {}  # Track Rhea pathway assignments
    
    for rxn_id, rxn_obj in RxnList_CL.items():
        if rxn_id not in [r.id for r in model.reactions]:
            # Check if this is a Rhea reaction (starts with RHEA)
            base_rxn_id = rxn_id.rsplit('_', 1)[0] if '_' in rxn_id else rxn_id
            if not base_rxn_id.startswith('RHEA'):
                continue
                
            try:
                # Detect if this is a transport reaction (format: RHEA12345_c1_c2)
                # Transport reactions have format RHEAxxxxx_e_c or similar (two compartment suffixes)
                parts = rxn_id.split('_')
                is_transport_rxn = (
                    len(parts) >= 3 and 
                    parts[0].startswith('RHEA') and
                    len(parts[-1]) <= 2 and  # compartment codes are short
                    len(parts[-2]) <= 2      # two compartment codes
                )
                
                # Get reaction compartment
                comp_name = rxn_id.split('_', 1)[1] if '_' in rxn_id else 'cytosol'
                comp_id = location_dict.get(comp_name.lower(), 'c')
                
                # Check if rxn_obj has the expected methods
                if hasattr(rxn_obj, 'Substrate') and callable(rxn_obj.Substrate):
                    # This is a proper reaction object
                    substrates = rxn_obj.Substrate()
                    products = rxn_obj.Product()
                    name = rxn_obj.Name() if hasattr(rxn_obj, 'Name') and callable(rxn_obj.Name) else rxn_id
                    ec_numbers = rxn_obj.EC() if hasattr(rxn_obj, 'EC') and callable(rxn_obj.EC) else []
                    annotation = rxn_obj.Annotation() if hasattr(rxn_obj, 'Annotation') and callable(rxn_obj.Annotation) else {}
                    
                    # Build reaction metabolites dict
                    metabolites = {}
                    for stoich, _, met_id_base in substrates:
                        # For transport reactions, metabolites already have compartment suffix
                        if is_transport_rxn:
                            met_id_full = met_id_base
                        else:
                            met_id_full = f"{met_id_base}_{comp_name}"
                        if met_id_full in model.metabolites:
                            metabolites[model.metabolites.get_by_id(met_id_full)] = -float(stoich)
                        else:
                            LOGGER.warning(f"Substrate {met_id_full} not found for reaction {rxn_id}")
                    
                    for stoich, _, met_id_base in products:
                        # For transport reactions, metabolites already have compartment suffix
                        if is_transport_rxn:
                            met_id_full = met_id_base
                        else:
                            met_id_full = f"{met_id_base}_{comp_name}"
                        if met_id_full in model.metabolites:
                            metabolites[model.metabolites.get_by_id(met_id_full)] = float(stoich)
                        else:
                            LOGGER.warning(f"Product {met_id_full} not found for reaction {rxn_id}")
                    
                    if not metabolites:
                        LOGGER.warning(f"Skipping reaction {rxn_id}: no valid metabolites")
                        continue
                    
                    # Create cobra reaction
                    cobra_rxn = cobra.Reaction(
                        id=rxn_id,
                        name=f"{name}_{comp_name}",
                    )
                    cobra_rxn.add_metabolites(metabolites)
                    
                    # Set bounds based on reversibility
                    if hasattr(rxn_obj, 'Termodyn') and callable(rxn_obj.Termodyn):
                        termodyn = rxn_obj.Termodyn()
                        if termodyn == '<=>':
                            cobra_rxn.lower_bound = -1000
                            cobra_rxn.upper_bound = 1000
                        else:
                            cobra_rxn.lower_bound = 0
                            cobra_rxn.upper_bound = 1000
                    else:
                        cobra_rxn.lower_bound = -1000
                        cobra_rxn.upper_bound = 1000
                    
                    # Add annotation
                    if annotation:
                        cobra_rxn.annotation = annotation
                    if ec_numbers:
                        cobra_rxn.annotation['ec-code'] = ec_numbers
                    
                    # Set GPR if available
                    if hasattr(rxn_obj, 'GPR') and callable(rxn_obj.GPR):
                        gpr_data = rxn_obj.GPR()
                        # GPR() may return string or tuple (sGPR, GPR)
                        if isinstance(gpr_data, tuple):
                            # Use the second element (full GPR) if available
                            gpr_str = gpr_data[1] if len(gpr_data) > 1 else gpr_data[0]
                        else:
                            gpr_str = gpr_data
                        if gpr_str and isinstance(gpr_str, str) and gpr_str.strip():
                            try:
                                cobra_rxn.gene_reaction_rule = gpr_str
                            except Exception as e:
                                LOGGER.warning(f"Could not set GPR for {rxn_id}: {e}")
                    
                    model.add_reactions([cobra_rxn])
                    new_rxns_added += 1
                    
                    # Track pathway for group assignment
                    if base_rxn_id in PathNameRxn:
                        pathway_name = PathNameRxn[base_rxn_id]
                        if pathway_name not in new_pathways:
                            new_pathways[pathway_name] = []
                        new_pathways[pathway_name].append(rxn_id)
                    
            except Exception as e:
                LOGGER.warning(f"Could not add reaction {rxn_id}: {e}")
                import traceback
                LOGGER.debug(traceback.format_exc())
    
    LOGGER.info(f"  Added {new_rxns_added} new Rhea reactions")
    
    # Add new pathway groups for Rhea reactions
    for pathway_name, rxn_ids in new_pathways.items():
        # Check if group already exists
        group_id = pathway_name.replace(" ", "_")
        existing_group = None
        for g in model.groups:
            if g.id == group_id or g.name == pathway_name:
                existing_group = g
                break
        
        if existing_group:
            # Add reactions to existing group
            for rxn_id in rxn_ids:
                if rxn_id in model.reactions:
                    existing_group.add_members([model.reactions.get_by_id(rxn_id)])
        else:
            # Create new group
            new_group = cobra.core.Group(group_id, pathway_name)
            for rxn_id in rxn_ids:
                if rxn_id in model.reactions:
                    new_group.add_members([model.reactions.get_by_id(rxn_id)])
            if new_group.members:
                model.add_groups([new_group])
    
    LOGGER.info(f"  Added/updated {len(new_pathways)} pathway groups")
    
    extended_model = model
    
    LOGGER.info(f"\n  Extended model: {extended_model.id}")
    LOGGER.info(f"  Reactions: {len(extended_model.reactions)}")
    LOGGER.info(f"  Metabolites: {len(extended_model.metabolites)}")
    LOGGER.info(f"  Genes: {len(extended_model.genes)}")
    
    # Step 7: Save extended model
    step_num += 1
    if output_model_path is None:
        if standalone:
            output_model_path = os.path.join(current_dir, f"{model_id}.xml")
        else:
            base, ext = os.path.splitext(input_model_path)
            output_model_path = f"{base}_rhea_extended{ext}"
    
    LOGGER.info(f"\n[Step {step_num}] Saving {'standalone' if standalone else 'extended'} model: {output_model_path}")
    write_sbml_model(extended_model, output_model_path)
    
    LOGGER.info("\n" + "=" * 80)
    if standalone:
        LOGGER.info("RHEA STANDALONE MODEL CREATION COMPLETE")
    else:
        LOGGER.info("RHEA EXTENSION COMPLETE")
    LOGGER.info("=" * 80)
    if not standalone:
        LOGGER.info(f"  Input:  {input_model_path}")
    LOGGER.info(f"  Output: {output_model_path}")
    LOGGER.info(f"  Reactions added:  {reactions_added}")
    
    return extended_model


# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Extend an existing KEGG-based model with Rhea reactions, or create a standalone Rhea-only model",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Extend with all human Rhea reactions
    python extend_model_with_rhea.py ../models/Human_database_20251201_fixed.xml
    
    # Test with specific Rhea IDs
    python extend_model_with_rhea.py ../models/Human_database_20251201_fixed.xml \\
        --rhea-ids 10040,10112,10116,10132,10164
    
    # Specify output path
    python extend_model_with_rhea.py ../models/Human_database_20251201_fixed.xml \\
        ../models/Human_database_extended.xml
    
    # Create a standalone Rhea-only model (no KEGG model needed)
    python extend_model_with_rhea.py --standalone --output ../models/Rhea_only_model.xml
    
    # Standalone with specific Rhea IDs
    python extend_model_with_rhea.py --standalone --rhea-ids 10040,10112 --model-id "My_Rhea_Model"
        """
    )
    
    parser.add_argument(
        "input_model",
        nargs="?",
        default=None,
        help="Path to existing SBML model file (not needed with --standalone)"
    )
    
    parser.add_argument(
        "--output", "-o",
        dest="output_model",
        default=None,
        help="Output path for model (default: adds _rhea_extended suffix or uses model ID)"
    )
    
    parser.add_argument(
        "--rhea-ids",
        type=str,
        default=None,
        help="Comma-separated list of Rhea master IDs to add (default: all human Rhea)"
    )
    
    parser.add_argument(
        "--impose-locations",
        type=int,
        choices=[0, 1],
        default=1,
        help="Compartmentalization mode: 1=restricted (default), 0=unrestricted"
    )
    
    parser.add_argument(
        "--standalone",
        action="store_true",
        help="Create a new Rhea-only model without requiring an existing KEGG model"
    )
    
    parser.add_argument(
        "--model-id",
        type=str,
        default="Rhea_Human_Model",
        help="Model ID for standalone mode (default: Rhea_Human_Model)"
    )
    
    args = parser.parse_args()
    
    # Validate arguments
    if not args.standalone and not args.input_model:
        parser.error("input_model is required unless --standalone is specified")
    
    # Parse Rhea IDs if provided
    rhea_ids = None
    if args.rhea_ids:
        rhea_ids = set(args.rhea_ids.split(','))
        LOGGER.info(f"Will process {len(rhea_ids)} specific Rhea IDs: {rhea_ids}")
    
    # Run extension or standalone creation
    extend_model_with_rhea(
        input_model_path=args.input_model,
        output_model_path=args.output_model,
        rhea_ids=rhea_ids,
        impose_locations=args.impose_locations,
        standalone=args.standalone,
        model_id=args.model_id,
    )


if __name__ == "__main__":
    main()
