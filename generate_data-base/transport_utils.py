"""
Transport Reaction Utilities for Rhea Extension

This module handles the special processing of transport reactions from Rhea.
Transport reactions are identified by (in)/(out) notation in the equation string
and require different compartmentalization logic than enzymatic reactions.

Two cases are handled:
1. Membrane transport (Case A): Gene is localized to a membrane compartment
   - Direct transport: substrates from Side1, products to Side2
   
2. Non-membrane transport (Case B): Gene is localized to non-membrane compartments
   - Match existing species in model
   - Create transport pairs based on EndoA adjacency matrix
   - Assign GPR based on genes localized to each compartment pair
"""

import logging
import os
import re
from typing import Dict, List, Optional, Set, Tuple, Any

import pandas as pd

LOGGER = logging.getLogger(__name__)

# Metabolites that should be omitted when matching transport reactions
# These are ions, cofactors, energy metabolites, and small molecules
TRANSPORT_OMIT_METABOLITES = {
    # Protons and water
    'C00080',  # H+
    'C00001',  # H2O
    # Ions
    'C00076',  # Ca2+
    'C01330',  # Na+
    'C00238',  # K+
    'C00305',  # Mg2+
    'C00698',  # Cl-
    'C14818',  # Fe2+
    'C14819',  # Fe3+
    'C00070',  # Cu2+
    'C00038',  # Zn2+
    # Cofactors
    'C00003',  # NAD+
    'C00004',  # NADH
    'C00005',  # NADPH
    'C00006',  # NADP+
    'C00016',  # FAD
    'C01352',  # FADH2
    'C00010',  # CoA
    # Energy metabolites (nucleotides)
    'C00002',  # ATP
    'C00008',  # ADP
    'C00020',  # AMP
    'C00044',  # GTP
    'C00035',  # GDP
    'C00075',  # UTP
    'C00015',  # UDP
    'C00063',  # CTP
    'C00112',  # CDP
    # Phosphate
    'C00009',  # Pi (orthophosphate)
    'C00013',  # PPi (pyrophosphate)
    # Other small molecules commonly omitted
    'C00007',  # O2
    'C00011',  # CO2
    'C00014',  # NH3/NH4+
    'C00059',  # SO4 2-
    'C00042',  # Succinate (sometimes)
}


def is_transport_reaction(equation: str) -> bool:
    """
    Detect if a Rhea reaction is a transport reaction.
    
    Transport reactions in Rhea are identified by (in) and (out) notation
    in the equation string, indicating metabolites on different sides
    of a membrane.
    
    Args:
        equation: The Rhea reaction equation string
        
    Returns:
        True if this is a transport reaction, False otherwise
    """
    if not equation:
        return False
    
    # Check for (in) or (out) pattern - case insensitive
    has_in = re.search(r'\(in\)', equation, re.IGNORECASE) is not None
    has_out = re.search(r'\(out\)', equation, re.IGNORECASE) is not None
    
    return has_in or has_out


def parse_transport_equation(equation: str) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
    """
    Parse a Rhea transport equation to extract metabolites with their side (in/out).
    
    Args:
        equation: The Rhea reaction equation string
        
    Returns:
        Tuple of (substrates, products) where each is a list of (metabolite_name, side)
        side is 'in', 'out', or 'both' (if no side specified)
    """
    # This is a simplified parser - the actual implementation would need to
    # properly parse ChEBI IDs and stoichiometry from Rhea
    
    substrates = []
    products = []
    
    # Split on reaction arrow
    if ' = ' in equation:
        left, right = equation.split(' = ', 1)
    elif ' => ' in equation:
        left, right = equation.split(' => ', 1)
    elif ' <=> ' in equation:
        left, right = equation.split(' <=> ', 1)
    else:
        LOGGER.warning(f"Could not parse equation: {equation}")
        return [], []
    
    # Parse each side
    for part, target_list in [(left, substrates), (right, products)]:
        # Split by + (but careful with chemical formulas)
        # For now, simple split
        compounds = re.split(r'\s*\+\s*', part)
        
        for compound in compounds:
            compound = compound.strip()
            if not compound:
                continue
                
            # Check for (in) or (out)
            if '(in)' in compound.lower():
                side = 'in'
                compound = re.sub(r'\s*\(in\)\s*', '', compound, flags=re.IGNORECASE)
            elif '(out)' in compound.lower():
                side = 'out'
                compound = re.sub(r'\s*\(out\)\s*', '', compound, flags=re.IGNORECASE)
            else:
                side = 'both'  # No side specified, present on both sides
            
            target_list.append((compound.strip(), side))
    
    return substrates, products


class TransportBoundaryManager:
    """
    Manages transport boundary definitions from the Excel file.
    
    This class loads and provides access to:
    - Transport_Boundaries: Membrane compartments and which spaces they connect
    - Def-Compartments: Mapping from UniProt locations to model compartments
    - EndoA: Adjacency matrix of allowed transport pairs
    - Endo1a_abb: Compartment name to abbreviation mapping
    """
    
    def __init__(self, excel_path: str):
        """
        Initialize the transport boundary manager.
        
        Args:
            excel_path: Path to the ListOfCompartments Excel file
        """
        self.excel_path = excel_path
        self._transport_boundaries: Optional[pd.DataFrame] = None
        self._def_compartments: Optional[pd.DataFrame] = None
        self._endo_a: Optional[pd.DataFrame] = None
        self._abb_to_name: Optional[Dict[str, str]] = None
        self._name_to_abb: Optional[Dict[str, str]] = None
        self._membrane_compartments: Optional[Set[str]] = None
        self._endo1a_values: Optional[Set[str]] = None  # Cache for known model compartment names
        
    def _load_data(self):
        """Load all required sheets from the Excel file."""
        if self._transport_boundaries is not None:
            return  # Already loaded
            
        LOGGER.info(f"Loading transport boundary data from {self.excel_path}")
        
        # Load Transport_Boundaries
        try:
            self._transport_boundaries = pd.read_excel(
                self.excel_path, 
                sheet_name='Transport_Boundaries'
            )
            LOGGER.info(f"Loaded {len(self._transport_boundaries)} transport boundary definitions")
        except Exception as e:
            LOGGER.error(f"Failed to load Transport_Boundaries sheet: {e}")
            self._transport_boundaries = pd.DataFrame()
            
        # Load Def-Compartments
        try:
            self._def_compartments = pd.read_excel(
                self.excel_path,
                sheet_name='Def-Compartments'
            )
            LOGGER.info(f"Loaded {len(self._def_compartments)} compartment definitions")
        except Exception as e:
            LOGGER.error(f"Failed to load Def-Compartments sheet: {e}")
            self._def_compartments = pd.DataFrame()
            
        # Load EndoA matrix
        try:
            df_endo = pd.read_excel(
                self.excel_path,
                sheet_name='EndoA',
                header=None
            )
            # Parse the matrix - row 1 has headers, column 0 has row labels
            headers = df_endo.iloc[1, 1:].tolist()  # Skip 'Endo1a' label
            
            # Build adjacency dict
            self._endo_a = {}
            for i in range(2, len(df_endo)):
                row_comp = df_endo.iloc[i, 0]
                if pd.isna(row_comp):
                    continue
                self._endo_a[row_comp] = set()
                for j, col_comp in enumerate(headers):
                    val = df_endo.iloc[i, j + 1]
                    if val == 1 or val == '1':
                        self._endo_a[row_comp].add(col_comp)
                        
            LOGGER.info(f"Loaded EndoA matrix with {len(self._endo_a)} compartments")
        except Exception as e:
            LOGGER.error(f"Failed to load EndoA sheet: {e}")
            self._endo_a = {}
            
        # Load abbreviation mapping
        try:
            df_abb = pd.read_excel(
                self.excel_path,
                sheet_name='Endo1a_abb'
            )
            self._abb_to_name = {}
            self._name_to_abb = {}
            for _, row in df_abb.iterrows():
                name = str(row.iloc[0]).lower()
                abb = str(row.iloc[1]).lower()
                self._abb_to_name[abb] = name
                self._name_to_abb[name] = abb
                
            # Add 'a' = 'cell membrane' (it's in EndoA but not in Endo1a_abb)
            if 'a' not in self._abb_to_name:
                self._abb_to_name['a'] = 'cell membrane'
                self._name_to_abb['cell membrane'] = 'a'
                
            LOGGER.info(f"Loaded {len(self._abb_to_name)} compartment abbreviations")
        except Exception as e:
            LOGGER.error(f"Failed to load Endo1a_abb sheet: {e}")
            self._abb_to_name = {}
            self._name_to_abb = {}
            
        # Build set of membrane compartments
        self._membrane_compartments = set()
        if self._transport_boundaries is not None and not self._transport_boundaries.empty:
            self._membrane_compartments = set(
                self._transport_boundaries['Membrane_Abb'].str.lower().tolist()
            )
            
    @property
    def transport_boundaries(self) -> pd.DataFrame:
        """Get the transport boundaries dataframe."""
        self._load_data()
        return self._transport_boundaries
    
    @property
    def membrane_compartments(self) -> Set[str]:
        """Get set of membrane compartment abbreviations."""
        self._load_data()
        return self._membrane_compartments
    
    def is_membrane_compartment(self, compartment_abb: str) -> bool:
        """Check if a compartment abbreviation is a membrane compartment."""
        self._load_data()
        return compartment_abb.lower() in self._membrane_compartments
    
    def get_membrane_sides(self, membrane_abb: str) -> Optional[Tuple[str, str]]:
        """
        Get the two sides (compartments) that a membrane connects.
        
        Args:
            membrane_abb: Membrane compartment abbreviation (e.g., 'a')
            
        Returns:
            Tuple of (side1_abb, side2_abb) or None if not a membrane
        """
        self._load_data()
        
        membrane_abb = membrane_abb.lower()
        if self._transport_boundaries is None or self._transport_boundaries.empty:
            return None
            
        row = self._transport_boundaries[
            self._transport_boundaries['Membrane_Abb'].str.lower() == membrane_abb
        ]
        
        if row.empty:
            return None
            
        return (
            row.iloc[0]['Side1_Abb'].lower(),
            row.iloc[0]['Side2_Abb'].lower()
        )
    
    def convert_location_to_compartment(self, uniprot_location: str) -> Optional[str]:
        """
        Convert a UniProt subcellular location to a model compartment.
        
        Also handles model compartment names directly (pass-through).
        
        Args:
            uniprot_location: The UniProt subcellular location string OR model compartment name
            
        Returns:
            Model compartment name (from Endo1-b column, with Endo1-a as fallback) or None if not found
        """
        self._load_data()
        
        if self._def_compartments is None or self._def_compartments.empty:
            return None
            
        # Normalize the location
        location_lower = uniprot_location.lower().strip()
        
        # FIRST: Check if this is already a model compartment name (from Endo1-b or Endo1-a column)
        # This handles cases where we're passed "inner mitochondria" or "cell membrane" directly
        if self._endo1a_values is None:
            # Build set of known model compartment names from both Endo1-b and Endo1-a
            endo1b_col = self._def_compartments['Endo1-b'].dropna().unique()
            endo1a_col = self._def_compartments['Endo1-a'].dropna().unique()
            self._endo1a_values = set(str(v).lower().strip() for v in endo1b_col)
            self._endo1a_values.update(str(v).lower().strip() for v in endo1a_col)
        
        if location_lower in self._endo1a_values:
            # Already a model compartment name, return as-is
            return location_lower
        
        # Look up in Def-Compartments (Compartment Name -> Endo1-b, fallback to Endo1-a)
        row = self._def_compartments[
            self._def_compartments['Compartment Name'].str.lower().str.strip() == location_lower
        ]
        
        if row.empty:
            # Try partial match
            row = self._def_compartments[
                self._def_compartments['Compartment Name'].str.lower().str.contains(
                    location_lower, regex=False, na=False
                )
            ]
            
        if row.empty:
            LOGGER.debug(f"Location '{uniprot_location}' not found in Def-Compartments")
            return None
            
        # Get the Endo1-b value first (primary), then Endo1-a as fallback
        endo1b = row.iloc[0].get('Endo1-b')
        if not pd.isna(endo1b):
            return str(endo1b).lower().strip()
            
        # Fallback to Endo1-a
        endo1a = row.iloc[0].get('Endo1-a')
        if pd.isna(endo1a):
            LOGGER.debug(f"Location '{uniprot_location}' has no Endo1-b or Endo1-a mapping")
            return None
            
        return str(endo1a).lower().strip()
    
    def get_compartment_abbreviation(self, compartment_name: str) -> Optional[str]:
        """Get the abbreviation for a compartment name."""
        self._load_data()
        return self._name_to_abb.get(compartment_name.lower().strip())
    
    def get_compartment_name(self, abbreviation: str) -> Optional[str]:
        """Get the full name for a compartment abbreviation."""
        self._load_data()
        return self._abb_to_name.get(abbreviation.lower().strip())
    
    def is_transport_allowed(self, comp1: str, comp2: str) -> bool:
        """
        Check if transport between two compartments is allowed (per EndoA matrix).
        
        Args:
            comp1: First compartment abbreviation
            comp2: Second compartment abbreviation
            
        Returns:
            True if transport is allowed, False otherwise
        """
        self._load_data()
        
        comp1 = comp1.lower()
        comp2 = comp2.lower()
        
        if self._endo_a is None:
            return False
            
        # Check both directions (matrix should be symmetric)
        allowed_from_1 = self._endo_a.get(comp1, set())
        allowed_from_2 = self._endo_a.get(comp2, set())
        
        return comp2 in allowed_from_1 or comp1 in allowed_from_2
    
    def get_allowed_transport_pairs(self, compartments: Set[str]) -> List[Tuple[str, str]]:
        """
        Get all allowed transport pairs from a set of compartments.
        
        Args:
            compartments: Set of compartment abbreviations
            
        Returns:
            List of (comp1, comp2) tuples representing allowed transport pairs
        """
        self._load_data()
        
        pairs = []
        comp_list = sorted(compartments)  # Sort for consistent ordering
        
        for i, comp1 in enumerate(comp_list):
            for comp2 in comp_list[i+1:]:
                if self.is_transport_allowed(comp1, comp2):
                    pairs.append((comp1, comp2))
                    
        return pairs


class TransportReactionProcessor:
    """
    Processes transport reactions from Rhea for the metabolic model.
    
    This class implements the two-case logic for transport reactions:
    - Case A: Membrane compartment → direct transport Side1→Side2
    - Case B: Non-membrane compartments → match existing species in model
    """
    
    def __init__(self, boundary_manager: TransportBoundaryManager):
        """
        Initialize the transport reaction processor.
        
        Args:
            boundary_manager: TransportBoundaryManager instance
        """
        self.boundary_manager = boundary_manager
        
    def get_transported_metabolite(
        self,
        substrates: List[Tuple[str, str, str]],  # (coef, url, kegg_id)
        products: List[Tuple[str, str, str]],
        omit_ids: Set[str] = None
    ) -> Optional[str]:
        """
        Identify the main transported metabolite in a transport reaction.
        
        The transported metabolite is typically one that appears on both sides
        of the reaction (substrate and product) with the same stoichiometry,
        excluding small molecules like H+, H2O, ATP, etc.
        
        Args:
            substrates: List of substrate tuples (coef, url, kegg_id)
            products: List of product tuples (coef, url, kegg_id)
            omit_ids: Set of KEGG IDs to omit (defaults to TRANSPORT_OMIT_METABOLITES)
            
        Returns:
            KEGG ID of the main transported metabolite, or None if not found
        """
        if omit_ids is None:
            omit_ids = TRANSPORT_OMIT_METABOLITES
            
        # Get metabolite IDs from substrates and products
        sub_ids = {s[2] for s in substrates if len(s) >= 3}
        prod_ids = {p[2] for p in products if len(p) >= 3}
        
        # Find metabolites that appear on both sides (transported)
        common = sub_ids & prod_ids
        
        # Remove omitted metabolites
        candidates = common - omit_ids
        
        if not candidates:
            # If no common metabolites after filtering, try to identify
            # by looking at what's different between sides
            LOGGER.debug("No common metabolites found, transport may be exchange type")
            return None
            
        # Return the first candidate (could be improved with stoichiometry check)
        return list(candidates)[0] if candidates else None
    
    def process_membrane_transport(
        self,
        reaction_id: str,
        membrane_comp: str,
        substrates: List[Tuple[str, str, str]],
        products: List[Tuple[str, str, str]],
        gene_ids: List[str],
        MetList_CL: Dict[str, Any],
        MetList: Dict[str, Any],
        MetEquiv: Dict[str, str]
    ) -> Optional[Dict[str, Any]]:
        """
        Process a transport reaction where the gene is in a membrane compartment (Case A).
        
        Creates direct transport: substrates from Side1, products to Side2.
        
        Args:
            reaction_id: The reaction identifier
            membrane_comp: The membrane compartment abbreviation (e.g., 'a')
            substrates: List of substrate tuples (coef, url, kegg_id)
            products: List of product tuples (coef, url, kegg_id)
            gene_ids: List of gene IDs for GPR
            MetList_CL: Compartment-specific metabolite dictionary
            MetList: Base metabolite dictionary
            MetEquiv: Metabolite equivalence dictionary
            
        Returns:
            Dictionary with reaction data for compartmentalized reaction, or None
        """
        # Get the sides this membrane connects
        sides = self.boundary_manager.get_membrane_sides(membrane_comp)
        if sides is None:
            LOGGER.warning(f"Membrane {membrane_comp} not found in Transport_Boundaries")
            return None
            
        side1, side2 = sides
        
        LOGGER.info(
            f"Processing membrane transport {reaction_id}: "
            f"{membrane_comp} connects {side1} <-> {side2}"
        )
        
        # Create compartmentalized substrates (Side1)
        subs_cl = []
        for coef, url, kegg_id in substrates:
            # Resolve equivalences
            resolved_id = MetEquiv.get(kegg_id, kegg_id)
            met_cl_id = f"{resolved_id}_{side1}"
            subs_cl.append((coef, url, met_cl_id))
            
        # Create compartmentalized products (Side2)
        prods_cl = []
        for coef, url, kegg_id in products:
            resolved_id = MetEquiv.get(kegg_id, kegg_id)
            met_cl_id = f"{resolved_id}_{side2}"
            prods_cl.append((coef, url, met_cl_id))
            
        # Build GPR (all genes contribute to this transport)
        gpr = ' or '.join(gene_ids) if gene_ids else ''
        
        return {
            'reaction_id': f"{reaction_id}_{side1}_{side2}",
            'substrates': subs_cl,
            'products': prods_cl,
            'gpr': gpr,
            'compartment_pair': (side1, side2),
            'membrane': membrane_comp
        }
    
    def process_non_membrane_transport(
        self,
        reaction_id: str,
        gene_compartments: Dict[str, Set[str]],  # gene_id -> set of compartments
        substrates: List[Tuple[str, str, str]],
        products: List[Tuple[str, str, str]],
        MetList_CL: Dict[str, Any],
        MetList: Dict[str, Any],
        MetEquiv: Dict[str, str]
    ) -> List[Dict[str, Any]]:
        """
        Process a transport reaction where genes are in non-membrane compartments (Case B).
        
        Creates transport reactions based on:
        1. Allowed pairs from EndoA matrix
        2. Existing species in model
        3. GPR based on genes in each compartment pair
        
        Args:
            reaction_id: The reaction identifier
            gene_compartments: Dict mapping gene_id to set of compartments
            substrates: List of substrate tuples (coef, url, kegg_id)
            products: List of product tuples (coef, url, kegg_id)
            MetList_CL: Compartment-specific metabolite dictionary
            MetList: Base metabolite dictionary
            MetEquiv: Metabolite equivalence dictionary
            
        Returns:
            List of reaction dictionaries for each valid transport pair
        """
        reactions = []
        
        # Get all unique compartments from genes
        all_compartments = set()
        for comps in gene_compartments.values():
            all_compartments.update(comps)
            
        LOGGER.info(
            f"Processing non-membrane transport {reaction_id}: "
            f"compartments from genes: {all_compartments}"
        )
        
        # Get allowed transport pairs
        allowed_pairs = self.boundary_manager.get_allowed_transport_pairs(all_compartments)
        
        if not allowed_pairs:
            LOGGER.warning(f"No allowed transport pairs for {reaction_id}")
            return []
            
        LOGGER.info(f"Allowed pairs for {reaction_id}: {allowed_pairs}")
        
        # Get the main transported metabolite
        transported_met = self.get_transported_metabolite(substrates, products)
        
        if transported_met:
            LOGGER.info(f"Main transported metabolite: {transported_met}")
        else:
            # If no clear transported metabolite, we still process but log warning
            LOGGER.warning(
                f"Could not identify main transported metabolite for {reaction_id}"
            )
        
        # For each allowed pair, create transport reactions
        # Note: For Rhea extension, we create transport reactions even if the metabolite
        # doesn't exist yet in the model, because the metabolite will be added with the reaction.
        # We log when metabolites already exist for informational purposes.
        for comp1, comp2 in allowed_pairs:
            if transported_met:
                resolved_met = MetEquiv.get(transported_met, transported_met)
                met_comp1 = f"{resolved_met}_{comp1}"
                met_comp2 = f"{resolved_met}_{comp2}"
                
                # Log whether metabolite exists (informational only, not a skip condition)
                exists_in_comp1 = met_comp1 in MetList_CL
                exists_in_comp2 = met_comp2 in MetList_CL
                
                if exists_in_comp1 or exists_in_comp2:
                    LOGGER.debug(
                        f"Metabolite {transported_met} exists in: "
                        f"{comp1}={exists_in_comp1}, {comp2}={exists_in_comp2}"
                    )
            
            # Build GPR for this pair: genes that are in comp1 OR comp2
            genes_for_pair = []
            for gene_id, gene_comps in gene_compartments.items():
                if comp1 in gene_comps or comp2 in gene_comps:
                    genes_for_pair.append(gene_id)
                    
            gpr = ' or '.join(genes_for_pair) if genes_for_pair else ''
            
            # Create compartmentalized reaction
            # Substrates in comp1, products in comp2
            subs_cl = []
            for coef, url, kegg_id in substrates:
                resolved_id = MetEquiv.get(kegg_id, kegg_id)
                met_cl_id = f"{resolved_id}_{comp1}"
                subs_cl.append((coef, url, met_cl_id))
                
            prods_cl = []
            for coef, url, kegg_id in products:
                resolved_id = MetEquiv.get(kegg_id, kegg_id)
                met_cl_id = f"{resolved_id}_{comp2}"
                prods_cl.append((coef, url, met_cl_id))
            
            reactions.append({
                'reaction_id': f"{reaction_id}_{comp1}_{comp2}",
                'substrates': subs_cl,
                'products': prods_cl,
                'gpr': gpr,
                'compartment_pair': (comp1, comp2),
                'genes': genes_for_pair
            })
            
            LOGGER.info(
                f"Created transport reaction {reaction_id}_{comp1}_{comp2} "
                f"with GPR: {gpr}"
            )
            
        return reactions
    
    def process_single_location_transport(
        self,
        reaction_id: str,
        single_comp: str,
        substrates: List[Tuple[str, str, str]],
        products: List[Tuple[str, str, str]],
        gene_ids: List[str],
        MetList_CL: Dict[str, Any],
        MetList: Dict[str, Any],
        MetEquiv: Dict[str, str]
    ) -> List[Dict[str, Any]]:
        """
        Process a transport reaction where genes have only ONE location.
        
        Logic:
        1. If the location is a membrane compartment:
           - Create transport between Side1 ↔ Side2 of that membrane
        2. If the location is NOT a membrane compartment:
           - Look up the EndoA adjacency matrix
           - Find all compartments that can connect to this location
           - Exclude membrane compartments as destinations
           - Only create transport if the metabolite already exists in the target compartment
        
        Args:
            reaction_id: The reaction identifier
            single_comp: The single compartment abbreviation (e.g., 'c', 'm', 'a')
            substrates: List of substrate tuples (coef, url, kegg_id)
            products: List of product tuples (coef, url, kegg_id)
            gene_ids: List of gene IDs for GPR
            MetList_CL: Compartment-specific metabolite dictionary
            MetList: Base metabolite dictionary
            MetEquiv: Metabolite equivalence dictionary
            
        Returns:
            List of reaction dictionaries for each valid transport
        """
        reactions = []
        
        # Check if this is a membrane compartment
        if self.boundary_manager.is_membrane_compartment(single_comp):
            # Case: Membrane compartment - create direct transport Side1 ↔ Side2
            result = self.process_membrane_transport(
                reaction_id=reaction_id,
                membrane_comp=single_comp,
                substrates=substrates,
                products=products,
                gene_ids=gene_ids,
                MetList_CL=MetList_CL,
                MetList=MetList,
                MetEquiv=MetEquiv
            )
            if result:
                reactions.append(result)
                
            LOGGER.info(
                f"Single-location transport {reaction_id}: "
                f"membrane {single_comp} -> processed as membrane transport"
            )
            return reactions
        
        # Case: Non-membrane compartment - use adjacency matrix
        LOGGER.info(
            f"Single-location transport {reaction_id}: "
            f"non-membrane {single_comp} -> checking adjacency matrix"
        )
        
        # Get the transported metabolite(s) to check existence
        transported_met = self.get_transported_metabolite(substrates, products)
        if transported_met:
            LOGGER.info(f"Main transported metabolite: {transported_met}")
        
        # Get all metabolites from the reaction (for existence check)
        all_met_ids = set()
        for _, _, kegg_id in substrates + products:
            resolved_id = MetEquiv.get(kegg_id, kegg_id)
            all_met_ids.add(resolved_id)
        
        # Get adjacent compartments from EndoA matrix
        self.boundary_manager._load_data()
        
        if single_comp not in self.boundary_manager._endo_a:
            LOGGER.warning(f"Compartment {single_comp} not found in EndoA matrix")
            return []
        
        # Find all compartments that this one can connect to
        adjacent_comps = list(self.boundary_manager._endo_a[single_comp])
        
        LOGGER.debug(f"Adjacent compartments to {single_comp}: {adjacent_comps}")
        
        # Build GPR (all genes contribute to all transports)
        gpr = ' or '.join(gene_ids) if gene_ids else ''
        
        # Process each adjacent compartment
        # Track which target compartments we've already processed to avoid duplicates
        processed_targets = set()
        
        for adj_comp in adjacent_comps:
            # Determine the actual target compartment
            if self.boundary_manager.is_membrane_compartment(adj_comp):
                # Adjacent is a membrane - transport to the OTHER side of the membrane
                sides = self.boundary_manager.get_membrane_sides(adj_comp)
                if sides is None:
                    LOGGER.warning(f"Membrane {adj_comp} not found in Transport_Boundaries")
                    continue
                    
                side1, side2 = sides
                
                # Skip if both sides are the same (no actual transport possible)
                if side1 == side2:
                    LOGGER.debug(
                        f"Skipping membrane {adj_comp}: both sides are {side1} (no transport)"
                    )
                    continue
                
                # Transport to the side that is NOT our current location
                if single_comp == side1:
                    target_comp = side2
                elif single_comp == side2:
                    target_comp = side1
                else:
                    # Current location is not directly on either side of the membrane
                    # This shouldn't happen with proper adjacency, but handle it
                    LOGGER.debug(
                        f"Compartment {single_comp} is adjacent to membrane {adj_comp} "
                        f"but not on Side1 ({side1}) or Side2 ({side2})"
                    )
                    continue
                    
                LOGGER.debug(
                    f"Adjacent membrane {adj_comp}: sides={side1},{side2}, "
                    f"current={single_comp} -> target={target_comp}"
                )
            else:
                # Adjacent is not a membrane - transport directly to it
                target_comp = adj_comp
            
            # Skip if we've already processed this target
            if target_comp in processed_targets:
                LOGGER.debug(f"Skipping duplicate target {target_comp}")
                continue
            processed_targets.add(target_comp)
            
            # Check if any transported metabolite exists in the target compartment
            exists_in_target = False
            for met_id in all_met_ids:
                met_cl_id = f"{met_id}_{target_comp}"
                if met_cl_id in MetList_CL:
                    exists_in_target = True
                    LOGGER.debug(f"Found {met_cl_id} in MetList_CL")
                    break
            
            if not exists_in_target:
                LOGGER.debug(
                    f"Skipping {single_comp} -> {target_comp}: "
                    f"no transported metabolites exist in {target_comp}"
                )
                continue
            
            # Create compartmentalized reaction
            # Substrates in single_comp, products in target_comp
            subs_cl = []
            for coef, url, kegg_id in substrates:
                resolved_id = MetEquiv.get(kegg_id, kegg_id)
                met_cl_id = f"{resolved_id}_{single_comp}"
                subs_cl.append((coef, url, met_cl_id))
                
            prods_cl = []
            for coef, url, kegg_id in products:
                resolved_id = MetEquiv.get(kegg_id, kegg_id)
                met_cl_id = f"{resolved_id}_{target_comp}"
                prods_cl.append((coef, url, met_cl_id))
            
            reactions.append({
                'reaction_id': f"{reaction_id}_{single_comp}_{target_comp}",
                'substrates': subs_cl,
                'products': prods_cl,
                'gpr': gpr,
                'compartment_pair': (single_comp, target_comp),
                'membrane': adj_comp if self.boundary_manager.is_membrane_compartment(adj_comp) else None
            })
            
            LOGGER.info(
                f"Created single-location transport {reaction_id}_{single_comp}_{target_comp} "
                f"(via membrane {adj_comp})" if self.boundary_manager.is_membrane_compartment(adj_comp) 
                else f"Created single-location transport {reaction_id}_{single_comp}_{target_comp}"
            )
        
        if not reactions:
            LOGGER.warning(
                f"No valid transport pairs for {reaction_id} from {single_comp}: "
                f"no metabolites exist in adjacent compartments"
            )
        
        return reactions

    def process_transport_reaction(
        self,
        reaction_id: str,
        equation: str,
        gene_locations: Dict[str, List[str]],  # gene_id -> list of UniProt locations
        substrates: List[Tuple[str, str, str]],
        products: List[Tuple[str, str, str]],
        MetList_CL: Dict[str, Any],
        MetList: Dict[str, Any],
        MetEquiv: Dict[str, str]
    ) -> List[Dict[str, Any]]:
        """
        Main entry point for processing a transport reaction.
        
        Determines whether to use Case A (membrane) or Case B (non-membrane)
        based on gene locations.
        
        Args:
            reaction_id: The reaction identifier
            equation: The reaction equation (for logging)
            gene_locations: Dict mapping gene_id to list of UniProt locations
            substrates: List of substrate tuples
            products: List of product tuples
            MetList_CL: Compartment-specific metabolite dictionary
            MetList: Base metabolite dictionary
            MetEquiv: Metabolite equivalence dictionary
            
        Returns:
            List of reaction dictionaries for compartmentalized reactions
        """
        LOGGER.info(f"Processing transport reaction {reaction_id}")
        LOGGER.debug(f"Equation: {equation}")
        LOGGER.debug(f"Gene locations received: {gene_locations}")
        
        # Filter out generic/ambiguous compartment assignments when more specific ones are available
        # The issue: "membrane" (generic) maps to "cell membrane" in the Excel, but when a gene
        # is also in specific mitochondrial locations, we should ignore the generic "cell membrane"
        # because it came from an ambiguous "membrane" annotation, not actual plasma membrane localization.
        # 
        # These are model compartment names (Endo1-b values) that can come from ambiguous source terms:
        AMBIGUOUS_MEMBRANE_COMPARTMENTS = {'cell membrane'}  # Can come from generic "membrane" term
        
        filtered_gene_locations: Dict[str, List[str]] = {}
        for gene_id, locations in gene_locations.items():
            # Check if we have both ambiguous membrane and specific non-membrane locations
            has_ambiguous = any(loc.lower() in AMBIGUOUS_MEMBRANE_COMPARTMENTS for loc in locations)
            specific_locations = [loc for loc in locations if loc.lower() not in AMBIGUOUS_MEMBRANE_COMPARTMENTS]
            
            if has_ambiguous and specific_locations:
                # Gene has both ambiguous membrane AND specific locations - keep only specific
                LOGGER.debug(f"Gene {gene_id}: filtering out ambiguous 'cell membrane', keeping: {specific_locations}")
                filtered_gene_locations[gene_id] = specific_locations
            else:
                # No ambiguity or only ambiguous - keep all
                filtered_gene_locations[gene_id] = locations
        
        # Convert gene locations to model compartments
        gene_compartments: Dict[str, Set[str]] = {}
        membrane_genes: Dict[str, str] = {}  # gene_id -> membrane_comp
        
        for gene_id, locations in filtered_gene_locations.items():
            gene_compartments[gene_id] = set()
            
            for location in locations:
                # Convert UniProt location to model compartment
                model_comp = self.boundary_manager.convert_location_to_compartment(location)
                
                if model_comp is None:
                    LOGGER.debug(f"Gene {gene_id} location '{location}' has no mapping")
                    continue
                    
                # Get the abbreviation
                abb = self.boundary_manager.get_compartment_abbreviation(model_comp)
                if abb is None:
                    abb = model_comp  # Use as-is if no abbreviation found
                    
                gene_compartments[gene_id].add(abb)
                
                # Check if this is a membrane compartment
                if self.boundary_manager.is_membrane_compartment(abb):
                    membrane_genes[gene_id] = abb
                    
        LOGGER.debug(f"Gene compartments (model): {gene_compartments}")
        LOGGER.debug(f"Membrane genes: {membrane_genes}")
        
        # Count total unique compartments across all genes
        all_compartments = set()
        for comps in gene_compartments.values():
            all_compartments.update(comps)
        
        LOGGER.debug(f"All unique compartments: {all_compartments}")
        
        # Determine which case to use
        if membrane_genes:
            # Case A: At least one gene is in a membrane compartment
            # Use the first membrane gene's location
            # (Could be improved to handle multiple membrane locations)
            gene_id = list(membrane_genes.keys())[0]
            membrane_comp = membrane_genes[gene_id]
            
            result = self.process_membrane_transport(
                reaction_id=reaction_id,
                membrane_comp=membrane_comp,
                substrates=substrates,
                products=products,
                gene_ids=list(gene_locations.keys()),
                MetList_CL=MetList_CL,
                MetList=MetList,
                MetEquiv=MetEquiv
            )
            
            return [result] if result else []
        
        elif len(all_compartments) == 1:
            # Case: Single location (non-membrane) - use adjacency matrix
            single_comp = list(all_compartments)[0]
            LOGGER.info(
                f"Transport {reaction_id}: Single non-membrane location '{single_comp}' "
                f"-> using adjacency matrix"
            )
            
            return self.process_single_location_transport(
                reaction_id=reaction_id,
                single_comp=single_comp,
                substrates=substrates,
                products=products,
                gene_ids=list(gene_locations.keys()),
                MetList_CL=MetList_CL,
                MetList=MetList,
                MetEquiv=MetEquiv
            )
        
        else:
            # Case B: Multiple non-membrane compartments - use existing logic
            return self.process_non_membrane_transport(
                reaction_id=reaction_id,
                gene_compartments=gene_compartments,
                substrates=substrates,
                products=products,
                MetList_CL=MetList_CL,
                MetList=MetList,
                MetEquiv=MetEquiv
            )


# Convenience function for module-level access
_boundary_manager_instance: Optional[TransportBoundaryManager] = None


def get_boundary_manager(excel_path: str = None) -> TransportBoundaryManager:
    """
    Get the singleton TransportBoundaryManager instance.
    
    Args:
        excel_path: Path to the Excel file (required on first call)
        
    Returns:
        TransportBoundaryManager instance
    """
    global _boundary_manager_instance
    
    if _boundary_manager_instance is None:
        if excel_path is None:
            raise ValueError("excel_path must be provided on first call")
        _boundary_manager_instance = TransportBoundaryManager(excel_path)
        
    return _boundary_manager_instance


def get_transport_processor(excel_path: str = None) -> TransportReactionProcessor:
    """
    Get a TransportReactionProcessor instance.
    
    Args:
        excel_path: Path to the Excel file
        
    Returns:
        TransportReactionProcessor instance
    """
    manager = get_boundary_manager(excel_path)
    return TransportReactionProcessor(manager)
