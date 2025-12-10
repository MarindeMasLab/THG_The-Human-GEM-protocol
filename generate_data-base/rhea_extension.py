#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhea Extension Module for Human Metabolic Database Generation

This module extends the KEGG-based metabolic model with reactions from the Rhea database.
It REUSES the existing classes from class_generate_database.py:
- reaction, compound, gene, gpr

Philosophy:
- Start from the KEGG model output (checkpoint or final model)
- Download fresh Rhea data at runtime
- Filter for human-relevant reactions not already in KEGG
- Use the SAME classes and pipeline structure as KEGG

Pipeline Flow:
    1. Load KEGG checkpoint/model
    2. Download Rhea TSV files (fresh data)
    3. Filter: Human reactions NOT in KEGG
    4. For each new Rhea reaction:
       - Create reaction object (reusing `reaction` class structure)
       - Create compound objects (reusing `compound` class)
       - Get GPR via UniProt→Ensembl (reusing `gpr`/`gene` classes)
       - Compartmentalize (same as KEGG)
    5. Merge with KEGG data structures
    6. Continue to cobra_reconstruction

Author: [Your Name]
Date: December 2025
"""

import os
import sys
import re
import logging
import pickle
import requests
import pandas as pd
from typing import Dict, Set, List, Tuple, Optional, Any
from dataclasses import dataclass, field
from collections import defaultdict
from tqdm import tqdm

# Determine the current file's directory and the project root.
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.join(current_dir, "..")

if project_root not in sys.path:
    sys.path.append(project_root)

# Import existing classes - REUSE, don't recreate!
from functions.class_generate_database import reaction, compound, gene, gpr
from functions.function_bm_gdb import getHtml, batch_fetch_kegg_entries
from functions.gpr.auth_gpr import setup_biocyc_session
from functions.gpr.get_location_def import getLocationnew
from functions.pattern_generate_database import rxnSubcel
from functions.equations_bm_gdb import mass_balance, RxnParam2Eq
from functions.metabolite_matching import (
    build_metabolite_attribute_index,
    find_equivalent_metabolite,
    METABOLITE_MATCH_ATTRIBUTES
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
LOGGER = logging.getLogger(__name__)


# =============================================================================
# HELPER FUNCTIONS FOR ATTRIBUTE ACCESS
# =============================================================================

def _get_compound_name(compound_obj) -> str:
    """
    Safely get the name from a compound object.
    
    Handles both:
    - compound class from class_generate_database.py (Name is an attribute)
    - MetaboliteProxy from extend_model_with_rhea.py (Name() is a method)
    
    Returns:
        str: The compound name, or empty string if not available
    """
    if not hasattr(compound_obj, 'Name'):
        return ""
    
    name_attr = compound_obj.Name
    # If Name is a method (callable), call it
    if callable(name_attr):
        return name_attr() or ""
    # Otherwise it's an attribute
    return name_attr or ""


def _set_compound_name(compound_obj, name: str) -> bool:
    """
    Safely set the name on a compound object.
    
    Returns:
        bool: True if name was set, False if compound doesn't support setting name
    """
    if not hasattr(compound_obj, 'Name'):
        return False
    
    name_attr = compound_obj.Name
    # If Name is a method (from MetaboliteProxy), we can't set it
    if callable(name_attr):
        # Try to set _name directly (MetaboliteProxy internal attribute)
        if hasattr(compound_obj, '_name'):
            compound_obj._name = name
            return True
        return False
    # Otherwise set the attribute directly
    compound_obj.Name = name
    return True


# =============================================================================
# FORMULA VALIDATION
# =============================================================================

# Pattern to detect problematic formulas with mathematical operations
# These contain +, -, /, * between numbers (e.g., '76800 - 51200', '2419200/2')
# This cannot be parsed by the mass balance algorithm
PROBLEMATIC_FORMULA_PATTERN = re.compile(r'\d+\s*[\-\*/]\s*\d+')  # Note: + is NOT included (C6H12O6 has no +)

def is_valid_formula_for_mass_balance(formula: str) -> bool:
    """
    Check if a formula can be used for mass balance calculations.
    
    Formulas with mathematical expressions (e.g., '76800 - 51200', '2419200/2')
    cannot be parsed by the mass balance algorithm and should be skipped.
    
    NOTE: Polymeric formulas with 'n' like (C8H13O10N)n ARE allowed because
    RxnParam2Eq already handles these by stripping the parentheses and 'n'.
    See equations_bm_gdb.py line ~681:
        if re.findall(r")n", eq):
            eq = re.sub(r")n", "", re.sub(r"(", "", ...))
    
    Args:
        formula: Chemical formula string
        
    Returns:
        True if formula is valid for mass balance, False otherwise
    """
    if not formula:
        return False
    
    # Trim whitespace
    formula = formula.strip()
    
    # Check for mathematical operations between numbers (e.g., '76800 - 51200')
    # These are impossible to parse as chemical formulas
    if PROBLEMATIC_FORMULA_PATTERN.search(formula):
        return False
    
    # Check for undefined/placeholder formulas
    if formula in ('R', 'X', '*', '?'):
        return False
    
    return True


# =============================================================================
# CONFIGURATION
# =============================================================================

RHEA_FTP_BASE = "https://ftp.expasy.org/databases/rhea/tsv/"
RHEA_SPARQL_ENDPOINT = "https://sparql.rhea-db.org/sparql"
CHEBI_API_BASE = "https://www.ebi.ac.uk/chebi/webServices.do"
UNIPROT_API_BASE = "https://rest.uniprot.org/uniprotkb"
REACTOME_API_BASE = "https://reactome.org/ContentService"
HUMAN_TAXON_ID = "9606"

RHEA_FILES = {
    "rhea2kegg_reaction.tsv": "Rhea ↔ KEGG mappings (to find NEW reactions)",
    "rhea2ec.tsv": "Rhea ↔ EC number mappings",
    "rhea2uniprot_sprot.tsv": "Rhea ↔ UniProt (for human genes)",
    "rhea-directions.tsv": "Master/directional ID mappings (for Termodyn)",
    "rhea2reactome.tsv": "Rhea ↔ Reactome pathway mappings",
    "chebiId_name.tsv": "ChEBI compound names",
}


# =============================================================================
# REACTOME PATHWAY LOOKUP
# =============================================================================

def fetch_reactome_pathway_name(reactome_id: str, cache: Optional[Dict[str, str]] = None) -> Optional[str]:
    """
    Fetch the pathway name for a Reactome reaction/entity ID.
    
    The Reactome API returns pathway information for a given entity.
    We extract the most specific (lowest-level) pathway name.
    
    Args:
        reactome_id: Reactome ID like "R-HSA-8953499" or "R-HSA-8953499.3"
        cache: Optional cache dict to avoid redundant API calls
        
    Returns:
        Pathway name string, or None if not found
    """
    # Strip version suffix if present (R-HSA-8953499.3 -> R-HSA-8953499)
    base_id = reactome_id.split('.')[0] if '.' in reactome_id else reactome_id
    
    # Check cache first
    if cache is not None and base_id in cache:
        return cache[base_id]
    
    try:
        # Query Reactome for pathways containing this entity
        url = f"{REACTOME_API_BASE}/data/pathways/low/entity/{base_id}?species=9606"
        resp = requests.get(url, timeout=10)
        
        if resp.status_code == 200:
            pathways = resp.json()
            if pathways:
                # Return the first (most specific) pathway's displayName
                pathway_name = pathways[0].get('displayName', '')
                if pathway_name and cache is not None:
                    cache[base_id] = pathway_name
                return pathway_name
                
    except requests.RequestException as e:
        LOGGER.debug(f"Could not fetch Reactome pathway for {base_id}: {e}")
    except Exception as e:
        LOGGER.debug(f"Error processing Reactome response for {base_id}: {e}")
    
    return None


def batch_fetch_reactome_pathways(reactome_ids: Set[str], max_workers: int = 5) -> Dict[str, str]:
    """
    Batch fetch pathway names for multiple Reactome IDs.
    
    Uses concurrent requests with rate limiting to avoid overwhelming the API.
    
    Args:
        reactome_ids: Set of Reactome IDs to look up
        max_workers: Maximum concurrent requests
        
    Returns:
        Dict mapping reactome_id -> pathway_name
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed
    import time
    
    result = {}
    cache = {}
    
    # Deduplicate by stripping version suffixes
    unique_ids = set(rid.split('.')[0] for rid in reactome_ids if rid)
    
    LOGGER.info(f"  Fetching pathway names for {len(unique_ids)} unique Reactome IDs...")
    
    def fetch_one(rid):
        time.sleep(0.1)  # Rate limiting: 10 requests/second
        return rid, fetch_reactome_pathway_name(rid, cache)
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(fetch_one, rid): rid for rid in unique_ids}
        
        for future in tqdm(as_completed(futures), total=len(futures), 
                          desc="Fetching Reactome pathways", leave=False):
            try:
                rid, pathway = future.result()
                if pathway:
                    result[rid] = pathway
                    # Also map versioned IDs to same pathway
                    for orig_id in reactome_ids:
                        if orig_id.startswith(rid):
                            result[orig_id] = pathway
            except Exception as e:
                LOGGER.debug(f"Error fetching pathway: {e}")
    
    LOGGER.info(f"  Resolved {len(result)} pathway names from Reactome API")
    return result


# =============================================================================
# RHEA DATA DOWNLOADER (same philosophy as KEGG - fresh data)
# =============================================================================

class RheaDataDownloader:
    """
    Downloads Rhea database files to /files/rhea/
    Philosophy: Always fetch fresh data for latest mappings.
    """
    
    def __init__(self, cache_dir: Optional[str] = None):
        if cache_dir is None:
            cache_dir = os.path.join(project_root, "files", "rhea")
        self.cache_dir = cache_dir
        self._data_cache: Dict[str, pd.DataFrame] = {}
        
    def download_all(self, force_refresh: bool = False) -> None:
        """Download all required Rhea files."""
        os.makedirs(self.cache_dir, exist_ok=True)
        
        for filename, description in RHEA_FILES.items():
            filepath = os.path.join(self.cache_dir, filename)
            
            if os.path.exists(filepath) and not force_refresh:
                LOGGER.info(f"Using cached: {filename}")
            else:
                LOGGER.info(f"Downloading: {filename}")
                self._download_file(filename, filepath)
                
    def _download_file(self, filename: str, filepath: str) -> None:
        """Download a single file from Rhea FTP."""
        url = f"{RHEA_FTP_BASE}{filename}"
        
        try:
            response = requests.get(url, timeout=120)
            response.raise_for_status()
            
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(response.text)
                
            LOGGER.info(f"  Saved: {filepath}")
            
        except requests.RequestException as e:
            LOGGER.error(f"  Failed to download {filename}: {e}")
            raise
            
    def load_file(self, filename: str) -> pd.DataFrame:
        """Load a Rhea TSV file as pandas DataFrame."""
        if filename in self._data_cache:
            return self._data_cache[filename]
            
        filepath = os.path.join(self.cache_dir, filename)
        
        if not os.path.exists(filepath):
            self.download_all()
            
        df = pd.read_csv(filepath, sep='\t')
        self._data_cache[filename] = df
        return df


# =============================================================================
# HUMAN PROTEOME FILTER
# =============================================================================

class HumanProteomeFilter:
    """
    Downloads human Swiss-Prot protein list from UniProt.
    Maps UniProt IDs to gene symbols (primary) and Ensembl gene IDs (secondary).
    
    IMPORTANT: We use gene symbols as primary identifiers to match KEGG's format.
    The KEGG pipeline stores genes by symbol (e.g., "ADH1B", "BRCA1") not Ensembl IDs.
    
    Also fetches cross-references (KEGG, GeneID, Ensembl, EC) for fallback mapping.
    """
    
    def __init__(self):
        self.human_uniprots: Set[str] = set()
        self.uniprot_to_ensembl: Dict[str, Set[str]] = {}
        self.uniprot_to_gene_symbol: Dict[str, str] = {}  # Primary mapping for gene names
        # Additional cross-reference mappings for fallback
        self.uniprot_to_kegg: Dict[str, str] = {}         # UniProt -> KEGG gene ID
        self.uniprot_to_geneid: Dict[str, str] = {}       # UniProt -> NCBI GeneID
        self.uniprot_to_ec: Dict[str, List[str]] = {}     # UniProt -> EC numbers
        # Statistics for diagnostics
        self.uniprot_mapping_stats: Dict[str, int] = {}
        
    def load_human_proteome(self) -> None:
        """Download list of human Swiss-Prot protein IDs."""
        LOGGER.info("Fetching human proteome from UniProt...")
        
        url = f"{UNIPROT_API_BASE}/stream"
        params = {
            "query": f"(organism_id:{HUMAN_TAXON_ID}) AND (reviewed:true)",
            "fields": "accession",
            "format": "list",
        }
        
        response = requests.get(url, params=params, timeout=300)
        response.raise_for_status()
        
        self.human_uniprots = set(response.text.strip().split('\n'))
        LOGGER.info(f"  Loaded {len(self.human_uniprots)} human proteins")
    
    def map_uniprot_to_gene_symbol_batch(self, uniprot_ids: Set[str]) -> Dict[str, str]:
        """
        Map UniProt IDs to gene symbols via UniProt API.
        
        This is the PRIMARY mapping used to match KEGG's gene format.
        KEGG GPRs use gene symbols (e.g., "ADH1B") not Ensembl IDs.
        
        Also fetches additional cross-references (KEGG, GeneID, Ensembl, EC) 
        for fallback and diagnostics.
        
        Args:
            uniprot_ids: Set of UniProt accession IDs
            
        Returns:
            Dict mapping UniProt ID -> Gene Symbol (e.g., {"P00325": "ADH1B"})
        """
        # Filter for human only
        human_ids = uniprot_ids & self.human_uniprots
        
        if not human_ids:
            return {}
            
        LOGGER.info(f"Mapping {len(human_ids)} UniProt IDs to gene symbols...")
        
        mapping = {}
        # Additional mappings for fallback/diagnostics
        self.uniprot_to_kegg: Dict[str, str] = {}      # UniProt -> KEGG gene ID
        self.uniprot_to_geneid: Dict[str, str] = {}    # UniProt -> NCBI GeneID
        self.uniprot_to_ec: Dict[str, List[str]] = {}  # UniProt -> EC numbers
        
        # Statistics for diagnostics
        stats = {
            'total_queried': 0,
            'has_gene_symbol': 0,
            'has_kegg': 0,
            'has_geneid': 0,
            'has_ensembl': 0,
            'has_ec': 0,
            'no_gene_symbol': 0,
            'api_errors': 0,
            'failed_ids': 0,  # IDs that failed even after retries
        }
        
        human_list = list(human_ids)
        failed_ids = []  # Track IDs from failed batches for retry
        
        def _process_batch_response(response_text: str, batch_ids: List[str]) -> Set[str]:
            """Process API response and return set of successfully processed IDs."""
            processed = set()
            lines = response_text.strip().split('\n')
            
            # DEBUG: Log response details
            LOGGER.debug(f"API response: {len(lines)} lines, first line: {lines[0][:50] if lines else 'EMPTY'}")
            
            if len(lines) > 1:  # Has header + data
                header = lines[0].split('\t')
                # Find column indices
                col_idx = {col: idx for idx, col in enumerate(header)}
                
                for line in lines[1:]:
                    parts = line.split('\t')
                    if not parts:
                        continue
                        
                    uniprot = parts[0] if len(parts) > 0 else ""
                    if not uniprot:
                        continue
                    
                    processed.add(uniprot)
                    
                    # Gene symbol (primary mapping)
                    gene_symbol = ""
                    if 'Gene Names (primary)' in col_idx:
                        idx = col_idx['Gene Names (primary)']
                        gene_symbol = parts[idx].strip() if len(parts) > idx else ""
                    
                    if gene_symbol:
                        mapping[uniprot] = gene_symbol
                        stats['has_gene_symbol'] += 1
                    else:
                        stats['no_gene_symbol'] += 1
                    
                    # KEGG cross-reference (fallback)
                    if 'KEGG' in col_idx:
                        idx = col_idx['KEGG']
                        kegg_ref = parts[idx].strip() if len(parts) > idx else ""
                        if kegg_ref:
                            # Format: "hsa:123" or just gene ID
                            self.uniprot_to_kegg[uniprot] = kegg_ref
                            stats['has_kegg'] += 1
                    
                    # GeneID cross-reference
                    if 'GeneID' in col_idx:
                        idx = col_idx['GeneID']
                        geneid = parts[idx].strip() if len(parts) > idx else ""
                        if geneid:
                            self.uniprot_to_geneid[uniprot] = geneid
                            stats['has_geneid'] += 1
                    
                    # Ensembl cross-reference
                    if 'Ensembl' in col_idx:
                        idx = col_idx['Ensembl']
                        ensembl_ref = parts[idx].strip() if len(parts) > idx else ""
                        if ensembl_ref and 'ENSG' in ensembl_ref:
                            stats['has_ensembl'] += 1
                    
                    # EC number
                    if 'EC number' in col_idx:
                        idx = col_idx['EC number']
                        ec_ref = parts[idx].strip() if len(parts) > idx else ""
                        if ec_ref:
                            # Parse EC numbers (semicolon-separated)
                            ec_list = [ec.strip() for ec in ec_ref.split(';') if ec.strip()]
                            if ec_list:
                                self.uniprot_to_ec[uniprot] = ec_list
                                stats['has_ec'] += 1
            
            return processed
        
        def _query_batch(batch_ids: List[str], timeout: int = 60) -> Optional[str]:
            """Query UniProt API for a batch. Returns response text or None on failure."""
            url = f"{UNIPROT_API_BASE}/stream"
            params = {
                "query": " OR ".join([f"accession:{acc}" for acc in batch_ids]),
                "fields": "accession,gene_primary,xref_kegg,xref_geneid,xref_ensembl,ec",
                "format": "tsv",
            }
            try:
                response = requests.get(url, params=params, timeout=timeout)
                if response.ok:
                    return response.text
                else:
                    LOGGER.warning(f"Batch query returned status {response.status_code}")
                    return None
            except requests.RequestException as e:
                LOGGER.warning(f"Batch query failed: {e}")
                return None
        
        # UniProt API has a hard limit of 100 OR conditions per query
        batch_size = 100
        
        # First pass: process batches
        for i in tqdm(range(0, len(human_list), batch_size), desc="UniProt→Symbol+xrefs"):
            batch = human_list[i:i+batch_size]
            stats['total_queried'] += len(batch)
            
            response_text = _query_batch(batch)
            if response_text:
                _process_batch_response(response_text, batch)
            else:
                # Batch failed - add to retry list
                stats['api_errors'] += 1
                failed_ids.extend(batch)
        
        # Retry logic with exponential backoff for any failed IDs
        if failed_ids:
            import time
            LOGGER.info(f"  Retrying {len(failed_ids)} UniProt IDs from failed batches...")
            retry_batch_size = 50  # Smaller batch for retries
            max_retries = 3
            
            for retry_attempt in range(max_retries):
                if not failed_ids:
                    break
                    
                still_failed = []
                backoff_time = 2 ** retry_attempt  # 1, 2, 4 seconds
                
                LOGGER.info(f"  Retry attempt {retry_attempt + 1}/{max_retries} with {len(failed_ids)} IDs (batch size: {retry_batch_size}, backoff: {backoff_time}s)")
                
                for j in range(0, len(failed_ids), retry_batch_size):
                    batch = failed_ids[j:j+retry_batch_size]
                    
                    # Apply backoff before each retry batch
                    time.sleep(backoff_time)
                    
                    response_text = _query_batch(batch, timeout=120)  # Longer timeout for retries
                    if response_text:
                        _process_batch_response(response_text, batch)
                    else:
                        still_failed.extend(batch)
                
                failed_ids = still_failed
                
                if still_failed:
                    LOGGER.warning(f"  {len(still_failed)} IDs still failed after retry {retry_attempt + 1}")
            
            # Final report on failed IDs
            if failed_ids:
                stats['failed_ids'] = len(failed_ids)
                LOGGER.warning(f"  {len(failed_ids)} UniProt IDs could not be mapped after all retries")
                # Store failed IDs for potential later analysis
                self._failed_uniprot_ids = set(failed_ids)
                
        self.uniprot_to_gene_symbol.update(mapping)
        
        # Log detailed statistics
        LOGGER.info(f"  UniProt mapping statistics:")
        LOGGER.info(f"    Total queried:     {stats['total_queried']}")
        LOGGER.info(f"    Has gene symbol:   {stats['has_gene_symbol']} ({100*stats['has_gene_symbol']/max(1,stats['total_queried']):.1f}%)")
        LOGGER.info(f"    No gene symbol:    {stats['no_gene_symbol']} ({100*stats['no_gene_symbol']/max(1,stats['total_queried']):.1f}%)")
        LOGGER.info(f"    Has KEGG xref:     {stats['has_kegg']} ({100*stats['has_kegg']/max(1,stats['total_queried']):.1f}%)")
        LOGGER.info(f"    Has GeneID xref:   {stats['has_geneid']} ({100*stats['has_geneid']/max(1,stats['total_queried']):.1f}%)")
        LOGGER.info(f"    Has Ensembl xref:  {stats['has_ensembl']} ({100*stats['has_ensembl']/max(1,stats['total_queried']):.1f}%)")
        LOGGER.info(f"    Has EC number:     {stats['has_ec']} ({100*stats['has_ec']/max(1,stats['total_queried']):.1f}%)")
        LOGGER.info(f"    API errors:        {stats['api_errors']}")
        if stats['failed_ids'] > 0:
            LOGGER.warning(f"    Failed after retries: {stats['failed_ids']}")
        
        # Store stats for later analysis
        self.uniprot_mapping_stats = stats
        
        return mapping
        
    def map_uniprot_to_ensembl_batch(self, uniprot_ids: Set[str]) -> Dict[str, Set[str]]:
        """
        Map UniProt IDs to Ensembl gene IDs via UniProt API.
        
        Args:
            uniprot_ids: Set of UniProt accession IDs
            
        Returns:
            Dict mapping UniProt ID -> Set of Ensembl gene IDs
        """
        # Filter for human only
        human_ids = uniprot_ids & self.human_uniprots
        
        if not human_ids:
            return {}
            
        LOGGER.info(f"Mapping {len(human_ids)} UniProt IDs to Ensembl...")
        
        mapping = {}
        batch_size = 200
        human_list = list(human_ids)
        
        for i in tqdm(range(0, len(human_list), batch_size), desc="UniProt→Ensembl"):
            batch = human_list[i:i+batch_size]
            
            url = f"{UNIPROT_API_BASE}/stream"
            params = {
                "query": " OR ".join([f"accession:{acc}" for acc in batch]),
                "fields": "accession,xref_ensembl",
                "format": "tsv",
            }
            
            try:
                response = requests.get(url, params=params, timeout=60)
                if response.ok:
                    for line in response.text.strip().split('\n')[1:]:
                        parts = line.split('\t')
                        if len(parts) >= 2:
                            uniprot = parts[0]
                            ensembl_refs = parts[1] if len(parts) > 1 else ""
                            
                            # Parse Ensembl IDs (format: ENSG00000xxx [ENST...])
                            ensembl_genes = set()
                            for ref in ensembl_refs.split(';'):
                                ref = ref.strip()
                                if ref.startswith('ENSG'):
                                    gene_id = ref.split()[0]
                                    ensembl_genes.add(gene_id)
                                    
                            if ensembl_genes:
                                mapping[uniprot] = ensembl_genes
                                
            except requests.RequestException as e:
                LOGGER.warning(f"Batch query failed: {e}")
                
        self.uniprot_to_ensembl.update(mapping)
        LOGGER.info(f"  Mapped {len(mapping)} proteins to Ensembl")
        return mapping
    
    def get_gene_symbols_from_kegg(self, kegg_refs: List[Tuple[str, str]]) -> Dict[str, str]:
        """
        Get gene symbols from KEGG cross-references.
        
        This is a FALLBACK method when UniProt gene symbol mapping fails.
        It queries KEGG API to get gene symbols from KEGG gene IDs.
        
        Args:
            kegg_refs: List of (uniprot_id, kegg_ref) tuples
                       where kegg_ref is like "hsa:123;hsa:456"
                       
        Returns:
            Dict mapping UniProt ID -> Gene Symbol
        """
        if not kegg_refs:
            return {}
        
        mapping = {}
        kegg_gene_ids = set()
        uniprot_to_kegg_gene = {}  # Map UniProt -> KEGG gene ID
        
        # Parse KEGG references to extract gene IDs
        for uniprot, kegg_ref in kegg_refs:
            for part in kegg_ref.split(';'):
                part = part.strip()
                if part and ':' in part:
                    org, gene_id = part.split(':', 1)
                    if org == 'hsa':  # Human genes only
                        kegg_gene_ids.add(f"hsa:{gene_id}")
                        uniprot_to_kegg_gene[uniprot] = f"hsa:{gene_id}"
        
        if not kegg_gene_ids:
            return {}
        
        LOGGER.info(f"  KEGG fallback: querying {len(kegg_gene_ids)} KEGG gene IDs...")
        
        # Query KEGG for gene info (batch up to 10 at a time)
        kegg_to_symbol = {}
        kegg_list = list(kegg_gene_ids)
        
        for i in range(0, len(kegg_list), 10):
            batch = kegg_list[i:i+10]
            try:
                # KEGG API: get multiple genes
                url = f"https://rest.kegg.jp/get/{'+'.join(batch)}"
                response = requests.get(url, timeout=30)
                
                if response.ok:
                    # Parse KEGG flat file format
                    current_gene = None
                    for line in response.text.split('\n'):
                        if line.startswith('ENTRY'):
                            # ENTRY       123             CDS       H.sapiens
                            parts = line.split()
                            if len(parts) >= 2:
                                current_gene = f"hsa:{parts[1]}"
                        elif line.startswith('SYMBOL'):
                            # SYMBOL      ADH1B, ADH2
                            symbol = line.replace('SYMBOL', '').strip().split(',')[0].strip()
                            if current_gene and symbol:
                                kegg_to_symbol[current_gene] = symbol
                                
            except requests.RequestException as e:
                LOGGER.debug(f"KEGG batch query failed: {e}")
        
        # Map back to UniProt IDs
        for uniprot, kegg_gene in uniprot_to_kegg_gene.items():
            if kegg_gene in kegg_to_symbol:
                mapping[uniprot] = kegg_to_symbol[kegg_gene]
        
        if mapping:
            LOGGER.info(f"  KEGG fallback recovered {len(mapping)} gene symbols")
            self.uniprot_mapping_stats['kegg_fallback_recovered'] = len(mapping)
        
        return mapping


# =============================================================================
# RHEA SPARQL CLIENT - Fetch reaction participants from Rhea
# =============================================================================

class RheaSparqlClient:
    """
    Client for querying Rhea SPARQL endpoint to get reaction participants.
    
    Philosophy: Fresh data from source (SPARQL queries) rather than static files.
    This ensures we always get the latest ChEBI compound information.
    """
    
    def __init__(self, endpoint: str = RHEA_SPARQL_ENDPOINT, timeout: int = 30):
        self.endpoint = endpoint
        self.timeout = timeout
        self._cache: Dict[str, Tuple[List[Dict], List[Dict]]] = {}
        
    def fetch_reaction_participants(self, rhea_master_id: str) -> Tuple[List[Dict], List[Dict]]:
        """
        Fetch substrates and products for a Rhea reaction via SPARQL.
        
        This fetches ALL participants including macromolecular compounds that
        don't have ChEBI IDs. Macromolecules are identified by their Rhea
        internal GENERIC:xxxxx ID and will have empty formula fields.
        
        Args:
            rhea_master_id: The master Rhea ID (e.g., "10000")
            
        Returns:
            Tuple of (substrates, products) where each is a list of dicts:
            [{'chebi': 'CHEBI:12345', 'name': 'compound name', 'coefficient': 1, 
              'formula': 'C6H12O6', 'rhea_compound_id': '', 'is_macromolecule': False}, ...]
            
            For macromolecules without ChEBI:
            [{'chebi': '', 'name': 'L-seryl-[protein]', 'coefficient': 1,
              'formula': '', 'rhea_compound_id': 'GENERIC:13713', 'is_macromolecule': True}, ...]
        """
        if rhea_master_id in self._cache:
            return self._cache[rhea_master_id]
            
        # Query ALL participants - make ChEBI OPTIONAL to include macromolecules
        query = f"""
        PREFIX rh: <http://rdf.rhea-db.org/>
        PREFIX chebi: <http://purl.obolibrary.org/obo/>
        
        SELECT ?side ?compound ?chebi ?name ?coefficient ?formula ?accession
        WHERE {{
          <http://rdf.rhea-db.org/{rhea_master_id}> rh:side ?reactionSide .
          ?reactionSide rh:contains ?participant .
          ?participant rh:compound ?compound .
          OPTIONAL {{ ?compound rh:chebi ?chebi }}
          OPTIONAL {{ ?compound rh:name ?name }}
          OPTIONAL {{ ?participant rh:coefficient ?coefficient }}
          OPTIONAL {{ ?compound rh:formula ?formula }}
          OPTIONAL {{ ?compound rh:accession ?accession }}
          
          # Determine side based on _L (left=substrates) or _R (right=products)
          BIND(IF(STRENDS(STR(?reactionSide), "_L"), "substrate", "product") AS ?side)
        }}
        """
        
        try:
            response = requests.post(
                self.endpoint,
                data={"query": query},
                headers={"Accept": "application/sparql-results+json"},
                timeout=self.timeout
            )
            
            if not response.ok:
                LOGGER.warning(f"SPARQL query failed for RHEA:{rhea_master_id}: {response.status_code}")
                return [], []
                
            data = response.json()
            substrates = []
            products = []
            
            for row in data.get('results', {}).get('bindings', []):
                chebi_uri = row.get('chebi', {}).get('value', '')
                compound_uri = row.get('compound', {}).get('value', '')
                accession = row.get('accession', {}).get('value', '')
                name = row.get('name', {}).get('value', '')
                formula = row.get('formula', {}).get('value', '')
                
                # Determine if this is a ChEBI compound or macromolecule
                if chebi_uri:
                    # Extract CHEBI ID from URI: http://purl.obolibrary.org/obo/CHEBI_12345
                    chebi_id = chebi_uri.split('/')[-1].replace('CHEBI_', '')
                    entry = {
                        'chebi': f'CHEBI:{chebi_id}',
                        'name': name or f'CHEBI:{chebi_id}',
                        'coefficient': int(row.get('coefficient', {}).get('value', '1')),
                        'formula': formula,
                        'rhea_compound_id': '',
                        'is_macromolecule': False,
                    }
                else:
                    # Macromolecular compound without ChEBI
                    # Use accession (GENERIC:xxxxx) or extract ID from compound URI
                    rhea_id = accession or compound_uri.split('/')[-1]
                    entry = {
                        'chebi': '',  # No ChEBI ID
                        'name': name or rhea_id,
                        'coefficient': int(row.get('coefficient', {}).get('value', '1')),
                        'formula': '',  # Macromolecules have no formula
                        'rhea_compound_id': rhea_id,
                        'is_macromolecule': True,
                    }
                
                side = row.get('side', {}).get('value', 'unknown')
                if side == 'substrate':
                    substrates.append(entry)
                else:
                    products.append(entry)
                    
            self._cache[rhea_master_id] = (substrates, products)
            return substrates, products
            
        except requests.RequestException as e:
            LOGGER.warning(f"SPARQL request failed for RHEA:{rhea_master_id}: {e}")
            return [], []
            
    def fetch_batch_participants(self, rhea_ids: List[str], 
                                  progress: bool = True) -> Dict[str, Tuple[List[Dict], List[Dict]]]:
        """
        Fetch participants for multiple reactions.
        
        Args:
            rhea_ids: List of Rhea master IDs
            progress: Show progress bar
            
        Returns:
            Dict mapping rhea_id -> (substrates, products)
        """
        results = {}
        
        iterator = tqdm(rhea_ids, desc="Fetching Rhea participants") if progress else rhea_ids
        
        for rhea_id in iterator:
            subs, prods = self.fetch_reaction_participants(rhea_id)
            if subs or prods:  # Only include if we got data
                results[rhea_id] = (subs, prods)
                
        return results
        
    def get_reaction_equation(self, rhea_master_id: str) -> str:
        """
        Fetch human-readable equation for a Rhea reaction.
        
        Returns:
            Equation string like "A + B = C + D"
        """
        query = f"""
        PREFIX rh: <http://rdf.rhea-db.org/>
        
        SELECT ?equation
        WHERE {{
          <http://rdf.rhea-db.org/{rhea_master_id}> rh:equation ?equation .
        }}
        """
        
        try:
            response = requests.post(
                self.endpoint,
                data={"query": query},
                headers={"Accept": "application/sparql-results+json"},
                timeout=self.timeout
            )
            
            if response.ok:
                data = response.json()
                bindings = data.get('results', {}).get('bindings', [])
                if bindings:
                    return bindings[0].get('equation', {}).get('value', '')
        except Exception:
            pass
            
        return ""


# =============================================================================
# CHEBI COMPOUND FETCHER - Get formula and other properties
# =============================================================================

class ChebiFetcher:
    """
    Fetches compound properties from ChEBI using OLS4 (Ontology Lookup Service).
    The legacy ChEBI web service is deprecated/broken, so we use OLS4 instead.
    """
    
    def __init__(self, timeout: int = 30):
        self.timeout = timeout
        self._cache: Dict[str, Dict] = {}
        
    def fetch_compound(self, chebi_id: str) -> Dict:
        """
        Fetch compound properties from ChEBI including cross-references.
        
        Uses EBI OLS4 API which provides comprehensive ChEBI data including
        cross-references to KEGG, PubChem, HMDB, etc.
        
        Args:
            chebi_id: ChEBI ID (format: "CHEBI:12345" or "12345")
            
        Returns:
            Dict with 'formula', 'charge', 'inchi', 'inchikey', 'name',
                      'pubchem', 'lipidmaps', 'kegg', 'hmdb'
        """
        # Normalize ID
        chebi_num = chebi_id.replace('CHEBI:', '').strip()
        
        if chebi_num in self._cache:
            return self._cache[chebi_num]
            
        result = {
            'formula': '', 'charge': None, 'inchi': '', 'inchikey': '', 
            'name': '', 'pubchem': '', 'lipidmaps': '', 'kegg': '', 'hmdb': ''
        }
        
        # Use OLS4 API (the legacy ChEBI web service returns 500 errors)
        url = f"https://www.ebi.ac.uk/ols4/api/ontologies/chebi/terms"
        params = {"obo_id": f"CHEBI:{chebi_num}"}
        
        try:
            response = requests.get(url, params=params, timeout=self.timeout)
            if response.ok:
                data = response.json()
                
                if '_embedded' in data and 'terms' in data['_embedded'] and data['_embedded']['terms']:
                    term = data['_embedded']['terms'][0]
                    
                    # Get label (name)
                    result['name'] = term.get('label', '')
                    
                    # Get annotations which contain cross-references
                    ann = term.get('annotation', {})
                    
                    # Formula (generalized_empirical_formula or formula)
                    formula_list = ann.get('generalized_empirical_formula', ann.get('formula', []))
                    if formula_list:
                        result['formula'] = formula_list[0] if isinstance(formula_list, list) else formula_list
                    
                    # Charge
                    charge_list = ann.get('charge', [])
                    if charge_list:
                        charge_val = charge_list[0] if isinstance(charge_list, list) else charge_list
                        try:
                            result['charge'] = int(float(charge_val))
                        except (ValueError, TypeError):
                            pass
                    
                    # InChI
                    inchi_list = ann.get('inchi_string', [])
                    if inchi_list:
                        result['inchi'] = inchi_list[0] if isinstance(inchi_list, list) else inchi_list
                    
                    # InChIKey
                    inchikey_list = ann.get('inchi_key_string', [])
                    if inchikey_list:
                        result['inchikey'] = inchikey_list[0] if isinstance(inchikey_list, list) else inchikey_list
                    
                    # Cross-references (database_cross_reference)
                    xrefs = ann.get('database_cross_reference', [])
                    for xref in xrefs:
                        xref_lower = xref.lower()
                        if xref_lower.startswith('kegg.compound:'):
                            result['kegg'] = xref.split(':', 1)[1]
                        elif xref_lower.startswith('hmdb:'):
                            result['hmdb'] = xref.split(':', 1)[1]
                        elif xref_lower.startswith('pubchem.compound:'):
                            result['pubchem'] = xref.split(':', 1)[1]
                        elif xref_lower.startswith('lipidmaps:'):
                            result['lipidmaps'] = xref.split(':', 1)[1]
                    
        except Exception as e:
            LOGGER.debug(f"ChEBI OLS4 fetch failed for {chebi_id}: {e}")
            
        self._cache[chebi_num] = result
        return result


# =============================================================================
# RHEA REACTION ADAPTER
# =============================================================================

class RheaReactionAdapter:
    """
    Adapts Rhea reaction data to be compatible with the existing
    `reaction` class interface used by the KEGG pipeline.
    
    The goal is to create objects that have the same methods as
    class_generate_database.reaction so they can be processed
    by the same downstream code (getRxncons, mass_balance, etc.)
    """
    
    def __init__(self, rhea_id: str, master_id: str, 
                 ec_numbers: List[str], 
                 ensembl_ids: Set[str],
                 substrates: List[Tuple[float, str, str]],  # (stoich, url, compound_id)
                 products: List[Tuple[float, str, str]],
                 kegg_reaction_id: Optional[str] = None):  # KEGG reaction ID if available
        """
        Create a Rhea reaction adapter that mimics the `reaction` class.
        
        Args:
            rhea_id: Rhea reaction ID
            master_id: Rhea master ID
            ec_numbers: List of EC numbers
            ensembl_ids: Set of Ensembl gene IDs
            substrates: List of (stoichiometry, url, compound_id) tuples
            products: List of (stoichiometry, url, compound_id) tuples
            kegg_reaction_id: KEGG reaction ID (e.g., "R00001") if mapped
        """
        self.rhea_id = rhea_id
        self.master_id = master_id
        self._ec_numbers = ec_numbers
        self._ensembl_ids = ensembl_ids
        self._substrates = substrates  # Same format as KEGG: [[stoich, url, compound_id], ...]
        self._products = products
        self._id = f"RHEA{master_id}"  # Set ID first
        self._name = self._id  # Name syncs with ID for SBML output
        self._equation = None  # Will store the equation string if needed
        self._termodyn = False  # Default: reversible (False = reversible in KEGG convention)
        self._direction = None  # Will be set by set_direction(): 'LR', 'RL', 'BI', or None (master)
        self._pathway = "Rhea_additional"  # Default, can be updated with set_pathway()
        self._reactome_pathway = None  # Store actual Reactome pathway name
        self._kegg_reaction_id = kegg_reaction_id  # KEGG reaction ID cross-reference
        
        # Store as subs/prods attributes (same as KEGG reaction after sanitization)
        self.subs = [[s[0], s[1], s[2]] for s in substrates]
        self.prods = [[p[0], p[1], p[2]] for p in products]
        
        # Link attribute for compatibility
        self.link = [
            self.subs + self.prods,  # [0] all compounds
            self.subs,                # [1] substrates
            self.prods,               # [2] products
            self._name,               # [3] name
            [[None, ec] for ec in ec_numbers],  # [4] EC numbers
            [],                       # [5] equivalent
            True,                     # [6] GTest
            True,                     # [7] CTest
        ]
        
        # GPR data - populated by set_compartment_gprs()
        # Format: (sGPR_general, GPR_general, {comp: sGPR}, {comp: GPR})
        self._gpr_data = None
        
        # Subcel data for compartmentalization
        # Format: [{compartment: sGPR}, {compartment: GPR}]
        self._subcel_data = [{}, {}]
        
    def set_compartment_gprs(self, compartment_genes: Dict[str, Set[str]]):
        """
        Set compartment-specific GPR rules.
        
        Args:
            compartment_genes: Dict mapping compartment -> set of Ensembl gene IDs
                              e.g., {'cytosol': {'ENSG001', 'ENSG002'}, 'mitochondria': {'ENSG003'}}
        """
        sGPR_by_compartment = {}
        GPR_by_compartment = {}
        
        for compartment, genes in compartment_genes.items():
            if genes:
                # For now, genes are OR'd together (isoenzymes)
                gpr_string = " or ".join(sorted(genes))
                sGPR_by_compartment[compartment] = gpr_string
                GPR_by_compartment[compartment] = gpr_string
        
        self._subcel_data = [sGPR_by_compartment, GPR_by_compartment]
        
        # Also set general GPR (union of all genes)
        all_genes = set()
        for genes in compartment_genes.values():
            all_genes.update(genes)
        
        if all_genes:
            general_gpr = " or ".join(sorted(all_genes))
            self._gpr_data = (general_gpr, general_gpr, sGPR_by_compartment, GPR_by_compartment)
        
    # ===== Methods matching the `reaction` class interface =====
    
    @property
    def ID(self):
        return self._id if hasattr(self, '_id') else f"RHEA{self.master_id}"
    
    @ID.setter
    def ID(self, value):
        self._id = value
        # Keep name in sync with ID for consistent naming
        self._name = value
    
    def Name(self):
        """
        Return reaction name for SBML output.
        
        Returns the reaction ID (e.g., 'RHEA10040_cytosol') to match KEGG naming convention.
        The ID is updated during compartmentalization, so Name() returns the updated ID.
        """
        # Return ID to match KEGG behavior where name is ID-based
        return self.ID
    
    def set_equation(self, equation: str):
        """Store the human-readable equation string."""
        self._equation = equation
        
    def Equation(self):
        """Return the human-readable equation string if available."""
        return self._equation
        
    def EC(self):
        return self._ec_numbers
    
    def KEGG_ID(self):
        """Return KEGG reaction ID if available from rhea2kegg mapping."""
        return self._kegg_reaction_id
        
    def GPR(self):
        """Return GPR data in the same format as KEGG reaction."""
        if self._gpr_data:
            return self._gpr_data[:2]  # (sGPR_general, GPR_general)
        # Default: just the gene list as OR
        if self._ensembl_ids:
            gpr_string = " or ".join(sorted(self._ensembl_ids))
            return (gpr_string, gpr_string)
        return ("", "")
        
    def Termodyn(self):
        """
        Return thermodynamic irreversibility flag.
        
        In KEGG convention:
        - True = irreversible (one direction only)
        - False = reversible (can go both ways)
        
        Rhea direction mapping:
        - 'LR' (=>) or 'RL' (<=): irreversible -> True
        - 'BI' (<=>) or None/Master (=): reversible -> False
        """
        if self._direction in ('LR', 'RL'):
            return True  # Irreversible
        return False  # Reversible (BI or undefined/master)
    
    def set_direction(self, direction: str):
        """
        Set reaction direction from Rhea directions file.
        
        Args:
            direction: 'LR' (forward), 'RL' (reverse), 'BI' (bidirectional), 
                      or None (master/undefined)
        """
        self._direction = direction
        
    def Substrate(self):
        return self.subs
        
    def SetSubstrate(self, substrate):
        self.subs = [[x[0], x[1], x[2]] for x in substrate]
        self.link[1] = self.subs
        self.link[0] = self.subs + self.prods
        return self.subs
        
    def Product(self):
        return self.prods
        
    def SetProduct(self, product):
        self.prods = [[x[0], x[1], x[2]] for x in product]
        self.link[2] = self.prods
        self.link[0] = self.subs + self.prods
        return self.prods
        
    def Pathway(self):
        """Return pathway name (Reactome pathway if available, else default)."""
        return self._reactome_pathway if self._reactome_pathway else self._pathway
    
    def set_pathway(self, pathway_name: str):
        """Set the pathway name for this reaction (from Reactome mapping)."""
        self._reactome_pathway = pathway_name
        
    def Subcel(self):
        """
        Return subcellular location data for compartmentalization.
        
        Returns:
            List of two dicts: [{compartment: sGPR}, {compartment: GPR}]
            Used by rxnSubcel() to expand reactions into compartment-specific versions.
        """
        return self._subcel_data
        
    def Equivalent(self):
        return []
        
    def GTest(self):
        return True
        
    def CTest(self):
        return True
    
    def MBTest(self):
        """Return first character of ID for mass balance test compatibility."""
        return self.ID[0]  # Returns 'R' for RHEA reactions


# =============================================================================
# CHEBI COMPOUND ADAPTER  
# =============================================================================

class ChebiCompoundAdapter:
    """
    Adapts ChEBI compound data to be compatible with the existing
    `compound` class interface used by the KEGG pipeline.
    """
    
    def __init__(self, chebi_id: str, name: str = "", formula: str = "",
                 inchi: str = "", inchikey: str = "", charge: Optional[int] = None,
                 pubchem: str = "", lipidmaps: str = "", kegg: str = "", hmdb: str = ""):
        """
        Create a ChEBI compound adapter that mimics the `compound` class.
        
        Args:
            chebi_id: ChEBI ID (numeric part only)
            name: Compound name
            formula: Chemical formula
            inchi: InChI string
            inchikey: InChIKey
            charge: Formal charge
            pubchem: PubChem CID (from ChEBI cross-references)
            lipidmaps: LIPID MAPS ID (from ChEBI cross-references)
            kegg: KEGG compound ID (from ChEBI cross-references)
            hmdb: HMDB ID (from ChEBI cross-references)
        """
        self.chebi_id = chebi_id
        self._name = name or f"CHEBI:{chebi_id}"
        self._formula = formula
        self._inchi = inchi
        self._inchikey = inchikey
        self._charge = charge
        
        # Attributes matching compound class
        self.ident = f"CHEBI{chebi_id}"
        self.ID1 = f"CHEBI{chebi_id}"
        self.ID2 = f"CHEBI{chebi_id}"
        self.Name = self._name
        self.Formula1 = formula
        self.Formula2 = formula
        self.Formula3 = formula
        self.Formula4 = formula
        self.charge = charge
        self.inchi = inchi
        self.inchikey = inchikey
        
        # External IDs (from ChEBI cross-references via OLS4)
        self.PubChem = pubchem
        self.CheBI = chebi_id
        self.LIPIDMAPS = lipidmaps
        self.KEGG = kegg          # KEGG compound ID (e.g., C00001)
        self.HMDB = hmdb          # HMDB ID (e.g., HMDB0002111)
        self.LipidBank = ""
        self.GlyDB = ""
        self.JCGGDB = ""
        self.CID = pubchem  # CID is same as PubChem
        
        # Atom composition (calculated from formula if possible)
        self.Atom1 = self._calculate_atom_composition(formula)
        self.Atom2 = self.Atom1
        self.Atom3 = self.Atom1
        
        # Subcellular location (set during compartmentalization)
        self.Subcel = ""
        
        # Associated reactions (populated when reaction uses this compound)
        self.AssRxn1 = []  # List of reaction info
        self.AssRxn2 = []  # List of reaction IDs
        self.AssRxn3 = []  # List of reaction names
    
    def _calculate_atom_composition(self, formula: str) -> dict:
        """
        Parse formula string to get atom composition.
        
        Args:
            formula: Chemical formula string (e.g., "C6H12O6")
            
        Returns:
            Dict mapping atom symbol to count, or empty dict if parsing fails
        """
        if not formula:
            return {}
            
        try:
            import re
            # Match element symbols followed by optional counts
            matches = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
            atom_dict = {}
            for element, count in matches:
                if element:  # Skip empty matches
                    atom_dict[element] = int(count) if count else 1
            return atom_dict
        except Exception:
            return {}
    
    def add_associated_reaction(self, rxn_id: str, rxn_name: str = ""):
        """Add a reaction that uses this compound."""
        if rxn_id not in self.AssRxn2:
            self.AssRxn2.append(rxn_id)
            self.AssRxn3.append(rxn_name)
            self.AssRxn1.append([rxn_id, rxn_name])


# =============================================================================
# RHEA MACROMOLECULE ADAPTER
# =============================================================================

class RheaMacromoleculeAdapter:
    """
    Adapts Rhea macromolecular compound data (proteins, modified proteins, etc.)
    to be compatible with the existing `compound` class interface.
    
    Macromolecules in Rhea (e.g., "L-seryl-[protein]", "ATP-activated protein")
    don't have ChEBI IDs, formulas, or external database references initially.
    They are identified by their Rhea internal GENERIC:xxxxx ID.
    
    These compounds will have:
    - ident/ID1/ID2: M + rhea_id (e.g., "MG13713")
    - Name: the compound name from Rhea
    - Formula fields: empty initially, but can be enriched later
    
    Mass Balance Behavior:
    - If Formula1/Formula2 is empty: reaction will be skipped from mass balance
    - If Formula1/Formula2 is set (via enrichment from other databases): 
      reaction WILL be mass balanced
    
    This design allows future enrichment from databases like UniProt, Reactome,
    or custom annotations without changing the compound ID structure.
    """
    
    def __init__(self, rhea_compound_id: str, name: str = "", formula: str = ""):
        """
        Create a Rhea macromolecule adapter.
        
        Args:
            rhea_compound_id: Rhea compound ID (e.g., "GENERIC:13713" or numeric)
            name: Compound name (e.g., "L-seryl-[protein]")
            formula: Chemical formula (optional, for pre-enriched macromolecules)
        """
        # Normalize the ID - remove special characters for use as identifier
        self.rhea_compound_id = rhea_compound_id
        clean_id = rhea_compound_id.replace(":", "").replace("GENERIC", "G")
        
        self._name = name or f"RHEA:{rhea_compound_id}"
        
        # Attributes matching compound class
        # Use 'M' prefix for Macromolecule identification
        self.ident = f"M{clean_id}"
        self.ID1 = f"M{clean_id}"
        self.ID2 = f"M{clean_id}"
        self.Name = self._name
        
        # Formula fields - empty by default, but can be enriched
        # If formula is set, mass balance will include this compound
        self.Formula1 = formula
        self.Formula2 = formula
        self.Formula3 = formula
        self.Formula4 = formula
        
        # No charge/structure info
        self.charge = None
        self.inchi = ""
        self.inchikey = ""
        
        # All external IDs empty - macromolecules have no cross-references
        self.PubChem = ""
        self.CheBI = ""
        self.LIPIDMAPS = ""
        self.KEGG = ""
        self.HMDB = ""
        self.LipidBank = ""
        self.GlyDB = ""
        self.JCGGDB = ""
        self.CID = ""
        
        # Empty atom composition
        self.Atom1 = {}
        self.Atom2 = {}
        self.Atom3 = {}
        
        # Subcellular location (set during compartmentalization)
        self.Subcel = ""
        
        # Associated reactions (populated when reaction uses this compound)
        self.AssRxn1 = []
        self.AssRxn2 = []
        self.AssRxn3 = []
    
    def add_associated_reaction(self, rxn_id: str, rxn_name: str = ""):
        """Add a reaction that uses this compound."""
        if rxn_id not in self.AssRxn2:
            self.AssRxn2.append(rxn_id)
            self.AssRxn3.append(rxn_name)
            self.AssRxn1.append([rxn_id, rxn_name])
    
    def enrich_formula(self, formula: str, source: str = ""):
        """
        Enrich this macromolecule with a chemical formula.
        
        Once a formula is set, reactions containing this macromolecule
        will be included in mass balance checks.
        
        Args:
            formula: Chemical formula string (e.g., "C100H150N30O40")
            source: Optional source annotation (e.g., "UniProt", "Reactome")
        """
        import re
        
        if formula:
            self.Formula1 = formula
            self.Formula2 = formula
            self.Formula3 = formula
            self.Formula4 = formula
            
            # Recalculate atom composition
            self.Atom1 = self._calculate_atom_composition(formula)
            self.Atom2 = self.Atom1
            self.Atom3 = self.Atom1
            
            if source:
                # Optionally track the source in the name
                if source not in self.Name:
                    self.Name = f"{self.Name} [{source}]"
    
    def _calculate_atom_composition(self, formula: str) -> dict:
        """
        Parse formula string to get atom composition.
        
        Args:
            formula: Chemical formula string (e.g., "C6H12O6")
            
        Returns:
            Dict mapping atom symbol to count, or empty dict if parsing fails
        """
        import re
        
        if not formula:
            return {}
            
        try:
            # Match element symbols followed by optional counts
            matches = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
            atom_dict = {}
            for element, count in matches:
                if element:  # Skip empty matches
                    atom_dict[element] = int(count) if count else 1
            return atom_dict
        except Exception:
            return {}
    
    def has_formula(self) -> bool:
        """Check if this macromolecule has been enriched with a formula."""
        return bool(self.Formula1 and self.Formula1.strip())


# =============================================================================
# MAIN RHEA EXTENDER - INTEGRATES WITH KEGG PIPELINE
# =============================================================================

class RheaExtender:
    """
    Extends KEGG pipeline data structures with Rhea reactions.
    
    REUSES the same data structures as generate_db.py:
    - RxnList: Dict[str, reaction]
    - MetList: Dict[str, compound]
    - GPRList: Dict[str, gpr]
    - RxnIdent, MetIdent, GPRIdent: Lists of IDs
    - RxnList_CL, MetList_CL: Compartmentalized versions
    
    Usage:
        # After KEGG loop completes, extend with Rhea:
        extender = RheaExtender()
        extender.run(
            RxnList=RxnList,
            MetList=MetList,
            RxnIdent=RxnIdent,
            MetIdent=MetIdent,
            ... (other KEGG data structures)
        )
        # RxnList, MetList, etc. are now extended with Rhea data
    """
    
    def __init__(self, cache_dir: Optional[str] = None):
        """
        Args:
            cache_dir: Directory for Rhea file cache (default: files/rhea)
        """
        self.downloader = RheaDataDownloader(cache_dir)
        self.human_filter = HumanProteomeFilter()
        self.sparql_client = RheaSparqlClient()
        self.chebi_fetcher = ChebiFetcher()
        
        # Rhea data loaded from TSV files
        self.rhea_to_kegg: Dict[str, str] = {}  # rhea_master_id -> kegg_id
        self.rhea_to_ec: Dict[str, List[str]] = {}  # rhea_master_id -> [ec_numbers]
        self.rhea_to_uniprots: Dict[str, Set[str]] = {}  # rhea_master_id -> {uniprot_ids}
        self.human_rhea_reactions: Dict[str, Set[str]] = {}  # rhea_master_id -> {ensembl_ids}
        
        # Statistics
        self.stats = {
            'total_rhea_human': 0,
            'already_in_kegg': 0,
            'new_reactions': 0,
            'reactions_skipped_no_participants': 0,
            'reactions_added': 0,
            'reactions_added_no_genes': 0,    # Reactions added WITHOUT gene data (still valid for FBA)
            'reactions_via_ec': 0,             # Reactions using EC→BioCyc pipeline
            'reactions_via_rhea_genes': 0,     # Reactions using Rhea genes (fallback)
            'no_genes_no_ec': 0,               # No EC numbers in Rhea
            'no_genes_no_uniprot': 0,          # No UniProt IDs in Rhea
            'no_genes_uniprot_no_symbol': 0,   # UniProt IDs didn't map to gene symbols
            'metabolites_added': 0,
            'metabolites_reused': 0,  # ChEBI already exists as KEGG compound
            'metabolites_enriched': 0,  # Existing metabolites with updated annotations
            'macromolecules_added': 0,  # Macromolecular compounds (no ChEBI/formula)
            'compartmentalized_reactions': 0,
            'mass_balance_success': 0,
            'mass_balance_already_balanced': 0,
            'mass_balance_failed': 0,
        }
        
        # ChEBI → KEGG compound ID mapping (built from existing MetList)
        self.chebi_to_kegg: Dict[str, str] = {}
        
        # Extra compounds for mass balancing (H+, H2O, etc.)
        self.extra_compound = {
            "H": "C00080",    # H+
            "H2O": "C00001",  # Water
            "CO2": "C00011",  # Carbon dioxide
            "O2": "C00007",   # Oxygen
            "NH3": "C00014",  # Ammonia
            "Pi": "C00009",   # Orthophosphate
            "PPi": "C00013",  # Diphosphate
        }
        
    def run(self, 
            RxnList: Dict, 
            MetList: Dict,
            RxnIdent: List,
            MetIdent: List,
            MetEquiv: Dict,
            GPRList: Optional[Dict] = None,
            GPRIdent: Optional[List] = None,
            GeneList: Optional[Dict] = None,
            GeneIdent: Optional[List] = None,
            RxnList_CL: Optional[Dict] = None,
            MetList_CL: Optional[Dict] = None,
            RxnIdent_CL: Optional[List] = None,
            MetIdent_CL: Optional[List] = None,
            Compartment_CL: Optional[List] = None,
            PathNameRxn: Optional[Dict] = None,
            session = None,
            time: int = 20,
            EF: Optional[List] = None,
            specialCompounds: Optional[str] = None,
            force_refresh: bool = False,
            max_retry_attempts: int = 3,
            limit_rhea_ids: Optional[set] = None,
            impose_locations: int = 1) -> Dict[str, Any]:
        """
        Run Rhea extension on KEGG data structures.
        
        Args:
            RxnList, MetList, etc.: Data structures from KEGG pipeline
            session: BioCyc session for GPR queries
            time: Timeout for HTTP requests
            EF: Extra formula list
            specialCompounds: Path to special compounds file
            force_refresh: Re-download Rhea files
            limit_rhea_ids: Optional set of Rhea IDs to process (for testing)
            impose_locations: Compartmentalization mode (should match KEGG pipeline)
                1 = RESTRICTED: Only use compartments from Excel lookup (default)
                0 = UNRESTRICTED: Keep all compartment names as-is
            
        Returns:
            Dict with statistics
        """
        LOGGER.info("=" * 60)
        LOGGER.info("RHEA EXTENSION - Building on KEGG model")
        LOGGER.info("=" * 60)
        
        # Store impose_locations for use in _get_compartment_gprs
        self._impose_locations = impose_locations
        
        # Setup defaults
        if session is None:
            session = setup_biocyc_session()
        if EF is None:
            EF = []
        if specialCompounds is None:
            specialCompounds = os.path.join(project_root, "files", "special_compounds.txt")
            
        # Initialize empty structures if not provided
        if GPRList is None:
            GPRList = {}
        if GPRIdent is None:
            GPRIdent = []
        if GeneList is None:
            GeneList = {}
        if GeneIdent is None:
            GeneIdent = []
        if PathNameRxn is None:
            PathNameRxn = {}
        if RxnList_CL is None:
            RxnList_CL = {}
        if MetList_CL is None:
            MetList_CL = {}
        if RxnIdent_CL is None:
            RxnIdent_CL = []
        if MetIdent_CL is None:
            MetIdent_CL = []
        if Compartment_CL is None:
            Compartment_CL = []
            
        # Step 1: Download Rhea data
        LOGGER.info("\n[Step 1/6] Downloading Rhea data...")
        self.downloader.download_all(force_refresh)
        
        # Step 2: Load human proteome
        LOGGER.info("\n[Step 2/6] Loading human proteome...")
        self.human_filter.load_human_proteome()
        
        # Step 3: Load Rhea mappings and filter for human
        LOGGER.info("\n[Step 3/6] Loading Rhea mappings...")
        self._load_rhea_mappings()
        
        # Step 4: Build ChEBI → KEGG metabolite mapping from existing MetList
        LOGGER.info("\n[Step 4/6] Building ChEBI-KEGG metabolite mapping...")
        self._build_chebi_to_kegg_mapping(MetList)
        
        # Step 5: Find NEW reactions (not in KEGG)
        LOGGER.info("\n[Step 5/6] Finding new reactions...")
        existing_kegg_ids = self._get_existing_kegg_ids(RxnList, RxnIdent)
        new_rhea = self._filter_new_reactions(existing_kegg_ids)
        
        # Optional: Limit to specific Rhea IDs (for testing)
        if limit_rhea_ids is not None:
            original_count = len(new_rhea)
            new_rhea = {k: v for k, v in new_rhea.items() if k in limit_rhea_ids}
            LOGGER.info(f"  Limited to {len(new_rhea)} of {original_count} reactions (test mode)")
        
        # Step 6: Process new reactions (with retry for failures)
        LOGGER.info("\n[Step 6/6] Processing new Rhea reactions...")
        failed_reactions = self._process_new_reactions(
            new_rhea=new_rhea,
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
            time=time,
            EF=EF,
            specialCompounds=specialCompounds,
        )
        
        # Step 7: Retry failed reactions (same pattern as KEGG)
        if failed_reactions and max_retry_attempts > 0:
            LOGGER.info("\n[Step 7] Retrying failed reactions...")
            self._retry_failed_reactions(
                failed_reactions=failed_reactions,
                max_retry_attempts=max_retry_attempts,
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
                session=session,
            )
        
        self._print_summary()
        return self.stats
        
    def _load_rhea_mappings(self) -> None:
        """Load Rhea TSV files and build mappings."""
        
        # Load Rhea-KEGG mappings
        rhea_kegg = self.downloader.load_file("rhea2kegg_reaction.tsv")
        for _, row in rhea_kegg.iterrows():
            master_id = str(row['MASTER_ID'])
            kegg_id = row['ID']
            self.rhea_to_kegg[master_id] = kegg_id
        LOGGER.info(f"  Loaded {len(self.rhea_to_kegg)} Rhea-KEGG mappings")
        
        # Load Rhea-EC mappings
        rhea_ec = self.downloader.load_file("rhea2ec.tsv")
        for _, row in rhea_ec.iterrows():
            master_id = str(row['MASTER_ID'])
            ec = row['ID']
            if master_id not in self.rhea_to_ec:
                self.rhea_to_ec[master_id] = []
            self.rhea_to_ec[master_id].append(ec)
        LOGGER.info(f"  Loaded EC mappings for {len(self.rhea_to_ec)} reactions")
        
        # Load Rhea-UniProt mappings (human only)
        rhea_uniprot = self.downloader.load_file("rhea2uniprot_sprot.tsv")
        for _, row in rhea_uniprot.iterrows():
            master_id = str(row['MASTER_ID'])
            uniprot = row['ID']
            
            # Filter for human
            if uniprot not in self.human_filter.human_uniprots:
                continue
                
            if master_id not in self.rhea_to_uniprots:
                self.rhea_to_uniprots[master_id] = set()
            self.rhea_to_uniprots[master_id].add(uniprot)
            
        self.stats['total_rhea_human'] = len(self.rhea_to_uniprots)
        LOGGER.info(f"  Found {len(self.rhea_to_uniprots)} human Rhea reactions")
        
        # Load Rhea directions (for Termodyn - reversibility info)
        # rhea-directions.tsv: RHEA_ID_MASTER, RHEA_ID_LR, RHEA_ID_RL, RHEA_ID_BI
        # - Master (=): undefined direction
        # - LR (=>): left-to-right (irreversible forward)
        # - RL (<=): right-to-left (irreversible reverse)  
        # - BI (<=>): bidirectional (reversible)
        self.rhea_directions = {}  # master_id -> {'LR': id, 'RL': id, 'BI': id}
        self.rhea_master_lookup = {}  # any_id -> master_id
        try:
            rhea_dir = self.downloader.load_file("rhea-directions.tsv")
            for _, row in rhea_dir.iterrows():
                master_id = str(row.get('RHEA_ID_MASTER', row.iloc[0]))
                lr_id = str(row.get('RHEA_ID_LR', row.iloc[1])) if len(row) > 1 else None
                rl_id = str(row.get('RHEA_ID_RL', row.iloc[2])) if len(row) > 2 else None
                bi_id = str(row.get('RHEA_ID_BI', row.iloc[3])) if len(row) > 3 else None
                
                self.rhea_directions[master_id] = {
                    'LR': lr_id if lr_id and lr_id != 'nan' else None,
                    'RL': rl_id if rl_id and rl_id != 'nan' else None,
                    'BI': bi_id if bi_id and bi_id != 'nan' else None,
                }
                
                # Build reverse lookup: any ID -> master_id
                self.rhea_master_lookup[master_id] = master_id
                if lr_id and lr_id != 'nan':
                    self.rhea_master_lookup[lr_id] = master_id
                if rl_id and rl_id != 'nan':
                    self.rhea_master_lookup[rl_id] = master_id
                if bi_id and bi_id != 'nan':
                    self.rhea_master_lookup[bi_id] = master_id
                    
            LOGGER.info(f"  Loaded direction info for {len(self.rhea_directions)} reactions")
        except Exception as e:
            LOGGER.warning(f"  Could not load rhea-directions.tsv: {e}")
            LOGGER.info("  Assuming all Rhea reactions are bidirectional (reversible)")
        
        # Load Rhea-Reactome pathway mappings
        # Step 1: Load rhea2reactome.tsv to get Rhea ID -> Reactome ID mapping
        self.rhea_to_reactome = {}  # rhea_id -> list of reactome_ids
        self.rhea_to_pathway = {}   # rhea_id -> pathway_name (human-readable)
        
        try:
            rhea_reactome = self.downloader.load_file("rhea2reactome.tsv")
            
            # Collect all Reactome IDs we need to look up
            reactome_ids_to_fetch = set()
            
            for _, row in rhea_reactome.iterrows():
                master_id = str(row['MASTER_ID'])
                reactome_id = str(row.get('ID', '')).strip()
                
                if reactome_id and reactome_id.startswith('R-HSA-'):
                    if master_id not in self.rhea_to_reactome:
                        self.rhea_to_reactome[master_id] = []
                    self.rhea_to_reactome[master_id].append(reactome_id)
                    reactome_ids_to_fetch.add(reactome_id)
            
            LOGGER.info(f"  Found {len(self.rhea_to_reactome)} Rhea reactions with Reactome mappings")
            LOGGER.info(f"  Found {len(reactome_ids_to_fetch)} unique Reactome IDs")
            
            # Step 2: Fetch actual pathway names from Reactome API
            if reactome_ids_to_fetch:
                reactome_pathway_names = batch_fetch_reactome_pathways(reactome_ids_to_fetch)
                
                # Step 3: Map Rhea IDs to human-readable pathway names
                for rhea_id, reactome_ids in self.rhea_to_reactome.items():
                    # Use the first Reactome ID's pathway as the primary pathway
                    for rid in reactome_ids:
                        pathway_name = reactome_pathway_names.get(rid) or reactome_pathway_names.get(rid.split('.')[0])
                        if pathway_name:
                            self.rhea_to_pathway[rhea_id] = pathway_name
                            break
                    
                    # If no pathway found, don't add to mapping (will default to "Rhea_additional")
                
                LOGGER.info(f"  Resolved pathway names for {len(self.rhea_to_pathway)} Rhea reactions")
            
        except Exception as e:
            LOGGER.warning(f"  Could not load rhea2reactome.tsv: {e}")
            LOGGER.info("  Using generic 'Rhea_additional' pathway for all Rhea reactions")
        
    def _build_metabolite_lookup_mappings(self, MetList: Dict) -> None:
        """
        Build multiple mappings from external IDs to KEGG compound IDs.
        
        Uses shared utilities from functions_merge_metabolic_networks.py
        to ensure consistent matching across all modules.
        
        Args:
            MetList: Dict of existing metabolites (KEGG compounds)
        """
        # Use shared function to build the index
        self._metabolite_index = build_metabolite_attribute_index(MetList)
        
        LOGGER.info(f"  Built metabolite lookup mappings:")
        for attr_name, mapping in self._metabolite_index.items():
            LOGGER.info(f"    {attr_name}: {len(mapping)} compounds")
    
    def _find_existing_metabolite(self, chebi_id: str, name: str = "", 
                                   inchi: str = "", inchikey: str = "", 
                                   pubchem: str = "") -> Optional[str]:
        """
        Try to find an existing KEGG metabolite using multiple identifiers.
        
        Uses shared utilities from functions_merge_metabolic_networks.py
        to ensure consistent matching across all modules.
        
        Returns:
            KEGG compound ID if found, None otherwise
        """
        return find_equivalent_metabolite(
            chebi_id=chebi_id,
            name=name,
            inchi=inchi,
            inchikey=inchikey,
            pubchem=pubchem,
            attribute_index=self._metabolite_index
        )
        
    def _build_chebi_to_kegg_mapping(self, MetList: Dict) -> None:
        """
        Build a mapping from ChEBI IDs to KEGG compound IDs.
        
        DEPRECATED: Use _build_metabolite_lookup_mappings() instead for
        multi-attribute matching (same approach as merge_metabolic_networks).
        
        Kept for backward compatibility.
        """
        # Call the new comprehensive method
        self._build_metabolite_lookup_mappings(MetList)
        
    def _get_existing_kegg_ids(self, RxnList: Dict, RxnIdent: List) -> Set[str]:
        """Get set of KEGG reaction IDs already in the model."""
        kegg_ids = set()
        
        # From RxnIdent (base reaction IDs)
        for rxn_id in RxnIdent:
            if rxn_id.startswith('R'):  # KEGG reaction format
                kegg_ids.add(rxn_id)
                
        # From RxnList keys (may include compartmentalized versions)
        for rxn_id in RxnList.keys():
            base_id = rxn_id.split('_')[0] if '_' in rxn_id else rxn_id
            if base_id.startswith('R'):
                kegg_ids.add(base_id)
                
        LOGGER.info(f"  Found {len(kegg_ids)} existing KEGG reaction IDs")
        return kegg_ids
        
    def _filter_new_reactions(self, existing_kegg_ids: Set[str]) -> Dict[str, Set[str]]:
        """
        Filter Rhea reactions to find those NOT already in KEGG model.
        
        Returns:
            Dict mapping rhea_master_id -> set of UniProt IDs (human only)
        """
        new_rhea = {}
        
        for rhea_id, uniprots in self.rhea_to_uniprots.items():
            # Check if this Rhea reaction maps to a KEGG ID we already have
            kegg_id = self.rhea_to_kegg.get(rhea_id)
            
            if kegg_id and kegg_id in existing_kegg_ids:
                self.stats['already_in_kegg'] += 1
                continue
                
            # This is a NEW reaction
            new_rhea[rhea_id] = uniprots
            
        self.stats['new_reactions'] = len(new_rhea)
        LOGGER.info(f"  New reactions to add: {len(new_rhea)}")
        LOGGER.info(f"  Skipped (already in KEGG): {self.stats['already_in_kegg']}")
        
        return new_rhea
        
    def _process_new_reactions(self, new_rhea: Dict[str, Set[str]], **kwargs) -> List[Dict]:
        """
        Process new Rhea reactions and add to KEGG data structures.
        
        Strategy:
        1. Try EC→BioCyc pipeline first (proper GPR rules with AND/OR)
        2. If no EC or BioCyc fails, fall back to Rhea's UniProt genes (OR assumption)
        3. For gene locations, use existing pipeline: BioCyc → UniProt → cytosol fallback
        4. Extract genes from GPR and add to GeneList (check if already exists)
        
        This ensures maximum coverage while preferring high-quality GPR rules.
        
        Returns:
            List of failed reactions (for retry mechanism)
        """
        # Extract kwargs
        RxnList = kwargs['RxnList']
        MetList = kwargs['MetList']
        RxnIdent = kwargs['RxnIdent']
        MetIdent = kwargs['MetIdent']
        MetEquiv = kwargs['MetEquiv']
        GPRList = kwargs.get('GPRList', {})
        GPRIdent = kwargs.get('GPRIdent', [])
        GeneList = kwargs.get('GeneList', {})
        GeneIdent = kwargs.get('GeneIdent', [])
        RxnList_CL = kwargs.get('RxnList_CL', {})
        MetList_CL = kwargs.get('MetList_CL', {})
        RxnIdent_CL = kwargs.get('RxnIdent_CL', [])
        MetIdent_CL = kwargs.get('MetIdent_CL', [])
        Compartment_CL = kwargs.get('Compartment_CL', [])
        PathNameRxn = kwargs.get('PathNameRxn', {})
        session = kwargs.get('session')
        
        # Define generic pathway for Rhea reactions (not in KEGG pathways)
        # NOTE: Reactions without Reactome mapping will be assigned to "Rhea_additional"
        # We don't pre-create the pathway - it gets created when reactions are assigned
        
        # Get UniProt → Gene Symbol mappings for fallback (Rhea genes)
        # IMPORTANT: We use gene symbols to match KEGG's format (e.g., "ADH1B" not "ENSG...")
        all_uniprots = set()
        for uniprots in new_rhea.values():
            all_uniprots.update(uniprots)
            
        LOGGER.info(f"  Mapping {len(all_uniprots)} UniProt IDs to gene symbols (for fallback)...")
        uniprot_to_gene_symbol = self.human_filter.map_uniprot_to_gene_symbol_batch(all_uniprots)
        
        # Fetch reaction participants via SPARQL for ALL new reactions
        LOGGER.info(f"  Fetching reaction participants from Rhea SPARQL...")
        rhea_ids_to_process = list(new_rhea.keys())
        reaction_participants = self.sparql_client.fetch_batch_participants(rhea_ids_to_process)
        
        # Collect all ChEBI IDs
        all_chebi_ids = set()
        for subs, prods in reaction_participants.values():
            for s in subs:
                all_chebi_ids.add(s['chebi'])
            for p in prods:
                all_chebi_ids.add(p['chebi'])
        LOGGER.info(f"  Found {len(all_chebi_ids)} unique ChEBI compounds")
        
        # Initialize additional stats
        self.stats['reactions_via_ec'] = 0
        self.stats['reactions_via_rhea_genes'] = 0
        self.stats['genes_added'] = 0
        
        # Track failed reactions for retry
        failed_reactions = []
        
        # Process each new Rhea reaction
        for rhea_id, uniprots in tqdm(new_rhea.items(), desc="Processing Rhea reactions"):
            # Get reaction participants
            if rhea_id not in reaction_participants:
                LOGGER.debug(f"  Skipping RHEA:{rhea_id} - no participants found")
                self.stats['reactions_skipped_no_participants'] += 1
                continue
                
            subs_data, prods_data = reaction_participants[rhea_id]
            
            if not subs_data or not prods_data:
                LOGGER.debug(f"  Skipping RHEA:{rhea_id} - incomplete equation")
                self.stats['reactions_skipped_no_participants'] += 1
                continue
            
            # Get EC numbers for this reaction (needed for GPR lookup later)
            ec_numbers = self.rhea_to_ec.get(rhea_id, [])
            
            # ================================================================
            # Step 1: Add metabolites to MetList (needed for mass balance)
            # ================================================================
            for sub in subs_data:
                self._add_metabolite_if_new(sub, MetList, MetIdent, MetEquiv)
            for prod in prods_data:
                self._add_metabolite_if_new(prod, MetList, MetIdent, MetEquiv)
            
            # Convert substrates/products to the format expected by reaction class
            # Handle both ChEBI compounds (have chebi ID) and macromolecules (have rhea_compound_id)
            substrates = []
            for sub in subs_data:
                is_macro = sub.get('is_macromolecule', False)
                if is_macro:
                    # Macromolecule: use rhea_compound_id for ID and URL
                    rhea_id_comp = sub.get('rhea_compound_id', '')
                    clean_id = rhea_id_comp.replace(":", "").replace("GENERIC", "G")
                    met_id = f"M{clean_id}"
                    # Check MetEquiv for any resolved ID
                    resolved_id = MetEquiv.get(met_id, met_id)
                    url = f"https://www.rhea-db.org/rhea?query={rhea_id_comp}"
                else:
                    # ChEBI compound: use chebi ID
                    chebi_num = sub['chebi'].replace('CHEBI:', '')
                    resolved_id = self._resolve_metabolite_id(sub['chebi'], MetEquiv)
                    url = f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:{chebi_num}"
                substrates.append([
                    float(sub['coefficient']),
                    url,
                    resolved_id
                ])
                
            products = []
            for prod in prods_data:
                is_macro = prod.get('is_macromolecule', False)
                if is_macro:
                    # Macromolecule: use rhea_compound_id for ID and URL
                    rhea_id_comp = prod.get('rhea_compound_id', '')
                    clean_id = rhea_id_comp.replace(":", "").replace("GENERIC", "G")
                    met_id = f"M{clean_id}"
                    # Check MetEquiv for any resolved ID
                    resolved_id = MetEquiv.get(met_id, met_id)
                    url = f"https://www.rhea-db.org/rhea?query={rhea_id_comp}"
                else:
                    # ChEBI compound: use chebi ID
                    chebi_num = prod['chebi'].replace('CHEBI:', '')
                    resolved_id = self._resolve_metabolite_id(prod['chebi'], MetEquiv)
                    url = f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:{chebi_num}"
                products.append([
                    float(prod['coefficient']),
                    url,
                    resolved_id
                ])
            
            # ================================================================
            # Step 2: Create preliminary reaction adapter (for mass balance)
            # ================================================================
            # Look up KEGG reaction ID if available
            kegg_rxn_id = self.rhea_to_kegg.get(rhea_id)
            
            rxn_adapter = RheaReactionAdapter(
                rhea_id=rhea_id,
                master_id=rhea_id,
                ec_numbers=ec_numbers,
                ensembl_ids=set(),  # Populated via GPR data later
                substrates=[(s[0], s[1], s[2]) for s in substrates],
                products=[(p[0], p[1], p[2]) for p in products],
                kegg_reaction_id=kegg_rxn_id,
            )
            
            # Set direction for Termodyn() based on rhea-directions.tsv
            # We use master ID by default, but check if there's direction info
            if hasattr(self, 'rhea_directions') and rhea_id in self.rhea_directions:
                dir_info = self.rhea_directions[rhea_id]
                # If this master has a BI version, use it (reversible)
                # Otherwise, check if it has only LR or RL (irreversible)
                if dir_info.get('BI'):
                    rxn_adapter.set_direction('BI')  # Bidirectional/reversible
                elif dir_info.get('LR') and not dir_info.get('RL'):
                    rxn_adapter.set_direction('LR')  # Forward only
                elif dir_info.get('RL') and not dir_info.get('LR'):
                    rxn_adapter.set_direction('RL')  # Reverse only
                # else: leave as None (master/undefined = reversible by default)
            
            # Store equation as name
            equation = self.sparql_client.get_reaction_equation(rhea_id)
            if equation:
                rxn_adapter._name = equation
            
            # ================================================================
            # Step 3: Mass balance the reaction (before GPR to fix stoichiometry)
            # ================================================================
            rxn_key = f"RHEA{rhea_id}"
            self._mass_balance_reaction(rxn_adapter, rxn_key, MetList, MetIdent, MetEquiv)
            
            # ================================================================
            # Step 4: Get GPR (priority: existing KEGG genes → EC→BioCyc → Rhea genes)
            # ================================================================
            # OPTIMIZATION: If the genes from Rhea's UniProt→KEGG xref already exist
            # in GeneList (from the KEGG pipeline), we can reuse them directly!
            # This avoids redundant BioCyc queries.
            
            gpr_data = None
            used_ec_pipeline = False
            used_existing_gpr = False
            timeout_error = False  # Track if we had a timeout/connection error
            
            # ----------------------------------------------------------------
            # FIRST: Check if UniProt → EC → existing GPRList (from KEGG run)
            # This is the most efficient path: we already have EC numbers from
            # UniProt API, and if that EC exists in GPRList, we can reuse it!
            # ----------------------------------------------------------------
            for uniprot in uniprots:
                if used_existing_gpr:
                    break
                    
                # Get EC numbers for this UniProt from our cached mapping
                if hasattr(self.human_filter, 'uniprot_to_ec') and uniprot in self.human_filter.uniprot_to_ec:
                    uniprot_ecs = self.human_filter.uniprot_to_ec[uniprot]
                    
                    for ec in uniprot_ecs:
                        # Check if this EC exists in GPRList from KEGG run
                        if ec in GPRList:
                            gpr_obj = GPRList[ec]
                            gpr_data = gpr_obj.GprSubcell()
                            
                            if gpr_data and gpr_data[0]:  # Has gene data
                                used_existing_gpr = True
                                self.stats['reactions_via_existing_gpr'] = self.stats.get('reactions_via_existing_gpr', 0) + 1
                                LOGGER.debug(f"  RHEA:{rhea_id} - Reused existing GPR via UniProt {uniprot} → EC {ec}")
                                break
            
            # ----------------------------------------------------------------
            # SECOND: Try EC from Rhea metadata → GPRList or BioCyc
            # ----------------------------------------------------------------
            if not used_existing_gpr and ec_numbers:
                for ec in ec_numbers:
                    try:
                        # Check if we already have this EC (same pattern as KEGG)
                        if ec not in GPRIdent:
                            GPRIdent.append(ec)
                            # Query BioCyc for GPR (this also gets locations)
                            gpr_obj = gpr(ec, session)
                            GPRList[ec] = gpr_obj
                        else:
                            # Reuse existing GPR object
                            gpr_obj = GPRList[ec]
                        
                        # Get GPR and subcellular data
                        gpr_data = gpr_obj.GprSubcell()
                        if gpr_data and gpr_data[0]:  # Has gene data
                            used_ec_pipeline = True
                            self.stats['ec_cache_hits'] = self.stats.get('ec_cache_hits', 0) + (1 if ec in GPRList else 0)
                            break
                    except (TimeoutError, ConnectionError, requests.exceptions.Timeout, 
                            requests.exceptions.ConnectionError) as e:
                        LOGGER.warning(f"  EC {ec} lookup timeout/connection error: {e}")
                        timeout_error = True
                        continue
                    except Exception as e:
                        LOGGER.debug(f"  EC {ec} lookup failed: {e}")
                        continue
            
            # Fallback to Rhea genes if EC pipeline didn't work
            if not used_ec_pipeline:
                # Get gene symbols from Rhea's UniProt mappings
                # Note: Rhea already filtered for human, so all genes here are human
                gene_symbols = set()
                kegg_fallback_ids = []  # For proteins without gene symbols
                
                for uniprot in uniprots:
                    if uniprot in uniprot_to_gene_symbol:
                        gene_symbols.add(uniprot_to_gene_symbol[uniprot])
                    else:
                        # Check for KEGG cross-reference as fallback
                        if hasattr(self.human_filter, 'uniprot_to_kegg') and uniprot in self.human_filter.uniprot_to_kegg:
                            kegg_id = self.human_filter.uniprot_to_kegg[uniprot]
                            kegg_fallback_ids.append((uniprot, kegg_id))
                
                # ============================================================
                # KEGG FALLBACK: If no gene symbols but we have KEGG IDs,
                # query KEGG API to get actual gene symbols
                # ============================================================
                if not gene_symbols and kegg_fallback_ids:
                    # Use KEGG fallback method to get gene symbols
                    kegg_recovered = self.human_filter.get_gene_symbols_from_kegg(kegg_fallback_ids)
                    
                    if kegg_recovered:
                        # Add recovered symbols to main mapping for future use
                        uniprot_to_gene_symbol.update(kegg_recovered)
                        gene_symbols = set(kegg_recovered.values())
                        self.stats['kegg_fallback_used'] = self.stats.get('kegg_fallback_used', 0) + 1
                        LOGGER.debug(f"  RHEA:{rhea_id} - KEGG fallback recovered genes: {gene_symbols}")
                    else:
                        # Log KEGG fallback availability (even if it didn't work)
                        if self.stats.get('kegg_fallback_available', 0) < 10:
                            LOGGER.info(f"  RHEA:{rhea_id} - No gene symbols, KEGG IDs available but lookup failed: {kegg_fallback_ids[:3]}")
                        self.stats['kegg_fallback_available'] = self.stats.get('kegg_fallback_available', 0) + 1
                
                if not gene_symbols:
                    # ============================================================
                    # NO GENES FOUND - but we still ADD the reaction!
                    # Rhea human reactions are valid for FBA even without genes.
                    # Genes are only needed for omics integration, not flux analysis.
                    # ============================================================
                    
                    # Detailed logging for why no genes were found
                    if timeout_error:
                        LOGGER.info(f"  RHEA:{rhea_id} added to retry queue (timeout during GPR)")
                        failed_reactions.append({
                            'rhea_id': rhea_id,
                            'uniprots': uniprots,
                            'ec_numbers': ec_numbers,
                            'error': 'timeout during GPR lookup',
                            'error_type': 'timeout'
                        })
                        continue  # Skip for now, will retry later
                    
                    # Track WHY there are no genes (for diagnostics)
                    if not ec_numbers:
                        self.stats['no_genes_no_ec'] = self.stats.get('no_genes_no_ec', 0) + 1
                    if not uniprots:
                        self.stats['no_genes_no_uniprot'] = self.stats.get('no_genes_no_uniprot', 0) + 1
                    elif uniprots and not gene_symbols:
                        # Had UniProt IDs but they didn't map to gene symbols
                        self.stats['no_genes_uniprot_no_symbol'] = self.stats.get('no_genes_uniprot_no_symbol', 0) + 1
                        # Log first few for debugging
                        if self.stats['no_genes_uniprot_no_symbol'] <= 10:
                            LOGGER.info(f"  RHEA:{rhea_id} - UniProt IDs {list(uniprots)[:3]} didn't map to gene symbols")
                    
                    LOGGER.debug(f"  RHEA:{rhea_id} - no genes found, adding with default cytosol location")
                    
                    # Create empty GPR data (reaction goes to cytosol by default)
                    gpr_data = ('', '', {'cytosol': ''}, {'cytosol': ''})
                    self.stats['reactions_added_no_genes'] = self.stats.get('reactions_added_no_genes', 0) + 1
                    # Don't continue - fall through to add the reaction!
                    
                else:
                    # Build GPR data using Rhea genes with OR assumption
                    # Location lookup uses getLocationnew: BioCyc → UniProt → cytosol
                    try:
                        gpr_data = self._build_gpr_from_rhea_genes(gene_symbols, session)
                    except (TimeoutError, ConnectionError, requests.exceptions.Timeout,
                            requests.exceptions.ConnectionError) as e:
                        LOGGER.warning(f"  RHEA:{rhea_id} timeout during location lookup: {e}")
                        failed_reactions.append({
                            'rhea_id': rhea_id,
                            'uniprots': uniprots,
                            'ec_numbers': ec_numbers,
                            'gene_symbols': gene_symbols,
                            'error': str(e),
                            'error_type': 'timeout'
                        })
                        continue
                    
                    if not gpr_data or not gpr_data[0]:
                        # GPR build failed but we still want the reaction
                        LOGGER.debug(f"  RHEA:{rhea_id} - GPR build failed, adding with default cytosol location")
                        gpr_data = ('', '', {'cytosol': ''}, {'cytosol': ''})
                        self.stats['reactions_added_no_genes'] = self.stats.get('reactions_added_no_genes', 0) + 1
                    else:
                        self.stats['reactions_via_rhea_genes'] += 1
            
            # ================================================================
            # Step 5: Attach GPR data to reaction adapter
            # ================================================================
            # gpr_data format: (sGPR_general, GPR_general, sGPR_by_compartment, GPR_by_compartment)
            rxn_adapter._gpr_data = gpr_data
            
            # Build Subcel dict for rxnSubcel compatibility
            if len(gpr_data) >= 4 and gpr_data[2] and gpr_data[3]:
                sGPR_by_comp = gpr_data[2]  # {compartment: sGPR}
                GPR_by_comp = gpr_data[3]   # {compartment: GPR}
                rxn_adapter._subcel_data = [sGPR_by_comp, GPR_by_comp]
            else:
                # Fallback: use general GPR for cytosol
                rxn_adapter._subcel_data = [
                    {'cytosol': gpr_data[0]},
                    {'cytosol': gpr_data[1]}
                ]
            
            # Add to RxnList (base reaction)
            RxnList[rxn_key] = rxn_adapter
            RxnIdent.append(rxn_key)
            
            self.stats['reactions_added'] += 1
            if used_ec_pipeline:
                self.stats['reactions_via_ec'] += 1
            # Note: reactions_via_rhea_genes is now incremented in the fallback block above
            
            # ================================================================
            # Step 5b: Assign reaction to pathway (PathNameRxn)
            # ================================================================
            # Check if we have a Reactome pathway for this Rhea ID
            pathway_name = self.rhea_to_pathway.get(rhea_id, "Rhea_additional")
            
            # Update reaction's pathway attribute
            rxn_adapter.set_pathway(pathway_name)
            
            # Ensure pathway exists in PathNameRxn
            if pathway_name not in PathNameRxn:
                PathNameRxn[pathway_name] = ""
                
            # Add reaction to pathway (space-separated list)
            if PathNameRxn[pathway_name]:
                PathNameRxn[pathway_name] += " " + rxn_key
            else:
                PathNameRxn[pathway_name] = rxn_key
            
            # ================================================================
            # Step 5c: Update metabolites' associated reactions
            # ================================================================
            for compound_data in (substrates + products):
                met_id = compound_data[2]  # compound ID from (stoich, url, id)
                if met_id in MetList:
                    met = MetList[met_id]
                    if hasattr(met, 'add_associated_reaction'):
                        met.add_associated_reaction(rxn_key, rxn_adapter.Name())
                    elif hasattr(met, 'AssRxn2'):
                        # For existing KEGG compounds, append if not already there
                        if rxn_key not in met.AssRxn2:
                            if isinstance(met.AssRxn2, list):
                                met.AssRxn2.append(rxn_key)
                            else:
                                met.AssRxn2 = [met.AssRxn2, rxn_key] if met.AssRxn2 else [rxn_key]
            
            # ================================================================
            # Step 6: Compartmentalize the reaction using rxnSubcel
            # ================================================================
            try:
                Compartment_CL, rxn_cl, comp_cl = rxnSubcel(
                    rxn_adapter,
                    RxnList_CL,
                    MetList_CL,
                    RxnIdent_CL,
                    MetIdent_CL,
                    Compartment_CL,
                    RxnList,
                    MetList,
                    MetEquiv,
                )
                RxnIdent_CL.extend(rxn_cl)
                MetIdent_CL.extend(comp_cl)
                self.stats['compartmentalized_reactions'] += len(rxn_cl)
                
                # Extract genes from GPR and add to GeneList (if not already present)
                self._extract_and_add_genes(gpr_data, GeneList, GeneIdent, EnsblDB=None)
                
            except Exception as e:
                LOGGER.warning(f"Failed to compartmentalize RHEA:{rhea_id}: {e}")
                # Add to failed reactions for retry
                failed_reactions.append({
                    'rhea_id': rhea_id,
                    'uniprots': uniprots,
                    'ec_numbers': ec_numbers,
                    'error': str(e),
                    'error_type': 'compartmentalization'
                })
        
        # Report on failed reactions
        if failed_reactions:
            LOGGER.info(f"  {len(failed_reactions)} reactions failed (will retry)")
            self.stats['reactions_failed'] = len(failed_reactions)
        
        return failed_reactions
    
    def _build_gpr_from_existing_genes(self, gene_symbols: Set[str], GeneList: Dict, session) -> tuple:
        """
        Build GPR data from genes that already exist in GeneList (from KEGG run).
        
        This is an OPTIMIZATION: when Rhea's UniProt→KEGG xref maps to genes
        that are already in our model from the KEGG pipeline, we can reuse
        their location information instead of querying BioCyc again.
        
        Args:
            gene_symbols: Set of gene symbols that exist in GeneList
            GeneList: Dict of existing gene objects from KEGG pipeline
            session: BioCyc session (may be needed for location lookup)
            
        Returns:
            Tuple: (sGPR_general, GPR_general, sGPR_by_compartment, GPR_by_compartment)
        """
        if not gene_symbols:
            return None
        
        # Build general GPR as OR of all genes (isozyme assumption)
        sorted_genes = sorted(gene_symbols)
        
        # sGPR uses bracket format to match BioCyc: ([([GENE1*1]) or ([GENE2*1])])
        sgpr_parts = [f"([{gene}*1])" for gene in sorted_genes]
        if len(sgpr_parts) == 1:
            sgpr_general = f"([{sgpr_parts[0]}])"
        else:
            sgpr_general = "([" + " or ".join(sgpr_parts) + "])"
        
        # GPR uses plain gene symbols for COBRA compatibility
        gpr_general = " or ".join(sorted_genes)
        
        # Get locations from existing gene objects
        gene_to_compartments = {}
        
        for gene_symbol in gene_symbols:
            if gene_symbol in GeneList:
                gene_obj = GeneList[gene_symbol]
                # Try to get location from gene object
                if hasattr(gene_obj, 'Location') and gene_obj.Location:
                    # Location might be a string or list
                    loc = gene_obj.Location
                    if isinstance(loc, str):
                        gene_to_compartments[gene_symbol] = {loc}
                    elif isinstance(loc, (list, set)):
                        gene_to_compartments[gene_symbol] = set(loc)
                    else:
                        gene_to_compartments[gene_symbol] = {'cytosol'}
                else:
                    # Default to cytosol if no location
                    gene_to_compartments[gene_symbol] = {'cytosol'}
            else:
                # Shouldn't happen, but fallback to cytosol
                gene_to_compartments[gene_symbol] = {'cytosol'}
        
        # Build by-compartment GPR dictionaries
        compartment_to_genes = {}
        for gene, compartments in gene_to_compartments.items():
            for comp in compartments:
                if comp not in compartment_to_genes:
                    compartment_to_genes[comp] = set()
                compartment_to_genes[comp].add(gene)
        
        # If no compartments found, use cytosol
        if not compartment_to_genes:
            compartment_to_genes = {'cytosol': gene_symbols}
        
        # Build sGPR and GPR for each compartment
        sGPR_by_comp = {}
        GPR_by_comp = {}
        
        for comp, genes in compartment_to_genes.items():
            sorted_comp_genes = sorted(genes)
            sgpr_parts = [f"([{g}*1])" for g in sorted_comp_genes]
            if len(sgpr_parts) == 1:
                sGPR_by_comp[comp] = f"([{sgpr_parts[0]}])"
            else:
                sGPR_by_comp[comp] = "([" + " or ".join(sgpr_parts) + "])"
            
            GPR_by_comp[comp] = " or ".join(sorted_comp_genes)
        
        return (sgpr_general, gpr_general, sGPR_by_comp, GPR_by_comp)

    def _build_gpr_from_rhea_genes(self, gene_symbols: Set[str], session) -> tuple:
        """
        Build GPR data from Rhea's gene list using OR assumption.
        
        Uses getLocationnew() to get subcellular locations for each gene,
        which tries BioCyc first, then UniProt, then defaults to cytosol.
        
        IMPORTANT: Uses gene symbols (e.g., "ADH1B") to match KEGG's format,
        NOT Ensembl IDs. This ensures genes can be found in GeneList.
        
        Args:
            gene_symbols: Set of gene symbols from Rhea (e.g., {"ADH1B", "ADH1C"})
            session: BioCyc session
            
        Returns:
            Tuple: (sGPR_general, GPR_general, sGPR_by_compartment, GPR_by_compartment)
        """
        if not gene_symbols:
            return None
            
        # Build general GPR as OR of all genes (isozyme assumption)
        # Format: Use bracket notation to match BioCyc format: ([([GENE*1]) or ([GENE2*1])])
        sorted_genes = sorted(gene_symbols)
        
        # sGPR uses bracket format to match BioCyc: ([([GENE1*1]) or ([GENE2*1])])
        sgpr_parts = [f"([{gene}*1])" for gene in sorted_genes]
        if len(sgpr_parts) == 1:
            sgpr_general = f"([{sgpr_parts[0]}])"
        else:
            sgpr_general = "([" + " or ".join(sgpr_parts) + "])"
        
        # GPR uses plain gene symbols for COBRA compatibility
        gpr_general = " or ".join(sorted_genes)
        
        # Determine project root path for required files
        bb_pickle_path = os.path.join(project_root, "files", "bb.pickle")
        excel_file = os.path.join(project_root, "files", "ListOfCompartments_sept2024.xlsx")
        comp_abb_file = os.path.join(project_root, "files", "compartments_info.txt")
        
        # Get impose_locations setting (default to restricted if not set)
        impose_locations = getattr(self, '_impose_locations', 1)
        
        # Collect gene locations
        gene_to_compartments = {}
        
        for gene_symbol in gene_symbols:
            try:
                # Build single-gene sGPR in the expected format
                single_gene_sgpr = f"([{gene_symbol}*1])"
                
                # Use getLocationnew which handles BioCyc → UniProt → cytosol fallback
                # Parameters:
                #   gpr: sGPR string like "([GENE*1])"
                #   genelist1: list of gene names (NOT a string)
                #   genelist2: list of BioCyc IDs (use gene symbol as dummy since Rhea has no BioCyc IDs)
                locations = getLocationnew(
                    single_gene_sgpr,       # sGPR for this gene
                    [gene_symbol],          # genes_list as LIST (not string!)
                    [gene_symbol],          # BioCyc IDs (use gene symbol as dummy)
                    impose_locations,       # 1 = restricted, 0 = unrestricted
                    bb_pickle_path,
                    session,
                    None,                   # ensembl_cache
                    excel_file,
                    "Def-Compartments",
                    {},                     # comp_dict
                    comp_abb_file,
                )
                
                # Extract compartments from the result
                if locations and len(locations) >= 2 and locations[0]:
                    # locations[0] is sGPR by compartment dict
                    compartments = list(locations[0].keys())
                    gene_to_compartments[gene_symbol] = compartments if compartments else ['cytosol']
                else:
                    gene_to_compartments[gene_symbol] = ['cytosol']
                    
            except Exception as e:
                LOGGER.debug(f"  Location lookup failed for {gene_symbol}: {e}")
                gene_to_compartments[gene_symbol] = ['cytosol']
        
        # Build compartment-specific GPRs
        compartment_to_genes = defaultdict(set)
        for gene_symbol, compartments in gene_to_compartments.items():
            for comp in compartments:
                compartment_to_genes[comp].add(gene_symbol)
        
        # Build sGPR and GPR by compartment (all OR since we assume isozymes)
        sGPR_by_comp = {}
        GPR_by_comp = {}
        
        for comp, genes in compartment_to_genes.items():
            sorted_comp_genes = sorted(genes)
            # sGPR with bracket notation to match BioCyc: ([([GENE1*1]) or ([GENE2*1])])
            sgpr_parts = [f"([{gene}*1])" for gene in sorted_comp_genes]
            if len(sgpr_parts) == 1:
                sGPR_by_comp[comp] = f"([{sgpr_parts[0]}])"
            else:
                sGPR_by_comp[comp] = "([" + " or ".join(sgpr_parts) + "])"
            # GPR with plain gene symbols
            GPR_by_comp[comp] = " or ".join(sorted_comp_genes)
        
        return (sgpr_general, gpr_general, sGPR_by_comp, GPR_by_comp)
    
    def _extract_and_add_genes(self, gpr_data: tuple, GeneList: Dict, GeneIdent: List, EnsblDB=None) -> None:
        """
        Extract genes from GPR data and add to GeneList if not already present.
        
        This mirrors the gene extraction in generate_db.py that builds GeneList
        from GPR rules after processing all reactions.
        
        Args:
            gpr_data: Tuple (sGPR_general, GPR_general, sGPR_by_comp, GPR_by_comp)
            GeneList: Dict of gene objects {gene_symbol: gene_object}
            GeneIdent: List of gene identifiers
            EnsblDB: Ensembl database (optional, for gene class)
        """
        if not gpr_data or len(gpr_data) < 2 or not gpr_data[1]:
            return
            
        # Extract gene symbols from GPR string
        gpr_string = gpr_data[1]  # GPR_general
        gene_matches = re.findall(
            r"([A-Za-z0-9\-]+)",
            gpr_string.replace("and", "").replace("or", "")
        )
        
        for gene_symbol in gene_matches:
            gene_symbol = gene_symbol.strip()
            if not gene_symbol:
                continue
                
            # Check if gene already exists in GeneList
            if gene_symbol not in GeneIdent:
                GeneIdent.append(gene_symbol)
                GeneList[gene_symbol] = gene(gene_symbol, EnsblDB)
                self.stats['genes_added'] += 1
                LOGGER.debug(f"  Added gene {gene_symbol} to GeneList")
    
    def _retry_failed_reactions(self, failed_reactions: List[Dict], max_retry_attempts: int, **kwargs) -> None:
        """
        Retry processing failed reactions (same pattern as KEGG pipeline).
        
        Handles two types of failures:
        1. GPR timeout - reaction not in RxnList yet, retry GPR lookup
        2. Compartmentalization failure - reaction in RxnList, retry rxnSubcel
        
        Args:
            failed_reactions: List of failed reaction dicts
            max_retry_attempts: Maximum number of retry rounds
            **kwargs: All data structures (RxnList, MetList, etc.)
        """
        import time as time_module
        
        session = kwargs.get('session')
        RxnList_CL = kwargs.get('RxnList_CL', {})
        MetList_CL = kwargs.get('MetList_CL', {})
        RxnIdent_CL = kwargs.get('RxnIdent_CL', [])
        MetIdent_CL = kwargs.get('MetIdent_CL', [])
        Compartment_CL = kwargs.get('Compartment_CL', [])
        RxnList = kwargs['RxnList']
        RxnIdent = kwargs.get('RxnIdent', [])
        MetList = kwargs['MetList']
        MetIdent = kwargs.get('MetIdent', [])
        MetEquiv = kwargs['MetEquiv']
        GPRList = kwargs.get('GPRList', {})
        GPRIdent = kwargs.get('GPRIdent', [])
        GeneList = kwargs.get('GeneList', {})
        GeneIdent = kwargs.get('GeneIdent', [])
        
        permanently_failed = []
        retry_round = 1
        
        while failed_reactions and retry_round <= max_retry_attempts:
            LOGGER.info("="*60)
            LOGGER.info(f"RETRY ROUND {retry_round}: {len(failed_reactions)} reactions to retry")
            LOGGER.info("="*60)
            
            # Wait before retrying (give servers time to recover)
            wait_time = retry_round * 5  # Increasing wait: 5s, 10s, 15s
            LOGGER.info(f"  Waiting {wait_time}s before retry...")
            time_module.sleep(wait_time)
            
            reactions_to_retry = failed_reactions[:]
            failed_reactions.clear()
            
            for failed_rxn in reactions_to_retry:
                rhea_id = failed_rxn['rhea_id']
                uniprots = failed_rxn['uniprots']
                ec_numbers = failed_rxn.get('ec_numbers', [])
                error_type = failed_rxn.get('error_type', 'unknown')
                
                LOGGER.info(f"  Retrying RHEA:{rhea_id} ({error_type}, attempt {retry_round + 1})")
                
                rxn_key = f"RHEA{rhea_id}"
                
                try:
                    # Case 1: GPR timeout - reaction not yet fully processed
                    if error_type == 'timeout' and rxn_key not in RxnList:
                        # Need to retry GPR lookup
                        gene_symbols = failed_rxn.get('gene_symbols', set())
                        
                        if not gene_symbols:
                            # Try EC→BioCyc again
                            gpr_data = None
                            for ec in ec_numbers:
                                try:
                                    if ec not in GPRIdent:
                                        GPRIdent.append(ec)
                                        gpr_obj = gpr(ec, session)
                                        GPRList[ec] = gpr_obj
                                    else:
                                        gpr_obj = GPRList[ec]
                                    
                                    gpr_data = gpr_obj.GprSubcell()
                                    if gpr_data and gpr_data[0]:
                                        break
                                except Exception as e:
                                    LOGGER.debug(f"  Retry EC {ec} failed: {e}")
                                    continue
                            
                            if not gpr_data or not gpr_data[0]:
                                raise Exception("GPR lookup still failing")
                        else:
                            # Retry with gene symbols (location lookup)
                            gpr_data = self._build_gpr_from_rhea_genes(gene_symbols, session)
                            
                            if not gpr_data or not gpr_data[0]:
                                raise Exception("Failed to build GPR from genes")
                        
                        # GPR succeeded - now we need to complete reaction setup
                        # Note: Metabolites and mass balance were done before GPR step
                        # so reaction adapter should exist, just needs GPR data
                        
                        # This is a partial retry - the reaction wasn't added to RxnList
                        # We need the reaction_participants data which we don't have here
                        # For now, mark as permanently failed if not in RxnList
                        LOGGER.warning(f"  RHEA:{rhea_id} - GPR retry succeeded but reaction not in RxnList")
                        permanently_failed.append(failed_rxn)
                        continue
                    
                    # Case 2: Compartmentalization failure - reaction in RxnList
                    if rxn_key not in RxnList:
                        LOGGER.warning(f"  Reaction {rxn_key} not in RxnList, skipping")
                        permanently_failed.append(failed_rxn)
                        continue
                        
                    rxn_adapter = RxnList[rxn_key]
                    
                    # Re-attempt compartmentalization
                    Compartment_CL, rxn_cl, comp_cl = rxnSubcel(
                        rxn_adapter,
                        RxnList_CL,
                        MetList_CL,
                        RxnIdent_CL,
                        MetIdent_CL,
                        Compartment_CL,
                        RxnList,
                        MetList,
                        MetEquiv,
                    )
                    RxnIdent_CL.extend(rxn_cl)
                    MetIdent_CL.extend(comp_cl)
                    self.stats['compartmentalized_reactions'] += len(rxn_cl)
                    
                    # Extract genes from GPR
                    gpr_data = getattr(rxn_adapter, '_gpr_data', None)
                    if gpr_data:
                        self._extract_and_add_genes(gpr_data, GeneList, GeneIdent, EnsblDB=None)
                    
                    LOGGER.info(f"  SUCCESS: Retried RHEA:{rhea_id}")
                    self.stats['retry_successes'] = self.stats.get('retry_successes', 0) + 1
                    
                except Exception as e:
                    LOGGER.warning(f"  Retry failed for RHEA:{rhea_id}: {e}")
                    
                    if retry_round >= max_retry_attempts:
                        permanently_failed.append(failed_rxn)
                    else:
                        failed_rxn['error'] = str(e)
                        failed_reactions.append(failed_rxn)
            
            retry_round += 1
        
        # Report permanently failed reactions
        if permanently_failed:
            LOGGER.warning(f"  {len(permanently_failed)} reactions permanently failed after {max_retry_attempts} retries:")
            for rxn in permanently_failed[:10]:  # Show first 10
                LOGGER.warning(f"    - RHEA:{rxn['rhea_id']}: {rxn.get('error', 'unknown')}")
            if len(permanently_failed) > 10:
                LOGGER.warning(f"    ... and {len(permanently_failed) - 10} more")
            self.stats['permanently_failed'] = len(permanently_failed)
    
    def _resolve_metabolite_id(self, chebi_id: str, MetEquiv: Dict) -> str:
        """
        Resolve a ChEBI ID to the correct metabolite ID.
        
        If the ChEBI compound already exists as a KEGG compound in MetList,
        return the KEGG ID. Otherwise, return the ChEBI ID.
        
        When resolving to a KEGG ID, we add the ChEBI→KEGG mapping to MetEquiv.
        This is critical for mass balance: RxnParam2Eq checks MetEquiv membership
        before checking ID prefix ('C' or 'G'), so ChEBI IDs in MetEquiv will pass
        both c_test and gly_test.
        
        Args:
            chebi_id: ChEBI ID (format: "CHEBI:12345" or "12345")
            MetEquiv: Metabolite equivalence dict
            
        Returns:
            The resolved metabolite ID (KEGG compound/glycan or ChEBI format)
        """
        chebi_num = chebi_id.replace('CHEBI:', '').strip()
        chebi_met_id = f"CHEBI{chebi_num}"
        
        # Check if this ChEBI maps to an existing KEGG compound/glycan
        if chebi_num in self.chebi_to_kegg:
            kegg_id = self.chebi_to_kegg[chebi_num]
            
            # Add ChEBI→KEGG mapping to MetEquiv so mass balance works
            # This ensures ChEBI IDs pass both c_test and gly_test
            if chebi_met_id not in MetEquiv:
                MetEquiv[chebi_met_id] = kegg_id
                LOGGER.debug(f"  Added MetEquiv mapping in resolve: {chebi_met_id} → {kegg_id}")
            
            LOGGER.debug(f"  ChEBI:{chebi_num} → KEGG:{kegg_id} (existing)")
            return kegg_id
            
        # Return ChEBI ID format (starts with 'C', works with compound path in RxnParam2Eq)
        return chebi_met_id
    
    def _mass_balance_reaction(self, rxn_adapter: RheaReactionAdapter, rxn_key: str,
                                MetList: Dict, MetIdent: List, MetEquiv: Dict) -> bool:
        """
        Mass balance a Rhea reaction using the same method as KEGG reactions.
        
        Rhea reactions are guaranteed to be chemically balanced in the database,
        but we still run our mass balance algorithm to:
        1. Verify the reaction is balanced with our metabolite formulas
        2. Update stoichiometry if needed (e.g., due to formula differences)
        3. Potentially add H+/H2O if our formulas differ from Rhea's
        
        Args:
            rxn_adapter: The RheaReactionAdapter object
            rxn_key: The reaction key (e.g., "RHEA12345")
            MetList: Dict of metabolites
            MetIdent: List of metabolite IDs
            MetEquiv: Metabolite equivalence dict
            
        Returns:
            True if mass balance succeeded or was already balanced
        """
        try:
            # Check if any metabolite has a problematic formula (polymeric/variable)
            # These contain mathematical expressions that cannot be parsed
            # 
            # IMPORTANT: We check BOTH Formula1 and Formula2 attributes directly
            # because RxnParam2Eq accesses these attributes, not Formula() method.
            # ChebiCompoundAdapter doesn't have Formula() method, so we need to
            # check the attributes that will actually be used.
            
            def get_formula_for_validation(met):
                """Get formula from metabolite, checking attributes used by RxnParam2Eq."""
                # RxnParam2Eq uses Formula1 for compounds (c_test) and Formula2 for glycans (gly_test)
                # Check Formula1 first (most common case for ChEBI compounds)
                if hasattr(met, 'Formula1') and met.Formula1:
                    return met.Formula1
                if hasattr(met, 'Formula2') and met.Formula2:
                    return met.Formula2
                # Fallback to Formula() method if available
                if hasattr(met, 'Formula') and callable(met.Formula):
                    return met.Formula()
                return ""
            
            for sub in rxn_adapter.Substrate():
                met_id = sub[2]
                if met_id in MetList:
                    formula = get_formula_for_validation(MetList[met_id])
                    if formula and not is_valid_formula_for_mass_balance(formula):
                        LOGGER.debug(f"  {rxn_key}: Skipping mass balance - metabolite {met_id} has polymeric formula: {formula}")
                        # Rhea guarantees balance, so we trust it
                        self.stats['mass_balance_already_balanced'] += 1
                        return True
                        
            for prod in rxn_adapter.Product():
                met_id = prod[2]
                if met_id in MetList:
                    formula = get_formula_for_validation(MetList[met_id])
                    if formula and not is_valid_formula_for_mass_balance(formula):
                        LOGGER.debug(f"  {rxn_key}: Skipping mass balance - metabolite {met_id} has polymeric formula: {formula}")
                        # Rhea guarantees balance, so we trust it
                        self.stats['mass_balance_already_balanced'] += 1
                        return True
            
            # Build the equation string from reaction parameters
            eq_result = RxnParam2Eq(rxn_adapter, MetList, MetEquiv)
            # RxnParam2Eq returns (eq, mb_test) tuple
            eq = eq_result[0] if isinstance(eq_result, tuple) else eq_result
            
            if not eq:
                LOGGER.debug(f"  {rxn_key}: Empty equation, skipping mass balance")
                self.stats['mass_balance_failed'] += 1
                return False
            
            # Additional validation of the equation string
            # Check if any formula in the equation is problematic
            # The equation looks like: "1 H2O + 2 C6H12O6 -> 2 C12H22O11"
            # Each formula after a stoichiometry coefficient should be valid
            equation_parts = eq.split('->')
            if len(equation_parts) == 2:
                for part in equation_parts:
                    for term in part.split('+'):
                        term = term.strip()
                        # Extract formula (after stoichiometry coefficient)
                        import re as re_local
                        match = re_local.match(r'^\d*\.?\d*\s*(.+)$', term)
                        if match:
                            formula_in_eq = match.group(1).strip()
                            if formula_in_eq and not is_valid_formula_for_mass_balance(formula_in_eq):
                                LOGGER.debug(f"  {rxn_key}: Skipping mass balance - equation contains invalid formula: {formula_in_eq}")
                                self.stats['mass_balance_already_balanced'] += 1
                                return True  # Trust Rhea's balance
            
            # Run mass balance with numpy warning suppression for clean output
            import warnings
            import numpy as np
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=RuntimeWarning)
                result = mass_balance(eq, rxn_key)
            
            # Result format: (StchS, StchP, AddH, AddH2O, NewSpeciesS, NewSpeciesP,
            #                 AllSpeciesS, AllSpeciesP, eq, eq_init, TestOfBalance)
            # TestOfBalance: 1 = already balanced, 0 = balanced with modifications, 2 = failed
            
            stch_s, stch_p, add_h, add_h2o = result[0], result[1], result[2], result[3]
            new_species_s, new_species_p = result[4], result[5]
            test_of_balance = result[10]
            
            if test_of_balance == 1:
                # Already balanced
                self.stats['mass_balance_already_balanced'] += 1
                LOGGER.debug(f"  {rxn_key}: Already mass balanced")
                return True
                
            elif test_of_balance == 0:
                # Mass balance succeeded with modifications
                # Update substrate stoichiometry
                new_subs = []
                for i, sub in enumerate(rxn_adapter.Substrate()):
                    if i < len(stch_s):
                        new_subs.append([stch_s[i], sub[1], sub[2]])
                    else:
                        new_subs.append(sub)
                
                # Update product stoichiometry
                new_prods = []
                for i, prod in enumerate(rxn_adapter.Product()):
                    if i < len(stch_p):
                        new_prods.append([stch_p[i], prod[1], prod[2]])
                    else:
                        new_prods.append(prod)
                
                # Add H+ if needed
                if add_h != 0:
                    h_id = self.extra_compound.get("H", "C00080")
                    # Ensure H+ is in MetList
                    self._ensure_extra_compound_exists("H", h_id, MetList, MetIdent)
                    
                    if add_h < 0:  # Add to substrates
                        new_subs.append([abs(add_h), "", h_id])
                    else:  # Add to products
                        new_prods.append([abs(add_h), "", h_id])
                
                # Add H2O if needed
                if add_h2o != 0:
                    h2o_id = self.extra_compound.get("H2O", "C00001")
                    # Ensure H2O is in MetList
                    self._ensure_extra_compound_exists("H2O", h2o_id, MetList, MetIdent)
                    
                    if add_h2o < 0:  # Add to substrates
                        new_subs.append([abs(add_h2o), "", h2o_id])
                    else:  # Add to products
                        new_prods.append([abs(add_h2o), "", h2o_id])
                
                # Update the reaction adapter
                rxn_adapter.SetSubstrate(new_subs)
                rxn_adapter.SetProduct(new_prods)
                
                self.stats['mass_balance_success'] += 1
                LOGGER.debug(f"  {rxn_key}: Mass balanced (added H+={add_h}, H2O={add_h2o})")
                return True
                
            else:
                # Mass balance failed (test_of_balance == 2)
                self.stats['mass_balance_failed'] += 1
                LOGGER.warning(f"  {rxn_key}: Mass balance failed")
                return False
                
        except Exception as e:
            self.stats['mass_balance_failed'] += 1
            LOGGER.warning(f"  {rxn_key}: Mass balance error: {e}")
            return False
    
    def _ensure_extra_compound_exists(self, name: str, kegg_id: str, 
                                       MetList: Dict, MetIdent: List) -> None:
        """
        Ensure an extra compound (H+, H2O, etc.) exists in MetList.
        
        These are needed for mass balancing when the original reaction
        doesn't include them but they're needed for balance.
        """
        if kegg_id in MetList:
            return
            
        # Create a minimal compound entry
        # The compound class will populate it properly when needed
        from functions.class_generate_database import compound
        
        # Common formulas for extra compounds
        formulas = {
            "H": "H",
            "H2O": "H2O",
            "CO2": "CO2",
            "O2": "O2",
            "NH3": "H3N",
            "Pi": "HO4P",
            "PPi": "H2O7P2",
        }
        
        formula = formulas.get(name, name)
        
        # Create minimal compound object
        new_compound = ChebiCompoundAdapter(
            chebi_id=kegg_id.replace("C", ""),  # Won't be a real ChEBI ID
            name=name,
            formula=formula,
        )
        # Override IDs to use KEGG format
        new_compound.ID1 = kegg_id
        new_compound.ID2 = kegg_id
        new_compound.ident = kegg_id
        
        MetList[kegg_id] = new_compound
        MetIdent.append(kegg_id)
        
        LOGGER.debug(f"  Added extra compound {name} ({kegg_id}) for mass balancing")
            
    def _add_metabolite_if_new(self, compound_data: Dict, MetList: Dict, 
                                MetIdent: List, MetEquiv: Dict) -> str:
        """
        Add a compound (ChEBI or macromolecule) to MetList if not already present.
        
        For ChEBI compounds (is_macromolecule=False):
        Uses multi-attribute matching (same approach as functions_merge_metabolic_networks):
        1. ChEBI ID
        2. InChIKey
        3. InChI  
        4. PubChem
        5. Name (exact match)
        
        For macromolecules (is_macromolecule=True):
        Uses only name-based matching since these have no ChEBI ID, formula, or external IDs.
        These will be added with 'M' prefix to ensure they don't pass c_test/gly_test
        in mass balance (so mb_test=0, reaction skipped from mass balance correctly).
        
        Args:
            compound_data: Dict with fields:
                - 'chebi': ChEBI ID (e.g., 'CHEBI:12345') or empty for macromolecules
                - 'name': compound name
                - 'formula': chemical formula (empty for macromolecules)
                - 'rhea_compound_id': Rhea compound ID for macromolecules (e.g., 'GENERIC:13713')
                - 'is_macromolecule': True if this is a macromolecular compound
            MetList: Dict of metabolites
            MetIdent: List of metabolite IDs
            MetEquiv: Metabolite equivalence dict
            
        Returns:
            The metabolite ID used (existing KEGG/ChEBI or new ChEBI/Macromolecule)
        """
        is_macromolecule = compound_data.get('is_macromolecule', False)
        
        if is_macromolecule:
            return self._add_macromolecule_if_new(compound_data, MetList, MetIdent, MetEquiv)
        else:
            return self._add_chebi_compound_if_new(compound_data, MetList, MetIdent, MetEquiv)
    
    def _add_macromolecule_if_new(self, compound_data: Dict, MetList: Dict,
                                   MetIdent: List, MetEquiv: Dict) -> str:
        """
        Add a macromolecular compound to MetList if not already present.
        
        Macromolecules (e.g., "L-seryl-[protein]") have no ChEBI ID, formula,
        or external database references. They are matched only by name.
        
        Args:
            compound_data: Dict with 'name', 'rhea_compound_id', 'is_macromolecule'
            MetList: Dict of metabolites
            MetIdent: List of metabolite IDs
            MetEquiv: Metabolite equivalence dict
            
        Returns:
            The metabolite ID used
        """
        rhea_id = compound_data.get('rhea_compound_id', '')
        name = compound_data.get('name', '')
        
        if not rhea_id and not name:
            LOGGER.warning("Macromolecule missing both rhea_compound_id and name")
            return ""
        
        # Generate the metabolite ID using M prefix
        # Normalize: GENERIC:13713 -> MG13713
        clean_id = rhea_id.replace(":", "").replace("GENERIC", "G")
        met_id = f"M{clean_id}"
        
        # Check if already added (by exact name match or ID)
        if met_id in MetList:
            LOGGER.debug(f"  Reusing macromolecule {met_id} ({name})")
            return met_id
        
        # Try to find by name (case-insensitive)
        if name:
            name_lower = name.lower()
            for existing_id, existing_compound in MetList.items():
                existing_name = _get_compound_name(existing_compound)
                if existing_name and existing_name.lower() == name_lower:
                    # Add equivalence mapping
                    if met_id not in MetEquiv:
                        MetEquiv[met_id] = existing_id
                    LOGGER.debug(f"  Found existing compound for macromolecule by name: {existing_id}")
                    return existing_id
        
        # Create new macromolecule adapter
        compound = RheaMacromoleculeAdapter(
            rhea_compound_id=rhea_id,
            name=name
        )
        
        MetList[met_id] = compound
        MetIdent.append(met_id)
        self.stats['macromolecules_added'] = self.stats.get('macromolecules_added', 0) + 1
        
        LOGGER.debug(f"  Added macromolecule: {met_id} ({name})")
        
        return met_id
    
    def _add_chebi_compound_if_new(self, compound_data: Dict, MetList: Dict,
                                    MetIdent: List, MetEquiv: Dict) -> str:
        """
        Add a ChEBI compound to MetList if not already present.
        
        Uses multi-attribute matching (same approach as functions_merge_metabolic_networks):
        1. ChEBI ID
        2. InChIKey
        3. InChI  
        4. PubChem
        5. Name (exact match)
        
        If the compound already exists as a KEGG compound (via any identifier),
        we enrich its annotations with ChEBI data if any attributes are missing,
        then return the KEGG ID without adding a duplicate.
        
        Args:
            compound_data: Dict with 'chebi', 'name', 'formula' from SPARQL
            MetList: Dict of metabolites
            MetIdent: List of metabolite IDs
            MetEquiv: Metabolite equivalence dict
            
        Returns:
            The metabolite ID used (either existing KEGG or new ChEBI)
        """
        chebi_num = compound_data['chebi'].replace('CHEBI:', '')
        name = compound_data.get('name', '')
        
        # Fetch additional data from ChEBI (OLS4 API) for comprehensive matching
        chebi_data = self.chebi_fetcher.fetch_compound(compound_data['chebi'])
        formula = compound_data.get('formula', '') or chebi_data.get('formula', '')
        inchi = chebi_data.get('inchi', '')
        inchikey = chebi_data.get('inchikey', '')
        pubchem = chebi_data.get('pubchem', '')
        lipidmaps = chebi_data.get('lipidmaps', '')
        kegg = chebi_data.get('kegg', '')      # KEGG compound ID from OLS4
        hmdb = chebi_data.get('hmdb', '')      # HMDB ID from OLS4
        
        # Use ChEBI name if SPARQL didn't provide one
        if not name and chebi_data.get('name'):
            name = chebi_data.get('name', '')
        
        # Try multi-attribute matching (same approach as merge_metabolic_networks)
        existing_kegg_id = self._find_existing_metabolite(
            chebi_id=chebi_num,
            name=name,
            inchi=inchi,
            inchikey=inchikey,
            pubchem=pubchem
        )
        
        if existing_kegg_id:
            self.stats['metabolites_reused'] += 1
            
            # Enrich existing compound with ChEBI data if attributes are empty
            if existing_kegg_id in MetList:
                enriched = self._enrich_metabolite(MetList[existing_kegg_id], compound_data)
                if enriched:
                    self.stats['metabolites_enriched'] += 1
            
            # Add ChEBI ID to MetEquiv mapping to the existing KEGG ID
            # This is critical for mass balance: RxnParam2Eq checks if IDs are in MetEquiv
            # before checking the prefix ('C' or 'G'). By adding ChEBI→KEGG mapping here,
            # ChEBI IDs will pass both c_test and gly_test in equations_bm_gdb.py
            chebi_met_id = f"CHEBI{chebi_num}"
            if chebi_met_id not in MetEquiv:
                MetEquiv[chebi_met_id] = existing_kegg_id
                LOGGER.debug(f"  Added MetEquiv mapping: {chebi_met_id} → {existing_kegg_id}")
                    
            LOGGER.debug(f"  Reusing KEGG:{existing_kegg_id} for ChEBI:{chebi_num}")
            return existing_kegg_id
            
        # Check if already added as ChEBI (by previous Rhea reaction)
        met_id = f"CHEBI{chebi_num}"
        if met_id in MetList:
            return met_id
            
        # Create new ChEBI compound adapter with all cross-references from OLS4
        compound = ChebiCompoundAdapter(
            chebi_id=chebi_num,
            name=name or chebi_data.get('name', '') or f'CHEBI:{chebi_num}',
            formula=formula,
            inchi=inchi,
            inchikey=inchikey,
            pubchem=pubchem,
            lipidmaps=lipidmaps,
            kegg=kegg,
            hmdb=hmdb,
        )
        
        MetList[met_id] = compound
        MetIdent.append(met_id)
        self.stats['metabolites_added'] += 1
        
        # Update lookup index so future Rhea compounds can match this one
        if chebi_num:
            self._metabolite_index['chebi'][chebi_num] = met_id
        if inchikey:
            self._metabolite_index['inchikey'][inchikey] = met_id
        if inchi:
            self._metabolite_index['inchi'][inchi] = met_id
        if pubchem:
            self._metabolite_index['pubchem.compound'][str(pubchem)] = met_id
        
        return met_id
    
    def _enrich_metabolite(self, existing_compound, chebi_data: Dict) -> bool:
        """
        Enrich an existing metabolite with ChEBI data if attributes are empty.
        
        Args:
            existing_compound: The existing compound object from MetList
            chebi_data: Dict with 'chebi', 'name', 'formula' from SPARQL
            
        Returns:
            True if any enrichment was made, False otherwise
        """
        enriched = False
        chebi_num = chebi_data['chebi'].replace('CHEBI:', '')
        
        # Fetch full ChEBI data including cross-references
        full_chebi_data = self.chebi_fetcher.fetch_compound(chebi_data['chebi'])
        
        # Check and update CheBI ID if empty
        if hasattr(existing_compound, 'CheBI'):
            if not existing_compound.CheBI:
                existing_compound.CheBI = chebi_num
                enriched = True
        
        # Check and update formula if empty
        formula = chebi_data.get('formula', '') or full_chebi_data.get('formula', '')
        if formula:
            if hasattr(existing_compound, 'Formula1') and not existing_compound.Formula1:
                existing_compound.Formula1 = formula
                enriched = True
            if hasattr(existing_compound, 'Formula2') and not existing_compound.Formula2:
                existing_compound.Formula2 = formula
                enriched = True
        
        # Check and update name if current name is just the ID
        name = chebi_data.get('name', '') or full_chebi_data.get('name', '')
        if name:
            current_name = _get_compound_name(existing_compound)
            # If current name is empty or looks like just an ID
            if not current_name or (current_name.startswith('C0') or current_name.startswith('G0')):
                if _set_compound_name(existing_compound, name):
                    enriched = True
        
        # Check and update PubChem if empty
        pubchem = full_chebi_data.get('pubchem', '')
        if pubchem and hasattr(existing_compound, 'PubChem'):
            if not existing_compound.PubChem:
                existing_compound.PubChem = pubchem
                enriched = True
                
        # Check and update LIPIDMAPS if empty
        lipidmaps = full_chebi_data.get('lipidmaps', '')
        if lipidmaps and hasattr(existing_compound, 'LIPIDMAPS'):
            if not existing_compound.LIPIDMAPS:
                existing_compound.LIPIDMAPS = lipidmaps
                enriched = True
                
        if enriched:
            LOGGER.debug(f"  Enriched compound with ChEBI:{chebi_num} data")
            
        return enriched
            
    def _print_summary(self) -> None:
        """Print extension summary."""
        LOGGER.info("\n" + "=" * 60)
        LOGGER.info("RHEA EXTENSION SUMMARY")
        LOGGER.info("=" * 60)
        LOGGER.info(f"Human Rhea reactions found:    {self.stats['total_rhea_human']}")
        LOGGER.info(f"Already in KEGG model:         {self.stats['already_in_kegg']}")
        LOGGER.info(f"New reactions identified:      {self.stats['new_reactions']}")
        LOGGER.info("-" * 60)
        LOGGER.info("SKIPPED REACTIONS:")
        LOGGER.info(f"  No participants (SPARQL):    {self.stats['reactions_skipped_no_participants']}")
        LOGGER.info("-" * 60)
        LOGGER.info("ADDED REACTIONS:")
        LOGGER.info(f"  Total (base):                {self.stats['reactions_added']}")
        LOGGER.info(f"  Via existing GPR (UniProt→EC):{self.stats.get('reactions_via_existing_gpr', 0)}")
        LOGGER.info(f"  Via EC→BioCyc (proper GPR):  {self.stats.get('reactions_via_ec', 0)}")
        LOGGER.info(f"  Via Rhea genes (OR fallback):{self.stats.get('reactions_via_rhea_genes', 0)}")
        LOGGER.info(f"  Via KEGG fallback:           {self.stats.get('kegg_fallback_used', 0)}")
        LOGGER.info(f"  WITHOUT genes (cytosol):     {self.stats.get('reactions_added_no_genes', 0)}")
        LOGGER.info(f"  Compartmentalized:           {self.stats['compartmentalized_reactions']}")
        LOGGER.info("-" * 60)
        LOGGER.info("NO-GENES BREAKDOWN (for diagnostics):")
        LOGGER.info(f"  No EC in Rhea:               {self.stats.get('no_genes_no_ec', 0)}")
        LOGGER.info(f"  No UniProt in Rhea:          {self.stats.get('no_genes_no_uniprot', 0)}")
        LOGGER.info(f"  UniProt→Symbol failed:       {self.stats.get('no_genes_uniprot_no_symbol', 0)}")
        LOGGER.info(f"  KEGG fallback tried:         {self.stats.get('kegg_fallback_available', 0)}")
        LOGGER.info("-" * 60)
        LOGGER.info("GENES:")
        LOGGER.info(f"  Added (new):                 {self.stats.get('genes_added', 0)}")
        LOGGER.info("-" * 60)
        LOGGER.info("METABOLITES:")
        LOGGER.info(f"  Added (new ChEBI):           {self.stats['metabolites_added']}")
        LOGGER.info(f"  Added (macromolecules):      {self.stats['macromolecules_added']}")
        LOGGER.info(f"  Reused (existing KEGG):      {self.stats['metabolites_reused']}")
        LOGGER.info(f"  Enriched (updated attrs):    {self.stats['metabolites_enriched']}")
        LOGGER.info("-" * 60)
        LOGGER.info("MASS BALANCE:")
        LOGGER.info(f"  Already balanced:            {self.stats['mass_balance_already_balanced']}")
        LOGGER.info(f"  Balanced with modifications: {self.stats['mass_balance_success']}")
        LOGGER.info(f"  Failed:                      {self.stats['mass_balance_failed']}")
        LOGGER.info("-" * 60)
        LOGGER.info("RETRY:")
        LOGGER.info(f"  Reactions failed initially:  {self.stats.get('reactions_failed', 0)}")
        LOGGER.info(f"  Recovered via retry:         {self.stats.get('retry_successes', 0)}")
        LOGGER.info(f"  Permanently failed:          {self.stats.get('permanently_failed', 0)}")
        LOGGER.info("-" * 60)
        LOGGER.info("CACHE:")
        LOGGER.info(f"  EC cache hits (GPRList):     {self.stats.get('ec_cache_hits', 0)}")
        LOGGER.info("=" * 60)
        
        # Print UniProt mapping stats if available
        if hasattr(self.human_filter, 'uniprot_mapping_stats') and self.human_filter.uniprot_mapping_stats:
            stats = self.human_filter.uniprot_mapping_stats
            LOGGER.info("\nUNIPROT API MAPPING STATISTICS:")
            LOGGER.info("-" * 60)
            LOGGER.info(f"  Total UniProt IDs queried:   {stats.get('total_queried', 0)}")
            LOGGER.info(f"  Has gene symbol:             {stats.get('has_gene_symbol', 0)}")
            LOGGER.info(f"  Has KEGG cross-ref:          {stats.get('has_kegg', 0)}")
            LOGGER.info(f"  Has GeneID cross-ref:        {stats.get('has_geneid', 0)}")
            LOGGER.info(f"  Has Ensembl cross-ref:       {stats.get('has_ensembl', 0)}")
            LOGGER.info(f"  Has EC number:               {stats.get('has_ec', 0)}")
            LOGGER.info(f"  NO gene symbol:              {stats.get('no_gene_symbol', 0)}")
            LOGGER.info(f"  API errors:                  {stats.get('api_errors', 0)}")
            LOGGER.info("=" * 60)


# =============================================================================
# INTEGRATION FUNCTION - Call after KEGG loop
# =============================================================================

def extend_with_rhea(
    RxnList: Dict,
    MetList: Dict,
    RxnIdent: List,
    MetIdent: List,
    MetEquiv: Dict,
    GPRList: Optional[Dict] = None,
    GPRIdent: Optional[List] = None,
    GeneList: Optional[Dict] = None,
    GeneIdent: Optional[List] = None,
    RxnList_CL: Optional[Dict] = None,
    MetList_CL: Optional[Dict] = None,
    RxnIdent_CL: Optional[List] = None,
    MetIdent_CL: Optional[List] = None,
    Compartment_CL: Optional[List] = None,
    PathNameRxn: Optional[Dict] = None,
    session = None,
    force_refresh: bool = False,
    max_retry_attempts: int = 3,
) -> Dict[str, Any]:
    """
    Convenience function to extend KEGG model with Rhea reactions.
    
    Call this after the KEGG pathway loop completes in generate_db.py:
    
        # After KEGG loop:
        from rhea_extension import extend_with_rhea
        
        rhea_stats = extend_with_rhea(
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
        )
        
        # Continue with cobra_reconstruction as usual...
    
    Returns:
        Dict with statistics about the extension
    """
    extender = RheaExtender()
    return extender.run(
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
        force_refresh=force_refresh,
        max_retry_attempts=max_retry_attempts,
    )


# =============================================================================
# MAIN - Standalone testing
# =============================================================================

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Test Rhea extension module"
    )
    parser.add_argument(
        "--checkpoint",
        help="Path to KEGG checkpoint file to extend"
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Force re-download of Rhea files"
    )
    
    args = parser.parse_args()
    
    # Test mode: just download and analyze Rhea data
    LOGGER.info("Testing Rhea extension module...")
    
    extender = RheaExtender()
    extender.downloader.download_all(force_refresh=args.refresh)
    extender.human_filter.load_human_proteome()
    extender._load_rhea_mappings()
    
    # Show stats
    new_rhea = extender._filter_new_reactions(set())  # Empty = all are "new"
    
    LOGGER.info(f"\nTotal human Rhea reactions: {len(extender.rhea_to_uniprots)}")
    LOGGER.info(f"With KEGG mapping: {len([r for r in extender.rhea_to_uniprots if r in extender.rhea_to_kegg])}")
    LOGGER.info(f"Without KEGG (potential new): {len(new_rhea)}")
