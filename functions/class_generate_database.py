#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import re
import urllib.request
from typing import TYPE_CHECKING

import pdb

from functions.gpr.auth_gpr import getGPR, setup_biocyc_session

if TYPE_CHECKING:  # Avoid circular import at runtime, keep type hints available
    from functions import function_bm_gdb as _bm_mod
    from functions.gpr import get_location_def as _loc_mod


def _bm():
    """Lazy import to avoid circular dependency with function_bm_gdb."""
    from functions import function_bm_gdb  # type: ignore

    return function_bm_gdb


def _loc():
    """Lazy import to avoid circular dependency with get_location_def."""
    from functions.gpr import get_location_def  # type: ignore

    return get_location_def


class pathway(object):
    def __init__(self, url, time, ID, urlReferer, PathName):
        bm = _bm()
        self.pagina = bm.getHtml(url, time, urlReferer)
        if type(self.pagina) == bytes:
            self.pagina = self.pagina.decode("utf-8")
        self.link = bm.getLinkPath(self.pagina)
        self.ID = ID
        self.PathName = PathName

    def ID(self):  # KEGG
        try:
            return self.ID
        except Exception:
            return ""

    def PathName(self):  # KEGG
        try:
            return self.PathName
        except Exception:
            return ""

    def Compounds(self):  # KEGG
        try:
            return [x for x in self.link[0]]
        except Exception:
            return ""

    def Reactions(self):  # KEGG
        try:
            return [x for x in self.link[1]]
        except Exception:
            return ""


class reaction(object):
    def __init__(self, url, time, ID, path, termdyn, *newparam):
        bm = _bm()
        self.pagina = bm.getHtml(url, time)
        self.link = bm.getReacParam(self.pagina, time)
        self.ID = ID
        self.path = path
        self.termdyn = termdyn
        self.newparam = newparam

    def ID(self):  # Patway-KEGG
        try:
            return self.ID
        except Exception:
            return ""

    def Name(self):  # KEGG
        try:
            if self.link[3]:
                N = self.link[3]
            else:
                N = self.ID
            return N
        except Exception:
            return ""

    def EC(self):  # KEGG
        try:
            # 			return self.link[4][0][1]
            return [x[1] for x in self.link[4]]
        except Exception:
            return ""

    def GPR(self):  # MetaCyc
        try:
            return self.newparam[0][0], self.newparam[0][1]
        except Exception:
            return ""

    # 	def GPR2(self): #MetaCyc
    # 		try:
    # 			return self.link[8]
    # 		except Exception:
    # 			return ''
    def Termodyn(self):  # Patway-KEGG
        try:
            return self.termdyn
        except Exception:
            return ""

    def Substrate(self):  # KEGG
        try:
            S = [[x[0]] + [x[1]] + [x[2]] for x in self.link[1]]
            return S  # link+ID+stoichometry
        except Exception:
            return ""

    def SetSubstrate(self, substrate):  # KEGG'
        self.link[1] = [[x[0]] + [x[1]] + [x[2]] for x in substrate]
        self.link[0] = self.link[1] + self.link[2]
        return self.link[1]

    def Product(self):  # KEGG
        try:
            P = [[x[0]] + [x[1]] + [x[2]] for x in self.link[2]]
            return P  # link+ID+stoichometry
        except Exception:
            return ""

    def SetProduct(self, product):  # KEGG
        self.link[2] = [[x[0]] + [x[1]] + [x[2]] for x in product]
        self.link[0] = self.link[1] + self.link[2]
        return self.link[2]

    def Pathway(self):  # Patway-KEGG
        try:
            return self.path
        except Exception:
            return ""

    def Subcel(self):  # Uniprot
        try:
            return self.newparam[0][2], self.newparam[0][3]
        except Exception:
            return ""

    def Equivalent(self):  # KEGG
        try:
            return self.link[5]
        except Exception:
            return ""

    def GTest(self):  # KEGG
        try:
            return self.link[6]
        except Exception:
            return ""

    def CTest(self):  # KEGG
        try:
            return self.link[7]
        except Exception:
            return ""

    def MBTest(self):  # KEGG
        try:
            return self.ID[0]
        except Exception:
            return ""


class gpr(object):
    def __init__(self, ec, session=None, impose_locations=1):
        """Initialize GPR object for an EC number.
        
        Args:
            ec: EC number to query
            session: BioCyc session (optional)
            impose_locations: Compartmentalization mode
                1 = RESTRICTED: Only use compartments from Excel lookup (default)
                0 = UNRESTRICTED: Keep all compartment names as-is
        """
        self.ec = ec
        # Use provided session or create new one
        if session is None:
            session = setup_biocyc_session()
        self.GPRPAss = getGPR(self.ec, session)
        location_module = _loc()
        
        # Determine project root path for bb.pickle file and Excel file
        current_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.join(current_dir, "..")
        bb_pickle_path = os.path.join(project_root, "files", "bb.pickle")
        excel_file = os.path.join(project_root, "files", "ListOfCompartments_sept2024.xlsx")
        comp_abb_file = os.path.join(project_root, "files", "compartments_info.txt")
        
        self.Subcell = location_module.getLocationnew(
            self.GPRPAss[3],
            self.GPRPAss[1],
            self.GPRPAss[2],
            impose_locations,  # 1 = restricted, 0 = unrestricted
            bb_pickle_path,
            session,
            None,  # ensembl_cache
            excel_file,
            "Def-Compartments",  # compartments_sheet_name
            {},  # comp_dict (empty, will be built by function)
            comp_abb_file,
        )

    def EC(self):  # Patway-KEGG
        try:
            return self.ec
        except Exception:
            return ""

    def GprSubcell(self):  # MetaCyc
        try:
            return (
                self.GPRPAss[3],
                self.GPRPAss[4],
                self.Subcell[0],
                self.Subcell[1],
            )  # General gene rule with stoichometry + General gene rule without stoichometry + Subcell-specific rule with stoichometry + Subcell-specific rule without stoichometry
        except Exception:
            return ""


class gene(object):
    # Class-level cache shared across all instances
    _location_cache = {}  # {gene_symbol: {'compartments': [...], 'biocyc_id': '...'}}
    _cache_stats = {'hits': 0, 'misses': 0, 'queries': 0}
    
    def __init__(self, gene, db):
        self.gene = gene
        self.db = db

    def Name(self):  # Patway-KEGG
        try:
            return self.gene
        except Exception:
            return ""

    @classmethod
    def get_locations(cls, gene_symbol, biocyc_id=None, session=None):
        """Get subcellular locations for a gene with caching.
        
        This method checks if the gene's locations are already cached.
        If not, it returns None to signal that the caller should perform
        the query and then cache the result using cache_location().
        
        Args:
            gene_symbol: Gene symbol (e.g., 'BRCA1')
            biocyc_id: BioCyc gene ID (optional, for cache metadata)
            session: requests.Session for BioCyc queries (unused in cache lookup)
            
        Returns:
            list or None: List of compartment names if cached, None if not cached
        """
        import logging
        LOGGER = logging.getLogger(__name__)
        
        if gene_symbol in cls._location_cache:
            cls._cache_stats['hits'] += 1
            cached_data = cls._location_cache[gene_symbol]
            LOGGER.debug(f"Cache HIT for {gene_symbol}: {cached_data['compartments']}")
            return cached_data['compartments']
        
        cls._cache_stats['misses'] += 1
        LOGGER.debug(f"Cache MISS for {gene_symbol} - needs querying")
        return None  # Signal to caller: not cached, perform query
    
    @classmethod
    def cache_location(cls, gene_symbol, compartments, biocyc_id=None):
        """Store gene location in cache.
        
        Args:
            gene_symbol: Gene symbol (e.g., 'BRCA1')
            compartments: List of compartment names (e.g., ['mitochondria', 'cytosol'])
            biocyc_id: BioCyc gene ID (optional, for metadata)
        """
        import logging
        LOGGER = logging.getLogger(__name__)
        
        cls._location_cache[gene_symbol] = {
            'compartments': compartments,
            'biocyc_id': biocyc_id
        }
        cls._cache_stats['queries'] += 1
        LOGGER.debug(f"Cached locations for {gene_symbol}: {compartments}")
    
    @classmethod
    def clear_cache(cls):
        """Clear location cache (useful for testing)."""
        cls._location_cache.clear()
        cls._cache_stats = {'hits': 0, 'misses': 0, 'queries': 0}
    
    @classmethod
    def get_cache_stats(cls):
        """Get cache performance statistics.
        
        Returns:
            dict: Statistics including hits, misses, queries, hit_rate, and cache_size
        """
        total = cls._cache_stats['hits'] + cls._cache_stats['misses']
        hit_rate = (cls._cache_stats['hits'] / total * 100) if total > 0 else 0
        return {
            **cls._cache_stats,
            'hit_rate': f"{hit_rate:.1f}%",
            'cache_size': len(cls._location_cache)
        }
    
    @classmethod
    def save_cache(cls, filepath):
        """Save location cache to pickle file.
        
        Args:
            filepath: Path to save the cache pickle file
        """
        import pickle
        import logging
        LOGGER = logging.getLogger(__name__)
        
        try:
            with open(filepath, 'wb') as f:
                pickle.dump(cls._location_cache, f)
            LOGGER.info(f"Saved gene location cache ({len(cls._location_cache)} genes) to {filepath}")
        except Exception as e:
            LOGGER.error(f"Failed to save gene location cache: {e}")
    
    @classmethod
    def load_cache(cls, filepath):
        """Load location cache from pickle file.
        
        Args:
            filepath: Path to load the cache pickle file from
            
        Returns:
            bool: True if cache was loaded successfully, False otherwise
        """
        import pickle
        import os
        import logging
        LOGGER = logging.getLogger(__name__)
        
        if os.path.exists(filepath):
            try:
                with open(filepath, 'rb') as f:
                    cls._location_cache = pickle.load(f)
                LOGGER.info(f"Loaded gene location cache ({len(cls._location_cache)} genes) from {filepath}")
                return True
            except Exception as e:
                LOGGER.error(f"Failed to load gene location cache: {e}")
                return False
        else:
            LOGGER.debug(f"Gene location cache file not found: {filepath}")
            return False

    def Ensg(self):  # MetaCyc
        try:
            iiii2 = re.search(
                r"gene=([A-Z0-9]+)",
                str(
                    urllib.request.urlopen(
                        "https://www.genome.jp/dbget-bin/www_bget?hsa+" + self.gene
                    ).read()
                ),
            )  # .group(1)
            if not iiii2:
                iiii2 = re.search(
                    r"(ENSG[0-9]+)",
                    str(
                        urllib.request.urlopen(
                            "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g="
                            + self.gene
                        ).read()
                    ),
                )
            if iiii2:
                EnsGene = iiii2.group(1)
            return EnsGene
        except Exception:
            return ""

    def Entrez(self):  # MetaCyc
        try:
            EntrezGene = re.findall(
                self.gene + r"_HUMAN[\S\s].*?\n", open(self.db).read()
            )[0].split("\t")[4]
            return EntrezGene
        except Exception:
            return ""

    def Uniprot(self):  # Uniprot
        try:
            UniProtGene = re.findall(
                self.gene + r"_HUMAN[\S\s].*?\n", open(self.db).read()
            )[0].split("\t")[0]
            return UniProtGene
        except Exception:
            return ""


class compound(object):
    def __init__(self, url, ident, time, EF, specialCompounds, *newparam):
        self.ident = ident
        bm = _bm()
        # Fetch data from REST API (ignoring legacy url parameter, always use REST)
        pagina_content = bm.getHtml(f"https://rest.kegg.jp/get/{self.ident}", time)
        self.pagina = (
            pagina_content.decode("utf-8")
            if isinstance(pagina_content, bytes)
            else pagina_content
        )
        # Use REST API parser (not HTML parser) since we're fetching from REST endpoint
        self.atributes = bm.getCompParamFromRestAPI(
            self.pagina, self.ident, time, EF, specialCompounds, RxnID=None
        )
        # store raw newparam then populate plain attributes
        self.newparam = newparam
        try:
            self._populate_attributes()
        except Exception:
            # be robust during parsing failures; attributes will default to empty
            pass

    @classmethod
    def from_batch_data(
        cls, ident, pagina_content, time, EF, specialCompounds, *newparam
    ):

        obj = cls.__new__(cls)  # Create instance without calling __init__
        obj.ident = ident
        obj.pagina = (
            pagina_content if isinstance(pagina_content, str) else str(pagina_content)
        )

        # Detect if this is KEGG flat file format (from REST API) or HTML format
        # Flat file format starts with "ENTRY" field
        is_flat_file = obj.pagina.strip().startswith("ENTRY")

        if is_flat_file:
            # Use new REST API parser
            bm = _bm()

            obj.atributes = bm.getCompParamFromRestAPI(
                obj.pagina, ident, time, EF, specialCompounds, RxnID=None
            )
        else:
            # Use legacy HTML parser
            bm = _bm()
            obj.atributes = bm.getCompParam(
                obj.pagina, ident, time, EF, specialCompounds, RxnID=None
            )

        obj.newparam = newparam
        # populate plain attributes for objects created via from_batch_data
        try:
            obj._populate_attributes()
        except Exception:
            pass

        return obj

    def _populate_attributes(self):
        """Populate plain attributes from self.atributes and self.newparam."""
        bm = _bm()

        # Helper function to safely get nested attributes
        def safe_get(data, *indices, default=""):
            try:
                result = data
                for idx in indices:
                    result = result[idx]
                return result
            except (IndexError, KeyError, TypeError):
                return default

        # Set all attributes with defaults
        self.ID1 = safe_get(self.atributes, 0, 0, 0)
        self.ID2 = safe_get(self.atributes, 0, 0, 1)

        self.Formula1 = safe_get(self.atributes, 1, 0, 0)
        self.Formula2 = safe_get(self.atributes, 1, 0, 1)
        self.Formula3 = self.Formula1  # Legacy mirror
        self.Formula4 = safe_get(self.atributes, 11)

        self.AssRxn1 = safe_get(self.atributes, 2, 0)
        try:
            self.AssRxn2 = (
                [x[0] for x in self.atributes[2][0]] if self.atributes[2][0] else ""
            )
            self.AssRxn3 = (
                [x[1] for x in self.atributes[2][0]] if self.atributes[2][0] else ""
            )
        except (IndexError, KeyError, TypeError):
            self.AssRxn2 = ""
            self.AssRxn3 = ""

        self.Name = safe_get(self.atributes, 3, 0)
        self.Subcel = self.newparam

        # External database identifiers
        self.PubChem = safe_get(self.atributes, 4, 0)
        self.CheBI = safe_get(self.atributes, 5, 0)
        self.LIPIDMAPS = safe_get(self.atributes, 6, 0)
        self.LipidBank = safe_get(self.atributes, 7, 0)
        self.GlyDB = safe_get(self.atributes, 8, 0)
        self.JCGGDB = safe_get(self.atributes, 9, 0)

        # Chemical properties
        self.charge = safe_get(self.atributes, 10)
        self.inchikey = safe_get(self.atributes, 12)
        self.inchi = safe_get(self.atributes, 13)
        self.CID = safe_get(self.atributes, 14)

        # Atom compositions
        self.Atom1 = ""
        self.Atom2 = ""
        self.Atom3 = ""

        try:
            id1 = safe_get(self.atributes, 0, 0, 0)
            if id1 and len(id1) > 0:
                if id1[0] == "C":
                    self.Atom1 = bm.atom(self.Formula1) if self.Formula1 else ""
                elif id1[0] == "G":
                    self.Atom1 = bm.glycan(self.Formula1) if self.Formula1 else ""

            id2 = safe_get(self.atributes, 0, 0, 1)
            if id2 and len(id2) > 0:
                if id2[0] == "C":
                    self.Atom2 = bm.atom(self.Formula2) if self.Formula2 else ""
                elif id2[0] == "G":
                    self.Atom2 = bm.glycan(self.Formula2) if self.Formula2 else ""

            self.Atom3 = self.Atom1  # Legacy mirror
        except Exception:
            pass
