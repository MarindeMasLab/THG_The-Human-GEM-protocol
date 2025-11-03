#!/usr/bin/python
# -*- coding: utf-8 -*-

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
    def __init__(self, ec, time):
        self.ec = ec
        session = setup_biocyc_session()
        self.GPRPAss = getGPR(self.ec, session)
        location_module = _loc()
        self.Subcell = location_module.getLocationnew(
            self.GPRPAss[3],
            self.GPRPAss[1],
            self.GPRPAss[2],
            1,
            "files/bb.pickle",
            session,
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
    def __init__(self, gene, db):
        self.gene = gene
        self.db = db

    def Name(self):  # Patway-KEGG
        try:
            return self.gene
        except Exception:
            return ""

    def Ensg(self):  # MetaCyc
        try:
            iiii2 = re.search(
                "gene=([A-Z0-9]+)",
                str(
                    urllib.request.urlopen(
                        "https://www.genome.jp/dbget-bin/www_bget?hsa+" + self.gene
                    ).read()
                ),
            )  # .group(1)
            if not iiii2:
                iiii2 = re.search(
                    "(ENSG[0-9]+)",
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
                self.gene + "_HUMAN[\S\s].*?\n", open(self.db).read()
            )[0].split("\t")[4]
            return EntrezGene
        except Exception:
            return ""

    def Uniprot(self):  # Uniprot
        try:
            UniProtGene = re.findall(
                self.gene + "_HUMAN[\S\s].*?\n", open(self.db).read()
            )[0].split("\t")[0]
            return UniProtGene
        except Exception:
            return ""


class compound(object):
    def __init__(self, url, ident, time, EF, specialCompounds, *newparam):
        self.ident = ident
        bm = _bm()
        pagina_content = bm.getHtml(f"https://rest.kegg.jp/get/{self.ident}", time)
        self.pagina = (
            pagina_content.decode("utf-8")
            if isinstance(pagina_content, bytes)
            else pagina_content
        )
        self.atributes = bm.getCompParam(
            self.pagina, self.ident, time, EF, specialCompounds, RxnID=None
        )
        self.newparam = newparam

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
        return obj

    def ID1(self):
        try:
            return self.atributes[0][0][0]
        # 			return getID(url)
        except Exception:
            return ""

    def ID2(self):
        try:
            return self.atributes[0][0][1]
        # 			return getID(url)
        except Exception:
            return ""

    def AssRxn1(self):
        try:
            return self.atributes[2][0]
        except Exception:
            return ""

    def AssRxn2(self):
        try:
            return [x[0] for x in self.atributes[2][0]]
        except Exception:
            return ""

    def AssRxn3(self):
        try:
            return [x[1] for x in self.atributes[2][0]]
        except Exception:
            return ""

    def Formula1(self):
        try:
            return self.atributes[1][0][0]
        except Exception:
            return ""

    def Formula2(self):
        try:
            return self.atributes[1][0][1]
        except Exception:
            return ""

    def Formula3(self):
        try:
            return self.atributes[1][0][0]
        except Exception:
            return ""

    def Atom1(self):
        try:
            bm = _bm()
            if self.atributes[0][0][0][0] == "C":
                composition = bm.atom(self.atributes[1][0][0])
            elif self.atributes[0][0][0][0] == "G":
                composition = bm.glycan(self.atributes[1][0][0])
            return composition
        except Exception:
            return ""

    def Atom2(self):
        try:
            bm = _bm()
            if self.atributes[0][0][1][0] == "C":
                composition = bm.atom(self.atributes[1][0][1])
            elif self.atributes[0][0][1][0] == "G":
                composition = bm.glycan(self.atributes[1][0][1])
            return composition
        except Exception:
            return ""

    def Atom3(self):
        try:
            bm = _bm()
            if self.atributes[0][0][0][0] == "C":
                composition = bm.atom(self.atributes[1][0][0])
            elif self.atributes[0][0][0][0] == "G":
                composition = bm.glycan(self.atributes[1][0][0])
            return composition
        except Exception:
            return ""

    def Name(self):
        try:
            return self.atributes[3][0]
        except Exception:
            return ""

    def Subcel(self):
        try:
            return self.newparam
        except Exception:
            return ""

    def PubChem(self):
        try:
            return self.atributes[4][0]
        except Exception:
            return ""

    def CheBI(self):
        try:
            return self.atributes[5][0]
        except Exception:
            return ""

    def LIPIDMAPS(self):
        try:
            return self.atributes[6][0]
        except Exception:
            return ""

    def LipidBank(self):
        try:
            return self.atributes[7][0]
        except Exception:
            return ""

    def GlyDB(self):
        try:
            return self.atributes[8][0]
        except Exception:
            return ""

    def JCGGDB(self):
        try:
            return self.atributes[9][0]
        except Exception:
            return ""

    def charge(self):
        try:
            return self.atributes[10][:]
        except Exception:
            return ""

    def Formula4(self):
        try:
            return self.atributes[11][:]
        except Exception:
            return ""

    def inchikey(self):
        try:
            return self.atributes[12][:]
        except Exception:
            return ""

    def inchi(self):
        try:
            return self.atributes[13][:]
        except Exception:
            return ""

    def CID(self):
        try:
            return self.atributes[14][:]
        except Exception as e:
            return ""
