import logging
import re
import socket
from io import StringIO
from typing import List, Optional, Tuple
from collections import defaultdict

import pandas as pd
import requests
import pdb
from functions.gpr.ast_gpr import sanitize_gpr
from functions.gpr.auth_gpr import setup_biocyc_session

# Set default timeout for all socket operations (including urllib)
# This prevents urllib.request.urlopen from hanging indefinitely
socket.setdefaulttimeout(30)

LOGGER = logging.getLogger(__name__)
GPRURL2 = "http://www.genome.jp/dbget-bin/www_bget?ec:"
REMOVE_PAT = re.compile(r"\)|\(|'| |\]|\[")
TRANSFERRED_PAT = re.compile(
    r'It is now listed as..\n+.+EC\-([0-9\.\-]+)" class\="EC\-NUMBER"', re.MULTILINE
)


def get_ecnumber_biocyc_html(
    ec_number: str, session: requests.Session, org: str = "META"
) -> str:
    """Get raw HTML representing a EC entry in BioCyc."""
    return session.get(
        f"https://websvc.biocyc.org/{org}/NEW-IMAGE?type=EC-NUMBER&object=EC-{ec_number}"
    ).text


def get_html(
    request_url: str, session: Optional[requests.Session] = None, timeout: int = 30
) -> str:
    """Fetch an html by perfoming a GET HTTPS request, maybe with session.

    Args:
        request_url: URL to fetch
        session: Optional requests session to use
        timeout: Timeout in seconds (default: 30)
    """
    if session is not None:
        return session.get(request_url, timeout=timeout).text
    else:
        return requests.get(request_url, timeout=timeout).text


def pattern_match_org(page: str, org: str = "Homo Sapiens") -> List[str]:
    return sorted(
        [
            x[0]
            for x in [
                re.findall(
                    r"\<b\>Gene:\</b\>..([a-zA-Z0-9]+)....([A-Za-z0-9:,-_\/]+)<br>.*"
                    + org,
                    str(x),
                )
                for x in page.splitlines()
            ]
            if x
        ]
    )


def pattern_match_org_3(page: str, org: str = "Homo Sapiens") -> List[Tuple[str, str]]:
    pattern = r"<b>Gene:</b>\s+([a-zA-Z0-9]+)\s+([A-Za-z0-9:,-_\/]+)<br>"
    matches = re.findall(pattern, page)

    result = [(gene.strip(), id.strip()) for gene, id in matches]

    return sorted(result)


def match_biocyc_page(page: str, humancyc: bool = False) -> List[str]:
    if humancyc:
        # since HumanCyc relates to human,
        # we don't need to look for a particular organism
        return pattern_match_org_3(page, "")
    # fetch in human
    urls0 = pattern_match_org_3(page, "Homo Sapiens")
    # if failed, failed in mouse
    if not urls0:
        urls0 = pattern_match_org_3(page, "Mus Musculus")
    return urls0


def getGPR(
    ec_number: str, session: Optional[requests.Session]
) -> Optional[Tuple[List[str], str, str, str, str]]:
    """Retrieve GPR given EC-number for BioCyc."""
    if session is None:
        session = setup_biocyc_session()
    page = get_ecnumber_biocyc_html(ec_number, session, org="HUMAN")
    LOGGER.info(f"Page length: {len( page )}")

    urls0 = match_biocyc_page(page, humancyc=True)
    LOGGER.info(f"Urls matched from matched from HUMAN: {len(urls0)}")
    # if failed, maybe the EC number was changed

    if not urls0:
        ec_number_match = TRANSFERRED_PAT.findall(page)
        if ec_number_match:
            LOGGER.info(f"Trying transferred EC number: {ec_number_match[0]}")
            gprs = getGPR(ec_number_match[0], session)
            if gprs is not None:
                return gprs
        else:
            pass
    # if failed, try MetaCyc
    if not urls0:
        page = get_ecnumber_biocyc_html(ec_number, session)
        urls0 = match_biocyc_page(page)
        LOGGER.info(f"Urls matched from changed MetaCyc: {len(urls0)}")
    if not urls0:

        LOGGER.warn("Trying Kegg")
        # return fetch_kegg_rest(ec_number)
        parsed = _fetch_kegg_from_ec_html(ec_number)
    else:
        parsed = parseGPRnewest(ec_number, urls0, page, session)
    if parsed[0]:
        a, b, c, d, gpr, gpr_dict = parsed
        gpr = sanitize_gpr(gpr)
        gpr = "([" + re.sub(r"([A-Z0-9\-\.]+)", r"([\1])", gpr) + "])"

        parsed = a, b, c, d, gpr, gpr_dict
    return parsed


def _fetch_kegg_from_ec_html(ec_number: str):
    # GPRPage2 = get_html(str(GPRURL2) + str(ec_number))
    GPRPage2 = get_html(str(GPRURL2) + str(ec_number))
    LOGGER.info(f"Page length from kegg: {len(GPRPage2)}")
    urls0 = sorted(
        set(
            REMOVE_PAT.sub(
                "",
                str(
                    re.findall(
                        r"(\([A-Za-z0-9]+\))",
                        str(
                            re.findall(
                                r"hsa:............................................",
                                GPRPage2,
                            )
                        ),
                    )
                ),
            ).split(",")
        )
    )
    LOGGER.info(f"Urls from HUMAN kegg: {urls0}")
    if not urls0[0]:  # if not for homo sapiens try with mus musculus
        urls0 = sorted(
            set(
                REMOVE_PAT.sub(
                    "",
                    str(
                        re.findall(
                            r"(\([A-Za-z0-9]+\))",
                            str(
                                re.findall(
                                    r"mmu:............................................",
                                    GPRPage2,
                                )
                            ),
                        )
                    ),
                ).split(",")
            )
        )
        LOGGER.warn(f"Urls from MOUSE kegg: {urls0}")
    if urls0[0]:
        urls1 = [x[0:] for x in urls0]
        urls2 = urls1
        urls3 = (
            "[(["
            + str(urls1)
            .replace("[", "")
            .replace("]", "")
            .replace(", ", "*1]) or ([")
            .replace("'", "")
            + "*1])]"
        )
        urls4 = re.sub(r"*[0-9]+", "", str(urls3))
    else:
        urls1 = ""
        urls2 = ""
        urls3 = ""
        urls4 = ""
    if not urls1:
        return ([], "", "", "", "")
    return urls0, urls1, urls2, urls3, urls4


def fetch_kegg_rest(ec_number: str) -> Optional[Tuple[List[str], str, str, str, str]]:
    """Fetch ec-number from kegg and link it to its genes.

    Perform the following requests:

    .. code-block::

        LINK ec-number -> KO (Kegg orthology)
        LINK KO (Kegg orthology) -> genes

    We map first to KO because LINK to genes only works in this way.

    There is also the possibility of linking from ec-number to HSA
    (human gene Kegg ids) the problem is that gene -> gene symbol
    requires additional parsing.
    """
    ec_to_ko = pd.read_csv(
        StringIO(get_html(f"https://rest.kegg.jp/link/ko/ec:{ec_number}")),
        sep="\t",
        names=["ec", "ko"],
    )
    LOGGER.info(f"{len(ec_to_ko)} Kegg orthologies fetched.")

    ko_to_genes = pd.read_csv(
        StringIO(
            get_html(
                f"https://rest.kegg.jp/link/genes/{'+'.join(ec_to_ko.ko.to_list())}"
            )
        ),
        sep="\t",
        names=["ko", "genes"],
    )

    genes = ko_to_genes.loc[
        ko_to_genes.genes.str.startswith("hsa:") & ~ko_to_genes.genes.isnull(), "genes"
    ].to_list()
    if not genes:
        return None
    LOGGER.info(f"{len(genes)} Kegg genes fetched.")
    genes_stripped = [gene.replace("hsa:", "") for gene in genes]
    return (
        genes,
        genes,
        genes_stripped,
        "[" + " or ".join(f"({genes_stripped}) *1") + "]",
        f"[{' or '.join(genes_stripped)}]",
    )


def multiple_replace(dict, text):
    """ "Multiple Replacement."""
    # Create a regular expression  from the dictionary keys
    regex = re.compile("(%s)" % "|".join(map(re.escape, dict.keys())))
    # For each match, look-up corresponding value in dictionary
    return regex.sub(lambda mo: dict[mo.string[mo.start() : mo.end()]], text)


def ParseNestedParen(string, level):
    """Generate strings contained in nested (), indexing i = level"""
    if len(re.findall(r"\[", string)) == len(re.findall(r"\]", string)):
        LeftRightIndex = [
            x
            for x in zip(
                [Left.start() + 1 for Left in re.finditer(r"\[", string)],
                reversed([Right.start() for Right in re.finditer(r"\]", string)]),
            )
        ]
    elif len(re.findall(r"\[", string)) > len(re.findall(r"\]", string)):
        return ParseNestedParen(string + "]", level)
    elif len(re.findall(r"\[", string)) < len(re.findall(r"\]", string)):
        return ParseNestedParen("[" + string, level)
    else:
        return "fail"
    return [string[LeftRightIndex[level][0] : LeftRightIndex[level][1]]]


def calculate_stoichiometry(input_str, gene_dict):
    gene_stoichiometry = defaultdict(int)

    def process_segment(segment, current_multiplier=1):
        matches = re.findall(r"(\[(?:[^\[\]]+|\[.*?\])*\]|\w+)\*?(\d*)", segment)
        for expr, multiplier in matches:
            multiplier = int(multiplier) if multiplier else 1
            if "[" in expr:  # Nested expression
                process_segment(expr[1:-1], current_multiplier * multiplier)
            else:  # Simple gene expression
                if expr in list(gene_dict.values()):  # Check if gene is in gene_dict
                    gene_stoichiometry[expr] += current_multiplier * multiplier

    process_segment(input_str)
    # Update gene_dict with calculated stoichiometry
    for gene in gene_stoichiometry:
        if gene not in list(gene_dict.values()):
            # If gene is not in gene_dict, remove it from gene_stoichiometry
            del gene_stoichiometry[gene]

    return gene_stoichiometry


def pattern_match_org_2(page: str, org: str = "Homo Sapiens") -> List[str]:
    org_escaped = re.escape(org)
    pattern = rf"<b>Gene:</b>\s*([a-zA-Z0-9]+)\s+([A-Za-z0-9]+)<br>"
    matched_genes = sorted(
        [match[0] for line in page.splitlines() for match in re.findall(pattern, line)]
    )

    return matched_genes


def extract_gene_pairs(gpr):
    # Extract gene pairs from gpr using regular expressions
    gene_pattern = re.compile(r"\[(\w+)\*\d+\*0\]")
    genes = gene_pattern.findall(gpr)
    pairs = [(genes[i], genes[i + 1]) for i in range(len(genes) - 1)]
    return pairs


def check_genes_in_d2(pair, d2):
    # Check if the gene pair appears in d2 separated by <br>
    gene1, gene2 = pair
    pattern = re.compile(rf"{gene1}.*?<br>.*?{gene2}", re.DOTALL)
    return bool(pattern.search(d2))


def replace_specific_and_with_or(gpr, pairs, d2):
    # Replace specific "and" with "or" in the gpr string
    for gene1, gene2 in pairs:
        if check_genes_in_d2((gene1, gene2), d2):
            gpr = re.sub(
                rf"\[{gene1}\*\d+\*0\] and \[{gene2}\*\d+\*0\]",
                rf"[{gene1}*48*0] or [{gene2}*12*0]",
                gpr,
            )
    return gpr


def extract_gene_pairs(gpr, gene_dict):
    # Use the gene dictionary to match gene names in the gpr string
    genes = list(gene_dict.values())
    # if the gene is not in the gpr string, remove it from the list
    genes = [gene for gene in genes if gene in gpr]
    # order these genes by the order they appear in the gpr string, if they appear
    genes = sorted(genes, key=lambda x: gpr.index(x))
    # Create a list of gene pairs
    pairs = [(genes[i], genes[i + 1]) for i in range(len(genes) - 1)]
    return pairs


def check_genes_in_d2(pair, d2):
    # Check if the gene pair appears in d2 separated by <br>
    gene1, gene2 = pair
    pattern = re.compile(rf"{gene1}.*?<br>.*?{gene2}", re.DOTALL)
    return bool(pattern.search(d2))


def replace_specific_and_with_or(gpr, pairs, d2):
    for gene1, gene2 in pairs:
        if check_genes_in_d2((gene1, gene2), d2):
            # Create a regex pattern to find "and" between gene1 and gene2 with wildcards around them
            pattern = re.compile(
                rf"{re.escape(gene1)}.*?and.*?{re.escape(gene2)}", re.DOTALL
            )
            # Find the match
            match = pattern.search(gpr)

            # Only replace if match was found
            if match:
                # Replace the "and" with "or" leaving the rest of the match unchanged
                gpr = (
                    gpr[: match.start()]
                    + gpr[match.start() : match.end()].replace("and", "or")
                    + gpr[match.end() :]
                )
    return gpr


def reorder_gpr_according_to_d2(gpr, d2, gene_dict):
    # Extract genes from d2 in the order they appear
    ordered_genes = [gene for gene in gene_dict.values() if gene in d2]
    ordered_genes = sorted(ordered_genes, key=lambda x: d2.index(x))

    # Extract the full gene entries from gpr based on the ordered genes
    gene_entries = []

    for gene in ordered_genes:
        # pattern should include the gene name, the square brackets just before and after the gene name but not the rest of genes
        pattern = re.compile(rf"\[{gene}.*?\]")

        # then the match will be only that relevant gene not the rest of them
        match = pattern.search(gpr)
        # out of the match we only want the current gene so keep everything after the last and in the match string
        if match:
            try:
                match = match.group(0).split("and")[-1]
            except:
                pass
        if match:
            gene_entries.append(match)

    # Join the gene entries with ' and ' to form the new gpr string
    new_gpr = " and ".join(gene_entries)

    return new_gpr


def fix_coefficients_in_gpr(gpr):

    # separate the gene names to go through them one by one, the gene names are separated by 'and' or 'or'
    genes = re.split(r"and|or", gpr)

    # if only one gene name is present, return the gpr as it is
    if len(genes) > 1:
        # detect by what the gene name is separated, 'and' or 'or' to later join them back together
        separator = re.search(r"and|or", gpr).group(0)

    # go through each gene name and fix the coefficients
    for i, gene in enumerate(genes):
        # if the gene name has a coefficient outside of the square brackets

        if "*" in gene:
            # consider there might be more than one coefficient in the gene name, separated by '*'
            gene_parts = gene.split("*")
            # count the number of closing square brackets in gene
            count = gene_parts[0].count("]")
            # join the actual gene name (no square brackets) with all the coefficients of that gene
            # in the gene name only add the letters, no square brackets or coefficients
            # check that the gene name is only composed of letters and remove any other characters
            gene = "".join([char for char in gene_parts[0] if char != "]"])

            gene += "*" + "*".join(gene_parts[1:])
            # now at the end of the gene name add the number of closing square brackets that are missing
            gene += "]" * count
            # replace the gene name in the list of genes with the fixed gene name
            genes[i] = gene

    # if there was more than one gene name, join them back together with the same separator that separated them before
    if len(genes) > 1:
        new_gpr = f" {separator} ".join(genes)
    else:
        new_gpr = genes[0]

    return new_gpr


def update_urls3_with_coefficients(urls3, w):
    # Create a dictionary from the list 'w' with gene name as key and coefficient as value
    gene_dict = {
        item.split("*")[0]: item.split("*")[1] for sublist in w for item in sublist
    }

    # Function to replace the coefficients in the urls3 string
    def replace_coefficients(match):
        gene = match.group(1)
        if gene in gene_dict:
            # Replace the coefficients with the corresponding value from gene_dict
            return f"{gene}*{gene_dict[gene]}"
        else:
            return match.group(0)

    # Regex pattern to match the gene names and coefficients
    pattern = re.compile(r"(\w+)\*\d+\*\d+")

    # Update the urls3 string
    updated_urls3 = pattern.sub(replace_coefficients, urls3)

    return updated_urls3


# New function (August 2024):


def parseGPRnewest(
    ec_number, urls0: List[str], page: str, session: requests.Session
) -> Tuple[List[str], str, str, str, str]:
    # Attach a short, per-invocation diagnostic id so we can trace repeated
    # calls or duplicated log lines more easily when running in threads.
    try:
        import threading

        _call_thread = threading.get_ident()
    except Exception:
        _call_thread = 0
    try:
        import uuid

        _call_uuid = uuid.uuid4().hex[:8]
    except Exception:
        _call_uuid = "-"
    call_id = f"{ec_number}:{_call_uuid}:{_call_thread}"
    LOGGER.debug(
        f"parseGPRnewest START for EC {ec_number}, urls0 length: {len(urls0)} call_id={call_id}"
    )
    urls1 = [
        i for i in reversed(sorted([x[0:][0] for x in urls0], key=len))
    ]  # Sort genes by name lenght
    LOGGER.debug(f"parseGPRnewest: sorted {len(urls1)} genes")
    Ls = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    urls2 = []
    # urls22 = []
    for x in urls1:
        for y in urls0:
            if re.findall(r"\('" + str(x) + "',", str(y)):
                urls2 = urls2 + [y[1]]
                # urls22 = urls22 + [(y[0], y[1], Ls.pop(0))]
    # LOGGER.debug(f"{urls22=}")
    a = page.replace("\n", "").replace("Enzymes and Genes:", "\nEnzymes and Genes:")
    a2 = re.findall(r"Enzymes and Genes:.*[\S\s]+", a)
    b = a2[0].replace("<br> <a href=", "\n<br> <a href=").replace("</a>)", "</a>)\n")
    c = re.findall(
        r"<br> <a href=.*[\S\s].*>Homo sapiens</i>", b
    )  # evaluate the ezymes active in human Extracted from Biocyc# c = re.findall(r'<br> <a href=.*[\S\s].*>Homo sapiens</a>',b)
    if not c:
        c = re.findall(
            r'href="/gene\?orgid=HUMAN&id=(.*?)</a>', b
        )  # evaluate the ezymes active in human extracted from HumanCyc# c = re.findall(r'<br> <a href=.*[\S\s].*>Mus musculus</a>',b)
    if not c:
        c = re.findall(
            r"<br> <a href=.*[\S\s].*>Mus musculus</i>", b
        )  # evaluate the ezymes active in human # c = re.findall(r'<br> <a href=.*[\S\s].*>Mus musculus</a>',b)

    dict1 = {
        "<SUB>": "*",
        "<@SUB>(": " and (",
        "<@SUB>[": " and [",
        "][": "] and [",
        ")(": " and ",
        "'": "",
        ",": "",
        ";": "",
        "<BR>": "",
    }
    dict2 = {
        ", ": "*",
        "), (": " and ",
        "(": "",
        ")": "",
        ",": "",
        ";": "",
        "'": "",
        "<b>": "",
        "</b>": "",
        "<br>": "",
        "</br>": "",
        "<i>": "",
        "</i>": "",
        "<a>": "",
        "</a>": "",
        ":": "",
        "=": ",",
    }
    gpr2 = []
    gpr2list = []

    # Create a dictionary to store the gene IDs as keys and the gene names as values

    gene_dict = {}

    for entry in c:
        gene_id_match = re.search(r"(HS\d+)", entry)
        gene_name_match = re.search(r"<b>Gene:</b>\s*(\w+)", entry)
        if gene_id_match and gene_name_match:
            gene_id = gene_id_match.group(1)
            gene_name = gene_name_match.group(1)
            gene_dict[gene_id] = gene_name

    # Next, for each of the ECs numbers, we will use them to define the variable d1, which is the page of the enzyme and from this variable we extract d2, which is the subunit composition of the enzyme. We will also extract the gene ID from the page and use it to get the gene name from the dictionary we created earlier.

    iters = 0
    # Create a dict to store the reaction id as key and the GPR as value
    rxn_gpr_dict = {}

    # From here :
    LOGGER.debug(
        f"parseGPRnewest: Starting to process {len(c)} isoforms call_id={call_id}"
    )
    while iters < len(c):  # isoforms
        LOGGER.debug(
            f"parseGPRnewest: Processing isoform {iters+1}/{len(c)} call_id={call_id}"
        )

        if re.findall(r'class="GENE"', c[iters]):

            c2 = c[iters].replace(" ", "\n")
            c3 = c2.replace('"\n', '" \n')
            isourl = re.findall(
                r"(http://biocyc.org/META/NEW-IMAGE\?type=ENZYME.*[\S\s])\" \n",
                c3.replace("/META", "http://biocyc.org/META"),
            )  # if it is a complex
            if not isourl:
                isourl = re.findall(
                    r"<br>\n<a\nhref=\"(http://biocyc.org/gene\?orgid.*)\"",
                    c3.replace("/gene?orgid", "http://biocyc.org/gene?orgid"),
                )  # if it is  not a complex
            # if not isourl: isourl = ['https://biocyc.org/gene?orgid=META&id='+re.findall('^(.*?)"',c3)[0].strip()] # if it is taken from HumanCyc
            if not isourl:
                if re.findall('^(.*?)"', c3):

                    isourl = [
                        "https://biocyc.org/gene?orgid=HUMAN&id="
                        + re.findall('^(.*?)"', c3)[0].strip()
                    ]

                else:
                    hs_match = re.findall("HS\w+", c3)
                    if hs_match:
                        isourl = [
                            "https://biocyc.org/gene?orgid=HUMAN&id="
                            + hs_match[0].strip()
                        ]
                    else:
                        # No valid gene ID found, skip this entry
                        LOGGER.debug("No valid gene ID found in c3, skipping isoform")
                        iters += 1
                        continue

            LOGGER.debug(f"parseGPRnewest: Fetching isopage from {isourl[0][:80]}...")
            isopage = str(
                get_html(isourl[0].replace(" ", "").replace('"', ""), session)
            )
            LOGGER.debug(f"parseGPRnewest: Got isopage, length={len(isopage)}")
            d1 = isopage.replace("\n", " ").replace(
                "</a>", "\n</a>"
            )  # d1 =  isopage.replace('\n',' ').replace('<br>','<br>\n')
            d2 = [
                x.strip() for x in re.findall(r"Subunit composition.*=(.*)\\n", d1)
            ]  # re.findall(r'</a> </p>  <p class=ecoparagraph>  Subunit composition .*=(.*)', d1.replace(" <br>",""))
            if not d2:
                d2 = re.findall(
                    r"</a> </p>  <p class=ecoparagraph>  Subunit Composition.*=(.*)",
                    d1.replace(" <br>", ""),
                )
            if not d2:
                d2 = re.findall(r"Subunit composition.*=(.*)", d1.replace(" <br>", ""))
                if d2:
                    d2 = d2[0]
            # if not d2:
            #    d2 = re.findall(r'Subunit Composition.*=(.*)', d1.replace(" <br>",""))
            #    if d2: d2 = d2[0]
            if d2:
                d2 = d2[
                    0
                ]  # if we found subunit composition then only use the first match
            if not d2:
                isourl2 = re.findall(r"/gene-tab.*META&tab=SUMMARY", str(isopage))
                if isourl2:
                    isourl2 = ["https://websvc.biocyc.org" + isourl2[0]]
                    isopage2 = str(get_html(isourl2[0], session))
                    d2 = re.findall(
                        r"Subunit Composition</td><td align=LEFT valign=TOP>(.*)",
                        isopage2.replace("<tr><td", "\n").replace("</td></tr>", "\n"),
                    )

            if not d2:
                d1 = str(isopage).replace("\n", " ").replace("</a>", "</a>\n")
                d2 = re.findall(
                    "/META/NEW-IMAGE\?type=REACTION.*\</a>",
                    d1.replace(" <br>", ""),
                )
            if not d2:
                # isourl = 'https://biocyc.org/gene?orgid=META&id='+re.findall('^(.*?)"',c3)[0].strip() # if it is taken from HumanCyc

                if re.findall('^(.*?)"', c3):

                    isourl = [
                        "https://biocyc.org/gene?orgid=HUMAN&id="
                        + re.findall('^(.*?)"', c3)[0].strip()
                    ]

                else:
                    hs_match = re.findall("HS\w+", c3)
                    if hs_match:
                        isourl = [
                            "https://biocyc.org/gene?orgid=HUMAN&id="
                            + hs_match[0].strip()
                        ]
                    else:
                        # No valid gene ID found, skip this entry
                        LOGGER.debug("No valid gene ID found in c3, skipping isoform")
                        iters += 1
                        continue
                isopage = str(
                    get_html(isourl[0].replace(" ", "").replace('"', ""), session)
                )

                d1 = str(isopage).replace("\n", " ").replace("</a>", "\n</a>")
                d2 = re.findall(
                    r">Subunit Composition</td><td.*?>(.*?)</td>", str(isopage)
                )
                if not d2:
                    ith_gene = re.findall(
                        r"gene\\n\\n<br><font class=\\\'header\\\'>(.*?)<",
                        isopage,
                    )
                    if ith_gene:
                        ith_gene = ith_gene[0].strip()
                        position = str(isopage).find("[" + ith_gene + "]<SUB>")
                        if position == -1:
                            d2 = ith_gene
                        else:
                            d2 = re.findall(
                                r"^(.*?</SUB>)",
                                isopage[position : position + 100],
                            )[0].strip()

            if not d2:

                d2 = []
                # Find all instances of the EC number in the HTML content​
                ec_pattern = re.escape(f"EC {ec_number}") + r"(?!\d)(.*?</a>)"
                ec_matches = re.findall(ec_pattern, d1, re.DOTALL)

                # Extract and Save Stoichiometry for each match​
                for match in ec_matches:
                    d2.append(match)
                    # Process d2 as needed​

            LOGGER.debug(f"{d2=}")

            # Extra condition if d2 is not found

            if not d2:

                # Check re.findall('^(.*?)"', c3) and if it is not empty, assign it to isourl
                if re.findall('^(.*?)"', c3):

                    isourl = [
                        "https://biocyc.org/gene?orgid=HUMAN&id="
                        + re.findall('^(.*?)"', c3)[0].strip()
                    ]

                else:

                    isourl = [
                        "https://biocyc.org/gene?orgid=HUMAN&id="
                        + re.findall("HS\w+", c3)[0].strip()
                    ]

                isopage = str(
                    get_html(isourl[0].replace(" ", "").replace('"', ""), session)
                )
                d1 = isopage.replace("\n", " ").replace("</a>", "\n</a>")

                d2 = []
                # Find all instances of the EC number in the HTML content​
                ec_pattern = re.escape(f"EC {ec_number}") + r"(?!\d)(.*?</a>)"
                ec_matches = re.findall(ec_pattern, d1, re.DOTALL)

                # Extract and Save Stoichiometry for each match​
                for match in ec_matches:
                    d2.append(match)

            LOGGER.debug(f"{d2=}")

            # Corrected regex pattern to capture text after Synonyms and before newline
            synonym_pattern = re.compile(r"<b>Synonyms:</b>\s*(.*?)\n", re.DOTALL)

            # Find and process synonyms
            synonym_matches = synonym_pattern.findall(d1)

            synonim = [
                multiple_replace(
                    dict2, multiple_replace(dict1, x.strip().replace("&nbsp;", ""))
                )
                for x in synonym_matches
            ]

            if not synonim and re.findall(r'top">\\nSynonyms</td>', d1):
                synonim = [
                    multiple_replace(dict2, multiple_replace(dict1, x))
                    for x in re.findall(r'Synonyms</td>.*top">(.*?)</td></tr>', d1)
                ]  # in case it comes from HumanCyc using getHtml2 function
                if synonim:
                    synonim = synonim[0].strip().split(";")
            if synonim:
                synonim = list(
                    filter(
                        None,
                        [re.sub(" $", "", re.sub("^ ", "", x)) for x in synonim],
                    )
                )

            # make d2 a list if it is not already a list
            if not isinstance(d2, list):
                d2 = [d2]

            if isinstance(d2, list):
                # call a variable d2_list to store the list of d2 in a way that the variable d2 can be used later on and not be modified
                d2_list = (
                    d2.copy()
                )  # any modifications to d2_list will not affect the original d2 object, and vice versa.

                d2_iter = -1

                for d2 in d2_list:

                    d2_iter += 1

                    # extract the Reaction ID from the d2 string
                    # Regular expression to find the reaction ID
                    rxn_id = re.search(r"ID:</b> (RXN-\d+|RXN\d+|RXN\d+-\d+)<br>", d2)

                    # Extracting the reaction ID
                    if rxn_id:
                        reaction_id = rxn_id.group(1)
                    else:
                        # build a reaction id based on the EC number and the d2 iteration number and the isoform iteration number
                        reaction_id = f"PUT_RXN_{ec_number}_{d2_iter}_{iters}"

                    d = list(
                        filter(
                            None,
                            re.findall(
                                "[\S]+",
                                str(d2)
                                .upper()
                                .replace("</SUB>", "<@SUB>")
                                .replace("/", " ")
                                .replace("(", "[")
                                .replace(")", "]"),
                            ),
                        )
                    )  # original without "list"
                    j = 0
                    e = [""] * len(d)

                    while j < len(d):
                        z = 0
                        while z < len(urls1):
                            IsVar = ""
                            if not re.findall(
                                "/NEW-IMAGE\?TYPE=REACTION", d[j].upper()
                            ):
                                if re.findall(urls1[z].upper(), d[j].upper()):
                                    IsVar = [d[j].upper()]
                                if not IsVar and synonim:
                                    dict3 = {}
                                    for x in synonim:
                                        dict3[x.upper()] = urls1[z].upper()
                                    IsVar = list(
                                        filter(
                                            None,
                                            [
                                                re.findall(
                                                    re.escape(x.upper()),
                                                    d[j].upper(),
                                                )
                                                for x in synonim
                                            ],
                                        )
                                    )  # original without list
                                    if IsVar:
                                        IsVar = [
                                            multiple_replace(dict3, x) for x in IsVar[0]
                                        ]
                                    del dict3
                                if not IsVar:
                                    IsVar = re.findall(urls1[z].upper(), d[j].upper())
                            if IsVar:
                                htry = (
                                    multiple_replace(dict1, IsVar[0])
                                    .replace("<@SUB>", "")
                                    .replace("] and [", "]*1 and [")
                                    .replace("] or [", "]*1 or [")
                                    .replace(")", "")
                                    .replace("(", "")
                                )

                                htry = re.sub("]$", "]*1", htry)

                                if not re.findall(r"]", htry) or not re.findall(
                                    r"\[", htry
                                ):
                                    htry = (
                                        "[" + str(htry) + "]*1"
                                    )  # to adapt the case of single gene reaction association

                                # If a ">" is found in htry, remove it and everything after it
                                if ">" in htry:
                                    htry = htry[: htry.index(">")]

                                gene_stoich = calculate_stoichiometry(htry, gene_dict)

                                if re.findall(
                                    r"Subunit composition.*=", d1
                                ) or re.findall(r"Subunit Composition:.*=", d1):
                                    z = len(
                                        urls1
                                    )  # if we are analyzing a complex, then if the gpr is found stop the iteration
                                    s = htry
                                else:
                                    NestLevel = re.findall(
                                        r"(?=(\]\*[0-9]+\]\*))", htry
                                    )
                                    # h0 = ParseNestedParen(str(htry), len(NestLevel))[0]
                                    h0 = htry
                                    Mult = re.findall(
                                        r"([\]\*0-9]+$)", htry
                                    )  # list of all sequences at the end of the string htry that consist of any combination of closing square brackets ], asterisks *, and digits 0-9

                                    if Mult and isinstance(Mult, list) and Mult[0]:

                                        parts = Mult[0].split("*")
                                        if len(parts) > 1:
                                            Mult = 1
                                            for part in parts:
                                                # Extract digits from Mult
                                                inner_part = re.findall(r"\d+", part)
                                                for inner in inner_part:
                                                    Mult *= int(str(inner))
                                        else:
                                            Mult = 1
                                            parts[0] = re.findall(r"\d+", parts[0])
                                            for part in parts[0]:
                                                Mult *= int(part)
                                    elif Mult:
                                        # Extract digits from Mult
                                        nums = re.findall(r"\d+", Mult)
                                        Mult = 1
                                        for part in nums:
                                            Mult *= int(part)
                                    else:
                                        Mult = 1

                                    htry = h0
                                    for x in urls1:
                                        htry = htry.replace(x, "@")

                                    htry = filter(None, htry.split("@"))
                                    dict4 = {
                                        "\\": "",
                                        '"': "",
                                        "'": "",
                                        "?": "",
                                        ")": "",
                                        "(": "",
                                        "]": "",
                                        "[": "",
                                        "*": "",
                                        "$": "",
                                        ".": "",
                                    }
                                    hh = []
                                    for x in htry:
                                        for w in urls1:
                                            for y in urls1:
                                                x1 = multiple_replace(dict4, x)
                                                m1 = str(y) + x1 + str(w)
                                                m2 = str(y) + x1
                                                if re.findall(
                                                    m1, multiple_replace(dict4, h0)
                                                ) or re.findall(
                                                    m2 + "$",
                                                    multiple_replace(dict4, h0),
                                                ):
                                                    t = re.findall(r"([0-9]+)", str(x))
                                                    if t:
                                                        t = str(
                                                            eval(
                                                                str(t)
                                                                .replace("'", "")
                                                                .replace('"', "")
                                                                .replace("]", "")
                                                                .replace("[", "")
                                                                .replace(", ", "*")
                                                            )
                                                        )
                                                    else:
                                                        t = "1"
                                                    hh = hh + [str(y) + "*" + t]

                                    hh = sorted(set(filter(None, hh)))
                                    # substitute the subunit composition in hh with the corresponding composition, stored in the gene_stoich dictionary, where the key is the gene name and the key is the stoichiometry
                                    for i, x in enumerate(hh):
                                        gene_name = x.split("*")[0]
                                        stoich = x.split("*")[1]
                                        if gene_name in gene_stoich:
                                            stoich = gene_stoich[gene_name]
                                        hh[i] = gene_name + "*" + str(stoich)

                                    if not hh:
                                        hh = [
                                            "["
                                            + str(multiple_replace(dict4, h0))
                                            + "]*1"
                                        ]  # no hh means that the string only contains the gene name # if not hh: hh = [str(h0)+"*1"]
                                    h22 = (
                                        str(hh).replace("'", "").replace(", ", " and ")
                                    )
                                    UniqGene = hh  # This line has been added and the next has been commented to express in a SGPR dimers, trimers, tetramers, ....
                                    # UniqGene = sorted(set([x.split("*")[0] for x in hh])) # if some of the subunits of the complex encoded by the same gene it is necessary to change"S and S" by 2*A # CHECK TO ADAPT TO THE GRASP CASE
                                    GeneCount = [
                                        (x.split("*")[0], x.split("*")[1])
                                        for x in h22.replace("[", "")
                                        .replace("]", "")
                                        .split(" and ")
                                    ]
                                    UniqCount = []
                                    for x in UniqGene:
                                        count = "0"
                                        for y in GeneCount:
                                            if re.findall(x, y[0]):
                                                count = count + "+" + y[1]
                                        UniqCount = UniqCount + [[x, str(eval(count))]]
                                    s = (
                                        "["
                                        + str(
                                            [
                                                str(x).replace("', '", "*")
                                                for x in UniqCount
                                            ]
                                        )
                                        .replace("'", "")
                                        .replace('"', "")
                                        .replace(", ", " and ")
                                        + "]"
                                    )
                                    z = len(urls1)

                            else:
                                s = ""
                                z += 1

                        e[j] = s  # multiple_replace(dict2, s)
                        j += 1
                    e = sorted(
                        set(filter(None, e))
                    )  # remove empty and duplicate elements
                    gpr = (
                        str(e)
                        .replace("'[", "")
                        .replace("]'", "")
                        .replace("'", "")
                        .replace(", ", " and ")
                    )

                    gpr = fix_coefficients_in_gpr(gpr)

                    gpr = reorder_gpr_according_to_d2(gpr, d2, gene_dict)

                    # Extract gene pairs from gpr
                    pairs = extract_gene_pairs(gpr, gene_dict)

                    # Replace specific "and" with "or" in gpr if condition is met
                    gpr = replace_specific_and_with_or(gpr, pairs, d2)

                    gpr2 = gpr2 + [
                        str(gpr).replace("'[", "").replace("]'", "").replace("'", "")
                    ]

                    gpr2list.extend(gpr2)

                    # dictionary with the reaction ID as key and the GPR as value

                    # check if the reaction ID is already in the dictionary
                    if reaction_id in rxn_gpr_dict:
                        # check if the rxn is already in the dictionary
                        if gpr not in rxn_gpr_dict[reaction_id]:

                            rxn_gpr_dict[reaction_id] = (
                                "[" + rxn_gpr_dict[reaction_id] + "] or [" + gpr + "]"
                            )
                    else:
                        rxn_gpr_dict[reaction_id] = gpr

        iters = iters + 1
        LOGGER.debug(
            f"parseGPRnewest: Completed isoform {iters}/{len(c)} call_id={call_id}"
        )

    urls3 = [
        str(sorted(set(gpr2list)))
        .replace("'[", "[")
        .replace("]'", "]")
        .replace(", ", "] or [")
        .replace("'", "")
        .replace("*0", "*1")
    ]

    # sort the dictionary by sorting its values
    rxn_gpr_dict = dict(sorted(rxn_gpr_dict.items(), key=lambda item: item[1]))

    # replace *0 with *1 in the GPRs dictionary
    # Function to replace *number*number*... with the product of those numbers using eval
    def replace_pattern(value):
        # This function will be used to replace the matched pattern
        def multiply_numbers(match):
            expression = match.group(1)
            product = eval(expression)
            return f"*{product}"

        # First replace *0 with *1
        value = value.replace("*0", "*1")

        # Regular expression to find patterns like *number*number*...
        return re.sub(r"\*(\d+(\*\d+)+)", multiply_numbers, value)

    # Update the dictionary values
    rxn_gpr_dict = {key: replace_pattern(value) for key, value in rxn_gpr_dict.items()}

    testIs = list(filter(None, [x.replace("[", "").replace("]", "") for x in urls3]))

    not_found = False

    if not urls3 or not testIs:
        urls3 = [str(x) + "*1" for x in urls1]

        # Save a record of the sGPR not found
        not_found = True

    g = []
    r = 0
    while r < len(urls3):
        GsprGenes = sorted(
            set(
                re.findall(
                    r"[A-Za-z0-9]+",
                    re.sub(
                        r"\*[0-9]+",
                        "",
                        str(urls3[r]).replace(" and ", " ").replace(" or ", " "),
                    ),
                )
            )
        )
        GsprGenes = [x.upper().replace("[", "").replace("]", "") for x in GsprGenes]
        intersection = len(GsprGenes) - len(
            set(GsprGenes).intersection([x.upper() for x in urls1])
        )
        if intersection < 1:
            g = g + [
                str(urls3[r].upper()).replace(" AND ", " and ").replace(" OR ", " or ")
            ]
        r += 1

    urls3 = (
        str(g)
        .replace("'', ", "")
        .replace("', ''", ")")
        .replace("', '", ") or (")
        .replace("'", "(")
        .replace("](]", "])]")
        .replace('"', "")
        .replace(" *", "*")
        .replace("[ ", "[")
        .replace(" ]", "]")
        .replace("(]", ")]")
        .replace("[)", "[(")
    )

    try:
        if re.findall(r"Subunit composition.*=", d1):
            urls3 = urls3
        else:
            w = []
    except NameError:
        w = []

    try:
        d1
    except NameError:
        d1 = ""

    if re.findall(r"Subunit composition.*=", d1):
        urls3 = urls3
    else:
        w = []

        for x in urls3.split(" or "):  # Evaluate Isoforms
            if re.findall(r"[A-Za-z0-9].*?", x):  # eliminate empty spots
                dict4 = {}
                for y in x.split(" and "):  # Evaluate complex
                    if re.findall(r"[A-Za-z0-9].*?", y):  # eliminate empty spots
                        # y2 = re.sub("(\(|\)|\]|\[)", "", y)
                        y2 = y
                        aa = re.findall(r"([A-Za-z0-9\-]+)\*([0-9\*]+)", y2)
                        LOGGER.debug(f"{aa=}")

                        if not aa[0][0] in dict4:
                            dict4[aa[0][0]] = eval(aa[0][1])
                        else:
                            dict4[aa[0][0]] = eval(
                                str(dict4[aa[0][0]]) + "+" + aa[0][1]
                            )
                        LOGGER.debug(f"{dict4=}")

                z = [(aa + "*" + str(dict4[aa])) for aa in dict4]

                w = w + [z]

    # Update urls3 with the coefficients from w
    urls3 = update_urls3_with_coefficients(urls3, w)

    urls4 = re.sub(r"\*[0-9]+", "", urls3)

    # If the sGPR was not found and a default one was set, add it to the dictionary so that the dictionary is not empty
    if not_found:
        rxn_gpr_dict["PUT_RXN_DEFAULT"] = urls3

    # Filter out key-value pairs where the key starts with 'PUT'
    put_entries = {
        key: value for key, value in rxn_gpr_dict.items() if key.startswith("PUT")
    }

    # Check if all 'PUT' values are the same
    if len(set(put_entries.values())) == 1:
        # Keep only one 'PUT' entry (choose the first one)
        key, value = next(iter(put_entries.items()))
        put_entries = {key: value}

    # Combine the reduced 'PUT' entries with the rest of the original dictionary
    rxn_gpr_dict = {
        **{
            key: value
            for key, value in rxn_gpr_dict.items()
            if not key.startswith("PUT")
        },
        **put_entries,
    }

    # If gpr has empty spots, remove them

    # Urls3:
    # Handle ors
    urls3 = urls3.replace("[] or ", "")
    urls3 = urls3.replace(" or []", "")
    # Handle ands
    urls3 = urls3.replace("[] and ", "")
    urls3 = urls3.replace(" and []", "")

    # Urls4:
    urls4 = urls4.replace("[] or ", "")
    urls4 = urls4.replace(" or []", "")
    # Handle ands
    urls4 = urls4.replace("[] and ", "")
    urls4 = urls4.replace(" and []", "")

    # rxn_gpr_dict: remove keys with empty values
    for key in list(rxn_gpr_dict.keys()):
        if not rxn_gpr_dict[key]:
            del rxn_gpr_dict[key]

    LOGGER.debug(f"parseGPRnewest: Finished processing, urls4 length={len(urls4)}")
    if not urls4:
        LOGGER.debug(f"parseGPRnewest: Returning empty result for EC {ec_number}")
        return ([], "", "", "", "", {})
    else:
        print("Full GPR: ", urls0, urls1, urls2, urls3, urls4, rxn_gpr_dict)
        LOGGER.debug(f"parseGPRnewest: Returning full result for EC {ec_number}")
        return urls0, urls1, urls2, urls3, urls4, rxn_gpr_dict
