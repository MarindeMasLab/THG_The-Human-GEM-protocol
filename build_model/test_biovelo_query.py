#!/usr/bin/env python3
"""
Quick test script to run a BioVelo xmlquery against the BioCyc web service
and extract the <string> values returned (e.g. subcellular-location).

Example (default):
    python build_model/test_biovelo_query.py --org HUMAN --gene TP53

This will issue a query equivalent to:
    [(g^product^subcellular-location) : g <- HUMAN^^genes, g^name = "TP53"]

The script URL-encodes the query and prints a short summary plus the
extracted string elements. Use --raw to show the raw XML response.
"""
import argparse
import sys
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path

# Add parent directory to path to import functions
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from functions.gpr.auth_gpr import setup_biocyc_session


def build_protein_query(org: str, gene: str) -> str:
    # Step 1: Get proteins associated with genes using the special function gene-to-proteins
    # Example: [gene-to-proteins(g) : g <- ECOLI^^genes, g^name = "fbaA"]
    q = f'[gene-to-proteins(g) : g <- {org}^^genes, g^name = "{gene}"]'
    return q


def build_protein_getxml_url(
    org: str, protein_frameid: str, detail: str = "full"
) -> str:
    # Step 2: Get protein object directly using getxml endpoint
    # Example: https://websvc.biocyc.org/getxml?ECOLI:FRUCTBISALD-CLASSII-MONOMER&detail=full
    url = f"https://websvc.biocyc.org/getxml?{org}:{protein_frameid}&detail={detail}"
    return url


def run_query(session, query: str, detail: str = "low", timeout: int = 30) -> str:
    base_url = "https://websvc.biocyc.org/xmlquery"
    encoded = urllib.parse.quote(query, safe="")
    url = f"{base_url}?query={encoded}&detail={detail}"
    resp = session.get(url, timeout=timeout)
    resp.raise_for_status()
    return resp.text


def extract_strings_from_xml(xml_text: str):
    """Return list of <string> element text values (anywhere in document).
    Falls back to empty list if parsing fails.
    """
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return []
    # collect text from any <string> tag
    strings = [el.text for el in root.findall(".//string") if el.text]
    return strings


def extract_genes_from_xml(xml_text: str):
    """Return list of gene dicts found in xmlquery gene responses.
    Each dict contains 'id' (e.g. ECOLI:EG10282), 'frameid', 'common_name', and a list of product frameids.
    """
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return []
    genes = []
    for g in root.findall(".//Gene"):
        gid = g.get("ID")
        frameid = g.get("frameid")
        common = None
        cn = g.find("common-name")
        if cn is not None and cn.text:
            common = cn.text
        # find product protein frameids
        products = []
        for prod in g.findall(".//product"):
            p = prod.find("Protein")
            if p is None:
                # sometimes product may be another tag; try any child
                for child in prod:
                    if child.tag and child.get("frameid"):
                        products.append(child.get("frameid"))
            else:
                if p.get("frameid"):
                    products.append(p.get("frameid"))
        genes.append(
            {"id": gid, "frameid": frameid, "common_name": common, "products": products}
        )
    return genes


def find_cco_resources(xml_text: str):
    """Return list of resource attribute values from <cco resource="..."/> tags inside the XML."""
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return []
    resources = []
    for cco in root.findall(".//cco"):
        res = cco.get("resource")
        if res:
            resources.append(res)
    return resources


def extract_locations_from_xml(xml_text: str):
    """Extract subcellular location information from <location> elements containing <cco> tags.
    Returns a list of dicts with 'resource', 'frameid', and 'orgid' for each location.
    """
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return []

    locations = []
    # Find all <location> elements
    for loc_elem in root.findall(".//location"):
        # Find <cco> child elements
        for cco in loc_elem.findall(".//cco"):
            location_info = {
                "resource": cco.get("resource"),
                "frameid": cco.get("frameid"),
                "orgid": cco.get("orgid"),
            }
            locations.append(location_info)

    return locations


def get_compartment_name(session, org: str, frameid: str) -> str:
    """Resolve a CCO/CCI frameid to its common-name using getxml.

    Args:
        session: Authenticated BioCyc session
        org: Organism ID (e.g., 'ECOLI', 'HUMAN')
        frameid: The CCO/CCI frameid (e.g., 'CCI-CYTOSOL-GN')

    Returns:
        The common-name of the compartment, or the frameid if not found
    """
    try:
        url = f"https://websvc.biocyc.org/getxml?{org}:{frameid}&detail=low"
        response = session.get(url, timeout=10)
        response.raise_for_status()
        xml_text = response.text

        # Parse and extract common-name
        root = ET.fromstring(xml_text)
        # Look for <common-name> element
        common_name_elem = root.find(".//common-name")
        if common_name_elem is not None and common_name_elem.text:
            return common_name_elem.text

        # Fallback: try to find it in string elements
        strings = [el.text for el in root.findall(".//string") if el.text]
        if strings:
            return strings[0]

        return frameid  # Return frameid if we can't resolve it
    except Exception:
        return frameid  # Return frameid on any error


def extract_proteins_from_xml(xml_text: str):
    """Extract protein information from query results.
    Returns a list of dicts with 'id', 'frameid', and 'common_name' for each protein.
    """
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return []

    proteins = []
    for prot in root.findall(".//Protein"):
        protein_info = {
            "id": prot.get("ID"),
            "frameid": prot.get("frameid"),
            "common_name": None,
        }
        # Try to get common name
        cn = prot.find("common-name")
        if cn is not None and cn.text:
            protein_info["common_name"] = cn.text

        proteins.append(protein_info)

    return proteins


def fetch_resource_and_extract(session, url_or_resource: str):
    """Fetch a resource (either full URL or relative 'getxml?...' resource) and extract any <string> values."""
    if url_or_resource.startswith("http"):
        url = url_or_resource
    else:
        # assume resource like getxml?ECOLI:FRAMEID or similar
        url = urllib.parse.urljoin("https://biocyc.org/", url_or_resource)
    try:
        r = session.get(url, timeout=30)
        r.raise_for_status()
        txt = r.text
    except Exception:
        return []
    # try to extract strings from returned XML
    vals = extract_strings_from_xml(txt)
    return vals


def main():
    p = argparse.ArgumentParser(
        description="Test BioVelo xmlquery for subcellular-location"
    )
    p.add_argument(
        "--org", default="HUMAN", help="organism DB id (e.g. ECOLI, HUMAN, META)"
    )
    p.add_argument(
        "--gene", default="TP53", help="gene name (exact match against g^name)"
    )
    p.add_argument(
        "--detail",
        default="low",
        choices=["none", "low", "full"],
        help="detail level for xmlquery",
    )
    p.add_argument("--raw", action="store_true", help="print raw XML response")
    args = p.parse_args()

    # Set up authenticated BioCyc session
    print("Setting up BioCyc session...")
    try:
        session = setup_biocyc_session()
        print("BioCyc session established.")
    except Exception as e:
        print(f"Failed to set up BioCyc session: {e}", file=sys.stderr)
        print(
            "Make sure BIOCYC_EMAIL and BIOCYC_PASSWORD are set in environment or .env file"
        )
        sys.exit(1)

    # Step 1: Query for proteins associated with the gene
    protein_query = build_protein_query(args.org, args.gene)
    print(f"\nStep 1 - Protein Query: {protein_query}")
    print("Retrieving proteins associated with gene...")
    try:
        protein_xml = run_query(session, protein_query, detail=args.detail)
    except Exception as e:
        print(f"Protein query failed: {e}", file=sys.stderr)
        sys.exit(2)

    if args.raw:
        print("--- RAW PROTEIN XML ---")
        print(protein_xml)

    # Extract proteins from the response
    proteins = extract_proteins_from_xml(protein_xml)
    if not proteins:
        print("No proteins found using proteins-of-gene function.")
        print("Response snippet:", protein_xml[:1000].replace("\n", " "))
        # Fall through to the existing fallback logic below
    else:
        print(f"\nFound {len(proteins)} protein(s):")
        for prot in proteins:
            print(
                f" - ID: {prot['id']}, frameid: {prot['frameid']}, common-name: {prot['common_name']}"
            )

        # Step 2: For each protein, retrieve it directly using getxml
        print("\nStep 2 - Retrieving protein objects using getxml...")
        for prot in proteins:
            frameid = prot["frameid"]
            print(f"\n  Protein: {frameid}")
            protein_url = build_protein_getxml_url(args.org, frameid, detail="full")
            print(f"  URL: {protein_url}")

            try:
                protein_xml = session.get(protein_url, timeout=30).text
            except Exception as e:
                print(f"  Protein retrieval failed: {e}")
                continue

            if args.raw:
                print("  --- RAW PROTEIN XML ---")
                print(protein_xml[:2000])  # Show first 2000 chars

            # Extract location information
            locations = extract_locations_from_xml(protein_xml)
            if locations:
                print(f"  Found {len(locations)} location(s):")
                for loc in locations:
                    frameid = loc["frameid"]
                    # Resolve frameid to common name
                    common_name = get_compartment_name(session, args.org, frameid)
                    print(f"    - {common_name}")
                    print(f"      (frameid: {frameid}, orgid: {loc['orgid']})")
            else:
                # Try extracting strings as fallback
                strings = extract_strings_from_xml(protein_xml)
                if strings:
                    print(f"  Found location string(s): {', '.join(strings)}")
                else:
                    print("  No locations found.")
                    print("  Checking for location elements in XML...")
                    # Show snippet to debug
                    if "location" in protein_xml.lower():
                        print("  XML contains 'location' - showing relevant snippet:")
                        for line in protein_xml.split("\n"):
                            if "location" in line.lower():
                                print(f"    {line.strip()[:200]}")
        return

    # No direct strings found. Try fallback: relaxed gene search (instringci),
    # then extract product frameids and query proteins for subcellular-location.
    print("No <string> elements found in direct query. Trying fallback gene search...")
    try:
        gene_xml = run_query(
            session,
            f'[g : g<- {args.org}^^genes, "{args.gene}" instringci g^name]',
            detail=args.detail,
        )
    except Exception as e:
        print(f"Relaxed gene search failed: {e}", file=sys.stderr)
        sys.exit(2)

    genes = extract_genes_from_xml(gene_xml)
    if not genes:
        print(
            "Relaxed gene search returned no genes. Response snippet:",
            gene_xml[:1000].replace("\n", " "),
        )
        sys.exit(1)

    print(
        f"Found {len(genes)} gene(s). Exploring product objects and querying protein locations..."
    )
    for gene in genes:
        print(
            "\nGene:",
            gene.get("id"),
            "frameid:",
            gene.get("frameid"),
            "common-name:",
            gene.get("common_name"),
        )
        prods = gene.get("products") or []
        if not prods:
            print(" - No product objects listed for this gene in the gene XML")
            continue
        for pf in prods:
            # construct full id like ORG:FRAMEID
            frameid = pf
            print(f" - Product frameid: {frameid}")

            # Retrieve protein using getxml
            protein_url = build_protein_getxml_url(args.org, frameid, detail="full")
            print(f"   URL: {protein_url}")
            try:
                prot_xml = session.get(protein_url, timeout=30).text
            except Exception as e:
                print(f"   Retrieval failed for product {frameid}: {e}")
                continue

            if args.raw:
                print("   --- RAW PROTEIN XML ---")
                print(prot_xml[:1000])

            # Extract locations
            locs = extract_locations_from_xml(prot_xml)
            if locs:
                print(f"   Found {len(locs)} location(s):")
                for loc in locs:
                    frameid = loc["frameid"]
                    # Resolve frameid to common name
                    common_name = get_compartment_name(session, args.org, frameid)
                    print(f"     - {common_name}")
                    print(f"       (frameid: {frameid}, orgid: {loc['orgid']})")
            else:
                # Try extracting strings as fallback
                strings = extract_strings_from_xml(prot_xml)
                if strings:
                    print(f"   Location string(s): {', '.join(strings)}")
                else:
                    print("   No locations found.")
                    if "location" in prot_xml.lower():
                        print("   XML contains 'location' - showing snippet:")
                        for line in prot_xml.split("\n"):
                            if "location" in line.lower():
                                print(f"     {line.strip()[:150]}")


if __name__ == "__main__":
    main()
