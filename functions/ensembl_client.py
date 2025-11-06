import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests
from tqdm import tqdm
import re

LOGGER = logging.getLogger(__name__)


def _make_requests_session(retries: int = 3, backoff: float = 0.3) -> requests.Session:
    """Create a requests.Session with simple retry/backoff behaviour."""
    s = requests.Session()
    adapter = requests.adapters.HTTPAdapter(max_retries=retries)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    s.headers.update(
        {
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "THG-ensembl-client/1.0",
        }
    )
    return s


def fetch_ensembl_annotations(
    gene_identifiers, batch_size: int = 50, max_workers: int = 8, timeout: int = 5
):
    """Fetch Ensembl annotations for a set/list of gene identifiers (symbols or Ensembl IDs).

    Tries to look up genes by symbol first (via /lookup/symbol), then falls back
    to direct ID lookup if the input looks like Ensembl IDs. Then fetches
    cross-references (/xrefs/id) in parallel for IDs that were found.

    Parameters:
    -----------
    gene_identifiers : list
        List of gene symbols (e.g., "BRCA1", "TP53") or Ensembl IDs (e.g., "ENSG00000139618")

    Returns a dict mapping identifier -> annotation dict with keys:
      - ensembl: Ensembl gene id
      - display_name: gene symbol (if available)
      - biotype
      - description
      - entrez: list of Entrez IDs (may be empty)
      - uniprot: list of UniProt accessions (may be empty)
    """
    if not gene_identifiers:
        return {}

    session = _make_requests_session()
    ids_unique = list({str(e).strip() for e in gene_identifiers if e})

    LOGGER.info(
        f"Starting Ensembl annotation fetch for {len(ids_unique)} unique identifiers"
    )
    LOGGER.info(
        f"Using batch_size={batch_size}, max_workers={max_workers}, timeout={timeout}"
    )

    def _chunks(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i : i + n]

    results = {}
    symbol_to_ensembl = {}  # Map original symbol to Ensembl ID

    # Determine if inputs look like Ensembl IDs or gene symbols
    ensembl_id_pattern = re.compile(r"^ENS[A-Z]*[0-9]+$")
    likely_ensembl_ids = [i for i in ids_unique if ensembl_id_pattern.match(i)]
    likely_symbols = [i for i in ids_unique if not ensembl_id_pattern.match(i)]

    LOGGER.info(
        f"Found {len(likely_ensembl_ids)} Ensembl IDs and {len(likely_symbols)} gene symbols"
    )

    # 1a) Look up gene symbols using /lookup/symbol endpoint (for human)
    if likely_symbols:
        import re as re_module

        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        symbol_batches = list(_chunks(likely_symbols, batch_size))
        LOGGER.info(f"Processing {len(symbol_batches)} symbol batches")

        for batch in tqdm(symbol_batches, desc="Fetching by symbol", unit="batch"):
            try:
                # The symbol lookup endpoint requires species
                payload = {"symbols": batch, "species": "human"}
                r = session.post(
                    "https://rest.ensembl.org/lookup/symbol/homo_sapiens",
                    json=payload,
                    headers=headers,
                    timeout=timeout,
                )
                if r.status_code != 200:
                    LOGGER.warning(
                        "Ensembl symbol lookup failed (status %s) for %d symbols",
                        r.status_code,
                        len(batch),
                    )
                    # Log response body for debugging when available
                    try:
                        LOGGER.debug("Ensembl response body: %s", r.text[:1000])
                    except Exception:
                        pass
                    # Fall back to per-symbol GET lookups below
                try:
                    j = r.json()
                except Exception:
                    j = {}
                found_in_batch = 0
                for symbol, obj in j.items():
                    if not obj:
                        continue
                    ens_id = obj.get("id")
                    if ens_id:
                        symbol_to_ensembl[symbol] = ens_id
                        results[symbol] = {  # Use original symbol as key
                            "ensembl": ens_id,
                            "display_name": obj.get("display_name") or symbol,
                            "biotype": obj.get("biotype"),
                            "description": obj.get("description"),
                            "entrez": [],
                            "uniprot": [],
                        }
                        found_in_batch += 1
                LOGGER.debug(
                    f"Found {found_in_batch}/{len(batch)} symbols in this batch"
                )
                # If nothing was found in the batch response, try per-symbol GET lookups
                if found_in_batch == 0:
                    LOGGER.debug(
                        "Batch lookup returned 0 results, trying per-symbol GET for batch"
                    )
                    for symbol in batch:
                        try:
                            url_sym = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}"
                            r2 = session.get(url_sym, headers=headers, timeout=timeout)
                            if r2.status_code != 200:
                                LOGGER.debug(
                                    "Per-symbol lookup failed for %s (status %s)",
                                    symbol,
                                    r2.status_code,
                                )
                                continue
                            obj = r2.json()
                            if not obj:
                                continue
                            ens_id = obj.get("id")
                            if ens_id:
                                symbol_to_ensembl[symbol] = ens_id
                                results[symbol] = {
                                    "ensembl": ens_id,
                                    "display_name": obj.get("display_name") or symbol,
                                    "biotype": obj.get("biotype"),
                                    "description": obj.get("description"),
                                    "entrez": [],
                                    "uniprot": [],
                                }
                                found_in_batch += 1
                        except Exception as e:
                            LOGGER.debug(
                                "Per-symbol lookup exception for %s: %s", symbol, e
                            )
            except Exception as e:
                LOGGER.warning(
                    "Failed Ensembl symbol lookup for %d symbols: %s", len(batch), e
                )

    # 1b) Use POST /lookup/id in batches for Ensembl IDs
    if likely_ensembl_ids:
        lookup_url = "https://rest.ensembl.org/lookup/id"
        headers = {"Content-Type": "application/json", "Accept": "application/json"}

        id_batches = list(_chunks(likely_ensembl_ids, batch_size))
        LOGGER.info(f"Processing {len(id_batches)} ID batches")

        for batch in tqdm(id_batches, desc="Fetching by ID", unit="batch"):
            try:
                payload = {"ids": batch}
                r = session.post(
                    lookup_url, json=payload, headers=headers, timeout=timeout
                )
                if r.status_code != 200:
                    LOGGER.warning(
                        "Ensembl batch lookup failed (status %s) for %d ids",
                        r.status_code,
                        len(batch),
                    )
                    continue
                j = r.json()
                found_in_batch = 0
                for ens_id, obj in j.items():
                    if not obj:
                        continue
                    results[ens_id] = {
                        "ensembl": obj.get("id"),
                        "display_name": obj.get("display_name"),
                        "biotype": obj.get("biotype"),
                        "description": obj.get("description"),
                        "entrez": [],
                        "uniprot": [],
                    }
                    found_in_batch += 1
                LOGGER.debug(f"Found {found_in_batch}/{len(batch)} IDs in this batch")
            except Exception as e:
                LOGGER.warning(
                    "Failed Ensembl batch lookup for %d ids: %s", len(batch), e
                )

    # 2) Fetch cross-references (/xrefs/id) in parallel for the found IDs
    # We need to get the Ensembl IDs from all results
    ensembl_ids_for_xrefs = []
    for key, val in results.items():
        ens_id = val.get("ensembl")
        if ens_id:
            ensembl_ids_for_xrefs.append((key, ens_id))  # Keep track of original key

    if not ensembl_ids_for_xrefs:
        LOGGER.warning("No results found from lookups")
        return {}

    LOGGER.info(
        f"Fetching cross-references for {len(ensembl_ids_for_xrefs)} IDs with {max_workers} workers"
    )

    def _fetch_xrefs(original_key_and_ens_id):
        original_key, ens_id = original_key_and_ens_id
        try:
            url_xrefs = f"https://rest.ensembl.org/xrefs/id/{ens_id}"
            r2 = session.get(url_xrefs, timeout=timeout)
            if r2.status_code != 200:
                return original_key, None
            entrez = []
            uniprot = []
            for x in r2.json():
                dbname = x.get("dbname") or x.get("db_display_name") or ""
                pid = x.get("primary_id")
                if not pid:
                    continue
                dn = dbname.lower()
                if "entrez" in dn or "ncbigene" in dn:
                    entrez.append(pid)
                if dn.startswith("uniprot") or "uniprot" in dn:
                    uniprot.append(pid)
            return original_key, {
                "entrez": sorted(set(entrez)),
                "uniprot": sorted(set(uniprot)),
            }
        except Exception as e:
            LOGGER.debug("Failed to fetch xrefs for %s: %s", ens_id, e)
            return original_key, None

    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = {
            ex.submit(_fetch_xrefs, item): item for item in ensembl_ids_for_xrefs
        }
        xrefs_found = 0
        for fut in tqdm(
            as_completed(futures),
            total=len(futures),
            desc="Fetching xrefs",
            unit="gene",
        ):
            original_key, xr = fut.result()
            if xr and original_key in results:
                results[original_key]["entrez"] = xr.get("entrez", [])
                results[original_key]["uniprot"] = xr.get("uniprot", [])
                if xr.get("entrez") or xr.get("uniprot"):
                    xrefs_found += 1

    LOGGER.info(
        f"Completed: {len(results)} genes annotated, {xrefs_found} with cross-references"
    )
    return results
