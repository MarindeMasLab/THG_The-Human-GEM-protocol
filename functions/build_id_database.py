#!/usr/bin/env python3
"""
Step 2: Build a comprehensive metabolite ID database.
Takes metabolite annotation data from analyze_annotations and creates a structured database.
"""

import json
from pathlib import Path

# Metabolite ID database - manually curated from multiple sources
# Structure: metabolite_name -> {database -> ID(s)}
METABOLITE_ID_DATABASE = {
    # Energy metabolites
    "ATP": {
        "kegg.compound": "C00002",
        "chebi": "CHEBI:15422",
        "bigg.metabolite": "atp",
        "pubchem.compound": "5957",
        "vmhmetabolite": "atp",
        "metanetx.chemical": "MNXM3",
        "inchikey": "ZKHQWZAMYRWXGA-KQYNXXCUSA-N"
    },
    "ADP": {
        "kegg.compound": "C00008",
        "chebi": "CHEBI:16761",
        "bigg.metabolite": "adp",
        "pubchem.compound": "6022",
        "vmhmetabolite": "adp",
        "metanetx.chemical": "MNXM7",
        "inchikey": "XTWYTFMLZFPYCI-KQYNXXCUSA-N"
    },
    "GTP": {
        "kegg.compound": "C00044",
        "chebi": "CHEBI:15996",
        "bigg.metabolite": "gtp",
        "pubchem.compound": "6830",
        "vmhmetabolite": "gtp",
        "metanetx.chemical": "MNXM68",
        "inchikey": "XKMLYUALXHKNFT-UUOKFMHZSA-N"
    },
    "GDP": {
        "kegg.compound": "C00035",
        "chebi": "CHEBI:17552",
        "bigg.metabolite": "gdp",
        "pubchem.compound": "8977",
        "vmhmetabolite": "gdp",
        "metanetx.chemical": "MNXM59",
        "inchikey": "QGWNDRXFNXRZMB-UUOKFMHZSA-N"
    },
    "H2O": {
        "kegg.compound": "C00001",
        "chebi": "CHEBI:15377",
        "bigg.metabolite": "h2o",
        "pubchem.compound": "962",
        "vmhmetabolite": "h2o",
        "metanetx.chemical": "MNXM2",
        "inchikey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
    },
    "phosphate": {
        "kegg.compound": "C00009",
        "chebi": "CHEBI:43474",
        "bigg.metabolite": "pi",
        "pubchem.compound": "1061",
        "vmhmetabolite": "pi",
        "metanetx.chemical": "MNXM9",
        "inchikey": "NBIIXXVUZAFLBC-UHFFFAOYSA-K"
    },
    "H+": {
        "kegg.compound": "C00080",
        "chebi": "CHEBI:15378",
        "bigg.metabolite": "h",
        "pubchem.compound": "1038",
        "vmhmetabolite": "h",
        "metanetx.chemical": "MNXM89",
        "inchikey": "GPRLSGONYQIRFK-UHFFFAOYSA-N"
    },
    
    # Amino acids
    "L-alanine": {
        "kegg.compound": "C00041",
        "chebi": "CHEBI:16977",
        "bigg.metabolite": "ala__L",
        "pubchem.compound": "5950",
        "vmhmetabolite": "ala_L",
        "metanetx.chemical": "MNXM66",
        "inchikey": "QNAYBMKLOCPYGJ-REOHCLBHSA-N"
    },
    "L-arginine": {
        "kegg.compound": "C00062",
        "chebi": "CHEBI:16467",
        "bigg.metabolite": "arg__L",
        "pubchem.compound": "6322",
        "vmhmetabolite": "arg_L",
        "metanetx.chemical": "MNXM83",
        "inchikey": "ODKSFYDXXFIFQN-BYPYZUCNSA-N"
    },
    "L-asparagine": {
        "kegg.compound": "C00152",
        "chebi": "CHEBI:17196",
        "bigg.metabolite": "asn__L",
        "pubchem.compound": "6267",
        "vmhmetabolite": "asn_L",
        "metanetx.chemical": "MNXM168",
        "inchikey": "DCXYFEDJOCDNAF-REOHCLBHSA-N"
    },
    "L-aspartate": {
        "kegg.compound": "C00049",
        "chebi": "CHEBI:17053",
        "bigg.metabolite": "asp__L",
        "pubchem.compound": "5960",
        "vmhmetabolite": "asp_L",
        "metanetx.chemical": "MNXM75",
        "inchikey": "CKLJMWTZIZZHCS-REOHCLBHSA-N"
    },
    "L-cysteine": {
        "kegg.compound": "C00097",
        "chebi": "CHEBI:17561",
        "bigg.metabolite": "cys__L",
        "pubchem.compound": "5862",
        "vmhmetabolite": "cys_L",
        "metanetx.chemical": "MNXM113",
        "inchikey": "XUJNEKJLAYXESH-REOHCLBHSA-N"
    },
    "L-glutamine": {
        "kegg.compound": "C00064",
        "chebi": "CHEBI:18050",
        "bigg.metabolite": "gln__L",
        "pubchem.compound": "5961",
        "vmhmetabolite": "gln_L",
        "metanetx.chemical": "MNXM85",
        "inchikey": "ZDXPYRJPNDTMRX-VKHMYHEASA-N"
    },
    "L-glutamate": {
        "kegg.compound": "C00025",
        "chebi": "CHEBI:16015",
        "bigg.metabolite": "glu__L",
        "pubchem.compound": "33032",
        "vmhmetabolite": "glu_L",
        "metanetx.chemical": "MNXM46",
        "inchikey": "WHUUTDBJXJRKMK-VKHMYHEASA-N"
    },
    "glycine": {
        "kegg.compound": "C00037",
        "chebi": "CHEBI:15428",
        "bigg.metabolite": "gly",
        "pubchem.compound": "750",
        "vmhmetabolite": "gly",
        "metanetx.chemical": "MNXM57",
        "inchikey": "DHMQDGOQFOQNFH-UHFFFAOYSA-N"
    },
    "L-histidine": {
        "kegg.compound": "C00135",
        "chebi": "CHEBI:15971",
        "bigg.metabolite": "his__L",
        "pubchem.compound": "6274",
        "vmhmetabolite": "his_L",
        "metanetx.chemical": "MNXM151",
        "inchikey": "HNDVDQJCIGZPNO-YFKPBYRVSA-N"
    },
    "L-isoleucine": {
        "kegg.compound": "C00407",
        "chebi": "CHEBI:17191",
        "bigg.metabolite": "ile__L",
        "pubchem.compound": "6306",
        "vmhmetabolite": "ile_L",
        "metanetx.chemical": "MNXM421",
        "inchikey": "AGPKZVBTJJNPAG-WHFBIAKZSA-N"
    },
    "L-leucine": {
        "kegg.compound": "C00123",
        "chebi": "CHEBI:15603",
        "bigg.metabolite": "leu__L",
        "pubchem.compound": "6106",
        "vmhmetabolite": "leu_L",
        "metanetx.chemical": "MNXM139",
        "inchikey": "ROHFNLRQFUQHCH-YFKPBYRVSA-N"
    },
    "L-lysine": {
        "kegg.compound": "C00047",
        "chebi": "CHEBI:18019",
        "bigg.metabolite": "lys__L",
        "pubchem.compound": "5962",
        "vmhmetabolite": "lys_L",
        "metanetx.chemical": "MNXM71",
        "inchikey": "KDXKERNSBIXSRK-YFKPBYRVSA-N"
    },
    "L-methionine": {
        "kegg.compound": "C00073",
        "chebi": "CHEBI:16643",
        "bigg.metabolite": "met__L",
        "pubchem.compound": "6137",
        "vmhmetabolite": "met_L",
        "metanetx.chemical": "MNXM96",
        "inchikey": "FFEARJCKVFRZRR-BYPYZUCNSA-N"
    },
    "L-phenylalanine": {
        "kegg.compound": "C00079",
        "chebi": "CHEBI:17295",
        "bigg.metabolite": "phe__L",
        "pubchem.compound": "6140",
        "vmhmetabolite": "phe_L",
        "metanetx.chemical": "MNXM102",
        "inchikey": "COLNVLDHVKWLRT-QMMMGPOBSA-N"
    },
    "L-proline": {
        "kegg.compound": "C00148",
        "chebi": "CHEBI:17203",
        "bigg.metabolite": "pro__L",
        "pubchem.compound": "145742",
        "vmhmetabolite": "pro_L",
        "metanetx.chemical": "MNXM164",
        "inchikey": "ONIBWKKTOPOVIA-BYPYZUCNSA-N"
    },
    "L-serine": {
        "kegg.compound": "C00065",
        "chebi": "CHEBI:17115",
        "bigg.metabolite": "ser__L",
        "pubchem.compound": "5951",
        "vmhmetabolite": "ser_L",
        "metanetx.chemical": "MNXM86",
        "inchikey": "MTCFGRXMJLQNBG-REOHCLBHSA-N"
    },
    "L-threonine": {
        "kegg.compound": "C00188",
        "chebi": "CHEBI:16857",
        "bigg.metabolite": "thr__L",
        "pubchem.compound": "6288",
        "vmhmetabolite": "thr_L",
        "metanetx.chemical": "MNXM203",
        "inchikey": "AYFVYJQAPQTCCC-GBXIJSLDSA-N"
    },
    "L-tryptophan": {
        "kegg.compound": "C00078",
        "chebi": "CHEBI:16828",
        "bigg.metabolite": "trp__L",
        "pubchem.compound": "6305",
        "vmhmetabolite": "trp_L",
        "metanetx.chemical": "MNXM101",
        "inchikey": "QIVBCDIJIAJPQS-VIFPVBQESA-N"
    },
    "L-tyrosine": {
        "kegg.compound": "C00082",
        "chebi": "CHEBI:17895",
        "bigg.metabolite": "tyr__L",
        "pubchem.compound": "6057",
        "vmhmetabolite": "tyr_L",
        "metanetx.chemical": "MNXM105",
        "inchikey": "OUYCCCASQSFEME-QMMMGPOBSA-N"
    },
    "L-valine": {
        "kegg.compound": "C00183",
        "chebi": "CHEBI:16414",
        "bigg.metabolite": "val__L",
        "pubchem.compound": "6287",
        "vmhmetabolite": "val_L",
        "metanetx.chemical": "MNXM198",
        "inchikey": "KZSNJWFQEVHDMF-BYPYZUCNSA-N"
    },
    
    # Regulatory metabolites
    "UDP-N-acetyl-alpha-D-glucosamine": {
        "kegg.compound": "C00043",
        "chebi": "CHEBI:16264",
        "bigg.metabolite": "uacgam",
        "pubchem.compound": "445675",
        "vmhmetabolite": "uacgam",
        "metanetx.chemical": "MNXM71",
        "inchikey": "LFTYTUAZOPRMMI-CFRASDGPSA-N"
    },
    "UDP": {
        "kegg.compound": "C00015",
        "chebi": "CHEBI:17659",
        "bigg.metabolite": "udp",
        "pubchem.compound": "6031",
        "vmhmetabolite": "udp",
        "metanetx.chemical": "MNXM26",
        "inchikey": "XCCTYIAWTASOJW-XVFCMESISA-N"
    },
    "1-phosphatidyl-1D-myo-inositol 4-phosphate": {
        "kegg.compound": "C04637",
        "chebi": "CHEBI:17526",
        "bigg.metabolite": "pi4p",
        "vmhmetabolite": "pi4p_hs",
        "metanetx.chemical": "MNXM1105",
        "inchikey": "MLRXMHJAPDAXBE-GVYPWVHMSA-N"
    }
}

def build_database_from_annotations(metabolite_annotations):
    """
    Build metabolite ID database from extracted annotations.
    
    Args:
        metabolite_annotations: Dict from analyze_annotations.extract_metabolite_annotations()
        
    Returns:
        dict: Standardized metabolite ID database
    """
    database = {}
    
    for met_name, annotations in metabolite_annotations.items():
        # Convert annotation keys to standard format
        standardized = {}
        
        for key, value in annotations.items():
            # Handle list values
            if isinstance(value, list):
                standardized[key] = value
            else:
                standardized[key] = value
        
        database[met_name] = standardized
    
    return database


def save_database(database, output_path):
    """Save the metabolite ID database to a JSON file."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump(database, f, indent=2)
    
    print(f"\n✓ Metabolite ID database saved to: {output_path}")
    print(f"  Total metabolites in database: {len(database)}")
    
    return output_path


def main():
    """Main function - connects to analyze_annotations and creates database."""
    print("=" * 70)
    print("STEP 2: BUILDING METABOLITE ID DATABASE")
    print("=" * 70)
    
    try:
        # Import and run annotation analysis
        try:
            from functions.analyze_annotations import extract_metabolite_annotations
            from functions.config import load_config, get_model_paths
        except ModuleNotFoundError:
            # Direct execution - adjust sys.path
            import sys
            from pathlib import Path
            sys.path.insert(0, str(Path(__file__).parent.parent))
            from functions.analyze_annotations import extract_metabolite_annotations
            from functions.config import load_config, get_model_paths
        
        # Load configuration
        config = load_config()
        model_path, _ = get_model_paths(config)
        targets = config.get('metabolites', {}).get('targets', [])
        output_path = config.get('metabolites', {}).get('database_file', 'data/metabolite_id_database.json')
        
        if not targets:
            print("✗ Error: No metabolite targets found in config.json")
            print("  Please add 'metabolites.targets' to config.json")
            return None
        
        print(f"\nExtracting annotations for {len(targets)} metabolites from model...")
        
        # Extract annotations
        metabolite_annotations, not_found = extract_metabolite_annotations(
            model_path, targets, verbose=False
        )
        
        if not_found:
            print(f"\n⚠ Warning: {len(not_found)} metabolites not found:")
            for met in not_found[:5]:  # Show first 5
                print(f"    - {met}")
            if len(not_found) > 5:
                print(f"    ... and {len(not_found) - 5} more")
        
        # Build database
        print(f"\nBuilding database from {len(metabolite_annotations)} annotations...")
        database = build_database_from_annotations(metabolite_annotations)
        
        # Save to file
        save_database(database, output_path)
        
        print("\n" + "=" * 70)
        print("Database creation complete!")
        print("=" * 70)
        
        return database
        
    except Exception as e:
        print(f"\n✗ Error: {e}")
        print("\nFalling back to hardcoded database...")
        print("(Consider fixing config.json and re-running)")
        
        # Fallback: save hardcoded database
        output_path = 'data/metabolite_id_database.json'
        save_database(METABOLITE_ID_DATABASE, output_path)
        return METABOLITE_ID_DATABASE

if __name__ == '__main__':
    main()
