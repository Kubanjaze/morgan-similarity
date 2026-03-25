"""
fetch_pubchem_smiles.py
-----------------------
Fetches canonical SMILES and properties from PubChem PUG REST API for a list
of KRAS G12C inhibitors. Writes a CSV you can copy into data/kras_g12c_inhibitors.csv.

Usage:
    python scripts/fetch_pubchem_smiles.py

Output:
    data/kras_g12c_inhibitors_pubchem.csv
"""

import csv
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

# ── Compounds to fetch ────────────────────────────────────────────────────────
# Each entry: (display_name, [query_strings_to_try_in_order])
COMPOUNDS = [
    ("sotorasib",   ["sotorasib", "AMG 510", "AMG510"]),
    ("adagrasib",   ["adagrasib", "MRTX849", "MRTX 849"]),
    ("divarasib",   ["divarasib", "GDC-6036", "GDC 6036"]),
    ("olomorasib",  ["olomorasib", "LY3537982", "LY 3537982"]),
    ("JDQ443",      ["JDQ443", "JDQ 443", "JDQ-443"]),
]

BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PROPERTIES = "CanonicalSMILES,IsomericSMILES,MolecularWeight,MolecularFormula"
SLEEP_BETWEEN = 0.4   # seconds — be polite to PubChem
TIMEOUT = 10          # seconds per request

OUTPUT_PATH = Path("data/kras_g12c_inhibitors_pubchem.csv")
OUTPUT_COLS = [
    "compound_name", "pubchem_query", "cid",
    "canonical_smiles", "isomeric_smiles",
    "molecular_weight", "molecular_formula", "status",
]


def fetch_json(url: str) -> dict | None:
    """GET url, return parsed JSON or None on error."""
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
            return json.loads(resp.read().decode())
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        print(f"    HTTP {e.code}: {url}")
        return None
    except Exception as e:
        print(f"    Error: {e}")
        return None


def get_cid(query: str) -> str | None:
    """Resolve a compound name/synonym to a PubChem CID."""
    url = f"{BASE_URL}/compound/name/{urllib.parse.quote(query)}/cids/JSON"
    data = fetch_json(url)
    if data and "IdentifierList" in data:
        cids = data["IdentifierList"].get("CID", [])
        if cids:
            return str(cids[0])
    return None


def get_properties(cid: str) -> dict | None:
    """Fetch compound properties by CID."""
    url = f"{BASE_URL}/compound/cid/{cid}/property/{PROPERTIES}/JSON"
    data = fetch_json(url)
    if data and "PropertyTable" in data:
        props = data["PropertyTable"].get("Properties", [])
        if props:
            return props[0]
    return None


def main() -> None:
    # urllib.parse not imported at top — fix:
    import urllib.parse  # noqa: F401 (used in get_cid)

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for display_name, queries in COMPOUNDS:
        print(f"{display_name}:")
        found = False
        for q in queries:
            time.sleep(SLEEP_BETWEEN)
            cid = get_cid(q)
            if cid:
                props = get_properties(cid)
                if props:
                    row = {
                        "compound_name":    display_name,
                        "pubchem_query":    q,
                        "cid":              cid,
                        "canonical_smiles": props.get("CanonicalSMILES", ""),
                        "isomeric_smiles":  props.get("IsomericSMILES", ""),
                        "molecular_weight": props.get("MolecularWeight", ""),
                        "molecular_formula":props.get("MolecularFormula", ""),
                        "status":           "ok",
                    }
                    rows.append(row)
                    print(f"  ok  query='{q}'  CID={cid}  MW={props.get('MolecularWeight','?')}")
                    found = True
                    break
        if not found:
            rows.append({
                "compound_name":    display_name,
                "pubchem_query":    "",
                "cid":              "",
                "canonical_smiles": "",
                "isomeric_smiles":  "",
                "molecular_weight": "",
                "molecular_formula":"",
                "status":           "not_found",
            })
            print(f"  not_found — tried: {queries}")

    with open(OUTPUT_PATH, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} rows to {OUTPUT_PATH}")
    print("Copy canonical_smiles values into data/kras_g12c_inhibitors.csv to populate missing SMILES.")


if __name__ == "__main__":
    import urllib.parse
    main()
