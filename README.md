# morgan-similarity — Phase 05

Compute pairwise Morgan fingerprint (ECFP) similarity between compounds using the Tanimoto coefficient and visualise the result as a clustered heatmap.

## What it does

Morgan fingerprints encode the circular neighbourhood of each atom up to a given radius. The Tanimoto coefficient measures the overlap between two bit vectors — the standard 2D similarity metric in drug discovery:

- Tanimoto ≥ 0.85 → near-identical structures
- Tanimoto 0.4–0.85 → structurally related (same scaffold family)
- Tanimoto < 0.4 → structurally distant

The clustered heatmap reorders compounds by structural similarity, surfacing scaffold families without any explicit scaffold perception.

## Inputs / Outputs

**Input CSV** (required columns):
```
compound_name   display label
smiles          SMILES string — invalid rows are skipped with a warning
```

**Outputs:**
```
output/tanimoto_matrix.csv     N×N similarity matrix (row/col = compound_name)
output/tanimoto_heatmap.png    clustered heatmap at 150 dpi
```

## Setup

### RDKit (recommended via conda)
```bash
conda create -n morgan-sim python=3.11 rdkit -c conda-forge
conda activate morgan-sim
pip install -r requirements.txt
```

### pip only
```bash
pip install rdkit
pip install -r requirements.txt
```

## Run

```bash
# Default (radius=2, nBits=2048)
python main.py --input data/kras_g12c_inhibitors.csv

# Custom fingerprint params
python main.py --input data/kras_g12c_inhibitors.csv --radius 3 --nbits 1024

# Include chirality
python main.py --input data/kras_g12c_inhibitors.csv --chirality

# All options
python main.py --help
```

## Populate missing SMILES

The seed dataset ships with 4 confirmed SMILES. JDQ443 is left blank.
Run the fetch script to pull all SMILES from PubChem and copy them in:

```bash
python scripts/fetch_pubchem_smiles.py
# writes data/kras_g12c_inhibitors_pubchem.csv
# copy canonical_smiles values into data/kras_g12c_inhibitors.csv
```

## Notes

- Invalid or empty SMILES are skipped — the run does not crash
- Heatmap requires ≥ 2 valid compounds; skipped silently if fewer
- Cell annotations are shown only when N ≤ 20
- Complexity is O(N²) — large libraries (N > 500) will be slow; use a batched approach
- Clustering uses `method='average'` linkage on `distance = 1 - similarity`
