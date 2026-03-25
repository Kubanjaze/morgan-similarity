# Phase 05 — Morgan Fingerprint Similarity Explorer

**Version:** 1.1 (final as-built)
**Author:** Kerwyn Medrano
**Date:** 2026-03-25
**Track:** Track 1 — Cheminformatics Core
**Tier:** Micro (2–3 hrs)
**API Cost:** $0.00 — pure RDKit + pandas + seaborn
**Status:** ✅ Complete
**Repo:** `morgan-similarity/`

## As-Built Results

| Pair | Tanimoto |
|---|---|
| sotorasib ↔ divarasib | 0.336 (highest) |
| divarasib ↔ adagrasib | 0.295 |
| divarasib ↔ olomorasib | 0.246 |
| sotorasib ↔ olomorasib | 0.233 |
| adagrasib ↔ olomorasib | 0.164 |
| sotorasib ↔ adagrasib | 0.153 (lowest) |

Mean similarity: 0.238. All four inhibitors are structurally distant despite targeting the same covalent pocket — diverse scaffolds converge on Cys12 alkylation.

**Deviation from plan:** `scipy` added to `requirements.txt` (required by `seaborn.clustermap`). JDQ443 SMILES left blank — PubChem fetch script provided.

---

---

## 1. Project Overview

### Goal

Given a CSV of compounds with SMILES, compute pairwise Morgan fingerprint similarity using the Tanimoto coefficient and visualise the result as a clustered heatmap.

A medicinal chemist should be able to run:

```bash
python main.py --input data/kras_g12c_inhibitors.csv --radius 2 --nbits 2048
```

...and get back:
- `output/tanimoto_matrix.csv` — N×N similarity matrix
- `output/tanimoto_heatmap.png` — clustered heatmap with compound names

### What This Phase Teaches

| Concept | Detail |
|---|---|
| Morgan fingerprints | Circular fingerprint parameterisation: radius, nBits, chirality |
| Tanimoto coefficient | Bit-vector Jaccard similarity via `DataStructs.TanimotoSimilarity` |
| Pairwise matrix | Efficient construction with `np.zeros` + symmetric fill |
| Clustered heatmap | `seaborn.clustermap` with linkage-based reordering |
| SMILES validation | Explicit failure counting, skip-and-warn pattern |

### Domain Context

Morgan fingerprints (ECFP-style) are the standard 2D similarity metric in drug discovery. Tanimoto > 0.4 is considered structurally related; > 0.85 is near-identical. For KRAS G12C inhibitors:
- Sotorasib and adagrasib share a piperazine-acrylamide warhead but differ significantly in the hinge-binding region — expect Tanimoto ~0.3–0.5
- Divarasib and sotorasib both hit the switch-II pocket covalently — expect moderate similarity
- Olomorasib has a distinct thienopyrimidine core — expect lower similarity to the others

The heatmap should surface these structural relationships visually and confirm or challenge prior SAR intuition.

---

## 2. Architecture

```
morgan-similarity/
├── main.py              ← CLI: load → fingerprint → matrix → visualise
├── data/
│   └── kras_g12c_inhibitors.csv   ← seed dataset (5 approved/Phase III)
├── output/              ← tanimoto_matrix.csv + tanimoto_heatmap.png
└── requirements.txt
```

Single-file design — no `src/` needed. All logic in `main.py` as clearly named functions.

---

## 3. Module Specification

### 3.1 `load_compounds(path) → pd.DataFrame`
- Read CSV, require `compound_name` + `smiles` columns
- Validate each SMILES with `Chem.MolFromSmiles()`
- Drop invalid rows, print warning with count
- Return DataFrame with added `mol` column

### 3.2 `compute_fingerprints(df, radius, nbits, chirality) → list[ExplicitBitVect]`
- Generate Morgan fingerprint for each `mol`
- Parameters explicit: `radius=2`, `nBits=2048`, `useChirality=False` (default)
- Return list aligned to DataFrame index

### 3.3 `tanimoto_matrix(fps) → np.ndarray`
- N×N symmetric matrix
- `DataStructs.TanimotoSimilarity(fps[i], fps[j])` for all pairs
- Diagonal = 1.0

### 3.4 `plot_heatmap(matrix, labels, output_path)`
- `seaborn.clustermap` with `method='average'`, `metric='euclidean'` on `1 - matrix`
- Annotate cells with similarity value (2 decimal places) if N ≤ 20
- Save to `output_path` at 150 dpi

### 3.5 `main()` — CLI
- `--input` (required): path to CSV
- `--radius` (default: 2): Morgan radius
- `--nbits` (default: 2048): fingerprint bit length
- `--chirality` (flag): enable chirality in fingerprint
- `--output-dir` (default: `output`): where to write artifacts
- Print summary: N compounds loaded, N valid, matrix shape, output paths

---

## 4. CLI Reference

```bash
# Default run
python main.py --input data/kras_g12c_inhibitors.csv

# Custom fingerprint params
python main.py --input data/kras_g12c_inhibitors.csv --radius 3 --nbits 1024

# With chirality
python main.py --input data/kras_g12c_inhibitors.csv --chirality
```

---

## 5. Input / Output Contract

### Input CSV (required columns)
```
compound_name  — display label in heatmap
smiles         — SMILES string (invalid rows skipped with warning)
```
All other columns are carried through but ignored.

### Outputs
```
output/tanimoto_matrix.csv    — N×N matrix, compound names as row/col headers
output/tanimoto_heatmap.png   — clustered heatmap
```

---

## 6. Key Patterns Being Learned

| Pattern | Location | What it teaches |
|---|---|---|
| `Chem.MolFromSmiles()` + None check | `load_compounds` | Standard RDKit SMILES validation idiom |
| `AllChem.GetMorganFingerprintAsBitVect()` | `compute_fingerprints` | Explicit radius/nBits/chirality params |
| `DataStructs.TanimotoSimilarity()` | `tanimoto_matrix` | Tanimoto on bit vectors |
| Symmetric matrix fill | `tanimoto_matrix` | Efficient pairwise computation |
| `seaborn.clustermap` on distance matrix | `plot_heatmap` | 1 - similarity = distance for clustering |

---

## 7. Verification Checklist

```bash
# Run on seed data
python main.py --input data/kras_g12c_inhibitors.csv

# Expected:
# - 5 compounds loaded, 5 valid
# - output/tanimoto_matrix.csv exists (5×5)
# - output/tanimoto_heatmap.png exists
# - Diagonal values = 1.0
# - Sotorasib-adagrasib Tanimoto in range 0.2–0.6
# - Olomorasib lowest similarity to the rest (distinct scaffold)

# Test with invalid SMILES row
echo "bad_compound,INVALIDSMILES,test,0,C" >> data/test_invalid.csv
python main.py --input data/test_invalid.csv
# Expected: warning printed, bad row skipped, run completes
```

---

## 8. Dependencies

```
rdkit
pandas
numpy
seaborn
matplotlib
```

All already present in the cheminformatics environment.

---

## 9. Risks / Assumptions / Next Step

**Risks:**
- RDKit import path differs across environments — use `from rdkit import Chem`
- `seaborn.clustermap` requires N ≥ 2 valid compounds — add guard
- Very large CSVs (N > 500): matrix becomes slow — note in CLI help, not in scope here

**Assumptions:**
- Python 3.10+, RDKit installed via conda or pip
- SMILES are canonical or at minimum parseable by RDKit

**Next step:** Phase 06 — Murcko scaffold extraction + frequency ranking. Uses the same compound CSV; asks what scaffolds appear most often in the KRAS G12C series and whether scaffold frequency correlates with potency.
