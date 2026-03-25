import sys
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


# ── Data loading ──────────────────────────────────────────────────────────────

def load_compounds(path: Path) -> pd.DataFrame:
    """
    Load CSV, validate SMILES, return DataFrame with compound_name, smiles, mol columns.
    Invalid SMILES rows are dropped with a warning.
    """
    df = pd.read_csv(path)
    for col in ("compound_name", "smiles"):
        if col not in df.columns:
            sys.exit(f"ERROR: CSV must contain a '{col}' column.")

    total = len(df)
    mols, valid_mask = [], []
    for _, row in df.iterrows():
        smi = str(row["smiles"]).strip() if pd.notna(row["smiles"]) else ""
        if not smi:
            mols.append(None)
            valid_mask.append(False)
            continue
        mol = Chem.MolFromSmiles(smi)
        mols.append(mol)
        valid_mask.append(mol is not None)

    df["mol"] = mols
    df = df[valid_mask][["compound_name", "smiles", "mol"]].reset_index(drop=True)
    n_valid = len(df)
    n_invalid = total - n_valid
    print(f"Loaded {total} rows, {n_valid} valid, {n_invalid} invalid skipped.")
    if n_invalid:
        print(f"  (Tip: run scripts/fetch_pubchem_smiles.py to populate missing SMILES)", file=sys.stderr)
    return df


# ── Fingerprints ──────────────────────────────────────────────────────────────

def compute_fingerprints(
    df: pd.DataFrame,
    radius: int = 2,
    nbits: int = 2048,
    chirality: bool = False,
) -> list:
    """Generate Morgan fingerprints for each compound row."""
    fps = []
    for mol in df["mol"]:
        fp = AllChem.GetMorganFingerprintAsBitVect(
            mol, radius=radius, nBits=nbits, useChirality=chirality
        )
        fps.append(fp)
    return fps


# ── Tanimoto matrix ───────────────────────────────────────────────────────────

def tanimoto_matrix(fps: list) -> np.ndarray:
    """Compute N×N pairwise Tanimoto similarity matrix (upper triangle + mirror)."""
    n = len(fps)
    matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        matrix[i, i] = 1.0
        for j in range(i + 1, n):
            sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            matrix[i, j] = sim
            matrix[j, i] = sim
    return matrix


# ── Heatmap ───────────────────────────────────────────────────────────────────

def plot_heatmap(
    matrix: np.ndarray,
    labels: list[str],
    output_path: Path,
) -> None:
    """Clustered heatmap of Tanimoto similarity. Clusters by distance = 1 - similarity."""
    n = len(labels)
    if n < 2:
        print("  Skipping heatmap: need at least 2 valid compounds.")
        return

    df_sim = pd.DataFrame(matrix, index=labels, columns=labels)
    annot = n <= 20
    figsize = max(6, n * 1.2)

    g = sns.clustermap(
        df_sim,
        method="average",
        metric="euclidean",
        cmap="YlOrRd",
        annot=annot,
        fmt=".2f" if annot else "",
        linewidths=0.5 if annot else 0,
        figsize=(figsize, figsize),
        vmin=0,
        vmax=1,
    )
    g.figure.suptitle("Tanimoto Similarity — Morgan Fingerprints (ECFP)", y=1.02, fontsize=13)
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")
    plt.tight_layout()
    g.figure.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Morgan fingerprint Tanimoto similarity matrix + clustered heatmap.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", required=True, help="CSV with compound_name + smiles columns")
    parser.add_argument("--radius", type=int, default=2, help="Morgan radius (e.g. 2 = ECFP4)")
    parser.add_argument("--nbits", type=int, default=2048, help="Fingerprint bit vector length")
    parser.add_argument("--chirality", action="store_true", help="Include chirality in fingerprint")
    parser.add_argument("--output-dir", default="output", help="Directory for output artifacts")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load + validate
    df = load_compounds(Path(args.input))
    n = len(df)
    if n < 2:
        sys.exit("ERROR: Need at least 2 valid compounds to compute a similarity matrix.")

    # Fingerprints
    print(f"Computing Morgan fingerprints (radius={args.radius}, nBits={args.nbits}, chirality={args.chirality})...")
    fps = compute_fingerprints(df, radius=args.radius, nbits=args.nbits, chirality=args.chirality)

    # Tanimoto matrix
    print(f"Computing {n}×{n} Tanimoto matrix...")
    matrix = tanimoto_matrix(fps)

    # Save matrix CSV
    labels = df["compound_name"].tolist()
    matrix_path = output_dir / "tanimoto_matrix.csv"
    pd.DataFrame(matrix, index=labels, columns=labels).to_csv(matrix_path, float_format="%.4f")
    print(f"  Saved: {matrix_path}")

    # Heatmap
    print("Plotting clustered heatmap...")
    plot_heatmap(matrix, labels, output_dir / "tanimoto_heatmap.png")

    # Summary stats
    upper = matrix[np.triu_indices(n, k=1)]
    print(f"\nSimilarity summary ({n}×{n}, {len(upper)} pairs):")
    print(f"  Mean : {upper.mean():.3f}")
    print(f"  Min  : {upper.min():.3f}")
    print(f"  Max  : {upper.max():.3f}")

    # Full matrix to stdout
    print(f"\nTanimoto matrix:")
    print(pd.DataFrame(matrix, index=labels, columns=labels).to_string(float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()
