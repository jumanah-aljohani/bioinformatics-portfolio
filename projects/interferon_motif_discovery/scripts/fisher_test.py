"""
Fisher's Exact Test for ISRE enrichment comparison
between ISG promoters and control promoters.
"""

import pandas as pd
from pathlib import Path
from scipy.stats import fisher_exact


# ==============================
# Configuration
# ==============================

ISG_FILE = "results/ISG_promoters_300bp_GRCh38_isre_scan.tsv"
CTRL_FILE = "results/control_gene_promoters_300bp_GRCh38_isre_scan.tsv"

# Define strong threshold (0–1 mismatches = strong ISRE)
STRONG_THRESHOLD = 1


# ==============================
# Count strong vs weak
# ==============================

def count_strong_weak(filepath):
    df = pd.read_csv(filepath, sep="\t")

    strong = (df["mismatches"] <= STRONG_THRESHOLD).sum()
    weak = (df["mismatches"] > STRONG_THRESHOLD).sum()

    return strong, weak


def main():

    isg_strong, isg_weak = count_strong_weak(ISG_FILE)
    ctrl_strong, ctrl_weak = count_strong_weak(CTRL_FILE)

    table = [
        [isg_strong, isg_weak],
        [ctrl_strong, ctrl_weak]
    ]

    oddsratio, p_value = fisher_exact(table)

    print("Contingency Table:")
    print(table)
    print("\nOdds ratio:", oddsratio)
    print("P-value:", p_value)

    # Save results
    out_dir = Path("results")
    out_dir.mkdir(parents=True, exist_ok=True)

    out_path = out_dir / "fisher_test_results.txt"

    out_path.write_text(
        f"Contingency Table:\n{table}\n\n"
        f"Odds ratio: {oddsratio}\n"
        f"P-value: {p_value}\n"
    )

    print(f"\nResults saved to: {out_path.resolve()}")


if __name__ == "__main__":
    main()
