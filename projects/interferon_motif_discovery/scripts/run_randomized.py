"""
Randomized Motif Search implementation.
Applied to ISG promoter sequences.
"""

import random
from pathlib import Path
from typing import Dict, List


# ==============================
# Configuration
# ==============================

FASTA_PATH = "data/ISG_promoters_300bp_GRCh38.fa"
K = 12
RESTARTS = 1000
SEED = 42


# ==============================
# FASTA reader
# ==============================

def read_fasta(fp: str) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    name = None
    parts: List[str] = []

    for line in Path(fp).read_text().splitlines():
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(parts).upper()
            name = line[1:].split()[0]
            parts = []
        else:
            parts.append(line)

    if name is not None:
        seqs[name] = "".join(parts).upper()

    return seqs


# ==============================
# Score function
# ==============================

def score_motifs(motifs: List[str]) -> int:
    t = len(motifs)
    k = len(motifs[0])
    total = 0

    for i in range(k):
        column = [m[i] for m in motifs]
        max_count = max(column.count(b) for b in "ACGT")
        total += t - max_count

    return total


# ==============================
# Build profile with pseudocounts
# ==============================

def build_profile(motifs: List[str]) -> Dict[str, List[float]]:
    k = len(motifs[0])
    profile = {b: [1] * k for b in "ACGT"}  # pseudocounts

    for m in motifs:
        for i, ch in enumerate(m):
            profile[ch][i] += 1

    t = len(motifs)
    for i in range(k):
        col_sum = sum(profile[b][i] for b in "ACGT")
        for b in "ACGT":
            profile[b][i] /= col_sum

    return profile


def most_probable_kmer(seq: str, k: int, profile: Dict[str, List[float]]) -> str:
    best_prob = -1
    best_kmer = seq[0:k]

    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        prob = 1.0
        for j, ch in enumerate(kmer):
            prob *= profile[ch][j]

        if prob > best_prob:
            best_prob = prob
            best_kmer = kmer

    return best_kmer


# ==============================
# Randomized Motif Search
# ==============================

def randomized_motif_search(dna: List[str], k: int) -> List[str]:
    motifs = []

    for seq in dna:
        start = random.randrange(len(seq) - k + 1)
        motifs.append(seq[start:start + k])

    best = motifs[:]
    best_score = score_motifs(best)

    while True:
        profile = build_profile(motifs)
        motifs = [most_probable_kmer(seq, k, profile) for seq in dna]
        current_score = score_motifs(motifs)

        if current_score < best_score:
            best = motifs[:]
            best_score = current_score
        else:
            return best


# ==============================
# Consensus calculation
# ==============================

def consensus_from_motifs(motifs: List[str]) -> str:
    k = len(motifs[0])
    bases = "ACGT"
    consensus = ""

    for i in range(k):
        column = [m[i] for m in motifs]
        consensus += max(bases, key=lambda b: column.count(b))

    return consensus


# ==============================
# Main execution
# ==============================

def main():
    if SEED is not None:
        random.seed(SEED)

    seqs = read_fasta(FASTA_PATH)
    genes = list(seqs.keys())
    dna = [seqs[g] for g in genes]

    global_best = None
    global_best_score = None

    for _ in range(RESTARTS):
        motifs = randomized_motif_search(dna, K)
        sc = score_motifs(motifs)

        if global_best is None or sc < global_best_score:
            global_best = motifs
            global_best_score = sc

    consensus = consensus_from_motifs(global_best)

    print("\nBest score:", global_best_score)
    print("Consensus:", consensus)
    print("\nGene\tMotif")

    for g, m in zip(genes, global_best):
        print(f"{g}\t{m}")

    # Save results
    out_dir = Path("results")
    out_dir.mkdir(parents=True, exist_ok=True)

    out_path = out_dir / f"randomized_k{K}_restarts{RESTARTS}.tsv"

    out_path.write_text(
        "gene\tmotif\n"
        + "\n".join(f"{g}\t{m}" for g, m in zip(genes, global_best))
        + f"\n\n# Best score: {global_best_score}"
        + f"\n# Consensus: {consensus}\n"
    )

    print(f"\nResults saved to: {out_path.resolve()}")


if __name__ == "__main__":
    main()
