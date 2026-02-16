"""
Gibbs Sampling implementation for de novo motif discovery.
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
N = 2000
RESTARTS = 50
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
    counts = {b: [1] * k for b in "ACGT"}  # pseudocounts

    for m in motifs:
        for i, ch in enumerate(m):
            counts[ch][i] += 1

    t = len(motifs)
    denom = t + 4

    profile = {b: [counts[b][i] / denom for i in range(k)] for b in "ACGT"}
    return profile


def kmer_probability(kmer: str, profile: Dict[str, List[float]]) -> float:
    p = 1.0
    for i, ch in enumerate(kmer):
        p *= profile[ch][i]
    return p


def weighted_choice(weights: List[float]) -> int:
    s = sum(weights)
    if s == 0:
        return random.randrange(len(weights))

    r = random.random() * s
    acc = 0.0
    for i, w in enumerate(weights):
        acc += w
        if acc >= r:
            return i

    return len(weights) - 1


def profile_random_kmer(seq: str, k: int, profile: Dict[str, List[float]]) -> str:
    kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]
    weights = [kmer_probability(km, profile) for km in kmers]
    idx = weighted_choice(weights)
    return kmers[idx]


# ==============================
# Gibbs Sampler (single run)
# ==============================

def gibbs_sampler(dna: List[str], k: int, n: int) -> List[str]:
    t = len(dna)

    motifs = []
    for seq in dna:
        start = random.randrange(0, len(seq) - k + 1)
        motifs.append(seq[start:start + k])

    best = motifs[:]
    best_score = score_motifs(best)

    for _ in range(n):
        i = random.randrange(t)
        motifs_excluding_i = motifs[:i] + motifs[i + 1:]

        profile = build_profile(motifs_excluding_i)
        motifs[i] = profile_random_kmer(dna[i], k, profile)

        sc = score_motifs(motifs)

        if sc < best_score:
            best = motifs[:]
            best_score = sc

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
        motifs = gibbs_sampler(dna, K, N)
        sc = score_motifs(motifs)

        if global_best is None or sc < global_best_score:
            global_best = motifs
            global_best_score = sc

    consensus = consensus_from_motifs(global_best)

    print(f"\nBest score: {global_best_score}")
    print("Consensus:", consensus)
    print("\nGene\tMotif")

    for g, m in zip(genes, global_best):
        print(f"{g}\t{m}")

    # Save results
    out_dir = Path("results")
    out_dir.mkdir(parents=True, exist_ok=True)

    out_path = out_dir / f"gibbs_k{K}_N{N}_restarts{RESTARTS}.tsv"

    out_path.write_text(
        "gene\tmotif\n"
        + "\n".join(f"{g}\t{m}" for g, m in zip(genes, global_best))
        + f"\n\n# Best score: {global_best_score}"
        + f"\n# Consensus: {consensus}\n"
    )

    print(f"\nResults saved to: {out_path.resolve()}")


if __name__ == "__main__":
    main()
