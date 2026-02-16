import random
from pathlib import Path
from typing import Dict, List, Tuple

FASTA_PATH = "ISG_promoters_300bp_GRCh38.fa"  # change to your actual path/name
K = 12                 # ISRE length
N = 2000               # iterations per run
RESTARTS = 50          # number of random starts
SEED = 42              # for reproducibility (set None for full randomness)


# -------------------------------
# FASTA reader
# -------------------------------
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


# -------------------------------
# Scoring
# Score = sum over columns of (t - max_base_count)
# Lower is better
# -------------------------------
def score_motifs(motifs: List[str]) -> int:
    t = len(motifs)
    k = len(motifs[0])
    total = 0
    for i in range(k):
        col = [m[i] for m in motifs]
        max_count = max(col.count(b) for b in "ACGT")
        total += t - max_count
    return total


# -------------------------------
# Profile with pseudocounts (+1 Laplace)
# Returns probabilities for A,C,G,T per position
# -------------------------------
def build_profile(motifs: List[str]) -> Dict[str, List[float]]:
    k = len(motifs[0])
    counts = {b: [1] * k for b in "ACGT"}  # pseudocounts
    for m in motifs:
        for i, ch in enumerate(m):
            counts[ch][i] += 1

    t = len(motifs)
    denom = t + 4  # because of +1 for each base
    profile = {b: [counts[b][i] / denom for i in range(k)] for b in "ACGT"}
    return profile


def kmer_probability(kmer: str, profile: Dict[str, List[float]]) -> float:
    p = 1.0
    for i, ch in enumerate(kmer):
        p *= profile[ch][i]
    return p


# -------------------------------
# Randomly sample an index proportional to weights
# -------------------------------
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


# -------------------------------
# Gibbs Sampler (single run)
# -------------------------------
def gibbs_sampler(dna: List[str], k: int, n: int) -> List[str]:
    t = len(dna)
    # Random initial motifs
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


# -------------------------------
# Main: multiple restarts
# -------------------------------
def main():
    if SEED is not None:
        random.seed(SEED)

    seqs = read_fasta(FASTA_PATH)
    genes = list(seqs.keys())
    dna = [seqs[g] for g in genes]

    global_best = None
    global_best_score = None

    for r in range(RESTARTS):
        motifs = gibbs_sampler(dna, K, N)
        sc = score_motifs(motifs)
        if global_best is None or sc < global_best_score:
            global_best = motifs
            global_best_score = sc

    # Print results aligned with gene names
    print(f"Best score: {global_best_score}")
    print("gene\tmotif")
    for g, m in zip(genes, global_best):
        print(f"{g}\t{m}")

    # Save to results file
    out_dir = Path("results")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"gibbs_k{K}_N{N}_restarts{RESTARTS}.tsv"
    out_path.write_text(
        "gene\tmotif\n" + "\n".join(f"{g}\t{m}" for g, m in zip(genes, global_best))
        + f"\n\n# Best score: {global_best_score}\n"
    )
    print(f"\nSaved: {out_path.resolve()}")


if __name__ == "__main__":
    main()
