from pathlib import Path

# Path to the 300bp promoter FASTA file
FASTA_PATH = "ISG_promoters_300bp_GRCh38.fa"  # Change if needed


# -------------------------------
# FASTA reader
# -------------------------------
def read_fasta(fp):
    sequences = {}
    name = None
    parts = []

    for line in Path(fp).read_text().splitlines():
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if name:
                sequences[name] = "".join(parts).upper()
            name = line[1:].split()[0]
            parts = []
        else:
            parts.append(line)

    if name:
        sequences[name] = "".join(parts).upper()

    return sequences


# -------------------------------
# Reverse complement function
# -------------------------------
_comp = str.maketrans("ACGT", "TGCA")

def reverse_complement(seq):
    return seq.translate(_comp)[::-1]


# -------------------------------
# ISRE consensus scoring
# Consensus: AGTTTCNN TTTC
# N positions are ignored in mismatch calculation
# -------------------------------
ISRE_LENGTH = 12

def count_mismatches_isre(window):
    mismatches = 0

    # Fixed part 1: AGTTTC (positions 0–5)
    fixed1 = "AGTTTC"
    for i, base in enumerate(fixed1):
        if window[i] != base:
            mismatches += 1

    # Fixed part 2: TTTC (positions 8–11)
    fixed2 = "TTTC"
    for j, base in enumerate(fixed2, start=8):
        if window[j] != base:
            mismatches += 1

    return mismatches


# -------------------------------
# Find best ISRE-like match
# Searches both forward and reverse strands
# -------------------------------
def find_best_isre(seq):
    best_mismatch = None
    best_hit = None  # (strand, position, sequence, mismatches)

    for strand, s in [("+", seq), ("-", reverse_complement(seq))]:
        for i in range(len(s) - ISRE_LENGTH + 1):
            window = s[i:i + ISRE_LENGTH]
            mismatches = count_mismatches_isre(window)

            if best_mismatch is None or mismatches < best_mismatch:
                best_mismatch = mismatches
                best_hit = (strand, i, window, mismatches)

    return best_hit


# -------------------------------
# Main execution
# -------------------------------
def main():
    sequences = read_fasta(FASTA_PATH)
    results = []

    for gene, seq in sequences.items():
        strand, position, match, mismatches = find_best_isre(seq)
        results.append((gene, mismatches, strand, position, match))

    # Sort by number of mismatches (ascending)
    results.sort(key=lambda x: x[1])

    print("gene\tmismatches\tstrand\tposition\tbest_match")
    for row in results:
        print("\t".join(map(str, row)))

    # Save results to file
    output_file = Path("isre_scan_results.tsv")
    output_file.write_text(
        "gene\tmismatches\tstrand\tposition\tbest_match\n" +
        "\n".join("\t".join(map(str, row)) for row in results)
    )

    print(f"\nResults saved to: {output_file.resolve()}")


if __name__ == "__main__":
    main()
