import requests

SERVER = "https://rest.ensembl.org"
HEADERS_JSON = {"Content-Type": "application/json"}
HEADERS_FASTA = {"Content-Type": "text/x-fasta"}


def get_gene_id_from_symbol(symbol: str) -> str:
    """
    Convert gene symbol to Ensembl Gene ID (Homo sapiens).
    """
    url = f"{SERVER}/xrefs/symbol/homo_sapiens/{symbol}?"
    r = requests.get(url, headers=HEADERS_JSON)
    r.raise_for_status()
    data = r.json()

    for entry in data:
        if entry["id"].startswith("ENSG"):
            return entry["id"]

    raise ValueError(f"No Ensembl Gene ID found for {symbol}")


def get_gene_info(gene_id: str) -> dict:
    """
    Retrieve gene coordinates and strand information from Ensembl.
    """
    url = f"{SERVER}/lookup/id/{gene_id}?expand=0"
    r = requests.get(url, headers=HEADERS_JSON)
    r.raise_for_status()
    return r.json()


def fetch_sequence(region: str, strand: int) -> str:
    """
    Retrieve genomic DNA sequence for a given region (GRCh38).
    """
    url = f"{SERVER}/sequence/region/human/{region}?coord_system_version=GRCh38"
    r = requests.get(url, headers=HEADERS_FASTA, params={"strand": strand})
    r.raise_for_status()
    return r.text


def define_promoter(start: int, end: int, strand: int, upstream: int):
    """
    Define promoter region based on strand orientation.
    """
    if strand == 1:  # Forward strand
        promoter_start = max(1, start - upstream)
        promoter_end = start - 1
        promoter_strand = 1
    else:  # Reverse strand
        promoter_start = end + 1
        promoter_end = end + upstream
        promoter_strand = -1

    return promoter_start, promoter_end, promoter_strand


# =====================================
# Your 15 Interferon-Stimulated Genes
# =====================================

gene_symbols = [
    "ISG15",
    "MX1",
    "MX2",
    "OAS1",
    "OAS2",
    "OAS3",
    "IFIT1",
    "IFIT2",
    "IFIT3",
    "RSAD2",
    "IFI6",
    "IFI27",
    "BST2",
    "IRF7",
    "STAT1"
]

UPSTREAM_LENGTH = 300  # Change to 2000 if needed

fasta_sequences = []

for symbol in gene_symbols:
    print(f"Processing {symbol}...")

    gene_id = get_gene_id_from_symbol(symbol)
    gene_info = get_gene_info(gene_id)

    chrom = gene_info["seq_region_name"]
    start = gene_info["start"]
    end = gene_info["end"]
    strand = gene_info["strand"]

    promoter_start, promoter_end, promoter_strand = define_promoter(
        start, end, strand, upstream=UPSTREAM_LENGTH
    )

    region = f"{chrom}:{promoter_start}..{promoter_end}"

    fasta = fetch_sequence(region, promoter_strand)

    # Replace FASTA header with gene symbol only
    lines = fasta.strip().split("\n")
    lines[0] = f">{symbol}"
    fasta_sequences.append("\n".join(lines))


# Save FASTA file
output_filename = "ISG_promoters_300bp_GRCh38.fa"

with open(output_filename, "w") as f:
    f.write("\n".join(fasta_sequences) + "\n")

print("\nDone.")
print("Promoter sequences saved to:", output_filename)
print("Total promoters extracted:", len(fasta_sequences))
