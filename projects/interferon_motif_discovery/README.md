🔬 Motif Discovery in Interferon-Stimulated Genes (ISGs)

📌 Project Overview

This project investigates the regulatory DNA motifs shared among interferon-stimulated genes (ISGs).

Using de novo motif discovery algorithms implemented from scratch, I aim to computationally infer common promoter elements that drive immune gene activation in response to interferon signaling.

This project is inspired by classical motif-finding studies such as the DosR regulatory analysis in Mycobacterium tuberculosis, but applied here to human immune response genes.

⸻

🧬 What Are Interferon-Stimulated Genes (ISGs)?

Interferon-stimulated genes (ISGs) are a set of genes that cells rapidly activate after receiving interferon signals, which typically occur during viral infection.

Instead of killing pathogens directly, interferons act as alarm signals that trigger a coordinated defensive gene program. The proteins encoded by ISGs implement antiviral actions such as:

- Inhibiting viral replication
- Degrading viral RNA
- Enhancing immune signaling
- Modulating cellular stress responses

In this project, I analyze promoter regions upstream of ISGs to discover shared regulatory DNA motifs (short recurring patterns) that transcription factors use to switch these genes on.

⸻

🧠 Biological Motivation

Interferon signaling activates transcription factors (such as STAT1, STAT2, and IRF family members) that bind to specific DNA motifs in promoter regions of ISGs.

Although certain consensus motifs (e.g., ISRE elements) have been experimentally characterized, regulatory sequences often vary between genes.

This project asks:

Can we computationally rediscover shared regulatory motifs across ISG promoters using motif-finding algorithms alone?

⸻

🔍 Computational Approach

Promoter regions upstream of selected ISGs are extracted and analyzed using:
- Median String algorithm
- Randomized Motif Search
- Gibbs Sampling
- k-mer frequency analysis
- Hamming distance & neighborhood generation
- Reverse complement scanning

The goal is to infer:

- Consensus motif sequences
- Profile matrices
- Motif length estimation
- Putative regulatory binding sites

⸻

📂 Project Structure

```text
interferon_motif_discovery/
├── data/
│   ├── gene_list.txt
│   └── promoters.fasta
│
├── scripts/
│   ├── run_median_string.py
│   ├── run_randomized.py
│   ├── run_gibbs.py
│   └── scan_genome.py
│
├── results/
│   ├── motifs_k8.txt
│   ├── motifs_k10.txt
│   ├── profile_matrix.txt
│   └── figures/
│
└── paper_summary.md
```

⸻

🎯 Project Goals

- Reproduce biologically meaningful regulatory motifs from ISG promoter data
- Compare deterministic and randomized motif-finding algorithms
- Evaluate motif conservation and variability
- Bridge computational genomics with immune biology

⸻

References

Schoggins et al., 2011, Nature.
Transcriptional regulation of antiviral interferon-stimulated genes, 2017.

⸻

👩‍🔬 Author

Jumanah Aljohani
B.Sc. in Medical Laboratory Sciences
Interested in genomics, computational bioinformatics, immune gene regulation, and translational biomedical research.
