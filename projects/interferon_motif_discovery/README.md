### 🔬 Motif Discovery and ISRE Enrichment in Interferon-Stimulated Genes (ISGs)

### 📌 Project Overview

This project investigates regulatory DNA motifs within proximal promoter regions of interferon-stimulated genes (ISGs).

Using both de novo motif discovery algorithms and targeted motif scanning approaches, the aim is to computationally identify and validate interferon-responsive elements in human ISG promoters.

### 🧬 Biological Background

Interferon-stimulated genes (ISGs) are a group of genes that are rapidly activated in response to interferon signaling, typically during viral infection.

Rather than directly eliminating pathogens, interferons act as signaling molecules that initiate a coordinated transcriptional defense program. The proteins encoded by ISGs contribute to antiviral immunity through mechanisms such as:
	
  - Inhibition of viral replication
  - Degradation of viral RNA
  - Amplification of immune signaling pathways
  - Regulation of cellular stress responses

At the molecular level, interferon signaling activates transcription factors such as STAT1, STAT2, and IRF9. These factors bind to specific regulatory DNA sequences within gene promoters, particularly the interferon-stimulated response element (ISRE).

The canonical ISRE consensus sequence is:

   ```text
AGTTTCNN TTTC
```

### 🔍 Dataset

- 15 curated ISGs from Schoggins et al. (2011)
- 15 randomly selected non-ISG control genes
- 300 bp upstream promoter regions extracted from GRCh38 (Ensembl REST API)

### 🧠 Computational Strategy

1.	De novo motif discovery
   
- MEME Suite
- Gibbs Sampling (k = 12)
- Randomized Motif Search (k = 12)

2.	Targeted ISRE motif scanning

- Consensus: AGTTTCNN TTTC
- Mismatch tolerance
- Bidirectional strand scanning

3.	Consensus sequence and PWM construction
     
4.	Enrichment analysis
    
    - ISGs vs control genes
    - Fisher’s Exact Test
   
 ### 🧪 Key Results
 
 - MEME identified an ISRE-like motif within proximal promoters.
   
 - Gibbs sampling converged to the consensus:
   
```text
AGTTTCAGTTTC
```
- Randomized Motif Search independently recovered motifs containing the conserved AGTTTC core.
   
 - ISRE scanning detected strong matches (≤2 mismatches) in:
      - 9/15 ISGs (60%)
      - 3/15 control genes (20%)
        
	- Fisher’s Exact Test:
      - Odds Ratio = 6.0
      - p-value = 0.060

Although statistical significance (α = 0.05) was not reached, a strong enrichment trend was observed.


### 📊 Interpretation

ISRE-like elements are enriched in proximal ISG promoters compared to random genes, consistent with interferon-mediated transcriptional regulation.

Multiple independent computational approaches converged on the same conserved regulatory motif.

⸻

### 📂 Project Structure

```text
interferon_motif_discovery/
│
├── data/
│   ├── ISG_promoters_300bp_GRCh38.fa
│   └── control_promoters_300bp.fa
│
├── scripts/
│   ├── run_gibbs.py
│   ├── run_randomized.py
│   ├── scan_isre.py
│   └── fisher_test.py
│
├── results/
│   ├── isre_scan_results.tsv
│   ├── control_isre_scan_results.tsv
│   ├── gibbs_results.tsv
│   ├── randomized_results.tsv
│   ├── motif_logo.png
│   └── project_summary.md
│
└── README.md
```

### References

1. Schoggins, J.W. et al. (2011). 
A diverse range of gene products are effectors of the type I interferon antiviral response. 
Nature, 472, 481–485. https://doi.org/10.1038/nature09907

### 👩‍🔬 Author

Jumanah Aljohani
B.Sc. in Medical Laboratory Sciences

Interested in genomics, computational bioinformatics, immune gene regulation, and translational biomedical research.
