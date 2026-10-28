# HIVphyloSeq V3

All analyses used publicly available sequences downloaded on Feb. 4, 2025, from the [LANL HIV Sequence Database](https://www.hiv.lanl.gov/): 518 near-full-length HIV-1 subtype B genomes (one per patient), sampled 2005–2022 in the United States, Canada, and European countries. See [`sequence_selection.md`](sequence_selection.md) for the full filtering workflow, including how to download the raw data from LANL. GenBank accessions for each region (near-full-length and the 12 subregions) are listed in [`accession_number/`](accession_number/).

---

## Pipeline Overview

```
LANL sequences (all fields, no filters, 14 batches by Se ID range)
    │
    ▼
Sequence filtering (→ 525 near-full-length subtype B genomes, one per patient)
    │
    ▼
Extract 12 subregions by pairwise alignment to HXB2
    │
    ▼
Multiple sequence alignment (DECIPHER)
    │
    ▼
Maximum-likelihood phylogeny (RAxML-NG)
    │
    ▼
Prune 7 duplicate patients (→ 518 analysis set)
    │
    ├─► ClusterPicker ──┐
    │                   ├─► Multi-metric concordance (compare_cluster.R)
    └─► HIV-TRACE ──────┘               │
                                        ▼
                      Dual_Ref_Results_full_new.tsv (fixed 3.0% reference)
                      Dual_Ref_Results_reference_grid.tsv (reference swept 0.5–4.5%)
```

---

## What changed from V2

| | V2 | V3 |
|---|---|---|
| Cohort | 496 | **518** |
| Subregions | 7 | **12** (adds *gag*, *pol*, *vpr*, *vpu*, *tat*) |
| Subregion extraction | substring by HXB2 offset | **pairwise alignment to HXB2** |
| *gag*-coverage filter | removed 114 genomes | **removed**; coverage handled per region |
| Cohort exclusions | — | LANL re-check (83), inverted segments (2), duplicate patients (7) |
| HIV-TRACE input | re-aligned to an HXB2 reference | **`--skip-alignment`** on the existing alignment |
| Reference threshold | fixed at 3.0% | fixed at 3.0% **plus a 0.5–4.5% sweep** |

The extraction change matters: V2 cut each subregion with `substr()` at fixed HXB2 offsets, which assumes every sequence is colinear with HXB2. Only 4 of 610 sequences were, so indels shifted the cut. V3 aligns each sequence to HXB2 first and cuts on the mapped coordinates.

The HIV-TRACE change matters too: without `--skip-alignment`, `bealign` trims every sequence to the span of the supplied `-r` reference. In V2 the near-full-length and *gag-pol* runs both used `HXB2_pol`, so both were silently reduced to 2,844 nt of *pol*, and the reported *gag-pol* F1 of 1.000 was *pol* compared against *pol*.

---

## Software Requirements

| Tool | Version | Notes |
|------|---------|-------|
| R | ≥ 4.5 | |
| DECIPHER | Bioconductor | alignment |
| Biostrings, pwalign | Bioconductor | FASTA I/O; HXB2 coordinate mapping |
| clevr | CRAN | V-measure |
| dplyr, readr, stringr, tibble, jsonlite, purrr, tidyr, ggplot2 | CRAN | data handling and figures |
| RAxML-NG | 1.2.2 | phylogeny |
| ClusterPicker | 1.2.5 | clustering |
| HIV-TRACE | 0.9.2 | clustering (run locally, from the command line) |

All tools were run inside a single Singularity/Apptainer container image, built to bundle every tool at the pinned versions above, so each step runs from the same reproducible environment:

```bash
singularity shell /path/to/multiple-alignment-container.sif
```

---

## Step 1 — Sequence Filtering

Filtering was run against a folder of raw LANL downloads (see [`sequence_selection.md`](sequence_selection.md) for the download procedure and full filtering criteria).

Output: one row per patient, tagged with HXB2 coordinates, ready for subregion extraction.

---

## Step 2 — Subregion Extraction

Each sequence is aligned pairwise to the HXB2 reference (K03455), and the twelve subregions plus the near-full-length genome are cut on the **mapped** coordinates (see the table in [`sequence_selection.md`](sequence_selection.md)), one FASTA file per region.

```r
library(Biostrings); library(pwalign)

aln <- pairwiseAlignment(subject = hxb2, pattern = seq, type = "global-local")
# map each HXB2 position to its position in the sequence, then cut on the mapped range
```

A sequence is written to a region's FASTA only if it fully spans that region; a genome lacking a region is dropped from that region alone, not from the cohort. *Tat* is spliced — its two exons are extracted separately and concatenated.

---

## Step 3 — Multiple Sequence Alignment (DECIPHER)

Alignments were generated using `DECIPHER::AlignSeqs()` with iterative refinement (5 iterations, 3 refinements). Run once per region FASTA file:

```r
library(DECIPHER)
library(Biostrings)

sequences <- readDNAStringSet("path/to/region_sequences.fasta")
aligned   <- AlignSeqs(sequences, iterations = 5, refinements = 3, processors = NULL)
writeXStringSet(aligned, "DECIPHER_aligned_region_sequences.fasta")
```

---

## Step 4 — Maximum-Likelihood Phylogeny (RAxML-NG)

Run once per region. The example below is for `Env`; repeat for each of the 13 regions (12 subregions + near-full-length), adjusting `--msa` and `--prefix` accordingly.

```bash
raxml-ng --all \
    --msa    path/to/DECIPHER_aligned_Env_sequences.fasta \
    --model  GTR+G \
    --prefix path/to/tree/Env/Env_mytree_bs \
    --seed   123 \
    --bs-trees 1000 \
    --threads auto{N}
```

**Model**: GTR + Gamma-distributed rate heterogeneity. GTR+G was prespecified, not selected by a model-selection procedure.
**Bootstrap**: 1,000 nonparametric replicates; support values mapped onto the best-scoring ML tree (`.raxml.support` used downstream).

For the two largest alignments (near-full-length and *tat*) the bootstrap replicates were computed as independent array jobs with distinct seeds and pooled to exactly 1,000; the ML search and support mapping are otherwise identical.

---

## Step 5 — Transmission Cluster Identification

Before clustering, the 7 duplicate patients are pruned from every tree and alignment, giving the 518-sequence analysis set.

### ClusterPicker

Run once per region × bootstrap threshold × genetic distance threshold. Bootstrap thresholds: 70, 90, 99. Genetic distance thresholds: 0.5%–4.5% in 0.5% steps.

| Argument | Value | Description |
|----------|-------|-------------|
| 1 | `DECIPHER_aligned_Env_sequences.fasta` | DECIPHER-aligned FASTA file |
| 2 | `Env_mytree_bs.raxml.support` | RAxML-NG bootstrap support tree |
| 3 | `90` | Initial bootstrap support threshold (%) |
| 4 | `90` | Final bootstrap support threshold (%) — set equal to argument 3 |
| 5 | `0.03` | Maximum pairwise genetic distance within a cluster (3%) |
| 6 | `2` | Minimum cluster size |
| 7 | `gap` | Distance method: gap positions ignored |

```bash
java -jar ClusterPicker_1.2.5.jar \
    DECIPHER_aligned_Env_sequences.fasta \
    Env_mytree_bs.raxml.support \
    90 90 \
    0.03 \
    2 \
    gap
```

### HIV-TRACE

Run once per region × genetic distance threshold, locally from the command line. **`--skip-alignment` is required**: the input is already a multiple sequence alignment, and without this flag HIV-TRACE re-aligns to the `-r` reference and silently truncates every sequence to that reference's span.

| Flag | Value | Description |
|------|-------|-------------|
| `-i` | `DECIPHER_aligned_Env_sequences.fasta` | DECIPHER-aligned FASTA file |
| `--skip-alignment` | — | treat the input as already aligned |
| `-a` | `resolve` | ambiguity handling |
| `--ambiguities-resolution` | `0.015` | sequences with >1.5% ambiguous bases are excluded |
| `-t` | `0.03` | genetic distance threshold for cluster membership |
| `-m` | `200` | minimum nucleotide overlap between two sequences |
| `-o` | `Env_hivtrace_GD030.json` | output JSON |

```bash
hivtrace \
    -i DECIPHER_aligned_Env_sequences.fasta \
    --skip-alignment \
    -a resolve --ambiguities-resolution 0.015 \
    -t 0.03 \
    -m 200 \
    -o Env_hivtrace_GD030.json
```

---

## Step 6 — Cluster Concordance Analysis

`compare_cluster.R` compares subregion clusters (ClusterPicker and HIV-TRACE, swept across bootstrap and genetic-distance thresholds) against the near-full-length reference clustering, method-matched. It reports:

| Metric | Meaning |
|--------|---------|
| `V_measure_refclustered` | **Headline metric.** V-measure scored only over sequences whose reference cluster has ≥ 2 members, so the large mutually-singleton background typical of HIV clustering doesn't saturate the score near 1. |
| `pw_recall` | Pairwise recall: fraction of reference-linked pairs the subregion also links. |
| `pw_precision` | Pairwise precision: fraction of subregion-linked pairs also linked in the reference. |
| `pw_f1` | Harmonic mean of `pw_precision` and `pw_recall`. |

For V-measure, every sequence not assigned to a cluster is given its own unique singleton label, rather than being pooled into a single "unclustered" category.

Edit the path variables at the top of `compare_cluster.R`, then:

```bash
Rscript compare_cluster.R regions.txt
```

Outputs:

- `Dual_Ref_Results_full_new.tsv` — one row per region × method × GD × bootstrap, scored against the reference held at GD 3.0%.
- `Dual_Ref_Results_reference_grid.tsv` — the same, repeated with the reference itself redefined at each threshold from 0.5% to 4.5% (the reference-threshold sensitivity analysis).

---

## File Structure

```
HIVphyloSeq_V3/
├── README.md                     # this file
├── sequence_selection.md         # LANL download instructions, filtering criteria, cohort description
├── compare_cluster.R             # multi-metric cluster concordance analysis (clevr R package)
└── accession_number/             # GenBank accessions (Se.ID, Accession), one CSV per region
    ├── accessions_fulllength.csv
    ├── accessions_Gag_Pol.csv
    ├── accessions_Pol.csv
    ├── accessions_Gag.csv
    ├── accessions_protease.csv
    ├── accessions_RT.csv
    ├── accessions_integrase.csv
    ├── accessions_Vif.csv
    ├── accessions_Vpr.csv
    ├── accessions_Vpu.csv
    ├── accessions_Tat.csv
    ├── accessions_Env.csv
    └── accessions_Nef.csv
```
