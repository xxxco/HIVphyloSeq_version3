# Sequence Selection

Sequences were obtained from the [Los Alamos National Laboratory (LANL) HIV Sequence Database](https://www.hiv.lanl.gov/), downloaded Feb. 4, 2025 via the [Advanced Search Interface](https://www.hiv.lanl.gov/components/sequence/HIV/asearch/map_db.comp). The filtering pipeline below reduced the full LANL dataset to the final analysis cohort of 518 near-full-length HIV-1 subtype B genomes (one per patient), sampled 2005–2022 in the US, Canada, and Europe.

## Data Download

LANL's search interface caps the number of records returned per query, so the full database was retrieved in several sequential batches by Se ID range, each saved as a tab-delimited `.txt` file in one folder.

### Step-by-step

1. Go to the [HIV Sequence Database](https://www.hiv.lanl.gov/), click **SEARCH** in the top menu, then **Advanced Search** (or go directly to the [Advanced Search Interface](https://www.hiv.lanl.gov/components/sequence/HIV/asearch/map_db.comp)). This is a different tool from the default "Search DB" page — it lets you choose exactly which database fields to return, which the plain search does not.
2. The page shows one box per LANL database table (sequence data, genome coordinates, publication linkage, cluster annotation, patient metadata, sample metadata, accession numbers, etc.), each with individual field checkboxes. Select the fields needed to reproduce the filtering criteria below: sequence and length, HXB2 genome-coordinate start/stop, subtype, country, sampling year, QC flag, patient-linkage ID, days-from-sample fields (for tie-breaking), and GenBank accession. Free-text fields not used by the filtering criteria (comments, journal titles, clinical fields like viral load or CD4/CD8 counts, etc.) aren't needed.
3. Leave all search-criteria fields blank so the query returns every record in the database — apply no subtype/country/year restriction at download time; all filtering happens afterward, per the criteria below.
4. Submit the search, then export/download the results as tab-delimited text. Because LANL caps records per query, repeat the search restricting by `Se ID` range (or another field) to pull the full database in batches, saving each batch as its own `.txt` file in one folder for the filtering scripts to read.

Each downloaded file has two header lines (`Number of records retrieved: N`, then a blank line) before the tab-delimited data, which the filtering scripts skip.

## Filtering Steps

Counts down to "One per patient" are of **database records** (LANL emits one row per sequence–publication link); from the de-duplication step onward they are of **unique sequences**.

| Step | Criterion | Remaining (n) | Excluded (n) |
|------|-----------|---------------|--------------|
| Full LANL dataset | — | 2,727,128 | — |
| Subtype filter | HIV-1 subtype B only | 1,361,074 | 1,366,054 |
| Geographic filter | US, Canada, or one of 19 European countries | 952,093 | 408,981 |
| Temporal filter | Sampling year 2005–2024 | 450,151 | 501,942 |
| Quality control | Passed LANL QC (`Problematic Sequence == "0"`) | 409,263 | 40,888 |
| Length filter | Sequence length > 7,000 bp | 23,670 | 385,593 |
| Sequence-ID dedup | One row per `Se ID` (removes duplicate publication rows) | 11,520 | 12,150 |
| Patient linkage | Sequence linked to a patient ID (`PAT id(SSAM)` non-missing) | 11,334 | 186 |
| One per patient | Earliest sample retained; ties broken by lowest `Days from first Sample`, then lowest `Accession` | 610 | 10,724 |
| LANL re-check | Not flagged problematic by LANL after the Feb. 2025 download (re-checked 2026-10-19) | 527 | 83 |
| Inversion check | No inverted (reverse-complemented) genome segment | 525 | 2 |
| Duplicate patients | Same virus deposited under two patient IDs removed | 518 | 7 |

**Final analysis set: n = 518**

The set of 525 is what the phylogenies were built from; the 7 duplicate patients are removed immediately before clustering, so the trees are pruned rather than rebuilt.

### Geographic filter — 21 eligible countries

United States, Canada, Portugal, Spain, France, Belgium, Luxembourg, Switzerland, Italy, Germany, Austria, Czech Republic, Slovenia, Denmark, Poland, Slovakia, Hungary, Croatia, Sweden, Norway, United Kingdom.

Eleven of these are represented in the final cohort: United States (340), Spain (64), Germany (35), France (26), Sweden (23), Canada (12), United Kingdom (10), Belgium (3), Slovenia (3), Italy (1), Switzerland (1).

## Rationale for Criteria

- **Subtype B**: Restricts analysis to a single subtype to reduce confounding in phylogenetic inference; subtype B predominates in North America and Western/Central Europe.
- **US, Canada & Europe**: A high-income-country sampling frame restricted to countries where subtype B circulates as the dominant clade, and where the primary mode of transmission is sexual rather than injecting drug use.
- **2005–2024**: The eligibility window. The resulting cohort spans 2005–2022, as no eligible sequence was sampled after 2022.
- **LANL QC**: Removes sequences flagged for excessive ambiguous nucleotides, frameshifts, hypermutation, or poor HXB2 alignment — artifacts that distort phylogenetic inference and genetic distance estimation.
- **Length > 7,000 bp**: Retains near-full-length genomes. Length is measured from the first to the last position covered after mapping to HXB2, so a genome meeting this criterion may still lack complete coverage of individual genes.
- **Sequence-ID dedup**: A single sequence can appear on multiple rows because LANL emits one row per linked publication; this collapses those to one row per `Se ID`.
- **Patient linkage**: Sequences that cannot be tied to a patient ID cannot be deduplicated at the patient level and are excluded.
- **One per patient**: Prevents over-representation of individuals with multiple longitudinal sequences; the earliest available sample is kept to approximate incident/baseline virus.
- **LANL re-check**: LANL's QC annotations are updated over time; 83 sequences that passed QC at download were flagged afterward and are excluded to match the current curation.
- **Inversion check**: Two genomes contain a reverse-complemented segment that LANL has not flagged; these distort both alignment and distance estimation.
- **Duplicate patients**: Seven sequences are re-deposits of a virus already present under a different LANL patient ID. Six were identified by a shared GenBank subject code, one by an identical genome from the same year and country with no subject code available.

## Subregions Analyzed

Twelve commonly deposited HIV-1 genomic subregions were extracted from the final 518 sequences. Extraction maps each sequence to HXB2 by pairwise alignment, so insertions and deletions are handled correctly; a sequence is retained for a subregion only if it fully spans that subregion's HXB2 coordinates.

| Subregion | HXB2 Start | HXB2 End | Length (bp) | Final N |
|-----------|-----------|---------|-------------|---------|
| Vpu | 6,062 | 6,310 | 249 | 518 |
| Vpr | 5,559 | 5,850 | 292 | 518 |
| Protease | 2,253 | 2,549 | 297 | 516 |
| Tat (exon 1 + exon 2) | 5,831 / 8,379 | 6,045 / 8,469 | 306 | 516 |
| Vif | 5,041 | 5,619 | 579 | 518 |
| Nef | 8,797 | 9,417 | 621 | 478 |
| Integrase (p31) | 4,230 | 5,096 | 867 | 518 |
| Reverse Transcriptase (p51/RT) | 2,550 | 3,869 | 1,320 | 517 |
| Gag | 790 | 2,292 | 1,503 | 434 |
| Env | 6,225 | 8,795 | 2,571 | 515 |
| Pol (CDS) | 2,085 | 5,096 | 3,012 | 513 |
| Gag-Pol | 790 | 5,096 | 4,307 | 434 |

Lengths are inclusive of both endpoints (`end − start + 1`); *tat* is spliced, so its length is the sum of its two exons.

These regions correspond to genomic portions frequently sequenced for drug resistance testing, vaccine studies, and immune-escape analyses. They were selected as the most commonly deposited regions in the LANL database: of seventeen candidate regions assessed, those shorter than 100 bp (*rev* exon 1, 76 bp) or rarely deposited (the 5′ LTR, spanned by only 1,659 database sequences) were not analyzed.

Because a genome lacking complete coverage of a region is dropped from that region's alignment rather than from the cohort, the final N differs between regions. Every comparison reported is therefore a pairwise comparison of complete regions within near-full-length genomes.

GenBank accessions for each region are listed in [`accession_number/`](accession_number/), one CSV per region.
