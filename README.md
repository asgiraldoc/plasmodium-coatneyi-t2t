# Plasmodium coatneyi T2T Genome Analysis Pipeline

This repository contains the scripts and pipelines used for the structural and comparative genomics analyses in the publication: **"A telomere-to-telomere Plasmodium coatneyi genome exhibits subtelomeric remodeling and dispersed gene duplications in the Plasmodium vivax clade"**

## Overview
The scripts provided here reproduce the core chromosomal architecture analyses described in the study, including:
1. Whole-genome synteny visualization.
2. Annotation and extraction of telomeric arrays and centromeric regions.
3. Karyotype feature visualization across the *P. vivax* clade.
4. Identification of major structural rearrangements (fusions, translocations, inversions) using contiguous ancestral regions (CARs).

## Dependencies
Ensure the following tools are installed and available in your `$PATH`:
- [MUMemto](https://github.com/MUMemto/mumemto)
- [BUSCO](https://busco.ezlab.org/) (v5+)
- [Telociraptor](https://github.com/slimsuite/telociraptor)
- [TIDK](https://github.com/tolkit/telomeric-identifier)
- [PGGB](https://github.com/pangenome/pggb)
- [ODGI](https://github.com/pangenome/odgi)
- [SeqKit](https://bioinf.shenwei.me/seqkit/)
- [ChromSyn](https://github.com/ChromSyn)
- Python 3.x (with standard libraries)
- R (for ChromSyn)

## Repository Structure

### `01_assembly_and_synteny/`
[cite_start]Contains scripts to generate whole-genome multi-MUM synteny plots comparing assemblies (e.g., *P. coatneyi* v1 vs v2)[cite: 53, 434].

### `02_chromosome_features/`
[cite_start]Pipelines to detect telomeric repeat arrays (`5′-GGGTT(T/C)A-3′`) using Telociraptor and TIDK [cite: 58, 267][cite_start], and to identify/extract centromeric regions (GC troughs) using reference-unbiased pangenome graphs constructed with PGGB[cite: 57, 65].

### `03_structural_rearrangements/`
[cite_start]Contains `identify_rearrangements.py` (an extension of the AGORA toolkit), which analyzes synteny blocks to quantify fusions, fissions, translocations, and inversions across lineages[cite: 77, 276].

### `04_karyotype_plotting/`
[cite_start]Wrappers to process feature lists (telomeres, gaps, BUSCO genes, TIDK windows) and plot the final chromosome-scale feature tracks using ChromSyn[cite: 55, 269].
