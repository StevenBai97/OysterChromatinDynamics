# Oyster Chromatin Dynamics

## Overview

This repository contains the analysis code associated with the research article:

**Evolutionary innovation within conserved gene regulatory networks underlying biomineralized skeletons in Bilateria**

This study integrates transcriptomics, chromatin accessibility profiling, histone modification analysis, and comparative evolutionary genomics to investigate the gene regulatory mechanisms underlying biomineralization in animals.

Using the bivalve mollusc *Crassostrea nippona* as a model, we explore how conserved transcription factors and dynamic chromatin remodeling contribute to the diversification of gene regulatory programs during larval and adult shell formation.

Our findings reveal a biphasic regulatory program in which stage-specific chromatin landscapes enable conserved upstream regulators to control distinct sets of downstream biomineralization genes, providing insights into the evolution of biomineralized skeletons across Bilateria.

## Main Analyses

The study includes the following computational analyses:

1. **Genome assembly and annotation improvement:** Refinement of the Crassostrea nippona reference genome assembly and improvement of structural and functional gene annotations to support downstream multi-omics analyses.

2. **Transcriptomic analysis:** Gene expression profiling across developmental stages and adult mantle tissue, including the identification of stage-specific biomineralization genes.

3. **Chromatin accessibility analysis:** ATAC-seq analysis to characterize accessible chromatin regions and their dynamic changes during shell development.

4. **Histone modification analysis:** CUT&Tag analysis of H3K4me3, H3K27ac, and H3K27me3 to investigate stage-specific epigenetic regulation.

5. **Chromatin state dynamics:** Integration of chromatin accessibility and histone modification profiles to characterize regulatory states and their transitions between larval and adult stages.

6. **Gene regulatory network analysis:** Reconstruction and comparison of regulatory networks associated with larval and adult shell formation.

7. **Comparative evolutionary analysis:** Cross-species comparisons of transcription factors and their expression dynamics to investigate the evolutionary conservation and divergence of biomineralization regulatory programs.

## Code and Reproducibility

These scripts support the computational analyses described in the manuscript.

Some analyses require external sequencing datasets, reference genomes, gene annotations, and software dependencies. Users should configure input and output paths and install the relevant dependencies before running the scripts.

For detailed information on datasets, analytical procedures, parameters, and software versions, please refer to the Materials and Methods section of the published article.

## Data Availability

The sequencing datasets and associated resources generated in this study are available through the following repositories:

- **NCBI BioProject:** PRJNA1267763 and PRJNA1268133
- **ProteomeXchange:** PXD064366
- **Improved genome annotation of *Crassostrea nippona*:** [Figshare](https://doi.org/10.6084/m9.figshare.29204657.v1)

Additional information on data availability is provided in the published article.

## Citation

If you use the code from this repository, please cite:

Bai Y, Min Y, Liu S, et al. (2026). **Evolutionary innovation within conserved gene regulatory networks underlying biomineralized skeletons in Bilateria.** *Molecular Biology and Evolution*, 43(2), msag019.

DOI: https://doi.org/10.1093/molbev/msag019
