# Albanian Population STR Profile Dataset

A dataset containing **STRAF 2000-compatible autosomal STR profiles** from individuals sampled in **Albania**.

## Overview

This repository contains a tab-delimited text file with autosomal STR genotypes formatted for use with population genetics and forensic analysis software.

The dataset includes:

- Sample identifiers
- Population information (Albania)
- Genotypes for multiple autosomal STR loci
- STRAF 2000 compatible formatting

Dataset to be published once the article gets published and gets approved for public availability in STRAF

Alleles are reported according to standard forensic STR nomenclature, including microvariants (e.g., 9.3, 15.2, 28.2).

## Repository Contents

* **`Albania_structure_full/`**
  Contains input file and parameters to run STRUCTURE software. In addition, I have included output and logfiles as well for reproducibility.

* **`NewPopulations/`**
  Contains allele frequency tables for some countries not available in STRidER

* **`Albania_Population_Genetics_Comparison_Analysis.R`**
  Performs population genetic comparisons between the Albanian population and populations from other European countries avaialable in STRidER Database, plus some data which we retrieved from literature and curated ourselves, such as Kosovo, Turkish Cyprus, and selected regions of Turkey.

* **`Albania_Post_Structure_Analysis.R`**
  Performs post-processing of STRUCTURE output, including summarization and visualization of the STRUCTURE results.

* **`Albania_Region_Level_Analysis.R`**
  Explores genetic heterogeneity within the Albanian population at the regional level and tests some hypothesis related to that

* **`Albania_Regions_AMOVA_Analysis.R`**
  Performs Analysis of Molecular Variance (AMOVA) using the available Albanian metadata to assess how genetic variation is partitioned among geographic and individual-level groupings.

* **`Albania_STRUCTURE_Data_Preparation.R`**
  Prepares input files for STRUCTURE analysis starting from the raw allele data.

## Software Compatibility

The file is suitable for use with:

- STRAF 2000
- R v4.6.1

## Data Format

- Text file (.txt)
- Tab-separated values (TSV)
- One individual per row
- Diploid genotype reported as two allele columns per locus

## Population

| Population | Country |
|------------|---------|
| ALBANIA | Albania |

## Citation: TBA

## Disclaimer

This repository is intended for research and educational purposes only. Users are responsible for ensuring compliance with local regulations and ethical requirements regarding the use of human genetic data.
