
# The Multi-Omics Landscape of *S. aureus* Colonisation in Atopic Dermatitis 
## MSc Research Project | Maisie Varcoe



This repository contains the analysis pipeline for a multi-omics study of atopic dermatitis (AD)
skin, investigating how *Staphylococcus aureus* (SA) abundance relates to the skin microbiome,
metabolome and single-cell transcriptome, and using those relationships to build patient and
cell-type-specific genome-scale metabolic models (GEMs).

The project integrates four data types from AD skin samples:
- **Shotgun metagenomics**: microbiome species composition and functional (MetaCyc) pathway abundance
- **NMR metabolomics**: peak-area metabolite quantification
- **Single-cell RNA-seq**: cell-type-resolved gene expression
- **Clinical/patient metadata**: EASI score, TEWL, pH, lesional status, treatment history

## Repository Structure

```
.
├── Scripts/
│   ├── Multi_Omics_Integration.Rmd   # Main integration pipeline (R)
│   ├── README.md                     # Details on the integration pipeline
│   └── GEM/
│       ├── GEM_Data_Loading.ipynb    # Loads & harmonises data for metabolic modelling
│       ├── GEM_Building.ipynb        # Builds & runs the genome-scale metabolic models
│       └── README.md                 # Details on the GEM pipeline
├── renv.lock                         # R package version lockfile
├── renv/                             # renv project infrastructure (activate.R, settings.json)
├── .Rprofile                         # Auto-activates renv for this project
├── .gitignore
├── .gitattributes
└── README.md                         # You are here
```

## Pipeline Overview

The analysis runs in two stages, each documented in more detail in its own README:

1. **[`Scripts/Multi_Omics_Integration.Rmd`](Scripts/README.md)** (R): Loads and harmonises
   the microbiome, metabolomics, pathway, and single-cell data; tests whether SA abundance
   structures the multi-omic profile (PERMANOVA, Mantel tests, diversity/dispersion analysis);
   identifies SA-associated metabolites and pathways (regression + mediation analysis); runs
   pseudobulk differential expression (DESeq2) and GSEA on the single-cell data; and exports
   cell-type-level expression matrices as inputs for the GEM pipeline.

2. **[`Scripts/GEM/`](Scripts/GEM/README.md)** (Python/Jupyter): Takes the outputs of the
   integration pipeline (expression matrices, NMR metabolite data, microbiome data) and the
   [Human-GEM](https://github.com/SysBioChalmers/Human-GEM) base model, applies
   expression- and metabolomics-derived constraints per cell type and SA group, and runs flux
   balance analysis (GIMME/E-Flux + pFBA/FVA) to compare predicted metabolic flux between
   High-SA and Low-SA conditions.

See each linked README for step-by-step breakdowns, required input files, and outputs.

## Setup

### R environment (`Scripts/Multi_Omics_Integration.Rmd`)

This project uses [`renv`](https://rstudio.github.io/renv/) to pin R package versions. After
cloning the repo, open it in RStudio (or run R from the repo root) and restore the environment:

```r
renv::restore()
```

This reads `renv.lock` and installs the exact package versions used in the analysis. See
`Scripts/README.md` for the full package list and any packages that require Bioconductor or
GitHub installation (e.g. `MOFA2`, `SeuratDisk`).

### Python environment (`Scripts/GEM/`)

The GEM notebooks require `cobra`, `mygene`, `optlang`, `pandas`, `numpy` and related
packages - see `Scripts/GEM/README.md` for the full list and install instructions.
`GEM_Building.ipynb` is designed to run in Google Colab with the Human-GEM model and data
loaded from Google Drive.

## Data

Raw data files are not tracked in this repository (patient metadata, sequencing/metabolomics
data and the single-cell RNA-seq object are kept outside version control). Scripts 
reference local absolute paths - see the individual READMEs for the expected file names/formats
and where to update paths before re-running.

## Caveats

This is a small-cohort exploratory study (18 samples with all three omics modalities across
patients; 6 single-cell samples from 3 patients), so several analyses, particularly the
single-cell differential expression and any patient-level comparisons, should be read as
hypothesis-generating. Sample-size and confounding caveats specific to
each stage are noted in the relevant sub-README.
