# Multi-Omics Integration: Skin Microbiome, Metabolomics & scRNA-seq

**File:** `Multi_Omics_Integration.Rmd`
**Author:** Maisie Varcoe

This R Markdown notebook integrates shotgun-metagenomic (microbiome + functional pathway),
NMR metabolomics, and single-cell RNA-seq data from atopic dermatitis skin samples. It tests
whether *Staphylococcus aureus* (SA) abundance structures the multi-omic profile, identifies
SA-associated molecular features, and prepares the pseudobulk expression inputs used for
downstream genome-scale metabolic modelling (GEM).

## Pipeline Overview

### Tier 1: Microbiome, Metabolomics & Pathways

| Step | What it does |
|---|---|
| **1. Data Preparation** | Loads the microbiome species table, NMR peak-area metabolomics, functional pathway abundance and the linking/metadata files (sequencing library, sample metadata, patient metadata, NMR–shotgun bridge). Extracts SA relative abundance from the microbiome table and joins everything into a single **master linking table** keyed by sample/patient IDs. |
| **2. Confounder Analysis** | Aggregates to one row per patient for visualisation, cleans clinical variables and checks whether SA abundance is confounded with age, EASI score, treatment history, TEWL or pH (scatter/box plots + Spearman correlations). Builds sample-matched, aligned matrices (NMR, microbiome, pathways, metadata) restricted to the **18 samples** present in all three modalities. Runs PERMANOVA (screening + a primary "SA adjusted for pH" model) on both pathway and metabolite composition, plus PCA/PCoA ordination coloured by SA status and pH. |
| **3. Confirm SA Structures Microbiome Composition** | Repeats PERMANOVA, partial Mantel tests (SA vs. composition, controlling for pH), beta-dispersion tests and richness/diversity LMMs across all three omics layers (microbiome, pathways, metabolites), each split into a High/Low SA group by median. Produces PCoA biplots (SA group centroids + dispersion, pH as an `envfit` vector) and a combined variance-partitioning summary plot. |
| **5. Linear Regressions - SA High/Low Associated Features** | Defines High (>35%) / Low (<10%) SA groups (mid-range excluded), then runs feature-wise linear regressions of SA status against NMR metabolites and microbial pathways (adjusting for lesional status and pH), a pathway→metabolite regression and Spearman sensitivity checks. Produces significance/high-confidence feature tables, boxplots, forest plots and a pathway–metabolite correlation heatmap. |
| **6. Mediation Analysis** | Tests whether significant pathways mediate the relationship between SA status and metabolite levels, for both the top metabolite hits and top 5 pathways, saving ACME/mediation summary tables. |

### Tier 2: scRNA-seq Mechanistic Subset

| Step | What it does |
|---|---|
| **Load & QC** | Loads the single-cell RNA-seq object, attaches metadata, computes % mitochondrial reads and filters to cells with >200 genes and <5% mitochondrial content. |
| **6. DEG for SA High vs. Low** *(Tier 2 numbering restarts from 6)* | Builds a patient-level pseudobulk (aggregating by true patient, not sample, since 3 of the 6 scRNA samples are lesional/non-lesional pairs from the same 3 individuals) and runs one DESeq2 model per cluster comparing High vs. Low SA. Includes a DEG-count vs. cluster-size diagnostic, volcano plots (patient-confounded "recurring" genes flagged separately) and curated dot plots of cluster-specific significant genes. |
| **6.5. GSEA** | Gene set enrichment analysis on the DE results, plus a pathway-recurrence diagnostic across clusters. |
| **7. GEM Input Preparation** | Normalises expression and builds mean-expression matrices per cell type, split by all-samples / High-SA / Low-SA, exporting `mean_expr_per_celltype.csv`, `mean_expr_high_sa.csv` and `mean_expr_low_sa.csv` for downstream GEM construction. |


## Requirements

R packages loaded at the top of the notebook:

```
readxl, dplyr, ggplot2, ggpubr, patchwork, phyloseq, vegan, tibble, MOFA2, ggrepel,
pheatmap, lme4, lmerTest, mediation, compositions, MultiAssayExperiment, Seurat,
reshape2, apeglm, magrittr, DESeq2, fgsea, msigdbr, tidyr, knitr, SeuratDisk
```

Most are on CRAN; the Bioconductor packages (`MOFA2`, `apeglm`, `phyloseq`, `DESeq2`, `fgsea`,
`msigdbr`, `MultiAssayExperiment`) need `BiocManager::install(...)` before first use.
`SeuratDisk` is used for reading/writing `.h5ad`/`.h5Seurat` single-cell files and is installed
from GitHub (`remotes::install_github("mojaveazure/seurat-disk")`) rather than CRAN/Bioconductor.

`MOFA2` is loaded and used to build the sample-matched objects referenced in Step 2, but no
`run_mofa()`/`create_mofa()` call appears later in this notebook - the matched matrices it
builds are reused for the PERMANOVA/regression analyses instead.

## Data Requirements

The notebook uses **hard-coded local paths** (partially redacted, e.g. `PXXXX` for the
sequencing library ID prefix) in the **Load Data** section. Before re-running, point these at
your local copies of:

- Microbiome species abundance table (CSV) and the corresponding `physeq` object (`.RData`)
- NMR peak-area metabolomics metadata (`.RData`)
- Cleaned functional pathway abundance table (CSV)
- Sequencing library file linking sequencing IDs to sample numbers (`.xlsx`)
- Shotgun sample metadata (eg_id, site, lesional status, local EASI) and patient metadata
  (age, sex, EASI, treatments) (`.xlsx`)
- NMR–shotgun sample-name linking file (CSV)
- Sample metadata with the `eg_id` linking factor (`.RData`)
- The single-cell RNA-seq object used in the Tier 2 section
- An output directory for the GEM-ready CSVs written in Step 7

## Outputs

- A results directory (set near Step 3) containing PERMANOVA/beta-dispersion/Mantel/LMM text
  reports, PCoA biplots, diversity plots, and the combined variance-partitioning plot.
- Regression/mediation result tables and figures for Steps 5–6 (significant feature lists,
  boxplots, forest plots, heatmaps).
- DESeq2 DEG tables, volcano plots, dot plots, and GSEA results for the Tier 2 Step 6 analysis.
- `mean_expr_per_celltype.csv`, `mean_expr_high_sa.csv`, `mean_expr_low_sa.csv` from Step 7,
  for use as GEM expression inputs.

## Caveats

- The scRNA-seq High/Low SA comparison is fundamentally a **3-patient** comparison (2 High,
  1 Low) once the two lesional/non-lesional biopsies per patient are correctly collapsed -
  DESeq2 can still fit a model, but there is no way to estimate true within-group variability
  for "Low" with a single biological replicate. Treat these DEGs as exploratory rather than
  validated.
- Sample sizes are small throughout (18 samples with all three Tier-1 modalities; 6 scRNA
  samples from 3 patients) — this pipeline is best read as hypothesis-generating.
