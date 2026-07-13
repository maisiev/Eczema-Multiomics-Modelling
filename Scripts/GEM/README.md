
# GEM: Genome-Scale Metabolic Modelling Pipeline

This folder contains the two notebooks used to build patient-derived, cell-type-specific
genome-scale metabolic models (GEMs) for skin keratinocytes (KER1) and T cells, constrained
by multi-omics data (NMR metabolomics, microbiome abundance, single-cell RNA-seq expression
and pathway/metadata tables) and grouped by *S. aureus* (SA) abundance status (High vs. Low).

The pipeline uses the [Human1/Human-GEM](https://github.com/SysBioChalmers/Human-GEM) model
as a base and applies **E-Flux** expression constraints plus NMR-derived
exchange constraints, then runs flux balance analysis (FBA/pFBA/FVA) to compare metabolic
flux between High-SA and Low-SA conditions.

## Notebooks

| Notebook | Purpose |
|---|---|
| `GEM_Data_Loading.ipynb` | Loads and harmonises all raw data sources (metabolomics, microbiome, single-cell expression, metadata, pathway scores), aligns samples across modalities, assigns SA status, maps metabolites to Human-GEM IDs and packages everything into `gem_input_data.pkl` file. |
| `GEM_Building.ipynb` | Loads the Human-GEM model and the packaged `.pkl` data, applies metabolite/expression constraints per cell type and SA group, runs GIMME/E-Flux + pFBA/FVA and produces pooled (High vs. Low) flux comparisons, pathway aggregation, visualisations and gene-driver/novelty tables. |

Run them in order: **Data Loading → Building**.

## Pipeline Overview

### 1. `GEM_Data_Loading.ipynb`
1. **Load Tier 1 data**: master table (metadata + pathway scores), single-cell bridge file,
   high-confidence metabolite/pathway lists, NMR metabolomics and microbiome species tables.
2. **Load single-cell expression**: mean expression per cell type (overall, High-SA, Low-SA),
   per-sample expression for KER1 and T-cell populations and matching cell counts.
3. **Harmonise samples**: cleans indices, filters to samples present in the scRNA-seq bridge
   and intersects all modalities down to a common, manually-curated sample list
   (6 samples across 3 patients in the current version).
4. **Assign SA status**: Low (< 20% abundance), High (> 49%), samples in between are excluded.
5. **Map metabolites**: links NMR metabolite names to Human-GEM metabolite IDs (`MAM...`) via
   a KEGG ID lookup against `metabolites.tsv`.
6. **Package and save**: bundles metadata, NMR, microbiome, pathways, expression matrices,
   metabolite mappings and sample lists into a single dictionary and pickles it to
   `gem_input_data.pkl`.

### 2. `GEM_Building.ipynb`
1. **Load the Human-GEM model** (SBML `.xml`) via COBRApy and the packaged input data.
2. **Define exchange reaction mapping**: links each measured metabolite to its Human-GEM
   exchange reaction ID and uptake/secretion direction.
3. **Map gene symbols → Ensembl IDs** via `mygene`, to align expression data with model GPRs.
4. **GPR-aware weighting**: recursively evaluates gene-protein-reaction (GPR) rules
   (AND → min, OR → max) to get a per-reaction expression value, used for GIMME weighting
   and E-Flux bound scaling.
5. **Apply constraints**:
   - Essential open exchanges (water, O2, ions, vitamins) with a small "trickle" flux for
     everything else.
   - NMR-derived uptake/secretion bounds, scaled by a glucose-calibrated scale factor.
   - Skin-pH context (differing proton exchange bounds for High vs. Low SA).
6. **Metabolic task validation**: checks that biomass, ATP, glucose uptake, etc. are
   achievable per cell type before proceeding.
7. **Run GIMME / E-Flux + pFBA** on pooled (group × cell-type) models, then FVA where needed.
8. **Aggregate to pathway (subsystem) level**, compute High-vs-Low flux differences, and
   export:
   - Flux matrices and summary tables (`FluxMatrix_*.csv`, `Summary_*.csv`)
   - Pathway-level heatmaps and polarisation plots
   - Per-pathway reaction detail CSVs and metabolic maps
   - Gene-driver tables and bar plots
   - A PubMed-based "novelty screen" cross-referencing driver genes against
     atopic dermatitis / *S. aureus* literature hits.

## Requirements

`GEM_Data_Loading.ipynb`:
- Python 3, `pandas`, `numpy`, `scanpy`, `pickle`

`GEM_Building.ipynb` (designed to run in Google Colab, with Google Drive mounted):
- `cobra`, `mygene`, `optlang`, `highspy`, `glpk-utils` (installed at the top of the notebook)
- `pandas`, `numpy`, `matplotlib`, `seaborn`, `networkx`, `biopython` (`Entrez`, for the PubMed
  novelty screen)
- A local/Drive copy of the Human-GEM SBML model (`Human-GEM.xml`)

Install locally with, e.g.:
```bash
pip install cobra mygene optlang highspy pandas numpy matplotlib seaborn networkx biopython scanpy
```
(`glpk-utils` also needs to be installed system-side, e.g. `apt-get install glpk-utils` on Debian/Ubuntu.)

## Data Requirements

Both notebooks currently reference **hard-coded absolute paths** (a local `Desktop/Research_Project/...`
path for Data Loading, and a Google Drive path for Building). Before re-running, update:
- `BASE` in `GEM_Data_Loading.ipynb` to point at your local `Data/` folder containing the CSV/`.h5ad` inputs.
- `GEM_PATH`, `project_folder`, and `DATA_PATH` in `GEM_Building.ipynb` to point at your Human-GEM
  model file, output folder, and the `.pkl` produced by the Data Loading notebook.

Expected raw inputs (Data Loading `Data/` folder):
`master_table_final.csv`, `scrna_bridge.csv`, `high_conf_metabolites.csv`, `high_conf_pathways.csv`,
`nmr_matched.csv`, `species_matched.csv`, `mean_expr_per_celltype.csv`, `mean_expr_high_sa.csv`,
`mean_expr_low_sa.csv`, `mean_expr_KER1_by_sample.csv`, `mean_expr_TCELL_by_sample.csv`,
`cell_counts_KER1.csv`, `cell_counts_TCELL.csv`, `adata_qc.h5ad`, `metabolites.tsv`.

## Output

Running `GEM_Building.ipynb` writes results under `<project_folder>/Pooled_PerCellType/`, including:
- `FluxMatrix_<CellType>.csv`, `Summary_<CellType>.csv`
- `Plots/` (pathway heatmaps, polarisation heatmaps)
- `PathwayDetails/` (per-reaction long-form CSVs)
- `PathwayMaps/` (metabolite–reaction network diagrams for top pathways)
- `PathwayGeneDrivers/` (gene-level driver tables and bar plots, plus PubMed novelty scores)

## Notes / Caveats
- The current pooled comparison averages 4 High-SA samples (2 patients) against 2 Low-SA
  samples (1 patient) - High vs. Low is therefore confounded with patient identity in the
  Low group, and results should be interpreted accordingly.
- Sample size is small (6 samples / 3 patients total); this pipeline is intended as a
  proof-of-concept / exploratory analysis rather than a definitive comparison.
