# Multi-Omics Data Integration 


Multi-omics data integration exploring metabolite-microbiome-host interactions in Atopic Dermatitis between lesional and non-lesional skin from a clinical dataset.

## Requirements

### Local Development
- macOS Sonoma 14.7.1 (Apple Silicon)
- R version 4.4.3

### HPC Analysis (SingleR annotation)
- Linux Ubuntu 22.04 (KCL CREATE HPC)
- R version 4.5.1

### Key R Packages
- Seurat (v5.4.0) — single cell analysis
- SingleR (v2.8.0) — cell type annotation
- zellkonverter (v1.16.0) — h5ad file conversion
- harmony (v1.2.4) — batch correction
- BiocParallel (v1.40.2) — parallelisation

### Environment
This project uses `renv` for package management. To restore the exact environment:
```r
renv::restore()
```
