[README.md](https://github.com/user-attachments/files/27042714/README.md)
# CMSubpathway Analysis Pipeline

Complete reproducible analysis pipeline for metabolic subpathway identification and validation in breast cancer.

## Software Environment

- **R version**: >= 4.2.0
- **Core packages**:
  - Bioconductor: `KEGGgraph`, `clusterProfiler`, `org.Hs.eg.db`, `GSVA`, `limma`, `survival`, `survminer`, `GEOquery`
  - CRAN: `igraph`, `e1071`, `pROC`, `reshape2`, `ggplot2`, `ggpubr`, `patchwork`, `data.table`, `dplyr`, `tidyr`, `openxlsx`, `stringr`, `pheatmap`
  - Single-cell: `Seurat`, `harmony`, `CytoTRACE` (optional)

## Installation

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("KEGGgraph", "clusterProfiler", "org.Hs.eg.db", "GSVA", "limma", "survival", "survminer", "GEOquery"))
install.packages(c("igraph", "e1071", "pROC", "reshape2", "ggplot2", "ggpubr", "patchwork", "data.table", "dplyr", "tidyr", "openxlsx", "stringr", "pheatmap", "Seurat", "harmony"))
```

## Execution Order

| Step | Script                                     | Description                                                                                                                                                                           |
| ---- | ------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 00   | `00_Download_Raw_Data.R`                   | Data source documentation and download instructions                                                                                                                                   |
| 01   | `01_Metabolic_Subpathway_Identification.R` | KEGG metabolic pathway download, network construction, subpathway identification (LTE algorithm)                                                                                      |
| 02   | `02_Differential_Expression_and_GSEA.R`    | Differential expression analysis and GSEA with 100 random resampling permutations                                                                                                     |
| 03   | `03_Prognostic_Assessment_GSVA.R`          | GSVA/ssGSEA scoring and survival analysis (log-rank test)                                                                                                                             |
| 04   | `04_Classification_Assessment_SVM.R`       | SVM classification and AUC evaluation for each subpathway                                                                                                                             |
| 05   | `05_BRCA_Validation_GSE42568.R`            | Independent BRCA dataset validation and Subpathway-CorSP benchmark                                                                                                                    |
| 06   | `06_Core_Module_MultiOmics_Validation.R`   | Core metabolic module (12 genes) validation across transcriptome, proteome, and metabolome                                                                                            |
| 07   | `07_Peripheral_Blood_Validation.R`         | Core module validation in 9 peripheral blood cohorts (GSE111842, GSE109761, GSE111065, GSE51827, GSE55807, GSE67939, GSE75367, GSE86978, GSE41245) with ssGSEA-tumor cell correlation |
| 08   | `08_Core_Module_Functional_Analysis.R`     | Immune infiltration analysis, chemotherapy subgroup survival, pathway enrichment                                                                                                      |
| 09   | `09_scRNA_Seq_Analysis.R`                  | scRNA-seq QC, clustering, cell annotation, SiPSiC scoring, CD8 T cell subtypes, CytoTRACE                                                                                             |
| 10   | `10_CRISPR_Analysis.R`                     | CRISPR-Cas9 dependency analysis for core module genes in breast cancer cell lines                                                                                                     |
| 11   | `11_Supplementary_Validation.R`            | Expression Atlas validation and fat interference experiment analysis                                                                                                                  |

<br />

# Reproducibility Notes

- All random processes use fixed seeds (`set.seed(123)`)
- All file paths are relative to project root directory

