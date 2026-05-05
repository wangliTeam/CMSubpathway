# CMSubpathway Analysis Pipeline Scripts

Identification and characterization of cancer-related risk metabolic subpathways reveal their functional significance in cancer

## Data Requirements

### Step 1: Metabolic Subpathway Identification
KEGG pathways: Auto-downloaded via `KEGGREST::keggGet()` (metabolic pathways) and TCGA expression

### Step 2: Independent Dataset Validation
TCGA-BRCA expression and GSE42568 dataset

### Step 3: Expression Atlas Validation
E-GEOD-45581 and E-MTAB-779 

### Step 4: Metabolic Reprogramming & CRISPR
GSE246231 dataset and DepMap CRISPR dataset

### Step 5: Functional Characteristics
TCGA expression and Clinical data

### Step 6: Single-Cell RNA-seq
GSE195861 raw data: (download from GEO, 10X matrix files)

### Step 7: CTC/Blood Cohort Validation
Bulk RNA-Seq of 9 breast cancer peripheral blood cohorts

## Datasets:
GSE111842; GSE109761; GSE111065; 
GSE51827; GSE55807 ; GSE67939; 
GSE75367; GSE86978; GSE41245; 
GSE51984  (download from GEO)

## Required R Packages

```r
# Core packages
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Bioinformatics packages
library(KEGGREST)
library(KEGGgraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)
library(limma)
library(survival)
library(survminer)
library(pROC)
library(e1071)
library(genefilter)
library(Biobase)

# Single-cell packages
library(Seurat)
library(harmony)
library(SingleCellExperiment)

# Visualization packages
library(ggpubr)
library(reshape2)
library(ggrepel)
library(rstatix)
library(UpSetR)
library(ggVennDiagram)

# Utility packages
library(openxlsx)
library(readxl)
library(stringr)
library(igraph)
```

## Usage Instructions

### Prerequisites
1. Install R >= 4.2.0
2. Install all required packages (see above)

### Key Results
1. Subpathway identification: Identified metabolic subpathways using LTE algorithm
2. GSEA screening: Subpathways screened by 100 random resampling (NES > 1, p < 0.05, fdr < 0.05)
3. Independent validation: Core module validated in TCGA-BRCA and GSE42568
4. Expression Atlas: Differential expression confirmed in E-GEOD-45581 and E-MTAB-779
5. Metabolic reprogramming: GSVA scores show downregulation in HFD/glucose conditions
6. CRISPR dependency: Gene essentiality scores in breast cancer cell lines
7. Survival analysis: Prognostic value across chemotherapy subgroups
8. Single-cell validation: Expression patterns across cell types in tumor microenvironment
9. Multi-cohort validation: Consistent downregulation across CTC datasets

