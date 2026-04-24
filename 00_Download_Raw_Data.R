# ==============================================================================
# 00_Download_Raw_Data.R
# Description: Document all raw data sources and download instructions
# ==============================================================================

# ==============================================================================
# DATA SOURCE DOCUMENTATION
# ==============================================================================
# This script documents ALL publicly available datasets used in this study.
# KEGG pathways are auto-downloaded via R packages (steps in 01_Metabolic_Subpathway_Identification.R).
# All other datasets must be manually downloaded as described below.
#
# After downloading, place files in the specified directories relative to project root.
# ==============================================================================

# --- 1. TCGA-BRCA Expression Data ---
# Source: The Cancer Genome Atlas (TCGA)
# Description: RNA-seq expression matrix for breast cancer (tumor + normal samples)
# Download: https://portal.gdc.cancer.gov/projects/TCGA-BRCA
#   - Use GDC Data Transfer Tool 
# Place at: ./data/raw/TCGA/TCGA_exp.rda
# Format: RData file containing 'exp' object (genes x samples matrix)
# Note: Sample barcodes follow TCGA convention (positions 14-15: 01=tumor, 11=normal)

# --- 2. TCGA-BRCA Clinical Data ---
# Source: TCGA via cBioPortal 
# Description: Clinical information including survival data and treatment history
# Download: https://www.cbioportal.org/study/summary?id=brca_tcga
#   - Download patient clinical data and treatment timeline
# Place at:
#   ./data/raw/TCGA/data_clinical_patient.txt
#   ./data/raw/TCGA/data_timeline_treatment.txt

# --- 3. GSE42568 (Independent BRCA Validation) ---
# Source: GEO (Gene Expression Omnibus)
# GEO Accession: GSE42568
# Description: 104 breast tumor + 17 normal samples with RNA-seq and survival data
# Download: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568
# Place at: ./data/raw/GSE42568/
# Files: expression matrix and clinical annotation

# --- 4. GSE195861 (Single-cell RNA-seq) ---
# Source: GEO
# GEO Accession: GSE195861
# Description: scRNA-seq data for breast cancer tumor microenvironment
# Download: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195861
#   - Download raw count matrices (10X Genomics format: barcodes.tsv, features.tsv, matrix.mtx)
# Place at: ./data/raw/GSE195861/

# --- 5. DepMap CRISPR Data ---
# Source: DepMap (Dependency Map) Portal
# Description: CRISPR-Cas9 gene dependency scores for breast cancer cell lines
# Download: https://depmap.org/portal/data_page/?tab=allData
#   - Download: CRISPRGeneEffect.csv (gene effect scores; Rows:ModelID, Columns:Gene) 
#   - Download: Model.csv (Metadata describing all cancer models/cell lines which are referenced by a dataset contained within the DepMap portal.)
# Place at:
#   ./data/raw/DepMap/CRISPRGeneEffect.csv
#   ./data/raw/DepMap/Model.csv

# --- 6. Expression Atlas: E-GEOD-45581 ---
# Source: EMBL-EBI Expression Atlas
# Accession: E-GEOD-45581
# Description: 40 untreated breast cancer + 5 normal control samples
# Download: https://www.ebi.ac.uk/gxa/experiments/E-GEOD-45581
# Place at: ./data/raw/ExpressionAtlas/E-GEOD-45581/

# --- 7. Expression Atlas: E-MTAB-779 ---
# Source: EMBL-EBI Expression Atlas
# Accession: E-MTAB-779
# Description: 20 untreated breast carcinoma + 22 normal control samples
# Download: https://www.ebi.ac.uk/gxa/experiments/E-MTAB-779
# Place at: ./data/raw/ExpressionAtlas/E-MTAB-779/

# --- 8. CPTAC Proteomics Data ---
# Source: CPTAC (Clinical Proteomic Tumor Analysis Consortium)
# Description: Proteomics data for breast cancer validation
# Download: https://proteomics.datacommons.cancer.gov/
#   - download from data portal
# Place at: ./data/raw/CPTAC/

# --- 9. HPA (Human Protein Atlas) ---
# Source: https://www.proteinatlas.org/
# Description: Protein expression validation for core module genes
# Download: Use HPA web interface
#   - Query genes: ADH1A, ADH1B, ADH1C, ADH4, ADH5, ADH6, ADH7, ALDH1B1, ALDH2, ALDH3A2, ALDH7A1, ALDH9A1
# Place at: ./data/raw/HPA/

# --- 10. GSE246231 (Fat Interference Experiment) ---
# Source: GEO
# GEO Accession: GSE246231
# Description: Gene expression data for fat/high-fat diet/glucose intervention experiment
# Download: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246231
# Place at: ./data/raw/GSE246231/

# --- 11. CTC/Peripheral Blood Datasets (9 BRCA cohorts) ---
# Source: GEO 
# Description: Bulk RNA-Seq of 9 breast cancer peripheral blood cohorts
# Datasets:
#   GSE111842; GSE109761; GSE111065; 
#   GSE51827; GSE55807 ; GSE67939; 
#   GSE75367; GSE86978; GSE41245; 
#   GSE51984  
# Download:
#   - GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXXXX
# Place at: ./data/raw/CTC/
# Files:
#   BR_GSE111842_gene_FPKM.xlsx
#   BR_GSE109761_gene_FPKM.xlsx
#   BR_GSE111065_gene_FPKM.xlsx
#   BR_GSE51827_gene_FPKM.xlsx
#   BR_GSE55807_gene_FPKM.xlsx
#   GSE67939_readCounts.txt.gz
#   GSE67939_annotation_file.txt.gz
#   BR_GSE75367_gene_FPKM.xlsx
#   BR_GSE86978_gene_FPKM.xlsx
#   GSE41245_processed_data.txt.gz
#   wbc_GSE51984_gene_FPKM.xlsx(control)
# Note: .xlsx files were converted from the original .txt.gz files downloaded from GEO and saved as Excel format for convenience.


# Data source documentation complete. Please download all datasets as described above.
