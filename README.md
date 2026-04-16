# CMSubpathway

## Overview
**CMSubpathway** is an integrative computational framework developed to identify and characterize cancer-specific risk metabolic subpathways. By coupling metabolic network topology with multi-omics data, the pipeline identifies localized functional modules that drive cancer progression and influence clinical outcomes.

## Methodology
The identification process follows a four-stage workflow:
1.  **Topological Extraction**: Partitioning 84 KEGG metabolic networks into sub-units using the **LTE (Local Tightness Expansion)** algorithm.
2.  **Dysregulation Screening**: Identifying significant expression perturbations between tumor and normal tissues via **GSEA** (|NES| > 1, FDR < 0.25).
3.  **Prognostic Analysis**: Stratifying patient survival based on subpathway activity calculated through **GSVA**.
4.  **Diagnostic Classification**: Evaluating the discrimination power of candidate subpathways using **Support Vector Machine (SVM)** models.

## Key Findings (Breast Cancer)
* Identified a **12-gene core metabolic module** (primarily ADH and ALDH families) bridging pyruvate and fatty acid metabolism.
* Activity of this module is inversely correlated with lactate accumulation and significantly associated with the differentiation and exhaustion of stem-like **CD8+ T cells**.

## File Structure
* `subpathway.R`: Core identification pipeline, including network community detection, subsampling GSEA, and SVM classification.
* `analysis.R`: Downstream validation scripts for BRCA-specific datasets, single-cell analysis (Seurat), and immune-metabolic correlation studies.

## Citation
> Zhao et al. (2026). **Identification and characterization of cancer-related risk metabolic subpathways reveal their functional significance in cancer**. Harbin Medical University.
