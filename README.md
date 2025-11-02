# Age-Associated Transcriptional Changes in Immune Cell Development and Function  
**Research Period:** July 16 2024 – February 20 2025  
**Authors:** [Anshuman Garg](https://github.com/AnshumanG99), Dr. Rachel Ringquist, Dr. Ankur Singh  
**Institutions:** Westwood High School (Austin TX) & Georgia Institute of Technology (Atlanta GA)

[![R](https://img.shields.io/badge/R-4.3-blue.svg)]()
[![Seurat](https://img.shields.io/badge/Seurat-v5.1.0-green.svg)]()
[![DESeq2](https://img.shields.io/badge/DESeq2-v3.20-lightgrey.svg)]()
[![clusterProfiler](https://img.shields.io/badge/clusterProfiler-v3.20-orange.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## Overview
Aging disrupts immune balance, increasing susceptibility to infections, cancers, and autoimmune disorders.  
This study investigates **transcriptional changes in B-cell subsets** across the human lifespan using **single-cell RNA-seq (scRNA-seq)** data from **317 samples / 166 healthy donors (ages 25 – 85)**.  

The work identifies how age affects gene expression, pathway regulation, and developmental trajectories in immune cells — focusing on **Plasma** and **Switched Memory B Cells**, which are crucial for antibody production and long-term immunity.

---

<img width="8800" height="11000" alt="Anshuman Garg - Age Associated Transcriptional Changes in B Cell Development and Function 2" src="https://github.com/user-attachments/assets/583ee9b0-2530-4dec-9e97-3e52e3940517" />

---

## Methodology

| Step | Tool / Version | Description |
|------|----------------|--------------|
| **Data Acquisition** | PubMed dataset (Terekhova et al., 2023) | 317 PBMC samples from 166 healthy donors, 6 age groups (25–85 yrs) |
| **Pre-processing** | `Seurat v5.1.0` in R 4.3 | Normalization · QC · UMAP visualization |
| **Differential Expression** | `DESeq2 v3.20`, `FindMarkers()` | Identified age-dependent genes across cell types |
| **Pathway Enrichment** | `clusterProfiler v3.20` (GSEA) | Mapped dysregulated immune pathways |
| **Trajectory Analysis** | `tradeSeq v3.20` + `Slingshot` | Modeled pseudotime-based B-cell development |
| **Cross-Referencing** | Custom R pipeline | Overlapped outputs to isolate 9 core genes |

All analyses were GPU-accelerated using multithreaded parallel computing.

---

## Key Findings

### Cell Composition Shifts
Older donors (65 +) exhibited marked declines in:
- **Plasma Cells:** ↓ 32.3 %  
- **Switched Memory B Cells:** ↓ 27.8 %  
- Reduced overall antibody-producing capacity  

### Dysregulated Pathways (GSEA Results)
| Pathway | NES | FDR / padj | Effect |
|----------|-----|-------------|--------|
| Immunoglobulin-mediated immune response | −2.98 | < 0.05 | Strongly downregulated |
| B cell-mediated immunity | −2.36 | 0.021 | Downregulated |
| Adaptive immune response | −2.21 | < 0.05 | Downregulated |

### Nine Consensus Genes Associated with Aging
| Gene | Expression Trend (Aging) | Function |
|------|---------------------------|-----------|
| **FOXP1**, **ARHGAP24** | ↑ Upregulated | Inhibit B-cell activation / class switching |
| **BACH2**, **CD83**, **CD69**, **REL**, **SELL**, **IGHM**, **HVCN1** | ↓ Downregulated | Reduce differentiation · antibody production · immune signaling |

These genes converge functionally within the **Germinal Center Light Zone**, the site of high-affinity antibody selection, suggesting a transcriptional bottleneck driving immunosenescence.

---

## Implications
- Establishes a **nine-gene signature** as potential biomarkers of immune aging.  
- Identifies **therapeutic targets** to improve vaccine response and resilience in older adults.  
- Proposes **gene therapy / CRISPR** and **organ-on-chip** validation to test causality.  
- Suggests modulation strategies (e.g., inhibit FOXP1 / ARHGAP24 while boosting REL / SELL) to restore adaptive immunity.

---

## Environment & Dependencies
```bash
R >= 4.3
Seurat >= 5.1.0
DESeq2 >= 3.20
clusterProfiler >= 3.20
tradeSeq >= 3.20
tidyverse >= 3.5.1
