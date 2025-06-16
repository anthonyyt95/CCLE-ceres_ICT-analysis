# CCLE-CERES ICT Analysis

A comprehensive cancer genomics analysis pipeline that integrates CERES dependency scores, CCLE expression data, and ICT (ion channel and transporter) gene lists to identify therapeutic targets and characterize cancer cell line relationships.

---

## Overview

This repository provides tools for analyzing cancer cell line dependencies using data from the Cancer Cell Line Encyclopedia (CCLE) and CERES (Cancer Dependency Map) databases. The pipeline focuses on ion channel and transporters (ICT) genes and their dependencies across different cancer types.

---

## 📂 Analysis Components

### MATLAB Pipeline

- **PCA Analysis**  
  Dimensionality reduction and visualization  
  `20190427_CCLE_PCAanalysis.m:223`

- **CRISPRi Analysis**  
  Gene dependency scoring and ranking  
  `20190427_CRISPRiAnalysis.m:72-73`

- **CERES Analysis**  
  Enhanced dependency analysis with z-score normalization  
  `20201002_CERESanalysis2.m:57-58`

### R Pipeline

- **Survival Analysis**  
  Clinical outcome modeling using TCGA data  
  `20200219_SurvivalCurve.Rmd:284-287`

- **Statistical Testing**  
  Survival curve generation and risk assessment

---

## 📊 Data Sources

- **CERES**: Gene dependency scores from genome-wide CRISPR screens  
- **CCLE**: Cancer cell line expression profiles and annotations  
- **TCGA**: Clinical survival data for validation  
- **ICT Gene Lists**: Curated ion channel and transporter genes

---
