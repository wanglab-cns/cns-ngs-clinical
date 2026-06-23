## Clinicogenomic Landscape of CNS Tumours: A Retrospective NGS-Based Study from a Canadian Quaternary Center

**Author:** [Farnoosh Abbas Aghababazadeh](https://github.com/RibaA)

**Contact:** farnoosh.abbasaghababazadeh@uhn.ca

---

## Description

This repository provides a reproducible computational framework for the analysis of clinicogenomic features in central nervous system (CNS) tumours using harmonized clinical and next-generation sequencing (NGS) data.

The project supports integrative analyses of molecular alterations, clinical characteristics, and survival outcomes in a retrospective CNS tumour cohort from a Canadian quaternary care center.

Clinical and molecular datasets are organized within a `MultiAssayExperiment` object to enable standardized, reproducible, and modular downstream analyses across mutation profiling, clinical association testing, visualization, and survival modeling workflows.

---

## Overview

The `cns-ngs-clinical` framework supports reproducible clinicogenomic analyses of CNS tumours using harmonized clinical and next-generation sequencing (NGS) data.

The repository includes workflows for:

- Clinical and molecular data integration
- Mutation profiling and visualization
- Clinicogenomic association analyses
- Survival and subgroup analyses
- Co-mutation analyses
- Generation of publication-ready figures and tables

The workflow is modular, allowing analyses to be run independently or as part of a complete analytical pipeline.

---

## Analytical Modules

- Data Integration
- Visualization
- Association and Survival Analyses (UV and MV)
- Co-Mutation Analysis

---

## Data Structure

Processed datasets are stored within a `MultiAssayExperiment` object containing:

- Binary mutation matrices
- OncoPrint-compatible mutation matrices
- Harmonized clinical metadata

This structure enables consistent patient alignment across molecular and clinical analyses.

---

## Set Up

### Prerequisites

This project uses **Pixi** for environment and dependency management.

Installation instructions are available at:

https://pixi.sh/latest/

---

### Installation

```bash
# Clone repository
git clone https://github.com/wanglab-cns/cns-ngs-clinical.git

# Navigate to repository
cd cns-ngs-clinical

# Install project environment
pixi install
```

---

## Repository Structure

```text
cns-ngs-clinical/
├── data/               # Input datasets and processed objects
├── result/             # Analysis outputs (tables and figures)
├── scripts/            # Modular R analysis scripts
├── docs/               # Documentation and project notes
├── pixi.toml           # Environment specification
├── mkdocs.yml          # Documentation configuration
└── README.md           # Project overview
```

---

## Notes

- Analyses are implemented in R/Bioconductor using standard statistical and bioinformatics approaches, with IDH status determined from combined IDH1 and IDH2 mutation data.
- The repository is designed to support reproducible and extensible clinicogenomic analyses in CNS tumours.
