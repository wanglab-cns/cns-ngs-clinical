## Clinicogenomic Landscape of CNS Tumours: A Retrospective NGS-Based Study from a Canadian Quaternary Center

**Author:** [Farnoosh Abbas Aghababazadeh](https://github.com/RibaA)

**Contact:** farnoosh.abbasaghababazadeh@uhn.ca

---

## Description

This repository provides a reproducible computational framework for investigating the clinicogenomic landscape of central nervous system (CNS) tumours using harmonized clinical and next-generation sequencing (NGS) data.

The project supports integrative analyses of molecular alterations, clinical characteristics, and survival outcomes in a retrospective CNS tumour cohort from a Canadian quaternary care center.

In addition to analyses of the overall cohort, patients are stratified into three clinically and molecularly relevant subgroups:

- IDH-wildtype glioblastoma 
- IDH-mutant CNS tumours 
- IDH-wildtype non-glioblastoma CNS tumours 

This subgroup-based framework enables the identification of genomic and clinical associations that may be obscured in analyses of a heterogeneous CNS tumour cohort.

Clinical and molecular datasets are organized within a `MultiAssayExperiment` object to support standardized, reproducible, and modular downstream analyses.

---

## Overview

The `cns-ngs-clinical` framework supports reproducible clinicogenomic analyses of CNS tumours using harmonized clinical and next-generation sequencing (NGS) data.

The repository includes workflows for:

- Clinical and molecular data integration
- Harmonization of histological and molecular annotations
- IDH-based molecular subgroup classification
- Mutation frequency profiling and visualization
- Clinicogenomic association analyses
- Overall and subgroup-specific survival analyses
- Co-mutation and mutual-exclusivity analyses
- Generation of publication-ready figures and tables

Analyses may be performed across the complete cohort or separately within subgroups.

The workflow is modular, allowing individual analyses to be run independently or as part of a complete analytical pipeline.

---

## Analytical Modules

- Data integration and harmonization
- Mutation profiling and visualization
- Clinical association analyses (UV and MV)
- Co-Mutation Analysis

---

## Data Structure

Processed datasets are stored within a `MultiAssayExperiment` object containing mutation matrices and harmonized clinical metadata, including histology, IDH status, subgroup classification, and survival information.

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

- IDH status is derived from available IDH1 and IDH2 mutation data.
- Histology is harmonized before subgroup classification.
- Analyses are performed for the overall cohort and for patients with IDH-wildtype glioblastoma, IDH-mutant tumours, and IDH-wildtype non-glioblastoma tumours.
- Subgroup results should be interpreted considering sample size and statistical power.
