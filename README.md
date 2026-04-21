# cns-ngs-clinical

**Authors:** [Farnoosh Abbas Aghababazadeh](https://github.com/RibaA)

**Contact:** [farnoosh.abbasaghababazadeh@uhn.ca](farnoosh.abbasaghababazadeh@uhn.ca)

**Description:** Clinicogenomic landscape of CNS tumours: a retrospective NGS-based study from a Canadian quaternary center

This repository provides a reproducible analysis framework for integrating clinical metadata and next-generation sequencing (NGS) mutation data in central nervous system (CNS) tumors.

The workflow supports mutation profiling, clinical association testing, and survival analyses using harmonized data structures and modular R scripts.

Data are organized within a `MultiAssayExperiment` object, enabling consistent alignment of mutation matrices and clinical variables across analyses.

----

## Overview

The `cns-ngs-clinical` framework supports:

- Integration of clinical and mutation data
- Construction of binary and oncoprint mutation matrices
- OncoPrint visualization 
- Mutation–clinical association analysis
- Clinical and mutation-based survival modelin
- IDH-stratified subgroup analyses
- IDH-stratified subgroup analyses
- Automated export of publication-ready tables and figures

The pipeline is modular, allowing individual analyses to be run independently or as part of a full workflow.

-----

**Analytical Modules**

The repository is organized into modular scripts:

- **Data Integration**
  - Clinical and mutation harmonization
  - Construction of `MultiAssayExperiment` objects

- **Visualization**
  - OncoPrint generation using `ComplexHeatmap`

- **Association Analyses**
  - Mutation–clinical associations (Fisher / logistic regression)
  - Binary outcome modeling (short vs long survival)

- **Survival Analysis**
  - Gene-level Cox proportional hazards models
  - Clinical variable survival associations
  - Multivariable and stratified models (IDH WT / Mut)

- **Co-Mutation Analysis**
  - Pairwise gene co-occurrence and mutual exclusivity

---

## Set Up

### Prerequisites

This project uses **Pixi** for environment and dependency management.

If you haven't installed it yet, [follow these instructions](https://pixi.sh/latest/)

### Installation

```bash
# Clone repository
git clone https://github.com/wanglab-cns/cns-ngs-clinical.git
cd cns-ngs-clinical

# Install environment
pixi install
```

## Repository Structure

```
cns-ngs-clinical/
├── data/               # Input datasets and intermediate objects
├── result/             # Analysis outputs (tables, figures)
├── scripts/            # Modular R analysis scripts
├── docs/               # Documentation (MkDocs)
├── pixi.toml           # Environment specification
├── mkdocs.yml          # Documentation configuration
└── README.md           # Project overview
```



