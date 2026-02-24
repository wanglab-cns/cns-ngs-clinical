# cns-ngs-clinical

**Authors:** [Farnoosh Abbas Aghababazadeh](https://github.com/RibaA)

**Contact:** [farnoosh.abbasaghababazadeh@uhn.ca](farnoosh.abbasaghababazadeh@uhn.ca)

**Description:** Clinical, Mutation, and Survival Analysis Framework for CNS Tumors

This repository provides a reproducible computational framework for analyzing somatic mutation profiles, clinical variables, and overall survival outcomes in central nervous system (CNS) tumors.

The workflow integrates mutation matrices and harmonized clinical metadata within structured `SummarizedExperiment` objects to enable modular, transparent, and publication-ready analyses.

--------------------------------------

## Project Overview

The `cns-ngs-clinical` framework supports:

- OncoPrint visualization of somatic alterations
- Mutation–clinical association testing
- Kaplan–Meier survival analysis
- IDH-stratified subgroup analyses
- Automated PDF figure generation
- Fully reproducible execution using **Pixi**

The design emphasizes clarity, reproducibility, and extensibility for translational cancer genomics research.

**Analytical Modules**

- Modular pipeline implemented in **R** 
- Dependency management via **Pixi**
- Configuration-driven execution to support multiple cohorts and analyses

---

## Set Up

### Prerequisites

Pixi is required to run this project.
If you haven't installed it yet, [follow these instructions](https://pixi.sh/latest/)

### Installation

```bash
# Clone the repository
git clone https://github.com/wanglab-cns/cns-ngs-clinical.git
cd cns-ngs-clinical

# Install dependencies via Pixi
pixi install
```

## Repository Structure

```
cns-ngs-clinical/
├── data/               # Processed input objects
├── result/             # Generated outputs
├── scripts/            # Modular analysis scripts
├── docs/               # MkDocs documentation source
├── pixi.toml           # Environment specification
├── mkdocs.yml          # Documentation configuration
└── README.md           # Project overview
```



