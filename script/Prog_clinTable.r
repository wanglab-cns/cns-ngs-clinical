##-------------------------------------------------------------------
## Script: CNS NGS OncoPrint Visualization (ComplexHeatmap)
##
## Purpose:
##   Generate publication-quality OncoPrint visualizations
##   from a curated CNS NGS mutation event matrix stored
##   in a MultiAssayExperiment object.
##
##   Standardize mutation alteration labels and enforce a
##   consistent within-cell priority ordering to ensure
##   reproducible visualization of multi-event mutations.
##
##   Integrate harmonized clinical metadata and generate
##   annotated OncoPrints for the full cohort and key
##   molecular subgroups.
##
##
## Input:
##   - result/data/mae_mut_clin.RData
##       MultiAssayExperiment object containing:
##         - mut_oncoprint assay:
##             gene × patient character matrix
##             ("" or semicolon-delimited alteration types)
##         - colData:
##             curated clinical metadata
##
##
## Output (4 figures):
##
##   1) result/oncoprint/fig1.pdf
##        OncoPrint (mutation events only; no annotations)
##
##   2) result/oncoprint/fig2.pdf
##        OncoPrint with clinical annotations (all patients)
##
##   3) result/oncoprint/fig3_wt.pdf
##        OncoPrint with annotations (IDH wild-type subset)
##
##   4) result/oncoprint/fig4_mut.pdf
##        OncoPrint with annotations (IDH mutant subset)
##
##
## Processing Overview:
##
##   1) Data Loading
##        The MultiAssayExperiment object is loaded and the
##        oncoprint assay (mutation matrix) and associated
##        clinical metadata are extracted.
##
##   2) Mutation Label Harmonization
##        Mutation labels are standardized to:
##           - SNV/Indel
##           - Amplification
##           - Fusion
##           - Deletion
##
##        Within each gene–patient cell, multiple alteration
##        types are:
##           - cleaned and deduplicated
##           - ordered using a fixed priority:
##                Fusion > Amplification > SNV/Indel > Deletion
##
##        Multi-hit events are preserved as semicolon-
##        separated entries.
##
##   3) OncoPrint Construction (ComplexHeatmap)
##        OncoPrints are generated using custom graphical
##        functions for each alteration type, enabling
##        stacked visual representation within each cell.
##
##        Features include:
##           - Removal of empty genes and samples
##           - Row-level mutation frequency barplots
##           - Consistent color mapping across alteration types
##
##   4) Clinical Annotation Processing
##        Clinical metadata are formatted for visualization:
##           - Age dichotomized (<40 vs ≥40)
##           - IDH status encoded as WT / Mut
##           - Histology harmonized and grouped
##           - WHO 2021 grade converted to I–IV
##
##        Annotation tracks are constructed and aligned
##        with the mutation matrix.
##
##   5) Annotated OncoPrint Generation
##        A full-cohort OncoPrint is generated with clinical
##        annotations displayed as top annotation tracks.
##
##   6) Subset Analyses
##        The dataset is stratified by IDH status:
##           - IDH wild-type patients
##           - IDH mutant patients
##
##        Independent OncoPrints are generated for each
##        subgroup with appropriately filtered mutation
##        matrices and clinical annotations.
##
##   7) Export
##        All OncoPrints are exported as PDF files with
##        consistent layout, legends, and formatting for
##        publication-ready visualization.
##-------------------------------------------------------------------
####################################################
## Load libraries
####################################################
library(MultiAssayExperiment)
library(ComplexHeatmap)
library(dplyr) 
library(tidyr)
library(stringr)

####################################################
## Setup directories
####################################################
dir_input <- 'result/data'
dir_output <- 'result/Table'

####################################################
## functions
####################################################
summarize_cont <- function(x) {
  x <- as.numeric(x)
  med <- median(x, na.rm = TRUE)
  q1  <- quantile(x, 0.25, na.rm = TRUE)
  q3  <- quantile(x, 0.75, na.rm = TRUE)
  sprintf("%.1f (%.1f–%.1f)", med, q1, q3)
}

summarize_cat <- function(x) {
  tbl <- table(x)
  pct <- prop.table(tbl) * 100
  paste0(names(tbl), ": ", tbl, " (", round(pct,1), "%)", collapse = "; ")
}

#################################################
## Load MultiAssayExperiment data
#################################################
load(file.path(dir_input, 'mae_mut_clin.RData'))
mut <- assay(mae[['mut_oncoprint']])
clin <- colData(mae[['mut_oncoprint']])

####################################################
## Clinbical metadata 
####################################################
# sample-level metadata
clin <- as.data.frame(clin)

# IDH status 
clin$IDH_status <- factor(clin$IDH_status, levels = c("WT", "Mut"))

# age
clin$Age <- ifelse(clin$Age >= 40, '>40', '<40')
clin$Age <- factor(clin$Age, levels = c(">40", "<40"))

# sex
clin$Sex <- factor(clin$Sex, levels = c("Male", "Female"))

# Grade
clin <- clin %>%
  mutate(
    Grade = case_when(
      WHO.2021.Grade == 1 ~ "I",
      WHO.2021.Grade == 2 ~ "II",
      WHO.2021.Grade == 3 ~ "III",
      WHO.2021.Grade == 4 ~ "IV",
      TRUE ~ NA_character_
    )
  )

clin$Grade <-  factor(clin$Grade, levels = c('I', 'II', 'III', 'IV'))

# Histology
clin$Histo <- clin$histo
clin <- clin %>%
  mutate(
    Histo = str_squish(Histo),
    Histo = str_to_lower(Histo),   # normalize first
    
    Histo = case_when(
      
      str_detect(Histo, "glioblastoma") ~ "Glioblastoma",
      str_detect(Histo, "oligodendroglioma") ~ "Oligodendroglioma",
      str_detect(Histo, "astrocytoma") ~ "Astrocytoma",
      str_detect(Histo, "ependymoma") ~ "Ependymoma",
      str_detect(Histo, "glioneuronal") ~ "Glioneuronal tumor",
      str_detect(Histo, "circumscribed glioma") ~ "Circumscribed glioma",
      str_detect(Histo, "pediatric-type high.?grade glioma") ~ "Pediatric-type HGG",
      str_detect(Histo, "pediatric-type low.?grade glioma") ~ "Pediatric-type LGG",
      TRUE ~ str_to_sentence(Histo)  # keep original but clean formatting
    )
  )

clin$Histo <- factor(clin$Histo, levels = c('Glioblastoma', 'Astrocytoma', 'Circumscribed glioma',
                                        'Glioneuronal tumor', 'Oligodendroglioma',  'Pediatric-type HGG',
                                        'Pediatric-type LGG', 'Ependymoma'))

# Location
clin <- clin %>%
  mutate(
    Primary.location = str_squish(Primary.location),
    Primary.location = str_to_title(Primary.location)
  )

clin <- clin %>%
  mutate(
    Location = case_when(
      Primary.location == "Lobar" ~ "Lobar",
      Primary.location == "Cerebellum" ~ "Cerebellum",
      Primary.location == "Thalamic" ~ "Thalamic",
      TRUE ~ "Other" # Spinal Cord, Intraventricular, Brainstem, Suprasellar
    )
  )

clin$Location <- factor(clin$Location, levels = c('Lobar', 'Cerebellum', 'Thalamic', 'Other'))

# ECOG
clin$ECOG <- clin$ECOG.Performance.Status
clin$ECOG <- substr(clin$ECOG, 1, 1)
clin$MGMT <- factor(clin$MGMT, levels = c("M", "U", "Unknown"))

# MGMT
clin$MGMT <- clin$MGMT.methylation
clin$MGMT <- case_when(
  str_detect(clin$MGMT, "Methylated") ~ "M",
  str_detect(clin$MGMT, "Unmethylated") ~ "U",
  TRUE ~ "Unknown"
)

clin$MGMT <- factor(clin$MGMT, levels = c("M", "U", "Unknown"))

# Recesion
clin$Resection <- clin$Extent.of.surgical.resection

#########################################################
## Summary statistics
#########################################################

vars_cat <- c('Age', 'Sex', 'Grade', 'Histo', 'Location', 'ECOG', 'MGMT', 'Resection', 'os.event')
vars_cont <- c("os.time")

clin_all <- clin
clin_wt  <- clin[clin$IDH_status == "WT", ]
clin_mut <- clin[clin$IDH_status == "Mut", ]


results <- data.frame(Variable = character(),
                      All = character(),
                      IDH_WT = character(),
                      IDH_Mut = character(),
                      stringsAsFactors = FALSE)

## Continuous
for (v in vars_cont) {
  results <- rbind(results, data.frame(
    Variable = v,
    All = summarize_cont(clin_all[[v]]),
    IDH_WT = summarize_cont(clin_wt[[v]]),
    IDH_Mut = summarize_cont(clin_mut[[v]])
  ))
}

## Categorical
for (v in vars_cat) {
  results <- rbind(results, data.frame(
    Variable = v,
    All = summarize_cat(clin_all[[v]]),
    IDH_WT = summarize_cat(clin_wt[[v]]),
    IDH_Mut = summarize_cat(clin_mut[[v]])
  ))
}

write.csv(results, file = file.path(dir_output, 'results.csv'))

