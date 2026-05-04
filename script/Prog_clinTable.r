##--------------------------------------------------------------------
## Description:
## This script processes curated clinical metadata from a CNS NGS
## cohort to generate baseline summary statistics stratified by IDH
## mutation status. The workflow integrates harmonized clinical
## annotations from a MultiAssayExperiment (MAE) object and produces
## a structured summary table for downstream reporting and
## publication-ready descriptive analyses.
##
## Workflow:
## 1. Load processed MultiAssayExperiment data
##      - Extract mutation oncoprint assay (for alignment)
##      - Extract associated clinical metadata (colData)
##
## 2. Curate and harmonize clinical metadata:
##      - Convert IDH status to categorical variable (WT vs Mut)
##      - Dichotomize age (<40 vs ≥40)
##      - Standardize sex, grade (WHO 2021), and histology labels
##      - Harmonize tumor location into grouped categories
##      - Clean ECOG performance status
##      - Encode MGMT methylation status (M / U / Unknown)
##      - Standardize surgical resection categories
##
## 3. Define analysis variables:
##      - Continuous variables:
##          * Overall survival time (os.time)
##      - Categorical variables:
##          * Age group, sex, grade, histology
##          * Tumor location, ECOG, MGMT, resection status
##          * Survival event (os.event)
##
## 4. Stratify cohort:
##      - All patients
##      - IDH wild-type (WT)
##      - IDH mutant (Mut)
##
## 5. Compute summary statistics:
##      - Continuous variables:
##          * Median and interquartile range (IQR)
##      - Categorical variables:
##          * Frequency counts and percentages (n, %)
##
## 6. Construct summary table:
##      - Rows correspond to clinical variables
##      - Columns represent:
##          * All patients
##          * IDH-WT subgroup
##          * IDH-Mut subgroup
##
## 7. Export results:
##      - Summary table saved as CSV file
##      - Suitable for Table 1 (baseline characteristics)
##
## Output:
## - Summary statistics table saved as:
##   result/Table/results.csv
##
## Notes:
## - Clinical data derived from curated MAE object
## - IDH status used as primary stratification variable
## - Age is dichotomized for categorical reporting
## - Continuous variables summarized using median (IQR)
## - Categorical variables summarized using n (%)
## - Output table represents descriptive statistics only
## - Table can be extended with statistical tests if required
##--------------------------------------------------------------------
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

