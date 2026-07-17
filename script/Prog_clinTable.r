##--------------------------------------------------------------------
## Script: CNS NGS Clinical Variables Table
## Description:
## This script generates baseline clinical summary statistics for
## a CNS NGS cohort stratified by IDH mutation status using
## harmonized clinical metadata from a MultiAssayExperiment object.
##
## Workflow:
##   1) Load processed clinical and molecular data
##   2) Harmonize and recode clinical variables including:
##        - Age
##        - Sex
##        - Grade
##        - Histology
##        - Location
##        - ECOG
##        - MGMT
##        - Resection
##   3) Stratify patients by:
##        - All patients
##        - IDH wild-type
##        - IDH mutant
##        - IDH wild-type glioblastoma subgroup
##        - IDH wild-type non-glioblastoma subgroup
##   4) Compute descriptive summary statistics:
##        - Median (IQR) for continuous variables
##        - Frequency and percentage for categorical variables
##   5) Export summary statistics table for downstream
##      reporting and visualization
##
## Output:
##   - result/Table/results.csv
##
## Notes:
##   - Clinical data are derived from a curated MAE object.
##   - IDH mutation status is used for cohort stratification.
##   - Output tables are intended for descriptive analyses
##     and manuscript reporting.
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
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(NA_character_)
  }

  med <- median(x)
  q1  <- quantile(x, 0.25)
  q3  <- quantile(x, 0.75)
  mn  <- min(x)
  mx  <- max(x)

  sprintf(
    "%.1f (%.1f–%.1f); range %.1f–%.1f",
    med, q1, q3, mn, mx
  )
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
clin_wt_gbm <- clin[clin$IDH_status == "WT" & clin$Histo == 'Glioblastoma', ]
clin_wt_nongbm <- clin[clin$IDH_status == "WT" & clin$Histo != 'Glioblastoma', ]

results <- data.frame(Variable = character(),
                      All = character(),
                      IDH_WT = character(),
                      IDH_Mut = character(),
                      IDH_WT_GBM = character(),
                      IDH_WT_nonGBM = character(),
                      stringsAsFactors = FALSE)

## Continuous
for (v in vars_cont) {
  results <- rbind(results, data.frame(
    Variable = v,
    All = summarize_cont(clin_all[[v]]),
    IDH_WT = summarize_cont(clin_wt[[v]]),
    IDH_Mut = summarize_cont(clin_mut[[v]]),
    IDH_WT_GBM = summarize_cont(clin_wt_gbm[[v]]),
    IDH_WT_nonGBM = summarize_cont(clin_wt_nongbm[[v]])
  ))
}

## Categorical
for (v in vars_cat) {
  results <- rbind(results, data.frame(
    Variable = v,
    All = summarize_cat(clin_all[[v]]),
    IDH_WT = summarize_cat(clin_wt[[v]]),
    IDH_Mut = summarize_cat(clin_mut[[v]]),
    IDH_WT_GBM = summarize_cat(clin_wt_gbm[[v]]),
    IDH_WT_nonGBM = summarize_cat(clin_wt_nongbm[[v]])
  ))
}

write.csv(results, file = file.path(dir_output, 'results.csv'))

