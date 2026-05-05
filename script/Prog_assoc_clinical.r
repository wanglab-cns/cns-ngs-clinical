##-------------------------------------------------------------------
## Script: Clinical Variable Association with Overall Survival (CNS Tumours)
##
## Purpose:
##   To evaluate the association between key clinical variables and
##   overall survival (OS) in CNS tumour patients using Cox
##   proportional hazards regression models.
##
##   Univariable survival analyses are performed across clinical
##   variables in the full cohort and within IDH-stratified subgroups.
##   Firth’s penalized Cox regression is applied to address small
##   sample size and sparse data issues.
##
##
## Input:
##   - result/data/mae_mut_clin.RData
##       A MultiAssayExperiment object containing:
##
##         • mut_binary assay:
##             Gene-level binary mutation matrix (used to derive IDH status)
##
##         • colData:
##             Harmonized clinical metadata including:
##               - os.time, os.event
##               - IDH_status
##               - Age, Sex
##               - WHO.2021.Grade
##               - Histology
##               - MGMT methylation
##               - ECOG performance status
##               - Extent of surgical resection
##
##
## Outputs:
##   CSV files summarizing univariable Cox regression results:
##
##     1) result/assoc/clinical/cox_os_clin_all.csv
##          All patients
##
##     2) result/assoc/clinical/cox_os_clin_wt.csv
##          IDH wild-type subset
##
##     3) result/assoc/clinical/cox_os_clin_mut.csv
##          IDH mutant subset
##
##     4) result/assoc/clinical/cox_os_clin_gbm.csv
##          Glioblastoma subset
##
##   Each output includes:
##     - variable name
##     - level (relative to reference)
##     - hazard ratio (HR)
##     - standard error (SE)
##     - sample size (N)
##     - 95% confidence interval (lower, upper)
##     - p-value
##
##
## Processing Overview:
##
##   1) Data Loading
##        Load MultiAssayExperiment object and extract clinical metadata.
##        IDH mutation status is re-derived from IDH1/IDH2 mutation data.
##
##   2) Clinical Variable Processing
##        Variables are harmonized and recoded:
##           - Age: <40 vs ≥40
##           - Grade: Low (I/II) vs High (III/IV)
##           - Histology: GBM vs Non-GBM
##           - MGMT: Methylated / Unmethylated / Unknown
##           - ECOG: categorical (0–3)
##           - Resection: Total / Subtotal / Biopsy
##
##        All variables are converted to factors with defined
##        reference levels for Cox model estimation.
##
##   3) Survival Preprocessing
##        Overall survival is administratively censored:
##           - All / IDH-WT / GBM: fixed time window (e.g., 36 months)
##           - IDH-Mut: extended follow-up window
##
##   4) Univariable Cox Modeling
##        For each clinical variable:
##           - Remove missing values
##           - Ensure ≥2 levels
##           - Fit penalized Cox model:
##                Surv(time, status) ~ variable
##
##   5) Stratified Analyses
##        Analyses are repeated for:
##           - Full cohort
##           - IDH wild-type
##           - IDH mutant
##           - GBM subset
##
##   6) Result Extraction
##        For each model:
##           - Estimate hazard ratios relative to reference levels
##           - Compute standard errors and confidence intervals
##
##   7) Export
##        Results are written to CSV files for downstream
##        analysis and reporting.
##
##
## Notes:
##   - Reference categories are defined by factor level ordering.
##   - Category collapsing is applied to reduce sparsity.
##   - Penalized Cox regression improves stability in small subgroups.
##   - Stratified analyses should be interpreted with caution due to
##     reduced sample size.
##-------------------------------------------------------------------
####################################################
## Load libraries
####################################################
library(MultiAssayExperiment)
library(ggplot2)
library(survival)
library(survminer)
library(coxphf)
library(dplyr) 
library(tidyr)
library(stringr)
library(paletteer)

####################################################
## Setup directories
####################################################
dir_input <- 'result/data'
dir_output <- 'result/assoc'

time.censor <- 36
time.censor.mut <- 120
n1.cutoff <- 3
n0.cutoff <- 3

#################################################
## Load data
#################################################
load(file.path(dir_input, 'mae_mut_clin.RData'))
mut <- assay(mae[['mut_binary']])
clin <- as.data.frame(colData(mae[['mut_binary']]))

## ---- Binary matrix
if(all(c("IDH1","IDH2") %in% rownames(mut))){

  idh_vec <- as.integer(
    mut["IDH1", ] == 1 | mut["IDH2", ] == 1
  )
  mut <- mut[!rownames(mut) %in% c("IDH1","IDH2"), , drop = FALSE]
  mut <- rbind(mut, IDH = idh_vec)

}

# Age
clin$Age <- ifelse(clin$Age >= 40, '>40', '<40')

# ECOG
clin$ECOG <- clin$ECOG.Performance.Status
clin$ECOG <- substr(clin$ECOG, 1, 1)

# MGMT
clin$MGMT <- clin$MGMT.methylation
clin$MGMT <- case_when(
  str_detect(clin$MGMT, "Methylated") ~ "M",
  str_detect(clin$MGMT, "Unmethylated") ~ "U",
  TRUE ~ "Unknown"
)

# Recesion
clin$Resection <- clin$Extent.of.surgical.resection

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

clin$Grade <- ifelse(clin$Grade %in% c('I', 'II'), 'Low', 'High')

# Histology
clin$Histo <- clin$histo
clin$Histo <- ifelse(
  clin$Histo == "Glioblastoma",
  "GBM",
  "Non-GBM"
)

mut <- mut[, colnames(mut) %in% clin$Study] # 199 patients
clin$IDH_status <- factor(
  clin$IDH_status,
  levels = c("WT", "Mut")
)

clin$Sex <- factor(
  clin$Sex,
  levels = c("Female", "Male")
)

clin$Age <- factor(
  clin$Age,
  levels = c(">40", "<40")
)

clin$Grade <-  factor(
  clin$Grade,
  levels = c("High", "Low")
)

clin$Histo <- factor(
  clin$Histo,
  levels = c("GBM", "Non-GBM")
)

clin$MGMT <- factor(
  clin$MGMT,
  levels = c("M", "U", "Unknown")
  )

clin$ECOG <- factor(
  clin$ECOG,
  levels = c("1", "0", "2", "3")
  )

clin$Resection <- factor(
  clin$Resection,
  levels = c("Subtotal", "Biopsy", "Total")
  ) 

####################################################
## OS association --> metadata clinical data
####################################################
## --- Step 1: all patients 
df <- clin[, c('Age', 'Sex', 'Grade', 'Histo', 'IDH_status', 'MGMT', 'ECOG', 'Resection')]

cox_res <- lapply(names(df), function(varname){

  data <- data.frame(
    status = clin$os.event,
    time   = as.numeric(clin$os.time),
    variable = df[[varname]]
  )

  # remove NA
  data <- data[!is.na(data$variable), ]

  # ensure factor (important)
  data$variable <- as.factor(data$variable)

  # censoring (vectorized)
  data$status[data$time > time.censor] <- 0
  data$time[data$time > time.censor] <- time.censor

  # skip if only one level
  if(nlevels(data$variable) < 2) return(NULL)

  # run Cox
  fit <- tryCatch(
    coxphf(Surv(time, status) ~ variable, data = data),
    error = function(e) return(NULL)
  )

  if(is.null(fit)) return(NULL)

  res <- data.frame(
    variable = varname,
    level = names(fit$coefficients),
    logHR = fit$coefficients,
    HR = exp(fit$coefficients),   
    se = sqrt(diag(fit$var)),
    n = nrow(data),
    low = fit$ci.lower,
    up  = fit$ci.upper,
    pval = fit$prob
  )

  res
})

cox_res <- do.call(rbind, cox_res)
cox_res <- cox_res[!is.na(cox_res$HR), ]

write.csv(cox_res,
          file = file.path(dir_output, 'clinical', 'cox_os_clin_all.csv'),
          row.names = FALSE)


## --- Step 2: WT patients 
clin_wt <- clin[clin$IDH_status == 'WT', ] # 152 cases
df <- clin_wt[, c('Age', 'Sex', 'Grade', 'Histo', 'MGMT', 'ECOG', 'Resection')]

cox_res <- lapply(names(df), function(varname){

  data <- data.frame(
    status = clin_wt$os.event,
    time   = as.numeric(clin_wt$os.time),
    variable = df[[varname]]
  )

  # remove NA
  data <- data[!is.na(data$variable), ]

  # ensure factor (important)
  data$variable <- as.factor(data$variable)

  # censoring (vectorized)
  data$status[data$time > time.censor] <- 0
  data$time[data$time > time.censor] <- time.censor

  # skip if only one level
  if(nlevels(data$variable) < 2) return(NULL)

  # run Cox
  fit <- tryCatch(
    coxphf(Surv(time, status) ~ variable, data = data),
    error = function(e) return(NULL)
  )

  if(is.null(fit)) return(NULL)

  res <- data.frame(
    variable = varname,
    level = names(fit$coefficients),
    logHR = fit$coefficients,
    HR = exp(fit$coefficients),   
    se = sqrt(diag(fit$var)),
    n = nrow(data),
    low = fit$ci.lower,
    up  = fit$ci.upper,
    pval = fit$prob
  )

  res
})

cox_res <- do.call(rbind, cox_res)
cox_res <- cox_res[!is.na(cox_res$HR), ]

write.csv(cox_res,
          file = file.path(dir_output, 'clinical', 'cox_os_clin_wt.csv'),
          row.names = FALSE)

## --- Step 3: Mut patients 
clin_mut <- clin[clin$IDH_status == 'Mut', ] # 47 cases
df <- clin_mut[, c('Age', 'Sex', 'Grade', 'MGMT', 'ECOG', 'Resection')]

cox_res <- lapply(names(df), function(varname){

  data <- data.frame(
    status = clin_mut$os.event,
    time   = as.numeric(clin_mut$os.time),
    variable = df[[varname]]
  )

  # remove NA
  data <- data[!is.na(data$variable), ]

  # ensure factor (important)
  data$variable <- as.factor(data$variable)

  # censoring (vectorized)
  data$status[data$time > time.censor.mut] <- 0
  data$time[data$time > time.censor.mut] <- time.censor.mut

  # skip if only one level
  if(nlevels(data$variable) < 2) return(NULL)

  # run Cox
  fit <- tryCatch(
    coxphf(Surv(time, status) ~ variable, data = data),
    error = function(e) return(NULL)
  )

  if(is.null(fit)) return(NULL)

  res <- data.frame(
    variable = varname,
    level = names(fit$coefficients),
    logHR = fit$coefficients,
    HR = exp(fit$coefficients),   
    se = sqrt(diag(fit$var)),
    n = nrow(data),
    low = fit$ci.lower,
    up  = fit$ci.upper,
    pval = fit$prob
  )

  res
})

cox_res <- do.call(rbind, cox_res)
cox_res <- cox_res[!is.na(cox_res$HR), ]

write.csv(cox_res,
          file = file.path(dir_output, 'clinical', 'cox_os_clin_mut.csv'),
          row.names = FALSE)
          

## --- Step 4: GBM patients 
clin_gbm <- clin[clin$Histo == 'GBM', ] # 105 cases
df <- clin_gbm[, c('Age', 'Sex', 'Grade', 'MGMT', 'ECOG', 'Resection')]

cox_res <- lapply(names(df), function(varname){

  data <- data.frame(
    status = clin_gbm$os.event,
    time   = as.numeric(clin_gbm$os.time),
    variable = df[[varname]]
  )

  # remove NA
  data <- data[!is.na(data$variable), ]

  # ensure factor (important)
  data$variable <- as.factor(data$variable)

  # censoring (vectorized)
  data$status[data$time > time.censor.mut] <- 0
  data$time[data$time > time.censor.mut] <- time.censor.mut

  # skip if only one level
  if(nlevels(data$variable) < 2) return(NULL)

  # run Cox
  fit <- tryCatch(
    coxphf(Surv(time, status) ~ variable, data = data),
    error = function(e) return(NULL)
  )

  if(is.null(fit)) return(NULL)

  res <- data.frame(
    variable = varname,
    level = names(fit$coefficients),
    logHR = fit$coefficients,
    HR = exp(fit$coefficients),   
    se = sqrt(diag(fit$var)),
    n = nrow(data),
    low = fit$ci.lower,
    up  = fit$ci.upper,
    pval = fit$prob
  )

  res
})

cox_res <- do.call(rbind, cox_res)
cox_res <- cox_res[!is.na(cox_res$HR), ]

write.csv(cox_res,
          file = file.path(dir_output, 'clinical', 'cox_os_clin_gbm.csv'),
          row.names = FALSE)
          
