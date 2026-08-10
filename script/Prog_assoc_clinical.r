##-------------------------------------------------------------------
## Script: Clinical Variable Associations with Overall Survival in
##         CNS Tumours
## Purpose:
##   To evaluate associations between key clinical variables and
##   overall survival (OS) in CNS tumour patients using Cox
##   proportional hazards regression models.
##
##   Univariable survival analyses are performed in the full cohort
##   and within IDH-stratified and GBM subgroups. Firth’s penalized
##   Cox regression is applied to improve model stability in small
##   or sparse subgroups.
##
## Input:
##   - result/data/mae_mut_clin.RData
##       MultiAssayExperiment object containing:
##         - binary mutation data
##         - harmonized clinical metadata
##
## Outputs:
##   CSV files summarizing univariable Cox regression results for:
##     - Full cohort
##     - IDH mutant subgroup
##     - IDH wild-type glioblastoma subgroup
##     - IDH wild-type non-glioblastoma subgroup
##
##   Output tables include:
##     - clinical variable
##     - hazard ratio (HR)
##     - standard error (SE)
##     - confidence intervals
##     - p-value
##
## Processing Overview:
##   1) Load mutation and clinical data
##   2) Harmonize and recode clinical variables including:
##        - Age
##        - Sex
##        - Grade
##        - Histology
##        - MGMT
##        - ECOG (0–1 vs 2–3)
##        - Resection
##   3) Apply administrative censoring to overall survival
##   4) Perform univariable penalized Cox regression analyses
##   5) Repeat analyses within clinically relevant subgroups
##   6) Export results for downstream analysis and reporting
##
## Notes:
##   - Factor reference levels are predefined for model estimation.
##   - Category collapsing is applied where appropriate to reduce sparsity.
##   - Stratified analyses should be interpreted cautiously due to reduced subgroup sample sizes.
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
clin$ECOG <- ifelse(clin$ECOG %in% c(0, 1), "ECOG 0–1", "ECOG 2–3") 

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
  levels = c("ECOG 0–1", "ECOG 2–3")
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
cox_res$study <- 'All'

write.csv(cox_res,
          file = file.path(dir_output, 'clinical', 'cox_os_clin_all.csv'),
          row.names = FALSE)

## --- Step 2: Mut patients 
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
cox_res$study <- 'Mut'

write.csv(cox_res,
          file = file.path(dir_output, 'clinical', 'cox_os_clin_mut.csv'),
          row.names = FALSE)          

## --- Step 3: WT GBM patients 
clin_wt_gbm <- clin[clin$Histo == 'GBM', ] # 105 cases
df <- clin_wt_gbm[, c('Age', 'Sex', 'Grade', 'MGMT', 'ECOG', 'Resection')]

cox_res <- lapply(names(df), function(varname){

  data <- data.frame(
    status = clin_wt_gbm$os.event,
    time   = as.numeric(clin_wt_gbm$os.time),
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
cox_res$study <- 'WT-GBM'

write.csv(cox_res,
          file = file.path(dir_output, 'clinical', 'cox_os_clin_wt_gbm.csv'),
          row.names = FALSE)          

## --- Step 4: nonGBM patients 
clin_wt_nongbm <- clin[clin$Histo != 'GBM' & clin$IDH_status == 'WT', ] # 43 cases
df <- clin_wt_nongbm[, c('Age', 'Sex', 'Grade', 'MGMT', 'ECOG', 'Resection')]

cox_res <- lapply(names(df), function(varname){

  data <- data.frame(
    status = clin_wt_nongbm$os.event,
    time   = as.numeric(clin_wt_nongbm$os.time),
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
cox_res$study <- 'WT-nonGBM'

write.csv(cox_res,
          file = file.path(dir_output, 'clinical', 'cox_os_clin_wt_nongbm.csv'),
          row.names = FALSE)
          
