##----------------------------------------------------------------------------
## Objective:
##   Assess associations between clinical variables and
##   overall survival (OS) using Cox proportional hazards
##   regression models.
##
##   Analyses performed:
##     1) Univariable Cox models (all patients)
##     2) IDH-stratified analyses (IDH-mutant and IDH-wildtype separately)
##
##   Overall survival is administratively censored at predefined
##   time horizons (time.censor, time.censor.mut; months) to ensure
##   comparable follow-up across analyses.
##
## Input:
##   - result/data/mae_mut_clin.RData
##       MultiAssayExperiment object containing:
##         • mut_binary: gene-level mutation matrix (not directly used here)
##         • colData: harmonized clinical metadata including:
##             - os.time   : overall survival time (months)
##             - os.event  : event indicator (1 = death, 0 = censored)
##             - IDH_status
##             - Age
##             - Sex
##             - WHO.2021.Grade
##             - Histology annotations
##             - MGMT methylation
##             - ECOG performance status
##             - Extent of surgical resection
##
## Output:
##   Per-variable Cox regression results (CSV files):
##
##     - cox_os_clin_all.csv
##         Univariable models (all patients)
##
##     - cox_os_clin_wt.csv
##         Univariable models within IDH-WT patients
##
##     - cox_os_clin_mut.csv
##         Univariable models within IDH-Mut patients
##
##   Each output includes:
##     variable, level, HR, SE, N,
##     95% CI (lower/upper), and Wald test p-value.
##
## Key Processing Steps:
##   1) Extract clinical variables from MultiAssayExperiment colData.
##   2) Recode and collapse variables to reduce sparsity:
##        - Age: <40 vs ≥40
##        - Grade: Low (I/II) vs High (III/IV)
##        - Histology: GBM vs Non-GBM
##        - MGMT: Methylated / Unmethylated / Unknown
##        - ECOG: simplified categories (0–3)
##        - Resection: Total / Subtotal / Biopsy
##   3) Remove samples with missing values per variable.
##   4) Apply administrative censoring at predefined time horizons.
##   5) For each variable:
##        - Ensure ≥2 levels after filtering
##        - Fit Cox proportional hazards model:
##              Surv(time, status) ~ variable
##   6) Extract hazard ratios relative to reference levels.
##
## Notes:
##   - Factor levels define reference categories for HR estimation.
##   - Collapsing categories is performed to avoid sparse groups
##     and unstable Cox model estimates (e.g., zero-event strata).
##   - Stratified analyses (IDH-WT and IDH-Mut) may have reduced
##     sample sizes and should be interpreted with caution.
##----------------------------------------------------------------------------
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
          file = file.path(dir_output, 'cox_os_clin_all.csv'),
          row.names = FALSE)


## --- Step 2: WT patients 
clin_wt <- clin[clin$IDH_status == 'WT', ]
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
          file = file.path(dir_output, 'cox_os_clin_wt.csv'),
          row.names = FALSE)

## --- Step 3: Mut patients 
clin_mut <- clin[clin$IDH_status == 'Mut', ]
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
          file = file.path(dir_output, 'cox_os_clin_mut.csv'),
          row.names = FALSE)
          
