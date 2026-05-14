##-------------------------------------------------------------------
## Script: Clinical Variables and Overall Survival Association
##         (CNS Tumours)
##
## Purpose:
##   To evaluate associations between clinical variables and
##   overall survival (OS) in CNS tumour patients using
##   penalized Cox proportional hazards regression models.
##
##   Analyses are performed in the full cohort and within
##   IDH- and histology-stratified subgroups.
##
##
## Input:
##   - result/data/mae_mut_clin.RData
##
##       MultiAssayExperiment object containing:
##         • binary mutation matrix including combined
##           IDH1/IDH2 status
##         • harmonized clinical metadata
##
##
## Outputs:
##   CSV files containing univariable Cox regression
##   results for:
##
##     - Full cohort
##     - IDH wild-type subgroup
##     - IDH mutant subgroup
##     - GBM subgroup
##
##
## Processing Overview:
##
##   1) Load mutation and clinical data
##
##   2) Harmonize and recode clinical variables including:
##        - Age
##        - Sex
##        - Grade
##        - Histology
##        - MGMT
##        - ECOG (0–1 vs 2–3)
##        - Resection
##
##   3) Apply administrative censoring to overall survival
##
##   4) Perform penalized Cox regression analyses
##
##   5) Repeat analyses within clinically relevant
##      subgroups
##
##   6) Export results for downstream analysis and
##      reporting
##
##
## Notes:
##   - Reference levels are predefined for model
##     estimation.
##
##   - Category collapsing is applied to reduce
##     sparsity.
##
##   - Stratified analyses should be interpreted
##     cautiously due to reduced subgroup sample sizes.
##-------------------------------------------------------------------
####################################################
## Load libraries
####################################################
library(MultiAssayExperiment)
library(ggplot2)
library(survival)
library(coxphf)
library(survminer)
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

########################################################################################################
########################################################################################################
############################### OS association --> no metadata adjustment ##############################
########################################################################################################
########################################################################################################
## --- Step 1: all patients 
df <- mut
cox_res <- lapply(1:nrow(df), function(k){

data <- data.frame( status=clin$os.event , 
                    time=clin$os.time , 
                    variable=as.numeric(unlist(df[k, ]) ) )

data <- data[!is.na(data$variable), ]
data$variable <- factor(data$variable)

data$status[data$time > time.censor] <- 0
data$time[data$time > time.censor] <- time.censor

n1 <- sum(data$variable == 1)
n0 <- sum(data$variable == 0)
e1 <- sum(data$status[data$variable == 1] == 1)
e0 <- sum(data$status[data$variable == 0] == 1)

  if( n1 >= n1.cutoff & n0 >= n0.cutoff & e1 >= 2 & e0 >= 2 ){
    
    fit <- coxphf(Surv(time, status) ~ variable, data = data)
    res <- data.frame(gene = rownames(df)[k],
                      logHR = fit$coefficients,
                      HR = exp(fit$coefficients),     
                      se = sqrt(diag(fit$var)),
                      n = nrow(data),
                      low = fit$ci.lower,
                      up  = fit$ci.upper,
                      pval = fit$prob)
                      

  } else{
    
   res <- data.frame(gene = rownames(df)[k],
                     logHR = NA,
                     HR = NA,
                     se = NA,
                     n = NA,
                     low = NA,
                     up = NA,
                     pval = NA)
    
  }
  
  res

})

cox_res <- do.call(rbind, cox_res)
cox_res <- cox_res[!is.na(cox_res$HR), ]
cox_res$fdr <- p.adjust(cox_res$pval, method = 'BH')
cox_res$study <- 'All'
write.csv(cox_res, file = file.path(dir_output, 'mutation', 'cox_os_mut_all.csv'), row.names=FALSE)

## --- Step 2: WT patients 
clin_wt <- clin[clin$IDH_status == 'WT', ] # 152
df <- mut[, colnames(mut) %in% clin_wt$Study]

cox_res <- lapply(1:nrow(df), function(k){

data <- data.frame( status=clin_wt$os.event , 
                    time=clin_wt$os.time , 
                    variable=as.numeric(unlist(df[k, ]) ) )

data <- data[!is.na(data$variable), ]
data$variable <- factor(data$variable)

data$status[data$time > time.censor] <- 0
data$time[data$time > time.censor] <- time.censor

n1 <- sum(data$variable == 1)
n0 <- sum(data$variable == 0)
e1 <- sum(data$status[data$variable == 1] == 1)
e0 <- sum(data$status[data$variable == 0] == 1)

  if( n1 >= n1.cutoff & n0 >= n0.cutoff & e1 >= 2 & e0 >= 2 ){
    
    fit <- coxphf(Surv(time, status) ~ variable, data = data)
    res <- data.frame(gene = rownames(df)[k],
                      logHR = fit$coefficients, 
                      HR = exp(fit$coefficients),   
                      se = sqrt(diag(fit$var)),
                      n = nrow(data),
                      low = fit$ci.lower,
                      up  = fit$ci.upper,
                      pval = fit$prob)
                      

  } else{
    
   res <- data.frame(gene = rownames(df)[k],
                     logHR = NA,
                     HR = NA,
                     se = NA,
                     n = NA,
                     low = NA,
                     up = NA,
                     pval = NA)
    
  }
  
  res

})

cox_res <- do.call(rbind, cox_res)
cox_res <- cox_res[!is.na(cox_res$HR), ]
cox_res$fdr <- p.adjust(cox_res$pval, method = 'BH')
cox_res$study <- 'IDH wilde-type'
write.csv(cox_res, file = file.path(dir_output, 'mutation', 'cox_os_mut_wt.csv'), row.names=FALSE)

## --- Step 3: Mut patients 
clin_mut <- clin[clin$IDH_status == 'Mut', ] # 47
df <- mut[, colnames(mut) %in% clin_mut$Study]

cox_res <- lapply(1:nrow(df), function(k){

data <- data.frame( status=clin_mut$os.event , 
                    time=clin_mut$os.time , 
                    variable=as.numeric(unlist(df[k, ]) ) )

data <- data[!is.na(data$variable), ]
data$variable <- factor(data$variable)

data$status[data$time > time.censor.mut] <- 0
data$time[data$time > time.censor.mut] <- time.censor.mut

n1 <- sum(data$variable == 1)
n0 <- sum(data$variable == 0)
e1 <- sum(data$status[data$variable == 1] == 1)
e0 <- sum(data$status[data$variable == 0] == 1)

  if( n1 >= n1.cutoff & n0 >= n0.cutoff & e1 >= 2 & e0 >= 2 ){
    
    fit <- coxphf(Surv(time, status) ~ variable, data = data)
    res <- data.frame(gene = rownames(df)[k],
                      logHR = fit$coefficients, 
                      HR = exp(fit$coefficients),   
                      se = sqrt(diag(fit$var)),
                      n = nrow(data),
                      low = fit$ci.lower,
                      up  = fit$ci.upper,
                      pval = fit$prob)
                      

  } else{
    
   res <- data.frame(gene = rownames(df)[k],
                     logHR = NA, 
                     HR = NA,
                     se = NA,
                     n = NA,
                     low = NA,
                     up = NA,
                     pval = NA)
    
  }
  
  res

})

cox_res <- do.call(rbind, cox_res)
cox_res <- cox_res[!is.na(cox_res$HR), ]
cox_res$fdr <- p.adjust(cox_res$pval, method = 'BH')
cox_res$study <- 'IDH mutant'
write.csv(cox_res, file = file.path(dir_output, 'mutation', 'cox_os_mut_mut.csv'), row.names=FALSE)

## --- Step 4: GBM patients 
clin_gbm <- clin[clin$Histo == 'GBM', ] # 105
df <- mut[, colnames(mut) %in% clin_gbm$Study]

cox_res <- lapply(1:nrow(df), function(k){

data <- data.frame( status=clin_gbm$os.event , 
                    time=clin_gbm$os.time , 
                    variable=as.numeric(unlist(df[k, ]) ) )

data <- data[!is.na(data$variable), ]
data$variable <- factor(data$variable)

data$status[data$time > time.censor] <- 0
data$time[data$time > time.censor] <- time.censor

n1 <- sum(data$variable == 1)
n0 <- sum(data$variable == 0)
e1 <- sum(data$status[data$variable == 1] == 1)
e0 <- sum(data$status[data$variable == 0] == 1)

  if( n1 >= n1.cutoff & n0 >= n0.cutoff & e1 >= 2 & e0 >= 2 ){
    
    fit <- coxphf(Surv(time, status) ~ variable, data = data)
    res <- data.frame(gene = rownames(df)[k],
                      logHR = fit$coefficients, 
                      HR = exp(fit$coefficients),   
                      se = sqrt(diag(fit$var)),
                      n = nrow(data),
                      low = fit$ci.lower,
                      up  = fit$ci.upper,
                      pval = fit$prob)
                      

  } else{
    
   res <- data.frame(gene = rownames(df)[k],
                     logHR = NA, 
                     HR = NA,
                     se = NA,
                     n = NA,
                     low = NA,
                     up = NA,
                     pval = NA)
    
  }
  
  res

})

cox_res <- do.call(rbind, cox_res)
cox_res <- cox_res[!is.na(cox_res$HR), ]
cox_res$fdr <- p.adjust(cox_res$pval, method = 'BH')
cox_res$study <- 'GBM'
write.csv(cox_res, file = file.path(dir_output, 'mutation', 'cox_os_mut_gbm.csv'), row.names=FALSE)

########################################################################################################
########################################################################################################
######################################### MV OS association (Part 1) ###################################
########################################################################################################
########################################################################################################
## --- Step 1: all patients 
res <- read.csv(file.path(dir_output, 'clinical', 'cox_os_clin_all.csv'))
varnames <- res[res$pval < 0.05, 'variable']

df <- mut
mut_freq <- rowMeans(df == 1)
genes <- names(mut_freq[mut_freq >= 0.1])
df <- df[rownames(df) %in% genes, ]

cox_res <- lapply(1:nrow(df), function(k){

data <- data.frame( status=clin$os.event , 
                    time=clin$os.time , 
                    variable=as.numeric(unlist(df[k, ]) ) , 
                    IDH =  clin$IDH_status, 
                    Age = clin$Age, 
                    Grade = clin$Grade, 
                    Histo = clin$Histo,
                    ECOG = clin$ECOG,
                    MGMT = clin$MGMT,
                    Resection = clin$Resection
                    )
  
data <- data[!is.na(data$variable), ]
data$variable <- factor(data$variable)

data$status[data$time > time.censor] <- 0
data$time[data$time > time.censor] <- time.censor
  
  n1 <- sum(data$variable == 1)
  n0 <- sum(data$variable == 0)
  e1 <- sum(data$status[data$variable == 1] == 1)
  e0 <- sum(data$status[data$variable == 0] == 1)

  if( n1 >= n1.cutoff & n0 >= n0.cutoff & e1 >= 2 & e0 >= 2 ){
    
    fit <- coxphf( formula= Surv( time , status ) ~ variable + Age    
                                                   + ECOG + Resection, data=data )
    res <- data.frame(gene = rownames(df)[k],
                      level = names(fit$coefficients),
                      logHR = fit$coefficients,
                      HR = exp(fit$coefficients),
                      se = sqrt(diag(fit$var)),
                      n = nrow(data),
                      low = fit$ci.lower,
                      up  = fit$ci.upper,
                      pval = fit$prob)
    
    res <- res[grep("^variable", res$level), ]

  } else{
    
   res <- data.frame(gene = rownames(df)[k],
                     level = NA,
                     logHR = NA,
                     HR = NA,
                     se = NA,
                     n = NA,
                     low = NA,
                     up = NA,
                     pval = NA)
    
  }
  
  res

})

cox_res <- do.call(rbind, cox_res)
cox_res <- cox_res[!is.na(cox_res$HR), ]
cox_res$fdr <- p.adjust(cox_res$pval, method = 'BH')
cox_res$study <- 'All'
write.csv(cox_res, file = file.path(dir_output, 'mutation & clinical', 'cox_os_all_mv_part1.csv'), row.names=FALSE)

## --- Step 2: WT patients
res <- read.csv(file.path(dir_output, 'clinical', 'cox_os_clin_wt.csv'))
varnames <- res[res$pval < 0.05, 'variable']

clin_wt <- clin[clin$IDH_status == 'WT', ]
df <- mut[, colnames(mut) %in% clin_wt$Study]
mut_freq <- rowMeans(df == 1)
genes <- names(mut_freq[mut_freq >= 0.1])
df <- df[rownames(df) %in% genes, ]

cox_res <- lapply(1:nrow(df), function(k){

data <- data.frame( status=clin_wt$os.event , 
                    time=clin_wt$os.time , 
                    variable=as.numeric(unlist(df[k, ]) ) , 
                    Age = clin_wt$Age, 
                    Grade = clin_wt$Grade,
                    Histo = clin_wt$Histo,
                    ECOG = clin_wt$ECOG,
                    MGMT = clin_wt$MGMT,
                    Resection = clin_wt$Resection
                    )

data <- data[!is.na(data$variable), ]
data$variable <- factor(data$variable)

data$status[data$time > time.censor] <- 0
data$time[data$time > time.censor] <- time.censor
  
  n1 <- sum(data$variable == 1)
  n0 <- sum(data$variable == 0)
  e1 <- sum(data$status[data$variable == 1] == 1)
  e0 <- sum(data$status[data$variable == 0] == 1)

  if( n1 >= n1.cutoff & n0 >= n0.cutoff & e1 >= 2 & e0 >= 2 ){
    
    fit <- coxphf( formula= Surv( time , status ) ~ variable + Age +   
                                                   + ECOG + Resection, data=data )
    res <- data.frame(gene = rownames(df)[k],
                      level = names(fit$coefficients),
                      logHR = fit$coefficients,
                      HR = exp(fit$coefficients),
                      se = sqrt(diag(fit$var)),
                      n = nrow(data),
                      low = fit$ci.lower,
                      up  = fit$ci.upper,
                      pval = fit$prob)
    
    res <- res[grep("^variable", res$level), ]

  } else{
    
   res <- data.frame(gene = rownames(df)[k],
                     level = NA,
                     logHR = NA,
                     HR = NA,
                     se = NA,
                     n = NA,
                     low = NA,
                     up = NA,
                     pval = NA)
    
  }
  
  res

})

cox_res <- do.call(rbind, cox_res)
cox_res <- cox_res[!is.na(cox_res$HR), ]
cox_res$fdr <- p.adjust(cox_res$pval, method = 'BH')
cox_res$study <- 'IDH wild-type'
write.csv(cox_res, file = file.path(dir_output, 'mutation & clinical', 'cox_os_wt_mv_part1.csv'), row.names=FALSE)

## --- Step 3: Mut patients
res <- read.csv(file.path(dir_output, 'clinical', 'cox_os_clin_mut.csv'))
varnames <- res[res$pval < 0.1, 'variable']

clin_mut <- clin[clin$IDH_status == 'Mut', ]
df <- mut[, colnames(mut) %in% clin_mut$Study]
mut_freq <- rowMeans(df == 1)
genes <- names(mut_freq[mut_freq >= 0.1])
df <- df[rownames(df) %in% genes, ]

cox_res <- lapply(1:nrow(df), function(k){

data <- data.frame( status=clin_mut$os.event , 
                    time=clin_mut$os.time , 
                    variable=as.numeric(unlist(df[k, ]) ), 
                    Age = clin_mut$Age,
                    Grade = clin_mut$Grade,
                    Histo = clin_mut$Histo,
                    ECOG = clin_mut$ECOG,
                    MGMT = clin_mut$MGMT,
                    Resection = clin_mut$Resection
                    )

data <- data[!is.na(data$variable), ]
data$variable <- factor(data$variable)

data$status[data$time > time.censor.mut] <- 0
data$time[data$time > time.censor.mut] <- time.censor.mut
  
  n1 <- sum(data$variable == 1)
  n0 <- sum(data$variable == 0)
  e1 <- sum(data$status[data$variable == 1] == 1)
  e0 <- sum(data$status[data$variable == 0] == 1)

  if( n1 >= n1.cutoff & n0 >= n0.cutoff & e1 >= 2 & e0 >= 2 ){
    
    fit <- coxphf( formula= Surv( time , status ) ~ variable +  Resection, data=data )
    res <- data.frame(gene = rownames(df)[k],
                      level = names(fit$coefficients),
                      logHR = fit$coefficients,
                      HR = exp(fit$coefficients),
                      se = sqrt(diag(fit$var)),
                      n = nrow(data),
                      low = fit$ci.lower,
                      up  = fit$ci.upper,
                      pval = fit$prob)
    
    res <- res[grep("^variable", res$level), ]

  } else{
    
   res <- data.frame(gene = rownames(df)[k],
                     level = NA,
                     logHR = NA,
                     HR = NA,
                     se = NA,
                     n = NA,
                     low = NA,
                     up = NA,
                     pval = NA)
    
  }
  
  res

})

cox_res <- do.call(rbind, cox_res)
cox_res <- cox_res[!is.na(cox_res$HR), ]
cox_res$fdr <- p.adjust(cox_res$pval, method = 'BH')
cox_res$study <- 'IDH mutant'
write.csv(cox_res, file = file.path(dir_output, 'mutation & clinical', 'cox_os_mut_mv_part1.csv'), row.names=FALSE)

## --- Step 3: Mut patients
res <- read.csv(file.path(dir_output, 'clinical', 'cox_os_clin_gbm.csv'))
varnames <- res[res$pval < 0.05, 'variable']

clin_gbm <- clin[clin$Histo == 'GBM', ]
df <- mut[, colnames(mut) %in% clin_gbm$Study]
mut_freq <- rowMeans(df == 1)
genes <- names(mut_freq[mut_freq >= 0.1])
df <- df[rownames(df) %in% genes, ]

cox_res <- lapply(1:nrow(df), function(k){

data <- data.frame( status=clin_gbm$os.event , 
                    time=clin_gbm$os.time , 
                    variable=as.numeric(unlist(df[k, ]) ), 
                    Age = clin_gbm$Age,
                    Grade = clin_gbm$Grade,
                    Histo = clin_gbm$Histo,
                    ECOG = clin_gbm$ECOG,
                    MGMT = clin_gbm$MGMT,
                    Resection = clin_gbm$Resection
                    )

data <- data[!is.na(data$variable), ]
data$variable <- factor(data$variable)

data$status[data$time > time.censor] <- 0
data$time[data$time > time.censor] <- time.censor
  
  n1 <- sum(data$variable == 1)
  n0 <- sum(data$variable == 0)
  e1 <- sum(data$status[data$variable == 1] == 1)
  e0 <- sum(data$status[data$variable == 0] == 1)

  if( n1 >= n1.cutoff & n0 >= n0.cutoff & e1 >= 2 & e0 >= 2 ){
    
    fit <- coxphf( formula= Surv( time , status ) ~ variable + Resection, data=data )
    res <- data.frame(gene = rownames(df)[k],
                      level = names(fit$coefficients),
                      logHR = fit$coefficients,
                      HR = exp(fit$coefficients),
                      se = sqrt(diag(fit$var)),
                      n = nrow(data),
                      low = fit$ci.lower,
                      up  = fit$ci.upper,
                      pval = fit$prob)
    
    res <- res[grep("^variable", res$level), ]

  } else{
    
   res <- data.frame(gene = rownames(df)[k],
                     level = NA,
                     logHR = NA,
                     HR = NA,
                     se = NA,
                     n = NA,
                     low = NA,
                     up = NA,
                     pval = NA)
    
  }
  
  res

})

cox_res <- do.call(rbind, cox_res)
cox_res <- cox_res[!is.na(cox_res$HR), ]
cox_res$fdr <- p.adjust(cox_res$pval, method = 'BH')
cox_res$study <- 'GBM'
write.csv(cox_res, file = file.path(dir_output, 'mutation & clinical', 'cox_os_gbm_mv_part1.csv'), row.names=FALSE)

