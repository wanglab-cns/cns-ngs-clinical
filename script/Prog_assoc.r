##-------------------------------------------------------------------
## Script: CNS NGS Gene-Level Mutation – Overall Survival Association
##
## Purpose:
##   Assess associations between binary gene-level mutation
##   status (0 = wild-type, 1 = mutated) and overall survival (OS)
##   in CNS tumor patients using Cox proportional hazards models.
##
##   Perform univariable and multivariable survival analyses,
##   including stratified analyses by IDH status, to identify
##   mutation-associated prognostic signals.
##
##   Apply administrative censoring to ensure comparable
##   follow-up across patients.
##
##
## Input:
##   - result/data/mae_mut_clin.RData
##       MultiAssayExperiment containing:
##         - mut_binary assay:
##             gene × patient binary mutation matrix (0/1)
##         - colData:
##             harmonized clinical metadata including:
##               • os.time, os.event
##               • IDH_status
##               • Age, Sex
##               • WHO.2021.Grade
##               • ECOG
##               • MGMT methylation
##               • Extent of resection
##               • Histology annotations
##
##
## Outputs:
##
##   Univariable Cox models:
##     - result/assoc/cox_os_all.csv
##     - result/assoc/cox_os_wt.csv
##     - result/assoc/cox_os_mut.csv
##
##   Multivariable Cox models (clinical-adjusted):
##
##     Part 1:
##       - cox_os_clin_mv_part1.csv
##       - cox_os_clin_wt_mv_part1.csv
##       - cox_os_clin_mut_mv_part1.csv
##
##       Covariates:
##         Age + ECOG + Resection
##
##     Part 2:
##       - cox_os_clin_mv_part2.csv
##       - cox_os_clin_wt_mv_part2.csv
##       - cox_os_clin_mut_mv_part2.csv
##
##       Covariates:
##         Age + ECOG + Resection + MGMT
##
##
##   Each output includes:
##     gene, HR, SE, N,
##     95% CI (lower/upper),
##     p-value, and BH-adjusted FDR
##
##
## Processing Overview:
##
##   1) Data Loading
##        The MultiAssayExperiment object is loaded and the
##        binary mutation matrix and clinical metadata are
##        extracted.
##
##        IDH mutation status is re-derived from IDH1/IDH2
##        and appended as a gene-level feature.
##
##   2) Clinical Variable Processing
##        Clinical variables are recoded and harmonized:
##           - Age dichotomized (>40 vs <40)
##           - ECOG simplified
##           - MGMT methylation recoded (M / U / Unknown)
##           - Grade grouped (Low vs High)
##           - Histology grouped (GBM vs Non-GBM)
##
##   3) Survival Preprocessing
##        Overall survival is administratively censored:
##           - All patients: time.censor (e.g., 36 months)
##           - IDH-mutant patients: extended censoring window
##
##   4) Gene Filtering per Model
##        For each gene:
##           - Require minimum mutated (n1.cutoff) and
##             wild-type (n0.cutoff) sample counts
##           - Require sufficient event counts in both groups
##
##   5) Univariable Cox Analysis
##        For each gene:
##           - Fit Cox proportional hazards model:
##                Surv(time, status) ~ mutation
##
##        Analyses performed for:
##           - All patients
##           - IDH wild-type subset
##           - IDH mutant subset
##
##   6) Multivariable Cox Analysis
##        Clinical-adjusted Cox models are fitted:
##
##           Part 1:
##             Surv ~ mutation + Age + ECOG + Resection
##
##           Part 2:
##             Surv ~ mutation + Age + ECOG + Resection + MGMT
##
##        Models are applied separately to:
##           - All patients
##           - IDH wild-type
##           - IDH mutant
##
##   7) Multiple Testing Correction
##        P-values are adjusted across genes using the
##        Benjamini–Hochberg procedure to control FDR.
##
##   8) Export
##        All Cox model results are written to CSV files
##        for downstream interpretation and reporting.
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
                      HR = exp(fit$coefficients),   
                      se = sqrt(diag(fit$var)),
                      n = nrow(data),
                      low = fit$ci.lower,
                      up  = fit$ci.upper,
                      pval = fit$prob)
                      

  } else{
    
   res <- data.frame(gene = rownames(df)[k],
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
write.csv(cox_res, file = file.path(dir_output, 'mutation', 'cox_os_mut_all.csv'), row.names=FALSE)

## --- Step 2: WT patients 
clin_wt <- clin[clin$IDH_status == 'WT', ]
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
                      HR = exp(fit$coefficients),   
                      se = sqrt(diag(fit$var)),
                      n = nrow(data),
                      low = fit$ci.lower,
                      up  = fit$ci.upper,
                      pval = fit$prob)
                      

  } else{
    
   res <- data.frame(gene = rownames(df)[k],
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
write.csv(cox_res, file = file.path(dir_output, 'mutation', 'cox_os_mut_wt.csv'), row.names=FALSE)

## --- Step 3: Mut patients 
clin_mut <- clin[clin$IDH_status == 'Mut', ]
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
                      HR = exp(fit$coefficients),   
                      se = sqrt(diag(fit$var)),
                      n = nrow(data),
                      low = fit$ci.lower,
                      up  = fit$ci.upper,
                      pval = fit$prob)
                      

  } else{
    
   res <- data.frame(gene = rownames(df)[k],
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
write.csv(cox_res, file = file.path(dir_output, 'mutation', 'cox_os_mut_mut.csv'), row.names=FALSE)
########################################################################################################
########################################################################################################
######################################### MV OS association (Part 1) ###################################
########################################################################################################
########################################################################################################
## --- Step 1: all patients 
res <- read.csv(file.path(dir_output, 'clinical', 'cox_os_clin_all.csv'))
varnames <- res[res$pval < 0.05, 'variable']

df <- mut
cox_res <- lapply(1:nrow(df), function(k){

data <- data.frame( status=clin$os.event , 
                    time=clin$os.time , 
                    variable=as.numeric(unlist(df[k, ]) ) , 
                    #IDH =  clin$IDH_status, 
                    Age = clin$Age, 
                    #Grade = clin$Grade, 
                    #Histo = clin$Histo,
                    ECOG = clin$ECOG,
                    #MGMT = clin$MGMT,
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
write.csv(cox_res, file = file.path(dir_output, 'mutation & clinical', 'cox_os_all_mv_part1.csv'), row.names=FALSE)

## --- Step 2: WT patients
res <- read.csv(file.path(dir_output, 'clinical', 'cox_os_clin_wt.csv'))
varnames <- res[res$pval < 0.05, 'variable']

clin_wt <- clin[clin$IDH_status == 'WT', ]
df <- mut[, colnames(mut) %in% clin_wt$Study]

cox_res <- lapply(1:nrow(df), function(k){

data <- data.frame( status=clin_wt$os.event , 
                    time=clin_wt$os.time , 
                    variable=as.numeric(unlist(df[k, ]) ) , 
                    Age = clin_wt$Age, 
                    #Grade = clin_wt$Grade,
                    #Histo = clin_wt$Histo,
                    ECOG = clin_wt$ECOG,
                    #MGMT = clin_wt$MGMT,
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
    
    fit <- coxphf( formula= Surv( time , status ) ~ variable + Age  
                                                   + ECOG + Resection, data=data )
    res <- data.frame(gene = rownames(df)[k],
                      level = names(fit$coefficients),
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
write.csv(cox_res, file = file.path(dir_output, 'mutation & clinical', 'cox_os_wt_mv_part1.csv'), row.names=FALSE)

## --- Step 3: Mut patients
res <- read.csv(file.path(dir_output, 'clinical', 'cox_os_clin_mut.csv'))
varnames <- res[res$pval < 0.1, 'variable']

clin_mut <- clin[clin$IDH_status == 'Mut', ]
df <- mut[, colnames(mut) %in% clin_mut$Study]

cox_res <- lapply(1:nrow(df), function(k){

data <- data.frame( status=clin_mut$os.event , 
                    time=clin_mut$os.time , 
                    variable=as.numeric(unlist(df[k, ]) ), 
                    #Age = clin_mut$Age
                    #Grade = clin_mut$Grade,
                    #Histo = clin_mut$Histo,
                    #ECOG = clin_mut$ECOG,
                    #MGMT = clin_mut$MGMT,
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
    
    fit <- coxphf( formula= Surv( time , status ) ~ variable + Resection, data=data )
    res <- data.frame(gene = rownames(df)[k],
                      level = names(fit$coefficients),
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
write.csv(cox_res, file = file.path(dir_output, 'mutation & clinical', 'cox_os_mut_mv_part1.csv'), row.names=FALSE)

########################################################################################################
########################################################################################################
######################################### MV OS association (Part 2) ###################################
########################################################################################################
########################################################################################################
## --- Step 1: all patients 
res <- read.csv(file.path(dir_output, 'clinical', 'cox_os_clin_all.csv'))
varnames <- res[res$pval < 0.05, 'variable']

df <- mut
cox_res <- lapply(1:nrow(df), function(k){

data <- data.frame( status=clin$os.event , 
                    time=clin$os.time , 
                    variable=as.numeric(unlist(df[k, ]) ) , 
                    #IDH =  clin$IDH_status, 
                    Age = clin$Age, 
                    #Grade = clin$Grade, 
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
                                                   + ECOG + Resection + MGMT, data=data )
    res <- data.frame(gene = rownames(df)[k],
                      level = names(fit$coefficients),
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
write.csv(cox_res, file = file.path(dir_output, 'mutation & clinical', 'cox_os_clin_mv_part2.csv'), row.names=FALSE)

## --- Step 2: WT patients
res <- read.csv(file.path(dir_output, 'clinical', 'cox_os_clin_wt.csv'))
varnames <- res[res$pval < 0.05, 'variable']

clin_wt <- clin[clin$IDH_status == 'WT', ]
df <- mut[, colnames(mut) %in% clin_wt$Study]

cox_res <- lapply(1:nrow(df), function(k){

data <- data.frame( status=clin_wt$os.event , 
                    time=clin_wt$os.time , 
                    variable=as.numeric(unlist(df[k, ]) ) , 
                    Age = clin_wt$Age, 
                    #Grade = clin_wt$Grade,
                    #Histo = clin_wt$Histo,
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
    
    fit <- coxphf( formula= Surv( time , status ) ~ variable + Age   
                                                   + ECOG + Resection + MGMT, data=data )
    res <- data.frame(gene = rownames(df)[k],
                      level = names(fit$coefficients),
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
write.csv(cox_res, file = file.path(dir_output, 'mutation & clinical', 'cox_os_clin_wt_mv_part2.csv'), row.names=FALSE)
