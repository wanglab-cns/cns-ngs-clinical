#####################################################
## Script: CNS NGS gene-level mutation – OS association
##
## Objective:
##   Assess associations between binary gene-level mutation
##   status (0 = wild-type, 1 = mutated) and overall survival (OS)
##   using Cox proportional hazards regression models.
##
##   Analyses performed:
##     1) Univariable Cox models (all patients)
##     2) Multivariable Cox models adjusted for clinical covariates
##     3) IDH-stratified analyses (IDH-mutant and IDH-wildtype separately)
##
##   Overall survival is administratively censored at a predefined
##   time horizon (time.censor, months) to ensure comparable follow-up.
##
## Input:
##   - result/data/se_mut_bin_clin.RData
##       SummarizedExperiment object containing:
##         • assay: binary mutation matrix (genes × patients; 0/1)
##         • colData: harmonized clinical metadata including:
##             - os.time   : overall survival time (months)
##             - os.event  : event indicator (1 = death, 0 = censored)
##             - IDH.status
##             - Age
##             - Sex
##             - WHO.2021.Grade
##             - Primary.location
##             - Histology annotations
##
## Output:
##   Per-gene Cox regression results (CSV files):
##
##     - cox_os_all.csv
##         Univariable models (all patients)
##
##     - cox_os_all_mv.csv
##         Multivariable models adjusted for:
##         IDH + Age + Sex + Grade + Location + Histology
##
##     - cox_os_wt.csv
##     - cox_os_wt_mv.csv
##         Univariable and multivariable models within IDH-WT patients
##
##     - cox_os_mut.csv
##     - cox_os_mut_mv.csv
##         Univariable and multivariable models within IDH-Mut patients
##
##   Each output includes:
##     gene, log(HR), SE, N, 95% CI (lower/upper),
##     Wald test p-value, and Benjamini–Hochberg FDR.
##
## Key Processing Steps:
##   1) Extract mutation matrix and harmonized clinical variables.
##   2) Recode and collapse clinical variables (age grouping,
##      grade, histology, location).
##   3) Apply administrative censoring at predefined time horizons.
##   4) For each gene:
##        - Require minimum mutated (n1.cutoff) and
##          wild-type (n0.cutoff) sample counts.
##        - Fit Cox proportional hazards model.
##   5) Adjust p-values across genes using BH FDR correction.
#####################################################
####################################################
## Load libraries
####################################################
library(MultiAssayExperiment)
library(ggplot2)
library(survival)
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
############################### OS association --> no metadata adjustment ##############################
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
write.csv(cox_res, file = file.path(dir_output, 'cox_os_all.csv'), row.names=FALSE)

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
write.csv(cox_res, file = file.path(dir_output, 'cox_os_wt.csv'), row.names=FALSE)

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
write.csv(cox_res, file = file.path(dir_output, 'cox_os_mut.csv'), row.names=FALSE)

##############################################################
## MV OS association --> metadata adjustment --> Part 1
##############################################################
## --- Step 1: all patients 
res <- read.csv(file.path(dir_output, 'cox_os_clin_all.csv'))
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
write.csv(cox_res, file = file.path(dir_output, 'cox_os_clin_mv_part1.csv'), row.names=FALSE)

## --- Step 2: WT patients
res <- read.csv(file.path(dir_output, 'cox_os_clin_wt.csv'))
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
write.csv(cox_res, file = file.path(dir_output, 'cox_os_clin_wt_mv_part1.csv'), row.names=FALSE)

## --- Step 3: Mut patients
res <- read.csv(file.path(dir_output, 'cox_os_clin_mut.csv'))
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
write.csv(cox_res, file = file.path(dir_output, 'cox_os_clin_mut_mv_part1.csv'), row.names=FALSE)

##############################################################
## MV OS association --> metadata adjustment --> Part 2
##############################################################
## --- Step 1: all patients 
res <- read.csv(file.path(dir_output, 'cox_os_clin_all.csv'))
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
write.csv(cox_res, file = file.path(dir_output, 'cox_os_clin_mv_part2.csv'), row.names=FALSE)

## --- Step 2: WT patients
res <- read.csv(file.path(dir_output, 'cox_os_clin_wt.csv'))
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
write.csv(cox_res, file = file.path(dir_output, 'cox_os_clin_wt_mv_part2.csv'), row.names=FALSE)

## --- Step 3: Mut patients
res <- read.csv(file.path(dir_output, 'cox_os_clin_mut.csv'))
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
write.csv(cox_res, file = file.path(dir_output, 'cox_os_clin_mut_mv_part2.csv'), row.names=FALSE)
