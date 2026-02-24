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
n1.cutoff <- 5
n0.cutoff <- 5

#################################################
## Load data
#################################################
load(file.path(dir_input, 'se_mut_bin_clin.RData'))
mut <- assay(eset)
clin <- as.data.frame(colData(eset))
clin <- clin[clin$IDH.status != 'NA', ]
clin$Age <- ifelse(clin$Age >= 40, '>40', '<40')
clin$IDH.status[clin$IDH.status == 'IDHmut'] <- "Mut"
clin$IDH.status[clin$IDH.status == 'IDHwt'] <- "WT"
clin$IDH.status <- factor(clin$IDH.status,
                          levels = c("WT", "Mut"))
# primary location
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
      TRUE ~ "Other"
    )
  )

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

# Histology
clin$Histo <- coalesce(clin$LGG, clin$HGG)
clin$Histo <- str_squish(clin$Histo)   # remove extra spaces
clin$Histo <- str_to_sentence(clin$Histo)
clin <- clin %>%
  mutate(
    Histo = str_squish(Histo),
    Histo = str_to_sentence(Histo),
    
    Histo = case_when(
      str_detect(Histo, "Glioblastoma") ~ "Glioblastoma",
      
      str_detect(Histo, "Diffuse hemispheric") ~ "Diffuse hemispheric",
      str_detect(Histo, "Diffuse midline") ~ "Diffuse midline",
      str_detect(Histo, "Diffuse high") ~ "Diffuse high-grade",
      str_detect(Histo, "Diffuse low") ~ "Diffuse low-grade",
      
      str_detect(Histo, "Astrocytoma") ~ "Astrocytoma",
      str_detect(Histo, "Oligodendroglioma") ~ "Oligodendroglioma",
      str_detect(Histo, "Glioneuronal") ~ "Glioneuronal",
      str_detect(Histo, "Pilocytic") ~ "Pilocytic astrocytoma",
      str_detect(Histo, "Pleomorphic") ~ "Pleomorphic xanthoastrocytoma",
      str_detect(Histo, "Ependymoma") ~ "Ependymoma",
      
      TRUE ~ Histo
    )
  )

histo_counts <- table(clin$Histo)

# Identify small groups (<10 patients)
small_groups <- names(histo_counts[histo_counts < 10])

clin <- clin %>%
  mutate(
    Histo = ifelse(
      Histo %in% small_groups,
      "Other",
      Histo
    )
  )

mut <- mut[, colnames(mut) %in% clin$Study] # 191 patients
clin$IDH.status <- factor(
  clin$IDH.status,
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
  levels = c("I", "II", "III", "IV")
)

clin$Location <-  factor(
  clin$Location,
  levels = c("Lobar", "Cerebellum", "Thalamic", "Other")
)

clin$Histo <- factor(
  clin$Histo,
  levels = c("Glioblastoma", "Astrocytoma", "Glioneuronal", "Oligodendroglioma", "Other")
)

####################################################
## OS association --> no metadata adjustment
####################################################
## --- Step 1: all patients 
df <- mut
cox_res <- lapply(1:nrow(df), function(k){

data <- data.frame( status=clin$os.event , time=clin$os.time , variable=df[k, ] )
data <- data[!is.na(data$variable), ]
data$time <- as.numeric(as.character(data$time))
  
for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor ){
      data[ i , "time" ] = time.censor
      data[ i , "status" ] = 0
      
    }
  }
  
  if( length( data$variable[ data$variable == 1 ] )>= n1.cutoff &
      length( data$variable[ data$variable == 0 ] ) >= n0.cutoff ){
    
    cox <- coxph( formula= Surv( time , status ) ~ variable , data=data )
    res <- data.frame(gene = rownames(df)[k],
                      hr = summary(cox)$coefficients[, "coef"],
                      se = summary(cox)$coefficients[, "se(coef)"],
                      n = round(summary(cox)$n),
                      low = summary(cox)$conf.int[, "lower .95"],
                      up = summary(cox)$conf.int[, "upper .95"],
                      pval = summary(cox)$coefficients[, "Pr(>|z|)"])

  } else{
    
   res <- data.frame(gene = rownames(df)[k],
                     hr = NA,
                     se = NA,
                     n = NA,
                     low = NA,
                     up = NA,
                     pval = NA)
    
  }
  
  res

})

cox_res <- do.call(rbind, cox_res)
cox_res <- cox_res[!is.na(cox_res$hr), ]
cox_res$fdr <- p.adjust(cox_res$pval, method = 'BH')
write.csv(cox_res, file = file.path(dir_output, 'cox_os_all.csv'), row.names=FALSE)

## --- Step 2: WT patients 
clin_wt <- clin[clin$IDH.status == 'WT', ]
df <- mut[, colnames(mut) %in% clin_wt$Study]

cox_res <- lapply(1:nrow(df), function(k){

data <- data.frame( status=clin_wt$os.event , time=clin_wt$os.time , variable=df[k, ] )
data <- data[!is.na(data$variable), ]
data$time <- as.numeric(as.character(data$time))
  
for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor ){
      data[ i , "time" ] = time.censor
      data[ i , "status" ] = 0
      
    }
  }
  
  if( length( data$variable[ data$variable == 1 ] )>= n1.cutoff &
      length( data$variable[ data$variable == 0 ] ) >= n0.cutoff ){
    
    cox <- coxph( formula= Surv( time , status ) ~ variable , data=data )
    res <- data.frame(gene = rownames(df)[k],
                      hr = summary(cox)$coefficients[, "coef"],
                      se = summary(cox)$coefficients[, "se(coef)"],
                      n = round(summary(cox)$n),
                      low = summary(cox)$conf.int[, "lower .95"],
                      up = summary(cox)$conf.int[, "upper .95"],
                      pval = summary(cox)$coefficients[, "Pr(>|z|)"])

  } else{
    
   res <- data.frame(gene = rownames(df)[k],
                     hr = NA,
                     se = NA,
                     n = NA,
                     low = NA,
                     up = NA,
                     pval = NA)
    
  }
  
  res

})

cox_res <- do.call(rbind, cox_res)
cox_res <- cox_res[!is.na(cox_res$hr), ]
cox_res$fdr <- p.adjust(cox_res$pval, method = 'BH')
write.csv(cox_res, file = file.path(dir_output, 'cox_os_wt.csv'), row.names=FALSE)

## --- Step 3: Mut patients 
clin_mut <- clin[clin$IDH.status == 'Mut', ]
df <- mut[, colnames(mut) %in% clin_mut$Study]

cox_res <- lapply(1:nrow(df), function(k){

data <- data.frame( status=clin_mut$os.event , time=clin_mut$os.time , variable=df[k, ] )
data <- data[!is.na(data$variable), ]
data$time <- as.numeric(as.character(data$time))
  
for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor.mut ){
      data[ i , "time" ] = time.censor.mut
      data[ i , "status" ] = 0
      
    }
  }
  
  if( length( data$variable[ data$variable == 1 ] )>= n1.cutoff &
      length( data$variable[ data$variable == 0 ] ) >= n0.cutoff ){
    
    cox <- coxph( formula= Surv( time , status ) ~ variable , data=data )
    res <- data.frame(gene = rownames(df)[k],
                      hr = summary(cox)$coefficients[, "coef"],
                      se = summary(cox)$coefficients[, "se(coef)"],
                      n = round(summary(cox)$n),
                      low = summary(cox)$conf.int[, "lower .95"],
                      up = summary(cox)$conf.int[, "upper .95"],
                      pval = summary(cox)$coefficients[, "Pr(>|z|)"])

  } else{
    
   res <- data.frame(gene = rownames(df)[k],
                     hr = NA,
                     se = NA,
                     n = NA,
                     low = NA,
                     up = NA,
                     pval = NA)
    
  }
  
  res

})

cox_res <- do.call(rbind, cox_res)
cox_res <- cox_res[!is.na(cox_res$hr), ]
cox_res$fdr <- p.adjust(cox_res$pval, method = 'BH')
write.csv(cox_res, file = file.path(dir_output, 'cox_os_mut.csv'), row.names=FALSE)

####################################################
## MV OS association --> metadata adjustment
####################################################
## --- Step 1: all patients 
df <- mut
cox_res <- lapply(1:nrow(df), function(k){

data <- data.frame( status=clin$os.event , time=clin$os.time , variable=df[k, ], 
                    IDH =  clin$IDH.status, 
                    Age = clin$Age, 
                    Sex = clin$Sex,
                    Grade = clin$Grade, 
                    Histo = clin$Histo, 
                    Location = clin$Location )
data <- data[!is.na(data$variable), ]
data$time <- as.numeric(as.character(data$time))
  
for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor ){
      data[ i , "time" ] = time.censor
      data[ i , "status" ] = 0
      
    }
  }
  
  if( length( data$variable[ data$variable == 1 ] )>= n1.cutoff &
      length( data$variable[ data$variable == 0 ] ) >= n0.cutoff ){
    
    cox <- coxph( formula= Surv( time , status ) ~ variable + IDH + Age + Sex 
                                                   + Grade + Location + Histo, data=data )
    res <- data.frame(gene = rownames(df)[k],
                      hr = summary(cox)$coefficients['variable', "coef"],
                      se = summary(cox)$coefficients['variable', "se(coef)"],
                      n = round(summary(cox)$n),
                      low = summary(cox)$conf.int['variable', "lower .95"],
                      up = summary(cox)$conf.int['variable', "upper .95"],
                      pval = summary(cox)$coefficients['variable', "Pr(>|z|)"])

  } else{
    
   res <- data.frame(gene = rownames(df)[k],
                     hr = NA,
                     se = NA,
                     n = NA,
                     low = NA,
                     up = NA,
                     pval = NA)
    
  }
  
  res

})

cox_res <- do.call(rbind, cox_res)
cox_res <- cox_res[!is.na(cox_res$hr), ]
cox_res$fdr <- p.adjust(cox_res$pval, method = 'BH')
write.csv(cox_res, file = file.path(dir_output, 'cox_os_all_mv.csv'), row.names=FALSE)

## --- Step 2: WT patients
clin_wt <- clin[clin$IDH.status == 'WT', ]
df <- mut[, colnames(mut) %in% clin_wt$Study]

cox_res <- lapply(1:nrow(df), function(k){

data <- data.frame( status=clin_wt$os.event , time=clin_wt$os.time , variable=df[k, ], 
                    Age = clin_wt$Age, Sex = clin_wt$Sex, Grade = clin_wt$Grade, 
                    Histo = clin_wt$Histo, Location = clin_wt$Location)
data <- data[!is.na(data$variable), ]
data$time <- as.numeric(as.character(data$time))
  
for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor ){
      data[ i , "time" ] = time.censor
      data[ i , "status" ] = 0
      
    }
  }
  
  if( length( data$variable[ data$variable == 1 ] )>= n1.cutoff &
      length( data$variable[ data$variable == 0 ] ) >= n0.cutoff ){
    
    cox <- coxph( formula= Surv( time , status ) ~ variable + Age + Sex 
                                                   + Grade + Location + Histo, data=data )
    res <- data.frame(gene = rownames(df)[k],
                      hr = summary(cox)$coefficients['variable', "coef"],
                      se = summary(cox)$coefficients['variable', "se(coef)"],
                      n = round(summary(cox)$n),
                      low = summary(cox)$conf.int['variable', "lower .95"],
                      up = summary(cox)$conf.int['variable', "upper .95"],
                      pval = summary(cox)$coefficients['variable', "Pr(>|z|)"])

  } else{
    
   res <- data.frame(gene = rownames(df)[k],
                     hr = NA,
                     se = NA,
                     n = NA,
                     low = NA,
                     up = NA,
                     pval = NA)
    
  }
  
  res

})

cox_res <- do.call(rbind, cox_res)
cox_res <- cox_res[!is.na(cox_res$hr), ]
cox_res$fdr <- p.adjust(cox_res$pval, method = 'BH')
write.csv(cox_res, file = file.path(dir_output, 'cox_os_wt_mv.csv'), row.names=FALSE)

## --- Step 3: Mut patients
clin_mut <- clin[clin$IDH.status == 'Mut', ]
df <- mut[, colnames(mut) %in% clin_mut$Study]

cox_res <- lapply(1:nrow(df), function(k){

data <- data.frame( status=clin_mut$os.event , time=clin_mut$os.time , variable=df[k, ], 
                    Age = clin_mut$Age, Sex = clin_mut$Sex, Grade = clin_mut$Grade, 
                    Histo = clin_mut$Histo, Location = clin_mut$Location)
data <- data[!is.na(data$variable), ]
data$time <- as.numeric(as.character(data$time))
  
for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor.mut ){
      data[ i , "time" ] = time.censor.mut
      data[ i , "status" ] = 0
      
    }
  }
  
  if( length( data$variable[ data$variable == 1 ] )>= n1.cutoff &
      length( data$variable[ data$variable == 0 ] ) >= n0.cutoff ){
    
    cox <- coxph( formula= Surv( time , status ) ~ variable + Age + Sex 
                                                   + Grade + Location + Histo, data=data )
    res <- data.frame(gene = rownames(df)[k],
                      hr = summary(cox)$coefficients['variable', "coef"],
                      se = summary(cox)$coefficients['variable', "se(coef)"],
                      n = round(summary(cox)$n),
                      low = summary(cox)$conf.int['variable', "lower .95"],
                      up = summary(cox)$conf.int['variable', "upper .95"],
                      pval = summary(cox)$coefficients['variable', "Pr(>|z|)"])

  } else{
    
   res <- data.frame(gene = rownames(df)[k],
                     hr = NA,
                     se = NA,
                     n = NA,
                     low = NA,
                     up = NA,
                     pval = NA)
    
  }
  
  res

})

cox_res <- do.call(rbind, cox_res)
cox_res <- cox_res[!is.na(cox_res$hr), ]
cox_res$fdr <- p.adjust(cox_res$pval, method = 'BH')
write.csv(cox_res, file = file.path(dir_output, 'cox_os_mut_mv.csv'), row.names=FALSE)
