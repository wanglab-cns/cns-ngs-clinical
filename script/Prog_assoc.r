#####################################################
## Script: CNS NGS mutation–OS association (binary 0/1)
##
## Purpose:
##   Test gene-level mutation status (0/1 per patient) for association with
##   overall survival (OS) using univariable Cox proportional hazards models.
##   OS is optionally administratively censored at a fixed time horizon.
##
## Input:
##   - result/data/se_mut_bin_clin.RData
##       SummarizedExperiment with:
##         assay: mat_bin (genes x patients; 0/1 indicating any alteration)
##         colData: clinical metadata including:
##           - os.time  (OS in months)
##           - os.event (1=death/event, 0=censored)
##
## Output:
##   - In-memory data.frame `cox_res` with per-gene Cox results:
##       gene, hr, se, n, 95% CI (low/up), pval, fdr
##
## Key steps:
##   1) Load SummarizedExperiment; extract binary mutation matrix and OS fields.
##   2) Apply administrative censoring at `time.censor` months (status -> 0 if
##      censored).
##   3) For each gene, fit Cox PH model: Surv(os.time, os.event) ~ mutation(0/1),
##      requiring minimum counts in both mutated and non-mutated groups.
##   4) Adjust p-values across genes using Benjamini–Hochberg FDR.
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
library(paletteer)

####################################################
## Setup directories
####################################################
dir_input <- 'result/data'
dir_output <- 'result/Fig'

#################################################
## Load data
#################################################
load(file.path(dir_input, 'se_mut_bin_clin.RData'))
mut <- assay(eset)
clin <- colData(eset)
clin <- clin[clin$IDH.status != 'NA', ]
clin$Age <- ifelse(clin$Age >= 40, '> 40', '< 40')
clin$IDH.status[clin$IDH.status == 'IDHmut'] <- "Mut"
clin$IDH.status[clin$IDH.status == 'IDHwt'] <- "WT"

mut <- mut[, colnames(mut) %in% clin$Study] # 191 patients

time.censor <- 36
n1.cutoff <- 5
n0.cutoff <- 5

####################################################
## KM figure
####################################################
surv_obj <- with(clin, Surv(os.time, os.event))

## IDH status
fit_idh <- survfit(surv_obj ~ IDH.status, data = clin)

p_idh <- ggsurvplot(
  fit_idh,
  data = clin,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "IDH status",
  legend.labs = c("WT", "Mut"),
  palette = c("#C44E52", "#55A868"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
)

pdf(file.path(dir_output, "KM_IDH_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

####################################################
## OS association --> no metadata adjustment
####################################################
## --- Step 1: all patients (Tier I and II ---> remove Tier III)
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
    
   res <- data.frame(gene = rownames(mut)[k],
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

####################################################
## OS association --> metadata adjustment
####################################################
## --- Step 2: all patients (Tier I and II ---> remove Tier III)
df <- mut
clin$IDH.status <- clin$IDH.status <- factor(
  clin$IDH.status,
  levels = c("WT", "Mut")
)

cox_res <- lapply(2:nrow(df), function(k){

data <- data.frame( status=clin$os.event , time=clin$os.time , variable=df[k, ], IDH =  clin$IDH.status)
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
    
    cox <- coxph( formula= Surv( time , status ) ~ variable + IDH , data=data )
    res <- data.frame(gene = rownames(df)[k],
                      hr = summary(cox)$coefficients['variable', "coef"],
                      se = summary(cox)$coefficients['variable', "se(coef)"],
                      n = round(summary(cox)$n),
                      low = summary(cox)$conf.int['variable', "lower .95"],
                      up = summary(cox)$conf.int['variable', "upper .95"],
                      pval = summary(cox)$coefficients['variable', "Pr(>|z|)"])

  } else{
    
   res <- data.frame(gene = rownames(mut)[k],
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
write.csv(cox_res, file = file.path(dir_output, 'cox_os_IDH.csv'), row.names=FALSE)

####################################################
## OS association --> stratify for IHD WT
####################################################
## --- Step 1: all patients (Tier I and II ---> remove Tier III)
df <- mut
clin$IDH.status <- clin$IDH.status <- factor(
  clin$IDH.status,
  levels = c("WT", "Mut")
)

clin_wt <- clin[clin$IDH.status == 'WT', ]
df <- mut[, colnames(mut) %in% clin_wt$Study] # 142 patients

cox_res <- lapply(2:nrow(df), function(k){

data <- data.frame( status=clin_wt$os.event , time=clin_wt$os.time , variable=df[k, ])
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
    
   res <- data.frame(gene = rownames(mut)[k],
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
write.csv(cox_res, file = file.path(dir_output, 'cox_os_IDH_patients.csv'), row.names=FALSE)
