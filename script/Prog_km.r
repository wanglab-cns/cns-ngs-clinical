#####################################################
## Script: CNS NGS Mutation and Clinical Variables –
##         Overall Survival (OS) Kaplan–Meier Analysis
##
## Purpose:
##   Perform univariate survival analyses to evaluate
##   associations between overall survival (OS) and:
##
##     (1) Clinical variables
##     (2) Gene-level binary mutation status (0/1 per patient)
##
##   Survival distributions are estimated using the
##   Kaplan–Meier method and compared using the log-rank test.
##
##   Administrative censoring is applied at predefined
##   time horizons to improve interpretability of survival curves:
##
##       - time.censor      (all patients / IDH-WT subgroup)
##       - time.censor.mut  (IDH-mutant subgroup)
##
##
## Input:
##   - result/data/se_mut_bin_clin.RData
##
##     SummarizedExperiment containing:
##       • assay:
##           Binary mutation matrix (genes × patients)
##           0 = Wild-type (WT), 1 = Mutant (Mut)
##       • colData:
##           Curated clinical metadata including:
##             - os.time        (overall survival time, months)
##             - os.event       (1 = death, 0 = censored)
##             - Age
##             - Sex
##             - IDH_status
##             - WHO.2021.Grade
##             - Primary.location
##             - Histology variables
##             - Therapy_status
##
##
## Output:
##
##   Kaplan–Meier survival plots (PDF format):
##
##   Clinical variables:
##       • Histology
##       • Tumor location
##       • WHO grade
##       • Age group
##       • Sex
##       • IDH mutation status
##       • Therapy status
##
##   Gene-level mutation status:
##       • All patients
##       • IDH-WT subgroup
##       • IDH-Mut subgroup
##
##   Output directories:
##       result/KM/clinical/all
##       result/KM/clinical/wt
##       result/KM/clinical/mut
##       result/KM/mutation/all
##       result/KM/mutation/wt
##       result/KM/mutation/mut
##
##
## Key Processing Steps:
##
##   1) Clinical Harmonization
##        - Exclude samples with missing IDH status
##        - Recode IDH (WT / Mut)
##        - Dichotomize age (<40 vs ≥40)
##        - Standardize tumor location and WHO grade
##        - Harmonize histology labels
##        - Collapse rare histology groups (<10 patients)
##
##   2) Administrative Censoring
##        - Truncate OS times at predefined thresholds
##        - Reset events beyond thresholds to censored
##        - Apply subgroup-specific censoring limits
##
##   3) Kaplan–Meier Estimation
##        - Surv(os.time, os.event) ~ variable
##        - Log-rank test for group comparison
##        - Risk tables displayed
##
##   4) Gene-Level Analysis
##        - Iterate across all genes
##        - Convert mutation status to WT vs Mut
##        - Skip genes with:
##              • No variation (all WT or all Mut)
##              • Group sizes below n1.cutoff / n0.cutoff
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
dir_output <- 'result/KM'

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

clin$Age <- ifelse(clin$Age >= 40, '>40', '<40')
clin$Therapy_status <- ifelse(clin$Therapy_status == 1, 'Yes', 'No')

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
# clin$Histo <- coalesce(clin$LGG, clin$HGG)
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

histo_counts <- table(clin$Histo)

# Identify small groups (<10 patients)
small_groups <- names(histo_counts[histo_counts <= 10])

clin <- clin %>%
  mutate(
    Histo = ifelse(
      Histo %in% small_groups,
      "Other",
      Histo
    )
  )

mut <- mut[, colnames(mut) %in% clin$Study] # 199 patients

####################################################
## KM figure --- clinical variables (all patients)
####################################################
data <- data.frame( status=clin$os.event , time=clin$os.time)
data$time <- as.numeric(as.character(data$time))
  
for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor ){
      data[ i , "time" ] = time.censor
      data[ i , "status" ] = 0
      
    }
  }  

## Histology status
data$variable <- clin$Histo
data$variable <- factor(data$variable,
                        levels = c("Glioblastoma", "Astrocytoma", "Oligodendroglioma", 
                                   "Glioneuronal tumor", "Circumscribed glioma", "Other"))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Histology status",
  legend.labs = c("Glioblastoma", "Astrocytoma", "Oligodendroglioma", 
                   "Glioneuronal tumor", "Circumscribed glioma", "Other"),
  palette = c("#4A7169FF", "#735231FF", "#E3CA97FF", "#99B6BDFF", "#E76254FF", "#96A5A5FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'clinical/all', "KM_Histo_OS.pdf"), width = 7, height = 8)
print(p_idh)
dev.off()

## Location status
data$variable <- clin$Location
data$variable <- factor(data$variable,
                        levels = c('Lobar', 'Cerebellum', 'Thalamic', 'Other'))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Location status",
  legend.labs = c("Lobar", "Cerebellum", "Thalamic", "Other"),
  palette = c("#A54B2DFF", "#577E2FFF", "#B49696FF", "#96A5A5FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'clinical/all', "KM_Location_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

## Grade status
data$variable <- clin$Grade
data$variable <- factor(data$variable,
                        levels = c('I', 'II', 'III', 'IV'))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Grade status",
  legend.labs = c("I", "II", "III", "IV"),
  palette = c("#A8C3A0FF", "#BC8E7DFF", "#FAE093FF", "#7C7189FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'clinical/all', "KM_Grade_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

## Age status
data$variable <- clin$Age
data$variable <- factor(data$variable,
                        levels = c('>40', '<40'))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Age status",
  legend.labs = c("younger", "older"),
  palette = c("#08306B", "#C6DBEF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'clinical/all', "KM_Age_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

## Sex status
data$variable <- clin$Sex
data$variable <- factor(data$variable,
                        levels = c('Female', 'Male'))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Sex status",
  legend.labs = c("Female", "Male"),
  palette = c("#4C72B0", "#DD8452"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
)

pdf(file.path(dir_output, 'clinical/all', "KM_Sex_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

## IDH status
data$variable <- clin$IDH_status
data$variable <- factor(data$variable,
                        levels = c('WT', 'Mut'))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "IDH status",
  legend.labs = c("WT", "Mut"),
  palette = c("#C44E52", "#55A868"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
)

pdf(file.path(dir_output, 'clinical/all', "KM_IDH_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

## Therapy status
data$variable <- clin$Therapy_status
data$variable <- factor(data$variable,
                          levels = c("Yes", "No"))

surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Therapy status",
  legend.labs = c("Yes", "No"),
  palette = c("#846D86FF", "#ABB2A5FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
)


pdf(file.path(dir_output, 'clinical/all', "KM_Therapy_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

####################################################
## KM figure --- clinical variables (WT patients)
####################################################
clin_wt <- clin[clin$IDH_status == 'WT', ]
data <- data.frame( status=clin_wt$os.event , time=clin_wt$os.time)
data$time <- as.numeric(as.character(data$time))
  
for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor ){
      data[ i , "time" ] = time.censor
      data[ i , "status" ] = 0
      
    }
  }  

## Histology status
data$variable <- clin_wt$Histo
data$variable <- factor(data$variable,
                        levels = c("Glioblastoma", "Astrocytoma", 
                                   "Glioneuronal tumor", "Circumscribed glioma", "Other"))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Histology status",
  legend.labs = c("Glioblastoma", "Astrocytoma",  
                   "Glioneuronal tumor", "Circumscribed glioma", "Other"),
  palette = c("#4A7169FF", "#735231FF", "#99B6BDFF", "#E76254FF", "#96A5A5FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'clinical/wt', "KM_Histo_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

## Location status
data$variable <- clin_wt$Location
data$variable <- factor(data$variable,
                        levels = c('Lobar', 'Cerebellum', 'Thalamic', 'Other'))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Location status",
  legend.labs = c("Lobar", "Cerebellum", "Thalamic", "Other"),
  palette = c("#A54B2DFF", "#577E2FFF", "#B49696FF", "#96A5A5FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'clinical/wt', "KM_Location_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

## Grade status
data$variable <- clin_wt$Grade
data$variable <- factor(data$variable,
                        levels = c('I', 'II', 'III', 'IV'))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Grade status",
  legend.labs = c("I", "II", "III", "IV"),
  palette = c("#A8C3A0FF", "#BC8E7DFF", "#FAE093FF", "#7C7189FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'clinical/wt', "KM_Grade_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

## Age status
data$variable <- clin_wt$Age
data$variable <- factor(data$variable,
                        levels = c('>40', '<40'))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Age status",
  legend.labs = c("younger", "older"),
  palette = c("#08306B", "#C6DBEF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'clinical/wt', "KM_Age_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

## Sex status
data$variable <- clin_wt$Sex
data$variable <- factor(data$variable,
                        levels = c('Female', 'Male'))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Sex status",
  legend.labs = c("Female", "Male"),
  palette = c("#4C72B0", "#DD8452"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
)

pdf(file.path(dir_output, 'clinical/wt', "KM_Sex_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

## Therapy status
data$variable <- clin_wt$Therapy_status
data$variable <- factor(data$variable,
                          levels = c("Yes", "No"))

surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Therapy status",
  legend.labs = c("Yes", "No"),
  palette = c("#846D86FF", "#ABB2A5FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
)

pdf(file.path(dir_output, 'clinical/wt', "KM_Therapy_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

####################################################
## KM figure --- clinical variables (Mut patients)
####################################################
clin_mut <- clin[clin$IDH_status == 'Mut', ]
data <- data.frame( status=clin_mut$os.event , time=clin_mut$os.time)
data$time <- as.numeric(as.character(data$time))
  
for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor.mut ){
      data[ i , "time" ] = time.censor.mut
      data[ i , "status" ] = 0
      
    }
  }  

## Histology status
data$variable <- clin_mut$Histo
data$variable <- factor(data$variable,
                        levels = c('Astrocytoma', 'Oligodendroglioma'))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Histology status",
  legend.labs = c("Astrocytoma", "Oligodendroglioma"),
  palette = c( "#735231FF", "#E3CA97FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'clinical/mut', "KM_Histo_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

## Location status
data$variable <- clin_mut$Location
data$variable <- factor(data$variable,
                        levels = c('Lobar', 'Cerebellum'))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Location status",
  legend.labs = c("Lobar", "Cerebellum"),
  palette = c("#A54B2DFF", "#577E2FFF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'clinical/mut', "KM_Location_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

## Grade status
data$variable <- clin_mut$Grade
data$variable <- factor(data$variable,
                        levels = c('II', 'III', 'IV'))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Grade status",
  legend.labs = c("II", "III", "IV"),
  palette = c( "#BC8E7DFF", "#FAE093FF", "#7C7189FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'clinical/mut', "KM_Grade_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

## Age status
data$variable <- clin_mut$Age
data$variable <- factor(data$variable,
                        levels = c('>40', '<40'))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Age status",
  legend.labs = c("younger", "older"),
  palette = c("#08306B", "#C6DBEF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'clinical/mut', "KM_Age_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

## Sex status
data$variable <- clin_mut$Sex
data$variable <- factor(data$variable,
                        levels = c('Female', 'Male'))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Sex status",
  legend.labs = c("Female", "Male"),
  palette = c("#4C72B0", "#DD8452"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
)

pdf(file.path(dir_output, 'clinical/mut', "KM_Sex_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

## Therapy status
data$variable <- clin_mut$Therapy_status
data$variable <- factor(data$variable,
                          levels = c("Yes", "No"))

surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Therapy status",
  legend.labs = c("Yes", "No"),
  palette = c("#846D86FF", "#ABB2A5FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
)


pdf(file.path(dir_output, 'clinical/mut', "KM_Therapy_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

####################################################
## KM figure --- mutation data (all patients)
####################################################
## ---- Binary matrix
if(all(c("IDH1","IDH2") %in% rownames(mut))){

  idh_vec <- as.integer(
    mut["IDH1", ] == 1 | mut["IDH2", ] == 1
  )
  mut <- mut[!rownames(mut) %in% c("IDH1","IDH2"), , drop = FALSE]
  mut <- rbind(mut, IDH = idh_vec)

}

data <- data.frame( status=clin$os.event , time=clin$os.time)
data$time <- as.numeric(as.character(data$time))
df <- t(mut)
  
for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor ){
      data[ i , "time" ] = time.censor.mut
      data[ i , "status" ] = 0
      
    }
  }  


for(i in 1:ncol(df)){

data$variable <- factor(
  as.numeric(df[,i]),
  levels = c(0,1),
  labels = c("WT", "Mut")
)

data$variable <- factor(data$variable,
                        levels = c('WT', 'Mut'))

grp_counts <- table(data$variable)
  
  if (length(grp_counts) < 2 || any(grp_counts == 0)) {
    message(paste("Skipping", i, "- only one group present"))
    next
  }

surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Mutation status",
  legend.labs = c("Mut", "WT"),
  palette = c("#A12A19FF", "#1E466EFF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'mutation/all', paste(paste('KM', colnames(df)[i], sep="_"), '.pdf', sep="")), width = 5, height = 7)
print(p_idh)
dev.off()

}

####################################################
## KM figure --- mutation data (WT patients)
####################################################
clin_wt <- clin[clin$IDH_status == 'WT', ]
data <- data.frame( status=clin_wt$os.event , time=clin_wt$os.time)
data$time <- as.numeric(as.character(data$time))
mut_wt <- mut[, colnames(mut) %in% clin_wt$Study]
df <- t(mut_wt)

for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor ){
      data[ i , "time" ] = time.censor
      data[ i , "status" ] = 0
      
    }
  }  


for(i in 1:ncol(df)){

print(i)
data$variable <- factor(
  as.numeric(df[,i]),
  levels = c(0,1),
  labels = c("WT", "Mut")
)

data$variable <- factor(data$variable,
                        levels = c('WT', 'Mut'))

grp_counts <- table(data$variable)
  
  if (length(grp_counts) < 2 || any(grp_counts == 0)) {
    message(paste("Skipping", i, "- only one group present"))
    next
  }
  
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Mutation status",
  legend.labs = c("Mut", "WT"),
  palette = c("#A12A19FF", "#1E466EFF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'mutation/wt', paste(paste('KM', colnames(df)[i], sep="_"), '.pdf', sep="")), width = 5, height = 7)
print(p_idh)
dev.off()

}


####################################################
## KM figure --- mutation data (Mut patients)
####################################################
clin_mut <- clin[clin$IDH_status == 'Mut', ]
data <- data.frame( status=clin_mut$os.event , time=clin_mut$os.time)
data$time <- as.numeric(as.character(data$time))
mut_mut <- mut[, colnames(mut) %in% clin_mut$Study]
df <- t(mut_mut)

for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor.mut ){
      data[ i , "time" ] = time.censor.mut
      data[ i , "status" ] = 0
      
    }
  }  


for(i in 1:ncol(df)){

print(i)
data$variable <- factor(
  as.numeric(df[,i]),
  levels = c(0,1),
  labels = c("WT", "Mut")
)

data$variable <- factor(data$variable,
                        levels = c('WT', 'Mut'))

grp_counts <- table(data$variable)
  
  if (length(grp_counts) < 2 || any(grp_counts == 0)) {
    message(paste("Skipping", i, "- only one group present"))
    next
  }
  
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Mutation status",
  legend.labs = c("Mut", "WT"),
  palette = c("#A12A19FF", "#1E466EFF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'mutation/mut', paste(paste('KM', colnames(df)[i], sep="_"), '.pdf', sep="")), width = 5, height = 7)
print(p_idh)
dev.off()

}
