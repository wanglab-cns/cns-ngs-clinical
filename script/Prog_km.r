#####################################################
## Script: CNS NGS mutation and clinical variables – 
##         Overall Survival (OS) Kaplan–Meier analysis
##
## Purpose:
##   To evaluate associations between:
##     (1) Clinical variables and OS
##     (2) Gene-level binary mutation status (0/1 per patient) and OS
##   Survival is analyzed using Kaplan–Meier estimation and 
##   log-rank testing.
##   OS is administratively censored at a fixed time horizon 
##   (time.censor, in months), with a separate censoring time 
##   applied for IDH-mutant patients when specified.
## Input:
##   - result/data/se_mut_bin_clin.RData
##     SummarizedExperiment object containing:
##       • assay: binary mutation matrix (genes × patients; 0 = WT, 1 = Mut)
##       • colData: clinical metadata including:
##           - os.time   (OS time in months)
##           - os.event  (1 = death/event, 0 = censored)
##           - Age
##           - Sex
##           - IDH.status
##           - WHO.2021.Grade
##           - Primary.location
##           - Histology variables
## Output:
##   - Kaplan–Meier survival plots (PDF format) for:
##     Clinical variables:
##       • Histology
##       • Tumor location
##       • WHO grade
##       • Age group
##       • Sex
##       • IDH mutation status
##     Gene-level mutation status:
##       • All patients
##       • IDH-WT subgroup
##       • IDH-Mut subgroup
## Key steps:
##   1) Load SummarizedExperiment object and extract:
##        - Binary mutation matrix
##        - Clinical metadata
##   2) Harmonize and recode clinical variables:
##        - Collapse small histology groups (<10 patients)
##        - Standardize location and grade
##        - Recode IDH and age groups
##   3) Apply administrative censoring:
##        - OS times truncated at time.censor (or time.censor.mut)
##        - Status reset to censored (0) beyond censoring threshold
##   4) Generate Kaplan–Meier curves:
##        Surv(os.time, os.event) ~ variable
##   5) For gene-level analyses:
##        - Convert mutation status to WT vs Mut
##        - Skip genes where one group has zero patients
##        - Apply minimum group size filtering if specified
##   6) Export survival plots with:
##        - Risk tables
##        - Log-rank p-values
##        - Consistent color palettes
## Notes:
##   - Mutation status is coded as:
##         0 = Wild-type (WT)
##         1 = Mutant (Mut)
##   - Genes with no variation (all WT or all Mut) are skipped.
##   - Interpretation should consider mutation frequency thresholds
##     (n1.cutoff, n0.cutoff) for stability of KM estimates.
##   - This analysis is univariate (log-rank only).
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
dir_output <- 'result/Fig'

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
                          levels = c("WT", "Mut")) # 191 patients
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
                        levels = c('Glioblastoma', 'Astrocytoma', 'Glioneuronal', 'Oligodendroglioma',  'Other'))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Histology status",
  legend.labs = c("Glioblastoma", "Astrocytoma", "Oligodendroglioma", "Glioneuronal",  "Other"),
  palette = c("#4A7169FF", "#735231FF", "#E3CA97FF", "#99B6BDFF", "#96A5A5FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'clin/all', "KM_Histo_OS.pdf"), width = 5, height = 7)
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

pdf(file.path(dir_output, 'clin/all', "KM_Location_OS.pdf"), width = 5, height = 7)
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

pdf(file.path(dir_output, 'clin/all', "KM_Grade_OS.pdf"), width = 5, height = 7)
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

pdf(file.path(dir_output, 'clin/all', "KM_Age_OS.pdf"), width = 5, height = 7)
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

pdf(file.path(dir_output, 'clin/all', "KM_Sex_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

## IDH status
data$variable <- clin$IDH.status
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

pdf(file.path(dir_output, 'clin/all', "KM_IDH_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()


####################################################
## KM figure --- clinical variables (WT patients)
####################################################
clin_wt <- clin[clin$IDH.status == 'WT', ]
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
                        levels = c('Glioblastoma', 'Astrocytoma', 'Glioneuronal', 'Other'))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Histology status",
  legend.labs = c("Glioblastoma", "Astrocytoma", "Glioneuronal",  "Other"),
  palette = c("#4A7169FF", "#735231FF",  "#99B6BDFF", "#96A5A5FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'clin/wt', "KM_Histo_OS.pdf"), width = 5, height = 7)
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

pdf(file.path(dir_output, 'clin/wt', "KM_Location_OS.pdf"), width = 5, height = 7)
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
  palette = c(,"#A8C3A0FF", "#BC8E7DFF", "#FAE093FF", "#7C7189FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'clin/wt', "KM_Grade_OS.pdf"), width = 5, height = 7)
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

pdf(file.path(dir_output, 'clin/wt', "KM_Age_OS.pdf"), width = 5, height = 7)
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

pdf(file.path(dir_output, 'clin/wt', "KM_Sex_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

####################################################
## KM figure --- clinical variables (Mut patients)
####################################################
clin_mut <- clin[clin$IDH.status == 'Mut', ]
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

pdf(file.path(dir_output, 'clin/mut', "KM_Histo_OS.pdf"), width = 5, height = 7)
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

pdf(file.path(dir_output, 'clin/mut', "KM_Location_OS.pdf"), width = 5, height = 7)
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

pdf(file.path(dir_output, 'clin/mut', "KM_Grade_OS.pdf"), width = 5, height = 7)
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

pdf(file.path(dir_output, 'clin/mut', "KM_Age_OS.pdf"), width = 5, height = 7)
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

pdf(file.path(dir_output, 'clin/mut', "KM_Sex_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

####################################################
## KM figure --- mutation data (all patients)
####################################################
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

pdf(file.path(dir_output, 'mut/all', paste(paste('KM', colnames(df)[i], sep="_"), '.pdf', sep="")), width = 5, height = 7)
print(p_idh)
dev.off()

}

####################################################
## KM figure --- mutation data (WT patients)
####################################################
clin_wt <- clin[clin$IDH.status == 'WT', ]
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

pdf(file.path(dir_output, 'mut/wt', paste(paste('KM', colnames(df)[i], sep="_"), '.pdf', sep="")), width = 5, height = 7)
print(p_idh)
dev.off()

}


####################################################
## KM figure --- mutation data (Mut patients)
####################################################
clin_mut <- clin[clin$IDH.status == 'Mut', ]
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

pdf(file.path(dir_output, 'mut/Mut', paste(paste('KM', colnames(df)[i], sep="_"), '.pdf', sep="")), width = 5, height = 7)
print(p_idh)
dev.off()

}
