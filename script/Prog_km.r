#####################################################
## Script: CNS NGS mutation and clinical variable – OS analysis
##
## Purpose:
##   Evaluate associations between clinical variables and 
##   gene-level binary mutation status (0/1 per patient) with 
##   overall survival (OS) using Kaplan–Meier analysis.
##
##   OS is administratively censored at a fixed time horizon 
##   (time.censor, in months).
##
## Input:
##   - result/data/se_mut_bin_clin.RData
##       SummarizedExperiment object containing:
##         • assay: binary mutation matrix (genes x patients; 0/1)
##         • colData: clinical metadata including:
##             - os.time   (OS time in months)
##             - os.event  (1 = death/event, 0 = censored)
##             - demographic and pathological variables
##
## Output:
##   - Kaplan–Meier survival plots (PDF) for:
##       • Clinical variables:
##           - Histology
##           - Tumor location
##           - WHO grade
##           - Age group
##           - Sex
##           - IDH mutation status
##       • Each gene-level mutation (0 = WT, 1 = Mut)
##
## Key steps:
##   1) Load SummarizedExperiment and extract:
##        - Binary mutation matrix
##        - Clinical metadata
##   2) Harmonize and recode clinical variables.
##   3) Apply administrative censoring at `time.censor` months.
##   4) Generate Kaplan–Meier survival curves using:
##        Surv(os.time, os.event) ~ variable
##   5) Export survival plots with risk tables and log-rank p-values.
##
## Notes:
##   - Mutation status is treated as:
##         0 = Wild-type (WT)
##         1 = Mutant (Mut)
##   - Minimum mutation frequency filtering should be applied
##     before interpretation if genes are rare.
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
      WHO.2021.Grade %in% c(1, 2) ~ "I/II",
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
## KM figure --- clinical variables
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

pdf(file.path(dir_output, 'clin', "KM_Histo_OS.pdf"), width = 5, height = 7)
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

pdf(file.path(dir_output, 'clin', "KM_Location_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

## Grade status
data$variable <- clin$Grade
data$variable <- factor(data$variable,
                        levels = c('I/II', 'III', 'IV'))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Grade status",
  legend.labs = c("I/II", "III", "IV"),
  palette = c("#BC8E7DFF", "#FAE093FF", "#7C7189FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'clin', "KM_Grade_OS.pdf"), width = 5, height = 7)
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

pdf(file.path(dir_output, 'clin', "KM_Age_OS.pdf"), width = 5, height = 7)
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

pdf(file.path(dir_output, 'clin', "KM_Sex_OS.pdf"), width = 5, height = 7)
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

pdf(file.path(dir_output, 'clin', "KM_IDH_OS.pdf"), width = 5, height = 7)
print(p_idh)
dev.off()

####################################################
## KM figure --- mutation data
####################################################
df <- t(mut)

for(i in 1:ncol(df)){

data$variable <- factor(
  as.numeric(df[,i]),
  levels = c(0,1),
  labels = c("WT", "Mut")
)

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
  legend.title = "Mutation status",
  legend.labs = c("Mut", "WT"),
  palette = c("#A12A19FF", "#1E466EFF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
) 

pdf(file.path(dir_output, 'mut', paste(paste('KM', colnames(df)[i], sep="_"), '.pdf', sep="")), width = 5, height = 7)
print(p_idh)
dev.off()

}
