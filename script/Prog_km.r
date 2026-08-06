##-------------------------------------------------------------------
## Script: CNS NGS Mutation and Clinical Variables –
##         Overall Survival Kaplan–Meier Analysis
## Purpose:
##   To evaluate associations between overall survival (OS)
##   and clinical or gene-level mutation variables in CNS
##   tumour patients using Kaplan–Meier survival analyses.
##
##   Survival analyses are performed across the full cohort
##   and IDH-stratified subgroups.
##
## Input:
##   - result/data/mae_mut_clin.RData
##       MultiAssayExperiment object containing:
##         - binary mutation matrix
##         - harmonized clinical metadata
##
## Outputs:
##   1) Kaplan–Meier survival plots for:
##        - Clinical variables
##        - Gene-level mutation status
##   2) Survival plots generated for:
##        - Full cohort
##        - IDH wild-type subgroup
##        - IDH mutant subgroup
##        - IDH wild-type glioblastoma subgroup
##        - IDH wild-type non-glioblastoma subgroup
##
## Processing Overview:
##   1) Load mutation and clinical data
##   2) Harmonize clinical variables including:
##        - Age
##        - Sex
##        - Grade
##        - Histology
##        - ECOG
##        - MGMT
##        - Resection
##   3) Apply administrative censoring to overall survival times
##   4) Perform Kaplan–Meier and log-rank analyses
##   5) Repeat analyses within IDH-stratified patient subgroups
##   6) Generate and export publication-quality survival figures
##
## Notes:
##   - Analyses are based on binary gene-level mutation matrices.
##   - Survival analyses are performed using the survival and survminer frameworks.
##   - Outputs are intended for downstream visualization and reporting.
##-------------------------------------------------------------------
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
library(readxl)

####################################################
## Setup directories
####################################################
dir_input <- 'result/data'
dir_output <- 'result/KM'

time.censor <- 36 # censoring for os including all patients or WT patients
time.censor.mut <- 120 # censoring for os including Mut patients
n1.cutoff <- 5
n0.cutoff <- 5

#################################################
## Load data
#################################################
load(file.path(dir_input, 'mae_mut_clin.RData'))
mut <- assay(mae[['mut_binary']])
clin <- as.data.frame(colData(mae[['mut_binary']]))

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

mut <- mut[, colnames(mut) %in% clin$Study] # 199 patients and 39 gene-mutations

###############################################################
## OS density
###############################################################
pdf(file=file.path(dir_output, "desnity_os.pdf"),
     width = 4, height = 3)

ggplot(clin) +
  geom_density(aes(x = os.time), fill = "#4e4a1c", alpha = 0.4) +
  labs(x = "Time (months)", y = "Density", fill = "") +
  theme_classic()

dev.off()

##############################################################################
## KM: OS
##############################################################################
# Fit models
fit_os  <- survfit(Surv(os.time, os.event) ~ 1, data = clin)

# Plots
pdf(file.path(dir_output, "KM_os.pdf"), width = 4, height = 4)
 ggsurvplot(
  fit_os,
  data = clin,
  conf.int = TRUE,
  palette = c("#4e4a1c"),
  ggtheme = theme_classic(),
  xlab = "Time (months)",
  ylab = "Overall Survival"
)
dev.off()

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
  ylab = "Overall survival probability",
 ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
) 

pdf(file.path(dir_output, 'clinical/all', "KM_Histo_OS.pdf"), width = 5.5, height = 6)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/all", "KM_Histo_OS.jpeg"), width = 5.5, height = 6,
  units = "in", res = 300, quality = 100)
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
  ylab = "Overall survival probability",
  ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
) 

pdf(file.path(dir_output, 'clinical/all', "KM_Grade_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/all", "KM_Grade_OS.jpeg"), width = 5, height = 6,
  units = "in", res = 300, quality = 100)
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
  legend.labs = c("\u200B>40", "<40"),
  palette = c("#08306B", "#C6DBEF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
  ggtheme = theme_bw(base_size = 10) +
    theme(
     panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
) 

pdf(file.path(dir_output, 'clinical/all', "KM_Age_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/all", "KM_Age_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
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
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/all', "KM_Sex_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/all", "KM_Sex_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
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
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/all', "KM_IDH_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/all", "KM_IDH_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

## MGMT
data$variable <- clin$MGMT
data$variable <- factor(
  data$variable,
  levels = c("M", "U", "Unknown")
   )
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "MGMT",
  legend.labs = c("M", "U", 'Unknown'),
  palette = c("#E76254FF", "#376795FF", "#697878FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/all', "KM_MGMT_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/all", "KM_MGMT_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

## ECOG
data$variable <- clin$ECOG
data$variable <- factor(
  data$variable,
  levels = c("0", "1", "2", "3")
   )
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "ECOG",
  legend.labs = c("0", "1", "2", "3"),
  palette = c("#5A4B3CFF", "#4B5A69FF", '#937D61FF', '#D88C27FF'),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/all', "KM_ECOG_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/all", "KM_ECOG_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

# Resection
data$variable <- clin$Resection
data$variable <- factor(
  data$variable,
  levels = c("Subtotal", "Biopsy", "Total")
   )

surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Resection",
  legend.labs = c("Subtotal", "Biopsy", "Total"),
  palette = c("#D04E59FF", "#2F3D70FF", "#BC8E7DFF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/all', "KM_Resection_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/all", "KM_Resection_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
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
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
) 

pdf(file.path(dir_output, 'clinical/wt', "KM_Histo_OS.pdf"), width = 6, height = 6)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/wt", "KM_Histo_OS.jpeg"), width = 6, height = 6,
  units = "in", res = 300, quality = 100)
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
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
) 

pdf(file.path(dir_output, 'clinical/wt', "KM_Grade_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/wt", "KM_Grade_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
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
  legend.labs = c("\u200B>40", "<40"),
  palette = c("#08306B", "#C6DBEF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability", ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
) 

pdf(file.path(dir_output, 'clinical/wt', "KM_Age_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/wt", "KM_Age_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
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
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/wt', "KM_Sex_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/wt", "KM_Sex_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

## MGMT
data$variable <- clin_wt$MGMT
data$variable <- factor(
  data$variable,
  levels = c("M", "U", "Unknown")
   )
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "MGMT",
  legend.labs = c("M", "U", 'Unknown'),
  palette = c("#E76254FF", "#376795FF", "#697878FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/wt', "KM_MGMT_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/wt", "KM_MGMT_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

## ECOG
data$variable <- clin_wt$ECOG
data$variable <- factor(
  data$variable,
  levels = c("0", "1", "2", "3")
   )
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "ECOG",
  legend.labs = c("0", "1", "2", "3"),
  palette = c("#5A4B3CFF", "#4B5A69FF", '#937D61FF', '#D88C27FF'),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/wt', "KM_ECOG_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/wt", "KM_ECOG_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

# Resection
data$variable <- clin_wt$Resection
data$variable <- factor(
  data$variable,
  levels = c("Subtotal", "Biopsy", "Total")
   )

surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Resection",
  legend.labs = c("Subtotal", "Biopsy", "Total"),
  palette = c("#D04E59FF", "#2F3D70FF", "#BC8E7DFF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/wt', "KM_Resection_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/wt", "KM_Resection_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
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
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
) 

pdf(file.path(dir_output, 'clinical/mut', "KM_Histo_OS.pdf"), width = 5, height = 5)
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
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
) 

pdf(file.path(dir_output, 'clinical/mut', "KM_Grade_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/mut", "KM_Grade_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
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
  legend.labs = c("\u200B>40", "<40"),
  palette = c("#08306B", "#C6DBEF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
) 

pdf(file.path(dir_output, 'clinical/mut', "KM_Age_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/mut", "KM_Age_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
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
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/mut', "KM_Sex_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/mut", "KM_Sex_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

## MGMT
data$variable <- clin_mut$MGMT
data$variable <- factor(
  data$variable,
  levels = c("M", "U", "Unknown")
   )
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "MGMT",
  legend.labs = c("M", "U", 'Unknown'),
  palette = c("#E76254FF", "#376795FF", "#697878FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/mut', "KM_MGMT_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/mut", "KM_MGMT_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

## ECOG
data$variable <- clin_mut$ECOG
data$variable <- factor(
  data$variable,
  levels = c("0", "1")
   )
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "ECOG",
  legend.labs = c("0", "1"),
  palette = c("#5A4B3CFF", "#4B5A69FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/mut', "KM_ECOG_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/mut", "KM_ECOG_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

# Resection
data$variable <- clin_mut$Resection
data$variable <- factor(
  data$variable,
  levels = c("Subtotal", "Biopsy", "Total")
   )

surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Resection",
  legend.labs = c("Subtotal", "Biopsy", "Total"),
  palette = c("#D04E59FF", "#2F3D70FF", "#BC8E7DFF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/mut', "KM_Resection_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/mut", "KM_Resection_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

######################################################
## KM figure --- clinical variables (WT GBM patients)
######################################################
clin_wt_gbm <- clin[clin$IDH_status == 'WT' & clin$Histo == "Glioblastoma", ]

data <- data.frame( status=clin_wt_gbm$os.event , time=clin_wt_gbm$os.time)
data$time <- as.numeric(as.character(data$time))
  
for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor ){
      data[ i , "time" ] = time.censor
      data[ i , "status" ] = 0
      
    }
  }  

## Grade status
data$variable <- clin_wt_gbm$Grade
data$variable <- factor(data$variable,
                        levels = c('I', 'II', 'IV'))
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Grade status",
  legend.labs = c("I", "II", "IV"),
  palette = c("#A8C3A0FF", "#BC8E7DFF", "#7C7189FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
) 

pdf(file.path(dir_output, 'clinical/wt_gbm', "KM_Grade_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/wt_gbm", "KM_Grade_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

## Age status
data$variable <- clin_wt_gbm$Age
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
  legend.labs = c("\u200B>40", "<40"),
  palette = c("#08306B", "#C6DBEF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability", ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
) 

pdf(file.path(dir_output, 'clinical/wt_gbm', "KM_Age_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/wt_gbm", "KM_Age_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

## Sex status
data$variable <- clin_wt_gbm$Sex
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
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/wt_gbm', "KM_Sex_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/wt_gbm", "KM_Sex_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

## MGMT
data$variable <- clin_wt_gbm$MGMT
data$variable <- factor(
  data$variable,
  levels = c("M", "U", "Unknown")
   )
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "MGMT",
  legend.labs = c("M", "U", 'Unknown'),
  palette = c("#E76254FF", "#376795FF", "#697878FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/wt_gbm', "KM_MGMT_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/wt_gbm", "KM_MGMT_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

## ECOG
data$variable <- clin_wt_gbm$ECOG
data$variable <- factor(
  data$variable,
  levels = c("0", "1", "2", "3")
   )
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "ECOG",
  legend.labs = c("0", "1", "2", "3"),
  palette = c("#5A4B3CFF", "#4B5A69FF", '#937D61FF', '#D88C27FF'),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/wt_gbm', "KM_ECOG_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/wt_gbm", "KM_ECOG_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

# Resection
data$variable <- clin_wt_gbm$Resection
data$variable <- factor(
  data$variable,
  levels = c("Subtotal", "Biopsy", "Total")
   )

surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Resection",
  legend.labs = c("Subtotal", "Biopsy", "Total"),
  palette = c("#D04E59FF", "#2F3D70FF", "#BC8E7DFF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/wt_gbm', "KM_Resection_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/wt_gbm", "KM_Resection_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()


######################################################
## KM figure --- clinical variables (WT nonGBM patients)
######################################################
clin_wt_nongbm <- clin[clin$IDH_status == 'WT' & clin$Histo != "Glioblastoma", ]

data <- data.frame( status=clin_wt_nongbm$os.event , time=clin_wt_nongbm$os.time)
data$time <- as.numeric(as.character(data$time))
  
for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor ){
      data[ i , "time" ] = time.censor
      data[ i , "status" ] = 0
      
    }
  }  

## Grade status
data$variable <- clin_wt_nongbm$Grade
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
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
) 

pdf(file.path(dir_output, 'clinical/wt_nongbm', "KM_Grade_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/wt_nongbm", "KM_Grade_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

## Age status
data$variable <- clin_wt_nongbm$Age
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
  legend.labs = c("\u200B>40", "<40"),
  palette = c("#08306B", "#C6DBEF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability", ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
) 

pdf(file.path(dir_output, 'clinical/wt_nongbm', "KM_Age_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/wt_nongbm", "KM_Age_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

## Sex status
data$variable <- clin_wt_nongbm$Sex
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
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/wt_nongbm', "KM_Sex_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/wt_nongbm", "KM_Sex_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

## MGMT
data$variable <- clin_wt_nongbm$MGMT
data$variable <- factor(
  data$variable,
  levels = c("M", "U", "Unknown")
   )
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "MGMT",
  legend.labs = c("M", "U", 'Unknown'),
  palette = c("#E76254FF", "#376795FF", "#697878FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/wt_nongbm', "KM_MGMT_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/wt_nongbm", "KM_MGMT_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

## ECOG
data$variable <- clin_wt_nongbm$ECOG
data$variable <- factor(
  data$variable,
  levels = c("0", "1", "2")
   )
surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "ECOG",
  legend.labs = c("0", "1", "2"),
  palette = c("#5A4B3CFF", "#4B5A69FF", '#937D61FF'),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/wt_nongbm', "KM_ECOG_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/wt_nongbm", "KM_ECOG_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

# Resection
data$variable <- clin_wt_nongbm$Resection
data$variable <- factor(
  data$variable,
  levels = c("Subtotal", "Biopsy", "Total")
   )

surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Resection",
  legend.labs = c("Subtotal", "Biopsy", "Total"),
  palette = c("#D04E59FF", "#2F3D70FF", "#BC8E7DFF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical/wt_nongbm', "KM_Resection_OS.pdf"), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical/wt_nongbm", "KM_Resection_OS.jpeg"), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
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
      data[ i , "time" ] = time.censor
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
  legend.labs = c("WT", "Mut"),
  palette = c("#A12A19FF", "#1E466EFF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
) 

pdf(file.path(dir_output, 'mutation/all', paste(paste('KM', colnames(df)[i], sep="_"), '.pdf', sep="")), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "mutation/all", paste(paste('KM', colnames(df)[i], sep="_"), '.jpeg', sep="")), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

}

####################################################
## KM figure --- mutation data (WT patients)
####################################################
clin_wt <- clin[clin$IDH_status == 'WT', ] # 152 WT patients
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
  legend.labs = c("WT", "Mut"),
  palette = c("#A12A19FF", "#1E466EFF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
) 

pdf(file.path(dir_output, 'mutation/wt', paste(paste('KM', colnames(df)[i], sep="_"), '.pdf', sep="")), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "mutation/wt", paste(paste('KM', colnames(df)[i], sep="_"), '.jpeg', sep="")), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
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
  legend.labs = c("WT", "Mut"),
  palette = c("#A12A19FF", "#1E466EFF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
) 

pdf(file.path(dir_output, 'mutation/mut', paste(paste('KM', colnames(df)[i], sep="_"), '.pdf', sep="")), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "mutation/mut", paste(paste('KM', colnames(df)[i], sep="_"), '.jpeg', sep="")), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

}

####################################################
## KM figure --- mutation data (WT GBM patients)
####################################################
clin_wt_gbm <- clin[clin$IDH_status == 'WT' & clin$Histo == 'Glioblastoma', ] # 105 WT GBM patients
data <- data.frame( status=clin_wt_gbm$os.event , time=clin_wt_gbm$os.time)
data$time <- as.numeric(as.character(data$time))
mut_wt_gbm <- mut[, colnames(mut) %in% clin_wt_gbm$Study]
df <- t(mut_wt_gbm)

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
  legend.labs = c("WT", "Mut"),
  palette = c("#A12A19FF", "#1E466EFF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
) 

pdf(file.path(dir_output, 'mutation/wt_gbm', paste(paste('KM', colnames(df)[i], sep="_"), '.pdf', sep="")), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "mutation/wt_gbm", paste(paste('KM', colnames(df)[i], sep="_"), '.jpeg', sep="")), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

}

####################################################
## KM figure --- mutation data (WT nonGBM patients)
####################################################
clin_wt_nongbm <- clin[clin$IDH_status == 'WT' & clin$Histo != 'Glioblastoma', ] # 43 WT nonGBM patients
data <- data.frame( status=clin_wt_nongbm$os.event , time=clin_wt_nongbm$os.time)
data$time <- as.numeric(as.character(data$time))
mut_wt_nongbm <- mut[, colnames(mut) %in% clin_wt_nongbm$Study]
df <- t(mut_wt_nongbm)

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
  legend.labs = c("WT", "Mut"),
  palette = c("#A12A19FF", "#1E466EFF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
   ggtheme = theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ),

  risk.table.height = 0.22,
  risk.table.fontsize = 3
) 

pdf(file.path(dir_output, 'mutation/wt_nongbm', paste(paste('KM', colnames(df)[i], sep="_"), '.pdf', sep="")), width = 5, height = 5)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "mutation/wt_nongbm", paste(paste('KM', colnames(df)[i], sep="_"), '.jpeg', sep="")), width = 5, height = 5,
  units = "in", res = 300, quality = 100)
print(p_idh)
dev.off()

}

#############################################################
## KM figure --- OS association with subtypes/treatment
#############################################################
## step 1 --- all patients: IDH-mut, IDH-wt-gbm and IDH-wt-other
clin$group <- NA
clin$group[clin$IDH_status == "WT" &
           clin$Histo == "Glioblastoma"] <- "IDH-WT-GBM"
clin$group[clin$IDH_status == "WT" &
           clin$Histo != "Glioblastoma"] <- "IDH-WT-nonGBM"
clin$group[clin$IDH_status == "Mut"] <- "IDH-Mut"

clin$group <- factor(
  clin$group,
  levels = c("IDH-WT-GBM",
             "IDH-WT-nonGBM",
             "IDH-Mut")
)

data <- data.frame( status=clin$os.event , time=clin$os.time)
data$time <- as.numeric(as.character(data$time))
    
for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor ){
      data[ i , "time" ] = time.censor
      data[ i , "status" ] = 0
      
    }
  }   

data$variable <- clin$group                      

surv_obj <- with(data, Surv(time, status))
fit_idh <- survfit(surv_obj ~ variable, data = data)
summary(fit_idh)

p_idh <- ggsurvplot(
  fit_idh,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = " ",
  legend.labs = c('IDH-WT-GBM', 'IDH-Mut', 'IDH-WT-nonGBM'),
  palette = c("#A12A19FF", "#1E466EFF", "#787878FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
  ggtheme = theme_bw(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),

      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 6),

      axis.text = element_text(size = 6),
      axis.title = element_text(size = 7),

      plot.title = element_text(size = 8)
    ),

  pval.size = 3,
  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

pdf(file.path(dir_output, 'clinical', "KM_OS_a.pdf"), width =4, height = 4)
print(p_idh)
dev.off()

jpeg(file.path(dir_output, "clinical", "KM_OS_a.jpeg"), width = 4, height = 4, units = "in", res = 300)
print(p_idh)
dev.off()



#############################################################
## KM figure --- OS association with grade
## IDH-WT-nonGBM only
#############################################################
## step 1 --- subset IDH-WT-nonGBM
df <- clin[clin$group == "IDH-WT-nonGBM", ]

## step 2 --- define grade groups
df$grade_group <- NA

df$grade_group[df$Grade %in% c("I", "II")] <- "Low grade"
df$grade_group[df$Grade %in% c("III", "IV")] <- "High grade"

df$grade_group <- factor(
  df$grade_group,
  levels = c("Low grade", "High grade")
)

## check
table(df$grade_group, useNA = "ifany")
table(df$histo, df$grade_group)

## step 3 --- prepare survival data
data <- data.frame(
  status = df$os.event,
  time   = df$os.time
)

data$time <- as.numeric(as.character(data$time))

## step 4 --- apply censoring
for(i in 1:nrow(data)){

  if(
    !is.na(as.numeric(as.character(data[i, "time"]))) &&
    as.numeric(as.character(data[i, "time"])) > time.censor
  ){
    data[i, "time"]   <- time.censor
    data[i, "status"] <- 0
  }
}

data$variable <- df$grade_group

## remove missing grade if any
data <- data[!is.na(data$variable), ]

## step 5 --- KM analysis
surv_obj <- with(data, Surv(time, status))

fit_grade <- survfit(
  surv_obj ~ variable,
  data = data
)

summary(fit_grade)

## step 6 --- plot
p_grade <- ggsurvplot(
  fit_grade,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = " ",
  legend.labs = c(
    "Low grade",
    "High grade"
  ),
  palette = c("#1E466EFF", "#A12A19FF"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",

  ggtheme = theme_bw(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),

      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 6),

      axis.text = element_text(size = 6),
      axis.title = element_text(size = 7),

      plot.title = element_text(size = 8)
    ),

  pval.size = 3,
  risk.table.height = 0.22,
  risk.table.fontsize = 3
)

## step 7 --- save
pdf(
  file.path(dir_output, "clinical", "KM_OS_WT_nonGBM_grade.pdf"),
  width = 4,
  height = 4
)
print(p_grade)
dev.off()

jpeg(
  file.path(dir_output, "clinical", "KM_OS_WT_nonGBM_grade.jpeg"),
  width = 4,
  height = 4,
  units = "in",
  res = 300
)
print(p_grade)
dev.off()
