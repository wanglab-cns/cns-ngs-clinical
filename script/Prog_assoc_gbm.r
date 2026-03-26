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
dir_output <- 'result/assoc_gbm'

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

## ---- Binary matrix
if(all(c("IDH1","IDH2") %in% rownames(mut))){

  idh_vec <- as.integer(
    mut["IDH1", ] == 1 | mut["IDH2", ] == 1
  )
  mut <- mut[!rownames(mut) %in% c("IDH1","IDH2"), , drop = FALSE]
  mut <- rbind(mut, IDH = idh_vec)

}

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
small_groups <- names(histo_counts[histo_counts < 10])

clin <- clin %>%
  mutate(
    Histo = ifelse(
      Histo %in% small_groups,
      "Other",
      Histo
    )
  )

clin$Sex <- factor(
  clin$Sex,
  levels = c("Female", "Male")
)

clin$Age <- factor(
  clin$Age,
  levels = c(">40", "<40")
)

clin <- clin[clin$Histo == 'Glioblastoma', ]
mut <- mut[, colnames(mut) %in% clin$Study] # 105 patients
clin$group <- case_when(
  clin$os.time < 6 ~ "Short",
  clin$os.time > 24 ~ "Long",
  TRUE ~ NA_character_
)

clin <- clin %>%
  filter(!is.na(group)) # 53

clin$group_bin <- ifelse(clin$group == "Long", 1, 0)
clin <- clin[clin$MGMT != 'Unknown', ]
mut <- mut[, colnames(mut) %in% clin$Study] # 105 patients

####################################################
## association --> clinical
####################################################
vars <- c("Age", "Sex", "ECOG", "MGMT", "Resection")

results <- data.frame()

for (v in vars) {
  
  fit <- glm(as.formula(paste("group_bin ~", v)),
             data = clin,
             family = "binomial")
  
  coef_table <- summary(fit)$coefficients
  
  # Skip intercept, take variable rows
  for (i in 2:nrow(coef_table)) {
    
    est <- coef_table[i, "Estimate"]
    se  <- coef_table[i, "Std. Error"]
    p   <- coef_table[i, "Pr(>|z|)"]
    
    results <- rbind(results, data.frame(
      variable = v,
      level = rownames(coef_table)[i],
      OR = exp(est),
      CI_low = exp(est - 1.96 * se),
      CI_high = exp(est + 1.96 * se),
      p_value = p
    ))
  }
}

write.csv(results, file = file.path(dir_output, 'clinical.csv'), row.names=FALSE)

####################################################
## association --> mutations
####################################################
genes <- rownames(mut)

mut_results <- lapply(genes, function(g) {
  print(g)
  df <- data.frame(
    group = clin$group_bin,
    mut_bin = mut[g, clin$Study]
  )
  
  # remove NAs
  df <- na.omit(df)
  
  # skip genes with no variation
  if (length(unique(df$mut_bin)) < 2) {
    return(NULL)
  }

  fit <- glm(group ~ mut_bin, data = df, family = "binomial")
  
  coef <- summary(fit)$coefficients
  
  data.frame(
    gene = g,
    OR = exp(coef[2,1]),
    pval = coef[2,4]
  )
})

mut_results <- do.call(rbind, mut_results)
mut_results$fdr <- p.adjust(mut_results$pval, method = "fdr")

write.csv(mut_results, file = file.path(dir_output, 'mutation.csv'), row.names=FALSE)

####################################################
## association --> mutations + clinical
####################################################
genes <- rownames(mut)

mut_results <- lapply(genes, function(g) {
  print(g)
  df <- data.frame(
    group = clin$group_bin,
    mut_bin = mut[g, clin$Study],
    Age = clin$Age,
    Sex = clin$Sex,
    ECOG = clin$ECOG,
    MGMT = clin$MGMT,
    Resection = clin$Resection
  )
  

  # remove NAs
  df <- na.omit(df)
  
  # skip genes with no variation
  if (length(unique(df$mut_bin)) < 2) {
    return(NULL)
  }
 
   # skip rare mutations
  if (sum(df$mut_bin == 1) < 3) {
    return(NULL)
  }

  fit <- glm(group ~ mut_bin + Age + Sex + ECOG + MGMT + Resection, data = df, family = "binomial")
  
  coef <- summary(fit)$coefficients
  
  data.frame(
    gene = g,
    OR = exp(coef[2,1]),
    pval = coef[2,4]
  )
})

mut_results <- do.call(rbind, mut_results)
mut_results$fdr <- p.adjust(mut_results$pval, method = "fdr")


write.csv(mut_results, file = file.path(dir_output, 'mut_clinical.csv'), row.names=FALSE)


####################################################
## association --> mutations + clinical
####################################################
genes <- rownames(mut)

mut_results <- lapply(genes, function(g) {
  print(g)

  df <- data.frame(
    group = clin$group_bin,
    mut_bin = mut[g, clin$Study],
    #Age = clin$Age,
    #Sex = clin$Sex,
    ECOG = clin$ECOG,
    #MGMT = clin$MGMT,
    Resection = clin$Resection
  )
  

  # remove NAs
  df <- na.omit(df)
  
  # skip genes with no variation
  if (length(unique(df$mut_bin)) < 2) {
    return(NULL)
  }
 
   # skip rare mutations
  if (sum(df$mut_bin == 1) < 3) {
    return(NULL)
  }

  fit <- glm(group ~ mut_bin + ECOG + Resection, data = df, family = "binomial")
  
  coef <- summary(fit)$coefficients
  
  data.frame(
    gene = g,
    OR = exp(coef[2,1]),
    pval = coef[2,4]
  )
})

mut_results <- do.call(rbind, mut_results)
mut_results$fdr <- p.adjust(mut_results$pval, method = "fdr")

write.csv(mut_results, file = file.path(dir_output, 'mut_sig_clinical.csv'), row.names=FALSE)

