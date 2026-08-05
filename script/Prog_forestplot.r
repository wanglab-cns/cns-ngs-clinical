##-----------------------------------------------------------------
## Script: Forest Plot Visualization of Survival Analyses
##
## Purpose:
##   To generate forest plots summarizing gene-level
##   overall survival associations in CNS tumours.
##
##   Forest plots compare mutation-only and mutation +
##   clinical survival models in IDH wild-type GBM and
##   IDH mutant subgroups.
##
## Input:
##
##   - Mutation survival analysis results
##   - Mutation + clinical survival analysis results
##
## Outputs:
##
##   - result/assoc/forestplot_WT_GBM.pdf
##   - result/assoc/forestplot_Mut.pdf
##
## Processing Overview:
##
##   1) Load survival analysis results
##   2) Prepare hazard ratio and confidence interval
##      matrices
##   3) Generate forest plot tables and visualizations
##   4) Export publication-quality PDF figures
##
## Notes:
##   - Forest plots display hazard ratios with
##     95% confidence intervals.
##   - Analyses are stratified by IDH status.
##-----------------------------------------------------------------
###########################################################
## library
###########################################################
library(PredictioR)
library(ggplot2)
library(ComplexHeatmap)
library(forestplot)
library(dplyr)
library(data.table)
library(MultiAssayExperiment)

###########################################################
## set up directory
###########################################################
dir_mut <-  "result/assoc/mutation"
dir_clin_mut <- "result/assoc/"
dir_output <- "result/assoc"

#########################################################
## Mutation: OS
#########################################################
os <- read.csv(file.path(dir_mut, "cox_os_mut_wt_gbm.csv"))
os_mut_wt_gbm_sig <- os[os$pval < 0.5, ]
os_mut_wt_gbm <- os

os <- read.csv(file.path(dir_mut, "cox_os_mut_mut.csv"))
os_mut_mut_sig <-  os[os$pval < 0.5, ]
#os_mut_mut_sig <- os
os_mut_mut <- os

#os <- read.csv(file.path(dir_mut, "cox_os_mut_wt_nongbm.csv"))
#os_mut_wt_nongbm_sig <-  os[os$pval < 0.5, ]
#os_mut_wt_nongbm <- os

#########################################################
## Mutation & clinical: OS
#########################################################
os <- read.csv(file.path(dir_clin_mut, 'mutation & clinical', "cox_os_wt_gbm_mv.csv"))
os_mut_clin_wt_gbm_sig <-  os[os$pval < 0.5, ]
os_mut_clin_wt_gbm <-  os

os <- read.csv(file.path(dir_clin_mut, 'mutation & clinical', "cox_os_mut_mv.csv"))
os_mut_clin_mut_sig <-  os[os$pval < 0.5, ]
os_mut_clin_mut <-  os

######################################################
## prepare data for forestplot:  WT GBM
######################################################
int <- union( os_mut_wt_gbm_sig$gene, os_mut_clin_wt_gbm_sig$gene)

os_mut_wt_gbm_upd <- lapply(1:length(int), function(k){
  
  res <- os_mut_wt_gbm[os_mut_wt_gbm$gene == int[k], ]
  if(nrow(res) == 0){
    
    res <- data.frame(gene = int[k],
                      logHR = NA,
                      HR = NA,
                      se = NA,
                      n = NA, 
                      low = NA,
                      up = NA,
                      pval = NA,
                      fdr = NA,
                      study = 'Mutation')
    
  }
  
  res
  
})

os_mut_wt_gbm_upd <- do.call(rbind, os_mut_wt_gbm_upd)
os_mut_wt_gbm_upd$Coef <- os_mut_wt_gbm_upd$HR
os_mut_wt_gbm_upd$Lower <- os_mut_wt_gbm_upd$low
os_mut_wt_gbm_upd$Upper <- os_mut_wt_gbm_upd$up
os_mut_wt_gbm_upd$N <- os_mut_wt_gbm_upd$n

os_mut_clin_wt_gbm_upd <- lapply(1:length(int), function(k){
  
  res <- os_mut_clin_wt_gbm[os_mut_clin_wt_gbm$gene == int[k], ]
  if(nrow(res) == 0){
    
 res <- data.frame(gene = int[k],
                   level = NA, 
                   logHR = NA,
                   HR = NA,
                   se = NA,
                   n = NA, 
                   low = NA,
                   up = NA,
                   pval = NA,
                   fdr = NA,
                   study = 'Mutation & Clinical')
    
  }
  
  res
  
})

os_mut_clin_wt_gbm_upd <- do.call(rbind, os_mut_clin_wt_gbm_upd)
os_mut_clin_wt_gbm_upd$Coef <- os_mut_clin_wt_gbm_upd$HR
os_mut_clin_wt_gbm_upd$Lower <- os_mut_clin_wt_gbm_upd$low
os_mut_clin_wt_gbm_upd$Upper <- os_mut_clin_wt_gbm_upd$up
os_mut_clin_wt_gbm_upd$N <- os_mut_clin_wt_gbm_upd$n

dat <- lapply(1:length(int), function(k){
  
  data.frame(gene = int[k],
             #N = paste(df_pan_upd$N[k], df_male_upd$N[k], df_female_upd$N[k], sep="/"),
             mut_Effect = os_mut_wt_gbm_upd$Coef[k],
             mut_Lower = os_mut_wt_gbm_upd$Lower[k],
             mut_Upper = os_mut_wt_gbm_upd$Upper[k],
             mut_clin_Effect = os_mut_clin_wt_gbm_upd$Coef[k],
             mut_clin_Lower = os_mut_clin_wt_gbm_upd$Lower[k],
             mut_clin_Upper = os_mut_clin_wt_gbm_upd$Upper[k])
})

dat <- do.call(rbind, dat)

## create format of forestplot
table_text <- rbind(
  c("Gene", "Mut", "Mut & Clin"),
  cbind(dat$gene, 
        formatC(dat$N, format="d"), 
        ifelse(is.na(dat$mut_Effect), " ", formatC(dat$mut_Effect, format="f", digits=2)),
        ifelse(is.na(dat$mut_clin_Effect), " ", formatC(dat$mut_clin_Effect, format="f", digits=2)))
)

# prepare effect size matrices
mean_values <- rbind(
  c(NA, NA),  # header
  cbind(dat$mut_Effect, dat$mut_clin_Effect)
)

lower_values <- rbind(
  c(NA, NA),
  cbind(dat$mut_Lower, dat$mut_clin_Lower)
)

upper_values <- rbind(
  c(NA, NA),
  cbind(dat$mut_Upper, dat$mut_clin_Upper)
)

# generate the forest plot

p <- forestplot(
  
  # Table labels
  labeltext = table_text,
  
  # Effect estimates
  mean  = mean_values,
  lower = lower_values,
  upper = upper_values,

  # Summary rows
  is.summary = c(TRUE, rep(FALSE, nrow(dat))),

  # Reference line
  zero = 1,
  lwd.zero = 1,

  # Axis
  xlab = "Hazard Ratio",
  xlog = FALSE,
  xticks = c(0, 0.5, 1, 2, 4),
  clip = c(0, 5),

  # Colors
  col = fpColors(
    box   = c("#526A83FF", "#C44E52"),      
    lines = c("#526A83FF", "#C44E52"),
    zero  = "#7C7C7CFF"
  ),

  # CI appearance
  lwd.ci = 2,
  ci.vertices = TRUE,
  ci.vertices.height = 0.08,

  # Box appearance
  boxsize = 0.18,

  # Spacing
  lineheight = unit(0.72, "cm"),
  mar = unit(c(3, 5, 2, 2), "mm"),

  # Text formatting
 txt_gp = fpTxtGp(
    label = gpar(fontsize = 8),
    ticks = gpar(fontsize = 10, fontface = "bold"),
    xlab  = gpar(fontsize = 12, fontface = "bold"),
    title = gpar(fontsize = 12, fontface = "bold")
  ),

  # Horizontal separator
  hrzl_lines = list(
    "2" = gpar(lwd = 1, col = "#272729FF")
  )

  # Title
  # title = "IDH Wild-Type Cohort"
)

pdf(file=file.path(dir_clin_mut, "forestplot_WT_GBM.pdf"), 
    width = 4.5, height = 3.5, useDingbats = FALSE)
print(p)
dev.off()


jpeg(file=file.path(dir_clin_mut, "forestplot_WT_GBM.jpeg"), 
    width = 4.5, height = 3.5, units = "in", res = 300, quality = 100)
print(p)
dev.off()

######################################################
## prepare data for forestplot: Mut
######################################################
int <- intersect(os_mut_mut$gene, union( os_mut_mut_sig$gene, os_mut_clin_mut_sig$gene))

os_mut_mut_upd <- lapply(1:length(int), function(k){
  
  res <- os_mut_mut[os_mut_mut$gene == int[k], ]
  if(nrow(res) == 0){
    
    res <- data.frame(gene = int[k],
                      logHR = NA,
                      HR = NA,
                      se = NA,
                      n = NA, 
                      low = NA,
                      up = NA,
                      pval = NA,
                      fdr = NA,
                      study = 'Mutation')
    
  }
  
  res
  
})

os_mut_mut_upd <- do.call(rbind, os_mut_mut_upd)
os_mut_mut_upd$Coef <- os_mut_mut_upd$HR
os_mut_mut_upd$Lower <- os_mut_mut_upd$low
os_mut_mut_upd$Upper <- os_mut_mut_upd$up
os_mut_mut_upd$N <- os_mut_mut_upd$n

os_mut_clin_mut_upd <- lapply(1:length(int), function(k){
  
  res <- os_mut_clin_mut[os_mut_clin_mut$gene == int[k], ]
  if(nrow(res) == 0){
    
 res <- data.frame(gene = int[k],
                   level = NA,
                   logHR = NA,
                   HR = NA,
                   se = NA,
                   n = NA, 
                   low = NA,
                   up = NA,
                   pval = NA,
                   fdr = NA,
                   study = 'Mutation & Clinical')
    
  }
  
  res
  
})

os_mut_clin_mut_upd <- do.call(rbind, os_mut_clin_mut_upd)
os_mut_clin_mut_upd$Coef <- os_mut_clin_mut_upd$HR
os_mut_clin_mut_upd$Lower <- os_mut_clin_mut_upd$low
os_mut_clin_mut_upd$Upper <- os_mut_clin_mut_upd$up
os_mut_clin_mut_upd$N <- os_mut_clin_mut_upd$n

dat <- lapply(1:length(int), function(k){
  
  data.frame(gene = int[k],
             #N = paste(df_pan_upd$N[k], df_male_upd$N[k], df_female_upd$N[k], sep="/"),
             mut_Effect = os_mut_mut_upd$Coef[k],
             mut_Lower = os_mut_mut_upd$Lower[k],
             mut_Upper = os_mut_mut_upd$Upper[k],
             mut_clin_Effect = os_mut_clin_mut_upd$Coef[k],
             mut_clin_Lower = os_mut_clin_mut_upd$Lower[k],
             mut_clin_Upper = os_mut_clin_mut_upd$Upper[k])
})

dat <- do.call(rbind, dat)

## create format of forestplot
table_text <- rbind(
  c("Gene", "Mut", "Mut & Clin"),
  cbind(dat$gene, 
        formatC(dat$N, format="d"), 
        ifelse(is.na(dat$mut_Effect), " ", formatC(dat$mut_Effect, format="f", digits=2)),
        ifelse(is.na(dat$mut_clin_Effect), " ", formatC(dat$mut_clin_Effect, format="f", digits=2)))
)

# prepare effect size matrices
mean_values <- rbind(
  c(NA, NA),  # header
  cbind(dat$mut_Effect, dat$mut_clin_Effect)
)

lower_values <- rbind(
  c(NA, NA),
  cbind(dat$mut_Lower, dat$mut_clin_Lower)
)

upper_values <- rbind(
  c(NA, NA),
  cbind(dat$mut_Upper, dat$mut_clin_Upper)
)


# generate the forest plot
p <- forestplot(
  
  # Table labels
  labeltext = table_text,
  
  # Effect estimates
  mean  = mean_values,
  lower = lower_values,
  upper = upper_values,

  # Summary rows
  is.summary = c(TRUE, rep(FALSE, nrow(dat))),

  # Reference line
  zero = 1,
  lwd.zero = 1,

  # Axis
  xlab = "Hazard Ratio",
  xlog = FALSE,
  xticks = c(0, 1, 2, 5, 10),
  clip = c(0, 10),

  # Colors
  col = fpColors(
    box   = c("#526A83FF", "#55A868"),      
    lines = c("#526A83FF", "#55A868"),
    zero  = "#7C7C7CFF"
  ),

  # CI appearance
  lwd.ci = 2,
  ci.vertices = TRUE,
  ci.vertices.height = 0.08,

  # Box appearance
  boxsize = 0.18,

  # Spacing
  lineheight = unit(0.72, "cm"),
  mar = unit(c(3, 5, 2, 2), "mm"),

  # Text formatting
  txt_gp = fpTxtGp(
    label = gpar(fontsize = 8),
    ticks = gpar(fontsize = 10, fontface = "bold"),
    xlab  = gpar(fontsize = 12, fontface = "bold"),
    title = gpar(fontsize = 12, fontface = "bold")
  ),

  # Horizontal separator
  hrzl_lines = list(
    "2" = gpar(lwd = 1, col = "#272729FF")
  )

  # Title
  # title = "IDH Wild-Type Cohort"
)

pdf(file=file.path(dir_clin_mut, "forestplot_Mut.pdf"), 
    width = 4.5, height = 3.5, useDingbats = FALSE)

print(p)

dev.off()


jpeg(file=file.path(dir_clin_mut, "forestplot_Mut.jpeg"), 
    width = 4.5, height = 3.5, units = "in", res = 300, quality = 100)
print(p)
dev.off()
