#####################################################
## Script: CNS NGS gene-level mutation – OS association
##
## Purpose:
##   Evaluate associations between binary gene-level mutation
##   status (0 = WT, 1 = Mut) and overall survival (OS) using
##   Cox proportional hazards models.
##
##   Analyses include:
##     1) Univariable Cox models (all patients)
##     2) Multivariable Cox models adjusted for clinical covariates
##     3) IDH-stratified Cox models (IDH Mut and IDH WT separately)
##
##   Overall survival is administratively censored at a fixed
##   time horizon (time.censor, in months).
##
## Input:
##   - result/data/se_mut_bin_clin.RData
##       SummarizedExperiment object containing:
##         • assay: binary mutation matrix (genes x patients; 0/1)
##         • colData: clinical metadata including:
##             - os.time   (OS time in months)
##             - os.event  (1 = death/event, 0 = censored)
##             - IDH status
##             - Age
##             - Sex
##             - WHO grade
##             - Tumor location
##             - Histology
##
## Output:
##   CSV files with per-gene Cox results:
##
##     - cox_os_all.csv
##         Univariable Cox (all patients)
##
##     - cox_os_all_mv.csv
##         Multivariable Cox adjusted for:
##         IDH + Age + Sex + Grade + Location + Histology
##
##     - cox_os_IDH_mut.csv
##         Univariable Cox within IDH Mut patients
##
##     - cox_os_IDH_wt.csv
##         Univariable Cox within IDH WT patients
##
##   Each output includes:
##     gene, log(HR), SE, N, 95% CI (lower/upper),
##     Wald p-value, and Benjamini–Hochberg FDR.
##
## Key Processing Steps:
##   1) Extract binary mutation matrix and harmonized clinical data.
##   2) Recode and collapse clinical variables where necessary.
##   3) Apply administrative censoring at `time.censor` months.
##   4) For each gene:
##        - Require minimum counts in mutated (n1.cutoff)
##          and wild-type (n0.cutoff) groups.
##        - Fit Cox PH model.
##   5) Adjust p-values across genes using BH FDR correction.
##
## Notes:
##   - Binary mutation coding: 0 = WT, 1 = Mut.
##   - Multivariable models may be unstable when mutation
##     status is strongly correlated with IDH status or when
##     event counts are low.
##   - Interpretation should consider event-per-variable
##     constraints and potential separation.
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
library(ggrepel)

####################################################
## Setup directories
####################################################
dir_input <- 'result/data'
dir_output <- 'result/assoc'

#################################################
## Load data
#################################################
load(file.path(dir_input, 'se_mut_bin_clin.RData'))
mut <- assay(eset)
clin <- as.data.frame(colData(eset))

##############################################################
## co-mutation analysis (pairwise Fisher exact tests)
##############################################################
## filter genes with sufficient mutation frequency
min_mut_freq <- 5

mut_counts <- rowSums(mut == 1)
genes_keep <- names(mut_counts[mut_counts >= min_mut_freq])

mut_filt <- mut[genes_keep, ]
gene_pairs <- combn(rownames(mut_filt), 2)


co_res <- lapply(1:ncol(gene_pairs), function(i) {

  g1 <- gene_pairs[1, i]
  g2 <- gene_pairs[2, i]

  vec1 <- as.numeric(mut_filt[g1, ])
  vec2 <- as.numeric(mut_filt[g2, ])

  tab <- table(vec1, vec2)

  if(all(dim(tab) == c(2,2))) {

    fit <- fisher.test(tab)

    data.frame(
      gene1 = g1,
      gene2 = g2,
      OR = as.numeric(fit$estimate),
      pval = fit$p.value,
      n11 = tab["1","1"],
      n10 = tab["1","0"],
      n01 = tab["0","1"],
      n00 = tab["0","0"]
    )

  } else {
    
     data.frame(
      gene1 = g1,
      gene2 = g2,
      OR = NA,
      pval = NA,
      n11 = NA,
      n10 =NA,
      n01 = NA,
      n00 = NA
    )

  }
})

co_res <- do.call(rbind, co_res)
co_res <- co_res %>%
  group_by(gene1) %>%
  mutate(fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()

write.csv(co_res, file = file.path(dir_output, 'coMut_res_all.csv'))

#############################################################
## Visualize
#############################################################

vol_idh1 <- co_res %>%
 # filter(gene1 == "IDH1") %>%
  mutate(
    pair = paste(gene1, gene2, sep = "–"),
    log2_or = log2(ifelse(OR == 0, 1e-6, OR)),
    neglog10 = -log10(pval),
    sig = ifelse(fdr < 0.05, "FDR < 0.05", "NS")
  )


pdf(file.path(dir_output, 'volcano_coMut.pdf'), width = 5, height = 4)

ggplot(vol_idh1, aes(x = log2_or, y = neglog10)) +
  geom_point(aes(color = sig), size = 2) +
  geom_text_repel(
    data = subset(vol_idh1, fdr < 0.05),
    aes(label = pair),
    size = 2.5
  ) +
  #geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    values = c("FDR < 0.05" = "#1d587a",   # light blue
               "NS" = "grey70")
  ) +
  theme_classic() +
  labs(
    title = " ",
    x = "log2(Odds Ratio)",
    y = "-log10(p-value)",
    color = ""
  )
dev.off()
