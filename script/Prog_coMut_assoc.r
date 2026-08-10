##-------------------------------------------------------------------
## Script: CNS NGS Gene-Level Co-Mutation Analysis
## Purpose:
##   To evaluate pairwise gene co-mutation and mutual
##   exclusivity patterns in CNS tumour patients using
##   curated binary mutation and clinical data.
##
##   Statistical associations between recurrent gene
##   alterations are evaluated in the full cohort and
##   IDH- and histology-stratified subgroups.
##
## Input:
##   - result/data/mae_mut_clin.RData
##       MultiAssayExperiment object containing:
##         • binary gene-level mutation matrix
##         • harmonized clinical metadata
##
## Outputs:
##   1) Pairwise co-mutation association tables
##   2) Volcano plot visualizations
##   3) UpSet plot visualizations
##
## Processing Overview:
##   1) Load mutation and clinical data
##   2) Combine IDH1 and IDH2 mutation calls into
##      a single IDH mutation indicator
##   3) Retain genes mutated in at least five
##      patients within each analysis cohort
##   4) Perform pairwise Fisher's exact tests and
##      estimate odds ratios
##   5) Adjust p-values for multiple testing using
##      the Benjamini-Hochberg FDR method
##   6) Repeat analyses within:
##        - Full cohort
##        - IDH wild-type subgroup
##        - IDH mutant subgroup
##        - IDH wild-type glioblastoma subgroup
##        - IDH wild-type non-glioblastoma subgroup
##   7) Generate volcano plots showing association
##      strength and statistical significance
##   8) Generate UpSet plots showing patient-level
##      mutation intersections among selected genes
##   9) Export statistical results and figures
##
## Notes:
##   - Analyses are based on binary gene-level
##     mutation matrices.
##   - IDH1 and IDH2 mutations are combined into
##     a single IDH mutation indicator.
##   - Fisher's exact test is used for pairwise
##     mutation association testing.
##   - Odds ratios greater than 1 indicate
##     co-occurrence, whereas odds ratios less than
##     1 indicate negative association or potential
##     mutual exclusivity.
##   - A Haldane-Anscombe correction of 0.5 is
##     applied when calculating odds ratios.
##   - Multiple testing correction is performed
##     using false discovery rate (FDR) adjustment.
##   - Volcano plots distinguish FDR-significant,
##     nominally significant, and non-significant
##     associations where applicable.
##   - UpSet gene selection is based primarily on
##     FDR-significant associations, with selected
##     nominal associations included for subgroups
##     with limited statistical power.
##-------------------------------------------------------------------
####################################################
## Load libraries
####################################################
library(MultiAssayExperiment)
library(ggplot2)
library(dplyr) 
library(tidyr)
library(stringr)
library(paletteer)
library(ggrepel)
library(UpSetR)
library(grid)

####################################################
## Setup directories
####################################################
dir_input <- 'result/data'
dir_output <- 'result/coMut'

#################################################
## Load data
#################################################
load(file.path(dir_input, 'mae_mut_clin.RData'))
mut <- assay(mae[['mut_binary']])
clin <- colData(mae[['mut_binary']])

## ---- Binary matrix
if(all(c("IDH1","IDH2") %in% rownames(mut))){

  idh_vec <- as.integer(
    mut["IDH1", ] == 1 | mut["IDH2", ] == 1
  )
  mut <- mut[!rownames(mut) %in% c("IDH1","IDH2"), , drop = FALSE]
  mut <- rbind(mut, IDH = idh_vec)

}

##############################################################
## co-mutation analysis (pairwise Fisher exact tests)
##############################################################
## Step 1 --- All patients
# filter genes with sufficient mutation frequency
min_mut_freq <- 5

mut_counts <- rowSums(mut == 1, na.rm = TRUE)
genes_keep <- names(mut_counts[mut_counts >= min_mut_freq])
mut_filt <- mut[genes_keep, ]

gene_pairs <- combn(rownames(mut_filt), 2)

co_res <- lapply(1:ncol(gene_pairs), function(i) {

  g1 <- gene_pairs[1, i]
  g2 <- gene_pairs[2, i]

  vec1 <- as.numeric(mut_filt[g1, ])
  vec2 <- as.numeric(mut_filt[g2, ])

  tab <- table(factor(vec1, levels = c(0,1)),
               factor(vec2, levels = c(0,1)))

  if(all(dim(tab) == c(2,2))) {
    
    n11 = tab["1","1"]
    n10 = tab["1","0"]
    n01 = tab["0","1"]
    n00 = tab["0","0"]

    fit <- fisher.test(tab)
    OR_corrected <- (n11 + 0.5)*(n00 + 0.5) / ((n10 + 0.5)*(n01 + 0.5))

    data.frame(
      gene1 = g1,
      gene2 = g2,
      OR_raw = as.numeric(fit$estimate),
      OR = OR_corrected,                 
      pval = fit$p.value,
      n11 = n11,
      n10 = n10,
      n01 = n01,
      n00 = n00
    )

  } else {
    
     data.frame(
      gene1 = g1,
      gene2 = g2,
      OR_raw = NA,
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
co_res <- co_res[!is.na(co_res$OR), ]
co_res$fdr <- p.adjust(co_res$pval, method = "BH")

write.csv(co_res, file = file.path(dir_output, 'coMut_res_all.csv'), row.names=FALSE)

## Step 2 --- WT patients
# filter genes with sufficient mutation frequency
min_mut_freq <- 5
clin_wt <- clin[clin$IDH_status == 'WT', ]
mut_wt <- mut[ , colnames(mut) %in% clin_wt$Study]

mut_counts <- rowSums(mut_wt == 1, na.rm = TRUE)
genes_keep <- names(mut_counts[mut_counts >= min_mut_freq])

mut_filt <- mut_wt[genes_keep, ]
gene_pairs <- combn(rownames(mut_filt), 2)

co_res <- lapply(1:ncol(gene_pairs), function(i) {
 
  g1 <- gene_pairs[1, i]
  g2 <- gene_pairs[2, i]

  vec1 <- as.numeric(mut_filt[g1, ])
  vec2 <- as.numeric(mut_filt[g2, ])

   tab <- table(factor(vec1, levels = c(0,1)),
               factor(vec2, levels = c(0,1)))

  if(all(dim(tab) == c(2,2))) {
    
    n11 = tab["1","1"]
    n10 = tab["1","0"]
    n01 = tab["0","1"]
    n00 = tab["0","0"]

    fit <- fisher.test(tab)
    OR_corrected <- (n11 + 0.5)*(n00 + 0.5) / ((n10 + 0.5)*(n01 + 0.5))

    data.frame(
      gene1 = g1,
      gene2 = g2,
      OR_raw = as.numeric(fit$estimate),
      OR = OR_corrected,                 
      pval = fit$p.value,
      n11 = n11,
      n10 = n10,
      n01 = n01,
      n00 = n00
    )

  } else {
    
     data.frame(
      gene1 = g1,
      gene2 = g2,
      OR_raw = NA,
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
co_res <- co_res[!is.na(co_res$OR), ]
co_res$fdr <- p.adjust(co_res$pval, method = "BH")

write.csv(co_res, file = file.path(dir_output, 'coMut_res_wt.csv'), row.names=FALSE)

## Step 3 --- MUT patients
# filter genes with sufficient mutation frequency
min_mut_freq <- 5
clin_mut <- clin[clin$IDH_status == 'Mut', ]
mut_mut <- mut[ , colnames(mut) %in% clin_mut$Study]

mut_counts <- rowSums(mut_mut == 1, na.rm = TRUE)
genes_keep <- names(mut_counts[mut_counts >= min_mut_freq])

mut_filt <- mut_mut[genes_keep, ]
gene_pairs <- combn(rownames(mut_filt), 2)

co_res <- lapply(1:ncol(gene_pairs), function(i) {
 
  g1 <- gene_pairs[1, i]
  g2 <- gene_pairs[2, i]

  vec1 <- as.numeric(mut_filt[g1, ])
  vec2 <- as.numeric(mut_filt[g2, ])

   tab <- table(factor(vec1, levels = c(0,1)),
               factor(vec2, levels = c(0,1)))

  if(all(dim(tab) == c(2,2))) {
    
    n11 = tab["1","1"]
    n10 = tab["1","0"]
    n01 = tab["0","1"]
    n00 = tab["0","0"]

    fit <- fisher.test(tab)
    OR_corrected <- (n11 + 0.5)*(n00 + 0.5) / ((n10 + 0.5)*(n01 + 0.5))

    data.frame(
      gene1 = g1,
      gene2 = g2,
      OR_raw = as.numeric(fit$estimate),
      OR = OR_corrected,                 
      pval = fit$p.value,
      n11 = n11,
      n10 = n10,
      n01 = n01,
      n00 = n00
    )

  } else {
    
     data.frame(
      gene1 = g1,
      gene2 = g2,
      OR_raw = NA,
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
co_res <- co_res[!is.na(co_res$OR), ]
co_res$fdr <- p.adjust(co_res$pval, method = "BH")

write.csv(co_res, file = file.path(dir_output, 'coMut_res_mut.csv'), row.names=FALSE)

## Step 4 --- WT GBM patients
# filter genes with sufficient mutation frequency
min_mut_freq <- 5
clin_wt_gbm <- clin[clin$IDH_status == 'WT' & clin$histo == 'Glioblastoma', ]
mut_wt_gbm <- mut[ , colnames(mut) %in% clin_wt_gbm$Study]

mut_counts <- rowSums(mut_wt_gbm == 1, na.rm = TRUE)
genes_keep <- names(mut_counts[mut_counts >= min_mut_freq])

mut_filt <- mut_wt_gbm[genes_keep, ]
gene_pairs <- combn(rownames(mut_filt), 2)

co_res <- lapply(1:ncol(gene_pairs), function(i) {
 
  g1 <- gene_pairs[1, i]
  g2 <- gene_pairs[2, i]

  vec1 <- as.numeric(mut_filt[g1, ])
  vec2 <- as.numeric(mut_filt[g2, ])

   tab <- table(factor(vec1, levels = c(0,1)),
               factor(vec2, levels = c(0,1)))

  if(all(dim(tab) == c(2,2))) {
    
    n11 = tab["1","1"]
    n10 = tab["1","0"]
    n01 = tab["0","1"]
    n00 = tab["0","0"]

    fit <- fisher.test(tab)
    OR_corrected <- (n11 + 0.5)*(n00 + 0.5) / ((n10 + 0.5)*(n01 + 0.5))

    data.frame(
      gene1 = g1,
      gene2 = g2,
      OR_raw = as.numeric(fit$estimate),
      OR = OR_corrected,                 
      pval = fit$p.value,
      n11 = n11,
      n10 = n10,
      n01 = n01,
      n00 = n00
    )

  } else {
    
     data.frame(
      gene1 = g1,
      gene2 = g2,
      OR_raw = NA,
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
co_res <- co_res[!is.na(co_res$OR), ]
co_res$fdr <- p.adjust(co_res$pval, method = "BH")

write.csv(co_res, file = file.path(dir_output, 'coMut_res_wt_gbm.csv'), row.names=FALSE)

## Step 5 --- WT nonGBM patients
# filter genes with sufficient mutation frequency
min_mut_freq <- 5
clin_wt_nongbm <- clin[clin$IDH_status == 'WT' & clin$histo != 'Glioblastoma', ]
mut_wt_nongbm <- mut[ , colnames(mut) %in% clin_wt_nongbm$Study]

mut_counts <- rowSums(mut_wt_nongbm == 1, na.rm = TRUE)
genes_keep <- names(mut_counts[mut_counts >= min_mut_freq])

mut_filt <- mut_wt_nongbm[genes_keep, ]
gene_pairs <- combn(rownames(mut_filt), 2)

co_res <- lapply(1:ncol(gene_pairs), function(i) {
 
  g1 <- gene_pairs[1, i]
  g2 <- gene_pairs[2, i]

  vec1 <- as.numeric(mut_filt[g1, ])
  vec2 <- as.numeric(mut_filt[g2, ])

   tab <- table(factor(vec1, levels = c(0,1)),
               factor(vec2, levels = c(0,1)))

  if(all(dim(tab) == c(2,2))) {
    
    n11 = tab["1","1"]
    n10 = tab["1","0"]
    n01 = tab["0","1"]
    n00 = tab["0","0"]

    fit <- fisher.test(tab)
    OR_corrected <- (n11 + 0.5)*(n00 + 0.5) / ((n10 + 0.5)*(n01 + 0.5))

    data.frame(
      gene1 = g1,
      gene2 = g2,
      OR_raw = as.numeric(fit$estimate),
      OR = OR_corrected,                 
      pval = fit$p.value,
      n11 = n11,
      n10 = n10,
      n01 = n01,
      n00 = n00
    )

  } else {
    
     data.frame(
      gene1 = g1,
      gene2 = g2,
      OR_raw = NA,
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
co_res <- co_res[!is.na(co_res$OR), ]
co_res$fdr <- p.adjust(co_res$pval, method = "BH")

write.csv(co_res, file = file.path(dir_output, 'coMut_res_wt_nongbm.csv'), row.names=FALSE)

#############################################################
## Visualize ---- volcano plot
#############################################################
## Step 1 --- all patients
co_res <- read.csv(file.path(dir_output, 'coMut_res_all.csv'))

vol_all <- co_res %>%
  mutate(
    pair = paste(gene1, gene2, sep = "–"),

    log2_or = log2(
      ifelse(is.na(OR) | OR <= 0, 1e-6, OR)
    ),

    neglog10 = -log10(
      ifelse(is.na(pval) | pval <= 0, 1e-300, pval)
    ),

    sig = case_when(
      fdr < 0.05 ~ "FDR < 0.05",
      pval < 0.05 ~ "p < 0.05",
      TRUE       ~ "NS"
    ),

    sig = factor(
      sig,
      levels = c(
        "FDR < 0.05",
        "p < 0.05",
        "NS"
      )
    )
  )

## Check the number of pairs in each group
table(vol_all$sig, useNA = "ifany")

p <- ggplot(
  vol_all,
  aes(x = log2_or, y = neglog10)
) +
  geom_point(
    aes(color = sig),
    size = 1.7,
    alpha = 0.9
  ) +

  ## Label pairs with FDR < 0.05
  geom_text_repel(
    data = vol_all %>%
      filter(fdr < 0.05),
    aes(label = pair),
    size = 1.7,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0
  ) +

  scale_color_manual(
    values = c(
      "FDR < 0.05"      = "#515260FF",
      "p < 0.05"        = "#A56A3EFF",  
      "NS"               = "#A5A6AEFF"
    ),
    drop = FALSE
  ) +

  theme_classic() +

  labs(
    title = "",
    x = "log2(Odds Ratio)",
    y = "-log10(p-value)",
    color = ""
  ) +

  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7)
  )

pdf(file.path(dir_output,  'volcano_coMut_all.pdf'), width = 4, height = 3)
print(p)
dev.off()

jpeg(file.path(dir_output,  'volcano_coMut_all.jpeg'), width = 4, height = 3,
     units = "in", res = 300, quality = 100)
print(p)
dev.off()

## Step 2 --- WT patients
co_res <- read.csv(file.path(dir_output, 'coMut_res_wt.csv'))

vol_idh1 <- co_res %>%
  mutate(
    pair = paste(gene1, gene2, sep = "–"),

    log2_or = log2(
      ifelse(is.na(OR) | OR <= 0, 1e-6, OR)
    ),

    neglog10 = -log10(
      ifelse(is.na(pval) | pval <= 0, 1e-300, pval)
    ),

    sig = case_when(
      fdr < 0.05 ~ "FDR < 0.05",
      pval < 0.05 ~ "p < 0.05",
      TRUE       ~ "NS"
    ),

    sig = factor(
      sig,
      levels = c(
        "FDR < 0.05",
        "p < 0.05",
        "NS"
      )
    )
  )

## Check the number of pairs in each group
table(vol_idh1$sig, useNA = "ifany")

p <- ggplot(
  vol_idh1,
  aes(x = log2_or, y = neglog10)
) +
  geom_point(
    aes(color = sig),
    size = 1.7,
    alpha = 0.9
  ) +

  ## Label pairs with FDR < 0.05
  geom_text_repel(
    data = vol_idh1 %>%
      filter(fdr < 0.05),
    aes(label = pair),
    size = 1.7,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0
  ) +

  scale_color_manual(
    values = c(
      "FDR < 0.05"      = "#515260FF",
      "p < 0.05"        = "#A56A3EFF",  
      "NS"               = "#A5A6AEFF"
    ),
    drop = FALSE
  ) +

  theme_classic() +

  labs(
    title = "",
    x = "log2(Odds Ratio)",
    y = "-log10(p-value)",
    color = ""
  ) +

  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7)
  )

pdf(file.path(dir_output,  'volcano_coMut_wt.pdf'), width = 3.5, height = 2.5)
print(p)
dev.off()

jpeg(file.path(dir_output,  'volcano_coMut_wt.jpeg'), width = 3.5, height = 2.5,
     units = "in", res = 300, quality = 100)
print(p)
dev.off()

## Step 3 --- Mut patients
co_res <- read.csv(file.path(dir_output, 'coMut_res_mut.csv'))

vol_idh1 <- co_res %>%
  mutate(
    pair = paste(gene1, gene2, sep = "–"),

    log2_or = log2(
      ifelse(is.na(OR) | OR <= 0, 1e-6, OR)
    ),

    neglog10 = -log10(
      ifelse(is.na(pval) | pval <= 0, 1e-300, pval)
    ),

    sig = case_when(
      fdr < 0.05 ~ "FDR < 0.05",
      pval < 0.05 ~ "p < 0.05",
      TRUE       ~ "NS"
    ),

    sig = factor(
      sig,
      levels = c(
        "FDR < 0.05",
        "p < 0.05",
        "NS"
      )
    )
  )

## Check the number of pairs in each group
table(vol_idh1$sig, useNA = "ifany")

p <- ggplot(
  vol_idh1,
  aes(x = log2_or, y = neglog10)
) +
  geom_point(
    aes(color = sig),
    size = 1.7,
    alpha = 0.9
  ) +

  ## Label pairs with pval < 0.05
  geom_text_repel(
    data = vol_idh1 %>%
      filter(pval < 0.05),
    aes(label = pair),
    size = 1.7,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0
  ) +

   scale_color_manual(
    values = c(
      "FDR < 0.05"      = "#515260FF",
      "p < 0.05"        = "#A56A3EFF",  
      "NS"               = "#A5A6AEFF"
    ),
    drop = FALSE
  ) +

  theme_classic() +

  labs(
    title = "",
    x = "log2(Odds Ratio)",
    y = "-log10(p-value)",
    color = ""
  ) +

  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7)
  )

pdf(file.path(dir_output,  'volcano_coMut_mut.pdf'), width = 3.5, height = 2.5)
print(p)
dev.off()

jpeg(file.path(dir_output,  'volcano_coMut_mut.jpeg'), width = 3.5, height = 2.5,
     units = "in", res = 300, quality = 100)
print(p)
dev.off()

## Step 4 --- WT GBM patients
co_res <- read.csv(file.path(dir_output, 'coMut_res_wt_gbm.csv'))

vol_idh1 <- co_res %>%
  mutate(
    pair = paste(gene1, gene2, sep = "–"),

    log2_or = log2(
      ifelse(is.na(OR) | OR <= 0, 1e-6, OR)
    ),

    neglog10 = -log10(
      ifelse(is.na(pval) | pval <= 0, 1e-300, pval)
    ),

    sig = case_when(
      fdr < 0.05 ~ "FDR < 0.05",
      pval < 0.05 ~ "p < 0.05",
      TRUE       ~ "NS"
    ),

    sig = factor(
      sig,
      levels = c(
        "FDR < 0.05",
        "p < 0.05",
        "NS"
      )
    )
  )

## Check the number of pairs in each group
table(vol_idh1$sig, useNA = "ifany")

## Select labels
label_pairs <- vol_idh1 %>%
  filter(fdr < 0.05) %>%
  bind_rows(
    vol_idh1 %>%
      filter(fdr >= 0.05, pval < 0.05) %>%
      arrange(pval) %>%
      slice_head(n = 6)
  ) %>%
  distinct(pair, .keep_all = TRUE)


p <- ggplot(
  vol_idh1,
  aes(x = log2_or, y = neglog10)
) +
  geom_point(
    aes(color = sig),
    size = 1.7,
    alpha = 0.9
  ) +

  ## Label pairs with pval < 0.05
  geom_text_repel(
  data = label_pairs,
  aes(label = pair),
  size = 1.7,
  max.overlaps = Inf,
  box.padding = 0.25,
  point.padding = 0.15,
  min.segment.length = 0
) +

   scale_color_manual(
    values = c(
      "FDR < 0.05"      = "#515260FF",
      "p < 0.05"        = "#A56A3EFF",  
      "NS"               = "#A5A6AEFF"
    ),
    drop = FALSE
  ) +

  theme_classic() +

  labs(
    title = "",
    x = "log2(Odds Ratio)",
    y = "-log10(p-value)",
    color = ""
  ) +

  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7)
  )

pdf(file.path(dir_output,  'volcano_coMut_wt_gbm.pdf'), width = 3.5, height = 2.5)
print(p)
dev.off()

jpeg(file.path(dir_output,  'volcano_coMut_wt_gbm.jpeg'), width = 3.5, height = 2.5,
     units = "in", res = 300, quality = 100)
print(p)
dev.off()

## Step 5 --- WT nonGBM patients
co_res <- read.csv(file.path(dir_output, 'coMut_res_wt_nongbm.csv'))

vol_idh1 <- co_res %>%
  mutate(
    pair = paste(gene1, gene2, sep = "–"),

    log2_or = log2(
      ifelse(is.na(OR) | OR <= 0, 1e-6, OR)
    ),

    neglog10 = -log10(
      ifelse(is.na(pval) | pval <= 0, 1e-300, pval)
    ),

    sig = case_when(
     # fdr < 0.05 ~ "FDR < 0.05",
      pval < 0.05 ~ "p < 0.05",
      TRUE       ~ "NS"
    ),

    sig = factor(
      sig,
      levels = c(
     #   "FDR < 0.05",
        "p < 0.05",
        "NS"
      )
    )
  )

## Check the number of pairs in each group
table(vol_idh1$sig, useNA = "ifany")

p <- ggplot(
  vol_idh1,
  aes(x = log2_or, y = neglog10)
) +
  geom_point(
    aes(color = sig),
    size = 1.7,
    alpha = 0.9
  ) +

  ## Label pairs with pval < 0.05
  geom_text_repel(
    data = vol_idh1 %>%
      filter(pval < 0.05),
    aes(label = pair),
    size = 1.7,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0
  ) +

   scale_color_manual(
    values = c(
    #  "FDR < 0.05"      = "#515260FF",
      "p < 0.05"        = "#A56A3EFF",  
      "NS"               = "#A5A6AEFF"
    ),
    drop = FALSE
  ) +

  theme_classic() +

  labs(
    title = "",
    x = "log2(Odds Ratio)",
    y = "-log10(p-value)",
    color = ""
  ) +

  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7)
  )

pdf(file.path(dir_output,  'volcano_coMut_wt_nongbm.pdf'), width = 3.5, height = 2.5)
print(p)
dev.off()

jpeg(file.path(dir_output,  'volcano_coMut_wt_nongbm.jpeg'), width = 3.5, height = 2.5,
     units = "in", res = 300, quality = 100)
print(p)
dev.off()

#############################################################
## Visualize ---- UpSet plot
#############################################################
gene_colors <- c(
  "TERTp"  = '#855C75FF',
  "IDH"    = '#99B6BDFF',
  "TP53"   = '#736F4CFF',
  "EGFR"   = '#4A7169FF',
  "CDK4"   = '#BF816BFF',
  "CDKN2A" = '#E3CA97FF',
  "CDKN2B" = '#7A8FA6FF',
  "PIK3CA" = '#B49696FF',
  "FGFR1"  = "#788787FF",
  "ATRX"   = "#E6BCACFF",
  "MDM4"   = '#9A7FA5FF',
  "MDM2"   = '#B5A56AFF',
  "PDGFRA" = '#6F8F8AFF'
)

## Step 1 --- all patients
co_res <- read.csv(file.path(dir_output, 'coMut_res_all.csv'))
sig_pairs <- co_res %>%
  filter(fdr < 0.05) %>%
  arrange(fdr)
selected_genes <- unique(c(sig_pairs$gene1, sig_pairs$gene2))
selected_genes <- intersect(selected_genes, rownames(mut))
gene_freq <- rowSums(mut, na.rm = TRUE)

selected_genes <- selected_genes[
  gene_freq[selected_genes] >= 5
]

gene_degree <- table(c(sig_pairs$gene1, sig_pairs$gene2))
gene_degree <- sort(gene_degree, decreasing = TRUE)
selected_genes <- intersect(names(gene_degree), selected_genes)

top_genes <- head(selected_genes, 8)
top_genes  # check what you selected

mut_top <- mut[top_genes, , drop = FALSE]
listInput <- lapply(rownames(mut_top), function(g) {
  colnames(mut_top)[mut_top[g, ] == 1]
})

names(listInput) <- rownames(mut_top)
col.id <- gene_colors[top_genes]
sets <- rownames(mut_top)
sets.bar.color <- col.id
mainbar.y.label <- "Number of Patients"
sets.x.label <- "Mutation Frequency"

pdf(file.path(dir_output,  'upset_coMut_all.pdf'), width = 7, height = 5)

suppressWarnings(
  upset.fun(fromList(listInput), 'Number of Patients')
)

dev.off()

jpeg(file.path(dir_output,  'upset_coMut_all.jpeg'), width = 7, height = 5,
     units = "in", res = 300, quality = 100)

suppressWarnings(
  upset.fun(fromList(listInput), 'Number of Patients')
)

dev.off()

## Step 2 --- WT patients
co_res <- read.csv(file.path(dir_output, 'coMut_res_wt.csv'))

sig_pairs <- co_res %>%
  filter(fdr < 0.05) %>%
  arrange(fdr)
selected_genes <- unique(c(sig_pairs$gene1, sig_pairs$gene2))
selected_genes <- intersect(selected_genes, rownames(mut))
gene_freq <- rowSums(mut_wt, na.rm = TRUE)

selected_genes <- selected_genes[
  gene_freq[selected_genes] >= 5
]

gene_degree <- table(c(sig_pairs$gene1, sig_pairs$gene2))
gene_degree <- sort(gene_degree, decreasing = TRUE)

selected_genes <- intersect(names(gene_degree), selected_genes)

top_genes <- head(selected_genes, 7)
top_genes  # check what you selected

mut_top <- mut_wt[top_genes, , drop = FALSE]
listInput <- lapply(rownames(mut_top), function(g) {
  colnames(mut_top)[mut_top[g, ] == 1]
})

names(listInput) <- rownames(mut_top)
col.id <- gene_colors[top_genes]
sets <- rownames(mut_top)
sets.bar.color <- col.id

pdf(file.path(dir_output,  'upset_coMut_wt.pdf'), width = 7, height = 5)

suppressWarnings(
  upset.fun(fromList(listInput), 'Number of Patients')
)

dev.off()

jpeg(file.path(dir_output,  'upset_coMut_wt.jpeg'), width = 7, height = 5,
     units = "in", res = 300, quality = 100)

suppressWarnings(
  upset.fun(fromList(listInput), 'Number of Patients')
)

dev.off()

## Step 3 --- Mut patients
co_res <- read.csv(file.path(dir_output, 'coMut_res_mut.csv'))

sig_pairs <- co_res %>%
  filter(pval < 0.05) %>%
  arrange(pval)
selected_genes <- unique(c(sig_pairs$gene1, sig_pairs$gene2))
selected_genes <- intersect(selected_genes, rownames(mut))
gene_freq <- rowSums(mut_mut, na.rm = TRUE)

selected_genes <- selected_genes[
  gene_freq[selected_genes] >= 5
]

gene_degree <- table(c(sig_pairs$gene1, sig_pairs$gene2))
gene_degree <- sort(gene_degree, decreasing = TRUE)

selected_genes <- intersect(names(gene_degree), selected_genes)

top_genes <- head(selected_genes, 7)
top_genes  # check what you selected

mut_top <- mut_mut[top_genes, , drop = FALSE]
listInput <- lapply(rownames(mut_top), function(g) {
  colnames(mut_top)[mut_top[g, ] == 1]
})

names(listInput) <- rownames(mut_top)
col.id <- gene_colors[top_genes]
sets <- rownames(mut_top)
sets.bar.color <- col.id

pdf(file.path(dir_output,  'upset_coMut_mut.pdf'), width = 7, height = 5)

suppressWarnings(
  upset.fun(fromList(listInput), 'Number of Patients')
)

dev.off()

jpeg(file.path(dir_output,  'upset_coMut_mut.jpeg'), width = 7, height = 5,
     units = "in", res = 300, quality = 100)

suppressWarnings(
  upset.fun(fromList(listInput), 'Number of Patients')
)

dev.off()

## Step 4 --- WT GBM patients
co_res <- read.csv(file.path(dir_output, 'coMut_res_wt_gbm.csv'))

sig_pairs <- co_res %>%
  filter(fdr < 0.05) %>%
  bind_rows(
    co_res %>%
      filter(fdr >= 0.05, pval < 0.05) %>%
      arrange(pval) %>%
      slice_head(n = 3)
  ) %>%
  arrange(fdr, pval)
selected_genes <- unique(c(sig_pairs$gene1, sig_pairs$gene2))
selected_genes <- intersect(selected_genes, rownames(mut))
gene_freq <- rowSums(mut_wt_gbm, na.rm = TRUE)

selected_genes <- selected_genes[
  gene_freq[selected_genes] >= 5
]

gene_degree <- table(c(sig_pairs$gene1, sig_pairs$gene2))
gene_degree <- sort(gene_degree, decreasing = TRUE)

selected_genes <- intersect(names(gene_degree), selected_genes)

top_genes <- head(selected_genes, 8)
top_genes  # check what you selected

mut_top <- mut_wt_gbm[top_genes, , drop = FALSE]
listInput <- lapply(rownames(mut_top), function(g) {
  colnames(mut_top)[mut_top[g, ] == 1]
})

names(listInput) <- rownames(mut_top)
col.id <- gene_colors[top_genes]
sets <- rownames(mut_top)
sets.bar.color <- col.id

pdf(file.path(dir_output,  'upset_coMut_wt_gbm.pdf'), width = 7, height = 5)

suppressWarnings(
  upset.fun(fromList(listInput), 'Number of Patients')
)

dev.off()

jpeg(file.path(dir_output,  'upset_coMut_wt_gbm.jpeg'), width = 7, height = 5,
     units = "in", res = 300, quality = 100)

suppressWarnings(
  upset.fun(fromList(listInput), 'Number of Patients')
)

dev.off()

## Step 5 --- WT nonGBM patients
co_res <- read.csv(file.path(dir_output, 'coMut_res_wt_nongbm.csv'))

sig_pairs <- co_res %>% 
           filter(pval < 0.05) %>% 
           arrange(pval)
selected_genes <- unique(c(sig_pairs$gene1, sig_pairs$gene2))
selected_genes <- intersect(selected_genes, rownames(mut))
gene_freq <- rowSums(mut_wt_nongbm, na.rm = TRUE)

selected_genes <- selected_genes[
  gene_freq[selected_genes] >= 5
]

gene_degree <- table(c(sig_pairs$gene1, sig_pairs$gene2))
gene_degree <- sort(gene_degree, decreasing = TRUE)

selected_genes <- intersect(names(gene_degree), selected_genes)

top_genes <- head(selected_genes, 8)
top_genes  # check what you selected

mut_top <- mut_wt_nongbm[top_genes, , drop = FALSE]
listInput <- lapply(rownames(mut_top), function(g) {
  colnames(mut_top)[mut_top[g, ] == 1]
})

names(listInput) <- rownames(mut_top)
col.id <- gene_colors[top_genes]
sets <- rownames(mut_top)
sets.bar.color <- col.id

pdf(file.path(dir_output,  'upset_coMut_wt_nongbm.pdf'), width = 7, height = 5)

suppressWarnings(
  upset.fun(fromList(listInput), 'Number of Patients')
)

dev.off()

jpeg(file.path(dir_output,  'upset_coMut_wt_nongbm.jpeg'), width = 7, height = 5,
     units = "in", res = 300, quality = 100)

suppressWarnings(
  upset.fun(fromList(listInput), 'Number of Patients')
)

dev.off()
