##############################################################
## Script: CNS NGS gene-level co-mutation analysis
##
## Purpose:
##   Assess pairwise gene co-mutation patterns using
##   Fisher’s exact tests applied to binary mutation data
##   (0 = WT, 1 = Mut).
##
##   The analysis quantifies mutation co-occurrence or
##   mutual exclusivity between gene pairs.
##
## Input:
##   - result/data/se_mut_bin_clin.RData
##       SummarizedExperiment containing:
##         • assay: binary mutation matrix (genes × patients; 0/1)
##         • colData: clinical metadata including IDH.status
##
## Analyses:
##
##   1) All patients
##      - Filter genes with mutation frequency ≥ min_mut_freq
##      - Construct pairwise 2×2 contingency tables
##      - Apply Fisher’s exact test
##
##   2) IDH-stratified analysis
##      - Repeat analysis separately in:
##          • IDH-WT patients
##          • IDH-Mut patients
##
## Statistical Framework:
##
##   For each gene pair (g1, g2), a 2×2 table is constructed:
##
##              g2 Mut   g2 WT
##     g1 Mut      n11     n10
##     g1 WT       n01     n00
##
##   Fisher’s exact test is applied to evaluate independence.
##
##   Odds Ratio (OR):
##       OR = (n11 × n00) / (n10 × n01)
##
##       OR > 1  → co-occurrence
##       OR < 1  → mutual exclusivity
##
##   Effect size visualization:
##       log2(OR)
##
##   Multiple testing correction:
##       Benjamini–Hochberg FDR
##
## Output:
##
##   CSV files:
##     - coMut_res_all.csv
##     - coMut_res_wt.csv
##     - coMut_res_mut.csv
##
##     Each file includes:
##       gene1, gene2,
##       OR, pval,
##       n11, n10, n01, n00,
##       FDR
##
##   Figures (PDF):
##     - volcano_coMut_all.pdf
##     - volcano_coMut_wt.pdf
##     - volcano_coMut_mut.pdf
##
##     Volcano plots display:
##       x-axis: log2(OR)
##       y-axis: −log10(p-value)
##       Significant pairs highlighted (FDR < 0.05)
##############################################################
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
clin <- clin[clin$IDH.status != 'NA', ]
clin$IDH.status <- ifelse(clin$IDH.status == 'IDHwt', 'WT', 'Mut') 

##############################################################
## co-mutation analysis (pairwise Fisher exact tests)
##############################################################
## Step 1 --- All patients
# filter genes with sufficient mutation frequency
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
co_res <- co_res[!is.na(co_res$OR), ]
co_res <- co_res %>%
  group_by(gene1) %>%
  mutate(fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()

write.csv(co_res, file = file.path(dir_output, 'coMut_res_all.csv'), row.names=FALSE)

## Step 2 --- WT patients
# filter genes with sufficient mutation frequency
min_mut_freq <- 5
clin_wt <- clin[clin$IDH.status == 'WT', ]
mut_wt <- mut[ , colnames(mut) %in% clin_wt$Study]

mut_counts <- rowSums(mut_wt == 1)
genes_keep <- names(mut_counts[mut_counts >= min_mut_freq])

mut_filt <- mut_wt[genes_keep, ]
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
co_res <- co_res[!is.na(co_res$OR), ]
co_res <- co_res %>%
  group_by(gene1) %>%
  mutate(fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()

write.csv(co_res, file = file.path(dir_output, 'coMut_res_wt.csv'), row.names=FALSE)

## Step 2 --- MUT patients
# filter genes with sufficient mutation frequency
min_mut_freq <- 5
clin_mut <- clin[clin$IDH.status == 'Mut', ]
mut_mut <- mut[ , colnames(mut) %in% clin_mut$Study]

mut_counts <- rowSums(mut_mut == 1)
genes_keep <- names(mut_counts[mut_counts >= min_mut_freq])

mut_filt <- mut_wt[genes_keep, ]
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
co_res <- co_res[!is.na(co_res$OR), ]
co_res <- co_res %>%
  group_by(gene1) %>%
  mutate(fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()

write.csv(co_res, file = file.path(dir_output, 'coMut_res_mut.csv'), row.names=FALSE)

#############################################################
## Visualize
#############################################################
## Step 1 --- all patients
co_res <- read.csv(file.path(dir_output, 'coMut_res_all.csv'))

vol_idh1 <- co_res %>%
  #filter(gene1 == "IDH1") %>%
  mutate(
    pair = paste(gene1, gene2, sep = "–"),
    log2_or = log2(ifelse(OR == 0, 1e-6, OR)),
    neglog10 = -log10(pval),
    sig = ifelse(fdr < 0.05, "FDR < 0.05", "NS")
  )

pdf(file.path(dir_output, 'Fig', 'volcano_coMut_all.pdf'), width = 5, height = 4)

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


## Step 2 --- WT patients
co_res <- read.csv(file.path(dir_output, 'coMut_res_wt.csv'))

vol_idh1 <- co_res %>%
  #filter(gene1 == "IDH1") %>%
  mutate(
    pair = paste(gene1, gene2, sep = "–"),
    log2_or = log2(ifelse(OR == 0, 1e-6, OR)),
    neglog10 = -log10(pval),
    sig = ifelse(fdr < 0.05, "FDR < 0.05", "NS")
  )

pdf(file.path(dir_output, 'Fig', 'volcano_coMut_wt.pdf'), width = 5, height = 4)

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


## Step 3 --- Mut patients
co_res <- read.csv(file.path(dir_output, 'coMut_res_mut.csv'))

vol_idh1 <- co_res %>%
  #filter(gene1 == "IDH1") %>%
  mutate(
    pair = paste(gene1, gene2, sep = "–"),
    log2_or = log2(ifelse(OR == 0, 1e-6, OR)),
    neglog10 = -log10(pval),
    sig = ifelse(fdr < 0.05, "FDR < 0.05", "NS")
  )

pdf(file.path(dir_output, 'Fig', 'volcano_coMut_mut.pdf'), width = 5, height = 4)

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
