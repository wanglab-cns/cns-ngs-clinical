##############################################################
## Script: CNS NGS Gene-Level Co-Mutation Analysis
##
## Purpose:
##   Quantify pairwise gene co-mutation and mutual exclusivity
##   patterns using Fisher’s exact test applied to binary
##   mutation data (0 = wild-type, 1 = mutant).
##
##   The analysis evaluates statistical dependence between
##   gene mutation events across patients.
##
##
## Input:
##   - result/data/se_mut_bin_clin.RData
##       SummarizedExperiment containing:
##         • assay: binary gene × patient mutation matrix (0/1)
##         • colData: harmonized clinical metadata
##             including IDH_status
##
##
## Analytical Overview:
##
##   1) Gene Filtering
##      - Genes with mutation frequency ≥ min_mut_freq
##        are retained to reduce instability from rare events.
##
##   2) Pairwise Testing
##      - All possible gene pairs are enumerated.
##      - For each pair, a 2×2 contingency table is constructed:
##
##              g2 Mut   g2 WT
##     g1 Mut      n11     n10
##     g1 WT       n01     n00
##
##      - Fisher’s exact test evaluates independence.
##
##   3) Effect Size Interpretation
##      - Odds Ratio (OR):
##            OR = (n11 × n00) / (n10 × n01)
##
##            OR > 1  → co-occurrence enrichment
##            OR < 1  → mutual exclusivity
##
##      - Effect size visualized as log2(OR).
##
##   4) Multiple Testing Correction
##      - Benjamini–Hochberg procedure applied
##        across all tested gene pairs.
##
##
## Stratified Analyses:
##
##   Analyses are repeated in:
##      • All patients
##      • IDH wild-type subgroup
##      • IDH mutant subgroup
##
##
## Output:
##
##   CSV files:
##     - coMut_res_all.csv
##     - coMut_res_wt.csv
##     - coMut_res_mut.csv
##
##     Each file contains:
##       gene1, gene2,
##       OR, pval,
##       n11, n10, n01, n00,
##       fdr
##
##
##   Figures (PDF):
##
##     Volcano plots:
##       - volcano_coMut_all.pdf
##       - volcano_coMut_wt.pdf
##       - volcano_coMut_mut.pdf
##
##         x-axis: log2(OR)
##         y-axis: −log10(p-value)
##         Significant pairs highlighted (FDR < 0.05)
##
##     UpSet plots:
##       - upset_coMut_all.pdf
##       - upset_coMut_wt.pdf
##       - upset_coMut_mut.pdf
##
##         UpSet plots display mutation combination patterns
##         across the top co-mutation network genes
##         (ranked by significant interaction degree).
##
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
load(file.path(dir_input, 'se_mut_bin_clin.RData'))
mut <- assay(eset)
clin <- as.data.frame(colData(eset))

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
clin_wt <- clin[clin$IDH_status == 'WT', ]
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
clin_mut <- clin[clin$IDH_status == 'Mut', ]
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
## Visualize ---- volcano plot
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

pdf(file.path(dir_output,  'volcano_coMut_all.pdf'), width = 5, height = 4)

ggplot(vol_idh1, aes(x = log2_or, y = neglog10)) +
  geom_point(aes(color = sig), size = 2) +
  geom_text_repel(
    data = subset(vol_idh1, fdr < 0.05),
    aes(label = pair),
    size = 2.5
  ) +
  scale_color_manual(
    values = c("FDR < 0.05" = "#1d587a",
               "NS" = "grey70")
  ) +
  theme_classic() +
  labs(
    title = " ",
    x = "log2(Odds Ratio)",
    y = "-log10(p-value)",
    color = ""
  ) +
  theme(
    axis.title = element_text(size = 10),
    axis.text  = element_text(size = 8),
   # plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
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

pdf(file.path(dir_output,  'volcano_coMut_wt.pdf'), width = 5, height = 4)

ggplot(vol_idh1, aes(x = log2_or, y = neglog10)) +
  geom_point(aes(color = sig), size = 2) +
  geom_text_repel(
    data = subset(vol_idh1, fdr < 0.05),
    aes(label = pair),
    size = 2.5
  ) +
  scale_color_manual(
    values = c("FDR < 0.05" = "#1d587a",
               "NS" = "grey70")
  ) +
  theme_classic() +
  labs(
    title = " ",
    x = "log2(Odds Ratio)",
    y = "-log10(p-value)",
    color = ""
  ) +
  theme(
    axis.title = element_text(size = 10),
    axis.text  = element_text(size = 8),
   # plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
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

pdf(file.path(dir_output,  'volcano_coMut_mut.pdf'), width = 5, height = 4)

ggplot(vol_idh1, aes(x = log2_or, y = neglog10)) +
  geom_point(aes(color = sig), size = 2) +
  geom_text_repel(
    data = subset(vol_idh1, fdr < 0.05),
    aes(label = pair),
    size = 2.5
  ) +
  scale_color_manual(
    values = c("FDR < 0.05" = "#1d587a",
               "NS" = "grey70")
  ) +
  theme_classic() +
  labs(
    title = " ",
    x = "log2(Odds Ratio)",
    y = "-log10(p-value)",
    color = ""
  ) +
  theme(
    axis.title = element_text(size = 10),
    axis.text  = element_text(size = 8),
   # plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  )

dev.off()

#############################################################
## Visualize ---- UpSet plot
#############################################################
## Step 1 --- all patients
co_res <- read.csv(file.path(dir_output, 'coMut_res_all.csv'))
sig_pairs <- subset(co_res, pval < 0.05)
selected_genes <- unique(c(sig_pairs$gene1, sig_pairs$gene2))
selected_genes <- intersect(selected_genes, rownames(mut))
gene_freq <- rowSums(mut)

selected_genes <- selected_genes[
  gene_freq[selected_genes] >= 5
]

gene_degree <- table(c(sig_pairs$gene1, sig_pairs$gene2))
gene_degree <- sort(gene_degree, decreasing = TRUE)

selected_genes <- intersect(names(gene_degree), selected_genes)

top5_genes <- head(selected_genes, 5)
top5_genes  # check what you selected

mut_top <- mut[top5_genes, , drop = FALSE]
listInput <- lapply(rownames(mut_top), function(g) {
  colnames(mut_top)[mut_top[g, ] == 1]
})

names(listInput) <- rownames(mut_top)
col.id <- c('#855C75FF', '#736F4CFF', "#99B6BDFF", '#BF816BFF', '#4A7169FF')
sets <- rownames(mut_top)
sets.bar.color <- col.id

pdf(file.path(dir_output,  'upset_coMut_all.pdf'), width = 8, height = 6)

suppressWarnings(
  upset.fun(fromList(listInput))
)

dev.off()


## Step 2 --- WT patients
co_res <- read.csv(file.path(dir_output, 'coMut_res_wt.csv'))

sig_pairs <- subset(co_res, pval < 0.05)
selected_genes <- unique(c(sig_pairs$gene1, sig_pairs$gene2))
selected_genes <- intersect(selected_genes, rownames(mut))
gene_freq <- rowSums(mut)

selected_genes <- selected_genes[
  gene_freq[selected_genes] >= 5
]

gene_degree <- table(c(sig_pairs$gene1, sig_pairs$gene2))
gene_degree <- sort(gene_degree, decreasing = TRUE)

selected_genes <- intersect(names(gene_degree), selected_genes)

top5_genes <- head(selected_genes, 5)
top5_genes  # check what you selected

mut_top <- mut[top5_genes, , drop = FALSE]
listInput <- lapply(rownames(mut_top), function(g) {
  colnames(mut_top)[mut_top[g, ] == 1]
})

names(listInput) <- rownames(mut_top)
col.id <- c('#855C75FF', '#736F4CFF', "#99B6BDFF", '#BF816BFF', '#4A7169FF')
sets <- rownames(mut_top)
sets.bar.color <- col.id

pdf(file.path(dir_output,  'upset_coMut_wt.pdf'), width = 8, height = 6)

suppressWarnings(
  upset.fun(fromList(listInput))
)

dev.off()

## Step 3 --- Mut patients
co_res <- read.csv(file.path(dir_output, 'coMut_res_mut.csv'))

sig_pairs <- subset(co_res, pval < 0.05)
selected_genes <- unique(c(sig_pairs$gene1, sig_pairs$gene2))
selected_genes <- intersect(selected_genes, rownames(mut))
gene_freq <- rowSums(mut)

selected_genes <- selected_genes[
  gene_freq[selected_genes] >= 5
]

gene_degree <- table(c(sig_pairs$gene1, sig_pairs$gene2))
gene_degree <- sort(gene_degree, decreasing = TRUE)

selected_genes <- intersect(names(gene_degree), selected_genes)

top5_genes <- head(selected_genes, 3)
top5_genes  # check what you selected

mut_top <- mut[top5_genes, , drop = FALSE]
listInput <- lapply(rownames(mut_top), function(g) {
  colnames(mut_top)[mut_top[g, ] == 1]
})

names(listInput) <- rownames(mut_top)
col.id <- c('#855C75FF', '#736F4CFF', "#99B6BDFF")
sets <- rownames(mut_top)
sets.bar.color <- col.id

pdf(file.path(dir_output,  'upset_coMut_mut.pdf'), width = 8, height = 6)

suppressWarnings(
  upset.fun(fromList(listInput))
)

dev.off()
