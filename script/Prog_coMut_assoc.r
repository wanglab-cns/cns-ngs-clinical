##-------------------------------------------------------------------
## Script: CNS NGS Gene-Level Co-Mutation Analysis
##
## Purpose:
##   Quantify pairwise gene co-mutation and mutual exclusivity
##   patterns in CNS tumor patients using binary mutation data.
##
##   Evaluate statistical dependence between gene mutation
##   events across patients and identify significant
##   co-occurring or mutually exclusive gene pairs.
##
##   Generate analysis-ready result tables and visualization
##   outputs for downstream interpretation.
##
##
## Input:
##   - result/data/mae_mut_clin.RData
##       MultiAssayExperiment containing:
##         - mut_binary assay:
##             gene × patient binary mutation matrix (0/1)
##         - colData:
##             harmonized clinical metadata (including IDH_status)
##
##
## Outputs:
##
##   1) result/coMut/coMut_res_all.csv
##        Pairwise co-mutation results (all patients)
##
##   2) result/coMut/coMut_res_wt.csv
##        Pairwise co-mutation results (IDH wild-type)
##
##   3) result/coMut/coMut_res_mut.csv
##        Pairwise co-mutation results (IDH mutant)
##
##        Each file contains:
##          - gene1, gene2
##          - OR (corrected odds ratio)
##          - pval, fdr
##          - n11, n10, n01, n00
##
##   4) Volcano plots (PDF):
##        - volcano_coMut_all_updated.pdf
##        - volcano_coMut_wt_updated.pdf
##        - volcano_coMut_mut_updated.pdf
##
##        Visualization:
##          x-axis: log2(OR)
##          y-axis: −log10(p-value)
##          significant pairs highlighted (FDR < 0.05)
##
##   5) UpSet plots (PDF):
##        - upset_coMut_all.pdf
##        - upset_coMut_wt.pdf
##        - upset_coMut_mut.pdf
##
##        Visualization:
##          mutation combination patterns across top
##          co-mutated genes ranked by interaction degree
##
##
## Processing Overview:
##
##   1) Data Loading
##        The MultiAssayExperiment object is loaded and the
##        binary mutation matrix and clinical metadata are
##        extracted.
##
##        IDH status is re-derived at the mutation level
##        (IDH1/IDH2) and appended as a gene-level feature.
##
##   2) Gene Filtering
##        Genes with mutation frequency ≥ min_mut_freq
##        (default = 5 patients) are retained to ensure
##        stability of statistical estimates.
##
##   3) Pairwise Co-Mutation Testing
##        All possible gene pairs are enumerated.
##
##        For each pair, a 2×2 contingency table is constructed:
##
##              g2 Mut   g2 WT
##     g1 Mut      n11     n10
##     g1 WT       n01     n00
##
##        Fisher’s exact test is used to evaluate independence.
##
##   4) Effect Size Estimation
##        Odds ratios (OR) are computed for each gene pair.
##        A continuity-corrected OR is applied to avoid
##        instability in sparse tables.
##
##           OR > 1  → co-occurrence enrichment
##           OR < 1  → mutual exclusivity
##
##        Effect sizes are visualized using log2(OR).
##
##   5) Multiple Testing Correction
##        P-values are adjusted using the Benjamini–Hochberg
##        procedure to control false discovery rate (FDR).
##
##   6) Stratified Analyses
##        Analyses are repeated for:
##           - All patients
##           - IDH wild-type subset
##           - IDH mutant subset
##
##        Mutation matrices are subset accordingly prior
##        to pairwise testing.
##
##   7) Visualization
##        Volcano plots are generated to display effect size
##        and statistical significance of gene pairs.
##
##        UpSet plots are constructed using top-ranked genes
##        to visualize higher-order mutation co-occurrence
##        patterns across patients.
##
##   8) Export
##        All results and figures are exported to the output
##        directory for downstream analysis and reporting.
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

mut_counts <- rowSums(mut == 1)
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

mut_counts <- rowSums(mut_wt == 1)
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

pdf(file.path(dir_output,  'volcano_coMut_all_updated.pdf'), width = 5, height = 4)

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
    axis.title = element_text(size = 9),
    axis.text  = element_text(size = 8),
   # plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
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

pdf(file.path(dir_output,  'volcano_coMut_wt_updated.pdf'), width = 5, height = 4)

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
    axis.title = element_text(size = 9),
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

pdf(file.path(dir_output,  'volcano_coMut_mut_updated.pdf'), width = 5, height = 4)

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
    axis.title = element_text(size = 9),
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
co_res <- co_res[order(co_res$pval), ]
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
mainbar.y.label <- "Number of Patients"
sets.x.label <- "Mutation Frequency"

pdf(file.path(dir_output,  'upset_coMut_all.pdf'), width = 7, height = 5)

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

pdf(file.path(dir_output,  'upset_coMut_wt.pdf'), width = 7, height = 5)

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

pdf(file.path(dir_output,  'upset_coMut_mut.pdf'), width = 7, height = 5)

suppressWarnings(
  upset.fun(fromList(listInput))
)

dev.off()
