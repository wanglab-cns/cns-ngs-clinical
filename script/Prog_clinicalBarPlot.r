##-------------------------------------------------------------------
## Script: CNS NGS Mutation–Clinical Association Analysis
##
## Purpose:
##   Evaluate associations between gene-level mutation
##   alterations and clinical variables in CNS tumor patients
##   using curated NGS and clinical data.
##
##   Identify statistically significant relationships between
##   mutation status and key clinical features, including
##   molecular subtype, histology, and demographic variables.
##
##   Generate analysis-ready summaries and results for
##   downstream interpretation and visualization.
##
##
## Input:
##   - result/data/mae_mut_clin.RData
##       MultiAssayExperiment containing:
##         - mut_binary assay:
##             gene × patient binary mutation matrix (0/1)
##         - colData:
##             curated clinical metadata
##
##
## Outputs:
##
##   1) result/association/mutation_clinical_association_results.csv
##        Table of statistical test results including:
##          - gene
##          - clinical variable
##          - test statistic
##          - p-value
##          - adjusted p-value
##
##   2) result/association/mutation_clinical_summary.csv
##        Summary table of mutation frequencies stratified
##        by clinical groups
##
##
## Processing Overview:
##
##   1) Data Loading
##        The MultiAssayExperiment object is loaded and the
##        binary mutation matrix and corresponding clinical
##        metadata are extracted.
##
##   2) Data Preparation
##        Mutation data are structured as gene × patient binary
##        indicators. Clinical variables of interest are selected
##        and formatted (e.g., categorical encoding, grouping).
##
##   3) Mutation Frequency Filtering
##        Genes with low mutation frequency may be filtered to
##        ensure sufficient statistical power for association
##        testing.
##
##   4) Association Testing
##        For each gene and clinical variable:
##           - appropriate statistical tests are applied
##             (e.g., Fisher’s exact test or chi-squared test)
##           - associations between mutation status and
##             clinical groups are evaluated
##
##   5) Multiple Testing Correction
##        P-values are adjusted across tests to control for
##        false discovery rate (e.g., Benjamini–Hochberg).
##
##   6) Result Summarization
##        Significant associations are summarized and mutation
##        frequencies are reported across clinical subgroups.
##
##   7) Export
##        Results and summary tables are written to output files
##        for downstream analysis and visualization.
##-------------------------------------------------------------------
####################################################
## Load libraries
####################################################
library(MultiAssayExperiment)
library(ComplexHeatmap)
library(dplyr) 
library(tidyr)
library(stringr)
library(paletteer)
library(ggplot2)

####################################################
## Setup directories
####################################################
dir_input <- 'result/data'
dir_output <- 'result/barPlot'

#################################################
## Load data
#################################################
load(file.path(dir_input, 'mae_mut_clin.RData'))
mut <- assay(mae[['mut_binary']])
clin <- as.data.frame(colData(mae[['mut_binary']]))
clin$Age <- ifelse(clin$Age >= 40, '>40', '<40')

## ---- Binary matrix
if(all(c("IDH1","IDH2") %in% rownames(mut))){

  idh_vec <- as.integer(
    mut["IDH1", ] == 1 | mut["IDH2", ] == 1
  )
  mut <- mut[!rownames(mut) %in% c("IDH1","IDH2"), , drop = FALSE]
  mut <- rbind(mut, IDH = idh_vec)

}

###################################################
## Bar plots: top mutations and sex/age
###################################################
## choose top mutations
mut_freq <- rowMeans(mut == 1)

genes <- names(mut_freq[mut_freq >= 0.05])
mut_filt <- mut[genes, ]

## convert to long format
mut_long <- as.data.frame(t(mut_filt)) %>%
  mutate(Sample = rownames(.)) %>%
  pivot_longer(
    cols = -Sample,
    names_to = "Gene",
    values_to = "Mutation"
  )

## merge with clinical age and sex variables
mut_long <- mut_long %>%
  left_join(
    clin %>% 
      mutate(Sample = Study) %>%
      select(Sample, Age, Sex, IDH_status),
    by = "Sample"
  )

###########################################################
## Fishe exact test --- Age and mutations
###########################################################
## Step 1 --- all patients
age_assoc <- mut_long %>%
  group_by(Gene) %>%
  summarise(
    pval = fisher.test(table(Mutation, Age))$p.value
  )

age_assoc$fdr <- p.adjust(age_assoc$pval, method = "BH")
write.csv(age_assoc, file = file.path(dir_output, 'mut_age.csv'), row.names=FALSE)

## bar plot age
gene_order <- mut_long %>%
  group_by(Gene) %>%
  summarise(overall = mean(Mutation)) %>%
  arrange(desc(overall)) %>%
  pull(Gene)

plot_age <- mut_long %>%
  group_by(Gene, Age) %>%
  summarise(freq = mean(Mutation) * 100, .groups = "drop")

age_assoc <- age_assoc %>%
  mutate(
    sig = case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE        ~ ""
    )
  )

plot_age$Gene <- factor(plot_age$Gene, levels = gene_order)
plot_age_sig <- plot_age %>%
  group_by(Gene) %>%
  summarise(y_pos = max(freq) + 5, .groups = "drop") %>%
  left_join(age_assoc %>% select(Gene, sig), by = "Gene")

pdf(file.path(dir_output, 'mut_age.pdf'), width = 6, height = 3.5)

ggplot(plot_age, aes(x = Gene, y = freq, fill = Age)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(
    values = c(">40" = "#08306B",
               "<40" = "#C6DBEF")
  ) +
  geom_text(
    data = plot_age_sig,
    aes(x = Gene, y = y_pos, label = sig),
    inherit.aes = FALSE,
    size = 5
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8)) +
  labs(y = "Mutation frequency (%)", x = "", 
       title = " ")

dev.off()

## Step 2 --- Mutant patients
mut_long_mut <- mut_long[mut_long$IDH_status == 'Mut', ]
age_assoc <- mut_long_mut %>%
  group_by(Gene) %>%
  summarise(
    pval = {
      tbl <- table(Mutation, Age)
      if (all(dim(tbl) == c(2, 2))) {
        fisher.test(tbl)$p.value
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  )

age_assoc <- age_assoc[!is.na(age_assoc$pval), ]
age_assoc$fdr <- p.adjust(age_assoc$pval, method = "BH")
write.csv(age_assoc, file = file.path(dir_output, 'mut_age_mut.csv'))

## bar plot age
gene_order <- mut_long_mut %>%
  group_by(Gene) %>%
  summarise(overall = mean(Mutation)) %>%
  arrange(desc(overall)) %>%
  pull(Gene)

plot_age <- mut_long_mut %>%
  group_by(Gene, Age) %>%
  summarise(freq = mean(Mutation) * 100, .groups = "drop")

age_assoc <- age_assoc %>%
  mutate(
    sig = case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE        ~ ""
    )
  )

plot_age$Gene <- factor(plot_age$Gene, levels = gene_order)
plot_age_sig <- plot_age %>%
  group_by(Gene) %>%
  summarise(y_pos = max(freq) + 5, .groups = "drop") %>%
  left_join(age_assoc %>% select(Gene, sig), by = "Gene")

pdf(file.path(dir_output, 'mut_age_mut.pdf'), width = 6, height = 3.5)

ggplot(plot_age, aes(x = Gene, y = freq, fill = Age)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(
    values = c(">40" = "#08306B",
               "<40" = "#C6DBEF")
  ) +
  geom_text(
    data = plot_age_sig,
    aes(x = Gene, y = y_pos, label = sig),
    inherit.aes = FALSE,
    size = 5
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8)) +
  labs(y = "Mutation frequency (%)", x = "", 
       title = " ")

dev.off()

## Step 3 --- Wild type patients
mut_long_wt <- mut_long[mut_long$IDH_status == 'WT', ]
age_assoc <- mut_long_wt %>%
  group_by(Gene) %>%
  summarise(
    pval = {
      tbl <- table(Mutation, Age)
      if (all(dim(tbl) == c(2, 2))) {
        fisher.test(tbl)$p.value
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  )

age_assoc <- age_assoc[!is.na(age_assoc$pval), ]
age_assoc$fdr <- p.adjust(age_assoc$pval, method = "BH")
write.csv(age_assoc, file = file.path(dir_output, 'mut_age_wt.csv'))

## bar plot age
gene_order <- mut_long_wt %>%
  group_by(Gene) %>%
  summarise(overall = mean(Mutation)) %>%
  arrange(desc(overall)) %>%
  pull(Gene)

plot_age <- mut_long_wt %>%
  group_by(Gene, Age) %>%
  summarise(freq = mean(Mutation) * 100, .groups = "drop")

age_assoc <- age_assoc %>%
  mutate(
    sig = case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE        ~ ""
    )
  )

plot_age$Gene <- factor(plot_age$Gene, levels = gene_order)
plot_age_sig <- plot_age %>%
  group_by(Gene) %>%
  summarise(y_pos = max(freq) + 5, .groups = "drop") %>%
  left_join(age_assoc %>% select(Gene, sig), by = "Gene")

pdf(file.path(dir_output, 'mut_age_wt.pdf'), width = 6, height = 3.5)

ggplot(plot_age, aes(x = Gene, y = freq, fill = Age)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(
    values = c(">40" = "#08306B",
               "<40" = "#C6DBEF")
  ) +
  geom_text(
    data = plot_age_sig,
    aes(x = Gene, y = y_pos, label = sig),
    inherit.aes = FALSE,
    size = 5
  ) +
  theme_classic() +
   theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8)) +
  labs(y = "Mutation frequency (%)", x = "", 
       title = " ")

dev.off()

###########################################################
## Fishe exact test --- Sex and mutations 
###########################################################
## Step 1 --- all patients
sex_assoc <- mut_long %>%
  group_by(Gene) %>%
  summarise(
    pval = fisher.test(table(Mutation, Sex))$p.value
  )

sex_assoc$fdr <- p.adjust(sex_assoc$pval, method = "BH")
write.csv(sex_assoc, file = file.path(dir_output, 'mut_sex.csv'))

## bar plot sex
plot_sex <- mut_long %>%
  group_by(Gene, Sex) %>%
  summarise(freq = mean(Mutation) * 100, .groups = "drop")

plot_sex$Gene <- factor(plot_sex$Gene, levels = gene_order)

pdf(file.path(dir_output, 'mut_sex.pdf'), width = 6, height = 3.5)

ggplot(plot_sex, aes(x = Gene, y = freq, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(
    values = c("Female" = "#DD8452",
               "Male"   = "#4C72B0")
  ) +
  theme_classic() +
   theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8)) +
  labs(y = "Mutation frequency (%)", x = "", 
       title = " ")

dev.off()

## Step 2 --- mutant patients
mut_long_mut <- mut_long[mut_long$IDH_status == 'Mut', ]
sex_assoc <- mut_long_mut %>%
  group_by(Gene) %>%
  summarise(
    pval = {
      tbl <- table(Mutation, Sex)
      if (all(dim(tbl) == c(2, 2))) {
        fisher.test(tbl)$p.value
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  )

sex_assoc <- sex_assoc[!is.na(sex_assoc$pval), ]
sex_assoc$fdr <- p.adjust(sex_assoc$pval, method = "BH")
write.csv(sex_assoc, file = file.path(dir_output, 'mut_sex_mut.csv'))

## bar plot sex
plot_sex <- mut_long_mut %>%
  group_by(Gene, Sex) %>%
  summarise(freq = mean(Mutation) * 100, .groups = "drop")

plot_sex$Gene <- factor(plot_sex$Gene, levels = gene_order)

pdf(file.path(dir_output, 'mut_sex_mut.pdf'), width = 6, height = 3.5)

ggplot(plot_sex, aes(x = Gene, y = freq, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(
    values = c("Female" = "#DD8452",
               "Male"   = "#4C72B0")
  ) +
  theme_classic() +
   theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8)) +
  labs(y = "Mutation frequency (%)", x = "", 
       title = " ")

dev.off()

## Step 2 --- mutant patients
mut_long_wt <- mut_long[mut_long$IDH_status == 'WT', ]
sex_assoc <- mut_long_wt %>%
  group_by(Gene) %>%
  summarise(
    pval = {
      tbl <- table(Mutation, Sex)
      if (all(dim(tbl) == c(2, 2))) {
        fisher.test(tbl)$p.value
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  )

sex_assoc <- sex_assoc[!is.na(sex_assoc$pval), ]
sex_assoc$fdr <- p.adjust(sex_assoc$pval, method = "BH")
write.csv(sex_assoc, file = file.path(dir_output, 'mut_sex_wt.csv'))

## bar plot sex
plot_sex <- mut_long_wt %>%
  group_by(Gene, Sex) %>%
  summarise(freq = mean(Mutation) * 100, .groups = "drop")

plot_sex$Gene <- factor(plot_sex$Gene, levels = gene_order)

pdf(file.path(dir_output, 'mut_sex_wt.pdf'), width = 6, height = 3.5)

ggplot(plot_sex, aes(x = Gene, y = freq, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(
    values = c("Female" = "#DD8452",
               "Male"   = "#4C72B0")
  ) +
  theme_classic() +
   theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8)) +
  labs(y = "Mutation frequency (%)", x = "", 
       title = " ")

dev.off()

