#####################################################
## Script: CNS NGS Mutation–Clinical Association
##         (Bar Plots + Fisher’s Exact Test)
##
## Purpose:
##   Evaluate associations between recurrent gene-level
##   mutation status and key clinical variables in a
##   CNS cohort, and generate publication-ready bar plots.
##
##   The script:
##     • Selects recurrent mutations (≥5% frequency)
##     • Computes mutation frequencies by clinical groups
##     • Performs Fisher’s exact tests
##     • Applies multiple-testing correction (BH FDR)
##     • Generates stratified visualizations
##
## Input:
##   - result/data/mae_mut_clin.RData
##       MultiAssayExperiment containing:
##         assay:
##           gene × sample binary mutation matrix (0/1)
##         colData:
##           harmonized clinical metadata
##
## Output:
##
##   Association tables (CSV):
##     - result/clinicalBarPlot/mut_age.csv
##     - result/clinicalBarPlot/mut_age_mut.csv
##     - result/clinicalBarPlot/mut_age_wt.csv
##     - result/clinicalBarPlot/mut_sex.csv
##     - result/clinicalBarPlot/mut_sex_mut.csv
##     - result/clinicalBarPlot/mut_sex_wt.csv
##
##   Bar plot figures (PDF):
##     - mut_age.pdf
##     - mut_age_mut.pdf
##     - mut_age_wt.pdf
##     - mut_sex.pdf
##     - mut_sex_mut.pdf
##     - mut_sex_wt.pdf
##
## Key Processing Steps:
##
##   1) Clinical Harmonization
##        - Age dichotomized (<40 vs ≥40)
##        - IDH status recoded (WT / Mut / Unknown)
##
##   2) Mutation Filtering
##        - Calculate gene-level mutation frequency
##        - Retain recurrent genes (≥5% prevalence)
##
##   3) Data Reshaping
##        - Convert mutation matrix to long format
##        - Merge with clinical metadata
##
##   4) Statistical Testing
##        - Perform Fisher’s exact test for:
##            • Age vs mutation status
##            • Sex vs mutation status
##        - Conduct analyses for:
##            • All patients
##            • IDH mutant subset
##            • IDH wild-type subset
##        - Adjust p-values using Benjamini–Hochberg FDR
##
##   5) Visualization
##        - Generate grouped bar plots showing mutation
##          frequency (%) by clinical category
##        - Annotate significance levels (*, **, ***)
##        - Export figures as PDF files
#####################################################
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

pdf(file.path(dir_output, 'mut_age.pdf'), width = 7, height = 4)

ggplot(plot_age, aes(x = Gene, y = freq, fill = Age)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
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
        axis.title = element_text(size = 10),
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

pdf(file.path(dir_output, 'mut_age_mut.pdf'), width = 7, height = 4)

ggplot(plot_age, aes(x = Gene, y = freq, fill = Age)) +
  geom_bar(stat = "identity", position = "dodge") +
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
        axis.title = element_text(size = 10),
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

pdf(file.path(dir_output, 'mut_age_wt.pdf'), width = 7, height = 4)

ggplot(plot_age, aes(x = Gene, y = freq, fill = Age)) +
  geom_bar(stat = "identity", position = "dodge") +
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
        axis.title = element_text(size = 10),
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

pdf(file.path(dir_output, 'mut_sex.pdf'), width = 7, height = 4)

ggplot(plot_sex, aes(x = Gene, y = freq, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c("Female" = "#DD8452",
               "Male"   = "#4C72B0")
  ) +
  theme_classic() +
   theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
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

pdf(file.path(dir_output, 'mut_sex_mut.pdf'), width = 7, height = 4)

ggplot(plot_sex, aes(x = Gene, y = freq, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c("Female" = "#DD8452",
               "Male"   = "#4C72B0")
  ) +
  theme_classic() +
   theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
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

pdf(file.path(dir_output, 'mut_sex_wt.pdf'), width = 7, height = 4)

ggplot(plot_sex, aes(x = Gene, y = freq, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c("Female" = "#DD8452",
               "Male"   = "#4C72B0")
  ) +
  theme_classic() +
   theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  labs(y = "Mutation frequency (%)", x = "", 
       title = " ")

dev.off()

