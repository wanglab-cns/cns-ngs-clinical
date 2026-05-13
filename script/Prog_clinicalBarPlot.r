##-------------------------------------------------------------------
## Script: CNS NGS Mutation–Clinical Association Analysis
##
## Purpose:
##   Evaluate associations between recurrent gene-level
##   mutation alterations and demographic or molecular
##   clinical variables in CNS tumor patients using
##   curated NGS and clinical datasets.
##
##   Perform statistical association testing between
##   mutation status and age or sex groups across the
##   full cohort and IDH-stratified patient subsets.
##
##   Generate publication-quality summary visualizations
##   and statistical result tables for downstream
##   interpretation and reporting.
##
##
## Input:
##
##   - result/data/mae_mut_clin.RData
##
##       MultiAssayExperiment object containing:
##
##         - mut_binary assay:
##             gene × patient binary mutation matrix (0/1)
##
##         - colData:
##             curated clinical metadata
##
##
## Outputs:
##
##   1) result/barPlot/csv/
##
##        Statistical association result tables:
##
##          - mut_age.csv
##          - mut_age_mut.csv
##          - mut_age_wt.csv
##          - mut_sex.csv
##          - mut_sex_mut.csv
##          - mut_sex_wt.csv
##
##
##   2) result/barPlot/
##
##        Publication-quality mutation frequency bar plots:
##
##          - mut_age.pdf
##          - mut_age_mut.pdf
##          - mut_age_wt.pdf
##          - mut_sex.pdf
##          - mut_sex_mut.pdf
##          - mut_sex_wt.pdf
##
##
## Processing Overview:
##
##   1) Data Loading
##
##        The MultiAssayExperiment object is loaded and the
##        binary mutation matrix and associated clinical
##        metadata are extracted.
##
##
##   2) Clinical Variable Processing
##
##        Clinical variables are standardized for analysis:
##
##           - Age dichotomized (<40 vs ≥40)
##           - Sex encoded as Male/Female
##           - IDH status extracted for subgroup analyses
##
##
##   3) Binary Mutation Matrix Processing
##
##        Gene-level mutation frequencies are calculated
##        from the binary mutation matrix.
##
##        IDH1 and IDH2 alterations are collapsed into a
##        unified patient-level IDH mutation variable.
##
##        Low-frequency genes (<5% mutation frequency)
##        are excluded from downstream visualization and
##        association analyses.
##
##
##   4) Data Restructuring
##
##        Mutation matrices are converted into long-format
##        patient-level tables and merged with clinical
##        metadata for downstream statistical testing and
##        visualization.
##
##
##   5) Association Testing
##
##        Gene-wise associations between mutation status
##        and clinical variables are evaluated using
##        Fisher’s exact tests.
##
##        Analyses are performed across:
##
##           - all patients
##           - IDH mutant patients
##           - IDH wild-type patients
##
##        Multiple-testing correction is performed using
##        the Benjamini–Hochberg false discovery rate method.
##
##
##   6) Mutation Frequency Visualization
##
##        Mutation frequencies are summarized as grouped
##        bar plots stratified by:
##
##           - Age group
##           - Sex
##
##        Features include:
##
##           - Consistent mutation ordering
##           - Cohort-specific subgroup analyses
##           - Significance annotations
##           - Harmonized color palettes
##           - Publication-quality formatting
##
##
##   7) Export
##
##        Statistical result tables and visualization
##        outputs are exported as CSV and PDF files for
##        downstream interpretation and reporting.
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
write.csv(age_assoc, file = file.path(dir_output, 'csv', 'mut_age.csv'), row.names=FALSE)

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

pdf(file.path(dir_output, 'mut_age.pdf'), width = 5, height = 3.5)

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
  guides(
  fill = guide_legend(
    keywidth  = unit(0.5, "cm"),
    keyheight = unit(0.5, "cm")
  )
) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 7, face= 'bold'),
        legend.text = element_text(size = 7)) +
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
write.csv(age_assoc, file = file.path(dir_output, 'csv', 'mut_age_mut.csv'))

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
   guides(
  fill = guide_legend(
    keywidth  = unit(0.5, "cm"),
    keyheight = unit(0.5, "cm")
  )
) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 7, face= 'bold'),
        legend.text = element_text(size = 7)) +
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
write.csv(age_assoc, file = file.path(dir_output, 'csv', 'mut_age_wt.csv'))

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
   guides(
  fill = guide_legend(
    keywidth  = unit(0.5, "cm"),
    keyheight = unit(0.5, "cm")
  )
) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 7, face= 'bold'),
        legend.text = element_text(size = 7)) +
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
write.csv(sex_assoc, file = file.path(dir_output, 'csv', 'mut_sex.csv'))

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
 guides(
  fill = guide_legend(
    keywidth  = unit(0.5, "cm"),
    keyheight = unit(0.5, "cm")
  )
) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 7, face= 'bold'),
        legend.text = element_text(size = 7)) +
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
write.csv(sex_assoc, file = file.path(dir_output, 'csv', 'mut_sex_mut.csv'))

## bar plot sex
plot_sex <- mut_long_mut %>%
  group_by(Gene, Sex) %>%
  summarise(freq = mean(Mutation) * 100, .groups = "drop")

plot_sex$Gene <- factor(plot_sex$Gene, levels = gene_order)

pdf(file.path(dir_output,  'mut_sex_mut.pdf'), width = 6, height = 3.5)

ggplot(plot_sex, aes(x = Gene, y = freq, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(
    values = c("Female" = "#DD8452",
               "Male"   = "#4C72B0")
  ) +
  guides(
  fill = guide_legend(
    keywidth  = unit(0.5, "cm"),
    keyheight = unit(0.5, "cm")
  )
) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 7, face= 'bold'),
        legend.text = element_text(size = 7)) +
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
write.csv(sex_assoc, file = file.path(dir_output, 'csv', 'mut_sex_wt.csv'))

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
 guides(
  fill = guide_legend(
    keywidth  = unit(0.5, "cm"),
    keyheight = unit(0.5, "cm")
  )
) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 7, face= 'bold'),
        legend.text = element_text(size = 7)) +
  labs(y = "Mutation frequency (%)", x = "", 
       title = " ")

dev.off()

