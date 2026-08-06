##-------------------------------------------------------------------
## Script: CNS NGS Mutation–Clinical Association Analysis
## Purpose:
##   To evaluate associations between recurrent gene-level
##   mutations and selected clinical variables in CNS
##   tumour patients using curated NGS and clinical data.
##
##   Statistical association analyses and mutation
##   frequency visualizations are generated for the
##   full cohort and IDH- and histology-stratified
##   subgroups.
##
## Input:
##   - result/data/mae_mut_clin.RData
##       MultiAssayExperiment object containing:
##         • binary mutation matrix
##         • harmonized clinical metadata
##
## Outputs:
##   1) Statistical association result tables
##   2) Mutation frequency bar plot visualizations
##
## Processing Overview:
##   1) Load mutation and clinical data
##   2) Harmonize clinical variables including:
##        - Age group
##        - Sex
##        - IDH status
##        - Histology
##   3) Combine IDH1 and IDH2 mutation calls into
##      a single IDH mutation indicator
##   4) Filter recurrent mutations (≥5% frequency)
##      and prepare patient-level datasets
##   5) Perform Fisher's exact tests between
##      mutation status and:
##        - Age group
##        - Sex
##   6) Repeat analyses within:
##        - Full cohort
##        - IDH wild-type subgroup
##        - IDH mutant subgroup
##        - IDH wild-type glioblastoma subgroup
##        - IDH wild-type non-glioblastoma subgroup
##   7) Apply false discovery rate (FDR)
##      adjustment
##   8) Generate mutation frequency bar plots
##   9) Export statistical results and figures
##
## Notes:
##   - Analyses are based on binary gene-level
##     mutation matrices.
##   - IDH1 and IDH2 mutations are combined into
##     a single IDH mutation indicator.
##   - Recurrent genes are defined using a minimum
##     mutation frequency threshold of 5%.
##   - Multiple testing correction is performed
##     using false discovery rate (FDR) adjustment.
##   - Outputs are intended for downstream
##     visualization and reporting.
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
      select(Sample, Age, Sex, IDH_status, histo),
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
age_assoc$pval <- round(age_assoc$pval, 1)

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
      pval < 0.001 ~ "****",
      pval < 0.01  ~ "***",
      pval < 0.05  ~ "**",
      pval <= 0.1  ~ "*",
      TRUE        ~ ""
    )
  )

plot_age$Gene <- factor(plot_age$Gene, levels = gene_order)
plot_age_sig <- plot_age %>%
  group_by(Gene) %>%
  summarise(y_pos = max(freq) + 5, .groups = "drop") %>%
  left_join(age_assoc %>% select(Gene, sig), by = "Gene")

p <- ggplot(plot_age, aes(x = Gene, y = freq, fill = Age)) +
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

pdf(file.path(dir_output, 'mut_age.pdf'), width = 5, height = 3.5)
print(p)
dev.off()

jpeg(file.path(dir_output, "mut_age.jpeg"), width = 5, height = 3.5, units = "in", 
     res = 300, quality = 100)
print(p)
dev.off()

## Step 2 --- Mutant patients
mut_long_mut <- mut_long[mut_long$IDH_status == 'Mut', ] %>%
  group_by(Gene) %>%
  filter(sum(Mutation) > 0) %>%
  ungroup()
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
age_assoc$pval <- round(age_assoc$pval, 1)

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
      pval < 0.001 ~ "****",
      pval < 0.01  ~ "***",
      pval < 0.05  ~ "**",
      pval <= 0.1  ~ "*",
      TRUE        ~ ""
    )
  )

plot_age$Gene <- factor(plot_age$Gene, levels = gene_order)
plot_age_sig <- plot_age %>%
  group_by(Gene) %>%
  summarise(y_pos = max(freq) + 5, .groups = "drop") %>%
  left_join(age_assoc %>% select(Gene, sig), by = "Gene")


p <- ggplot(plot_age, aes(x = Gene, y = freq, fill = Age)) +
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

pdf(file.path(dir_output, 'mut_age_mut.pdf'), width = 5.5, height = 3.5)
print(p)
dev.off()

jpeg(file.path(dir_output, "mut_age_mut.jpeg"), width = 5.5, height = 3.5, units = "in", 
     res = 300, quality = 100)
print(p)
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
age_assoc$pval <- round(age_assoc$pval, 1)

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
      pval < 0.001 ~ "****",
      pval < 0.01  ~ "***",
      pval < 0.05  ~ "**",
      pval <= 0.1  ~ "*",
      TRUE        ~ ""
    )
  )

plot_age$Gene <- factor(plot_age$Gene, levels = gene_order)
plot_age_sig <- plot_age %>%
  group_by(Gene) %>%
  summarise(y_pos = max(freq) + 5, .groups = "drop") %>%
  left_join(age_assoc %>% select(Gene, sig), by = "Gene")

p <- ggplot(plot_age, aes(x = Gene, y = freq, fill = Age)) +
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

pdf(file.path(dir_output, 'mut_age_wt.pdf'), width = 6, height = 3.5)
print(p)
dev.off()

jpeg(file.path(dir_output, "mut_age_wt.jpeg"), width = 6, height = 3.5, units = "in", 
     res = 300, quality = 100)
print(p)
dev.off()

## Step 4 --- Wild type GBM patients
mut_long_wt_gbm <- mut_long[
  mut_long$IDH_status == 'WT' & mut_long$histo == 'Glioblastoma',
] %>%
  filter(Gene != "IDH")
age_assoc <- mut_long_wt_gbm %>%
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
write.csv(age_assoc, file = file.path(dir_output, 'csv', 'mut_age_wt_gbm.csv'))
age_assoc$pval <- round(age_assoc$pval, 1)

## bar plot age
gene_order <- mut_long_wt_gbm %>%
  group_by(Gene) %>%
  summarise(overall = mean(Mutation)) %>%
  arrange(desc(overall)) %>%
  pull(Gene)

plot_age <- mut_long_wt_gbm %>%
  group_by(Gene, Age) %>%
  summarise(freq = mean(Mutation) * 100, .groups = "drop")

age_assoc <- age_assoc %>%
  mutate(
    sig = case_when(
      pval < 0.001 ~ "****",
      pval < 0.01  ~ "***",
      pval < 0.05  ~ "**",
      pval <= 0.1  ~ "*",
      TRUE        ~ ""
    )
  )

plot_age$Gene <- factor(plot_age$Gene, levels = gene_order)
plot_age_sig <- plot_age %>%
  group_by(Gene) %>%
  summarise(y_pos = max(freq) + 5, .groups = "drop") %>%
  left_join(age_assoc %>% select(Gene, sig), by = "Gene")

p <- ggplot(plot_age, aes(x = Gene, y = freq, fill = Age)) +
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

pdf(file.path(dir_output, 'mut_age_wt_gbm.pdf'), width = 6, height = 3.5)
print(p)
dev.off()

jpeg(file.path(dir_output, "mut_age_wt_gbm.jpeg"), width = 6, height = 3.5, units = "in", 
     res = 300, quality = 100)
print(p)
dev.off()

## Step 5 --- Wild type nonGBM patients
mut_long_wt_nongbm <- mut_long[
  mut_long$IDH_status == 'WT' & mut_long$histo != 'Glioblastoma',
] %>%
  group_by(Gene) %>%
  filter(sum(Mutation) > 0) %>%
  ungroup()
age_assoc <- mut_long_wt_nongbm %>%
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
write.csv(age_assoc, file = file.path(dir_output, 'csv', 'mut_age_wt_nongbm.csv'))
age_assoc$pval <- round(age_assoc$pval, 1)

## bar plot age
gene_order <- mut_long_wt_nongbm %>%
  group_by(Gene) %>%
  summarise(overall = mean(Mutation)) %>%
  arrange(desc(overall)) %>%
  pull(Gene)

plot_age <- mut_long_wt_nongbm %>%
  group_by(Gene, Age) %>%
  summarise(freq = mean(Mutation) * 100, .groups = "drop")

age_assoc <- age_assoc %>%
  mutate(
    sig = case_when(
      pval < 0.001 ~ "****",
      pval < 0.01  ~ "***",
      pval < 0.05  ~ "**",
      pval <= 0.1  ~ "*",
      TRUE        ~ ""
    )
  )

plot_age$Gene <- factor(plot_age$Gene, levels = gene_order)
plot_age_sig <- plot_age %>%
  group_by(Gene) %>%
  summarise(y_pos = max(freq) + 5, .groups = "drop") %>%
  left_join(age_assoc %>% select(Gene, sig), by = "Gene")

p <- ggplot(plot_age, aes(x = Gene, y = freq, fill = Age)) +
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

pdf(file.path(dir_output, 'mut_age_wt_nongbm.pdf'), width = 5, height = 3.5)
print(p)
dev.off()

jpeg(file.path(dir_output, "mut_age_wt_nongbm.jpeg"), width = 5, height = 3.5, units = "in", 
     res = 300, quality = 100)
print(p)
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
sex_assoc$pval <- round(sex_assoc$pval, 1)

## bar plot sex
## bar plot sex
plot_sex <- mut_long %>%
  group_by(Gene, Sex) %>%
  summarise(freq = mean(Mutation) * 100, .groups = "drop")


sex_assoc <- sex_assoc %>%
  mutate(
    sig = case_when(
      pval < 0.001 ~ "****",
      pval < 0.01  ~ "***",
      pval < 0.05  ~ "**",
      pval <= 0.1  ~ "*",
      TRUE        ~ ""
    )
  )

plot_sex$Gene <- factor(plot_sex$Gene, levels = gene_order)
plot_sex_sig <- plot_sex %>%
  group_by(Gene) %>%
  summarise(y_pos = max(freq) + 5, .groups = "drop") %>%
  left_join(sex_assoc %>% select(Gene, sig), by = "Gene")


p <- ggplot(plot_sex, aes(x = Gene, y = freq, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(
    values = c("Female" = "#DD8452",
               "Male"   = "#4C72B0")
  ) +
  geom_text(
    data = plot_sex_sig,
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

pdf(file.path(dir_output, 'mut_sex.pdf'), width = 6, height = 3.5)
print(p)
dev.off()

jpeg(file.path(dir_output, "mut_sex.jpeg"), width = 6, height = 3.5, units = "in", 
     res = 300, quality = 100)
print(p)
dev.off()

## Step 2 --- Mutant patients
mut_long_mut <- mut_long[mut_long$IDH_status == 'Mut', ] %>%
  group_by(Gene) %>%
  filter(sum(Mutation) > 0) %>%
  ungroup()
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
sex_assoc$pval <- round(sex_assoc$pval, 1)

## bar plot sex
## bar plot sex
plot_sex <- mut_long_mut %>%
  group_by(Gene, Sex) %>%
  summarise(freq = mean(Mutation) * 100, .groups = "drop")

sex_assoc <- sex_assoc %>%
  mutate(
    sig = case_when(
      pval < 0.001 ~ "****",
      pval < 0.01  ~ "***",
      pval < 0.05  ~ "**",
      pval <= 0.1  ~ "*",
      TRUE        ~ ""
    )
  )

plot_sex$Gene <- factor(plot_sex$Gene, levels = gene_order)
plot_sex_sig <- plot_sex %>%
  group_by(Gene) %>%
  summarise(y_pos = max(freq) + 5, .groups = "drop") %>%
  left_join(sex_assoc %>% select(Gene, sig), by = "Gene")


p <- ggplot(plot_sex, aes(x = Gene, y = freq, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(
    values = c("Female" = "#DD8452",
               "Male"   = "#4C72B0")
  ) +
  geom_text(
    data = plot_sex_sig,
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

pdf(file.path(dir_output,  'mut_sex_mut.pdf'), width = 5.5, height = 3.5)
print(p)
dev.off()

jpeg(file.path(dir_output, "mut_sex_mut.jpeg"), width = 5.5, height = 3.5, units = "in", 
     res = 300, quality = 100)
print(p)
dev.off()

## Step 3 --- Wild type patients
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
sex_assoc$pval <- round(sex_assoc$pval, 1)

## bar plot sex
## bar plot sex
plot_sex <- mut_long_wt %>%
  group_by(Gene, Sex) %>%
  summarise(freq = mean(Mutation) * 100, .groups = "drop")

sex_assoc <- sex_assoc %>%
  mutate(
    sig = case_when(
      pval < 0.001 ~ "****",
      pval < 0.01  ~ "***",
      pval < 0.05  ~ "**",
      pval <= 0.1  ~ "*",
      TRUE        ~ ""
    )
  )

plot_sex$Gene <- factor(plot_sex$Gene, levels = gene_order)
plot_sex_sig <- plot_sex %>%
  group_by(Gene) %>%
  summarise(y_pos = max(freq) + 5, .groups = "drop") %>%
  left_join(sex_assoc %>% select(Gene, sig), by = "Gene")

p <- ggplot(plot_sex, aes(x = Gene, y = freq, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(
    values = c("Female" = "#DD8452",
               "Male"   = "#4C72B0")
  ) +
  geom_text(
    data = plot_sex_sig,
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

pdf(file.path(dir_output, 'mut_sex_wt.pdf'), width = 6, height = 3.5)
print(p)
dev.off()

jpeg(file.path(dir_output, "mut_sex_wt.jpeg"), width = 6, height = 3.5, units = "in", 
     res = 300, quality = 100)
print(p)
dev.off()

## Step 4 --- Wild type GBM patients
mut_long_wt_gbm <- mut_long[
  mut_long$IDH_status == 'WT' & mut_long$histo == 'Glioblastoma',
] %>%
  filter(Gene != "IDH")
sex_assoc <- mut_long_wt_gbm %>%
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
write.csv(sex_assoc, file = file.path(dir_output, 'csv', 'mut_sex_wt_gbm.csv'))
sex_assoc$pval <- round(sex_assoc$pval, 1)

## bar plot sex
plot_sex <- mut_long_wt_gbm %>%
  group_by(Gene, Sex) %>%
  summarise(freq = mean(Mutation) * 100, .groups = "drop")

sex_assoc <- sex_assoc %>%
  mutate(
    sig = case_when(
      pval < 0.001 ~ "****",
      pval < 0.01  ~ "***",
      pval < 0.05  ~ "**",
      pval <= 0.1  ~ "*",
      TRUE        ~ ""
    )
  )

plot_sex$Gene <- factor(plot_sex$Gene, levels = gene_order)
plot_sex_sig <- plot_sex %>%
  group_by(Gene) %>%
  summarise(y_pos = max(freq) + 5, .groups = "drop") %>%
  left_join(sex_assoc %>% select(Gene, sig), by = "Gene")


p <- ggplot(plot_sex, aes(x = Gene, y = freq, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(
    values = c("Female" = "#DD8452",
               "Male"   = "#4C72B0")
  ) +
  geom_text(
    data = plot_sex_sig,
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

pdf(file.path(dir_output, 'mut_sex_wt_gbm.pdf'), width = 6, height = 3.5)
print(p)
dev.off()

jpeg(file.path(dir_output, "mut_sex_wt_gbm.jpeg"), width = 6, height = 3.5, units = "in", 
     res = 300, quality = 100)
print(p)
dev.off()

## Step 5 --- Wild type nonGBM patients
mut_long_wt_nongbm <- mut_long[
  mut_long$IDH_status == 'WT' & mut_long$histo != 'Glioblastoma',
] %>%
  group_by(Gene) %>%
  filter(sum(Mutation) > 0) %>%
  ungroup()
sex_assoc <- mut_long_wt_nongbm %>%
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
write.csv(sex_assoc, file = file.path(dir_output, 'csv', 'mut_sex_wt_nongbm.csv'))
sex_assoc$pval <- round(sex_assoc$pval, 1)

## bar plot sex
## bar plot sex
plot_sex <- mut_long_wt_nongbm %>%
  group_by(Gene, Sex) %>%
  summarise(freq = mean(Mutation) * 100, .groups = "drop")

sex_assoc <- sex_assoc %>%
  mutate(
    sig = case_when(
      pval < 0.001 ~ "****",
      pval < 0.01  ~ "***",
      pval < 0.05  ~ "**",
      pval <= 0.1  ~ "*",
      TRUE        ~ ""
    )
  )

plot_sex$Gene <- factor(plot_sex$Gene, levels = gene_order)
plot_sex_sig <- plot_sex %>%
  group_by(Gene) %>%
  summarise(y_pos = max(freq) + 5, .groups = "drop") %>%
  left_join(sex_assoc %>% select(Gene, sig), by = "Gene")


p <- ggplot(plot_sex, aes(x = Gene, y = freq, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(
    values = c("Female" = "#DD8452",
               "Male"   = "#4C72B0")
  ) +
  geom_text(
    data = plot_sex_sig,
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


pdf(file.path(dir_output, 'mut_sex_wt_nongbm.pdf'), width = 5, height = 3.5)
print(p)
dev.off()

jpeg(file.path(dir_output, "mut_sex_wt_nongbm.jpeg"), width = 5, height = 3.5, units = "in", 
     res = 300, quality = 100)
print(p)
dev.off()
