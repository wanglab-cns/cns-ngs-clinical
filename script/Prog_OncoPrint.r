#####################################################
## Script: CNS NGS OncoPrint (ComplexHeatmap)
##
## Purpose:
##   Generate publication-quality OncoPrints from a curated CNS NGS
##   mutation event matrix stored in a SummarizedExperiment object.
##   The script standardizes alteration labels, enforces within-cell
##   alteration priority, optionally collapses multi-event cells into
##   a single "Other" category, and visualizes mutation frequencies
##   with and without clinical annotations.
##
## Input:
##   - result/data/se_mut_onco_clin.RData
##       SummarizedExperiment containing:
##         assay: gene × sample character matrix
##                ("" or semicolon-delimited alteration types)
##         colData: sample-level clinical metadata
##
## Output:
##   - result/Fig/fig1.pdf
##       OncoPrint of mutation events only (no clinical metadata)
##   - result/Fig/fig2.pdf
##       OncoPrint with clinical annotations (Sex, Age group, IDH status)
##
## Key steps:
##   1) Load the SummarizedExperiment and extract the mutation matrix
##      and clinical metadata.
##   2) Standardize alteration labels to a controlled vocabulary:
##        snv/indel            → SNV/Indel
##        copy number variant → CNV
##        fusion              → Fusion
##      and enforce a consistent within-cell priority:
##        Fusion > CNV > SNV/Indel.
##   3) Optionally collapse combined alteration labels
##      (e.g., "Fusion;CNV", "CNV;SNV/Indel") into a single "Other"
##      category to simplify visualization and legends.
##   4) Define custom graphical functions for alteration rendering
##      and generate OncoPrints using ComplexHeatmap.
##   5) Add sample-level clinical annotations (Sex, Age group, IDH
##      status) and alteration frequency barplots, and export figures
##      as PDF files.
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

####################################################
## Setup directories
####################################################
dir_input <- 'result/data'
dir_output <- 'result/Fig'

#################################################
## Load data
#################################################
load(file.path(dir_input, 'se_mut_onco_clin.RData'))
mut <- assay(eset)
clin <- colData(eset)

# Fix oncoprint labels
label_map <- c(
  "SNV/Indel"     = "SNV/Indel",
  "Amplification" = "Amplification",
  "Fusion"        = "Fusion",
  "Deletion"      = "Deletion"
)

priority <- c("Fusion", "Amplification", "SNV/Indel", "Deletion")

fix_cell <- function(x) {
  if (is.na(x) || x == "") return("")

  parts <- strsplit(x, ";", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts <- parts[parts != ""]
  parts <- gsub("\\s+", " ", parts)

  # map known labels
  parts <- ifelse(parts %in% names(label_map),
                  unname(label_map[parts]),
                  parts)

  # keep unique + consistent ordering
  parts <- unique(parts)
  ord <- match(parts, priority)
  parts <- parts[order(is.na(ord), ord, parts)]

  paste(parts, collapse = ";")
}

mut_fixed <- mut
mut_fixed[] <- vapply(mut, fix_cell, character(1))

## check
mut_fixed[1:4, 1:6]
table(mut_fixed)

####################################################
## OncoPrint ---> no clinbical metadata
####################################################
combo_levels <- c("Amplification;SNV/Indel", "Fusion;Amplification;SNV/Indel", "Fusion;Amplification")

mut_other <- mut_fixed
mut_other[mut_other %in% combo_levels] <- "Multi-hit"

col <- c(
  "SNV/Indel"     = "#68855CFF",
  "Amplification" = "#A06177FF",
  "Fusion"        = "#A5693CFF",
  "Deletion"      = "#526A83FF", 
  "Multi-hit"     = "#D9AF6BFF"
)

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(2, "pt"), h - unit(2, "pt"),
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  "SNV/Indel" = function(x, y, w, h) {
  grid.rect(x, y, w - unit(2, "pt"), h - unit(2, "pt"),
            gp = gpar(fill = col["SNV/Indel"], col = NA))
},
  "Amplification" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(2, "pt"), h - unit(2, "pt"),
              gp = gpar(fill = col["Amplification"], col = NA))
  },
  "Fusion" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(2, "pt"), h - unit(2, "pt"),
              gp = gpar(fill = col["Fusion"], col = NA))
  },
  "Deletion" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(2, "pt"), h - unit(2, "pt"),
              gp = gpar(fill = col["Deletion"], col = NA))
  },
  "Multi-hit" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(2, "pt"), h - unit(2, "pt"),
              gp = gpar(fill = col["Multi-hit"], col = NA))
  }
)

heatmap_legend_param <- list(
  title  = "Alterations",
  at     = c("SNV/Indel", "Amplification", "Fusion", "Deletion", "Multi-hit"),
  labels = c("SNV/Indel", "Amplification", "Fusion", "Deletion", "Multi-hit")
)

p <- oncoPrint(
  mut_other,
  alter_fun = alter_fun,
  col = col,
  remove_empty_rows = TRUE,
  remove_empty_columns = TRUE,
  heatmap_legend_param = heatmap_legend_param,
    row_names_gp = grid::gpar(fontsize = 7),
  left_annotation = rowAnnotation(
    rbar = anno_oncoprint_barplot(
      axis_param = list(direction = "reverse")
    )
  ),
  right_annotation = NULL
)

pdf(file.path(dir_output, 'fig1.pdf'),  width = 5, height = 7)
p
dev.off()

####################################################
## OncoPrint ---> clinbical metadata (all patients)
####################################################
# sample-level metadata
clin <- as.data.frame(clin)
clin$Age <- ifelse(clin$Age >= 40, '>40', '<40')
clin$IDH.status[clin$IDH.status == 'NA'] <- "Unknown"
clin$IDH.status[clin$IDH.status == 'IDHmut'] <- "Mut"
clin$IDH.status[clin$IDH.status == 'IDHwt'] <- "WT"

# primary location
clin <- clin %>%
  mutate(
    Primary.location = str_squish(Primary.location),
    Primary.location = str_to_title(Primary.location)
  )

clin <- clin %>%
  mutate(
    Location = case_when(
      Primary.location == "Lobar" ~ "Lobar",
      Primary.location == "Cerebellum" ~ "Cerebellum",
      Primary.location == "Thalamic" ~ "Thalamic",
      TRUE ~ "Other"
    )
  )

# Grade
clin <- clin %>%
  mutate(
    Grade = case_when(
      WHO.2021.Grade %in% c(1, 2) ~ "I/II",
      WHO.2021.Grade == 3 ~ "III",
      WHO.2021.Grade == 4 ~ "IV",
      TRUE ~ NA_character_
    )
  )

# Histology
clin$Histo <- coalesce(clin$LGG, clin$HGG)
clin$Histo <- str_squish(clin$Histo)   # remove extra spaces
clin$Histo <- str_to_sentence(clin$Histo)
clin <- clin %>%
  mutate(
    Histo = str_squish(Histo),
    Histo = str_to_sentence(Histo),
    
    Histo = case_when(
      str_detect(Histo, "Glioblastoma") ~ "Glioblastoma",
      
      str_detect(Histo, "Diffuse hemispheric") ~ "Diffuse hemispheric",
      str_detect(Histo, "Diffuse midline") ~ "Diffuse midline",
      str_detect(Histo, "Diffuse high") ~ "Diffuse high-grade",
      str_detect(Histo, "Diffuse low") ~ "Diffuse low-grade",
      
      str_detect(Histo, "Astrocytoma") ~ "Astrocytoma",
      str_detect(Histo, "Oligodendroglioma") ~ "Oligodendroglioma",
      str_detect(Histo, "Glioneuronal") ~ "Glioneuronal",
      str_detect(Histo, "Pilocytic") ~ "Pilocytic astrocytoma",
      str_detect(Histo, "Pleomorphic") ~ "Pleomorphic xanthoastrocytoma",
      str_detect(Histo, "Ependymoma") ~ "Ependymoma",
      
      TRUE ~ Histo
    )
  )

histo_counts <- table(clin$Histo)

# Identify small groups (<10 patients)
small_groups <- names(histo_counts[histo_counts < 10])

clin <- clin %>%
  mutate(
    Histo = ifelse(
      Histo %in% small_groups,
      "Other",
      Histo
    )
  )

anno_df <- data.frame(
  Sex = factor(clin$Sex, levels = c("Male", "Female")),
  Age = factor(clin$Age, levels = c(">40", "<40")),
  IDH = factor(clin$IDH.status, levels = c("WT", "Mut", "Unknown")),
  Location = factor(clin$Location, levels = c('Lobar', 'Cerebellum', 'Thalamic', 'Other')),
  Histo = factor(clin$Histo, levels = c('Glioblastoma', 'Astrocytoma', 'Oligodendroglioma', 'Glioneuronal')),
  Grade = factor(clin$Grade, levels = c('I/II', 'III', 'IV'))
)

# make sure rownames match sample IDs
rownames(anno_df) <- colnames(mut_other)

# colors
anno_col <- list(
  Sex = c(
    Male   = "#4C72B0",
    Female = "#DD8452"
  ),
  IDH = c(
    WT  = "#C44E52",
    Mut = "#55A868",
    Unknown = "#96A5A5FF"
  ),
  Age = c(
    '>40'  = "#08306B",
    '<40' = "#C6DBEF"
  ),
  Grade = c(
    'I/II'  = "#BC8E7DFF",
    'III' = "#FAE093FF",
    'IV' = "#7C7189FF"
  ),
  Location = c(
    'Lobar'  = "#A54B2DFF",
    'Cerebellum' = "#577E2FFF",
    'Thalamic'  = "#B49696FF",
    'Other' = "#96A5A5FF"
  ),
  Histo = c(
    'Glioblastoma'       = "#4A7169FF",
    'Astrocytoma'        = "#735231FF",
    'Oligodendroglioma'  = "#E3CA97FF",
    'Glioneuronal'       = "#99B6BDFF",
    'Other'              = "#96A5A5FF"
  )
)

# create annotation
ha_meta <- HeatmapAnnotation(
  df = anno_df,
  col = anno_col,
  annotation_name_side = "left"
)

ha_bar <- HeatmapAnnotation(
  "Alterations" = anno_oncoprint_barplot()
)

top_anno <- c(ha_bar, ha_meta)
p <- oncoPrint(
  mut_other,
  alter_fun = alter_fun,
  col = col,
  remove_empty_rows = TRUE,
  remove_empty_columns = TRUE,
  heatmap_legend_param = heatmap_legend_param,
  row_names_gp = grid::gpar(fontsize = 7),
  top_annotation = top_anno,
  left_annotation = rowAnnotation(
    rbar = anno_oncoprint_barplot(axis_param = list(direction = "reverse"))
  ),
  right_annotation = NULL
)

pdf(file.path(dir_output, 'fig2.pdf'), width = 8, height = 9)
p
dev.off()


####################################################
## OncoPrint ---> clinbical metadata (IDH WT)
####################################################
# sample-level metadata
clin <- as.data.frame(clin)
clin_wt <- clin[clin$IDH.status == 'IDHwt', ]
clin_wt$Age <- ifelse(clin_wt$Age >= 40, '>40', '<40')

# primary location
clin_wt <- clin_wt %>%
  mutate(
    Primary.location = str_squish(Primary.location),
    Primary.location = str_to_title(Primary.location)
  )

clin_wt <- clin_wt %>%
  mutate(
    Location = case_when(
      Primary.location == "Lobar" ~ "Lobar",
      Primary.location == "Cerebellum" ~ "Cerebellum",
      Primary.location == "Thalamic" ~ "Thalamic",
      TRUE ~ "Other"
    )
  )

# Grade
clin_wt <- clin_wt %>%
  mutate(
    Grade = case_when(
      WHO.2021.Grade %in% c(1, 2) ~ "I/II",
      WHO.2021.Grade == 3 ~ "III",
      WHO.2021.Grade == 4 ~ "IV",
      TRUE ~ NA_character_
    )
  )

# Histology
clin_wt$Histo <- coalesce(clin_wt$LGG, clin_wt$HGG)
clin_wt$Histo <- str_squish(clin_wt$Histo)   # remove extra spaces
clin_wt$Histo <- str_to_sentence(clin_wt$Histo)
clin_wt <- clin_wt %>%
  mutate(
    Histo = str_squish(Histo),
    Histo = str_to_sentence(Histo),
    
    Histo = case_when(
      str_detect(Histo, "Glioblastoma") ~ "Glioblastoma",
      
      str_detect(Histo, "Diffuse hemispheric") ~ "Diffuse hemispheric",
      str_detect(Histo, "Diffuse midline") ~ "Diffuse midline",
      str_detect(Histo, "Diffuse high") ~ "Diffuse high-grade",
      str_detect(Histo, "Diffuse low") ~ "Diffuse low-grade",
      
      str_detect(Histo, "Astrocytoma") ~ "Astrocytoma",
      str_detect(Histo, "Oligodendroglioma") ~ "Oligodendroglioma",
      str_detect(Histo, "Glioneuronal") ~ "Glioneuronal",
      str_detect(Histo, "Pilocytic") ~ "Pilocytic astrocytoma",
      str_detect(Histo, "Pleomorphic") ~ "Pleomorphic xanthoastrocytoma",
      str_detect(Histo, "Ependymoma") ~ "Ependymoma",
      
      TRUE ~ Histo
    )
  )

histo_counts <- table(clin_wt$Histo)

# Identify small groups (<10 patients)
small_groups <- names(histo_counts[histo_counts < 10])

clin_wt <- clin_wt %>%
  mutate(
    Histo = ifelse(
      Histo %in% small_groups,
      "Other",
      Histo
    )
  )

anno_df <- data.frame(
  Sex = factor(clin_wt$Sex, levels = c("Male", "Female")),
  Age = factor(clin_wt$Age, levels = c(">40", "<40")),
  #IDH = factor(clin_wt$IDH.status, levels = c("WT", "Mut", "Unknown")),
  Location = factor(clin_wt$Location, levels = c('Lobar', 'Cerebellum', 'Thalamic', 'Other')),
  Histo = factor(clin_wt$Histo, levels = c('Glioblastoma', 'Astrocytoma', 'Oligodendroglioma', 'Glioneuronal')),
  Grade = factor(clin_wt$Grade, levels = c('I/II', 'III', 'IV'))
)

# make sure rownames match sample IDs
mut_wt <- mut_other[, colnames(mut_other) %in% clin_wt$Study]
rownames(anno_df) <- colnames(mut_wt)

# colors
anno_col <- list(
  Sex = c(
    Male   = "#4C72B0",
    Female = "#DD8452"
  ),
 # IDH = c(
 #   WT  = "#C44E52",
 #   Mut = "#55A868",
 #   Unknown = "#96A5A5FF"
 # ),
  Age = c(
    '>40'  = "#08306B",
    '<40' = "#C6DBEF"
  ),
  Grade = c(
    'I/II'  = "#BC8E7DFF",
    'III' = "#FAE093FF",
    'IV' = "#7C7189FF"
  ),
  Location = c(
    'Lobar'  = "#A54B2DFF",
    'Cerebellum' = "#577E2FFF",
    'Thalamic'  = "#B49696FF",
    'Other' = "#96A5A5FF"
  ),
  Histo = c(
    'Glioblastoma'       = "#4A7169FF",
    'Astrocytoma'        = "#735231FF",
    'Oligodendroglioma'  = "#E3CA97FF",
    'Glioneuronal'       = "#99B6BDFF",
    'Other'              = "#96A5A5FF"
  )
)

# create annotation
ha_meta <- HeatmapAnnotation(
  df = anno_df,
  col = anno_col,
  annotation_name_side = "left"
)

ha_bar <- HeatmapAnnotation(
  "Alterations" = anno_oncoprint_barplot()
)

top_anno <- c(ha_bar, ha_meta)
p <- oncoPrint(
  mut_wt,
  alter_fun = alter_fun,
  col = col,
  remove_empty_rows = TRUE,
  remove_empty_columns = TRUE,
  heatmap_legend_param = heatmap_legend_param,
  row_names_gp = grid::gpar(fontsize = 7),
  top_annotation = top_anno,
  left_annotation = rowAnnotation(
    rbar = anno_oncoprint_barplot(axis_param = list(direction = "reverse"))
  ),
  right_annotation = NULL
)

pdf(file.path(dir_output, 'fig3_wt.pdf'), width = 7, height = 9)
p
dev.off()


####################################################
## OncoPrint ---> clinbical metadata (IDH MUT)
####################################################
# sample-level metadata
clin <- as.data.frame(clin)
clin_mut <- clin[clin$IDH.status == 'IDHmut', ]
clin_mut$Age <- ifelse(clin_mut$Age >= 40, '>40', '<40')

# primary location
clin_mut <- clin_mut %>%
  mutate(
    Primary.location = str_squish(Primary.location),
    Primary.location = str_to_title(Primary.location)
  )

clin_mut <- clin_mut %>%
  mutate(
    Location = case_when(
      Primary.location == "Lobar" ~ "Lobar",
      Primary.location == "Cerebellum" ~ "Cerebellum",
      TRUE ~ Location
    )
  )

# Grade
clin_mut <- clin_mut %>%
  mutate(
    Grade = case_when(
      WHO.2021.Grade %in% c(1, 2) ~ "I/II",
      WHO.2021.Grade == 3 ~ "III",
      WHO.2021.Grade == 4 ~ "IV",
      TRUE ~ NA_character_
    )
  )

# Histology
clin_mut$Histo <- coalesce(clin_mut$LGG, clin_mut$HGG)
clin_mut$Histo <- str_squish(clin_mut$Histo)   # remove extra spaces
clin_mut$Histo <- str_to_sentence(clin_mut$Histo)
clin_mut <- clin_mut %>%
  mutate(
    Histo = str_squish(Histo),
    Histo = str_to_sentence(Histo),
    Histo = case_when(
      str_detect(Histo, "Astrocytoma") ~ "Astrocytoma",
      str_detect(Histo, "Oligodendroglioma") ~ "Oligodendroglioma",
      TRUE ~ Histo
    )
  )


anno_df <- data.frame(
  Sex = factor(clin_mut$Sex, levels = c("Male", "Female")),
  Age = factor(clin_mut$Age, levels = c(">40", "<40")),
  #IDH = factor(clin_wt$IDH.status, levels = c("WT", "Mut", "Unknown")),
  Location = factor(clin_mut$Location, levels = c('Lobar', 'Cerebellum')),
  Histo = factor(clin_mut$Histo, levels = c('Astrocytoma', 'Oligodendroglioma')),
  Grade = factor(clin_mut$Grade, levels = c('I/II', 'III', 'IV'))
)

# make sure rownames match sample IDs
mut_mut <- mut_other[, colnames(mut_other) %in% clin_mut$Study]
rownames(anno_df) <- colnames(mut_mut)

# colors
anno_col <- list(
  Sex = c(
    Male   = "#4C72B0",
    Female = "#DD8452"
  ),
 # IDH = c(
 #   WT  = "#C44E52",
 #   Mut = "#55A868",
 #   Unknown = "#96A5A5FF"
 # ),
  Age = c(
    '>40'  = "#08306B",
    '<40' = "#C6DBEF"
  ),
  Grade = c(
    'I/II'  = "#BC8E7DFF",
    'III' = "#FAE093FF",
    'IV' = "#7C7189FF"
  ),
  Location = c(
    'Lobar'  = "#A54B2DFF",
    'Cerebellum' = "#577E2FFF"
  ),
  Histo = c(
    'Astrocytoma'        = "#735231FF",
    'Oligodendroglioma'  = "#E3CA97FF"
  )
)

# create annotation
ha_meta <- HeatmapAnnotation(
  df = anno_df,
  col = anno_col,
  annotation_name_side = "left"
)

ha_bar <- HeatmapAnnotation(
  "Alterations" = anno_oncoprint_barplot()
)

top_anno <- c(ha_bar, ha_meta)
p <- oncoPrint(
  mut_mut,
  alter_fun = alter_fun,
  col = col,
  remove_empty_rows = TRUE,
  remove_empty_columns = TRUE,
  heatmap_legend_param = heatmap_legend_param,
  row_names_gp = grid::gpar(fontsize = 7),
  top_annotation = top_anno,
  left_annotation = rowAnnotation(
    rbar = anno_oncoprint_barplot(axis_param = list(direction = "reverse"))
  ),
  right_annotation = NULL
)

pdf(file.path(dir_output, 'fig4_mut.pdf'), width = 6, height = 7)
p
dev.off()

###################################################
## Bar plots: top mutations and sex/age
###################################################
load(file.path(dir_input, 'se_mut_bin_clin.RData'))
mut <- assay(eset)

## choose top mutations
mut_freq <- rowMeans(mut == 1)

genes <- names(mut_freq[mut_freq >= 0.10])
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
      select(Sample, Age, Sex),
    by = "Sample"
  )

###########################################################
## Fishe exact test --- Age/sex and mutations
###########################################################
age_assoc <- mut_long %>%
  group_by(Gene) %>%
  summarise(
    pval = fisher.test(table(Mutation, Age))$p.value
  )

age_assoc$fdr <- p.adjust(age_assoc$pval, method = "BH")
write.csv(age_assoc, file = file.path(dir_output, 'mut_age.csv'))

## Fishe exact test --- Sex and mutations
sex_assoc <- mut_long %>%
  group_by(Gene) %>%
  summarise(
    pval = fisher.test(table(Mutation, Sex))$p.value
  )

sex_assoc$fdr <- p.adjust(sex_assoc$pval, method = "BH")
write.csv(sex_assoc, file = file.path(dir_output, 'mut_sex.csv'))

## bar plot age
gene_order <- mut_long %>%
  group_by(Gene) %>%
  summarise(overall = mean(Mutation)) %>%
  arrange(desc(overall)) %>%
  pull(Gene)

plot_age <- mut_long %>%
  group_by(Gene, Age) %>%
  summarise(freq = mean(Mutation) * 100, .groups = "drop")

plot_age$Gene <- factor(plot_age$Gene, levels = gene_order)

pdf(file.path(dir_output, 'mut_age.pdf'), width = 7, height = 6)

ggplot(plot_age, aes(x = Gene, y = freq, fill = Age)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c(">40" = "#08306B",
               "<40" = "#C6DBEF")
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Mutation frequency (%)", x = "", 
       title = " ")

dev.off()

## bar plot sex
plot_sex <- mut_long %>%
  group_by(Gene, Sex) %>%
  summarise(freq = mean(Mutation) * 100, .groups = "drop")

plot_sex$Gene <- factor(plot_sex$Gene, levels = gene_order)

pdf(file.path(dir_output, 'mut_sex.pdf'), width = 7, height = 6)

ggplot(plot_sex, aes(x = Gene, y = freq, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c("Female" = "#DD8452",
               "Male"   = "#4C72B0")
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Mutation frequency (%)", x = "", 
       title = " ")

dev.off()

