#####################################################
## Script: CNS NGS OncoPrint (ComplexHeatmap)
##
## Purpose:
##   Generate publication-quality OncoPrint visualizations
##   from a curated CNS NGS mutation event matrix stored
##   in a SummarizedExperiment object.
##
##   The script:
##     • Standardizes alteration labels
##     • Enforces consistent within-cell alteration priority
##     • Preserves multi-hit events (no collapsing)
##     • Integrates clinical metadata for annotation
##     • Produces OncoPrints for:
##         - All patients
##         - IDH wild-type subset
##         - IDH mutant subset
##
##
## Input:
##   - result/data/se_mut_onco_clin.RData
##       SummarizedExperiment containing:
##         assay:
##           gene × sample character matrix
##           ("" or semicolon-delimited alteration types)
##         colData:
##           harmonized clinical metadata
##
##
## Output (4 figures):
##
##   - result/oncoprint/fig1.pdf
##       OncoPrint (mutations only; no clinical annotations)
##
##   - result/oncoprint/fig2.pdf
##       OncoPrint with clinical annotations (all patients)
##
##   - result/oncoprint/fig3_wt.pdf
##       OncoPrint with annotations (IDH wild-type subset)
##
##   - result/oncoprint/fig4_mut.pdf
##       OncoPrint with annotations (IDH mutant subset)
##
##
## Key Processing Steps:
##
##   1) Mutation Label Harmonization
##        - Standardize alteration labels to:
##            SNV/Indel
##            Amplification
##            Fusion
##            Deletion
##
##        - Enforce within-cell priority ordering:
##            Fusion > Amplification > SNV/Indel > Deletion
##
##        - Preserve multiple alteration types per gene/sample
##          as semicolon-separated entries (multi-hit retained).
##
##   2) OncoPrint Construction (ComplexHeatmap)
##        - Custom graphical rendering for each alteration type
##        - Row-level mutation frequency barplots
##        - Removal of empty rows and columns
##
##   3) Clinical Annotation Harmonization
##        - Age dichotomized (<40 vs ≥40)
##        - IDH status recoded (WT / Mut / Unknown)
##        - Tumor location grouped into major categories
##        - Histology harmonized and rare groups collapsed
##        - WHO 2021 grade formatted (I–IV)
##
##   4) Subset Analyses
##        - Independent OncoPrints generated for:
##            • IDH wild-type patients
##            • IDH mutant patients
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
dir_output <- 'result/oncoprint'

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
#combo_levels <- c("Amplification;SNV/Indel", "Fusion;Amplification;SNV/Indel", "Fusion;Amplification")

#mut_other <- mut_fixed
#mut_other[mut_other %in% combo_levels] <- "Multi-hit"

col <- c(
  "SNV/Indel"     = "#68855CFF",
  "Amplification" = "#A06177FF",
  "Fusion"        = "#A5693CFF",
  "Deletion"      = "#526A83FF"
 # "Multi-hit"     = "#D9AF6BFF"
)

alter_fun <- list(
  
  background = function(x, y, w, h) {
  grid.rect(x, y, w, h,
            gp = gpar(fill =  "#f5f5f5", col = NA))
  },

  "SNV/Indel" = function(x, y, w, h) {
    grid.rect(x, y + h*0.30, w*0.95, h*0.25,
              gp = gpar(fill = col["SNV/Indel"], col = NA))
  },

  "Amplification" = function(x, y, w, h) {
    grid.rect(x, y + h*0.05, w*0.95, h*0.25,
              gp = gpar(fill = col["Amplification"], col = NA))
  },

  "Fusion" = function(x, y, w, h) {
    grid.rect(x, y - h*0.20, w*0.95, h*0.25,
              gp = gpar(fill = col["Fusion"], col = NA))
  },

  "Deletion" = function(x, y, w, h) {
    grid.rect(x, y - h*0.45, w*0.95, h*0.25,
              gp = gpar(fill = col["Deletion"], col = NA))
  }
)

left_annotation = rowAnnotation(
  rbar = anno_oncoprint_barplot(
    axis_param = list(direction = "reverse"),
    border = FALSE
  )
)

heatmap_legend_param <- list(
  title = "Alterations",
  at = names(col),
  labels = names(col),
  nrow = 4,                             
  grid_width = unit(3, "mm"),             
  grid_height = unit(3, "mm"),
  title_gp = gpar(fontsize = 10, face='bold'),          
  labels_gp = gpar(fontsize = 8)          
)

p <- oncoPrint(
  mut_fixed,
  get_type = function(x) {
  if(x == "") return(NULL)
  strsplit(x, ";", fixed = TRUE)[[1]]
},
  alter_fun = alter_fun,
  col = col,
  remove_empty_rows = TRUE,
  remove_empty_columns = TRUE,
  heatmap_legend_param = heatmap_legend_param,
  row_names_gp = gpar(fontsize = 8),
  left_annotation = left_annotation,
  right_annotation = NULL
)

pdf(file.path(dir_output, 'fig1.pdf'),  width = 7, height = 9)

draw(
  p,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

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
      WHO.2021.Grade == 1 ~ "I",
      WHO.2021.Grade == 2 ~ "II",
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
  Grade = factor(clin$Grade, levels = c('I', 'II', 'III', 'IV'))
)

# make sure rownames match sample IDs
rownames(anno_df) <- colnames(mut_fixed)

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
    'I'  = "#A8C3A0FF",
    'II'  = "#BC8E7DFF",
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
top_anno <- HeatmapAnnotation(
  
  # ---- alteration barplot ----
  Alterations = anno_oncoprint_barplot(
    border = FALSE,
    height = unit(1.2, "cm")
  ),

  # ---- metadata ----
  df = anno_df,
  col = anno_col,
  simple_anno_size = unit(3, "mm"), 
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 7) 
)

heatmap_legend_param <- list(
  title = "Alterations",
  at = names(col),
  labels = names(col),
  nrow = 4,                             
  grid_width = unit(3, "mm"),             
  grid_height = unit(3, "mm"),
  title_gp = gpar(fontsize = 10, face='bold'),          
  labels_gp = gpar(fontsize = 8)          
)

p <- oncoPrint(
  mut_fixed,
  get_type = function(x) {
  if(x == "") return(NULL)
  strsplit(x, ";", fixed = TRUE)[[1]]
},
  alter_fun = alter_fun,
  col = col,
  remove_empty_rows = TRUE,
  remove_empty_columns = TRUE,
  heatmap_legend_param = heatmap_legend_param,
  row_names_gp = grid::gpar(fontsize = 8),
  top_annotation = top_anno,
  left_annotation = rowAnnotation(
    rbar = anno_oncoprint_barplot(axis_param = list(direction = "reverse"))
  ),
  right_annotation = NULL
)

pdf(file.path(dir_output, 'fig2.pdf'), width = 7, height = 8)

draw(
  p,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

####################################################
## OncoPrint ---> clinbical metadata (IDH WT)
####################################################
# sample-level metadata
clin <- as.data.frame(clin)
clin_wt <- clin[clin$IDH.status == 'WT', ]

anno_df <- data.frame(
  Sex = factor(clin_wt$Sex, levels = c("Male", "Female")),
  Age = factor(clin_wt$Age, levels = c(">40", "<40")),
  #IDH = factor(clin_wt$IDH.status, levels = c("WT", "Mut", "Unknown")),
  Location = factor(clin_wt$Location, levels = c('Lobar', 'Cerebellum', 'Thalamic', 'Other')),
  Histo = factor(clin_wt$Histo, levels = c('Glioblastoma', 'Astrocytoma', 'Oligodendroglioma', 'Glioneuronal')),
  Grade = factor(clin_wt$Grade, levels = c('I', 'II', 'III', 'IV'))
)

# make sure rownames match sample IDs
mut_wt <- mut_fixed[, colnames(mut_fixed) %in% clin_wt$Study]
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
    'I'  = "#A8C3A0FF",
    'II'  = "#BC8E7DFF",
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
top_anno <- HeatmapAnnotation(
  
  # ---- alteration barplot ----
  Alterations = anno_oncoprint_barplot(
    border = FALSE,
    height = unit(1.2, "cm")
  ),

  # ---- metadata ----
  df = anno_df,
  col = anno_col,
  simple_anno_size = unit(3, "mm"), 
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 7) 
)

heatmap_legend_param <- list(
  title = "Alterations",
  at = names(col),
  labels = names(col),
  nrow = 4,                             
  grid_width = unit(3, "mm"),             
  grid_height = unit(3, "mm"),
  title_gp = gpar(fontsize = 10, face='bold'),          
  labels_gp = gpar(fontsize = 8)          
)

p <- oncoPrint(
  mut_wt,
  get_type = function(x) strsplit(x, ";")[[1]],
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

pdf(file.path(dir_output, 'fig3_wt.pdf'), width = 7, height = 8)

draw(
  p,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

####################################################
## OncoPrint ---> clinbical metadata (IDH MUT)
####################################################
# sample-level metadata
clin <- as.data.frame(clin)
clin_mut <- clin[clin$IDH.status == 'Mut', ]

anno_df <- data.frame(
  Sex = factor(clin_mut$Sex, levels = c("Male", "Female")),
  Age = factor(clin_mut$Age, levels = c(">40", "<40")),
  #IDH = factor(clin_wt$IDH.status, levels = c("WT", "Mut", "Unknown")),
  Location = factor(clin_mut$Location, levels = c('Lobar', 'Cerebellum')),
  Histo = factor(clin_mut$Histo, levels = c('Astrocytoma', 'Oligodendroglioma')),
  Grade = factor(clin_mut$Grade, levels = c('I', 'II', 'III', 'IV'))
)

# make sure rownames match sample IDs
mut_mut <- mut_fixed[, colnames(mut_fixed) %in% clin_mut$Study]
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
    'I'  = "#A8C3A0FF",
    'II'  = "#BC8E7DFF",
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
top_anno <- HeatmapAnnotation(
  
  # ---- alteration barplot ----
  Alterations = anno_oncoprint_barplot(
    border = FALSE,
    height = unit(1.2, "cm")
  ),

  # ---- metadata ----
  df = anno_df,
  col = anno_col,
  simple_anno_size = unit(3, "mm"), 
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 7) 
)

heatmap_legend_param <- list(
  title = "Alterations",
  at = names(col),
  labels = names(col),
  nrow = 4,                             
  grid_width = unit(3, "mm"),             
  grid_height = unit(3, "mm"),
  title_gp = gpar(fontsize = 10, face='bold'),          
  labels_gp = gpar(fontsize = 8)          
)

p <- oncoPrint(
  mut_mut,
  get_type = function(x) strsplit(x, ";")[[1]],
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

pdf(file.path(dir_output, 'fig4_mut.pdf'), width = 6, height = 5)

draw(
  p,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

###################################################
## Bar plots: top mutations and sex/age
###################################################
load(file.path(dir_input, 'se_mut_bin_clin.RData'))
mut <- assay(eset)
clin <- as.data.frame(colData(eset))
clin$Age <- ifelse(clin$Age >= 40, '>40', '<40')
clin$IDH.status[clin$IDH.status == 'NA'] <- "Unknown"
clin$IDH.status[clin$IDH.status == 'IDHmut'] <- "Mut"
clin$IDH.status[clin$IDH.status == 'IDHwt'] <- "WT"

## choose top mutations
mut_freq <- rowMeans(mut == 1)

genes <- names(mut_freq[mut_freq >= 0.1])
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
      select(Sample, Age, Sex, IDH.status),
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
write.csv(age_assoc, file = file.path(dir_output, 'mut_age.csv'))

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

pdf(file.path(dir_output, 'mut_age.pdf'), width = 7, height = 6)

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
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Mutation frequency (%)", x = "", 
       title = " ")

dev.off()

## Step 2 --- Mutant patients
mut_long_mut <- mut_long[mut_long$IDH.status == 'Mut', ]
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

pdf(file.path(dir_output, 'mut_age_mut.pdf'), width = 7, height = 6)

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
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Mutation frequency (%)", x = "", 
       title = " ")

dev.off()

## Step 3 --- Wild type patients
mut_long_wt <- mut_long[mut_long$IDH.status == 'WT', ]
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

pdf(file.path(dir_output, 'mut_age_wt.pdf'), width = 7, height = 6)

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
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
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

## Step 2 --- mutant patients
mut_long_mut <- mut_long[mut_long$IDH.status == 'Mut', ]
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

pdf(file.path(dir_output, 'mut_sex_mut.pdf'), width = 7, height = 6)

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

## Step 2 --- mutant patients
mut_long_wt <- mut_long[mut_long$IDH.status == 'WT', ]
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

pdf(file.path(dir_output, 'mut_sex_wt.pdf'), width = 7, height = 6)

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

