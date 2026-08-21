##-------------------------------------------------------------------
## Script: CNS NGS OncoPrint Visualization
## Purpose:
##   To generate publication-quality OncoPrint
##   visualizations for CNS tumour patients using
##   curated NGS and clinical data.
##
##   Annotated OncoPrint visualizations are generated
##   for the full cohort and IDH-stratified
##   subgroups.
##
## Input:
##   - result/data/mae_mut_clin.RData
##       MultiAssayExperiment object containing:
##         • OncoPrint mutation matrix
##         • harmonized clinical metadata
##
## Outputs:
##   1) Annotated OncoPrint visualizations
##   2) Publication-quality PDF and JPEG figures
##
## Processing Overview:
##   1) Load mutation and clinical data
##   2) Standardize mutation alteration labels
##   3) Harmonize clinical variables including:
##        - Age group
##        - Sex
##        - IDH status
##        - Histology
##        - WHO grade
##   4) Generate annotated OncoPrints for:
##        - Full cohort
##        - IDH wild-type subgroup
##        - IDH mutant subgroup
##        - IDH wild-type glioblastoma subgroup
##        - IDH wild-type non-glioblastoma subgroup
##   5) Filter low-frequency alterations for
##      subgroup-specific visualizations
##   6) Export publication-quality PDF and
##      JPEG figures
##
## Notes:
##   - Analyses are based on categorical
##     gene-level alteration matrices.
##   - Alteration classes include:
##        - SNV/Indel
##        - Amplification
##        - Fusion
##        - Deletion
##   - Clinical annotations include demographic,
##     histological, and molecular subgroup
##     information.
##   - The IDH wild-type cohort is further
##     divided into glioblastoma and
##     non-glioblastoma subgroups.
##   - Visualizations are generated using the
##     ComplexHeatmap framework.
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

#################################################
## functions
#################################################

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

####################################################
## Setup directories
####################################################
dir_input <- 'result/data'
dir_output <- 'result/oncoprint'

#################################################
## Load MultiAssayExperiment data
#################################################
load(file.path(dir_input, 'mae_mut_clin.RData'))
mut <- assay(mae[['mut_oncoprint']])
clin <- colData(mae[['mut_oncoprint']])

# Fix oncoprint labels
label_map <- c(
  "SNV/Indel"     = "SNV/Indel",
  "Amplification" = "Amplification",
  "Fusion"        = "Fusion",
  "Deletion"      = "Deletion"
)

priority <- c("Fusion", "Amplification", "SNV/Indel", "Deletion")

mut_fixed <- mut
mut_fixed[] <- vapply(mut, fix_cell, character(1))

## check
mut_fixed[1:4, 1:6]
table(mut_fixed)

####################################################
## OncoPrint settings
####################################################
col <- c(
  "SNV/Indel"     = "#A06177FF",
  "Amplification" = "#68855CFF",
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

annotation_legend_param <- list(
  Sex = list(
    title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6)
  ),
  
  Age = list(
    title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6)
  ),
  
  IDH = list(
    title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6)
  ),
  
  Histo = list(
    title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 5.5)
  ),
  
  Grade = list(
    title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6)
  )
)

####################################################
## OncoPrint ---> clinbical metadata (all patients)
####################################################
# sample-level metadata
clin <- as.data.frame(clin)
clin$Age <- ifelse(clin$Age >= 40, '>40', '<40')

# Grade
clin <- clin %>%
  mutate(
    Grade = case_when(
      WHO.2021.Grade == 1 ~ "1",
      WHO.2021.Grade == 2 ~ "2",
      WHO.2021.Grade == 3 ~ "3",
      WHO.2021.Grade == 4 ~ "4",
      TRUE ~ NA_character_
    )
  )

# Histology
clin$Histo <- clin$histo
clin <- clin %>%
  mutate(
    Histo = str_squish(Histo),
    Histo = str_to_lower(Histo),   # normalize first
    
    Histo = case_when(
      
      str_detect(Histo, "glioblastoma") ~ "Glioblastoma",
      str_detect(Histo, "oligodendroglioma") ~ "Oligodendroglioma",
      str_detect(Histo, "astrocytoma") ~ "Astrocytoma",
      str_detect(Histo, "ependymoma") ~ "Ependymoma",
      str_detect(Histo, "glioneuronal") ~ "Glioneuronal tumor",
      str_detect(Histo, "circumscribed glioma") ~ "Circumscribed glioma",
      str_detect(Histo, "pediatric-type high.?grade glioma") ~ "Pediatric-type HGG",
      str_detect(Histo, "pediatric-type low.?grade glioma") ~ "Pediatric-type LGG",
      TRUE ~ str_to_sentence(Histo)  # keep original but clean formatting
    )
  )

anno_df <- data.frame(
  Sex = factor(clin$Sex, levels = c("Male", "Female")),
  Age = factor(clin$Age, levels = c(">40", "<40")),
  IDH = factor(clin$IDH_status, levels = c("WT", "Mut")),
  Histo = factor(clin$Histo, levels = c('Glioblastoma', 'Astrocytoma', 'Circumscribed glioma',
                                        'Glioneuronal tumor', 'Oligodendroglioma',  'Pediatric-type HGG',
                                        'Pediatric-type LGG', 'Ependymoma')),
  Grade = factor(clin$Grade, levels = c('1', '2', '3', '4'))
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
    Mut = "#55A868"
  ),
  Age = c(
    '>40'  = "#08306B",
    '<40' = "#C6DBEF"
  ),
  Grade = c(
    '1'  = "#A8C3A0FF",
    '2'  = "#BC8E7DFF",
    '3' = "#FAE093FF",
    '4' = "#7C7189FF"
  ),
  Histo = c(
    'Glioblastoma'          = "#4A7169FF",
    'Astrocytoma'           = "#735231FF",
    'Circumscribed glioma'  = "#E76254FF", 
    'Glioneuronal tumor'    = "#99B6BDFF",
    'Oligodendroglioma'     = "#E3CA97FF",
    'Pediatric-type HGG'    = "#B49696FF",
    'Pediatric-type LGG'    = "#AAC197FF",
    'Ependymoma'            = "#96A5A5FF" 
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
  annotation_name_gp = gpar(fontsize = 7), 

  annotation_legend_param = annotation_legend_param
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

pdf(file.path(dir_output, 'fig1.pdf'), width = 7, height = 8)

draw(
  p,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

jpeg(filename = file.path(dir_output, "fig1.jpeg"), width = 7, height = 8, units = "in", 
    res = 300, quality = 100)

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
clin_wt <- clin[clin$IDH_status == 'WT', ]

anno_df <- data.frame(
  Sex = factor(clin_wt$Sex, levels = c("Male", "Female")),
  Age = factor(clin_wt$Age, levels = c(">40", "<40")),
  Histo = factor(clin_wt$Histo, c('Glioblastoma', 'Circumscribed glioma',
                                  'Glioneuronal tumor', 'Pediatric-type HGG', 
                                  'Pediatric-type LGG', 'Ependymoma')),
  Grade = factor(clin_wt$Grade, levels = c('1', '2', '3', '4'))
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
  Age = c(
    '>40'  = "#08306B",
    '<40' = "#C6DBEF"
  ),
  Grade = c(
    '1'  = "#A8C3A0FF",
    '2'  = "#BC8E7DFF",
    '3' = "#FAE093FF",
    '4' = "#7C7189FF"
  ),
  Histo = c(
    'Glioblastoma'          = "#4A7169FF",
    'Circumscribed glioma'  = "#E76254FF", 
    'Glioneuronal tumor'    = "#99B6BDFF",
    'Pediatric-type HGG'    = "#B49696FF",
    'Pediatric-type LGG'    = "#AAC197FF",
    'Ependymoma'            = "#96A5A5FF"
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

freq <- rowMeans(mut_wt != '')
keep <- freq[freq > 0.05]
mut_wt_filtered <- mut_wt[rownames(mut_wt) %in% names(keep), ]

p <- oncoPrint(
  mut_wt_filtered,
  get_type = function(x) strsplit(x, ";")[[1]],
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

pdf(file.path(dir_output, 'fig2_wt.pdf'), width = 6.5, height = 5)

draw(
  p,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

jpeg(filename = file.path(dir_output, "fig2_wt.jpeg"), width = 6.5, height = 5, units = "in", 
    res = 300, quality = 100)

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
clin_mut <- clin[clin$IDH_status == 'Mut', ]

anno_df <- data.frame(
  Sex = factor(clin_mut$Sex, levels = c("Male", "Female")),
  Age = factor(clin_mut$Age, levels = c(">40", "<40")),
  Histo = factor(clin_mut$Histo, levels = c('Astrocytoma', 'Oligodendroglioma')),
  Grade = factor(clin_mut$Grade, levels = c('2', '3', '4'))
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
  Age = c(
    '>40'  = "#08306B",
    '<40' = "#C6DBEF"
  ),
  Grade = c(
    '2'  = "#BC8E7DFF",
    '3' = "#FAE093FF",
    '4' = "#7C7189FF"
  ),
  Histo = c(
    'Astrocytoma'           = "#735231FF",
    'Oligodendroglioma'     = "#E3CA97FF"
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

freq <- rowMeans(mut_mut != '')
keep <- freq[freq > 0.05]
mut_mut_filtered <- mut_mut[rownames(mut_mut) %in% names(keep), ]

p <- oncoPrint(
  mut_mut_filtered,
  get_type = function(x) strsplit(x, ";")[[1]],
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

pdf(file.path(dir_output, 'fig3_mut.pdf'), width = 6.5, height = 5)

draw(
  p,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

jpeg(filename = file.path(dir_output, "fig3_mut.jpeg"), width = 6.5, height = 5, units = "in", 
    res = 300, quality = 100)

draw(
  p,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

####################################################
## OncoPrint ---> clinbical metadata (IDH WT-GBM)
####################################################
# sample-level metadata
clin <- as.data.frame(clin)
clin_wt_gbm <- clin[clin$IDH_status == 'WT' & clin$Histo == 'Glioblastoma', ]

anno_df <- data.frame(
  Sex = factor(clin_wt_gbm$Sex, levels = c("Male", "Female")),
  Age = factor(clin_wt_gbm$Age, levels = c(">40", "<40")),
  Grade = factor(clin_wt_gbm$Grade, levels = c('1', '2', '4'))
)

# make sure rownames match sample IDs
mut_wt_gbm <- mut_fixed[, colnames(mut_fixed) %in% clin_wt_gbm$Study]
rownames(anno_df) <- colnames(mut_wt_gbm)

# colors
anno_col <- list(
  Sex = c(
    Male   = "#4C72B0",
    Female = "#DD8452"
  ),
  Age = c(
    '>40'  = "#08306B",
    '<40' = "#C6DBEF"
  ),
  Grade = c(
    '1'  = "#A8C3A0FF",
    '2'  = "#BC8E7DFF",
    '4' = "#7C7189FF"
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

freq <- rowMeans(mut_wt_gbm != '')
keep <- freq[freq > 0.05]
mut_wt_gbm_filtered <- mut_wt_gbm[rownames(mut_wt_gbm) %in% names(keep), ]

p <- oncoPrint(
  mut_wt_gbm_filtered,
  get_type = function(x) strsplit(x, ";")[[1]],
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

pdf(file.path(dir_output, 'fig4_wt_gbm.pdf'), width = 6.5, height = 5)

draw(
  p,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

jpeg(filename = file.path(dir_output, "fig4_wt_gbm.jpeg"), width = 6.5, height = 5, units = "in", 
    res = 300, quality = 100)

draw(
  p,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

####################################################
## OncoPrint ---> clinbical metadata (IDH WT-nonGBM)
####################################################
# sample-level metadata
clin <- as.data.frame(clin)
clin_wt_nongbm <- clin[clin$IDH_status == 'WT' & clin$Histo != 'Glioblastoma', ]

anno_df <- data.frame(
  Sex = factor(clin_wt_nongbm$Sex, levels = c("Male", "Female")),
  Age = factor(clin_wt_nongbm$Age, levels = c(">40", "<40")),
  Histo = factor(clin_wt_nongbm$Histo, c('Circumscribed glioma', 'Glioneuronal tumor', 
                                         'Pediatric-type HGG', 'Pediatric-type LGG', 
                                         'Ependymoma')),
  Grade = factor(clin_wt_nongbm$Grade, levels = c('1', '2', '3', '4'))
)

# make sure rownames match sample IDs
mut_wt_nongbm <- mut_fixed[, colnames(mut_fixed) %in% clin_wt_nongbm$Study]
rownames(anno_df) <- colnames(mut_wt_nongbm)

# colors
anno_col <- list(
  Sex = c(
    Male   = "#4C72B0",
    Female = "#DD8452"
  ),
  Age = c(
    '>40'  = "#08306B",
    '<40' = "#C6DBEF"
  ),
  Grade = c(
    '1'  = "#A8C3A0FF",
    '2'  = "#BC8E7DFF",
    '3' = "#FAE093FF",
    '4' = "#7C7189FF"
  ),
  Histo = c(
    'Circumscribed glioma'  = "#E76254FF", 
    'Glioneuronal tumor'    = "#99B6BDFF",
    'Pediatric-type HGG'    = "#B49696FF",
    'Pediatric-type LGG'    = "#AAC197FF",
    'Ependymoma'            = "#96A5A5FF"
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

freq <- rowMeans(mut_wt_nongbm != '')
keep <- freq[freq > 0.05]
mut_wt_nongbm_filtered <- mut_wt_nongbm[rownames(mut_wt_nongbm) %in% names(keep), ]

p <- oncoPrint(
  mut_wt_nongbm_filtered,
  get_type = function(x) strsplit(x, ";")[[1]],
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

pdf(file.path(dir_output, 'fig5_wt_nongbm.pdf'), width = 6.5, height = 5)

draw(
  p,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

jpeg(filename = file.path(dir_output, "fig5_wt_nongbm.jpeg"), width = 6.5, height = 5, units = "in", 
    res = 300, quality = 100)

draw(
  p,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()