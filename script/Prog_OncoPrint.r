##-------------------------------------------------------------------
## Script: CNS NGS OncoPrint Visualization (ComplexHeatmap)
##
## Purpose:
##   Generate publication-quality OncoPrint visualizations
##   from a curated CNS NGS mutation event matrix stored
##   in a MultiAssayExperiment object.
##
##   Standardize mutation alteration labels and enforce a
##   consistent within-cell priority ordering to ensure
##   reproducible visualization of multi-event mutations.
##
##   Integrate harmonized clinical metadata and generate
##   annotated OncoPrints for the full cohort and key
##   molecular subgroups.
##
##
## Input:
##   - result/data/mae_mut_clin.RData
##       MultiAssayExperiment object containing:
##         - mut_oncoprint assay:
##             gene × patient character matrix
##             ("" or semicolon-delimited alteration types)
##         - colData:
##             curated clinical metadata
##
##
## Output (4 figures):
##
##   1) result/oncoprint/fig1.pdf
##        OncoPrint (mutation events only; no annotations)
##
##   2) result/oncoprint/fig2.pdf
##        OncoPrint with clinical annotations (all patients)
##
##   3) result/oncoprint/fig3_wt.pdf
##        OncoPrint with annotations (IDH wild-type subset)
##
##   4) result/oncoprint/fig4_mut.pdf
##        OncoPrint with annotations (IDH mutant subset)
##
##
## Processing Overview:
##
##   1) Data Loading
##        The MultiAssayExperiment object is loaded and the
##        oncoprint assay (mutation matrix) and associated
##        clinical metadata are extracted.
##
##   2) Mutation Label Harmonization
##        Mutation labels are standardized to:
##           - SNV/Indel
##           - Amplification
##           - Fusion
##           - Deletion
##
##        Within each gene–patient cell, multiple alteration
##        types are:
##           - cleaned and deduplicated
##           - ordered using a fixed priority:
##                Fusion > Amplification > SNV/Indel > Deletion
##
##        Multi-hit events are preserved as semicolon-
##        separated entries.
##
##   3) OncoPrint Construction (ComplexHeatmap)
##        OncoPrints are generated using custom graphical
##        functions for each alteration type, enabling
##        stacked visual representation within each cell.
##
##        Features include:
##           - Removal of empty genes and samples
##           - Row-level mutation frequency barplots
##           - Consistent color mapping across alteration types
##
##   4) Clinical Annotation Processing
##        Clinical metadata are formatted for visualization:
##           - Age dichotomized (<40 vs ≥40)
##           - IDH status encoded as WT / Mut
##           - Histology harmonized and grouped
##           - WHO 2021 grade converted to I–IV
##
##        Annotation tracks are constructed and aligned
##        with the mutation matrix.
##
##   5) Annotated OncoPrint Generation
##        A full-cohort OncoPrint is generated with clinical
##        annotations displayed as top annotation tracks.
##
##   6) Subset Analyses
##        The dataset is stratified by IDH status:
##           - IDH wild-type patients
##           - IDH mutant patients
##
##        Independent OncoPrints are generated for each
##        subgroup with appropriately filtered mutation
##        matrices and clinical annotations.
##
##   7) Export
##        All OncoPrints are exported as PDF files with
##        consistent layout, legends, and formatting for
##        publication-ready visualization.
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
## OncoPrint ---> no clinbical metadata
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

pdf(file.path(dir_output, 'fig1.pdf'),  width = 7, height = 8)

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
#  Therapy =  factor(clin$Therapy_status, levels = c("Yes", "No")),
#  Location = factor(clin$Location, levels = c('Lobar', 'Cerebellum', 'Thalamic', 'Other')),
  Histo = factor(clin$Histo, levels = c('Glioblastoma', 'Astrocytoma', 'Circumscribed glioma',
                                        'Glioneuronal tumor', 'Oligodendroglioma',  'Pediatric-type HGG',
                                        'Pediatric-type LGG', 'Ependymoma')),
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
    Mut = "#55A868"
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
clin_wt <- clin[clin$IDH_status == 'WT', ]

anno_df <- data.frame(
  Sex = factor(clin_wt$Sex, levels = c("Male", "Female")),
  Age = factor(clin_wt$Age, levels = c(">40", "<40")),
  #IDH = factor(clin_wt$IDH_status, levels = c("WT", "Mut")),
  #Therapy =  factor(clin_wt$Therapy_status, levels = c("Yes", "No")),
  #Location = factor(clin_wt$Location, levels = c('Lobar', 'Cerebellum', 'Thalamic', 'Other')),
  Histo = factor(clin_wt$Histo, c('Glioblastoma', 'Astrocytoma', 'Circumscribed glioma',
                                  'Glioneuronal tumor', 'Pediatric-type HGG', 'Pediatric-type LGG', 
                                  'Ependymoma')),
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
 #   Mut = "#55A868"
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
 # Location = c(
 #   'Lobar'  = "#A54B2DFF",
 #   'Cerebellum' = "#577E2FFF",
 #   'Thalamic'  = "#B49696FF",
 #   'Other' = "#96A5A5FF"
 # ),
  Histo = c(
    'Glioblastoma'          = "#4A7169FF",
    'Astrocytoma'           = "#735231FF",
    'Circumscribed glioma'  = "#E76254FF", 
    'Glioneuronal tumor'    = "#99B6BDFF",
  #  'Oligodendroglioma'     = "#E3CA97FF",
    'Pediatric-type HGG'    = "#B49696FF",
    'Pediatric-type LGG'    = "#AAC197FF",
    'Ependymoma'            = "#96A5A5FF"
  )
  # Therapy = c(
  #  'Yes'       = "#846D86FF",
  #  'No'        = "#ABB2A5FF"
  #)
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

freq <- rowMeans(mut_wt != '')
keep <- freq[freq >= 0.03]
mut_wt_filtered <- mut_wt[rownames(mut_wt) %in% names(keep), ]

p <- oncoPrint(
  mut_wt_filtered,
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

pdf(file.path(dir_output, 'fig3_wt.pdf'), width = 6.5, height = 5)

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
 #IDH = factor(clin_mut$IDH_status, levels = c("WT", "Mut")),
 # Therapy =  factor(clin_mut$Therapy_status, levels = c("Yes", "No")),
 # Location = factor(clin_mut$Location, levels = c('Lobar', 'Cerebellum')),
  Histo = factor(clin_mut$Histo, levels = c('Astrocytoma', 'Oligodendroglioma')),
  Grade = factor(clin_mut$Grade, levels = c('II', 'III', 'IV'))
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
    'II'  = "#BC8E7DFF",
    'III' = "#FAE093FF",
    'IV' = "#7C7189FF"
  ),
 # Location = c(
 #   'Lobar'  = "#A54B2DFF",
 #   'Cerebellum' = "#577E2FFF"
 # ),
  Histo = c(
    'Astrocytoma'           = "#735231FF",
    'Oligodendroglioma'     = "#E3CA97FF"
  )
  # Therapy = c(
  #  'Yes'       = "#846D86FF",
  #  'No'        = "#ABB2A5FF"
  #)
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

pdf(file.path(dir_output, 'fig4_mut.pdf'), width = 6.5, height = 5)

draw(
  p,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

