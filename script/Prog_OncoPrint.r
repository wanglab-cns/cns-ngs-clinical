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
  "snv/indel"           = "SNV/Indel",
  "copy number variant" = "CNV",
  "fusion"              = "Fusion"
)

priority <- c("Fusion", "CNV", "SNV/Indel")

fix_cell <- function(x) {
  if (is.na(x) || x == "") return("")

  parts <- strsplit(x, ";", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts <- parts[parts != ""]
  parts <- tolower(parts)
  parts <- gsub("\\s+", " ", parts)  

  # map known labels
  parts <- ifelse(parts %in% names(label_map), unname(label_map[parts]), parts)

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
combo_levels <- c("CNV;SNV/Indel", "Fusion;CNV", "Fusion;CNV;SNV/Indel")

mut_other <- mut_fixed
mut_other[mut_other %in% combo_levels] <- "Other"

col <- c(
  "SNV/Indel" = "#68855CFF",
  "CNV"       = "#A06177FF",
  "Fusion"    = "#526A83FF",
  "Other"     = "#D9AF6BFF"
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
  "CNV" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(2, "pt"), h - unit(2, "pt"),
              gp = gpar(fill = col["CNV"], col = NA))
  },
  "Fusion" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(2, "pt"), h - unit(2, "pt"),
              gp = gpar(fill = col["Fusion"], col = NA))
  },
  "Other" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(2, "pt"), h - unit(2, "pt"),
              gp = gpar(fill = col["Other"], col = NA))
  }
)

heatmap_legend_param <- list(
  title  = "Alterations",
  at     = c("SNV/Indel", "CNV", "Fusion", "Other"),
  labels = c("SNV/Indel", "CNV", "Fusion", "Other")
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
## OncoPrint ---> clinbical metadata
####################################################
# sample-level metadata
clin$Age <- ifelse(clin$Age >= 40, '> 40', '< 40')
clin$IDH.status[clin$IDH.status == 'NA'] <- "Unknown"
clin$IDH.status[clin$IDH.status == 'IDHmut'] <- "Mut"
clin$IDH.status[clin$IDH.status == 'IDHwt'] <- "WT"

anno_df <- data.frame(
  Sex = factor(clin$Sex, levels = c("Male", "Female")),
  Age = factor(clin$Age, levels = c("> 40", "< 40")),
  IDH = factor(clin$IDH.status, levels = c("WT", "Mut", "Unknown"))
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
    Unknown = "grey"
  ),
  Age = c(
    '> 40'  = "#08306B",
    '< 40' = "#C6DBEF"
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

pdf(file.path(dir_output, 'fig2.pdf'), width = 6, height = 7)
p
dev.off()
