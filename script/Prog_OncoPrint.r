#####################################################
## Script: CNS NGS OncoPrint (ComplexHeatmap)
##
## Purpose:
##   Generate an OncoPrint from a curated CNS NGS mutation event matrix stored
##   in a SummarizedExperiment object. Standardize alteration labels, optionally
##   collapse multi-event cells into "Other", and export the figure as a PDF.
##
## Input:
##   - result/data/se_mut_onco_clin.RData
##       SummarizedExperiment with:
##         assay: gene x patient character matrix ("" or "type1;type2")
##
## Output:
##   - result/Fig/fig1.pdf
##
## Key steps:
##   1) Load SummarizedExperiment and extract the oncoprint matrix.
##   2) Relabel alteration types to a consistent vocabulary:
##        snv/indel -> SNV/Indel
##        copy number variant -> CNV
##        fusion -> Fusion
##      and enforce within-cell ordering (Fusion > CNV > SNV/Indel).
##   3) Collapse combined event labels (e.g., "Fusion;CNV") into "Other" for a
##      simplified legend (optional).
##   4) Draw and export OncoPrint using ComplexHeatmap.
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
## OncoPrint ---> ComplexHeatmap
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
    grid.rect(x, y, w - unit(2, "pt"), h * 0.33,
              gp = gpar(fill = col["SNV/Indel"], col = NA))
  },
  "CNV" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(2, "pt"), h - unit(2, "pt"),
              gp = gpar(fill = col["CNV"], col = NA))
  },
  "Fusion" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(2, "pt"), h * 0.33,
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
left_annotation =  rowAnnotation(
        rbar = anno_oncoprint_barplot(
            axis_param = list(direction = "reverse")
    )),
    right_annotation = NULL)

pdf(file.path(dir_output, 'fig1.pdf'), width = 5, height = 7)
p
dev.off()

## ADD ANNOTATIONS