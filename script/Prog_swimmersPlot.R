##-------------------------------------------------------------------
## Script: CNS NGS Targeted Therapy Swimmer Plot (done by Steven Cook)
##
## Purpose:
##   To visualize progression-free survival (PFS) outcomes
##   following molecularly guided targeted therapy in
##   patients with glioma using a swimmer plot framework.
##
##   The plot summarizes treatment duration, radiographic
##   response, treatment status, molecular alterations,
##   and clinical subgroup classification at the individual
##   patient level.
##
##
## Input:
##
##   - data/CNS_NGS_Swimmers_Plot_Data.csv
##
##       Patient-level dataset containing:
##         • clinical characteristics
##         • targeted therapy information
##         • actionable molecular alterations
##         • treatment dates
##         • response assessments
##         • follow-up information
##
##
## Outputs:
##
##   1) Publication-quality swimmer plot PDF
##
##        - swimmers_plot.pdf
##
##   2) Visualization of:
##
##        - Progression-free survival duration
##        - Best radiographic response
##        - Treatment status
##        - Molecular alterations
##        - Targeted therapies administered
##        - Clinical subgroup membership
##
##
## Processing Overview:
##
##   1) Import patient-level targeted therapy dataset
##
##   2) Standardize and parse treatment and follow-up dates
##
##   3) Derive progression-free survival (PFS) from
##        targeted therapy initiation
##
##   4) Annotate actionable molecular alterations using
##        simplified gene-level labels
##
##   5) Classify patients into molecular subgroups:
##
##        - Glioblastoma (GBM)
##        - BRAF-driven non-GBM glioma
##        - IDH-mutant glioma
##
##   6) Order patients within each subgroup by decreasing
##        PFS duration
##
##   7) Generate swimmer plot displaying:
##
##        - PFS duration as horizontal bars
##        - Best response by bar colour
##        - Treatment status by endpoint symbol
##        - Targeted therapy administered
##        - Actionable molecular alteration
##        - GBM status
##        - Patient identifier
##
##   8) Add subgroup annotations, legends, and
##        publication-ready formatting
##
##   9) Export final figure as PDF
##
##
## Notes:
##
##   - PFS is calculated from targeted therapy initiation
##     to progression, treatment discontinuation, or last
##     follow-up for ongoing therapy.
##
##   - Bars exceeding 25 months are truncated and
##     indicated with arrowheads.
##
##   - Best response categories include:
##        • CR = Complete response
##        • PR = Partial response
##        • SD = Stable disease
##        • PD = Progressive disease
##
##   - Treatment status is displayed using endpoint
##     symbols indicating ongoing treatment, completion
##     without progression, or progression while on
##     therapy.
##
##   - Patients are ranked within each molecular subgroup
##     according to decreasing PFS to facilitate visual
##     comparison of treatment outcomes.
##-------------------------------------------------------------------
###############################################################
## load libraries
###############################################################
library(ggplot2)
library(dplyr)
library(stringr)
library(paletteer)

###############################################################
## setup directory
###############################################################
dir_input <- 'data'
dir_output <- 'result/data'

###############################################################
## define functins
###############################################################

parse_date <- function(x) {
  if (inherits(x, "Date"))    return(x)
  if (inherits(x, "POSIXct")) return(as.Date(x))

  x <- as.character(x)

  out <- as.Date(rep(NA_real_, length(x)), origin = "1970-01-01")

  num <- suppressWarnings(as.numeric(x))
  is_num <- !is.na(num)

  out[is_num] <- as.Date(num[is_num], origin = "1899-12-30")

  is_chr <- !is_num & !is.na(x) & !(trimws(x) %in% c("NA", "", "Ongoing"))

  if (any(is_chr)) {
    out[is_chr] <- as.Date(
      suppressWarnings(
        as.POSIXct(
          trimws(x[is_chr]),
          tryFormats = c(
            "%Y-%m-%d",
            "%Y/%m/%d",
            "%d/%m/%Y",
            "%m/%d/%Y",
            "%B %d, %Y",
            "%b %d, %Y"
          )
        )
      )
    )
  }

  out
}

short_mutation <- function(s) {
  s <- as.character(s)
  dplyr::case_when(
    str_detect(s, "Val600Glu")               ~ "BRAF V600E",
    str_detect(s, "Asn581His")               ~ "BRAF N581H",
    str_detect(s, "FGFR3") & str_detect(s, "TACC3") ~ "FGFR3-TACC3 fusion",
    str_detect(s, "GOPC")  & str_detect(s, "ROS1")  ~ "GOPC-ROS1 fusion",
    str_detect(s, "KIAA1549")                ~ "KIAA1549-BRAF fusion",
    str_detect(s, "BCAS1")                   ~ "BCAS1-BRAF fusion",
    str_detect(s, "NF1")                     ~ "NF1 mut",
    str_detect(s, "Arg132Gly")               ~ "IDH1 R132G",
    str_detect(s, "Arg132His")               ~ "IDH1 R132H",
    str_detect(s, "Arg172Lys")               ~ "IDH2 R172K",
    TRUE                                      ~ s
  )
}


get_col <- function(pattern) nm[str_detect(nm, regex(pattern, ignore_case = TRUE))][1]


add_legend_block <- function(p, title, y0, labels, shapes, colours,
                             fills = NA, line_gap = 0.85) {
  ys <- y0 - seq_along(labels) * line_gap
  p +
    annotate("text", x = x_leg_g - 0.15, y = y0, label = title,
             hjust = 0, fontface = "bold", size = 4.5) +
    annotate("point", x = rep(x_leg_g, length(labels)), y = ys,
             shape = shapes, colour = colours,
             fill = if (length(fills) == 1) rep(fills, length(labels)) else fills,
             size = 4, stroke = 1.2) +
    annotate("text", x = rep(x_leg_t, length(labels)), y = ys,
             label = labels, hjust = 0, size = 4)
}

###############################################################
## data preparation
###############################################################

Swimmersplot <- read.csv(
  file.path(dir_input, "CNS_NGS_Swimmers_Plot_Data.csv"),
  fileEncoding = "Windows-1252",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
df <- Swimmersplot
nm <- names(df)

col_ptnum   <- get_col("^Patient Number")
col_ptid    <- get_col("^Patient Identifier")
col_therapy <- get_col("^Targeted Therapy$")
col_mut     <- get_col("Actionable Mutation")
col_start   <- get_col("Targeted Therapy Start Date")
col_end     <- get_col("Targeted Therapy End Date")
col_status  <- get_col("Treatment Status")
col_reason  <- get_col("Reason Stopped")
col_resp    <- get_col("Best Response")
col_gbm     <- get_col("^GBM")
col_lfu     <- get_col("Last Follow Up Date")

dat <- df %>%
  transmute(
    pt_num    = suppressWarnings(as.integer(.data[[col_ptnum]])),
    pt_id     = .data[[col_ptid]],
    therapy   = trimws(.data[[col_therapy]]),
    mut_raw   = .data[[col_mut]],
    start_dt  = parse_date(.data[[col_start]]),
    end_dt    = parse_date(.data[[col_end]]),
    lfu_dt    = parse_date(.data[[col_lfu]]),
    status    = trimws(.data[[col_status]]),
    reason    = trimws(.data[[col_reason]]),
    response  = trimws(.data[[col_resp]]),
    gbm       = trimws(.data[[col_gbm]])
  ) %>%
  # Drop legend / note rows at the bottom of the sheet
  filter(!is.na(pt_num), !is.na(pt_id), !is.na(therapy), therapy != "") %>%
  arrange(pt_num)

# End date used for PFS:
#   Ongoing                       -> Last Follow Up Date
#   Reason == Disease Progression -> Targeted Therapy End Date
#   otherwise                     -> Targeted Therapy End Date
dat <- dat %>%
  mutate(
    end_for_pfs = dplyr::case_when(
      status == "Ongoing"                  ~ lfu_dt,
      reason == "Disease Progression"      ~ end_dt,
      TRUE                                 ~ end_dt
    ),
    pfs_months = as.numeric(end_for_pfs - start_dt) / 30.44,
    prog_on_tx = reason == "Disease Progression",
    mutation   = short_mutation(mut_raw),
    is_gbm     = str_to_lower(gbm) == "yes",
    response   = factor(response, levels = c("PD", "SD", "PR", "CR")),
    # End-of-bar marker category
    marker = dplyr::case_when(
      status == "Ongoing" ~ "Ongoing",
      prog_on_tx          ~ "Progressed on therapy",
      TRUE                ~ "Finished (no progression)"
    ),
    # Clinical group (defined by patient number, as before)
    group = dplyr::case_when(
      pt_num <= 5  ~ "GBM",
      pt_num <= 11 ~ "BRAF",
      TRUE         ~ "IDH"
    ),
    group = factor(group, levels = c("GBM", "BRAF", "IDH")),
    # Cap bars at 25 months; flag bars that hit the cap (need an arrow)
    pfs_capped = pmin(pfs_months, 25),
    hit_cap    = pfs_months >= 25
  )

n_pt <- nrow(dat)   # 28

# Order rows: keep the groups top->bottom (GBM, BRAF, IDH); within each group
# sort by PFS so the LONGEST lane is on top and the shortest at the bottom.
dat <- dat %>%
  arrange(group, desc(pfs_months)) %>%
  mutate(
    row_order = row_number(),
    y_pos     = n_pt - row_order + 1   # first row at the top
  )

###############################################################
## layout constants
###############################################################
# Colours for Best Response
resp_cols <- c(PD = "#A5693CFF", SD = "#526A83FF", PR = "#A5872DFF", CR = "#68855CFF")
resp_labs <- c(PD = "Progressive disease", SD = "Stable disease",
               PR = "Partial response",    CR = "Complete response")

# x anchors for the left-hand annotation columns (text is right-aligned)
x_therapy <- -14.5   # was -9.0
x_gene    <- -7.5    # was -3.4
x_gbm     <- -5.2    # was -2.2
x_pt      <- -2.8    # was -0.6

# Vertical separators sit in the GAPS between columns (never over text):
#   -2.8  : between gene alteration and GBM   (right of gene's right edge)
#   -1.4  : between GBM and Pt
#   -0.05 : between Pt and the bars
sep_x <- c(-6.35, -4, -1.4)

# Far-right legend column (single tidy stack)
x_leg_g <- 33.0   # glyph x
x_leg_t <- 33.7   # label x

grey_line <- "grey80"
header_y  <- n_pt + 1.4

# Group section boundaries, computed from the (possibly reordered) rows so
# they always track the actual group blocks rather than fixed patient numbers.
grp_sizes <- as.numeric(table(factor(dat$group, levels = c("GBM", "BRAF", "IDH"))))
grp_cum   <- cumsum(grp_sizes)
# Divider lines sit between the last row of one group and the first of the next
div_y <- n_pt - grp_cum[-length(grp_cum)] + 0.5

sections <- dat %>%
  group_by(group) %>%
  summarise(ymid = mean(y_pos), .groups = "drop") %>%
  mutate(label = c("GBM", "BRAF-driven", "IDH-mutant")[as.integer(group)]) %>%
  arrange(group)

# ---------------------------------------------------------------------
# 4. Build the plot
# ---------------------------------------------------------------------

p <- ggplot(dat) +
  
  # --- dashed vertical grid lines at 5,10,15,20,25 ---
  geom_vline(xintercept = c(5, 10, 15, 20, 25),
             linetype = "dashed", colour = grey_line, linewidth = 0.3) +
  
  # --- horizontal group divider lines ---
  geom_hline(yintercept = div_y, colour = grey_line, linewidth = 0.3) +
  
  # --- bars (length = PFS, colour = Best Response) ---
  geom_segment(aes(x = 0, xend = pfs_capped, y = y_pos, yend = y_pos,
                   colour = response), linewidth = 4.2, lineend = "butt") +
  
  # --- arrow for any bar that reaches the 25-month cap ---
  geom_segment(data = filter(dat, hit_cap),
               aes(x = pfs_capped - 0.01, xend = 25.6, y = y_pos, yend = y_pos,
                   colour = response), linewidth = 4.2,
               arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  
  # --- end-of-bar markers ---
  geom_point(data = filter(dat, marker == "Ongoing"),
             aes(x = pfs_capped, y = y_pos), shape = 17, size = 2.6,
             colour = "black") +
  geom_point(data = filter(dat, marker == "Finished (no progression)"),
             aes(x = pfs_capped, y = y_pos), shape = 15, size = 2.4,
             colour = "black") +
  geom_point(data = filter(dat, marker == "Progressed on therapy"),
             aes(x = pfs_capped, y = y_pos), shape = 4, size = 2.6,
             stroke = 1.3, colour = "black") +
  
  # --- PFS value to the right of each bar end (only for bars < 25) ---
  geom_text(data = filter(dat, !hit_cap),
            aes(x = pfs_capped + 0.45, y = y_pos,
                label = sprintf("%.1f", pfs_months)),
            hjust = 0, size = 4, colour = "grey20") +
  
  # ----- Left annotation columns -----
geom_text(aes(x = x_therapy, y = y_pos, label = therapy),
          hjust = 1, size = 4.5) + # fontface = "bold"
  geom_text(aes(x = x_gene, y = y_pos, label = mutation),
            hjust = 1, size = 4, colour = "grey15") +
  # GBM dot: filled navy = GBM, open = non-GBM
  geom_point(aes(x = x_gbm, y = y_pos, shape = is_gbm),
             size = 3.8, colour = "#4B5A69FF", fill = "#4B5A69FF", stroke = 0.9) +
  geom_text(aes(x = x_pt, y = y_pos, label = pt_num),
            hjust = 1, size = 4, colour = "grey15") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  
  # ----- Column headers + subtle underline -----
annotate("text", x = x_therapy, y = header_y, label = "Targeted therapy",
         hjust = 1, fontface = "bold", size = 5) +
  annotate("text", x = x_gene, y = header_y, label = "Gene alteration",
           hjust = 1, fontface = "bold", size = 5) +
  annotate("text", x = x_gbm, y = header_y, label = "GBM",
           hjust = 0.5, fontface = "bold", size = 5) +
  annotate("text", x = x_pt, y = header_y, label = "Pt.",
           hjust = 1, fontface = "bold", size = 5) +
  annotate("segment", x = x_therapy - 4.0, xend = x_pt + 0.3,
           y = header_y - 0.6, yend = header_y - 0.6,
           colour = "grey50", linewidth = 0.4) +
  
  # ----- Thin vertical separators between annotation columns -----
annotate("segment", x = sep_x, xend = sep_x,
         y = 0.4, yend = header_y - 0.2,
         colour = grey_line, linewidth = 0.3) +
  
  # ----- Italic grey section labels (left-aligned, just right of the bars) -----
annotate("text", x = 25.5, y = sections$ymid, label = sections$label,
         hjust = 0, fontface = "italic", colour = "grey55", size = 5.2) +
  
  # ----- Bar colours (legend drawn manually below, so no auto guide) -----
scale_colour_manual(values = resp_cols, drop = FALSE, guide = "none") +
  
  # ----- Scales / coords -----
scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25),
                   limits = c(x_therapy - 5, 41),
                   expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, header_y + 0.6), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "Progression-free survival from therapy start (months)", y = NULL) +
  
  theme_classic(base_size = 14) +
  theme(
plot.background = element_rect(fill = "white", colour = NA),
panel.background = element_rect(fill = "white", colour = NA),
axis.line.y = element_blank(),
axis.ticks.y = element_blank(),
axis.text.y = element_blank(),
axis.line.x = element_line(colour = "grey40"),
axis.title.x = element_text(size = 16),
axis.text.x = element_text(size = 14),
legend.position = "none",
plot.margin = margin(15, 15, 15, 15)
)

# ---------------------------------------------------------------------
# 5. Legends
# ---------------------------------------------------------------------

# --- Group (top) ---
p <- add_legend_block(
  p, "Group", y0 = 27,
  labels  = c("GBM", "Non-GBM"),
  shapes  = c(16, 1),
  colours = c("#4B5A69FF", "#4B5A69FF"))

# --- Treatment status (middle) ---
p <- add_legend_block(
  p, "Treatment status", y0 = 23,
  labels  = c("Ongoing", "Finished (no progression)", "Progressed on therapy"),
  shapes  = c(17, 15, 4),
  colours = c("black", "black", "black"))

# --- Best response (bottom) — filled squares in the response colours ---
p <- add_legend_block(
  p, "Best response", y0 = 18,
  labels  = c("Progressive disease", "Stable disease",
              "Partial response", "Complete response"),
  shapes  = c(15, 15, 15, 15),
  colours = unname(resp_cols[c("PD", "SD", "PR", "CR")]))

# ---------------------------------------------------------------------
# 6. Display
# ---------------------------------------------------------------------
dev.new(width = 22, height = 11)

print(p)

ggsave(file.path(dir_output, "swimmers_plot.pdf"), p, width = 22, height = 11)
