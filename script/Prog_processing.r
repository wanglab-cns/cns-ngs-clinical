##----------------------------------------------------------------------------------------------------
## Script: CNS NGS Clinical–Mutation Integration
##         Binary and Oncoprint Mutation Matrices
##
## Purpose:
##   To integrate curated clinical metadata with CNS tumour
##   next-generation sequencing (NGS) alteration data and
##   generate harmonized molecular and clinical datasets for
##   downstream analyses.
##
##   The workflow constructs:
##     (1) A binary gene-level mutation matrix
##     (2) An oncoprint-compatible mutation matrix
##
##   Clinical and molecular data are harmonized and exported
##   within a unified MultiAssayExperiment object for
##   downstream survival analysis, biomarker discovery,
##   clustering, and visualization.
##
##
## Inputs:
##
##   1) data/Feb 18 2026 CNS Tumours NGS & Clinical Data Lock.xlsx
##        • Clinical metadata
##        • Curated Tier I/II mutation calls
##
##   2) data/clin.xlsx
##        • Updated clinical annotations
##
##
## Outputs:
##
##   1) result/data/mae_mut_clin.RData
##
##        MultiAssayExperiment containing:
##
##          - Binary mutation matrix
##          - Oncoprint mutation matrix
##          - Harmonized clinical metadata
##
##
## Processing Overview:
##
##   1) Import and standardize clinical and mutation data
##
##   2) Curate and harmonize clinical variables including:
##        - Overall survival
##        - Age
##        - Treatment status
##        - Histology
##
##   3) Filter and recode clinically curated mutation events
##
##   4) Generate:
##        - Binary gene-level mutation matrix
##        - Oncoprint-compatible mutation matrix
##
##   5) Derive patient-level IDH mutation status
##
##   6) Harmonize patients across datasets
##
##   7) Construct and export a MultiAssayExperiment object
##
##
## Notes:
##   - Only Tier I and Tier II clinically curated variants are
##     retained.
##
##   - Mutation classes include:
##        SNV/Indel, Fusion, Amplification, and Deletion.
##
##   - Clinical and molecular datasets are aligned using
##     harmonized patient identifiers.
##
##   - The exported object is designed for downstream analyses
##     within the Bioconductor ecosystem.
##----------------------------------------------------------------------------------------------------
####################################################
## Load libraries
####################################################
library(MultiAssayExperiment)
library(dplyr)
library(readxl)  
library(tidyr)
library(stringr)
library(tibble)

####################################################
## Setup directories
####################################################
dir_input <- 'data'
dir_output <- 'result/data'

####################################################
## Load clinical CNS NGS data
####################################################
## -- Step 1: clinical data
clin <- readxl::read_xlsx(file.path(dir_input, 'Feb 18 2026 CNS Tumours NGS & Clinical Data Lock.xlsx'), 
                          sheet = 2, .name_repair = "minimal")
colnames(clin)[1] <- 'Study'
clin <- clin[!is.na(clin$Study), ] ## remove NAs for patients
length(unique(clin$Study))   # 279 samples

## -- Step 2: mutation data
mut_dat <- read_xlsx(file.path(dir_input, 'Feb 18 2026 CNS Tumours NGS & Clinical Data Lock.xlsx'),
                          sheet = 3, .name_repair = "minimal") 

colnames(mut_dat)[1] <- 'Study'
mut_dat <- mut_dat[!is.na(mut_dat$Study), ]
mut_dat <- mut_dat[mut_dat$Tier %in% c('I', 'II'), ]
length(unique(mut_dat$Study)) # 199 patients

## remove NA columns in clinical data
bad <- which(is.na(names(clin)) | names(clin) == "")
bad
names(clin)[bad]

# 2) Repair names (NA/blank -> unnamed_#, then make all names unique)
nm <- names(clin)
nm[is.na(nm) | nm == ""] <- paste0("unnamed_", which(is.na(nm) | nm == ""))
names(clin) <- make.unique(nm)

# 3) Confirm it's clean
any(is.na(names(clin)) | names(clin) == "")

#################################################
## Curate clinical data
#################################################
## OS event
clin <- clin %>%
  mutate(`Survival Status` = recode(`Survival Status`,
                                   "Ailve" = "Alive"))

clin <- clin %>%
  mutate(
    os.event = case_when(
      is.na(`Survival Status`) ~ NA_integer_,
      `Survival Status` == "Dead" ~ 1L,
      `Survival Status` %in% c("Alive", "LTF") ~ 0L,
      TRUE ~ NA_integer_  # catches unexpected values
    )
  )

## rename variable 
clin$os.time <- clin$'Survival (Months)'
clin$Age <- clin$'Age at diagnosis'
clin <- clin[order(clin$Study), ]

## therapy status 
#clin <- clin %>%
#  mutate(
#    RT_Treated = case_when(
#      `Was patient treated with radiotherapy?.1` == "Yes" ~ 1L,
#      `Was patient treated with radiotherapy?.1` == "No"  ~ 0L,
#      TRUE ~ NA_integer_
#    ),
    
#    Concurrent_TMZ = case_when(
#      `Was temozolomide given concurrent to radiotherapy?` == "Yes" ~ 1L,
#      `Was temozolomide given concurrent to radiotherapy?` == "No"  ~ 0L,
#      TRUE ~ NA_integer_
#    ),
    
#    Therapy_status = case_when(
#      RT_Treated == 1 | Concurrent_TMZ == 1 ~ 1L,
#      RT_Treated == 0 & Concurrent_TMZ == 0 ~ 0L,
#      TRUE ~ NA_integer_
#    )
#  )

#################################################
## Create mutation matrix
#################################################
mut_dat <- mut_dat[!is.na(mut_dat$Study), ]
mut_dat <- mut_dat[mut_dat$'Mutation Type' != 'Unclassified', ] # 1 patient  
mut_dat <- mut_dat[order(mut_dat$Study), ]

## extract IDH status
#IDH_status <- mut_dat %>%
#  group_by(`Study`) %>%
#  summarise(
#    IDH_status = ifelse(
#      any(Gene %in% c("IDH1", "IDH2") & `Mutation Type` == "SNV/Indel"), 
#      "Mut",
#      "WT"
#    ),
#    .groups = "drop"
#  )

## table(IDH_status$IDH_status)

## extract alterations
mut_dat2 <- mut_dat %>%
  mutate(
    Mut_Class = case_when(
      `Mutation Type` == "SNV/Indel" ~ "SNV/Indel",
      `Mutation Type` == "Fusion" ~ "Fusion",
      `Mutation Type` == "Copy Number Variant" &
        str_detect(`Full Mutation Description`, regex("amplification|gain", ignore_case = TRUE)) ~ "Amplification",
      `Mutation Type` == "Copy Number Variant" &
        str_detect(`Full Mutation Description`, regex("deletion|loss", ignore_case = TRUE)) ~ "Deletion",
      TRUE ~ NA_character_
    )
  )

## Option A: binary gene-level matrix
mat_bin <- mut_dat2 %>%
  mutate(
    Study = trimws(Study),
    Gene  = trimws(Gene),
    Gene  = dplyr::na_if(Gene, "")
  ) %>%
  filter(!is.na(Study), Study != "", !is.na(Gene), !is.na(Mut_Class)) %>%
  transmute(patient = Study, gene = Gene) %>%
  distinct() %>%
  mutate(value = 1L) %>%
  pivot_wider(names_from = patient, values_from = value, values_fill = 0L) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

mat_bin <- mat_bin[, order(colnames(mat_bin))]

## Option B: oncoprint event-type matrix
mut_order <- c("SNV/Indel", "Fusion", "Amplification", "Deletion")
mat_onco <- mut_dat2 %>%
  mutate(
    Study     = trimws(Study),
    Gene      = trimws(Gene),
    Mut_Class = trimws(Mut_Class),
    Mut_Class = na_if(Mut_Class, "")
  ) %>%
  filter(
    !is.na(Study), Study != "",
    !is.na(Gene),  Gene  != "",
    !is.na(Mut_Class)
  ) %>%
  group_by(Gene, Study) %>%
  summarise(
    value = paste(
      mut_order[mut_order %in% unique(Mut_Class)],
      collapse = ";"
    ),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = Study,
    values_from = value,
    values_fill = list(value = "")
  )
  
mat_onco <- as.data.frame(mat_onco)
rownames(mat_onco) <- mat_onco$Gene
mat_onco$Gene <- NULL
mat_onco <- as.matrix(mat_onco)

mat_onco <- mat_onco[, sort(colnames(mat_onco)), drop = FALSE]

###############################################
## Create SE objects
###############################################
## Common patients  
clin_updated <- read.csv(file.path(dir_input, 'clin.csv'))
clin_updated <- clin_updated[order(clin_updated$Study), ]

int <- intersect(clin$Study, colnames(mat_bin)) # 213 patients (all) & 199 for Tier 1 & II
clin <- clin[clin$Study %in% int, ]
clin <- as.data.frame(clin)
rownames(clin) <- clin$Study
mat_bin <- mat_bin[, colnames(mat_bin) %in% int]
mat_onco <- mat_onco[, colnames(mat_onco) %in% int]

clin <- clin %>%
  left_join(
    clin_updated %>%
      select(
        Study,
        histo = Histo_updated,
        IDH_status = IDH_status
      ),
    by = "Study"
  )

rownames(clin) <- clin$Study
###################################################################
## MAE object
###################################################################
# Create SE objects
se_bin <- SummarizedExperiment(
  assays = list(binary = mat_bin),
  colData = clin
)

se_onco <- SummarizedExperiment(
  assays = list(oncoprint = mat_onco),
  colData = clin
)

# Combine into MAE
mae <- MultiAssayExperiment(
  experiments = list(
    mut_binary = se_bin,
    mut_oncoprint = se_onco
  ),
  colData = clin
)

# Save
save(mae, file = file.path(dir_output, "mae_mut_clin.RData"))
