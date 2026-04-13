############################################################
## Script: CNS NGS Clinical–Mutation Integration
##         Binary and Oncoprint Mutation Matrices
##         SummarizedExperiment Export
##
## Purpose:
##   Integrate curated clinical metadata with Tier I/II
##   next-generation sequencing (NGS) mutation calls for
##   central nervous system (CNS) tumor patients.
##
##   Construct two gene-level mutation representations:
##     (1) A binary mutation matrix encoding the presence
##         or absence of any alteration per gene per patient.
##     (2) An oncoprint-ready matrix preserving alteration
##         classes and multi-hit events.
##
##   Harmonize patient identifiers across clinical and
##   molecular data and export analysis-ready
##   SummarizedExperiment objects for downstream modeling,
##   survival analysis, clustering, and visualization.
##
##
## Inputs:
##   - data/Feb 18 2026 CNS Tumours NGS & Clinical Data Lock.xlsx
##       • Sheet 2: Clinical metadata
##       • Sheet 3: Curated mutation calls
##
##   - data/clin.xlsx
##       • Curated clinical updates including refined histology
##         annotations (Histo_updated)
##
##
## Outputs:
##
##   1) result/data/se_mut_bin_clin.RData
##        SummarizedExperiment containing:
##          - assay: binary gene-level mutation matrix (0/1)
##          - colData: curated clinical metadata
##
##   2) result/data/se_mut_onco_clin.RData
##        SummarizedExperiment containing:
##          - assay: oncoprint event-type mutation matrix
##                   ("" or "type1;type2;...")
##          - colData: curated clinical metadata
##
##
## Processing Overview:
##
##   1) Data Import and Standardization
##        Clinical and mutation sheets are imported from the
##        NGS & Clinical Data Lock workbook. The patient
##        identifier column is standardized as 'Study'.
##        Invalid or blank column names in the clinical sheet
##        are repaired to ensure compatibility with dplyr
##        operations.
##
##   2) Clinical Data Curation
##        Survival variables are harmonized:
##           - os.event: Dead = 1; Alive/LTF = 0; NA preserved
##           - os.time:  Survival (Months)
##
##        Key analysis variables are standardized,
##        including histology, IDH status (derived from SNV/Indel
##        mutations in IDH1/IDH2), age at diagnosis, and therapy
##        status (binary indicator for RT and/or concurrent
##        temozolomide).
##
##   3) External Clinical Annotation Integration
##        Updated histology annotations are imported from an
##        external clinical dataset (clin.xlsx). The variable
##        Histo_updated is merged into the main clinical dataset
##        using the patient identifier ('Study') and stored as
##        'histo' for downstream analysis.
##
##   4) Mutation Filtering and Recoding
##        Mutation calls are restricted to Tier I and Tier II
##        events. Records labeled as "Unclassified" are removed.
##
##        Mutation types are recoded into harmonized classes:
##           - SNV/Indel
##           - Fusion
##           - Amplification (copy number gains)
##           - Deletion (copy number losses)
##
##   5) Binary Gene-Level Mutation Matrix (mat_bin)
##        A gene × patient matrix is constructed where:
##           - 1 indicates ≥1 alteration in the gene
##           - 0 indicates no alteration
##        Multiple alterations in the same gene/patient pair
##        are collapsed to a single binary indicator.
##
##   6) Oncoprint Event-Type Matrix (mat_onco)
##        A gene × patient matrix is constructed preserving
##        alteration classes. Multiple events per gene/patient
##        are concatenated as semicolon-separated strings
##        (e.g., "Amplification;SNV/Indel").
##
##   7) IDH Derivation
##        Patient-level IDH status is derived based on the
##        presence of IDH1 or IDH2 SNV/Indel alterations.
##        Status is encoded as Mut or WT and appended to
##        clinical metadata.
##
##   8) Patient Harmonization
##        Only patients present in both curated clinical data
##        and mutation matrices are retained. Clinical metadata
##        and mutation matrices are aligned and ordered
##        consistently.
##
##   9) MultiAssayExperiment Export
##        Two SummarizedExperiment objects (binary and oncoprint) are combined into a single
##        MultiAssayExperiment with shared clinical metadata.
############################################################
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
clin <- clin %>%
  mutate(
    RT_Treated = case_when(
      `Was patient treated with radiotherapy?.1` == "Yes" ~ 1L,
      `Was patient treated with radiotherapy?.1` == "No"  ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    Concurrent_TMZ = case_when(
      `Was temozolomide given concurrent to radiotherapy?` == "Yes" ~ 1L,
      `Was temozolomide given concurrent to radiotherapy?` == "No"  ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    Therapy_status = case_when(
      RT_Treated == 1 | Concurrent_TMZ == 1 ~ 1L,
      RT_Treated == 0 & Concurrent_TMZ == 0 ~ 0L,
      TRUE ~ NA_integer_
    )
  )

#################################################
## Create mutation matrix
#################################################
mut_dat <- mut_dat[!is.na(mut_dat$Study), ]
mut_dat <- mut_dat[mut_dat$'Mutation Type' != 'Unclassified', ] # 1 patient  
mut_dat <- mut_dat[order(mut_dat$Study), ]

## extract IDH status
IDH_status <- mut_dat %>%
  group_by(`Study`) %>%
  summarise(
    IDH_status = ifelse(
      any(Gene %in% c("IDH1", "IDH2") & `Mutation Type` == "SNV/Indel"), 
      "Mut",
      "WT"
    ),
    .groups = "drop"
  )

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
clin_updated <- readxl::read_xlsx(file.path(dir_input, 'clin.xlsx'))
clin_updated <- clin_updated[order(clin_updated$Study), ]

int <- intersect(clin$Study, colnames(mat_bin)) # 213 patients (all) & 199 for Tier 1 & II
clin <- clin[clin$Study %in% int, ]
clin <- as.data.frame(clin)
rownames(clin) <- clin$Study
mat_bin <- mat_bin[, colnames(mat_bin) %in% int]
mat_onco <- mat_onco[, colnames(mat_onco) %in% int]

clin <- clin %>%
  left_join(IDH_status, by = "Study")

clin <- clin %>%
  left_join(
    clin_updated %>% 
      select(Study, histo = Histo_updated),
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
