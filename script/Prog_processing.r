#####################################################
## Script: CNS NGS clinical + mutation integration,
##         mutation matrices (binary + oncoprint),
##         and SummarizedExperiment export
##
## Purpose:
##   Load CNS NGS clinical metadata and curated mutation calls from the
##   CNS Tumours NGS & Clinical Data Lock workbook; derive (i) a binary
##   gene-level mutation matrix for downstream modeling and (ii) an
##   oncoprint-ready event-type matrix that preserves alteration classes;
##   harmonize common patients and export SummarizedExperiment objects.
##
## Inputs:
##   - data/Jan 21 2026 CNS Tumours NGS & Clinical Data Lock.xlsx
##       * Sheet 2: Clinical data
##       * Sheet 3: Mutation data
##
## Outputs:
##   - result/data/se_mut_bin_clin.RData
##       SummarizedExperiment with:
##         assay: mat_bin (genes x patients; 0/1 any alteration)
##         colData: curated clinical metadata (OS time/event, histology, IDH, age)
##   - result/data/se_mut_onco_clin.RData
##       SummarizedExperiment with:
##         assay: mat_onco (genes x patients; "" or "type1;type2" per cell)
##         colData: curated clinical metadata (same as above)
##
## Key steps:
##   1) Load clinical + mutation sheets; set patient ID column to 'Study'.
##   2) Repair invalid clinical column names (NA/blank) to enable dplyr verbs.
##   3) Curate OS variables:
##        - os.event: Dead=1; Alive/LTF=0; NA preserved
##        - os.time: Survival (Months)
##      and create analysis-friendly fields (histo, IDH.status, Age).
##   4) Filter mutation records:
##        - remove Unclassified mutation types
##        - remove 'Multi-Gene' entries
##   5) Build mutation matrices:
##        - Option A (mat_bin): binary gene-level (any alteration per gene/patient)
##        - Option B (mat_onco): oncoprint-style strings (semicolon-separated types)
##   6) Harmonize common patients between clinical metadata and matrices; ensure
##      consistent ordering; export SummarizedExperiment objects.
#####################################################
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
clin <- readxl::read_xlsx(file.path(dir_input, 'Jan 21 2026 CNS Tumours NGS & Clinical Data Lock.xlsx'), 
                          sheet = 2, .name_repair = "minimal")
colnames(clin)[1] <- 'Study'
clin <- clin[!is.na(clin$Study), ] ## remove NAs for patients
length(unique(clin$Study))   # 279 samples

## -- Step 2: mutation data
mut_dat <- read_xlsx(file.path(dir_input, 'Jan 21 2026 CNS Tumours NGS & Clinical Data Lock.xlsx'),, 
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
  mutate(
    os.event = case_when(
      is.na(`Survival Status`) ~ NA_integer_,
      `Survival Status` == "Dead" ~ 1L,
      TRUE ~ 0L
    )
  )

## rename variable 
clin$os.time <- clin$'Survival (Months)'
clin$histo <- clin$"Histology"
clin$IDH.status <- clin$'IDH status'
clin$Age <- clin$'Age at diagnosis'
clin <- clin[order(clin$Study), ]

#################################################
## Create mutation matrix
#################################################
mut_dat <- mut_dat[!is.na(mut_dat$Study), ]
mut_dat <- mut_dat[mut_dat$'Mutation Type' != 'Unclassified', ] # 1 patient  
mut_dat <- mut_dat[mut_dat$Gene != 'Multi-Gene', ] 
mut_dat <- mut_dat[order(mut_dat$Study), ]

mut_dat2 <- mut_dat %>%
  mutate(
    Mut_Class = case_when(
      `Mutation Type` == "SNV/Indel" ~ "SNV/Indel",
      `Mutation Type` == "Fusion" ~ "Fusion",
      `Mutation Type` == "Copy Number Variant" &
        str_detect(`Full Mutation Description`, "Amplification") ~ "Amplification",
      `Mutation Type` == "Copy Number Variant" &
        str_detect(`Full Mutation Description`, "Deletion") ~ "Deletion",
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
  filter(!is.na(Study), Study != "", !is.na(Gene)) %>%
  transmute(patient = Study, gene = Gene) %>%
  distinct() %>%
  mutate(value = 1L) %>%
  pivot_wider(names_from = patient, values_from = value, values_fill = 0L) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

mat_bin <- mat_bin[, order(colnames(mat_bin))]

## Option B: oncoprint event-type matrix
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
    value = paste(sort(unique(Mut_Class)), collapse = ";"),
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

which(mat_onco == "Amplification;Fusion;SNV/Indel", arr.ind = TRUE)
mut_dat2 %>%
  filter(Gene == rownames(mat_onco)[7]) %>%
  filter(Study %in% c("C22-081","C23-131","C24-035","C24-045")) %>%
  arrange(Study)

## Common patients  
int <- intersect(clin$Study, colnames(mat_bin)) # 213 patients (all) & 199 for Tier 1 & II
clin <- clin[clin$Study %in% int, ]
clin <- as.data.frame(clin)
rownames(clin) <- clin$Study

mat_bin <- mat_bin[, colnames(mat_bin) %in% int]
mat_onco <- mat_onco[, colnames(mat_onco) %in% int]

## SE object
eset <- SummarizedExperiment(assay= list("gene_expression"=mat_bin),    
                            colData=clin)

save(eset, file = file.path(dir_output, 'se_mut_bin_clin.RData'))

eset <- SummarizedExperiment(assay= list("gene_expression"=mat_onco),    
                            colData=clin)

save(eset, file = file.path(dir_output, 'se_mut_onco_clin.RData'))
