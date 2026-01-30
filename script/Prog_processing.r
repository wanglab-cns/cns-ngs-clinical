####################################################
## Load libraries
####################################################
library(MultiAssayExperiment)
library(dplyr)
library(readxl)  
library(tidyr)

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
length(unique(clin$Study))   

## -- Step 2: mutation data
mut_dat <- read_xlsx(file.path(dir_input, 'Jan 21 2026 CNS Tumours NGS & Clinical Data Lock.xlsx'),, 
                          sheet = 3, .name_repair = "minimal") 
colnames(mut_dat)[1] <- 'Study'
mut_dat <- mut_dat[!is.na(mut_dat$Study), ]
length(unique(mut_dat$Study))


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
mut_dat <- mut_dat[mut_dat$'Mutation Type' != 'Unclassified', ]
mut_dat <- mut_dat[mut_dat$Gene != 'Multi-Gene', ]
mut_dat <- mut_dat[, order(colnames(mut_dat))]

## Option A: binary gene-level matrix
mat_bin <- mut_dat %>%
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

## Option B: oncoprint event-type matrix
mat_onco <- mut_dat %>%
  mutate(
    Study = trimws(Study),
    Gene  = trimws(Gene),
    Type  = tolower(`Mutation Type`)   # optional: make labels cleaner
  ) %>%
  filter(!is.na(Gene), Gene != "") %>%
  group_by(Gene, Study) %>%
  summarise(
    value = paste(sort(unique(Type)), collapse = ";"),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = Study,
    values_from = value,
    values_fill = ""
  ) %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()

## Common patients  
int <- intersect(clin$Study, colnames(mat_bin)) # 213 patients
clin <- clin[clin$Study %in% int, ]
clin <- as.data.frame(clin)
rownames(clin) <- clin$Study

mat_bin <- mat_bin[, colnames(mat_bin) %in% int]
mat_bin <- mat_bin[, order(colnames(mat_bin))]
mat_onco <- mat_onco[, colnames(mat_onco) %in% int]
mat_onco <- mat_onco[, order(colnames(mat_onco))]

## SE object
eset <- SummarizedExperiment(assay= list("gene_expression"=mat_bin),    
                            colData=clin)

save(eset, file = file.path(dir_output, 'se_mut_bin_clin.RData'))

eset <- SummarizedExperiment(assay= list("gene_expression"=mat_onco),    
                            colData=clin)

save(eset, file = file.path(dir_output, 'se_mut_onco_clin.RData'))


# Questions: 
- what does 'LTF' mean in survival status? ---> Lost To Follow up ---> censoring
- Any censoring for time-to-event outcome (OS)?
- Histology type? 
- Tier variable (Tier 1? or 1 and 2?) 
- Analysis based on IDH status?  (let's consider both stratification or even as not stratification, IDH wild type GBM)
- Cut-off for age variable? --> median of age at diagnosis? ~ 40 yrs
- 'C23-098' & 'C24-002' with mutation type, Unclassified
