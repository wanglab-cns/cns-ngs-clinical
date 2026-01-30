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

## Common patients ---> Not necessary as MAE object 
int <- intersect(clin$Study, mut_dat$Study) # 213 patients
clin <- clin[clin$Study %in% int, ]
mut_dat <- mut_dat[mut_dat$Study %in% int, ]

#################################################
## Create mutation matrix
#################################################





## ADD SCRIPT OR STEPS FOR MAE OBJECT
table(clin$'Survival Status')
#Alive  Dead   LTF
#  126    59    28

# any censoring for survival? 
summary(clin$'Survival (Months)')
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.230   5.388  13.568  33.312  25.263 321.651       2

table(clin$'Histology...46')                   
#      Astrocytoma       Diffuse HGG      Glioblastoma            Glioma
#               51                 1               105                23
# Oligodendroglioma             Other
#               13                20

table(clin$'IDH status...59') 
#IDHmut  IDHwt     NA
#    50    152     10

table(mut_dat$'Mutation Type')
#Copy Number Variant              Fusion           SNV/Indel        Unclassified 
#                149                  50                 545                   2

mut_dat[mut_dat$'Mutation Type' == 'Unclassified', ]

# for patient C24-002
mut_dat[mut_dat$Gene == 'Multi-Gene', ]

# WHO Grade
table(clin$'WHO 2021 Grade')

# 1   2   3   4  NA
# 22  23  23 136   9

df <- mut_dat[!duplicated(mut_dat$Study), ]
df <- df[order(df$Study), ]
df_clin <- clin[order(clin$Study), ]

table(df$Tier, df_clin$'WHO 2021 Grade') 

#        1   2   3   4  NA
#  I     9  20  17 114   5
#  II   10   0   4  19   1
#  III   0   1   1   0   2

summary(clin$"Age at diagnosis")
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6516 33.5858 46.8648 45.7562 57.4057 82.8956

# Questions: 
- what does 'LTF' mean in survival status? ---> Lost To Follow up ---> censoring
- Any censoring for time-to-event outcome (OS)?
- Histology type? 
- Tier variable (Tier 1? or 1 and 2?) 
- Analysis based on IDH status?  (let's consider both stratification or even as not stratification, IDH wild type GBM)
- Cut-off for age variable? --> median of age at diagnosis? ~ 40 yrs
- 'C23-098' & 'C24-002' with mutation type, Unclassified
