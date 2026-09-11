library(dplyr)
library(data.table)
library(PheWAS)

# To create phecodes
setwd("~/Documents/PRS_CR/data")
ori_ICD <- readRDS("ICD.rds") # Load the ICD data frame; in our case, it contains SUBJID, AGE_AT_EVENT, ICD_CODE, ICD_FLAG, and site
ex_ICD <- fread(file = "exclusion_ICD.csv") # Load predefined ICD codes for exclusion due to their potential effects on lowering WBC
ex_ICD <- ex_ICD[-c(1:2),]
colnames(ex_ICD) <- c("Type", "ICD_CODE", "Descr")  

df_chemo <- ori_ICD %>% filter(ICD_CODE %in% unique(ex_ICD$ICD_CODE)) # Extract SUBJIDs of individuals with exclusion ICD codes
ICD_clean <- ori_ICD %>% filter(!(SUBJID %in% unique(df_chemo$SUBJID))) # Exclude individuals with predefined exclusion ICD codes

ICD_clean <- ICD_clean %>% filter(AGE_AT_EVENT != ".") %>% 
              mutate(AGE_AT_EVENT = as.numeric(AGE_AT_EVENT)) %>% 
              filter(!is.na(AGE_AT_EVENT))

# Create a specific data frame form to apply "createPhenotypes" function in PheWAS package
ICD_count <- ICD_clean %>% select(-c(site)) %>%
            dplyr::group_by(SUBJID, ICD_CODE) %>%
            dplyr::mutate(count = n()) %>%
            ungroup()

uni_ICD_count <- ICD_count %>%
              group_by(SUBJID, ICD_CODE, count) %>%
              arrange(AGE_AT_EVENT) %>%
              slice(1) %>%
              ungroup()

# extract sex column
df_race <- fread(file = paste0("sample_manifest.csv")); 
ind_id <- uni_ICD_count %>% pull(SUBJID) %>% unique()
df_race2 <- df_race %>% select(IID, sex2) %>% rename(SUBJID = IID)
sex.dt <- uni_ICD_count %>% left_join(df_race2, by = "SUBJID") %>% select(SUBJID, sex2)
sex.dt$SUBJID <- as.character(sex.dt$SUBJID)
sex.dt <- sex.dt %>% mutate(sex = case_when(
            sex2 == "male" ~ "M",
            sex2 == "female" ~ "F"))
sex.dt <- as.data.table(sex.dt)
sex.dt <- sex.dt %>% rename(id = SUBJID) %>% select(-sex2)

# When we have the count column which is manually counted, use aggregate.fun = sum. The count indicates how many times a patient had a particular ICD code.
# The createPhenotypes function allows to have long format as well (don't necessarily need to manually create a count column.)

foo <- uni_ICD_count %>% rename(id = SUBJID, code = ICD_CODE) %>%
        mutate(vocabulary_id = case_when(
          ICD_FLAG == "9"  ~ "ICD9CM",
          ICD_FLAG == "10" ~ "ICD10CM",
          TRUE           ~ NA_character_
        )) %>% select(-ICD_FLAG)

foo <- foo %>% mutate(id = as.character(id))
foo <- foo %>% select(-AGE_AT_EVENT)
foo <- as.data.table(foo)
foo <- foo %>% select(id, vocabulary_id, code, count)

# Note that column in the data table should be ordered id, vocab, code, count for "createPhenotypes" function.
phenotypes = createPhenotypes(foo, 
                              min.code.count = 2, 
                              aggregate.fun = sum,
                              id.sex = sex.dt)

pheno_dt <- data.table(phenotypes)

# Remove duplicated rows
# Returns a wide matrix where each row is a unique patient ID
# Each column represents a aggregated Phecode count.

pheno_dt_unique <- unique(pheno_dt, by = "id")
#saveRDS(pheno_dt, "/home/sum/pheno_dt.rds")
#saveRDS(pheno_dt_unique, "/home/sum/pheno_dt_unique.rds")

