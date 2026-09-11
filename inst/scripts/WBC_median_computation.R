# data frame for WBC median 
library(data.table)
library(dplyr)
setwd("~/data")
df_demo <- fread("GWAS.csv")
df_WBC <- fread("Leukocyte.csv", fill=Inf)
df_ID_leuk <- readRDS("id_ex_wbc.rds") 
df_WBC <- df_WBC %>% filter(!(EMERGEID %in% df_ID_leuk) )
df_WBC <- df_WBC %>% filter(EMERGEID %in% df_demo$SUBJID)
df_WBC <- df_WBC %>% mutate(Leukocyte = case_when(Unit_Leuk == "/CU M" ~ Leukocyte/1000,
                                                  Unit_Leuk == ".CMM" ~ Leukocyte/1000,
                                                  TRUE ~ Leukocyte))
df_WBC <- df_WBC %>% filter(!(is.na(Leukocyte))) 
df_WBC <- df_WBC %>% filter(!(Unit_Leuk %in% c(".HPF", ".hpf", "hpf")))
df_WBC <- df_WBC %>% select(-V8, -V9)
df_WBC <- df_WBC %>% filter(AGE >= 18 & AGE < 90)
df_WBC <- df_WBC %>% filter(Leukocyte < 200)

# multiple measurements withing a week including 0
week_yr <- 7 / 365.25
df_WBC <- df_WBC %>% arrange(EMERGEID, AGE) %>%
          group_by(EMERGEID) %>%
          mutate(week_grp = cumsum(c(0, diff(AGE) > week_yr))) %>%
          ungroup()

df_WBC <- df_WBC %>%
          group_by(EMERGEID, week_grp) %>%
          mutate(
            week_has_zero    = any(Leukocyte == 0, na.rm = TRUE),
            week_has_nonzero = any(Leukocyte != 0, na.rm = TRUE)) %>%
          ungroup()

df_WBC <- df_WBC %>% filter(!(Leukocyte == 0 & week_has_nonzero))
df_WBC <- df_WBC %>% group_by(EMERGEID) %>% mutate(med_age = median(AGE))
df_WBC <- df_WBC %>% select(-AGE)
df_WBC <- df_WBC %>% rename(AGE = med_age)

temp <- df_WBC %>%
        select(EMERGEID, AGE, Leukocyte) %>% 
        group_by(EMERGEID) %>%
        mutate(
          subj_median = Leukocyte[which.min(abs(Leukocyte - median(Leukocyte, na.rm = TRUE)))]
        ) %>% ungroup() %>%
        mutate(logic = Leukocyte == subj_median)

temp_short <- temp %>% filter(logic == TRUE) %>% arrange(EMERGEID, Leukocyte, AGE) %>% 
              distinct(EMERGEID, .keep_all = TRUE) %>%  
              rename(IID = EMERGEID)

temp_short <- temp_short %>% mutate(subj_median_st = as.numeric(scale(subj_median)))

saveRDS(temp_short, file = paste0("/WBC_median_df.rds"))
