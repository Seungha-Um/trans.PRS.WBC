# data frame for WBC median 
library(data.table)
library(dplyr)
setwd("~/Documents/PRS_CR/data")
df_demo <- fread("person_GWAS.csv")
df_WBC <- fread("Leukocyte-EMERGE-2020-05-13.csv", fill=Inf)
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
df_WBC <- df_WBC %>% distinct(EMERGEID, AGE, Leukocyte, Unit_Leuk, .keep_all = TRUE)

wbc_summary <- df_WBC %>%
                group_by(EMERGEID) %>%
                summarise(
                  n_measurements = n(),
                  followup_years = max(AGE) - min(AGE),
                  median_WBC = median(Leukocyte),
                  variation_WBC = sd(Leukocyte, na.rm = TRUE),
                  .groups = "drop"
                )

summary_table <- tibble(
  metric = c(
    "Total participants",
    "Total number of repeated WBC measurements",
    "Median number of measurements per participant (IQR)",
    "Median observation duration, years (IQR)",
    "Median of within-participant WBC standard deviation (IQR)",
    "Median of participant-level median WBC counts (IQR)"
  ),
  value = c(
    n_distinct(df_WBC$EMERGEID),
    nrow(df_WBC),
    sprintf(
      "%.0f (%.0f–%.0f)",
      median(wbc_summary$n_measurements),
      quantile(wbc_summary$n_measurements, 0.25),
      quantile(wbc_summary$n_measurements, 0.75)
    ),
    sprintf(
      "%.2f (%.2f–%.2f)",
      median(wbc_summary$followup_years),
      quantile(wbc_summary$followup_years, 0.25),
      quantile(wbc_summary$followup_years, 0.75)
    ),
    sprintf(
      "%.2f (%.2f–%.2f)",
      median(wbc_summary$variation_WBC, na.rm = TRUE),
      quantile(wbc_summary$variation_WBC, 0.25, na.rm = TRUE),
      quantile(wbc_summary$variation_WBC, 0.75, na.rm = TRUE)
    ),
    sprintf(
      "%.2f (%.2f–%.2f)",
      median(wbc_summary$median_WBC),
      quantile(wbc_summary$median_WBC, 0.25),
      quantile(wbc_summary$median_WBC, 0.75)
    )
  )
)

summary_table

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
saveRDS(temp_short, file = paste0("WBC_median_df.rds"))
