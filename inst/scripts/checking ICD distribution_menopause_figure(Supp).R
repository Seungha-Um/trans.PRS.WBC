# data frame for WBC median 
rm(list=ls())
library(data.table)
library(dplyr)
library(ggsci)
library(ggpubr)
library(ggsci)

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
              distinct(EMERGEID, .keep_all = TRUE) 

df_race <- fread("/manifest.csv"); 
df_race <- df_race %>% rename(EMERGEID = IID) 
df_race <- df_race %>% select(EMERGEID, sex2)

df_fig <- temp_short %>% select(- c(Leukocyte, logic)) %>% 
          left_join(df_race, by = "EMERGEID") %>% filter(sex2 != "unknown")

p1 <- df_fig %>% ggplot(aes(x = AGE)) +
          geom_histogram(aes(fill = sex2, color = sex2),
                          bins = 30, color = "black", 
                          linewidth = 0.3, alpha =0.8) + 
          scale_fill_npg() +
          theme_bw() + theme(legend.position = "bottom") + 
          labs(fill = NULL)

# To creat ICD distribution
ori_ICD <- readRDS("~/ICD.rds")
icd <- fread("/ICD_GWAS.csv", nThread = parallel::detectCores())
#head(icd); tail(icd); dim(icd) # 20333657

ex_ICD <- fread(file = "~/exclusion_ICD.csv")
ex_ICD <- ex_ICD[-c(1:2),]
colnames(ex_ICD) <- c("Type", "ICD_CODE", "Descr")  

df_chemo <- ori_ICD %>% filter(ICD_CODE %in% unique(ex_ICD$ICD_CODE)) 
ICD_clean <- ori_ICD %>% filter(!(SUBJID %in% unique(df_chemo$SUBJID)))

ICD_clean <- ICD_clean %>% filter(AGE_AT_EVENT != ".") %>% 
              mutate(AGE_AT_EVENT = as.numeric(AGE_AT_EVENT)) %>% 
              filter(!is.na(AGE_AT_EVENT))

# individuals with WBC
ind_WBC <- df_fig %>% pull(EMERGEID)
ICD_WBC <- ICD_clean %>% filter(SUBJID %in% ind_WBC) %>% select(-site)
df_fig2 <- df_fig %>% select(EMERGEID, sex2) #df from lab test, need it for gender

ICD_WBC <- ICD_WBC %>% rename(EMERGEID = SUBJID) %>% left_join(df_fig2, by = "EMERGEID") 

ICD_count1 <- ICD_WBC %>%
            group_by(ICD_CODE) %>%
            summarize(count = n_distinct(EMERGEID)) %>%
            ungroup()

ICD_count1 %>% arrange(desc(count))

ICD_count <- ICD_WBC %>% group_by(ICD_CODE, sex2) %>%
            summarize(count = n_distinct(EMERGEID)) %>%
            ungroup()

# ICD count distn with lab result
ICD_count %>% filter(sex2 =="female") %>% arrange(desc(count)) %>% select(-sex2)
ICD_count %>% filter(sex2 =="male") %>% arrange(desc(count)) %>% select(-sex2)

# Define menopause group
ICD9_mena <- unique(c("627.1", "627.2", "627.3", "627.4", "627.8", "627.9", 
                          "V49.81", "256.31", "256.39"))

ICD10_mena <- unique(c("N95.0", "N95.1", "N95.2", "N95.3",	"N95.8", "N95.9", "Z78.0", "E28.310", 
                          "E28.319", "E28.39", "E89.40", "E89.41"))
  
ICD_menopause <- c(ICD9_mena, ICD10_mena)  
ID_menopause <- ICD_WBC %>% filter(ICD_CODE %in% ICD_menopause) %>% pull(EMERGEID)

temp <- ICD_WBC %>% mutate(menopause = case_when(EMERGEID %in% ID_menopause ~ "menopause",
                                      TRUE ~"non-menopause"))

# this figure represents the density plot with participants with lab test 
# but with age at ICD code observation not age at the lab test

p2 <- temp %>% filter(sex2 == "female") %>% ggplot() + 
      geom_density(aes(x = AGE_AT_EVENT, color=menopause, fill=menopause), alpha= 0.6) + 
      theme(legend.position = "bottom")

########## pop is restricted with having lab test results.
# let's define menopause group more rigoursly. 
df_menopause <- ICD_WBC %>% filter(ICD_CODE %in% ICD_menopause) %>% select(-sex2) %>% 
               distinct(EMERGEID, .keep_all = TRUE) 
                
df_female <- df_fig %>% filter(sex2 == "female") 
df_menopause2 <- df_menopause %>% filter(EMERGEID %in% df_female$EMERGEID)

dim(df_menopause2)
dim(df_female)

df_menopause_WBC <- df_female %>% 
                   left_join(df_menopause2, by = "EMERGEID") %>%
                   mutate(menopause = case_when(EMERGEID %in% ID_menopause ~ "menopause ICD = Yes",
                                                  TRUE ~"menopause ICD = NO")) 

p11 <- df_menopause_WBC %>% ggplot(aes(x = AGE)) +
            geom_histogram(aes(fill = menopause, color = menopause),
                           bins = 30, color = "black", 
                           linewidth = 0.3, alpha =0.8) + 
            scale_fill_npg() +
            theme_bw() + theme(legend.position = "bottom") + 
            labs(fill = NULL, y="Count", x="Age")

pal <- pal_npg("nrc")(10)
#scales::show_col(pal)
p22 <- ggplot(df_menopause_WBC,
               aes(x = AGE,
                   y = subj_median,
                   color = menopause,
                   group = menopause)) +
        geom_smooth(method = "loess", se = FALSE) +
        geom_rug(data = df_menopause_WBC %>% filter(menopause == "menopause ICD = Yes"),
          aes(x = AGE), inherit.aes = FALSE, sides = "b", color = pal[5], alpha=0.8) +
        scale_color_npg() + theme_bw() +
        theme(legend.position = "bottom") + 
        labs(x = "Age", y = "Median WBC", color = NULL) + 
        coord_cartesian(ylim = c(6.7, 9))

plot2 <- ggarrange(p11, p22)
