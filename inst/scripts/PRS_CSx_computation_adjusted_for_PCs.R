rm(list=ls())
library(purrr)
library(dplyr)
library(data.table)
library(ggplot2)
library(data.table)
library(dplyr)
library(class)

# For post-hoc adjustment, we also need to compute PRS of 1KG.
# Load PCA data
pca_ref <- readRDS("/PCA/pca_ref.rds")   #reference 1KG PCA
snp_load <- readRDS("/PCA/snp_load.rds")   #loading from 1KG
samp_load <- readRDS("/samp_load.rds")   #projected eMERGE on 1KG
select <- dplyr::select

# option 1 for PRS_CS (sperated PRS for each ancestry group)
# method <- "PRS_CS/"
# anc_group_vec <- c("afr", "EA", "eas")
# CS_logic = TRUE

# option 2 for PRS_CSx
method <- "PRS_CSx/"
anc_group_vec <- c("AFR", "EAS", "EUR")
CS_logic = FALSE # is this PRS-CS = FALSE 

emerge_list <- KG_list <- list()
for(t in anc_group_vec){
  anc_group <- t
  # raw PRS of 1KG
  if(CS_logic == TRUE){
    data_use <- "1KG"
    setwd(paste0("/data/", method, data_use))
    prefix1 <- "CS_sep"
    prefix2 <- "_1KG_"
  }else{
    data_use <- "1KG" 
    setwd(paste0("/data/", method, data_use))
    prefix1 <- "CSx_"
    prefix2 <- "_1KG_"
  }
  
  prs_list <- list()
  for(j in 1:22){
    file_name <- paste0(prefix1, anc_group, prefix2, j, ".csv")
    prs_list[[j]] <- read.csv(file_name, header = TRUE) 
    colnames(prs_list[[j]])[2] <- paste0("PRS", j)
  }
  
  joined_prs <- reduce(prs_list, left_join, by = "ID")
  joined_prs$PRS <- rowSums(joined_prs %>% dplyr::select(-ID) )
  joined_prs <- joined_prs %>% rename(IID = ID) %>% select(IID, PRS)
  
  KG_list[[anc_group]] <- joined_prs
  
  # raw PRS of eMERGE
  data_use <- "EMERGE" 
  setwd(paste0("/data/", method, data_use))
  
  if(CS_logic == TRUE){
    prefix1 <- "CS_sep"
    prefix2 <- "_emerge_chr"
  }else{
    prefix1 <- "CS_sep_"
    prefix2 <- "_emerge_chr_"
  }  
  
  prs_list_emerge <- list()
  for(j in 1:22){
    file_name <- paste0(prefix1, anc_group, prefix2, j, ".csv")
    prs_list_emerge[[j]] <- read.csv(file_name, header = TRUE) 
    colnames(prs_list_emerge[[j]])[2] <- paste0("PRS", j)
  }
  
  joined_prs_emerge <- reduce(prs_list_emerge, left_join, by = "ID")
  joined_prs_emerge$PRS <- rowSums(joined_prs_emerge %>% dplyr::select(-ID) )
  joined_prs_emerge <- joined_prs_emerge %>% rename(IID = ID) %>% select(IID, PRS)
  
  emerge_list[[anc_group]] <- joined_prs_emerge
}

emerge_prs_raw <- Reduce(function(x, y) full_join(x, y, by = "IID"), emerge_list)
names(emerge_prs_raw)[-1] <- names(emerge_list)
emerge_prs_raw <- emerge_prs_raw %>% mutate(AFR_st = scale(AFR),
                                            EAS_st = scale(EAS),
                                            EUR_st = scale(EUR))

KG_prs_raw <- Reduce(function(x, y) full_join(x, y, by = "IID"), KG_list)
names(KG_prs_raw)[-1] <- names(KG_list)

KG_prs_raw <- KG_prs_raw %>% mutate(AFR_st = scale(AFR),
                                    EAS_st = scale(EAS),
                                    EUR_st = scale(EUR))

df_WBC <- fread("/Leukocyte.csv", fill=Inf)
df_WBC <- df_WBC %>% mutate(Leukocyte = case_when(Unit_Leuk == "/CU M" ~ Leukocyte/1000,
                                                  Unit_Leuk == ".CMM" ~ Leukocyte/1000,
                                                  TRUE ~ Leukocyte))
temp <- df_WBC %>% filter(!(is.na(Leukocyte))) 
temp <- temp %>% filter(!(Unit_Leuk %in% c(".HPF", ".hpf", "hpf")))
temp <- temp %>% filter(AGE >= 18 & AGE < 90)

df_ID_leuk <- readRDS("id_ex_wbc.rds") 
temp <- temp %>% filter(!(EMERGEID %in% df_ID_leuk) )
temp <- temp %>% filter(Leukocyte < 200)

temp_clean <- temp %>% group_by(EMERGEID) %>%
              dplyr::mutate(subj_median = median(Leukocyte, na.rm = TRUE)) %>%
              ungroup()

temp_short <- temp_clean %>% select(EMERGEID, subj_median) %>% 
              distinct(EMERGEID, .keep_all = TRUE) %>%  #keep only one row for each patient to work with median WBC
              rename(IID = EMERGEID)

# Filtered by WBC IDs 
df_PRS_lm <- emerge_prs_raw %>% filter(IID %in% temp_short$IID) %>%
             left_join(temp_short, by = "IID") 

# compute PRS-CSx
cov <- Hmisc::Cs(AFR_st, EUR_st, EAS_st)
formula <- reformulate(cov, response = "subj_median")
fit_WBC <- lm(formula, data = df_PRS_lm)
fitted <- as.matrix((emerge_prs_raw %>% select(AFR_st, EUR_st, EAS_st))) %*% coef(fit_WBC)[-1] 
emerge_prs_raw <- emerge_prs_raw %>% 
                  mutate(PRS_CSx = fitted)

fitted_KG <- as.matrix((KG_prs_raw %>% select(AFR_st, EUR_st, EAS_st))) %*% coef(fit_WBC)[-1] 
KG_prs_raw <- KG_prs_raw %>% mutate(PRS_CSx = fitted_KG)

## adjustment for PCs 
pc_num <- 5
PC_1KG <- data.frame(IID = pca_ref$sample.id, pca_ref$eigenvect[,1:pc_num]) 
colnames(PC_1KG)[-1] <- paste0("PC", 1:pc_num)
df_1KG <- KG_prs_raw %>% filter(IID %in% PC_1KG$IID) %>% left_join(PC_1KG, by = "IID")

cov <- paste0("PC", 1:pc_num)  
formula <- reformulate(cov, response = "PRS_CSx")
fit_mu <- lm(formula, data = df_1KG)
df_1KG <- df_1KG %>% mutate(mu_hat = predict(fit_mu, newdata = df_1KG),
                            resid_sq = resid(fit_mu)^2)

formula_var <- reformulate(cov, response = "resid_sq")
fit_var <- lm(formula_var, data = df_1KG)
df_1KG <- df_1KG %>% mutate(sig_hat = predict(fit_var, newdata = df_1KG))
df_1KG <- df_1KG %>% mutate(adj_PRS = (PRS_CSx - mu_hat) / sqrt(sig_hat))
  
df_1KG$log_resid2  <- log(pmax(df_1KG$resid_sq, .Machine$double.eps))
formula_var        <- reformulate(cov, response = "log_resid2")
fit_var2           <- lm(formula_var, data = df_1KG)
df_1KG <- df_1KG %>% mutate(sig_hat_log = exp(0.5 * predict(fit_var2, newdata=df_1KG)))
df_1KG <- df_1KG %>% mutate(adj_PRS_log = (PRS_CSx - mu_hat) / sig_hat_log)

# 1KG
KG_panel <- read.table("ALL.panel",
                       header = TRUE, stringsAsFactors = FALSE)
KG_panel2 <- KG_panel %>% rename(IID = sample)
df_1KG2 <- df_1KG %>% left_join(KG_panel2, by = "IID")

## EMERGE (validation data)
PC_emerge <- data.frame(IID = samp_load$sample.id, samp_load$eigenvect[,1:pc_num]) 
colnames(PC_emerge)[-1] <- paste0("PC", 1:pc_num)
emerge_prs_raw$IID <- as.character(emerge_prs_raw$IID)
df_emerge <- emerge_prs_raw %>% filter(IID %in% PC_emerge$IID) %>% left_join(PC_emerge, by = "IID")
df_emerge <- df_emerge %>% mutate(mu_hat = predict(fit_mu, newdata = df_emerge),
                                  sig_hat = predict(fit_var, newdata = df_emerge),
                                  sig_hat2 = exp(0.5 * predict(fit_var2, newdata = df_emerge)))

df_emerge <- df_emerge  %>% mutate(adj_PRS = (PRS_CSx - mu_hat) / sqrt(sig_hat),
                                   adj_PRS_log = (PRS_CSx - mu_hat) / sig_hat2)

# saveRDS(df_emerge, file = paste0("/PRS_emerge/PRS_CSx_new.rds") )

# density plot 
KG_panel <- read.table("ALL.panel", header = TRUE, stringsAsFactors = FALSE)

panel_sub <- KG_panel %>% filter(sample %in% pca_ref$sample.id) %>% 
             rename(sample.id = sample)

emerge_pop <- knn(train = pca_ref$eigenvect[, 1:2],
                  test = samp_load$eigenvect[, 1:2],
                  cl = panel_sub$super_pop,
                  k = 5)

knn_df <- data.frame(IID = samp_load$sample.id, race = emerge_pop)
df_emerge <- df_emerge %>% left_join(knn_df, by= "IID")

df_emerge %>% ggplot() + geom_point(aes(x=PC1, y=PC2, color = race))

foo1 <- df_emerge %>% select(PRS_CSx, race) %>% mutate(group = "Raw PRS") %>%
        rename(PRS = PRS_CSx)

foo2 <- df_emerge %>% select(adj_PRS_log, race) %>% mutate(group = "Adjusted PRS") %>%
        rename(PRS = adj_PRS_log)

foo2 <- df_emerge %>% select(adj_PRS, race) %>% mutate(group = "Adjusted PRS") %>%
  rename(PRS = adj_PRS)


df_graph <- rbind(foo1, foo2)
df_graph$group <- factor(df_graph$group, levels = c("Raw PRS", "Adjusted PRS"))
df_graph$race <- factor(df_graph$race, levels = c("EUR", "AFR", "EAS", "SAS", "AMR"))

library(ggsci)
p1 <- df_graph %>% ggplot() + geom_density(aes(x = PRS, 
                                          color = race, fill = race), alpha = 0.2) + 
                          theme_bw() +
                          theme(legend.position = "bottom",
                                legend.title = element_blank(), 
                                legend.box.spacing = unit(-2, "pt")) +
                          labs(y = NULL, x = NULL) +
                          # scale_color_d3() +
                          # scale_fill_d3() +
                          facet_wrap(~ group, scales = "free")

ggsave("/density_PRS.pdf", p1, device = cairo_pdf, width = 6, height = 3)

