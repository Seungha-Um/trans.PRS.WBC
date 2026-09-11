rm(list = ls())
library(data.table)
library(dplyr)
library(ggpubr)
library(BART)
library(SoftBart)
library(scales)
library(stringr)
library(forcats)

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
df_WBC <- df_WBC %>% group_by(EMERGEID) %>% mutate(U_age = max(AGE),
                                                   L_age = min(AGE))

WBC_range <- df_WBC %>% select(EMERGEID, U_age, L_age) %>% distinct()
# Repeated measurements within a week which including 0 values
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

# Exclude measurements during pregnancy
ori_ICD <- fread("/ICD_GWAS.csv")
# Exclude observations from children hospitals
# Exclude infants who have mother's record
df_pre <- ori_ICD %>% filter(!(site %in% c("ccmc", "chop", "cchmc"))) %>%
  filter(AGE_AT_EVENT != ".") %>%
  mutate(AGE_AT_EVENT = as.numeric(AGE_AT_EVENT)) %>%
  rename(EMERGEID = SUBJID) %>%
  filter(EMERGEID %in% unique(df_WBC$EMERGEID)) %>%
  filter(AGE_AT_EVENT != ".") %>%
  filter(AGE_AT_EVENT >= 18) %>%
  filter((ICD_FLAG == 10 & grepl("^(O|Z33|Z34|Z3A)", ICD_CODE)) |
           (ICD_FLAG == 9  & grepl("^(63[0-9]|64[0-9]|65[0-9]|66[0-9]|67[0-9]|V22|V23)", ICD_CODE)))

df_pre <- df_pre %>% left_join(WBC_range, by = "EMERGEID") %>%
  filter((AGE_AT_EVENT <= U_age & AGE_AT_EVENT >= L_age))

#“Pregnancy was identified using ICD-10 codes O*, Z33*, Z34*, and Z3A*, and ICD-9 codes 630–679 and V22–V23.”
df_pre <- df_pre %>% mutate(AGE_AT_EVENT = as.numeric(AGE_AT_EVENT)) %>%
  filter(AGE_AT_EVENT < 60) %>%
  group_by(EMERGEID) %>%
  summarize(UB_pre = max(AGE_AT_EVENT, na.rm = TRUE) + 280/365.25,
            LB_pre = min(AGE_AT_EVENT, na.rm = TRUE) - 280/365.25) %>%
  ungroup()

df_WBC <- df_WBC %>% left_join(df_pre, by = "EMERGEID") %>%
  mutate(preg_ix = case_when(AGE >= LB_pre & AGE <= UB_pre ~ "pregnant",
                             TRUE ~ "non pregnant")) %>%
  filter(preg_ix != "pregnant")

df_race <- fread("/manifest.csv");
df_race <- df_race %>% rename(EMERGEID = IID)
df_race <- df_race %>% select(EMERGEID, sex2, site)

#load PCA data
samp_load <- readRDS("/samp_load.rds")   #projected eMERGE on 1KG
pc_num <- 10
df_pca <- data.frame(EMERGEID = as.integer(samp_load$sample.id), samp_load$eigenvect[,1:pc_num])
colnames(df_pca)[-1] <- paste0("PC", 1:pc_num)


#df_PRS <- readRDS("/Users/ums3/Documents/PRS_CR/data/PRS_emerge/PRS_CSx_new.rds")
df_PRS <- readRDS("/PRS_CSx_after_adj.rds")

df_PRS <- df_PRS %>% mutate(PRS_CSx_st = scale(adj_PRS_log))
df_PRS <- df_PRS %>% select(PRS_CSx_st, IID) %>%
  dplyr::rename(EMERGEID = IID) %>%
  mutate(EMERGEID = as.integer(EMERGEID))



# BMI
df_bmi_ori <- fread("/EMERGE_201907_BMI_GWAS.csv") #AGE and number of repetition does not match with WBC df
df_bmi_ori <- df_bmi_ori %>% filter(!is.na(BODY_MASS_INDEX)) %>%
  rename(EMERGEID = SUBJECT_ID)
df_bmi_ori <- df_bmi_ori %>% filter(EMERGEID %in% unique(df_WBC$EMERGEID))
no_wei_in <- df_bmi_ori %>% filter(is.na(WEIGHT)) %>% pull(EMERGEID)
no_hei_in <- df_bmi_ori %>% filter(is.na(HEIGHT))%>% pull(EMERGEID)
only_bmi_in <- unique(c(no_hei_in, no_wei_in))
only_bmi_df <- df_bmi_ori %>% filter(EMERGEID %in% only_bmi_in) %>%
  arrange(EMERGEID) %>%
  select(EMERGEID, BODY_MASS_INDEX)

# Define ranges to detect outliers for BMI
temp <- only_bmi_df %>% group_by(EMERGEID) %>%
  mutate(BMI_med = median(BODY_MASS_INDEX, na.rm = TRUE),
         BMI_sd  = sd(BODY_MASS_INDEX, na.rm = TRUE),
         outlier = case_when(is.na(BMI_sd) ~ "non-outlier",
                             abs(BODY_MASS_INDEX - BMI_med) > 5 * BMI_sd ~ "outlier",
                             TRUE ~ "non-outlier")) %>% ungroup()

out_df <- temp %>% filter(outlier == "outlier") %>% select(EMERGEID, BODY_MASS_INDEX)
df_bmi_clean <- df_bmi_ori %>% anti_join(out_df, by = c("EMERGEID", "BODY_MASS_INDEX"))
df_bmi_clean <- df_bmi_clean %>% filter(BODY_MASS_INDEX < 100)

# Compute median BMI within the age window corresponding to WBC measurements
df_bmi_clean <- df_bmi_clean %>% left_join(WBC_range, by = "EMERGEID") %>%
  filter((AGE_AT_OBSERVATION <= U_age & AGE_AT_OBSERVATION >= L_age))


# Exclude BMI during pregnancy
df_bmi_clean <- df_bmi_clean %>% left_join(df_pre, by = "EMERGEID") %>%
  mutate(preg_ix = case_when(AGE_AT_OBSERVATION >= LB_pre & AGE_AT_OBSERVATION <= UB_pre ~ "pregnant",
                             TRUE ~ "non pregnant")) %>%
  filter(preg_ix != "pregnant")

# Compute median of BMI with clean BMI dataset
df_bmi_clean <- df_bmi_clean %>% group_by(EMERGEID) %>%
  summarize(BMI_med = median(BODY_MASS_INDEX, na.rm = TRUE))

df_bmi <- df_bmi_clean %>% mutate(BMI_cat = case_when(BMI_med >= 30 ~ "High",
                                                      TRUE ~ "Low")) %>% select(-BMI_med)

# WBC
temp <- df_WBC %>%
  group_by(EMERGEID) %>%
  summarize(WBC_med = median(Leukocyte),
            Age_med = median(AGE)) %>%
  ungroup()
temp <- temp %>% rename(AGE = Age_med)

temp_clean <- temp %>% left_join(df_PRS, by = "EMERGEID")
temp_clean <- temp_clean %>% left_join(df_race, by = "EMERGEID")
temp_clean <- temp_clean %>% left_join(df_bmi, by = "EMERGEID")
temp_clean <- temp_clean %>% left_join(df_pca, by = "EMERGEID")
temp_clean <- temp_clean %>% filter(!is.na(PRS_CSx_st))
temp_clean <- temp_clean %>% filter(!is.na(AGE))
temp_clean <- temp_clean %>% filter(!is.na(BMI_cat))
#temp_clean <- temp_clean %>% filter(!is.na(BMI_med))
temp_clean <- temp_clean %>% filter(!(site %in% c("ccmc", "chop", "cchmc")))

X_mat <- temp_clean %>% ungroup() %>%
  select(PRS_CSx_st, num_range("PC", 1:10), sex2, AGE, BMI_cat)

X_mat <- X_mat %>% mutate(PRS_CSx_st = as.numeric(PRS_CSx_st))

X_mat_site <- temp_clean %>% ungroup() %>% pull(site)
num_site <- length(unique(X_mat_site))
name_site <- unique(X_mat_site)
factor_site <- factor(X_mat_site, levels = name_site)

########################## Implement riBART (BART with random intercept)

library(zeallot)
Y <- temp_clean$WBC_med
Y_scaled <- scale(Y)
c(Y_dims, center_Y, scale_Y) %<-% attributes(Y_scaled)

X_mat <- X_mat %>%
  mutate(
    sex2 = if_else(sex2 == "female", 1, 0),
    BMI_cat = if_else(BMI_cat == "High", 1, 0)
  )

source("/normalized_funs.R")
df_train <- as.matrix(quantile_normalize_bart(X_mat))

###
PRS_seq <- as.data.frame(df_train) %>% mutate(decile_temp = ntile(PRS_CSx_st, 4)) %>%
            group_by(decile_temp) %>% summarize(Q_mean = mean(PRS_CSx_st)) %>% pull(Q_mean)

X_mat_test <- as.data.frame(df_train) %>% select(-PRS_CSx_st)

X_list <- lapply(PRS_seq, function(a){cbind(PRS_CSx_st = a, X_mat_test)})
X_prs_mat <- do.call(rbind, X_list)
df_test <- data.matrix(X_prs_mat)
###

library(SoftBart)
hypers <- Hypers(X = df_train, Y = Y_scaled, normalize_Y = FALSE)
num_tree <- 200
hypers$num_tree <- num_tree
opts <- Opts(update_s = FALSE)
forest <- MakeForest(hypers, opts, warn = FALSE)
#step_size <- 100
sim_num <- 5000
NN <- dim(X_mat)[1]
N_test <- dim(df_test)[1]

mu_hat <- matrix(NA, sim_num, NN)
mu_test <- matrix(NA, sim_num, N_test)
alpha_mat <- matrix(NA, sim_num, num_site)

#mu_test_hat <- matrix(NA, sim_num/step_size, df_train)
ave_mu_test <- list()
ave_mu_test_all <- list()

beta_mat <- rep(1, num_site)
ran_int <- model.matrix(~ factor_site - 1) %*% beta_mat

################### random intercept BART
Sigma_vec <- rep(NA, sim_num)
Z_beta <- list()
var_count <- matrix(NA, sim_num, dim(df_train)[2])

for(s in 1:sim_num){
  R = Y_scaled - ran_int
  mu_hat[s,] <- forest$do_gibbs(df_train, R, df_train, 1)
  mu_test[s,] <- forest$do_predict(df_test)
  Sigma <- forest$get_sigma()
  var_count[s,] <- forest$get_counts()
  Sigma_vec[s] <- Sigma
  delta <- Y_scaled - mu_hat[s,]
  alpha_vec <- rep(NA, num_site)
  ran_int_new <- rep(NA, NN)

  for(g in name_site){
    ind_sites <- X_mat_site == g
    num_obs <- sum(ind_sites)
    Cov <- Sigma^2/num_obs
    mu_vec <- sum(delta[ind_sites])/num_obs
    alpha_g <- rnorm(1, mu_vec, Cov)
    alpha_vec[which(g == name_site)] <- alpha_g
    ran_int_new[ind_sites] <- alpha_g
  }
  alpha_mat[s,] <- alpha_vec
  #if(s %% step_size == 0) mu_test_hat[(s/step_size),] <- EY_est
  ran_int <- ran_int_new
  if(s %% 100 == 0) cat("MCMC iter :", s, "\n")
}

bin_mat <- model.matrix(~ factor_site - 1)
rand_ind_mat <- alpha_mat %*% t(bin_mat)
temp_f <- mu_hat + rand_ind_mat
temp_f <- temp_f * scale_Y + center_Y

mu_test_hat <- mu_test * scale_Y + center_Y

ind_Y <- Y < 20
var_mat <- var_count[(sim_num/2):sim_num, ]
var_im <- var_mat / rowSums(var_mat)

##### variable importance
temp_var <- as.data.frame(rbind(mean = colMeans(var_im),
                                apply(var_im, 2, quantile, probs = c(0.025, 0.975))))

name_fac <- names(as.data.frame(df_train))
name_fac[name_fac=="PRS_CSx_st"] <- "PRS"
name_fac[name_fac=="sex2"] <- "Sex"
name_fac[name_fac=="BMI_cat"] <- "BMI"
colnames(temp_var ) <- name_fac

library(tidyr)
library(dplyr)
library(Polychrome)
library(ggsci)

temp_long <- temp_var %>%
              tibble::rownames_to_column("stat") %>%
              pivot_longer(-stat, names_to = "variable", values_to = "value") %>%
              pivot_wider(names_from = stat, values_from = value)

npg_cols <- pal_npg("nrc")(10)
extra_cols <- c("#8C8C8C", "#6B8FA3", "#DD8452","#55A868")
cols14 <- c(npg_cols, extra_cols)

p1 <- ggplot(temp_long, aes(x = reorder(variable, mean), y = mean, fill = variable)) +
  geom_col() +
  scale_color_manual(values = cols14) +
  scale_fill_manual(values = cols14) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.4, alpha=0.4) +
  coord_flip() +  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
  xlab("Variable Importance") + labs(y = NULL) + guides(fill = guide_legend(nrow = 2))

#ggsave("/Users/ums3/Desktop/VI.pdf", p1, width = 6, height = 4)

y_hat_train_mean <- colMeans(temp_f[(sim_num/2):sim_num,])

X_pred <- as.data.frame(X_mat) %>% mutate(fitted = y_hat_train_mean, Y = Y)
X_pred <- X_pred %>% mutate(decile = ntile(PRS_CSx_st, 4))
X_pred <- X_pred %>% mutate(decile2 = case_when(decile == 4 ~ "Q4(75–100%)",
                                                decile == 3 ~ "Q3(50–75%)",
                                                decile == 2 ~ "Q2(25–50%)",
                                                decile == 1 ~ "Q1(0–25%)"))

X_pred <- X_pred %>% mutate(BMI_cat3 = case_when(BMI_cat == 1 ~ "BMI>30 (High)",
                                                 BMI_cat == 0 ~ "BMI<30 (Low)"),
                            sex3 = case_when(sex2 == 1 ~ "Female",
                                             sex2 == 0 ~ "Male") )

X_pred$BMI_cat3 <- factor(X_pred$BMI_cat3, levels = c( "BMI<30 (Low)", "BMI>30 (High)"))

p2 <- X_pred %>% ggplot(aes(x = AGE, y = fitted, color = as.factor(decile2), linetype = sex3)) +
  geom_smooth(se = FALSE, method = "loess", span = 0.7) +
  geom_rug(
    aes(x = AGE),
    sides = "b",
    inherit.aes = FALSE,
    color = "gray30",
    alpha = 0.1
  ) +
  scale_color_npg() +
  scale_fill_npg() +
  facet_grid(cols = vars(BMI_cat3), switch = "y") +
  theme_bw() +
  theme(legend.position = "bottom") + xlab("Age") +
  labs(y = "WBC", color = "PRS Quartiles", linetype = "sex") +
  guides(linetype = guide_legend(override.aes = list(color = "gray35")))

X_rpart <- as.data.frame(X_mat)
X_rpart <- X_rpart %>% rename(PRS = PRS_CSx_st) %>%
  rename(age = AGE) %>%
  rename(BMI = BMI_cat)

X_rpart <- X_rpart %>% mutate(sex = case_when(sex2 == 1 ~ "female",
                                              TRUE ~ "male")) %>% select(-sex2)
library(rpart.plot)
tree_fit <- rpart(y_hat_train_mean ~ ., data = X_rpart,
                  control = rpart.control(cp = 0.005, maxdepth = 4))

npg_cols <- pal_npg("nrc")(10)
extra_cols <- c("#8C8C8C", "#6B8FA3", "#DD8452","#55A868")

base <- npg_cols[4]
my_pal <- colorRampPalette(c("white", base))
rpart.plot(tree_fit, box.palette = my_pal(30))   # 30 = smooth gradient

