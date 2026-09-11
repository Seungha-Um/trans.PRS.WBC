rm(list = ls())
library(data.table)
library(dplyr)
library(ggpubr)
library(BART)
library(Hmisc)

setwd("~/Documents/PRS_CR/data")
df_demo <- fread("person_GWAS.csv")
df_WBC <- fread("Leukocyte-EMERGE-2020-05-13.csv", fill=Inf)
df_ID_leuk <- readRDS("/id_ex_wbc.rds") 
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

df_race <- fread("/manifest.csv"); 
df_race <- df_race %>% rename(EMERGEID = IID) 
pc_num <- 10
df_race <- df_race %>% select(EMERGEID, sex2, site, num_range("PC", 1:pc_num))

df_PRS <- readRDS("/PRS_CSx_new.rds")
df_PRS <- df_PRS %>% mutate(PRS_CSx_st = scale(adj_PRS))
df_PRS <- df_PRS %>% select(PRS_CSx_st, IID) %>%
              dplyr::rename(EMERGEID = IID) %>%
              mutate(EMERGEID = as.integer(EMERGEID))

df_bmi_ori <- fread("/bmi_GWAS.csv") #AGE and number of repeatation does not match, 

df_bmi <- df_bmi_ori %>% filter(MEASUREMENT_CONCEPT_ID == "3038553") %>% 
          select(SUBJID, VALUE_AS_NUMBER) %>%
          rename(EMERGEID = SUBJID) %>%
          mutate(VALUE_AS_NUMBER = as.numeric(VALUE_AS_NUMBER)) %>% 
          filter(!is.na(VALUE_AS_NUMBER))

df_bmi <- df_bmi %>% group_by(EMERGEID) %>% 
          summarize(BMI_med = median(VALUE_AS_NUMBER, rm.na=TRUE))

df_bmi <- df_bmi %>% filter(BMI_med < 100)

df_bmi %>% ggplot() + geom_density(aes(x = BMI_med)) 

# Option2: consider median of WBC for each patient
temp <- df_WBC %>%
        select(EMERGEID, AGE, Leukocyte) %>% 
        group_by(EMERGEID) %>%
        mutate(
          subj_median = Leukocyte[which.min(abs(Leukocyte - median(Leukocyte, na.rm = TRUE)))]
        ) %>%
        ungroup() %>%
        mutate(logic = Leukocyte == subj_median)

temp <- temp %>% filter(logic == TRUE) %>% arrange(EMERGEID, Leukocyte, AGE) %>% 
        distinct(EMERGEID, .keep_all = TRUE) 

temp_clean <- temp %>% left_join(df_PRS, by = "EMERGEID")
temp_clean <- temp_clean %>% left_join(df_race, by = "EMERGEID")
temp_clean <- temp_clean %>% left_join(df_bmi, by = "EMERGEID")
temp_clean <- temp_clean %>% filter(!is.na(PRS_CSx_st))
temp_clean <- temp_clean %>% filter(!is.na(AGE))
temp_clean <- temp_clean %>% filter(!is.na(BMI_med))
temp_clean <- temp_clean %>% filter(!(site %in% c("ccmc", "chop")))

temp_clean <- temp_clean %>% mutate(BMI_ind = case_when(BMI_med < 30 ~ "Low",
                                              BMI_med >= 30 ~ "High"))

X_mat <- temp_clean %>% ungroup() %>% select(PRS_CSx_st, num_range("PC", 1:pc_num), 
                                             sex2, AGE, BMI_med) 

X_mat <- X_mat %>% mutate(PRS_CSx_st = as.numeric(PRS_CSx_st))
Y <- temp_clean$Leukocyte
X_mat_test <- X_mat %>% select(-AGE)

# Prepare the test set with equally spaced age intervals
age_seq <- seq(round(min(X_mat$AGE)), round(max(X_mat$AGE)), 5)
X_list <- lapply(age_seq, function(a){cbind(X_mat_test, AGE = a)})
X_age_mat <- do.call(rbind, X_list)
X_age_mat <- data.matrix(X_age_mat)
X_mat <- data.matrix(X_mat)

#fit BART model 
fit <- wbart(x.train = X_mat, y.train = Y, x.test = X_age_mat, 
             ndpost=2500L, nskip=2500L)

fitted_train <- fit$yhat.train.mean
fitted_test <- fit$yhat.test.mean

X_pred <- as.data.frame(X_mat) %>% mutate(fitted = fitted_train, Y = Y)
X_pred <- X_pred %>% mutate(decile = ntile(PRS_CSx_st, 10))
X_pred %>% ggplot() + geom_line(aes(x = AGE, y = fitted, color = as.factor(decile)))

# figures with predicted values 
X_pred %>%
  ggplot(aes(x = AGE, y = fitted, color = as.factor(decile))) +
  geom_smooth(se = FALSE, method = "loess", span = 0.7) +
  theme_minimal()

X_pred %>%
  ggplot(aes(x = AGE, y = fitted, color = as.factor(decile))) +
  geom_smooth(se = FALSE, method = "loess", span = 0.7) + 
  theme(legend.position = "bottom") + 
  facet_grid(~ sex2)

X_pred %>%
  ggplot(aes(x = AGE, y = fitted, color = as.factor(decile))) +
  geom_smooth(se = FALSE, method = "loess", span = 0.7) + 
  theme(legend.position = "bottom") + 
  facet_grid(~ BMI_ind)

X_pred2 <- as.data.frame(X_age_mat) %>% mutate(fitted = fitted_test,
                                               decile = ntile(PRS_CSx_st, 10))  

X_pred2 <- X_pred2 %>% group_by(AGE, decile) %>% summarize(fitted_mean = mean(fitted))

X_pred2 %>% ggplot() + 
            geom_line(aes(x = AGE, y = fitted_mean, 
                          color = as.factor(decile)), linewidth=0.6) + 
            labs(color = "decile") +
            ylab("WBC count") + 
            theme_bw() +
            theme( legend.position = "top") 

# Representative decision tree with split rules and fitted values
library(rpart)

X_rpart <- as.data.frame(X_mat)
X_rpart <- X_rpart %>% rename(PRS = PRS_CSx_st) %>%
            rename(age = AGE) %>% 
            rename(BMI = BMI_med)
X_rpart <- X_rpart %>% mutate(sex = case_when(sex2 == 1 ~ "female",
                              TRUE ~ "male")) %>% select(-sex2) 

tree_fit <- rpart(fit$yhat.train.mean ~ ., data = X_rpart,
                  control = rpart.control(cp = 0.005, maxdepth = 4))

tree_fit$variable.importance
colMeans(fit$varcount)

stat_df <- t(rbind(fit$varcount %>% apply(2, mean), fit$varcount %>% apply(2, quantile, probs = c(0.025, 0.975))))
stat_df

cor(fit$yhat.train.mean, Y)

library(rpart.plot)

rpart.plot(tree_fit, box.palette = "Reds")
rpart.plot(tree_fit, box.palette = "Greens")
rpart.plot(tree_fit, box.palette = "Blues")

setEPS()                     # tell R to create an EPS device
postscript("tree_plot.eps",width = 5, height = 4)
rpart.plot(tree_fit)
dev.off()

temp_test2 <- temp_test %>% ungroup() %>%
                filter(!is.na(PRS_CSx_st)) %>%
                mutate(PRS_decile = ntile(PRS_CSx_st, 10)) 

temp_clean2 <- temp_test2 %>%
               group_by(PRS_decile, age) %>%
               mutate(fitted_median = median(fitted_train))

temp_clean2 %>% ggplot() + geom_line(aes(x = age, y = fitted_median, 
                                         color = as.factor(PRS_decile)))

temp_clean2 <- temp_clean %>% ungroup() %>%
              filter(!is.na(PRS_CSx_st)) %>%
              mutate(PRS_decile = ntile(PRS_CSx_st, 10)) 

temp_clean2 <- temp_clean2 %>%
                group_by(PRS_decile, AGE) %>%
                mutate(fitted_median = median(fitted_train))

temp_clean2 %>% ggplot() + geom_line(aes(x = AGE, y = fitted_median, 
                                         color = as.factor(PRS_decile)))

p1 <- temp_clean2 %>%
      ggplot(aes(x = AGE, y = fitted_median, color = as.factor(PRS_decile))) +
      geom_point(alpha = 0.3) +          # optional scatter points
      #geom_smooth(se = FALSE, method = "loess", span = 0.6, size = 1) + 
      theme( legend.position = "bottom") +
      scale_color_discrete(name = "PRS Decile") +
      ylab("predicted WBC")

p2 <- temp_clean2 %>%
      ggplot(aes(x = AGE, y = fitted_median, color = as.factor(PRS_decile))) +
      #geom_point(alpha = 0.3) +          # optional scatter points
      geom_smooth(se = FALSE, method = "loess", span = 0.6, size = 1) + 
      theme(legend.position = "bottom") +
      scale_color_discrete(name = "PRS Decile") +
      ylab("predicted WBC")

p3 <- temp_clean2 %>% filter(AGE > 18 & AGE < 90) %>%
      ggplot(aes(x = AGE, y = fitted_median, color = as.factor(PRS_decile))) +
      #geom_point(alpha = 0.3) +          # optional scatter points
      geom_smooth(se = FALSE, method = "loess", span = 0.6, size = 1) + 
      scale_color_discrete(name = "PRS Decile") +
      theme(legend.position = "bottom") + ylab("predicted WBC")

ggsave("regression.pdf", plot = p2, device = cairo_pdf, width = 8, height = 4)
ggsave("regression2.pdf", plot = p3, device = cairo_pdf, width = 8, height = 4)


temp_clean2 %>% ggplot() + geom_boxplot(aes(x = as.factor(PRS_decile), y = fitted_median, 
                              color = as.factor(PRS_decile)),outlier.shape = NA) + 
                coord_cartesian(ylim = c(3, 10)) +
                theme(legend.position = "bottom") + 
                facet_grid(~ age_grp)

