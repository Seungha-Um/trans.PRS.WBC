rm(list = ls())
library(data.table)
library(dplyr)
library(tidyr)

# Load PRSs from multiple methods
df_EUR_Cs <- readRDS("PRS_CS_EA.rds") %>%
              rename(PRS_EUR = adj_PRS_log) %>% select(IID, PRS_EUR) 

df_Multi_Cs <- readRDS("PRS_CS_multi.rds") %>%
              rename(PRS_multi = adj_PRS_log) %>% select(IID, PRS_multi) 

df_Meta <- readRDS("PRS_CSx_META.rds") %>%
            select(IID, adj_PRS_log) %>% rename(PRS_META = adj_PRS_log) 

df_CSx <- readRDS( "PRS_CSx_after_adj.rds") %>%
            select(IID, adj_PRS_log) %>% rename(PRS_CSx = adj_PRS_log) 

df_PRS_all <- df_EUR_Cs %>% left_join(df_Multi_Cs, by = "IID") %>% 
              left_join(df_Meta, by = "IID") %>%
              left_join(df_CSx, by = "IID")

# Load pre-calculated median WBC counts for each subject
temp_short <- readRDS("WBC_median_df.rds")
df_WBC <- temp_short %>% select(IID, subj_median_st) %>% 
          mutate(IID = as.character(IID))  %>%
          rename(subj_median = subj_median_st) %>% as.data.frame()

df_race <- fread("manifest.csv"); 
df_race <- df_race %>% select(IID, self_reported_race, genetic_ancestry_race, ethnicity2, site)

# Exclude children's hospital observations
df_gene_race <- df_race %>% filter(!(site %in% c("ccmc", "chop"))) %>% 
                select(IID, genetic_ancestry_race, ethnicity2) %>%
                mutate(IID = as.character(IID))
CH_ID <- df_race %>% filter(site %in% c("ccmc", "chop")) %>% pull(IID)
df_PRS_all <- df_PRS_all %>% filter(!(IID %in% CH_ID)) 
df_PRS_all <- df_PRS_all %>% left_join(df_gene_race, by = "IID")
df_PRS_all <- df_PRS_all %>% filter(IID %in% df_WBC$IID) %>% left_join(df_WBC, by = "IID")

# R square computation between WBC and PRS 
rsq_func <- function(ancestry, data_df){
  ancestry <- ancestry
  multi_sum <- lm(subj_median ~ PRS_multi, data = data_df)
  rsq_multi <- summary(multi_sum)$r.squared
  
  csx_sum <- lm(subj_median ~ PRS_CSx, data = data_df)
  rsq_csx <- summary(csx_sum)$r.squared
  
  meta_sum <- lm(subj_median ~ PRS_META, data = data_df)
  rsq_meta <- summary(meta_sum)$r.squared
  
  eur_sum <- lm(subj_median ~ PRS_EUR, data = data_df)
  rsq_eur <- summary(eur_sum)$r.squared
  
  output = data.frame(group = ancestry, multi = rsq_multi, 
                      csx = rsq_csx, meta = rsq_meta, eur = rsq_eur)
  return(output)
}

anc_group <- unique(df_race$genetic_ancestry_race)
df_list <- list()

for(k in anc_group){
  print(k)
  temp_df <- df_PRS_all %>% filter(genetic_ancestry_race == k)    
  print(dim(temp_df))
  ind <- which(anc_group == k)
  df_list[[ind]] <- rsq_func(ancestry = k, data_df = temp_df)
}

df_rsq_all <- rsq_func(ancestry = "all", data_df = df_PRS_all)
df_hisp <- df_PRS_all %>% filter(ethnicity2 == "hispanic or latino")    
df_rsq_hisp <- rsq_func(ancestry = "Hispanic", data_df = df_hisp)
df_rsq_final <- rbind(df_rsq_all, do.call(rbind, df_list), df_rsq_hisp)

df_long <- pivot_longer(df_rsq_final, cols = -group, names_to = "method", values_to = "R2")
df_long$method <- factor(df_long$method, levels = c("csx","meta", "multi", "eur"))

df_long$method <- factor(
  df_long$method,
  levels = c("csx", "meta", "multi", "eur"),
  labels = c("PRS-CSx", "PRS-CS meta", "PRS-CS multi", "PRS-CS EUR")
)
df_long$group <- factor(df_long$group, levels = c("all", "european", "asian", "african", "Hispanic"))

p1 <- ggplot(df_long, aes(x = method, y = R2, fill = method)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_wrap(~ group, nrow = 1, scales = "free_y") + 
      labs(y = expression(R^2), x = NULL) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "bottom",
        legend.title = element_blank()
      )

#ggsave("comp.eps", plot=p1, device=cairo_ps, width = 8, height = 4)
################################################################################

df_emerge_WBC <- df_PRS_all %>% mutate(Median = subj_median)
df_emerge_WBC <- df_emerge_WBC %>% mutate(PRS_decile = as.factor(ntile(PRS_CSx, 10)))
df_no_outliers2 <- df_emerge_WBC %>%
                  group_by(PRS_decile) %>%
                  mutate(LB = quantile(Median, 0.025),
                         UB = quantile(Median, 0.975)) %>%
                  filter(Median >= LB,
                         Median <= UB) %>%
                  ungroup()

df_no_outliers2$genetic_ancestry_race <- factor(df_no_outliers2$genetic_ancestry_race, levels = c("european", "african", "asian"))

p1 <- df_no_outliers2 %>% ggplot(aes(x = PRS_decile, y = Median, color = PRS_decile)) +
      theme_bw() +
      geom_boxplot() + labs(x = "PRS Decile", y = "Median WBC") + 
      theme(legend.title = element_blank(), legend.position = "bottom") +
      guides(color = guide_legend(nrow = 1))
#ggsave("comp2.eps", plot=p1, device=cairo_ps, width = 8, height = 6)

p2 <- df_no_outliers2 %>% ggplot(aes(x = PRS_decile, y = Median, color = PRS_decile)) +
      theme_bw() +
      geom_boxplot() + labs(x = "PRS Decile", y = "Median WBC") + 
      theme(legend.title = element_blank(), legend.position = "bottom") +
      facet_wrap(~ genetic_ancestry_race, nrow = 1) + 
      guides(color = guide_legend(nrow = 1))

#ggsave("comp3.eps", plot=p2, device=cairo_ps, width = 8, height = 5)

df %>% ggpadj_PRS2df %>% ggplot(aes(x = PRS, col = race)) +
  geom_density(alpha = 0.3) +
  #xlim(-0.5,0.5) + 
  labs(
    x = "Polygenic Risk Score (PRS)",
    y = "Density", color = "Race", fill  = "Race"
  ) + theme_minimal()

cov <- paste0("PC.", 1:pc_num)  
formula <- reformulate(cov, response = "PRS")
fit <- lm(formula, data = df)

df <- df %>% mutate(mu_hat = predict(fit, newdata=df),
                    resid_sq = resid(fit)^2)

formula_var <- reformulate(cov, response = "resid_sq")
fit_var <- lm(formula_var, data = df)
df <- df %>% mutate(sig_pre= predict(fit_var, newdata=df))
df <- df %>% mutate(adj_PRS = (df$PRS - df$mu_hat) / sqrt(df$sig_pre))

df2 %>% ggplot(aes(x = adj_PRS, col = race, fill = race)) +
  geom_density(alpha = 0.2) +
  #xlim(-3, 3) + 
  labs(
    x = "Polygenic Risk Score (PRS)",
    y = "Density", color = "Race", fill  = "Race"
  ) + theme_minimal()

# 2) Variance model on log scale: log(resid^2) ~ PCs
df$log_resid2  <- log(pmax(df$resid_sq, .Machine$double.eps))
form_var       <- reformulate(cov, response = "log_resid2")
fit_var        <- lm(form_var, data = df)
# sigma_hat = exp(0.5 * E[log(resid^2)|PCs])
sigma_hat      <- exp(0.5 * predict(fit_var, newdata=df))
df$adj_PRS_log <- (df$PRS - df$mu_hat) / sigma_hat


df %>% ggplot(aes(x = adj_PRS_log, col = race, fill = race)) +
  geom_density(alpha = 0.2) +
  #xlim(-3, 3) + 
  labs(
    x = "Polygenic Risk Score (PRS)",
    y = "Density", color = "Race", fill  = "Race"
  ) + theme_minimal()


df1 <- df %>% mutate(group = "before ADJ", PRS_2 = PRS)
df2 <- df %>% mutate(group = "after ADJ", PRS_2 = adj_PRS_log)
dff <- rbind(df1, df2)
dff$group <- factor(dff$group, levels =c("before ADJ", "after ADJ"))

p1<- dff %>%
  ggplot(aes(x = PRS_2, col = race)) +
  geom_density(alpha = 0.2, linewidth=0.6) +
  #xlim(-3, 3) + 
  labs(
    x = "Polygenic Risk Score (PRS)",
    y = "Density", color = "Race", fill  = "Race"
  ) + theme_bw() + facet_wrap(~group, scales = "free") + 
  theme(legend.title = element_blank())


p2 <- dff %>% filter(group == "before ADJ") %>%
  ggplot(aes(x = PRS_2, col = race)) +
  geom_density(alpha = 0.2, linewidth=0.6) +
  #xlim(-3, 3) + 
  labs(
    x = "Polygenic Risk Score (PRS)",
    y = "Density", color = "Race", fill  = "Race"
  ) + theme_bw() + 
  theme(legend.title = element_blank())

ggsave("both.eps", plot=p1, device=cairo_ps, width = 7, height = 4)

