rm(list=ls())
library(data.table)
library(PheWAS)
library(dplyr)
library(logistf)
library(purrr)
library(dplyr)
library(purrr)
library(logistf)
library(ggrepel)
library(ggplot2)

group_by <- dplyr::group_by
mutate <- dplyr::mutate
rename <- dplyr::rename
summarize <- dplyr::summarize

df_PRS_all <- readRDS("PRS_CSx_new.rds") %>%
              select(IID, adj_PRS) %>% rename(PRS_CSx_st = adj_PRS) %>%
              mutate(PRS_CSx_st  = scale(PRS_CSx_st))

pheno_dt_unique <- readRDS("pheno_dt_unique.rds")
df_race <- fread("manifest.csv"); 
df_race <- df_race %>% rename(id = IID) %>% mutate(id = as.character(id))
df_race <- df_race %>% mutate(
            sex = case_when(
            sex2 == "male" ~ "M",
            sex2 == "female" ~ "F")) %>% select(-sex2)
df_bcancer <- df_race %>% filter(id %in% pheno_dt_unique$id) %>% select(id, sex)
  
# Check male breast cancer cases 
males <- df_bcancer %>% filter(sex == "M") %>% pull(id) 
pheno_dt_unique <- data.table(pheno_dt_unique)
br.cancer.cols <- grep("174", names(pheno_dt_unique), value = TRUE)
phenotypes <- pheno_dt_unique %>%
                mutate(across(all_of(br.cancer.cols), ~ if_else(id %in% males, NA, .))) %>%
                arrange(id) %>%
                as.data.frame()

# following code is for checking how many males are incldued in breast cancer cases. 
# pheno_dt_unique %>% select(all_of(br.cancer.cols)) %>% Hmisc::describe()
# pheno.dt.males <- pheno_dt_unique %>% filter(id %in% males)
# pheno.dt.females <- pheno_dt_unique %>% filter(!(id %in% males))
# pheno.dt.females %>% select(all_of(br.cancer.cols)) %>% skim() # # Check frequency table 
# pheno.dt.males %>% select(all_of(br.cancer.cols)) %>% Hmisc::describe()
# pheno.dt.males <- pheno.dt.males %>% mutate(across(all_of(br.cancer.cols), ~ NA))
# output <- rbind(pheno.dt.males, pheno.dt.females)
# setkey(output, "id") ; describe(output[,.SD,.SDcols=br.cancer.cols])
# phenotypes <- as.data.frame(output)

pheno_with_cov <- phenotypes %>% left_join(df_race, by = "id")
df_demo <- fread("person_GWAS.csv")
df_demo <- df_demo %>% select(SUBJID, YEAR_OF_BIRTH) %>% 
           mutate(id = as.character(SUBJID)) %>% select(-SUBJID)

df_demo <- df_demo %>% mutate(age = as.integer(format(Sys.Date(), "%Y")) - YEAR_OF_BIRTH) %>% 
            select(-YEAR_OF_BIRTH)
pheno_with_cov <- pheno_with_cov %>% left_join(df_demo, by = "id")

PRS_obj <- "PRS_CSx_st"
#PRS_obj <- "PRS_st_EUR"

df_PRS_temp <- df_PRS_all %>% select(IID, all_of(PRS_obj)) %>% 
               rename(id = IID)
pheno_with_cov <- pheno_with_cov %>% left_join(df_PRS_temp, by = "id")
phecode_sex_specific <- c()

phewas_fun <- function(target_phecode, covariates, genotypes, phecode_input, min_cases = 20) {
  df <- target_phecode %>%
        left_join(genotypes, by = "id") %>%
        left_join(covariates, by = "id")
  
  n_cases <- sum(df$phecode == 1, na.rm = TRUE)
  
  if (n_cases == 0) {
    message("n_cases are zero — skipping regression.")
    return(NULL)
  }
  
  sex_vec <- df %>% group_by(phecode, sex) %>% summarize(n=n(), .groups = "drop") %>% 
            filter(!is.na(phecode)) %>% pull(sex) %>% unique()
  # not sex-specific diseases
  if(length(sex_vec) == 2){
    cov_vec <- setdiff(names(df), c("id", "phecode"))
    form <- as.formula(
      paste("phecode", "~", paste(cov_vec, collapse = " + "))
    )
  }else{
    cov_vec <- setdiff(names(df), c("id", "phecode", "sex"))
    form <- as.formula(
      paste("phecode", "~", paste(cov_vec, collapse = " + "))
    )
    message(paste0("this phecode(", phecode_input, ") is a sex-specific code in eMERGE."))
  }
  phecode_sex_specific <- c(phecode_sex_specific, phecode_input)
  if (n_cases >= min_cases) {
    tryCatch({
      fit <- glm(form, data = df, family = binomial)
      summary_fit <- summary(fit)
      coef_mat <- coef(summary_fit)
      beta <- coef_mat[PRS_obj, "Estimate"]
      se   <- coef_mat[PRS_obj, "Std. Error"]
      pval <- coef_mat[PRS_obj, "Pr(>|z|)"]
      result_tab <- tibble(
        phecode = phecode,
        n_cases = n_cases,
        method = "glm",
        beta = beta,
        se = se,
        p = pval
      )
    }, error = function(e) {
      return(NULL)
    })
  } else {
    tryCatch({
      # Firth’s logistic regression for rare cases
      fit <- logistf(formula = form, data = df,
                     control = logistf.control(maxit = 100, maxstep = 0.5))
      result_tab <- tibble(
        phecode = phecode,
        n_cases = n_cases,
        method = "firth",
        beta = coef(fit)[PRS_obj],
        se   = sqrt(diag(vcov(fit)))[PRS_obj],
        p = fit$prob[PRS_obj]
      )
    }, error = function(e) {
      return(NULL)
    })
  }
  return(result_tab)
}

ex_names <- names(df_race)[-1]
ex_names[ex_names=="sex2"] <- "sex"
ex_names <- c(ex_names, "age", PRS_obj)
pheno_in <- pheno_with_cov %>% select(-all_of(ex_names))
cov_in <- pheno_with_cov %>% select(id, age, sex, paste0("PC",1:5))
geno_in  <- pheno_with_cov %>% select(id, all_of(PRS_obj))

phecodes_vec <- setdiff(names(pheno_in), "id")
tab_list <- list()
phecodes_vec <- setdiff(phecodes_vec, "175") # male breast cancer. we don't have male breast cancer in eMERGE

for(k in 1:length(phecodes_vec)){
  phecode <- phecodes_vec[k]
  phecode_col <- pheno_in[colnames(pheno_in) == phecode]
  temp_phe <- data.frame(id = pheno_in$id, phecode = phecode_col)
  names(temp_phe) <- c("id", "phecode")
  tab <- phewas_fun(target_phecode = temp_phe, covariates = cov_in, genotypes = geno_in, phecode_input = phecode)
  cat(paste0("phecode: ", phecode, ", index:", k, "\n"))
  if (is.null(tab)) next 
  tab_list[[k]] <- tab
}
  
combined <- bind_rows(tab_list)
#saveRDS(combined, file = "combined.rds")

setwd("/data")
#combined <- readRDS("combined_prs_eur.rds")
combined <- readRDS("combined.rds")

# load phecode map 
data("phecode_map")
data("pheinfo")
ref_phecode <- phecode_map %>% left_join(pheinfo, by = "phecode") 
results <- combined %>% inner_join(ref_phecode %>% distinct(phecode, .keep_all = TRUE), by = "phecode") # keep rows that show up in both dfs
results <- results %>% rename(phenotype = phecode)  
results <- results %>% mutate(sign = case_when(sign(beta) == 1 ~ "Pos", TRUE ~"Neg"))

results <- results %>%
            arrange(group, phenotype) %>%
            mutate(phenotype = factor(phenotype, levels = unique(phenotype)))

log10_bonf <- -log10(0.05 / nrow(results))

results <- results %>% mutate(phenotype_idx = as.numeric(phenotype))

# x position for x-label
group_positions <- results %>%
                    group_by(group) %>%
                    summarize(x_start = min(phenotype_idx),
                              x_end   = max(phenotype_idx),
                              x_mid   = mean(phenotype_idx))
x_labels <- rep("", length(levels(results$phenotype)))
x_labels[round(group_positions$x_mid)] <- group_positions$group

# PheWAS plot
p1 <- results %>% filter(!is.na(group)) %>% ggplot() +
      geom_point(aes(x = phenotype, y = -log10(p), color = group, shape = sign),
                 alpha = 0.6, size = 1.3) +
        annotate("segment",
               x = group_positions$x_end[-nrow(group_positions)] + 0.5,
               xend = group_positions$x_end[-nrow(group_positions)] + 0.5,
               y = -0, yend = 0.7,
               color = "black", linewidth = 0.3) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(-0, NA)) +
      scale_x_discrete(labels = x_labels) +
      geom_text_repel(
        data = results %>% filter(-log10(p) > log10_bonf),
        aes(x = phenotype, y = -log10(p), label = description),
        size = 2.7, max.overlaps = 100
      ) + 
      theme_minimal(base_size = 12) +
      geom_hline(yintercept = log10_bonf, linetype = "dashed", color = "navy", linewidth = 0.4) + 
      ylab(expression(-log[10](p))) + 
      theme(
        axis.text.x = element_text(angle = 40, hjust = 1),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray85", linewidth = 0.3),
        legend.position = "none"
      ) +
      scale_shape_manual(values = c("Pos" = 17, "Neg" = 25)) 
#ggsave("/fig/phewas2.pdf", p1, device = cairo_pdf, width = 5.7, height = 3)

### Volcano plot 
p2 <- results %>% filter(!is.na(group)) %>% ggplot() +
            geom_point(aes(x = beta, y = -log10(p), 
                           size = n_cases, 
                           color = (beta > 0)), 
                       alpha = 0.6) +
            geom_hline(yintercept = log10_bonf, 
                       linetype = "dashed", color = "gray20", linewidth = 0.4) +
            geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "darkgreen", linewidth = 0.4) +
            ylab(expression(-log[10](p))) + 
            xlab("log(OR) per 1 SD increase in PRS") +
            theme_minimal() +
            guides(size = "none") +
            scale_color_manual(
              name = NULL,
              values = c("TRUE" = "firebrick2", "FALSE" = "steelblue"),
              labels = c("Negative", "Positive")) +
            geom_text_repel(data = subset(results, -log10(p) > log10_bonf & abs(beta) > 0.2),
                            aes(x = beta, y = -log10(p), label = description), size = 3) + 
            theme(legend.title = element_blank(), legend.position = "bottom")

#ggsave("/fig/phe_vol2.pdf", p2, device = cairo_pdf, width = 6, height = 3)
