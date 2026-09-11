# module load R/4.4.1
# options(repos = c(CRAN = "https://cran.r-project.org"))
# install.packages(c("remotes", "withr", "desc"))

lib_path <- "/Rlib"
.libPaths(c(lib_path, .libPaths()))
# Package installation check
my.packages <- c("data.table", "Hmisc", "tidyverse", "gdsfmt", "SNPRelate",
                 "bigparallelr", "ggplot2","BiocManager", "bigstatsr", "MASS",
                 "bigsnpr")
my.packages2 <- c("GENESIS", "SeqArray", "SeqVarTools", "GWASTools", "GENESIS")

lapply(c(my.packages, my.packages2), library, character.only = TRUE)
data_path <- "/final_merged_pca/"
num_threads <- parallel::detectCores() # For parallel computing in PCA

# 1KG data load
bed.fn <- paste0(data_path, "megred_pruned_all.bed")
fam.fn <- paste0(data_path, "megred_pruned_all.fam")
bim.fn <- paste0(data_path, "megred_pruned_all.bim")
gds.fn <- paste0(data_path, "megred_pruned_all.gds")

# Converts PLINK binary format -> GDS.
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, gds.fn)

# Summary
genofile <- snpgdsOpen(gds.fn)

# To check gds file
# head(read.gdsn(index.gdsn(genofile, "sample.id")))
# read.gdsn(index.gdsn(genofile, "snp.id"))[1:5]
# read.gdsn(index.gdsn(genofile, "genotype"))[1:5, 1:5]
# exclude kinship > 0.1

kin_df <- read.table(paste0(data_path, "kinship_filtered.txt"), header = TRUE, sep = "\t")

sample_ids_ori <- read.gdsn(index.gdsn(genofile, "sample.id"))
index <- sample_ids_ori %in% unique(kin_df$ID1)
sample_ids <- sample_ids_ori[!index]

KG_samples <- sample_ids[grepl("^(NA|HG)", sample_ids)]
emerge_id <- sample_ids[!(sample_ids %in% KG_samples)]

# (same scale of standardization of 1KG should be applied to validation data)

# PCA with 1KG
pca_ref <- snpgdsPCA(genofile, sample.id = KG_samples, num.thread = num_threads, eigen.cnt = 30)
#str(pca_ref); pc_percent <- pca_ref$varprop * 100; head(round(pc_percent, 2))

# loadings/eigenvectors from reference PCA (1KG)
snp_load <- snpgdsPCASNPLoading(pca_ref, gdsobj = genofile)

# Project eMERGE into the reference 1KG space
samp_load <- snpgdsPCASampLoading(snp_load, sample.id = emerge_id, gdsobj = genofile, num.thread = num_threads)

#all.equal(samp_load$sample.id, emerge_id)

#saveRDS(pca_ref, file = paste0(data_path, "pca_ref.rds"))
#saveRDS(snp_load, file = paste0(data_path, "snp_load.rds"))
#saveRDS(samp_load, file = paste0(data_path, "samp_load.rds"))



####################################################################
# PCA figure

library(MASS)
library(dplyr)
library(ggplot2)
library(ggsci)
library(scales)
library(ggplot2)
library(ggpubr)

pca_ref <- readRDS("/PCA/pca_ref.rds")   #reference 1KG PCA
snp_load <- readRDS("/PCA/snp_load.rds")   #loading from 1KG
samp_load <- readRDS("/PCA/samp_load.rds")   #projected eMERGE on 1KG

# ancestry from 1KG
KG_panel <- read.table("ALL.panel", header = TRUE, stringsAsFactors = FALSE)

panel_sub <- KG_panel %>% filter(sample %in% pca_ref$sample.id) %>%
              rename(sample.id = sample)

## 1KG
tab <- data.frame(sample.id = pca_ref$sample.id,
                  EV1 = pca_ref$eigenvect[,1],
                  EV2 = pca_ref$eigenvect[,2])

tab <- tab %>% left_join(panel_sub, by = "sample.id")

plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
pc_percent <- (pca_ref$varprop * 100)[1:20]
pc_percent

df_var <- data.frame(PC = 1:length(pc_percent), prop = pc_percent) %>% mutate(cumsum = cumsum(prop))
p1 <- df_var %>% ggplot(aes(x = PC, y = prop, linetype = 1)) + geom_line() + geom_point() + geom_bar()
p1 %>% ggplot() + geom_line(aes(x=PC, y=cumsum, linetype = 2 ))

# pal <- pal_npg()(10)
# show_col(pal)

pal <- pal_lancet()(9)
show_col(pal)

p1 <- ggplot(df_var, aes(x = PC)) +
      geom_col(aes(y = prop), fill = pal[4], alpha = 0.8) +
      geom_line(aes(y = prop, color = "Proportion", linetype = "Proportion"), linewidth = 0.5) +
      geom_point(aes(y = prop, color = "Proportion"), size = 1.2, alpha=0.7) +
      geom_line(aes(y = cumsum, color = "Cumulative", linetype = "Cumulative"), linewidth = 0.5) +
      geom_point(aes(y = cumsum, color = "Cumulative"), size = 1.2, alpha=0.7) +
      scale_color_manual(values = c("Proportion" = pal[1], "Cumulative" = pal[7])) +
      scale_linetype_manual(values = c("Proportion" = 1, "Cumulative" = 5)) +
      labs(x = "Principal Component", y = "Variance Explained (%)", color = "",
        linetype = "") + theme_bw() +
      theme(legend.position = "bottom")

#ggsave("~/pc_var_fig.pdf", plot = p1, width = 7, height = 3.8, device = cairo_pdf)

lbls <- paste("PC", 1:4, "\n", format(pc_percent[1:4], digits=2), "%", sep="")
pop <- factor(panel_sub$super_pop)
pop_colors <- rainbow(length(levels(pop)))[pop]
pairs(pca_ref$eigenvect[,1:4], col = pop, labels=lbls)

## eMERGE group assignment

pc_all <- rbind(pca_ref$eigenvect[, 1:30], samp_load$eigenvect[, 1:30])
group_label <- c(panel_sub$super_pop, rep(NA, length(samp_load$sample.id)))

library(class)
emerge_pop <- knn(train = pca_ref$eigenvect[, 1:2],
                    test = samp_load$eigenvect[, 1:2],
                    cl = panel_sub$super_pop,
                    k = 5)

KG_df <- data.frame(score = pca_ref$eigenvect[, 1:2],
                    race = panel_sub$super_pop,
                    group = "1KG")

emerge_df <- data.frame(score = samp_load$eigenvect[, 1:2],
                    race = emerge_pop,
                    group = "eMERGE")

df <- rbind(KG_df, emerge_df) %>%
      mutate(alpha_val = ifelse(group == "1KG", 0.5, 1))

p11 <- df %>% ggplot(aes(x = score.1, y = score.2, col = race)) +
      geom_point(alpha = 0.4, size = 0.7) + facet_grid(~ group) +
      scale_color_brewer(palette = "Dark2") +
      theme_bw() +
      theme(legend.title = element_blank()) + xlab("PC1") + ylab("PC2")


p1 <- df %>% filter(group == "eMERGE") %>% ggplot(aes(x = score.1, y = score.2)) +
  geom_point(alpha = 0.4, size = 0.7, color = "gray") +
  theme_bw() +
  theme(legend.title = element_blank()) + xlab("PC1") + ylab("PC2")

p2 <- p1 + geom_point(data = df %>% filter(group != "eMERGE"),
                aes(x = score.1, y = score.2, col = race), alpha = 0.4, size = 0.7 ) +
      scale_color_brewer(palette = "Dark2") +
      theme_bw() + theme(legend.title = element_blank()) +
      xlab("PC1") + ylab("PC2")

plot <- ggarrange(p11, p2, ncol = 2, widths = c(1.7, 1), nrow = 1,
          common.legend = TRUE, legend = "bottom")
