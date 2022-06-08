# Author: Zheng Li
# Date: 2022-04-11
# Purpose:
# Plot results from analysis where genes were randomly
# excluded from the spatial transcriptomics data set

setwd("~/BASS_pjt/11_revision/4_exclude_genes/")
library(tidyverse)
library(cowplot)

main_files <- list.files("results", pattern = "scenario", full.names = T)
spagcn_files <- list.files("results/SpaGCN", full.names = T)
proc_HMRF_out <- function(df)
{
  df <- df %>% select(contains("HMRF"))
  ari <- apply(df, 2, median, na.rm = T)
  best_idx <- which.max(ari)
  data.frame(HMRF = df[, best_idx])
}

dat <- lapply(main_files, function(f){
  dat.i <- read.table(f, header = T, na.strings = "error")
  hmrf.i <- proc_HMRF_out(dat.i)
  dat.i <- dat.i %>% select(!contains("HMRF"))
  dat.i <- cbind(dat.i, hmrf.i)
  spagcn.i <- read.table(grep(strsplit(f, "/")[[1]][2], spagcn_files, 
    value = T), header = T, na.strings = "error")
  dat.i <- merge(dat.i, spagcn.i, by = "seed")
  dat.i$scenario <- strsplit(f, "scenario|_C")[[1]][2]
  dat.i$I <- strsplit(f, "_I|_prop")[[1]][2]
  dat.i
})
dat <- do.call(rbind, dat)


# impact on cell type clustering
df <- dat %>% 
  select(c("scenario", "I", "BASS_c_ari", "seu_ari", "sc3_ari", "fict_ari")) %>%
  gather(Method, ARI, BASS_c_ari:fict_ari) %>%
  filter(complete.cases(.)) %>%
  mutate(Method = recode(Method,
    "BASS_c_ari" = "BASS",
    "fict_ari" = "FICT",
    "sc3_ari" = "SC3",
    "seu_ari" = "Seurat")) %>%
  mutate(Method = factor(Method, levels = c("BASS", "SC3", "Seurat", "FICT")))
scenarios <- c(paste("Scenario", c("I", "II", "III", "IV")))
names(scenarios) <- 1:4

c_cols <- c("firebrick", "#377eb8", "#984ea3", "#e9c46a")
p1 <- ggplot(df, aes(x = I, y = ARI, fill = Method)) +
  facet_grid(Method ~ scenario, labeller = labeller(scenario = scenarios)) +
  geom_boxplot(lwd = 0.3, outlier.size = 1) +
  geom_boxplot(aes(color = Method), fatten = NULL, fill = NA, coef = 0, 
    outlier.alpha = 0, show.legend = F, lwd = 0.5) +
  xlab("Number of subsampled genes (original = 1,000)") +
  ylab("ARI of cell type clustering") +
  theme_bw() + 
  theme(legend.position = "none",
    strip.text = element_text(face = "bold")) +
  scale_fill_manual(values = c_cols) +
  scale_color_manual(values = c_cols)


# impact on spatial domain detection
df <- dat %>% 
  select(c("scenario", "I", "BASS_z_ari", "HMRF", "BayesSpace", "SpaGCN")) %>%
  gather(Method, ARI, BASS_z_ari:SpaGCN) %>%
  filter(complete.cases(.)) %>%
  mutate(Method = recode(Method, "BASS_z_ari" = "BASS")) %>%
  mutate(Method = factor(Method, levels = c("BASS", "HMRF", "BayesSpace", "SpaGCN")))
scenarios <- c(paste("Scenario", c("I", "II", "III", "IV")))
names(scenarios) <- 1:4

z_cols <- c("firebrick", "#3182bd", "#FA9F42", "#6A994E")
p2 <- ggplot(df, aes(x = I, y = ARI, fill = Method)) +
  facet_grid(Method ~ scenario, labeller = labeller(scenario = scenarios)) +
  geom_boxplot(lwd = 0.3, outlier.size = 1) +
  geom_boxplot(aes(color = Method), fatten = NULL, fill = NA, coef = 0, 
    outlier.alpha = 0, show.legend = F, lwd = 0.5) +
  xlab("Number of subsampled genes (original = 1,000)") +
  ylab("ARI of spatial domain detection") +
  theme_bw() + 
  theme(legend.position = "none",
    strip.text = element_text(face = "bold")) +
  scale_fill_manual(values = z_cols) +
  scale_color_manual(values = z_cols)

pdf("exclude_genes.pdf", width = 6.7, height = 8)
plot_grid(p1, p2, ncol = 1, labels = c("A", "B"))
dev.off()

