# Author: Zheng Li
# Date: 2022-04-09
# Plot results of sensitivity analysis on C and R

library(tidyverse)
library(cowplot)
# impact of C on C and C on R
setwd("~/BASS_pjt/11_revision/3_hyperParam_CnR/")
files <- list.files("results/C", full.names = T)
dat <- lapply(files, function(f){
  dat.i <- read.table(f, header = T, na.string = "error")
  C_truth <- strsplit(f, split = "_C|_R")[[1]][2]
  dat.i$C <- C_truth
  dat.i
})
dat <- do.call(rbind, dat)
## impact on ARI
df <- dat %>% 
  select("C", "BASS_all_ari", "sc3_all_ari", "seu_all_ari", "fict_all_ari") %>%
  gather(Method, ARI, BASS_all_ari:fict_all_ari) %>%
  filter(complete.cases(.)) %>%
  mutate(C = factor(C, levels = c(2, 4, 6, 8, 10))) %>%
  mutate(Method = recode(Method,
    "BASS_all_ari" = "BASS",
    "sc3_all_ari" = "SC3",
    "seu_all_ari" = "Seurat",
    "fict_all_ari" = "FICT")) %>%
  mutate(Method = factor(Method, levels = c("BASS", "SC3", "Seurat", "FICT")))

p1 <- ggplot(df, aes(x = C, y = ARI)) +
  geom_boxplot(fill = "#6c757d", outlier.size = 1) +
  geom_boxplot(color = "#6c757d", fatten = NULL, fill = NA, coef = 0, 
    outlier.alpha = 0, show.legend = F) +
  facet_grid(~ Method) +
  theme_bw() +
  ylim(c(0, 1)) +
  ylab("ARI of cell\ntype clustering") +
  theme(axis.title.x = element_blank())

## impact on estimated number of cell types
df <- dat %>% 
  select("C", "BASS_C_est", "sc3_C_est", "seu_C_est", "fict_C_est") %>%
  gather(Method, C_est, BASS_C_est:fict_C_est) %>%
  filter(complete.cases(.)) %>%
  mutate(C = factor(C, levels = c(2, 4, 6, 8, 10))) %>%
  mutate(Method = recode(Method,
    "BASS_C_est" = "BASS",
    "sc3_C_est" = "SC3",
    "seu_C_est" = "Seurat",
    "fict_C_est" = "FICT")) %>%
  mutate(Method = factor(Method, levels = c("BASS", "SC3", "Seurat", "FICT")))

p2 <- ggplot(df, aes(x = C, y = C_est)) +
  geom_boxplot(fill = "#6c757d", outlier.size = 1) +
  geom_boxplot(color = "#6c757d", fatten = NULL, fill = NA, coef = 0, 
    outlier.alpha = 0, show.legend = F) +
  facet_grid(~ Method) +
  theme_bw() +
  ylab("Estimated number\nof cell types") +
  theme(axis.title.x = element_blank()) +
  scale_y_continuous(breaks = seq(0, 10, 2), limits = c(0, 10))

## impact on the proportion of each cell type cluster
df <- dat %>% 
  select("C", "BASS_c_clust_prop", "sc3_c_clust_prop", 
    "seu_c_clust_prop", "fict_c_clust_prop")
colnames(df) <- c("C", "BASS", "SC3", "Seurat", "FICT")
df <- lapply(c("BASS", "SC3", "Seurat", "FICT"), function(method){
  df.i <- lapply(seq(2, 10, 2), function(c){
    df.ii <- df %>% select(C, method) %>% filter(C == c)
    props <- strsplit(df.ii[, 2], ",")
    max_nClust <- max(sapply(props, length))
    props <- lapply(props, function(props.i){
      as.numeric(c(props.i, rep(0, max_nClust - length(props.i))))
    })
    props <- apply(do.call(rbind, props), 2, mean, na.rm = T)
    data.frame(Method = method, C = c, 
      Cluster = paste("Cluster", 1:max_nClust),
      Proportion = props)
  })
  do.call(rbind, df.i)
})
df <- do.call(rbind, df)
df$C <- factor(df$C, levels = seq(2, 10, 2))
df$Method <- factor(df$Method, levels = c("BASS", "SC3", "Seurat", "FICT"))
df$Cluster <- factor(df$Cluster, levels = paste("Cluster", 1:10))

p3 <- ggplot(df, aes(x = C, y = Proportion, fill = Cluster)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(~ Method) +
  theme_bw() +
  xlab("Specified number of cell types (Truth = 4)") +
  ylab("Proportion of total\ncell population") +
  scale_fill_manual(values = c("#001219","#005f73","#0a9396","#94d2bd","#e9d8a6",
    "#ee9b00","#ca6702","#bb3e03","#ae2012","#9b2226")) +
  theme(legend.position = "bottom",
    legend.key.size = unit(0.2, 'cm'))

## impact on R
df <- dat %>% select("C", "BASS_z", "BASS_R_est", "BASS_z_clust_prop")
df$C <- factor(df$C, levels = seq(0, 10, 2))
p4 <- ggplot(df, aes(x = C, y = BASS_z)) +
  geom_boxplot(fill = "#6c757d", outlier.size = 1) +
  geom_boxplot(color = "#6c757d", fatten = NULL, fill = NA, coef = 0, 
    outlier.alpha = 0, show.legend = F) +
  theme_bw() +
  ylim(c(0, 1)) +
  ylab("ARI of spatial\ndomain detection") +
  theme(axis.title.x = element_blank())

p5 <- ggplot(df, aes(x = C, y = BASS_R_est)) +
  geom_boxplot(fill = "#6c757d", outlier.size = 1) +
  geom_boxplot(color = "#6c757d", fatten = NULL, fill = NA, coef = 0, 
    outlier.alpha = 0, show.legend = F) +
  theme_bw() +
  ylab("Estimated number\nof spatial domains") +
  scale_y_continuous(breaks = 2:5, limits = c(2, 5)) +
  theme(axis.title.x = element_blank())

df <- dat %>% select("C", "BASS_z", "BASS_R_est", "BASS_z_clust_prop")
df <- lapply(seq(2, 10, 2), function(c){
  df.i <- df %>% filter(C == c)
  props <- strsplit(df.i$BASS_z_clust_prop, ",")
  max_nClust <- max(sapply(props, length))
  props <- lapply(props, function(props.i){
    as.numeric(c(props.i, rep(0, max_nClust - length(props.i))))
  })
  props <- apply(do.call(rbind, props), 2, mean, na.rm = T)
  data.frame(C = c, Domain = paste("Domain", 1:max_nClust), 
    Proportion = props)
})
df <- do.call(rbind ,df)
df$C <- factor(df$C, levels = seq(2, 10, 2))

p6 <- ggplot(df, aes(x = C, y = Proportion, fill = Domain)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  ylab("Proportion of total\ncell population") +
  scale_fill_manual(values = c("#001219","#005f73","#0a9396","#94d2bd")) +
  theme(axis.title.x = element_blank(),
    legend.key.size = unit(0.3, "cm"))

xlab <- ggdraw() + 
  draw_label("Specified number of cell types (Truth = 4)", size = 12) +
  theme(plot.margin = margin(10, 0, 20, 0))

p123 <- plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(0.6, 0.6, 1))
p456 <- plot_grid(p4, p5, p6, ncol = 3, rel_widths = c(0.6, 0.6, 1))
pdf("impact_C.pdf", width = 6.7, height = 7)
plot_grid(p123, p456, xlab, ncol = 1, rel_heights = c(3, 0.8, 0.1))
dev.off()

#############################
# impact of R on C and R on R
proc_HMRF_out <- function(df)
{
  ari <- df %>% 
    select(contains("HMRF")) %>% 
    select(contains("ari")) %>%
    apply(2, median, na.rm = T)
  best_beta <- strsplit(names(ari)[which.max(ari)], "_")[[1]][1]
  df %>% select(contains(best_beta))
}

files <- list.files("results/R", full.names = T, pattern = "scenario")
spagcn_files <-  list.files("results/R/SpaGCN", full.names = T)
dat <- lapply(files, function(f){
  dat.i <- read.table(f, header = T, na.string = "error")
  hmrf.i <- proc_HMRF_out(dat.i)
  names(hmrf.i) <- c("HMRF_ari", "HMRF_R_est", "HMRF_z_clust_prop")
  dat.i <- dat.i %>% select(contains(c("BASS", "BayesSpace", "seed")))
  dat.i <- cbind(dat.i, hmrf.i)
  sg.i <- read.table(grep(strsplit(f, "/")[[1]][3], spagcn_files, 
    value = T), header = T)
  dat.i <- merge(dat.i, sg.i, by = "seed")
  R_truth <- strsplit(f, split = "_R|_J")[[1]][2]
  dat.i$R <- R_truth
  dat.i
})
dat <- do.call(rbind, dat)

## impact on ARI
df <- dat %>% 
  select("R", "BASS_z", "HMRF_ari", "BayesSpace_ari", "SpaGCN_ari") %>%
  gather(Method, ARI, BASS_z:SpaGCN_ari) %>%
  filter(complete.cases(.)) %>%
  mutate(R = factor(R, levels = c(2, 4, 6, 8, 10))) %>%
  mutate(Method = recode(Method,
    "BASS_z" = "BASS",
    "HMRF_ari" = "HMRF",
    "BayesSpace_ari" = "BayesSpace",
    "SpaGCN_ari" = "SpaGCN")) %>%
  mutate(Method = factor(Method, levels = c("BASS", "HMRF", "BayesSpace", "SpaGCN")))

p1 <- ggplot(df, aes(x = R, y = ARI)) +
  geom_boxplot(fill = "#6c757d", outlier.size = 1) +
  geom_boxplot(color = "#6c757d", fatten = NULL, fill = NA, coef = 0, 
    outlier.alpha = 0, show.legend = F) +
  facet_grid(~ Method) +
  theme_bw() +
  ylim(c(0, 1)) +
  ylab("ARI of spatial\ndomain detection") +
  theme(axis.title.x = element_blank())

## impact on estimated number of cell types
df <- dat %>% 
  select("R", "BASS_R_est", "HMRF_R_est", "BayesSpace_R_est", "SpaGCN_R_est") %>%
  gather(Method, R_est, BASS_R_est:SpaGCN_R_est) %>%
  filter(complete.cases(.)) %>%
  mutate(R = factor(R, levels = c(2, 4, 6, 8, 10))) %>%
  mutate(Method = recode(Method,
    "BASS_R_est" = "BASS",
    "HMRF_R_est" = "HMRF",
    "BayesSpace_R_est" = "BayesSpace",
    "SpaGCN_R_est" = "SpaGCN")) %>%
  mutate(Method = factor(Method, levels = c("BASS", "HMRF", "BayesSpace", "SpaGCN")))

p2 <- ggplot(df, aes(x = R, y = R_est)) +
  geom_boxplot(fill = "#6c757d", outlier.size = 1) +
  geom_boxplot(color = "#6c757d", fatten = NULL, fill = NA, coef = 0, 
    outlier.alpha = 0, show.legend = F) +
  facet_grid(~ Method) +
  theme_bw() +
  ylab("Estimated number\nof spatial domains") +
  theme(axis.title.x = element_blank()) +
  scale_y_continuous(breaks = seq(0, 12, 2), limits = c(0, 12))

## impact on the proportion of each cell type cluster
df <- dat %>% 
  select("R", "BASS_z_clust_prop", "HMRF_z_clust_prop", 
    "BayesSpace_z_clust_prop", "SpaGCN_z_clust_prop")
colnames(df) <- c("R", "BASS", "HMRF", "BayesSpace", "SpaGCN")
df <- lapply(c("BASS", "HMRF", "BayesSpace", "SpaGCN"), function(method){
  df.i <- lapply(seq(2, 10, 2), function(r){
    df.ii <- df %>% select(R, method) %>% filter(R == r)
    props <- strsplit(df.ii[, 2], ",")
    max_nClust <- max(sapply(props, length))
    props <- lapply(props, function(props.i){
      as.numeric(c(props.i, rep(0, max_nClust - length(props.i))))
    })
    props <- apply(do.call(rbind, props), 2, mean, na.rm = T)
    data.frame(Method = method, R = r, 
      Domain = paste("Domain", 1:max_nClust),
      Proportion = props)
  })
  do.call(rbind, df.i)
})
df <- do.call(rbind, df)
df$R <- factor(df$R, levels = seq(2, 10, 2))
df$Method <- factor(df$Method, levels = c("BASS", "HMRF", "BayesSpace", "SpaGCN"))
df$Domain <- factor(df$Domain, levels = paste("Domain", 1:11))

p3 <- ggplot(df, aes(x = R, y = Proportion, fill = Domain)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(~ Method) +
  theme_bw() +
  xlab("Specified number of spatial domains (Truth = 4)") +
  ylab("Proportion of total\ncell population") +
  scale_fill_manual(values = c("#001219","#005f73","#0a9396","#94d2bd","#e9d8a6",
    "#ee9b00","#ca6702","#bb3e03","#ae2012","#9b2226", "#6a040f")) +
  theme(legend.position = "bottom",
    legend.key.size = unit(0.2, 'cm')) +
  guides(fill = guide_legend(nrow = 2, byrow = T))

## impact on of R on C
df <- dat %>% select("R", "BASS_all_ari", "BASS_C_est", "BASS_c_clust_prop")
df$R <- factor(df$R, levels = seq(0, 10, 2))
p4 <- ggplot(df, aes(x = R, y = BASS_all_ari)) +
  geom_boxplot(fill = "#6c757d", outlier.size = 1) +
  geom_boxplot(color = "#6c757d", fatten = NULL, fill = NA, coef = 0, 
    outlier.alpha = 0, show.legend = F) +
  theme_bw() +
  ylim(c(.75, 1)) +
  ylab("ARI of cell\ntype clustering") +
  theme(axis.title.x = element_blank())

p5 <- ggplot(df, aes(x = R, y = BASS_C_est)) +
  geom_boxplot(fill = "#6c757d", outlier.size = 1) +
  geom_boxplot(color = "#6c757d", fatten = NULL, fill = NA, coef = 0, 
    outlier.alpha = 0, show.legend = F) +
  theme_bw() +
  ylab("Estimated number\nof cell types") +
  scale_y_continuous(breaks = 3:5, limits = c(3, 5)) +
  theme(axis.title.x = element_blank())

df <- dat %>% select("R", "BASS_all_ari", "BASS_C_est", "BASS_c_clust_prop")
df <- lapply(seq(2, 10, 2), function(r){
  df.i <- df %>% filter(R == r)
  props <- strsplit(df.i$BASS_c_clust_prop, ",")
  max_nClust <- max(sapply(props, length))
  props <- lapply(props, function(props.i){
    as.numeric(c(props.i, rep(0, max_nClust - length(props.i))))
  })
  props <- apply(do.call(rbind, props), 2, mean, na.rm = T)
  data.frame(R = r, Cluster = paste("Cluster", 1:max_nClust), 
    Proportion = props)
})
df <- do.call(rbind ,df)
df$R <- factor(df$R, levels = seq(2, 10, 2))

p6 <- ggplot(df, aes(x = R, y = Proportion, fill = Cluster)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  ylab("Proportion of total\ncell population") +
  scale_fill_manual(values = c("#001219","#005f73","#0a9396","#94d2bd")) +
  theme(axis.title.x = element_blank(),
    legend.key.size = unit(0.3, "cm"))

xlab <- ggdraw() + 
  draw_label("Specified number of spatial domains (Truth = 4)", size = 12) +
  theme(plot.margin = margin(10, 0, 20, 0))

p123 <- plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(0.6, 0.6, 1))
p456 <- plot_grid(p4, p5, p6, ncol = 3, rel_widths = c(0.6, 0.6, 1))
pdf("impact_R.pdf", width = 6.7, height = 7)
plot_grid(p123, p456, xlab, ncol = 1, rel_heights = c(3, 0.8, 0.1))
dev.off()


