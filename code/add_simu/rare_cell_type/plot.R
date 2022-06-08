# Author: Zheng Li
# Date: 2022-04-04
# Plot results

library(tidyverse)
library(ggpattern)
library(gtable)
library(grid)
setwd("~/BASS_pjt/11_revision/2_rare_type/")
files <- list.files("results", full.names = T, pattern = "scenario5")

# cell type clustering
vars <- expand.grid(c("ari", "F1", "MCC"), "all", 
                    c("BASS", "sc3", "seu", "fict"))[, 3:1]
vars <- apply(vars, 1, paste, collapse = "_")
dat <- lapply(files, function(f){
  dat.i <- read.table(f, header = T, na.strings = "error")
  dat.i <- dat.i[, vars]
  dat.i$C <- strsplit(f, split = "_C|_R")[[1]][2]
  dat.i$dist <- strsplit(f, split = "_C|rareDist")[[1]][2]
  dat.i
})
dat <- do.call(rbind, dat)
dat <- gather(dat, case, value, BASS_all_ari:fict_all_MCC)
dat <- dat[complete.cases(dat), ]
dat$Metric <- sapply(strsplit(dat$case, "_"), `[`, 3)
dat$Method <- sapply(strsplit(dat$case, "_"), `[`, 1)
dat$Metric <- factor(dat$Metric)
levels(dat$Metric) <- c("ARI", "F1", "MCC")
dat$Method <- factor(dat$Method, levels = c("BASS", "sc3", "seu", "fict"))
levels(dat$Method) <- c("BASS", "SC3", "Seurat", "FICT")
dat$case <- factor(dat$case, levels = vars)

dist_labels <- c("Random", "Spatial")
names(dist_labels) <- c("random", "spatial")
C_labels <- c("10 (4 major + 6 rare types)", "14 (4 major + 10 rare types)")
names(C_labels) <- c("10", "14")
c_cols <- c("firebrick", "#377eb8", "#984ea3", "#e9c46a")

p <- ggplot(dat, aes(x = case, y = value, fill = Method, pattern = Metric)) +
  geom_boxplot_pattern(
    pattern_fill = 'black', pattern_colour = 'black', 
    pattern_density = .1, pattern_spacing = 0.02,
    lwd = 0.4, outlier.size = 1) +
  facet_grid(C ~ dist, labeller = labeller(C = C_labels, dist = dist_labels)) + 
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom") +
  ylim(c(0, 1)) +
  ylab("Accuracy") +
  xlab("") +
  scale_fill_manual(values = c_cols) +
  scale_pattern_manual(values = c('none', 'stripe', 'circle')) +
  guides(fill = guide_legend(override.aes = list(pattern = "none", order = 1)))

z <- ggplotGrob(p)
# add top title
z <- gtable_add_rows(z, z$height[7], pos = 6)
z <- gtable_add_grob(z, list(
  rectGrob(gp = gpar(col = "black", fill = "gray85", size = .5)),
  textGrob("Spatial distribution of rare cell types", gp = gpar(
    cex = 1, fontface = 'bold', col = "black"))), 
  t = 7, l = 5, b = 7, r = 7, name = c("a", "b"))
z <- gtable_add_rows(z, unit(1/10, "line"), 7)
# add right title
z <- gtable_add_cols(z, z$widths[8], pos = 11)
z <- gtable_add_grob(z, list(
  rectGrob(gp = gpar(col = "black", fill = "gray85", size = .5)),
  textGrob("Number of cell types (C)", rot = -90,  gp = gpar(
    cex = .75, fontface = 'bold', col = "black"))), 
  t = 10, l = 12, b = 12, r = 12, name = c("a", "b"))
z <- gtable_add_cols(z, unit(1/10, "line"), 11)
pdf("rare_cell_type.pdf", width = 6.7, height = 5)
grid.draw(z)
dev.off()


# spatial domain detection
proc_HMRF_out <- function(df)
{
  df <- df %>% select(contains("HMRF"))
  ari <- apply(df, 2, median, na.rm = T)
  worst_idx <- which.min(ari)
  med_idx <- which.min(abs(ari - median(ari)))
  best_idx <- which.max(ari)
  data.frame(HMRF_worst = df[, worst_idx], HMRF_med = df[, med_idx],
    HMRF_best = df[, best_idx])
}

spagcn_files <-  list.files("results/SpaGCN/", full.names = T)
dat <- lapply(files, function(f){
  dat.i <- read.table(f, header = T, na.strings = "error")
  hmrf.i <- proc_HMRF_out(dat.i)
  dat.i <- cbind(dat.i[, c("BASS_z", "BayesSpace", "seed")], hmrf.i)
  # spagnc
  sg.i <- read.table(grep(strsplit(f, "/")[[1]][2], spagcn_files, value = T), 
    header = T)
  dat.i <- merge(dat.i, sg.i, by = "seed")
  dat.i$C <- strsplit(f, split = "_C|_R")[[1]][2]
  dat.i$dist <- strsplit(f, split = "_C|rareDist")[[1]][2]
  dat.i
})

dat <- do.call(rbind, dat)
dat <- gather(dat, Method, ARI, BASS_z:SpaGCN)
dat <- dat[complete.cases(dat), ]
dat$Method <- factor(dat$Method, levels = c("BASS_z", "HMRF_worst", 
  "HMRF_med", "HMRF_best", "BayesSpace", "SpaGCN"))

z_cols <- c("firebrick", "#FA9F42", "#6A994E", "#deebf7", "#9ecae1", "#3182bd")
p <- ggplot(dat, aes(x = Method, y = ARI)) +
  geom_boxplot(aes(fill = Method), lwd = 0.3, outlier.size = 1) +
  geom_boxplot(aes(color = Method), fatten = NULL, fill = NA, coef = 0, 
    outlier.alpha = 0, show.legend = F, lwd = 0.5) +
  facet_grid(C ~ dist, labeller = labeller(C = C_labels, dist = dist_labels)) + 
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "bottom") +
  ylim(c(min(dat$ARI), 1)) +
  ylab("ARI") +
  xlab("") +
  scale_fill_manual(values = z_cols, 
    breaks=c("BASS_z", "BayesSpace", "SpaGCN", 
      "HMRF_worst", "HMRF_med", "HMRF_best"),
    labels = c("BASS", "BayesSpace", "SpaGCN", 
      "HMRF\n(worst)", "HMRF\n(median)", "HMRF\n(best)")) +
  scale_color_manual(values = z_cols[c(1, 4:6, 2:3)]) +
  guides(fill = guide_legend(byrow = T))

z <- ggplotGrob(p)
# add top title
z <- gtable_add_rows(z, z$height[7], pos = 6)
z <- gtable_add_grob(z, list(
  rectGrob(gp = gpar(col = "black", fill = "gray85", size = .5)),
  textGrob("Spatial distribution of rare cell types", gp = gpar(
    cex = 1, fontface = 'bold', col = "black"))), 
  t = 7, l = 5, b = 7, r = 7, name = c("a", "b"))
z <- gtable_add_rows(z, unit(1/10, "line"), 7)
# add right title
z <- gtable_add_cols(z, z$widths[8], pos = 11)
z <- gtable_add_grob(z, list(
  rectGrob(gp = gpar(col = "black", fill = "gray85", size = .5)),
  textGrob("Number of cell types (C)", rot = -90,  gp = gpar(
    cex = .75, fontface = 'bold', col = "black"))), 
  t = 10, l = 12, b = 12, r = 12, name = c("a", "b"))
z <- gtable_add_cols(z, unit(1/10, "line"), 11)
pdf("rare_cell_type_spatial_domian_detection.pdf", width = 6.7, height = 5)
grid.draw(z)
dev.off()

