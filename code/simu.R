# Author: Zheng Li
# Date: 2020-11-15
# Purpose:
# Simulation analysis for BASS propject with a 
# single slice as input
# Log:
# 2020-12-29: For unknown resons, when submitting jobs to other nodes, only the first replication
#             has a unique seed. All the other replications have the same seed used to generate
#             count data. Therefore, modified code to sample NREP seeds at the beginning and use 
#             those seeds to generate count data. Not a problem when running in local.
#             Parallel the simulation.
library(BASS)
library(pbmcapply)
wkdir <- "~/BASS_pjt/"
setwd(wkdir)
source("2_simu/simu_utils.R")

# data for inferring parameters for splatter
cnts <- readRDS("data/STARmap/20180417_BZ5_control/processed/cnts.RDS")
info <- readRDS("data/STARmap/20180417_BZ5_control/processed/info_mult.RDS")
xy <- as.matrix(info[, c("x", "y")])
starmap <- list(cnts = cnts, info = info)

# main
run_main <- function(
  scenario = args[1],
  C = as.numeric(args[2]),
  R = as.numeric(args[3]),
  J = as.numeric(args[4]),
  de_prop = as.numeric(args[5]),
  de_facLoc = as.numeric(args[6]),
  case = args[7],
  NREP = 50,
  outdir = "./",
  seed = 0,
  cores = 5
  )
{
  set.seed(seed)
  seeds <- sample(20201230, NREP) # 20201230
  filename <- paste0("scenario", scenario, "_C", C, "_R", R, "_J", J,
    "_prop", de_prop, "_de", de_facLoc)
  outfile <- file.path(outdir, filename)
  cat(filename, "\n")
  resolutions <- seq(0.3, 1.5, by = 0.1)
  header <- c("BASS_z", "BASS_c", "HMRF" %&% seq(0, by = 2, length = 26),
    "BayesSpace", "seu" %&% resolutions, "km", "fict", "mse_pi",
    "BASS_pi" %&% apply(expand.grid(1:C,1:R), 1, function(x) paste0(x[1], x[2])), 
    "seed")
  write.table(t(header), file = outfile, quote = F, col.names = F, row.names = F)

  verbose <- pbmclapply(1:NREP, function(rep){
      sim_dat <- simu(
        starmap = starmap,
        scenario = scenario, 
        C = C, 
        J = J, 
        L = 1,
        batch_facLoc = 0,
        de_prop = de_prop, 
        de_facLoc = de_facLoc, 
        de_facScale = 0.4,
        sim_seed = seeds[rep],
        debug = FALSE
        )
      BASS_out <- run_BASS(sim_dat, xy, "ACCUR_EST", beta = 1, C, R)
      HMRF_out <- run_HMRF(sim_dat, xy, info$z, R, case, rep,
        usePCs = F, dosearchSEgenes = T)
      BayesSpace_out <- run_BayesSpace(sim_dat, xy, info$z, R)
      seu_out <- seu_cluster(sim_dat, resolutions)
      km_out <- km_cluster(sim_dat, C)
      fict_out <- fict_cluster(sim_dat, xy, C, case, rep)

      res_out <- c(BASS_out$z_ari, BASS_out$c_ari, HMRF_out, 
        BayesSpace_out, seu_out, km_out, fict_out, BASS_out$mse_pi, 
        as.numeric(BASS_out$pi_est), sim_dat[[3]])
      write.table(t(res_out), file = outfile, quote = F, 
        col.names = F, row.names = F, append = T)
    }, mc.cores = cores)
}

args <- commandArgs(trailingOnly = T)
run_main(outdir = "2_simu/results/results_single/20211028", cores = 10)

