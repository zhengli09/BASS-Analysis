# Author: Zheng Li
# Date: 2021-04-18
# Main simulation analyses for BASS project with multiple
# tissue sections data as input

library(pbmcapply)
wkdir <- "~/BASS_pjt/1_code/BASS-analysis/"
setwd(wkdir)
source("code/simu_utils.R")

# data for inferring parameters for splatter
cnts <- readRDS("data/simu_cnts.RDS")
info <- readRDS("data/simu_info.RDS")
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
  L = as.numeric(args[7]),
  batch_facLoc = as.numeric(args[8]),
  NREP = 50,
  outdir = "./",
  seed = 0,
  cores = 5
  )
{
  set.seed(seed)
  seeds <- sample(20201230, NREP)
  filename <- paste0("scenario", scenario, "_C", C, "_R", R, "_J", J,
    "_prop", de_prop, "_de", de_facLoc, "_L", L, "_batch", batch_facLoc)
  outfile <- file.path(outdir, filename)
  header <- c("BASS_z_ari", "BASS_c_ari", "BASS_beta", "BASS_pi_mse",
    "BASS_pi" %&% apply(expand.grid(1:C, 1:R), 1, paste0, collapse = ""),
    "seu_ari", "sc3_ari", "seed")
  write.table(t(header), file = outfile, quote = F, col.names = F, 
    row.names = F)

  verbose <- pbmclapply(1:NREP, function(rep){
      sim_dat <- simu(
        starmap = starmap,
        scenario = scenario, 
        C = C, 
        J = J, 
        L = L,
        batch_facLoc = batch_facLoc,
        de_prop = de_prop, 
        de_facLoc = de_facLoc, 
        de_facScale = 0.4,
        sim_seed = seeds[rep],
        debug = FALSE
        )
      BASS_out <- run_BASS(sim_dat, xy, "SW", beta = 1, C, R,
        init_method = "kmeans")
      seu_out <- seu_cluster(sim_dat, C)
      sc3_out <- sc3_cluster(sim_dat, C)

      res_out <- c(BASS_out$z_ari, BASS_out$c_ari[1], BASS_out$beta,
        BASS_out$mse_pi, as.numeric(BASS_out$pi_est), 
        seu_out$metric[1], sc3_out$metric[1], sim_dat[[3]])
      write.table(t(res_out), file = outfile, quote = F, 
        col.names = F, row.names = F, append = T)
    }, mc.cores = cores)
}

args <- commandArgs(trailingOnly = T)
verbose <- run_main(outdir = "~/BASS_pjt/2_simu/results/results_mult/20220421", 
  cores = 10)
print(verbose) # debug purpose

