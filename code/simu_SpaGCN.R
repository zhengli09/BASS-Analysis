# Author: Zheng Li
# 2021-10-28
# Compare with SpaGCN

library(pbmcapply)
wkdir <- "~/BASS_pjt/1_code/BASS-analysis/"
setwd(wkdir)
source("code/simu_utils.R")
source_python("code/run_SpaGCN.py")

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
  case = args[7],
  NREP = 50,
  outdir = "./",
  seed = 0,
  cores = 5
  )
{
  set.seed(seed)
  seeds <- sample(20201230, NREP)
  filename <- paste0("scenario", scenario, "_C", C, "_R", R, "_J", J,
    "_prop", de_prop, "_de", de_facLoc)
  outfile <- file.path(outdir, filename)
  header <- c("SpaGCN", "seed")
  write.table(t(header), file = outfile, quote = F, 
    col.names = F, row.names = F)

  # SpaGCN has errors if running in Parallel in pbmclapply
  SpaGCN_res <- sapply(1:NREP, function(rep){
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
      SpaGCN_out <- run_SpaGCN(sim_dat, xy, info$z, R)$ari
    })
  res_out <- cbind(SpaGCN_res, seeds)
  write.table(res_out, file = outfile, quote = F, 
    col.names = F, row.names = F, append = T)
}

args <- commandArgs(trailingOnly = T)
run_main(outdir = "~/BASS_pjt/2_simu/results/results_single/20220421/SpaGCN", 
  cores = 1)

