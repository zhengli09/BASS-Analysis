# Author: Zheng Li
# Date: 2021-04-06
# Purpose:
# Set the number of cell types to be
# 2, 4, 6, 8, 10 while the number
# of spatial domains is set to be 4

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
  case = args[7],
  NREP = 50,
  outdir = "./",
  seed = 0,
  cores = 5
)
{
  set.seed(seed)
  seeds <- sample(20201230, NREP)
  filename <- paste0("scenario", scenario, "_C", C, "_R", R, 
    "_J", J, "_prop", de_prop, "_de", de_facLoc)
  outfile <- file.path(outdir, filename)
  header <- expand.grid(c("all_ari", "C_est", "c_clust_prop"), 
    c("seu", "sc3", "fict"))[2:1]
  header <- apply(header, 1, paste, collapse = "_")
  header <- c("BASS_z", "BASS_all_ari", "BASS_C_est", "BASS_R_est",
    "BASS_c_clust_prop", "BASS_z_clust_prop", header, "seed")
  write.table(t(header), file = outfile, quote = F, col.names = F,
    row.names = F)
  
  verbose <- pbmclapply(1:NREP, function(rep){
    sim_dat <- simu(
      starmap = starmap,
      scenario = scenario, 
      C = 4, 
      J = J, 
      L = 1,
      batch_facLoc = 0,
      de_prop = de_prop, 
      de_facLoc = de_facLoc,
      de_facScale = 0.4,
      sim_seed = seeds[rep],
      debug = FALSE
    )
    
    # run all the methods
    BASS_out <- run_BASS(sim_dat, xy, "SW", beta = 1, 
      C, R, init_method = "kmeans")
    seu_out <- seu_cluster(sim_dat, C)
    sc3_out <- sc3_cluster(sim_dat, C)
    fict_out <- fict_cluster(sim_dat, xy, C, case, rep)

    # output
    res_out <- c(BASS_out$z_ari, BASS_out$c_ari[1], 
      unlist(BASS_out[c("C_est", "R_est", "c_clust_prop", "z_clust_prop")]),
      seu_out$metric[1], seu_out$C_est, seu_out$clust_props,
      sc3_out$metric[1], sc3_out$C_est, sc3_out$clust_props,
      fict_out$metric[1], fict_out$C_est, fict_out$clust_props,
      sim_dat[[3]])
    write.table(t(res_out), file = outfile, quote = F, col.names = F, 
      row.names = F, append = T)
  }, mc.cores = cores)
}

args <- commandArgs(trailingOnly = T)
verbose <- run_main(outdir = "~/BASS_pjt/11_revision/3_hyperParam_CnR/results/C", 
         cores = 10)
print(verbose) # debugging
