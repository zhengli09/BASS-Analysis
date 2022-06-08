# Author: Zheng Li
# Date: 2021-04-03
# Purpose:
# New simulation scenario (scenario 5) that have simulated rare cell types

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
  rare_dist = args[2],
  C = as.numeric(args[3]),
  R = as.numeric(args[4]),
  J = as.numeric(args[5]),
  de_prop = as.numeric(args[6]),
  de_facLoc = as.numeric(args[7]),
  case = args[8],
  NREP = 50,
  outdir = "./",
  seed = 0,
  cores = 5
)
{
  set.seed(seed)
  seeds <- sample(20201230, NREP)
  filename <- paste0("scenario", scenario, "_rareDist", rare_dist,  
                     "_C", C, "_R", R, "_J", J, "_prop", de_prop, 
                     "_de", de_facLoc)
  outfile <- file.path(outdir, filename)
  header <- expand.grid(c("all", "major", "rare"), c("ari", "F1", "MCC"),
    c("BASS", "seu", "sc3", "fict"))[, c(3, 1, 2)]
  header <- apply(header, 1, paste, collapse = "_")
  header <- c("BASS_z", header, "HMRF" %&% seq(0, by = 2, length = 26), 
              "BayesSpace", "seed")
  write.table(t(header), file = outfile, quote = F, col.names = F,
    row.names = F)
  
  verbose <- pbmclapply(1:NREP, function(rep){
    sim_dat <- simu(
      starmap = starmap,
      scenario = scenario, 
      rare_dist = rare_dist,
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
    
    # run all the methods
    BASS_out <- run_BASS(sim_dat, xy, "SW", beta = 1, C, R,
      init_method = "kmeans")
    seu_out <- seu_cluster(sim_dat, C)
    sc3_out <- sc3_cluster(sim_dat, C)
    fict_out <- fict_cluster(sim_dat, xy, C, case, rep)
    HMRF_out <- run_HMRF(sim_dat, xy, info$z, R, case, rep, 
                         usePCs = F, dosearchSEgenes = T)
    BayesSpace_out <- run_BayesSpace(sim_dat, xy, info$z, R)
    
    # output
    res_out <- c(BASS_out$z_ari, unlist(BASS_out[1:3]), seu_out$metric, 
      sc3_out$metric, fict_out$metric, HMRF_out$ari, BayesSpace_out$ari, 
      sim_dat[[3]])
    write.table(t(res_out), file = outfile, quote = F, col.names = F, 
      row.names = F, append = T)
  }, mc.cores = cores)
}

args <- commandArgs(trailingOnly = T)
verbose <- run_main(outdir = "~/BASS_pjt/11_revision/2_rare_type/results", cores = 10)
print(verbose) # debugging

