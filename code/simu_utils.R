# Author: Zheng Li
# Date: 2021-05-19
# Purpose:
# This script includes all functions to conduct the simulation study for BASS

suppressPackageStartupMessages({
  library(BASS)
  library(Giotto)
  library(BayesSpace)
  library(Seurat)
  library(SC3)
  library(mclust)
  library(tidyverse)
  library(scater)
  library(gtools)
  library(splatter)
  library(reticulate)
})


#' Concatenate two strings
"%&%" <- function(x, y) paste0(x, y)


#' Cluster cells with Seurat
seu_cluster <- function(sim_dat, resolutions)
{
  cnts <- do.call(cbind, sim_dat[[1]])
  colnames(cnts) <- "Cell" %&% 1:ncol(cnts)
  seu <- CreateSeuratObject(counts = cnts, min.cells = 1, min.features = 1)
  seu <- NormalizeData(seu, verbose = F)
  seu <- ScaleData(seu, features = rownames(seu), verbose = F)
  seu <- RunPCA(seu, features = rownames(seu), verbose = F)
  seu <- FindNeighbors(seu, dims = 1:20, verbose = F)
  seu <- FindClusters(seu, resolution = resolutions, verbose = F)
  ari <- sapply("RNA_snn_res." %&% resolutions, function(res){
    adjustedRandIndex(seu[[res, drop = T]], unlist(sim_dat[[2]]))
    })
  kest <- sapply("RNA_snn_res." %&% resolutions, function(res){
    length(unique(seu[[res, drop = T]]))
    })
  seu_out <- paste(ari, kest, sep = ",")
}


#' Cluster cells with Gaussian mixture models (GMM)
gmm_cluster <- function(sim_dat, C, model)
{
  cnts <- do.call(cbind, sim_dat[[1]])
  x <- scater::normalizeCounts(cnts, log = TRUE)
  x <- scater::calculatePCA(x, scale = T)
  gmm <- Mclust(x[, 1:20], G = C, modelNames = model, verbose = T)
  ari <- adjustedRandIndex(gmm$classification, unlist(sim_dat[[2]]))
}


#' Cluster cells with k-means algorithm
km_cluster <- function(sim_dat, C)
{
  cnts <- do.call(cbind, sim_dat[[1]])
  x <- scater::normalizeCounts(cnts, log = TRUE)
  x <- scater::calculatePCA(x, scale = T)
  km <- kmeans(x[, 1:20], C, nstart = 25)
  c_est <- km$cluster
  ari <- adjustedRandIndex(c_est, unlist(sim_dat[[2]]))
}


#' Cluster cells with SC3
sc3_cluster <- function(sim_dat, C, seed = 0)
{
  set.seed(seed)
  cnts <- do.call(cbind, sim_dat[[1]])
  x <- scater::normalizeCounts(cnts, log = TRUE)
  sce <- SingleCellExperiment(assays = list(counts = cnts, logcounts = x))
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- tryCatch({
    sc3(sce, ks = C, n_cores = 1)}, 
    error = function(e){
      print(e)
      return("error")}
    )
  if(is.character(sce)){
    return(sce)
  } else{
    if(ncol(x) > 5000) sce <- sc3_run_svm(sce, ks = C)
    c_est <- colData(sce)[, "sc3_" %&% C %&% "_clusters"]
    ari <- adjustedRandIndex(c_est, unlist(sim_dat[[2]]))
    return(ari)
  }
}

#' Cluster cells with FICT
#' FICT input files:
#' 1.expression, NxJ
#' 2.coordinates,
fict_cluster <- function(sim_dat, xy, C, case, rep)
{
  # The following environmental variables need to be set 
  # to run FICT successfully
  PYTHONPATH <- paste(
    "/net/mulan/home/zlisph/BASS_pjt/9_FICT/FICT-SAMPLE/FICT/",
    "/net/mulan/home/zlisph/BASS_pjt/9_FICT/FICT-SAMPLE/GECT/", sep = ":")
  Sys.setenv(PYTHONPATH = PYTHONPATH)
  Sys.setenv(MPLBACKEND='Agg')
  fict <- "~/softwares/miniconda3/bin/fict"
  tmp_fder <- "~/BASS_pjt/2_simu/fict_tmp/case" %&% case %&% "/rep" %&% rep
  if(!file.exists(tmp_fder)) dir.create(tmp_fder, recursive = T)
  # prepare input
  # 1.coordinates
  write.table(xy, file = tmp_fder %&% "/simu.coordinates", 
    col.names = F, quote = F)
  # 2.expression
  cnts <- do.call(cbind, sim_dat[[1]])
  ge <- scater::normalizeCounts(cnts, log = TRUE)
  ge <- t(ge) # ncell x ngene
  write.table(ge, file = tmp_fder %&% "/simu.expression",
    col.names = F, quote = F)

  # run FICT
  # Note:
  # by default hidden = 20, size of denoise auto-encoder
  cmd <- paste(fict, "-p", tmp_fder %&% "/simu", "-o", 
    tmp_fder %&% "/out/", "--n_type", C)
  system(cmd, ignore.stdout = T, ignore.stderr = T)
  
  # evaluate
  # occasionally have error: numpy.linalg.LinAlgError: Internal Error.
  c_est <- tryCatch({
    read.table(tmp_fder %&% "/out/cluster_result.csv")
    }, error = function(e) e)
  if(!is.data.frame(c_est)){
    return("error")
  } else{
    ari <- adjustedRandIndex(c_est$V1, unlist(sim_dat[[2]]))
    return(ari)
  }
  # remove tmp files
  unlink(tmp_fder, recursive = T)
}


#' Match spatial domain labels and cell type labels 
#' in the estimated proportion matrix (pi_est) with 
#' labels in the true proportion matrix
match_pi <- function(pi_est, pi_true, z_est, z_true)
{
  C <- nrow(pi_true)
  R <- ncol(pi_true)
  # 1. Find tissue structure label correspondence
  perms_R <- permutations(n = R, r = R)
  accur <- rep(NA, nrow(perms_R))
  for(i in 1:nrow(perms_R))
  {
    perm <- perms_R[i, ]
    z_true_perm <- case_when(
      z_true == 1 ~ perm[1], 
      z_true == 2 ~ perm[2], 
      z_true == 3 ~ perm[3], 
      z_true == 4 ~ perm[4])
    accur[i] <- mean(z_est == z_true_perm)
  }
  perm_final_R <- perms_R[which.max(accur), ]

  # 2. Match spatial domain label in pi_est
  pi_est_match <- pi_est[, perm_final_R]

  # 3. Match cell type label
  perms_C <- permutations(n = C, r = C)
  sse <- rep(NA, nrow(perms_C))
  for(i in 1:nrow(perms_C))
  {
    perm <- perms_C[i, ]
    sse[i] <- sum((pi_est_match[perm, ] - pi_true)^2)
  }
  perm_final_C <- perms_C[which.min(sse), ]
  pi_est_match <- pi_est_match[perm_final_C, ]
  out <- list(
    pi = pi_est_match, 
    perm_K = perm_final_R, 
    perm_C = perm_final_C)
}


#' Map tissue structure label to cell types in that region
#' 1 -> (1, 2, 3), 2 -> (2, 3, 4)
#' 3 -> (3, 4, 1), 4 -> (4, 1, 2)
map_z2c <- function(z)
{
  case_when(
    z == 1 ~ c(1, 2, 3),
    z == 2 ~ c(2, 3, 4),
    z == 3 ~ c(3, 4, 1),
    z == 4 ~ c(4, 1, 2)
    )
}


#' Generate simulated spatial transcriptomic data with splatter package
#'
#' Spatial domain labels are based on STARmap (20180417_BZ5_control)
#' and annotated based on marker genes. Spatial domain labels for section 1 
#' to section 9 are generated based on the annotated labels.
#'
#' @param starmap starmap data for inferring simulation parameters, providing
#'                realistic cell coordinates and spatial domains (eL2/3 cells)
#' @param scenario Simulation scenario (refer to the manuscript for details)
#' @param C Number of cell types
#' @param J Number of genes
#' @param L Number of tissue sections
#' @param batch_facLoc Batch factor location (refer to splatter package). Zero 
#'                     if there is no batch effect
#' @param de_prop Probability that a gene will be selected to be differentially 
#'                expressed for each cell type (refer to splatter package)
#' @param de_facLoc DE factor location (refer to splatter package)
#' @param de_facScale DE factor scale (refer to splatter package)
#' @param sim_seed Random seed
#' @param debug Output DE genes for each cell type
simu <- function(
  starmap,
  scenario, 
  C, 
  J, 
  L, 
  batch_facLoc, 
  de_prop, 
  de_facLoc, 
  de_facScale,
  sim_seed, 
  debug = FALSE
  )
{
  cnts <- starmap$cnts
  info <- starmap$info
  N <- nrow(info)
  init_params <- splatEstimate(cnts[, which(info$c == "eL2/3")])

  # reproduce randomness in cell type assignment and 
  # count data selection.
  set.seed(sim_seed)

  # 1.simulate count data
  noBatch <- ifelse(batch_facLoc == 0, TRUE, FALSE)
  params <- setParams(
    init_params, 
    batchCells = rep(3 * N, L),
    batch.rmEffect = noBatch,
    batch.facLoc = batch_facLoc,
    nGenes = J,
    group.prob = rep(1, C) / C,
    out.prob = 0,
    de.prob = de_prop,
    de.facLoc = de_facLoc,
    de.facScale = de_facScale,
    seed = sim_seed)
  sim_groups <- splatSimulate(
    params = params, 
    method = "groups", 
    verbose = FALSE)
  # remove cells having no expressed genes
  idx_zerosize <- apply(counts(sim_groups), MARGIN = 2, sum) == 0
  sim_groups <- sim_groups[, !idx_zerosize]

  # 2.parse proportion of cells types
  if(scenario == 2){
    prop <- c(0.8, 0.1, 0.1)
  } else if(scenario == 3){
    prop <- c(0.5, 0.25, 0.25)
  } else if(scenario == 4){
    prop <- rep(1, 3) / 3
  }

  c_simu <- list()
  sim_cnt <- list()
  for(l in 1:L)
  {
    # 3.generate cell types
    if(l == 1){
      ztrue <- info$z
    } else{
      ztrue <- info[["slice" %&% (l-1)]]
    }
    c_simu[[l]] <- rep(NA, length(ztrue))
    if(scenario == 1){
      c_simu[[l]] <- ztrue
    } else if(scenario %in% c(2, 3, 4)){
      for(z in unique(info$z))
      {
        zi_idx <- ztrue == z
        c_simu[[l]][zi_idx] <- sample(map_z2c(z), sum(zi_idx), prob = prop, replace = T)
      }
    }
  
    # 4.assign count data
    groups <- as.data.frame(colData(sim_groups)) %>% filter(Batch == "Batch" %&% l)
    sim_cnt[[l]] <- array(NA, c(J, N))
    for(c in 1:C)
    {
      c_size <- sum(c_simu[[l]] == c)
      c_cells <- groups$Cell[grepl(c, groups$Group)]
      cells_select <- sample(as.character(c_cells), c_size, replace = F)
      sim_cnt[[l]][, c_simu[[l]] == c] <- as.matrix(counts(sim_groups)[, cells_select])
    }
    colnames(sim_cnt[[l]]) <- "Cell" %&% 1:N
    rownames(sim_cnt[[l]]) <- "Gene" %&% 1:J
  }

  # return DE genes for each group
  if(debug){
    de_genes <- list()
    for(c in 1:C)
    {
      de_genes[[c]] <- rownames(rowData(sim_groups))[rowData(sim_groups)[, "DEFacGroup" %&% c] != 1]
    }
    return(list(sim_cnt, c_simu, sim_seed, de_genes))
  }

  return(list(sim_cnt, c_simu, sim_seed))
}


#' Run BASS algorithm
run_BASS <- function(sim_dat, xy, beta_est_approach, beta, C, R,
  init_method, cov_struc)
{
  L <- length(sim_dat[[1]])
  xys <- lapply(1:L, function(x) xy)
  sim_cnt <- sim_dat[[1]]
  for(l in 1:L)
  {
    colnames(sim_cnt[[l]]) <- rownames(xys[[l]])
  }

  # run algorithms
  BASS <- createBASSObject(sim_cnt, xys, C = C, R = R, 
    beta_est_approach = beta_est_approach, beta = beta,
    init_method = init_method, cov_struc = cov_struc,
    beta_tol = 0.01, burn_in = 5000, samples = 5000)
  BASS <- BASS.preprocess(BASS, doLogNormalize = T, doPCA = T)
  BASS <- BASS.run(BASS)
  BASS <- BASS.postprocess(BASS)

  ctrue <- unlist(sim_dat[[2]])
  ztrue <- unlist(info[, 4:(L+3)])
  pi_true <- table(ctrue, ztrue)
  pi_true <- pi_true %*% diag(1 / apply(pi_true, 2, sum))
  c_est <- BASS@res_postprocess$c_ls
  z_est <- BASS@res_postprocess$z_ls
  pi_est <- BASS@res_postprocess$pi_ls

  c_ari <- adjustedRandIndex(ctrue, unlist(c_est))
  z_ari <- adjustedRandIndex(ztrue, unlist(z_est))
  pi_est <- match_pi(pi_est, pi_true, unlist(z_est), ztrue)$pi
  mse_pi <- sqrt(sum((pi_est - pi_true)^2))

  output <- list(
    c_ari = c_ari, 
    z_ari = z_ari, 
    pi_est = pi_est, 
    mse_pi = mse_pi)
  return(output)
}


#' Run HMRF algorithm
run_HMRF <- function(sim_dat, xy, ztrue, R, case, rep, 
  usePCs = F, SEgenes = NULL, dosearchSEgenes = T)
{
  sim_cnt <- sim_dat[[1]][[1]]
  J <- nrow(sim_cnt)

  # 0.prepare
  my_python_path <- "/net/mulan/home/zlisph//softwares/miniconda3/bin/python"
  results_folder <- "/net/mulan/home/zlisph/BASS_pjt/2_simu/HMRF_tmp/case" %&% 
    case %&% "/rep" %&% rep
  if(!file.exists(results_folder)) dir.create(results_folder, recursive = T)
  instrs <- createGiottoInstructions(python_path = my_python_path)

  # 1.run HMRF
  starmap <- createGiottoObject(
    raw_exprs = as.data.frame(sim_cnt), 
    spatial_locs = xy, 
    instructions = instrs)
  starmap <- filterGiotto(
    gobject = starmap, 
    expression_threshold = 0.5,
    gene_det_in_min_cells = 20, 
    min_det_genes_per_cell = 0)
  starmap <- normalizeGiotto(starmap)
  starmap <- addStatistics(starmap)
  starmap <- createSpatialNetwork(starmap, minimum_k = 2)

  # select SE genes using
  if(dosearchSEgenes){
    genes <- binSpect(starmap, bin_method = 'kmeans')$genes
    genes <- if(length(genes) < 100) genes else genes[1:100]
  } else{
    genes <- rownames((starmap@norm_scaled_expr))
  }

  if(usePCs){
    starmap <- Giotto::runPCA(
      gobject = starmap, 
      genes_to_use = genes, 
      ncp = 20)
    betas <- c(0, 1, 21)
    output_folder <- results_folder %&% "/result_k" %&% R %&% "_pca"
    HMRF_out <- doHMRF(
      gobject = starmap, 
      dim_reduction_to_use = "pca", 
      dimensions_to_use = 1:20, 
      k = R, 
      betas = betas, 
      output_folder = output_folder, 
      overwrite_output = T)
  } else{
    betas <- c(0, 2, 26)
    output_folder <- results_folder %&% "/result_k" %&% R
    HMRF_out <- doHMRF(
      gobject = starmap, 
      expression_values = "scaled", 
      spatial_genes = genes,
      k = R, 
      betas = betas, 
      output_folder = output_folder, 
      overwrite_output = T)
  }

  # 2.evaluate results
  res_dir <- output_folder %&% "/result.spatial.zscore/k_" %&% R
  allbetas <- seq(betas[1], by = betas[2], length = betas[3])
  ari <- list()
  for(beta in allbetas)
  {
    beta <- beta %&% ".0"
    probs <- tryCatch(
      read.table(res_dir %&% "/ftest.beta." %&% beta %&% ".unnormprob.txt", 
        row.names = 1),
      error = function(e) "error",
      warning = function(w) "warning"
      )

    if(!is.data.frame(probs)){
      return(rep("error", length(allbetas)))
    } else{
      z_est <- apply(probs, MARGIN = 1, which.max)
      ari[[beta]] <- adjustedRandIndex(ztrue, z_est)
    }
  }
  unlink(results_folder, recursive = T)
  return(unlist(ari))
}


#' Run BayesSpace
run_BayesSpace <- function(sim_dat, xy, ztrue, R, seed = 0)
{
  set.seed(seed)
  sim_cnt <- sim_dat[[1]][[1]]
  colnames(xy) <- c("row", "col")
  info <- data.frame(xy)
  sce <- SingleCellExperiment(assays = list(counts = sim_cnt), colData = info)
  sce <- spatialPreprocess(sce, n.PCs = 15, n.HVGs = nrow(sce), log.normalize = T)
  sce <- spatialCluster(sce, q = R, d = 15, 
    init.method = "mclust", model = "t",
    nrep = 10000, burn.in = 1000)
  z_est <- colData(sce)$spatial.cluster
  return(adjustedRandIndex(ztrue, z_est))
}


#' run_SpGCN
run_SpaGCN <- function(sim_dat, xy, ztrue, R, p = 0.5)
{
  sim_cnt <- sim_dat[[1]][[1]]
  xy <- data.frame(xy)
  z_est <- run_SpaGCN_py(t(sim_cnt), xy, R, p)$refined_pred
  return(adjustedRandIndex(ztrue, z_est))
}