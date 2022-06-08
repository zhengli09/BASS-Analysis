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
  library(GapClust)
  library(mclust)
  library(tidyverse)
  library(scater)
  library(gtools)
  library(splatter)
  library(reticulate)
  library(CVXR)
})


#' Concatenate two strings
"%&%" <- function(x, y) paste0(x, y)


#' Cluster cells with Seurat
seu_cluster <- function(sim_dat, C, resolutions = seq(0.1, 4, by = 0.1))
{
  cnts <- do.call(cbind, sim_dat[[1]])
  colnames(cnts) <- "Cell" %&% 1:ncol(cnts)
  seu <- CreateSeuratObject(counts = cnts, min.cells = 1, min.features = 1)
  seu <- NormalizeData(seu, verbose = F)
  seu <- ScaleData(seu, features = rownames(seu), verbose = F)
  seu <- RunPCA(seu, features = rownames(seu), verbose = F)
  seu <- FindNeighbors(seu, dims = 1:20, verbose = F)
  seu <- FindClusters(seu, resolution = resolutions, verbose = F)
  kest <- sapply("RNA_snn_res." %&% resolutions, function(res){
    length(unique(seu[[res, drop = T]]))
  })
  out_res <- names(which.min(abs(kest - C)))
  c_est <- as.numeric(as.character(seu[[out_res, drop = T]])) + 1
  # evaluation
  ctrue <-  unlist(sim_dat[[2]])
  ari <- eval_ARI(c_est, ctrue)
  F1 <- eval_F1(c_est, ctrue)
  MCC <- eval_MCC(c_est, ctrue)
  C_est <- length(unique(c_est))
  clust_props <- calc_clust_prop(c_est, ctrue)
  
  list(metric = c(ari, F1, MCC), C_est = C_est, clust_props = clust_props)
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
    out <- list(metric = rep("error", 9), C_est = "error", 
      clust_props = "error")
  } else{
    if(ncol(x) > 5000) sce <- sc3_run_svm(sce, ks = C)
    c_est <- colData(sce)[, "sc3_" %&% C %&% "_clusters"]
    ctrue <- unlist(sim_dat[[2]])
    ari <- eval_ARI(c_est, ctrue)
    F1 <- eval_F1(c_est, ctrue)
    MCC <- eval_MCC(c_est, ctrue)
    C_est <- length(unique(c_est))
    clust_props <- calc_clust_prop(c_est, ctrue)
    out <- list(metric = c(ari, F1, MCC), C_est = C_est, 
                clust_props = clust_props)
  }
  return(out)
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
  Sys.setenv(MPLBACKEND = 'Agg')
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
    out <- list(metric = rep("error", 9), C_est = "error", 
                clust_props = "error")
  } else{
    c_est <- c_est$V1+1
    ctrue <- unlist(sim_dat[[2]])
    ari <- eval_ARI(c_est, ctrue)
    F1 <- eval_F1(c_est, ctrue)
    MCC <- eval_MCC(c_est, ctrue)
    C_est <- length(unique(c_est))
    clust_props <- calc_clust_prop(c_est, ctrue)
    out <- list(metric = c(ari, F1, MCC), C_est = C_est, 
                clust_props = clust_props)
  }
  # remove tmp files
  unlink(tmp_fder, recursive = T)
  return(out)
}


#' Match spatial domain labels and cell type labels 
#' in the estimated proportion matrix (pi_est) with 
#' labels in the true proportion matrix
match_pi <- function(pi_est, pi_true, z_est, z_true)
{
  C <- nrow(pi_true)
  R <- ncol(pi_true)
  # 1. Find spatial domain label correspondence
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


#' Map spatial domain label to cell types in that domain
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


#' Produce the confusion matrix
#' One difficulty for constructing the confusion matrix is that 
#' the labeling from the clustering and truth may not match.
#' We get around this difficulty by identifying the confusion matrix 
#' such that the number of true positives (sum of diagonal values)
#' is maximized.
#' Note the number of clusters in est_labels needs to be smaller than
#' or equal to the number of clusters in true_labels
confusion <- function(est_labels, true_labels)
{
  C <- length(unique(true_labels))
  est_labels <- factor(est_labels, 1:C)
  A <- matrix(table(est_labels, true_labels), C, C)
  P <- Variable(C, C, boolean = T) # row permutation matrix
  Y <- Variable(C, C)
  
  problem <- Problem(Maximize(matrix_trace(Y)), list(Y == P %*% A, 
    sum_entries(P, axis = 1) == 1, sum_entries(P, axis = 2) == 1))
  result <- solve(problem)
  P <- result$getValue(P)
  confusion_mtx <- P %*% A
}


#' Evaluate the F1-score
#' Refer to paper Chicco and Jurman, BMC Genomics, 2020
eval_F1 <- function(est_labels, true_labels)
{
  C <- length(unique(true_labels))
  confusion_mtx <- suppressMessages(confusion(est_labels, true_labels))
  TP <- diag(confusion_mtx)
  FP <- apply(confusion_mtx, 1, sum) - diag(confusion_mtx)
  FN <- apply(confusion_mtx, 2, sum) - diag(confusion_mtx)
  F1 <- 2 * TP / (2 * TP + FP + FN)
  c(all_F1 = mean(F1), major_F1 = mean(F1[1:4]), rare_F1 = mean(F1[5:C]))
}


#' Evaluate Matthew correlation coefficient (MCC) score
#' Refer to paper Chicco and Jurman, BMC Genomics, 2020
eval_MCC <- function(est_labels, true_labels)
{
  C <- length(unique(true_labels))
  confusion_mtx <- suppressMessages(confusion(est_labels, true_labels))
  TP <- diag(confusion_mtx)
  FP <- apply(confusion_mtx, 1, sum) - diag(confusion_mtx)
  FN <- apply(confusion_mtx, 2, sum) - diag(confusion_mtx)
  TN <- sum(confusion_mtx) - TP - FP - FN
  MCC <- (TP * TN - FP * FN) / 
    sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  MCC[is.na(MCC)] <- 0
  c(all_MCC = mean(MCC), major_MCC = mean(MCC[1:4]), rare_MCC = mean(MCC[5:C]))
}


#' Evaluate the adjusted Random Index(ARI)
eval_ARI <- function(est_labels, true_labels)
{
  major_cells <- true_labels %in% 1:4
  rare_cells <- !major_cells
  all_ari <- adjustedRandIndex(est_labels, true_labels)
  major_ari <- adjustedRandIndex(
    est_labels[major_cells], true_labels[major_cells])
  rare_ari <- adjustedRandIndex(
    est_labels[rare_cells], true_labels[rare_cells])
  c(all_ari = all_ari, major_ari = major_ari, rare_ari = rare_ari)
}


#' Reorder the labeling of cell type clusters/spatial domains to
#' cell types/spatial domains 1, 2, 3, 4 and then order the remaining
#' redundant clusters/domains based on their cluster/domain size from
#' large to small. Finally, calculate the proportion. Note this quantity
#' is specifically used for evaluation in the case of a mis-specified
#' number of cell types/spatial domains.
calc_clust_prop <- function(est_labels, true_labels)
{
  C <- length(unique(true_labels))
  C_est <- length(unique(est_labels))
  C_min <- min(C, C_est)
  C_max <- max(C, C_est)
  mtx <- matrix(table(est_labels, true_labels), C_est, C)
  P <- Variable(C_est, C_est, boolean = T) # row permutation matrix
  Y <- Variable(C_est, C)
  problem <- suppressMessages(Problem(
    Maximize(matrix_trace(Y[1:C_min, 1:C_min])), 
    list(Y == P %*% mtx, sum_entries(P, axis = 1) == 1, 
      sum_entries(P, axis = 2) == 1)))
  result <- suppressMessages(solve(problem))
  P <- result$getValue(P)
  mtx <- P %*% mtx
  order_idx <- order(apply(mtx, 1, sum)[-c(1:C_min)], decreasing = T)
  mtx <- mtx[c(1:C_min, order_idx+C_min), , drop = F]
  props <- apply(mtx, 1, sum) / sum(mtx)
  names(props) <- 1:C_est
  paste(props, collapse = ",")
}


#' Generate simulated spatial transcriptomic data with splatter package
#'
#' Spatial domain labels are based on STARmap (20180417_BZ5_control) data
#' and annotated based on the marker genes. Spatial domain labels for 
#' section 1 to section 9 are generated based on the annotated labels.
#'
#' @param starmap STARmap data for inferring simulation parameters and for 
#'  providing realistic cell coordinates and spatial domains.
#' @param scenario Simulation scenario (refer to the manuscript for details)
#' @param rare_dist Rare cell types either randomly ditributed across the entire
#'  tissue section or located in certain spatial domians. 
#' @param C Number of cell types that consist of 4 major cell types and (C-4)
#'  rare cell types. Rare cell types together consist of 30% of the total cell 
#'  population and randomly distribued across the entire tissue region. 
#'  Each rare cell type consist of 30%/(C-4) of the total cell population.
#' @param J Number of genes.
#' @param I Number of randomly selected genes among the J genes to remain in
#'   the final dataset.
#' @param L Number of tissue sections.
#' @param batch_facLoc Batch factor location (refer to the splatter package). 
#'  Zero if there is no batch effect.
#' @param de_prop Probability that a gene will be selected to be differentially 
#'  expressed for each cell type (refer to splatter package).
#' @param de_facLoc DE factor location (refer to splatter package).
#' @param de_facScale DE factor scale (refer to splatter package).
#' @param sim_seed Random seed.
#' @param debug Output DE genes for each cell type.
simu <- function(
  starmap,
  scenario, 
  rare_dist = c("random", "spatial"),
  C, 
  J,
  I = NULL,
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
  group_prob <- if(C == 4){
    rep(0.25, 4)
  } else if(C > 4){
    c(rep(0.7 / 4, 4), rep(1 / (C - 4), C - 4) * 0.3)
  }
  params <- setParams(
    init_params, 
    batchCells = rep(3 * N, L), # 3N here represents a large number such that
                                # we have sufficient cells of each type to be 
                                # allocated to the spatial transcriptomics data
    batch.rmEffect = noBatch,
    batch.facLoc = batch_facLoc,
    nGenes = J,
    group.prob = group_prob,
    out.prob = 0,
    de.prob = de_prop,
    de.facLoc = de_facLoc,
    de.facScale = de_facScale,
    seed = sim_seed)
  sim_groups <- splatSimulate(
    params = params, 
    method = "groups", 
    verbose = FALSE)
  if(!is.null(I)){
    # down-sample J genes to I genes
    keep_idx <- sort(sample(1:J, I, replace = F))
    sim_groups <- sim_groups[keep_idx, ]
  } else{
    I <- J
  }
  # remove cells having no expressed genes
  idx_zerosize <- apply(counts(sim_groups), MARGIN = 2, sum) == 0
  sim_groups <- sim_groups[, !idx_zerosize]

  # 2.parse proportion of cells types in each spatial domain
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
    # 3.set spatial domain labels
    if(l == 1){
      ztrue <- info$z
    } else{
      ztrue <- info[["slice" %&% (l-1)]]
    }
    
    # 4.generate cell types
    c_simu[[l]] <- rep(NA, length(ztrue))
    if(scenario == 1){
      c_simu[[l]] <- ztrue
    } else if(scenario %in% c(2, 3, 4)){
      for(z in unique(info$z))
      {
        zi_idx <- ztrue == z
        c_simu[[l]][zi_idx] <- sample(map_z2c(z), sum(zi_idx), prob = prop, 
          replace = T)
      }
    } else if(scenario == 5){ # more challenging scenario with rare cell types
      # assign rare cell types
      if(rare_dist[1] == "spatial"){
        rare_idx <- rep(NA, length(ztrue))
        nrare_each <- round(0.3 / (C - 4) * length(ztrue))
        for(c in 5:C)
        {
          repeat{
            zi_idx <- sample(unique(info$z), 1) == ztrue
            zi_noAlloc_idx <- zi_idx & is.na(rare_idx)
            if(sum(zi_noAlloc_idx) >= nrare_each) break
          }
          ci_idx <- sample(which(zi_noAlloc_idx), nrare_each, replace = F)
          c_simu[[l]][ci_idx] <- c
          rare_idx[ci_idx] <- TRUE
        }
        rare_idx[is.na(rare_idx)] <- FALSE
      } else if(rare_dist[1] == "random"){
        rare_idx <- sample(c(TRUE, FALSE), length(ztrue), prob = c(0.3, 0.7), 
                           replace = T)
        c_simu[[l]][rare_idx] <- sample(5:C, sum(rare_idx), replace = T)
      }
      major_idx <- !rare_idx
      
      # assign major cell types using proportion of 
      # cell types in scenario 3
      for(z in unique(info$z))
      {
        zi_idx <- ztrue == z
        zi_major_idx <- zi_idx & major_idx
        c_simu[[l]][zi_major_idx] <- sample(map_z2c(z), sum(zi_major_idx), 
          prob = c(0.5, 0.25, 0.25), replace = T)
      }
    }
  
    # 5.assign count data
    groups <- as.data.frame(colData(sim_groups)) %>% 
      filter(Batch == "Batch" %&% l)
    sim_cnt[[l]] <- array(NA, c(I, N))
    for(c in 1:C)
    {
      c_size <- sum(c_simu[[l]] == c)
      c_cells <- groups$Cell[grepl("Group" %&% c %&% "$", groups$Group)]
      cells_select <- sample(as.character(c_cells), c_size, replace = F)
      sim_cnt[[l]][, c_simu[[l]] == c] <- 
        as.matrix(counts(sim_groups)[, cells_select])
    }
    colnames(sim_cnt[[l]]) <- "Cell" %&% 1:N
    rownames(sim_cnt[[l]]) <- "Gene" %&% 1:I
  }

  # return DE genes for each group
  if(debug){
    de_genes <- list()
    for(c in 1:C)
    {
      de_genes[[c]] <- rownames(rowData(sim_groups))[
        rowData(sim_groups)[, "DEFacGroup" %&% c] != 1]
    }
    return(list(sim_cnt, c_simu, sim_seed, de_genes))
  }

  return(list(sim_cnt, c_simu, sim_seed))
}


#' Run BASS algorithm
run_BASS <- function(sim_dat, xy, beta_method, beta, C, R,
  init_method)
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
    init_method = init_method, beta_method = "SW", 
    beta = beta, tol = 1e-4, burnin = 2000, nsample = 2000)
  BASS <- BASS.preprocess(BASS)
  BASS <- BASS.run(BASS)
  BASS <- BASS.postprocess(BASS)

  ctrue <- unlist(sim_dat[[2]])
  ztrue <- unlist(info[, 4:(L+3)])
  pi_true <- table(ctrue, ztrue)
  pi_true <- pi_true %*% diag(1 / apply(pi_true, 2, sum))
  c_est <- unlist(BASS@results$c)
  z_est <- unlist(BASS@results$z)
  pi_est <- BASS@results$pi

  # evaluation
  c_ari <- eval_ARI(c_est, ctrue)
  F1 <- eval_F1(c_est, ctrue)
  MCC <- eval_MCC(c_est, ctrue)
  z_ari <- adjustedRandIndex(ztrue, z_est)
  if(C != 4 | R != 4){ # only evaluate in the main simulation
    pi_est <- NA
    mse_pi <- NA
  } else{
    pi_est <- match_pi(pi_est, pi_true, z_est, ztrue)$pi
    mse_pi <- sqrt(sum((pi_est - pi_true)^2))
  }
  C_est <- length(unique(c_est))
  R_est <- length(unique(z_est))
  c_clust_prop <- calc_clust_prop(c_est, ctrue)
  z_clust_prop <- calc_clust_prop(z_est, ztrue)
  
  output <- list(
    c_ari = c_ari,
    c_F1 = F1, 
    c_MCC = MCC,
    C_est = C_est,
    c_clust_prop = c_clust_prop,
    z_ari = z_ari,
    R_est = R_est,
    z_clust_prop = z_clust_prop,
    pi_est = pi_est, 
    mse_pi = mse_pi,
    beta = BASS@results$beta)
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
  R_est <- list()
  z_clust_prop <- list()
  for(beta in allbetas)
  {
    beta <- beta %&% ".0"
    probs <- tryCatch(
      read.table(res_dir %&% "/ftest.beta." %&% beta %&% ".unnormprob.txt", 
        row.names = 1),
      error = function(e) "error",
      warning = function(w) "warning")

    if(!is.data.frame(probs)){
      out <- list(ari = rep("error", length(allbetas)), 
        R_est = rep("error", length(allbetas)),
        z_clust_prop = rep("error", length(allbetas)))
      return(out)
    } else{
      z_est <- apply(probs, MARGIN = 1, which.max)
      ari[[beta]] <- adjustedRandIndex(ztrue, z_est)
      R_est[[beta]] <- length(unique(z_est))
      z_clust_prop[[beta]] <- calc_clust_prop(z_est, ztrue)
    }
  }
  unlink(results_folder, recursive = T)
  out <- list(ari = ari, R_est = R_est, z_clust_prop = z_clust_prop)
  return(out)
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
  ari <- adjustedRandIndex(ztrue, z_est)
  R_est <- length(unique(z_est))
  z_clust_prop <- calc_clust_prop(z_est, ztrue)
  out <- list(ari = ari, R_est = R_est, z_clust_prop = z_clust_prop)
  return(out)
}


#' run_SpGCN
run_SpaGCN <- function(sim_dat, xy, ztrue, R, p = 0.5)
{
  sim_cnt <- sim_dat[[1]][[1]]
  xy <- data.frame(xy)
  z_est <- run_SpaGCN_py(t(sim_cnt), xy, R, p)$refined_pred
  ari <- adjustedRandIndex(ztrue, z_est)
  R_est <- length(unique(z_est))
  z_clust_prop <- calc_clust_prop(z_est, ztrue)
  out <- list(ari = ari, R_est = R_est, z_clust_prop = z_clust_prop)
  return(out)
}
