# Author: Zheng Li
# Date: 2021-07-28
# SpaGCN function for simulation
import anndata as ad
import SpaGCN as spg
import pandas as pd
import numpy as np
import scanpy as sc
import random
import torch
import os

def run_SpaGCN_py(sim_cnt, info, R, p = 0.5):
  adata = ad.AnnData(sim_cnt, obs = info)
  # calculate adjacency matrix
  x_array = adata.obs["x"].tolist()
  y_array = adata.obs["y"].tolist()
  adj = spg.calculate_adj_matrix(x = x_array, y = y_array, histology = False)

  # Expression data preprocessing
  adata.var_names_make_unique()
  spg.prefilter_genes(adata)
  spg.prefilter_specialgenes(adata)
  sc.pp.normalize_per_cell(adata)
  sc.pp.log1p(adata)

  # set hyper-parameters
  l = spg.search_l(p, adj, start = 0.01, end = 1000, tol = 0.01, max_run = 100)
  n_clusters = R
  r_seed = t_seed = n_seed = 0
  res = spg.search_res(adata, adj, l, n_clusters, 
    start = 0.7, step = 0.1, tol = 5e-3, lr = 0.05, 
    max_epochs = 20, r_seed = r_seed, t_seed = t_seed, 
    n_seed = n_seed)

  # run spaGCN
  clf = spg.SpaGCN()
  clf.set_l(l)
  random.seed(r_seed)
  torch.manual_seed(t_seed)
  np.random.seed(n_seed)
  clf.train(adata, adj, init_spa = True, init = "louvain", res = res, 
    tol = 5e-3, lr = 0.05, max_epochs = 200)
  y_pred, prob = clf.predict()
  adata.obs["pred"]= y_pred
  adata.obs["pred"] = adata.obs["pred"].astype('category')

  # cluster refinement
  adj_2d = spg.calculate_adj_matrix(x = x_array, y = y_array, histology = False)
  refined_pred = spg.refine(sample_id = adata.obs.index.tolist(), 
    pred = adata.obs["pred"].tolist(), dis = adj_2d, shape = "hexagon")
  adata.obs["refined_pred"] = refined_pred
  adata.obs["refined_pred"] = adata.obs["refined_pred"].astype('category')

  return(adata.obs[["refined_pred"]])