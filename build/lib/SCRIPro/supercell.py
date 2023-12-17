import math
import multiprocessing as mp
import random
import warnings
from collections import Counter
from multiprocessing import Pool
from scipy.spatial.distance import cdist
import anndata as ad
import h5py
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from lisa import FromGenes

def preprocess_data(sub_adata, n_genes=2000, npcs=40, percent_cells=0.7):
    sc.pp.normalize_total(sub_adata, target_sum=1e4)
    sc.pp.log1p(sub_adata)
    sc.pp.highly_variable_genes(sub_adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=n_genes)
    sub_adata = sub_adata[:, sub_adata.var.highly_variable]
    ad_sub = sc.AnnData(sub_adata.X, obs=sub_adata.obs, var=sub_adata.var,dtype='float32')
    
    saved_log_mat = ad_sub.X.copy()
    sc.pp.scale(ad_sub)
    sc.tl.pca(ad_sub, svd_solver='arpack', n_comps=min(npcs, sub_adata.shape[0] - 1))
    sc.pp.neighbors(ad_sub, n_neighbors=10, n_pcs=min(npcs, sub_adata.shape[0] - 1))
    sc.tl.umap(ad_sub)
    ad_sub.X = saved_log_mat

    return ad_sub
def merge_cells(ads, variable_to_merge='leiden'):
    clust_assignment = ads.obs.loc[:, variable_to_merge].to_numpy()
    clusts = np.unique(clust_assignment)
    ncells = len(clusts)
    new_mat = np.zeros((ncells, ads.shape[1]))
    for i, c in enumerate(clusts):
        locs = np.where(clust_assignment == c)[0]
        new_mat[i, :] = np.mean(ads.X[locs, :], axis=0)
    merged_cells = sc.AnnData(new_mat, var=ads.var,dtype='float32')
    return merged_cells
def get_leiden_based_on_ncell(ad_sub, num_cells, verbose):
    resolutions = np.arange(0.1, 1000, 0.1)
    vec_length = len(resolutions)
    iter_ = int(vec_length / 2)
    last_iter = 0
    if ad_sub.shape[0] > num_cells:

        # binary search to find optimal resolution
        while True:
            vec_length = abs(iter_ - last_iter)
            sc.tl.leiden(ad_sub, resolution=resolutions[iter_])

            if abs(iter_ - last_iter) <= 1:
                break
	
            ncells_in_merged = len(np.unique(ad_sub.obs.leiden.to_numpy().ravel()))

            last_iter = iter_

            if ncells_in_merged < num_cells:
                iter_ += int(vec_length / 2)
            elif ncells_in_merged >= num_cells:
                iter_ -= int(vec_length / 2)
	
        
        sc.tl.leiden(ad_sub, resolution=resolutions[iter_])
    else:
        sc.tl.leiden(ad_sub,resolution=0.1)
    if verbose:
        print('There are ',
              len(np.unique(ad_sub.obs.leiden.to_numpy().ravel())),
              ' supercells')

    return ad_sub
def get_merged_dataset(adata_all, obs):
    
    num_super_cells = 0
    
    num_super_cells += len(np.unique(obs.leiden))

    
    new_data = np.zeros((num_super_cells, adata_all.X.shape[1]))

    new_celltypes = []
    current_loc = 0
    sub_cell = adata_all[np.isin(adata_all.obs.index,obs.index)]
    num_current_cell = len(np.unique(obs.leiden))
    sub_cell.obs = obs
    merged = merge_cells(sub_cell)
    new_data[current_loc:(current_loc + num_current_cell), :] = merged.X
    current_loc += num_current_cell

    new_celltypes = np.array(new_celltypes)

    all_merged = sc.AnnData(new_data, obs=new_celltypes, var=adata_all.var,dtype='float32')

    all_merged.X = np.round(all_merged.X)

    return all_merged
def supercell_pipeline(adata, ngenes=2000, npcs=40,cell_num=50,min_cell=30,verbose=True):
    ncell=adata.shape[0]/cell_num
    saved_counts = adata.X.copy()
    # Run PCA and find nearest neighbors
    if verbose: print('preprocessing data...')
    sub_cells = preprocess_data(adata, n_genes=ngenes, npcs=npcs)
    
    if verbose: print('finding optimal resolution...')
    temp = get_leiden_based_on_ncell(sub_cells, ncell, verbose)

    adata.X = saved_counts
    if temp.obs.leiden.value_counts()[0]>min_cell and temp.obs.leiden.value_counts().shape[0]>2:
        clusters = temp.obs.leiden.cat.categories
        centroids = pd.DataFrame({cluster: np.mean(temp[temp.obs.index[np.where(temp.obs.leiden == cluster)]].obsm['X_umap'],axis=0) for cluster in clusters}).T
        small_clusters = temp.obs.leiden.value_counts()[temp.obs.leiden.value_counts() < min_cell].index
        replace_dict= {}
        for cluster in small_clusters:
            distances = cdist(centroids.loc[[cluster]], centroids)
            distances = distances[0][0:len(clusters)-len(small_clusters)]
            nearest_cluster = clusters[np.argmin(distances)]
            replace_dict[cluster] = nearest_cluster
        temp.obs.leiden = temp.obs.leiden.replace(replace_dict)
        
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    if verbose: print('merging cells...')
    merged_data = get_merged_dataset(adata, temp.obs)
    
    return merged_data,temp.obs

def get_supercell(adata,adata_super):
    adList = []
    adList_obs = {}
    for i in set(adata.obs.leiden):
        leiden_index = adata.obs.loc[adata.obs.leiden == i].index
        sub_cluster = adata_super[leiden_index]
        merged_data,obs = supercell_pipeline(sub_cluster,cell_num=50,verbose=False)
        
        merged_data_index = [i + "_" + str(j) for j in range(0,merged_data.shape[0])]
        merged_data.obs = pd.DataFrame(index = merged_data_index)
        adList.append(merged_data)
        adList_obs[i] = obs
    ad_all = ad.concat(adList,join='outer')
    ad_all = ad_all.to_df()
    adata.obs.insert(loc=6, column='new_leiden', value=3)
    for i in adList_obs.keys():
        for k in adList_obs[i].index:
            adata.obs.loc[k,'new_leiden'] = i+"_"+adList_obs[i].loc[k,'leiden']  
    return ad_all
    