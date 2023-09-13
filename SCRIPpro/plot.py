# -*- coding: utf-8 -*-
# @Author: Zhanhe Chang

import scanpy as sc
import numpy as np
import pandas as pd
import anndata
import math
import matplotlib
import matplotlib.pyplot as plt


def get_cell_barcode_score(super_leiden,super_leiden_value,adList_obs):
    s = []
    for i in super_leiden:
        origin_leiden = adList_obs[i.split(sep = '_')[0]]
        num = len(origin_leiden.loc[origin_leiden.leiden ==i.split(sep = '_')[1]].index)
        s.extend(list(zip([pd.DataFrame(super_leiden_value).loc[i,:]]*num,origin_leiden.loc[origin_leiden.leiden ==i.split(sep = '_')[1]].index)))
    return s

def plot_spatial(TF, TF_matrix, adList_obs, adata, spot_size, fig_path):
    plot_basis='spatial'
    s = get_cell_barcode_score(TF_matrix.loc[:,TF].index,TF_matrix.loc[:,TF],adList_obs)
    plt.figure(figsize=(10, 7), dpi=80)
    f, ax = plt.subplots()

    spot = sorted(s,key = lambda x: x[0][0],reverse=False)

    init_points = []
    value = []
    for _, cell in spot:
        init_points.append(cell)
        value.append(_)
        
    points = ax.scatter(adata[init_points].obsm[plot_basis][:, 0], adata[init_points].obsm[plot_basis][:, 1], c=value, s = spot_size,cmap="PuRd",edgecolors = None)
    plt.scatter(adata[init_points].obsm[plot_basis][:, 0], adata[init_points].obsm[plot_basis][:, 1], c=value, s = spot_size,cmap="PuRd",edgecolors = None)
    f.colorbar(points)
    ax.grid(False)
    # # ax = plt.gca()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    plt.savefig(fig_path)
