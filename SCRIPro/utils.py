import math
import random
import os
import warnings
from collections import Counter
from multiprocessing import Pool
from itertools import islice
import anndata as ad
import math
import h5py
import shutil
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
import scanpy as sc
import multiprocessing
import seaborn as sns
from lisa import FromCoverage, FromGenes, FromRegions
from lisa.core.assays import LISA_RP_Assay, delta_RP_wrapper, transform_RP,get_delta_RP,get_deltaRP_activation
from functools import partial
from scipy import stats, sparse
from scipy.stats import kendalltau, pearsonr, spearmanr

def write_dataframe_to_tsv(df_subset, filename):
    series = df_subset.stack().reset_index()
    series = series[series.iloc[:, 2] != 0]
    with open(filename, 'w') as tsvfile:
        for _, row in series.iterrows():
            barcode = row['barcode']
            position_parts = row['level_1'].split('_')
            chr_name = position_parts[0]
            start_pos = position_parts[1]
            end_pos = position_parts[2]
            open_value = row[0]
            tsvfile.write(f"{chr_name}\t{start_pos}\t{end_pos}\t{barcode}\t{open_value}\n")
            

def get_marker_for_group(adata, i, log2fc, pval, gene_list_len):
    tem_gene = list(sc.get.rank_genes_groups_df(adata, group=i, log2fc_min=log2fc, pval_cutoff=pval).sort_values(by = 'logfoldchanges',ascending = False)[0:500].names)
    if len(tem_gene) < gene_list_len:
        return i, None
    else:
        return i, tem_gene

def process_and_merge(df, folder_path, n_num=20):
    num_partitions = min(n_num, len(df))
    dfs = np.array_split(df, num_partitions)
    pool = multiprocessing.Pool(processes=num_partitions)

    
    for i, df_subset in enumerate(dfs):
        filename = f"{folder_path}/output_{i}.tsv"
        pool.apply_async(write_dataframe_to_tsv, args=(df_subset, filename))

    
    pool.close()
    pool.join()

    
    merge_and_delete_tsv(folder_path)

def merge_and_delete_tsv(folder_path):
    
    subprocess.run(f'cat {folder_path}/*.tsv > {folder_path}/merge.tsv', shell=True)
    
    subprocess.run(f'rm {folder_path}/output*.tsv', shell=True)


def cal_tf(input_mat,species1, assays1):
    all_tf_result = {}
    lisa_test2 = FromGenes(species=species1,rp_map='enhanced_10K', assays=assays1,isd_method='chipseq', verbose = 0)  
    for k,v in input_mat.items():
        results, metadata = lisa_test2.predict(v, num_background_genes = 3000, background_strategy = 'regulatory')
        results = pd.DataFrame(results.to_dict())
        all_tf_result[k] = results
    return all_tf_result 


def process_group(group,adata,log,pval):
    return group, list(sc.get.rank_genes_groups_df(adata, group=group, log2fc_min=log, pval_cutoff=pval).sort_values(by='logfoldchanges', ascending=False)[0:500].names)

def get_supercell_fragment(leiden_clusters,base_dir,fragment_file,chunksize = 10000000):
    folder_name = base_dir+"/supercell_fragment"
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)
    else:
        print("文件夹已经存在：", folder_name)
    for i,chunk in enumerate(pd.read_csv(fragment_file, delimiter='\t',comment='#', chunksize=chunksize,header=None, names=['chrom', 'start', 'end', 'cell', 'fragment'])):
        chunk.set_index('cell', inplace=True)
        groups = chunk.groupby(leiden_clusters['new_leiden'])
        for group_name, group_data in groups:
            file_path = os.path.join(base_dir, 'supercell_fragment', f'{group_name}.tsv')
            group_data = group_data[group_data['chrom'].str.startswith('chr')]
            group_data.to_csv(file_path, sep='\t',mode='a', header=False, index=False)
        print(f"Processed chunk {i+1}")
    print('final')

def process_tsv(working_directory, bedGraph,species):
    
    sort_command = 'ls ./*.tsv | xargs -I {} sh -c \'bedtools sort -i "$1" > ./sort_tsv/"$(basename -- "$1")"\' -- {}'
    sort_folder = os.path.join(working_directory, "sort_tsv")
    if not os.path.exists(sort_folder):
        os.mkdir(sort_folder)
    p = subprocess.Popen(sort_command, cwd=working_directory, shell=True)
    p.wait()

    
    merge_command = 'ls *.tsv | xargs -P 4 -n 1 -I {} sh -c \'bedtools merge -d 1000 -c 4 -o sum -i "{}" > "./merge_tsv/{}"\''
    merge_folder = os.path.join(sort_folder, "merge_tsv")
    if not os.path.exists(merge_folder):
        os.mkdir(merge_folder)
    p = subprocess.Popen(merge_command, cwd=sort_folder, shell=True)
    p.wait()

    
    bigwig_command = 'find . -name \'*.tsv\' -type f -print0 | xargs -0 -P 4 -n 1 sh -c \'file=\"$1\"; {0} \"$file\" {1} \"./bigwig/${{file%.*}}.bw\"\' sh'.format(bedGraph,species)
    bigwig_folder = os.path.join(merge_folder, "bigwig")
    if not os.path.exists(bigwig_folder):
        os.mkdir(bigwig_folder)
    p = subprocess.Popen(bigwig_command, cwd=merge_folder, shell=True)
    p.wait()

    
    shutil.move(bigwig_folder, os.path.join(os.getcwd(), "bigwig"))
    shutil.rmtree(working_directory)
    
def process_marker(i,lisa_info,bw_path,rpmap_enhanced,factor_binging,factor_metadata,species):
    cell_info = lisa_info[i]
    gene_mask = cell_info[0]
    label_vector = cell_info[1]
    bigwig_path = bw_path + str(i) + '.bw'
    testbw = FromCoverage.convert_bigwig(bigwig_path, species)

    test_profile = testbw[:, np.newaxis]
    test_profile = test_profile / test_profile.sum() * 1e5
     
    rp_matrix = rpmap_enhanced.dot(test_profile)
    subset_rp_matrix = rp_matrix[gene_mask, :]
    
    bin_mask = np.squeeze(np.array(rpmap_enhanced[gene_mask, :].tocsc().sum(axis=0) > 0))
    
    subset_factor_binding = factor_binging[bin_mask, :]
    subset_rp_map = rpmap_enhanced[gene_mask, :][:, bin_mask]
    subset_accessibility = test_profile[bin_mask, :]
    
    rp_knockout = get_delta_RP(subset_accessibility, subset_factor_binding, subset_rp_map)[:, np.newaxis, :]
    deltaX = get_deltaRP_activation(subset_rp_matrix[:, :, np.newaxis], rp_knockout)
    p_vals = get_delta_RP_p_value(deltaX.transpose(0, 2, 1)[:, :, 0], label_vector)
    factor_metadata_pd = pd.DataFrame(factor_metadata)
    factor_metadata_pd['p_vals'] = p_vals
    return i, factor_metadata_pd


def glue_supercell(combined):
    for i in set(combined.obs.leiden):
        leiden_index = combined.obs.loc[combined.obs.leiden == i].index
        sub_test = combined[combined.obs.leiden == i]
        get_leiden_based_on_ncell(sub_test,50,verbose=False)
        if sub_test.obs.leiden.value_counts()[0]>30 and sub_test.obs.leiden.value_counts().shape[0]>2:
            clusters = sub_test.obs.leiden.cat.categories
            centroids = pd.DataFrame({cluster: np.mean(sub_test[sub_test.obs.index[np.where(sub_test.obs.leiden == cluster)]].obsm['X_glue'],axis=0) for cluster in clusters}).T
            small_clusters = sub_test.obs.leiden.value_counts()[sub_test.obs.leiden.value_counts() < 30].index
            replace_dict= {}
            for cluster in small_clusters:
                distances = cdist(centroids.loc[[cluster]], centroids)
                distances = distances[0][0:len(clusters)-len(small_clusters)]
                nearest_cluster = clusters[np.argmin(distances)]
                replace_dict[cluster] = nearest_cluster
            sub_test.obs.leiden = sub_test.obs.leiden.replace(replace_dict)
            sub_test.obs['leiden'] = i + "_" + sub_test.obs['leiden'].astype(str)
        else:
            sub_test.obs['leiden'] = i + "_0"
        combined.obs.loc[leiden_index,'new_leiden'] = sub_test.obs['leiden']
def get_delta_RP_p_value(gene_TF_scores, label_vector):
        '''
        gene_TF_scores: gene x TF, model output of delta-RP matrix. more purturbation of genes of interest correspond with higher delta regulation score
        '''
        
        query_delta = gene_TF_scores[label_vector.astype(bool)]
        background_delta = gene_TF_scores[~label_vector.astype(bool)]

        
        test_parameters = list(zip(query_delta.T, background_delta.T))

        p_vals = [
            mannu_test_function((q,b)) for q,b in test_parameters
        ]

        _, p_values = list(zip(*p_vals))

        return p_values
    
def mannu_test_function(x):
    query, background = x
    try:
        return stats.mannwhitneyu(query, background, alternative = 'greater')
    
    except ValueError:
        
        return (None, 1.0)

def get_quert_gene(chip_matrix,target_h5,gene_name):
    tem = {}
    for i in set(chip_matrix.values.flatten()):
        tem[i]=list(set(target_h5[i][target_h5[i]>2].index).intersection(set(gene_name)))
    return tem

def maxmin(dataframe):
    return (dataframe - dataframe.min())/(dataframe.max() - dataframe.min())

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