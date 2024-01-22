import math
import random
import warnings
from collections import Counter
from multiprocessing import Pool
from itertools import islice
import anndata as ad
import math
import h5py
from tqdm import tqdm
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from lisa import FromGenes
from tqdm.contrib.concurrent import process_map
from .supercell import *
from .utils import *
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from scipy.stats import kendalltau, pearsonr, spearmanr
class Ori_Data():
    def __init__(self,adata,Cell_num,use_glue = False,cluster_method='leiden'):
        self.adata = adata[:, ~adata.var_names.str.contains('\.|\-')]
        self.Cell_num = Cell_num
        self.cluster_method = cluster_method
        if use_glue == False:
            #self._fliter()
            self._supercell(cluster_method=self.cluster_method)
            self.super_gene_exp = self.ad_all.copy()
            self.super_gene_mean = self.super_gene_exp.mean()
            self.super_gene_std = self.super_gene_exp.std()
            sc.tl.rank_genes_groups(self.adata,'new_leiden', method='t-test')
            self.cellgroup = pd.DataFrame(self.adata.obs.iloc[:,-1]) 
    def get_glue_cluster(self,rna_leiden_clusters):
        self.adata.obs['new_leiden'] = rna_leiden_clusters[self.adata.obs.index]
        self.ad_all=pd.DataFrame(self.adata.raw.X,index=self.adata.raw.obs_names,columns=self.adata.raw.var_names).groupby(self.adata.obs['new_leiden']).mean()
        self.super_gene_exp = self.ad_all.copy()
        self.super_gene_mean = self.super_gene_exp.mean()
        self.super_gene_std = self.super_gene_exp.std()
        sc.tl.rank_genes_groups(self.adata,'new_leiden', method='t-test')
        #self.cellgroup = pd.DataFrame(self.adata.obs.iloc[:,-1])

    def _spatialcluster(self,Target_sum=1e4,N_top_genes=3000,Rad_cutoff=50,Resolution=0.8,Expression_alpha=0.5):
        sc.pp.normalize_total(self.adata, target_sum=Target_sum)
        sc.pp.log1p(self.adata)
        sc.pp.highly_variable_genes(self.adata,n_top_genes=N_top_genes)
        STAGATE.Cal_Spatial_Net(self.adata, rad_cutoff=Rad_cutoff)
        STAGATE.Stats_Spatial_Net(adata)
        self.adata = STAGATE.train_STAGATE(self.adata, alpha=Expression_alpha)
        sc.pp.neighbors(self.adata, use_rep='SCRIPro')
        sc.tl.umap(self.adata)
        sc.tl.leiden(self.adata, resolution=Resolution)
        sc.pl.embedding(self.adata, basis="spatial", color=['leiden'],s=6, show=False, title='SCRIPro',save="Spatial_cluster.pdf")
        sc.pl.umap(self.adata, color=['leiden'],show=False,save="Spatial_umap_leiden.pdf")
    
    def _fliter(self,Target_sum=1e4,N_top_genes=3000,Max_value=10,N_neighbors=10,N_pcs=40,Resolution=0.6):
        sc.pp.normalize_total(self.adata, target_sum=Target_sum)
        sc.pp.log1p(self.adata)
        sc.pp.highly_variable_genes(self.adata,n_top_genes=N_top_genes)
        self.adata.raw = self.adata
        self.adata = self.adata[:, self.adata.var.highly_variable]
        #sc.pp.regress_out(self.adata, ['total_counts', 'pct_counts_mt'])
        sc.pp.scale(self.adata, max_value=Max_value)
        sc.tl.pca(self.adata,svd_solver='arpack')
        sc.pp.neighbors(self.adata, n_neighbors=N_neighbors, n_pcs=N_pcs)
        sc.tl.umap(self.adata)
        sc.tl.leiden(self.adata,resolution=Resolution)
 
    def get_positive_marker_gene_parallel(self, gene_list_len=30,cores=2):
        pool = mp.Pool(cores)
        groups = set(self.adata.obs['new_leiden'])
        func = partial(get_marker_for_group, self.adata, gene_list_len=gene_list_len)
        results = pool.map(func, groups)
        pool.close()
        self.marker_list = {i: tem_gene for i, tem_gene in results if tem_gene is not None}
        
    def get_largedata_markergene(self,topgene=1500,quantile=0.6):
        sc.tl.rank_genes_groups(self.adata,self.cluster_method, method='t-test')
        all_marker_gene_list = {}
        for i in set(self.adata.obs[self.cluster_method]):
            all_marker_gene_list[i] = set(self.adata.uns['rank_genes_groups']['names'][i][0:topgene])
        marker_gene_list = {}
        candidate_gene = self.ad_all[self.ad_all.apply(lambda x: x > x.quantile(0.6),axis=0)]
        for j in candidate_gene.index:
            leiden_cluster = j[0]
            df = candidate_gene.filter(items=all_marker_gene_list[leiden_cluster])
            row_df = df.loc[j,:]
            meta_marker = row_df[~(np.isnan(row_df))].index
            if len(meta_marker) < 35:
                continue
            else:
                marker_gene_list[j] = meta_marker
        self.marker_list = marker_gene_list
        
    def _supercell(self,cluster_method='leiden'):
        adList = []
        adList_obs = {}
        for i in set(self.adata.obs[cluster_method]):
            leiden_index = self.adata.obs.loc[self.adata.obs[cluster_method] == i].index
            sub_cluster = self.adata.raw.to_adata()[leiden_index]
            merged_data,obs = supercell_pipeline(sub_cluster,cell_num=self.Cell_num,verbose=False)
            merged_data_index = [i + "_" + str(j) for j in range(0,merged_data.shape[0])]
            merged_data.obs = pd.DataFrame(index = merged_data_index)
            adList.append(merged_data)
            adList_obs[i] = obs
        ad_all = ad.concat(adList,join='outer')
        sc.pp.normalize_total(ad_all, target_sum=1e4)
        sc.pp.log1p(ad_all)
        ad_all = ad_all.to_df()
        self.adata.obs['new_leiden'] = np.zeros(len(self.adata.obs))
        for i in adList_obs.keys():
            for k in adList_obs[i].index:
                self.adata.obs.loc[k,'new_leiden'] = i+"_"+adList_obs[i].loc[k,'leiden']  
        self.ad_all = ad_all
        
    def get_supercell_exp(self,pvalue_matrix):
        super_exp_tra = self.ad_all.loc[pvalue_matrix.index,pvalue_matrix.columns]
        super_exp_tra = self._replace_all_outliers(super_exp_tra)
        super_exp_tra = super_exp_tra - np.min(super_exp_tra)/np.max(super_exp_tra) - np.min(super_exp_tra)
        super_exp_tra = super_exp_tra.replace(np.nan,0)
        self.super_exp_tra = super_exp_tra       

class SCRIPro_Multiome():
    def __init__(self,cores,species,Ori_Data):
        self.cores = cores
        self.Ori_Data = Ori_Data
        self.species = species
        '''
    def cal_ISD_parallel(self, bw_path):
        lisa_test = FromGenes(self.species, rp_map='enhanced_10K', assays=['Direct'], isd_method='chipseq', verbose=0)
        datainterface = lisa_test.data_interface
        rpmap_enhanced = datainterface.get_rp_map('enhanced_10K')
        factor_binging, factor_dataset_ids, factor_metadata = datainterface.get_binding_data(technology='ChIP-seq')
        final_dict = {}
        lisa_info = {}
        for i in self.Ori_Data.marker_list.keys():
            query_list = self.Ori_Data.marker_list[i]
            query_genes, background_genes = lisa_test._get_query_and_background_genes(query_list)
            gene_mask, label_vector, gene_info_dict = lisa_test._make_gene_mask(query_genes, background_genes)
            lisa_info[i] = [gene_mask, label_vector, gene_info_dict, self.species]
        results = process_map(process_marker, [(i, lisa_info, bw_path, rpmap_enhanced, factor_binging, factor_metadata, self.species) for i in lisa_info.keys()], max_workers=self.cores)
        for i, factor_metadata_pd in results:
            final_dict[i] = factor_metadata_pd
        self.results = final_dict
        '''
    '''
    def cal_ISD_parallel(self, bw_path):
        lisa_test = FromGenes(self.species, rp_map='enhanced_10K', assays=['Direct'], isd_method='chipseq', verbose=0)
        datainterface = lisa_test.data_interface
        rpmap_enhanced = datainterface.get_rp_map('enhanced_10K')
        factor_binging, factor_dataset_ids, factor_metadata = datainterface.get_binding_data(technology='ChIP-seq')
        final_dict = {}
        lisa_info = {}
        for i in self.Ori_Data.marker_list.keys():
            query_list = self.Ori_Data.marker_list[i]
            query_genes, background_genes = lisa_test._get_query_and_background_genes(query_list)
            gene_mask, label_vector, gene_info_dict = lisa_test._make_gene_mask(query_genes, background_genes)
            lisa_info[i] = [gene_mask, label_vector, gene_info_dict, self.species]
        #tasks = [(i, lisa_info[i], bw_path, rpmap_enhanced, factor_binging, factor_metadata, self.species) for i in lisa_info.keys()]
        tasks = [(i,) + lisa_info[i] + (bw_path, rpmap_enhanced, factor_binging, factor_metadata, self.species) for i in lisa_info.keys()]
        results = process_map(process_marker, tasks, max_workers=self.cores)
        for i, factor_metadata_pd in results:
            final_dict[i] = factor_metadata_pd
        self.results = final_dict

    def cal_ISD_parallel(self, bw_path):
        lisa_test = FromGenes(self.species, rp_map='enhanced_10K', assays=['Direct'], isd_method='chipseq', verbose=0)
        datainterface = lisa_test.data_interface
        rpmap_enhanced = datainterface.get_rp_map('enhanced_10K')
        factor_binging, factor_dataset_ids, factor_metadata = datainterface.get_binding_data(technology='ChIP-seq')
        final_dict = {}
        lisa_info = {}
        for i in self.Ori_Data.marker_list.keys():
            query_list = self.Ori_Data.marker_list[i]
            query_genes, background_genes = lisa_test._get_query_and_background_genes(query_list)
            gene_mask, label_vector, gene_info_dict = lisa_test._make_gene_mask(query_genes, background_genes)
            lisa_info[i] = [gene_mask, label_vector, gene_info_dict, self.species]

        tasks = [(i, lisa_info, bw_path, rpmap_enhanced, factor_binging, factor_metadata, self.species) for i in lisa_info.keys()]
        # 使用 Pool.starmap 在多进程中运行 process_marker 函数
        with Pool(processes=self.cores) as pool:
            # 这里我们使用 tqdm 来创建进度条
            results = []
            for result in tqdm(pool.starmap(process_marker, tasks), total=len(tasks)):
                results.append(result)
        # 收集结果
        for i, factor_metadata_pd in results:
            final_dict[i] = factor_metadata_pd

        self.results = final_dict
'''     
    def cal_ISD_parallel(self, bw_path):
        lisa_test = FromGenes(self.species, rp_map='enhanced_10K', assays=['Direct'], isd_method='chipseq', verbose=0)
        datainterface = lisa_test.data_interface
        rpmap_enhanced = datainterface.get_rp_map('enhanced_10K')
        factor_binging, factor_dataset_ids, factor_metadata = datainterface.get_binding_data(technology='ChIP-seq')
        final_dict = {}
        lisa_info = {}
        for i in self.Ori_Data.marker_list.keys():
            query_list = self.Ori_Data.marker_list[i]
            query_genes, background_genes = lisa_test._get_query_and_background_genes(query_list)
            gene_mask, label_vector, gene_info_dict = lisa_test._make_gene_mask(query_genes, background_genes)
            lisa_info[i] = [gene_mask, label_vector, gene_info_dict, self.species]

        # 使用ProcessPoolExecutor替换原有的Pool
        with ProcessPoolExecutor(max_workers=self.cores) as executor:
            # 准备任务列表
            tasks = {executor.submit(process_marker, i, lisa_info, bw_path, rpmap_enhanced, factor_binging, factor_metadata, self.species): i for i in lisa_info}

            # 使用tqdm创建进度条
            for future in tqdm(as_completed(tasks), total=len(tasks), desc="Processing markers"):
                i, factor_metadata_pd = future.result()
                final_dict[tasks[future]] = factor_metadata_pd

        self.results = final_dict
    
    def get_P_value_matrix(self):
        results_combine = []
        for i in self.results.keys():
            test = self.results[i].sort_values(by = 'p_vals').drop_duplicates(subset='factor', keep='first').loc[:,['factor','p_vals']]
            test[i] = -np.log(test['p_vals'])
            test.index = test['factor']
            test = test.loc[:,i]
            results_combine.append(test)
        all_result_combine = pd.concat([i for i in results_combine],axis=1)
        all_result_combine = (all_result_combine - all_result_combine.mean())/(all_result_combine.std())
        all_result_combine.clip(upper=5, inplace=True)
        all_result_combine = (all_result_combine - all_result_combine.min())/(all_result_combine.max() - all_result_combine.min())
        all_result_combine = all_result_combine.T
        P_value_matrix = all_result_combine
        self.P_value_matrix = P_value_matrix
    
    def get_chip_matrix(self):
        results_combine = []
        for i in self.results.keys():
            test = self.results[i].sort_values(by = 'p_vals').drop_duplicates(subset='factor', keep='first').loc[:,['sample_id','factor']]
            test.index = test['factor']
            test=test.iloc[:,0]
            test.name = i
            results_combine.append(test)
        all_result_combine = pd.concat([i for i in results_combine],axis=1)
        all_result_combine = all_result_combine.T
        chip_matrix = all_result_combine
        self.chip_matrix = chip_matrix
    
    
    def get_tf(self):
        target_h5_path=pkg_resources.resource_filename('scripro', 'data/TF_target_RP.h5')
        target_h5 = pd.HDFStore(target_h5_path)
        super_gene_exp = self.Ori_Data.super_gene_exp
        super_gene_mean = self.Ori_Data.super_gene_mean
        super_gene_std = self.Ori_Data.super_gene_std
        chip_matrix2 = self.chip_matrix.copy()
        chip_matrix2 = chip_matrix2.loc[:,list(set(chip_matrix2.columns).intersection(set(super_gene_exp.columns)))]
        
        querygene_dict = get_quert_gene(chip_matrix2,target_h5,super_gene_exp.columns)
        chip_matrix2 = chip_matrix2.applymap(lambda x: querygene_dict[x])
        
        chip_matrix2=chip_matrix2.reset_index().rename(columns={'index': 'row'})
        chip_matrix2_melt = chip_matrix2.melt(id_vars='row', var_name='column', value_name='value')
        self.chip_matrix_melt = chip_matrix2_melt
        def cal_tf_score(row):
            tf_name = row['column']
            supercell_name = row['row']
            query_gene = row['value']
            return ((super_gene_exp.loc[supercell_name, query_gene] - super_gene_mean[query_gene])/super_gene_std[query_gene]).mean() + (super_gene_exp.loc[supercell_name, tf_name] - super_gene_mean[tf_name])/ super_gene_std[tf_name]
        
        
        chip_matrix2_melt['value2'] = chip_matrix2_melt.apply(cal_tf_score, axis=1)
        chip_matrix2_new = chip_matrix2_melt.pivot(index='row', columns='column', values='value2')
        chip_matrix2_new = chip_matrix2_new.reset_index()
        chip_matrix2_new= chip_matrix2_new.rename_axis(None, axis=1)
        chip_matrix2_new = chip_matrix2_new.set_index('row')
        
        chip_matrix2_new.clip(upper=4,lower=-4,inplace=True)
        chip_matrix2_new = (chip_matrix2_new - chip_matrix2_new.min())/(chip_matrix2_new.max() - chip_matrix2_new.min())
        chip_matrix2_new = chip_matrix2_new.replace(np.nan,0)
        
        self.P_value_matrix = self.P_value_matrix.loc[chip_matrix2_new.index,chip_matrix2_new.columns]
        test_mat = chip_matrix2_new.mul(self.P_value_matrix)
        
        
        self.tf_score = test_mat
    
    def get_tf_only_target(self):
        target_h5=pkg_resources.resource_filename('scripro', 'data/TF_target_RP.h5')
        super_gene_exp = self.Ori_Data.super_gene_exp
        super_gene_mean = self.Ori_Data.super_gene_mean
        super_gene_std = self.Ori_Data.super_gene_std
        chip_matrix2 = self.chip_matrix.copy()
        chip_matrix2 = chip_matrix2.loc[:,list(set(chip_matrix2.columns).intersection(set(super_gene_exp.columns)))]
        
        querygene_dict = get_quert_gene(chip_matrix2,target_h5,super_gene_exp.columns)
        chip_matrix2 = chip_matrix2.applymap(lambda x: querygene_dict[x])
        
        chip_matrix2=chip_matrix2.reset_index().rename(columns={'index': 'row'})
        chip_matrix2_melt = chip_matrix2.melt(id_vars='row', var_name='column', value_name='value')
        self.chip_matrix_melt = chip_matrix2_melt
        def cal_tf_score(row):
            tf_name = row['column']
            supercell_name = row['row']
            query_gene = row['value']
            return ((super_gene_exp.loc[supercell_name, query_gene] - super_gene_mean[query_gene])/super_gene_std[query_gene]).mean()
        
        
        chip_matrix2_melt['value2'] = chip_matrix2_melt.apply(cal_tf_score, axis=1)
        chip_matrix2_new = chip_matrix2_melt.pivot(index='row', columns='column', values='value2')
        chip_matrix2_new = chip_matrix2_new.reset_index()
        chip_matrix2_new= chip_matrix2_new.rename_axis(None, axis=1)
        chip_matrix2_new = chip_matrix2_new.set_index('row')
        
        chip_matrix2_new.clip(upper=4,lower=-4,inplace=True)
        chip_matrix2_new = (chip_matrix2_new - chip_matrix2_new.min())/(chip_matrix2_new.max() - chip_matrix2_new.min())
        chip_matrix2_new = chip_matrix2_new.replace(np.nan,0)
        
        self.P_value_matrix = self.P_value_matrix.loc[chip_matrix2_new.index,chip_matrix2_new.columns]
        self.chip_matrix2_new = chip_matrix2_new
        test_mat = self.chip_matrix2_new.mul(self.P_value_matrix)
        
        
        self.tf_score = test_mat
    
    def get_tf_target(self,TF):
        C =self.chip_matrix_melt[self.chip_matrix_melt.column == TF]
        
        merged = pd.merge(self.Ori_Data.adata.obs, self.Ori_Data.adata.raw.to_adata().to_df(), left_index=True, right_index=True)
        results = []
        for index, row in C.iterrows():
            cluster = row['row']
            target_genes = row['value']
            cluster_data = merged[merged['new_leiden'] == cluster]
            tf_results = pd.DataFrame(index=[cluster], columns=target_genes)
            for gene in target_genes:
                tf_expression = cluster_data[TF]
                gene_expression = cluster_data[gene]
                corr, _ = pearsonr(tf_expression, gene_expression)
                tf_results.loc[cluster, gene] = corr
            results.append(tf_results*self.tf_score.loc[cluster,TF])
        return pd.concat(results).clip(lower=0).replace(np.nan,0)
        
class SCRIPro_RNA(SCRIPro_Multiome):  
    def __init__(self, cores, species, Ori_Data, assays):
        super().__init__(cores, species, Ori_Data)
        self.assays = assays
    def _chunks(self,data):
        SIZE=math.ceil(len(data)/self.cores)
        it = iter(data)
        for i in range(0, len(data), SIZE):
            yield {k:data[k] for k in islice(it, SIZE)}
            
    def cal_ISD_cistrome(self):
        input_list = self.Ori_Data.marker_list
        cal = partial(cal_tf, species1=self.species, assays1=self.assays)
        pbar = tqdm(total=len(input_list))
        def update_pbar(*args):
            pbar.update()
        all_tf_result = {}
        with ProcessPoolExecutor(max_workers=self.cores) as executor:
            futures = [executor.submit(cal, input_data) for input_data in input_list.items()]
            for future in futures:
                future.add_done_callback(update_pbar)
                result = future.result()
                all_tf_result.update(result)
        pbar.close()
        self.results = all_tf_result
        
        
    def get_P_value_matrix(self):
        results_combine = []
        for i in self.results.keys():
            test = self.results[i].sort_values(by = 'summary_p_value').drop_duplicates(subset='factor', keep='first').loc[:,['factor','summary_p_value']]
            test[i] = -np.log(test['summary_p_value'])
            test.index = test['factor']
            test = test.loc[:,i]
            results_combine.append(test)
        all_result_combine = pd.concat([i for i in results_combine],axis=1)
        all_result_combine = (all_result_combine - all_result_combine.mean())/(all_result_combine.std())
        all_result_combine.clip(upper=5, inplace=True)
        all_result_combine = (all_result_combine - all_result_combine.min())/(all_result_combine.max() - all_result_combine.min())
        all_result_combine = all_result_combine.T
        P_value_matrix = all_result_combine
        self.P_value_matrix = P_value_matrix 
    def get_chip_matrix(self):
        results_combine = []
        for i in self.results.keys():
            test = self.results[i].sort_values(by = 'summary_p_value').drop_duplicates(subset='factor', keep='first').loc[:,['sample_id','factor']]
            test.index = test['factor']
            test=test.iloc[:,0]
            test.name = i
            results_combine.append(test)
        all_result_combine = pd.concat([i for i in results_combine],axis=1)
        all_result_combine = all_result_combine.T
        chip_matrix = all_result_combine
        self.chip_matrix = chip_matrix
    
        

          

