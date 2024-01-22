import argparse
import pickle
from .supercell import *
from .utils import *
from .Ori_data import *
import sys
from numpy import require
import warnings
warnings.filterwarnings("ignore")

def run_enrich_only_rna(args):
    feature_matrix_path = args.feature_matrix
    num = args.cell_num
    species = args.species
    project = args.project
    n_cores = args.n_cores
    filename=args.project
    print("Data preprocessing...")
    if feature_matrix_path.endswith('.h5'):
        feature_matrix = sc.read_10x_h5(feature_matrix_path, gex_only=False)
    elif feature_matrix_path.endswith('.h5ad'):
        feature_matrix = sc.read_h5ad(feature_matrix_path)
    else:
        feature_matrix = sc.read_10x_mtx(feature_matrix_path)
    feature_matrix.var_names_make_unique()
    sc.pp.normalize_total(feature_matrix, target_sum=1e4)
    sc.pp.log1p(feature_matrix)
    sc.pp.highly_variable_genes(feature_matrix,n_top_genes=3000)
    feature_matrix.raw = feature_matrixq
    feature_matrix = feature_matrix[:, feature_matrix.var.highly_variable]
    #sc.pp.regress_out(self.adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(feature_matrix, max_value=10)
    sc.tl.pca(feature_matrix,svd_solver='arpack')
    sc.pp.neighbors(feature_matrix, n_neighbors=10, n_pcs=40)
    sc.tl.umap(feature_matrix)
    sc.tl.leiden(feature_matrix,resolution=0.6)
    print("Calculating supercell...")
    test_data = Ori_Data(feature_matrix,Cell_num=num)
    print("Calculating markergene...")
    test_data.get_positive_marker_gene_parallel(cores=n_cores)
    rna_seq_data = SCRIPro_RNA(n_cores,species,test_data,assays=['Direct','DNase','H3K27ac'])
    print("Executing In Silico Deletion")
    rna_seq_data.cal_ISD_cistrome()
    rna_seq_data.get_P_value_matrix()
    rna_seq_data.get_chip_matrix()
    rna_seq_data.P_value_matrix
    rna_seq_data.get_tf()
    rna_seq_data.tf_score.to_csv(filename+'.csv', index=False, sep=',')
    with open('scripro_result.pkl', 'wb') as file:
        pickle.dump(rna_seq_data, file)
    return 


def run_enrich_multiome(args):
    feature_matrix_path = args.feature_matrix
    num = args.cell_num
    species = args.species
    project = args.project
    n_cores = args.n_cores
    atac_path = args.atac_path
    filename=args.project
    barcodes = args.barcodes
    atac_file_type=args.atac_file_type
    glue_annotation = args.glue_annotation
    
    
    
    if feature_matrix_path.endswith('.h5'):
        feature_matrix = sc.read_10x_h5(feature_matrix_path, gex_only=False)
    elif feature_matrix_path.endswith('.h5ad'):
        feature_matrix = sc.read_h5ad(feature_matrix_path)
    else:
        feature_matrix = sc.read_10x_mtx(feature_matrix_path)
    feature_matrix.var_names_make_unique()
    
    
    
    if atac_path.endswith('.tsv'):
        pass
    elif feature_matrix_path.endswith('.h5'):
        atac_matrix = sc.read_10x_h5(atac_path, gex_only=False)
        atac_matrix.var_names_make_unique()
    elif feature_matrix_path.endswith('.h5ad'):
        atac_matrix = sc.read_h5ad(atac_path)
        atac_matrix.var_names_make_unique()
    else:
        atac_matrix = sc.read_10x_mtx(atac_path)
        atac_matrix.var_names_make_unique()
    print("Data preprocessing...")
    feature_matrix.layers["counts"] = feature_matrix.X.copy()
    sc.pp.normalize_total(feature_matrix, target_sum=1e4)
    sc.pp.log1p(feature_matrix)
    sc.pp.highly_variable_genes(feature_matrix,n_top_genes=3000)
    feature_matrix.raw = feature_matrix
    feature_matrix = feature_matrix[:, feature_matrix.var.highly_variable]
    #sc.pp.regress_out(self.adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(feature_matrix, max_value=10)
    sc.tl.pca(feature_matrix,svd_solver='arpack')
    sc.pp.neighbors(feature_matrix, n_neighbors=10, n_pcs=40)
    sc.tl.umap(feature_matrix)
    sc.tl.leiden(feature_matrix,resolution=0.6)
    
    print("Calculating supercell...")
    test_data = Ori_Data(feature_matrix,Cell_num=num)
    
    if barcodes == '0':
        print("Calculating Chromatin landscape")
        cellgroup = test_data.adata.obs.loc[:,['new_leiden']]
        test_data.get_positive_marker_gene_parallel(cores=n_cores)
        if atac_file_type=='fragment':
            get_supercell_fragment(cellgroup,'.',atac_path,chunksize = 10000000)
            process_tsv('./supercell_fragment/', species)
        else: 
            dataframe_to_sparse_tsv(atac.to_df(), 'test.tsv')
            get_supercell_fragment(cellgroup,'.','test.tsv',chunksize = 10000000)
            process_tsv('./supercell_fragment/', species)
    else:
        print("USE GLue to calculate barcode")
        scglue.data.lsi(atac_matrix, n_components=5, n_iter=15)
        sc.pp.neighbors(atac_matrix, use_rep="X_lsi", metric="cosine")
        sc.tl.umap(atac_matrix)
        scglue.data.get_gene_annotation(feature_matrix, gtf=glue_annotation,gtf_by="gene_name")
        genes_to_remove = feature_matrix.var[~(feature_matrix.var.loc[:,"chromStart"]>0)].index
        feature_matrix = feature_matrix[:, ~feature_matrix.var.index.isin(genes_to_remove)]
        split = atac_matrix.var_names.str.split(r"[_]")
        atac_matrix.var["chrom"] = split.map(lambda x: x[0])
        atac_matrix.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
        atac_matrix.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)
        guidance = scglue.genomics.rna_anchored_guidance_graph(feature_matrix, atac_matrix)
        scglue.graph.check_graph(guidance, [feature_matrix, atac_matrix])
        scglue.models.configure_dataset(feature_matrix, "NB", use_highy_variable=True,use_layer="counts", use_rep="X_pca")
        scglue.models.configure_dataset(atac_matrix, "NB", use_highly_variable=True,use_rep="X_lsi")
        guidance_hvf = guidance.subgraph(chain(feature_matrix.var.query("highly_variable").index,atac_matrix.var.query("highly_variable").index)).copy()
        glue = scglue.models.fit_SCGLUE({"rna": feature_matrix, "atac": atac_matrix}, guidance_hvf,fit_kws={"directory": "glue"})
        dx = scglue.models.integration_consistency(glue, {"rna": feature_matrix, "atac": atac_matrix}, guidance_hvf)
        _ = sns.lineplot(x="n_meta", y="consistency", data=dx).axhline(y=0.05, c="darkred", ls="--")
        feature_matrix.obsm["X_glue"] = glue.encode_data("rna", feature_matrix)
        atac_matrix.obsm["X_glue"] = glue.encode_data("atac", atac_matrix)
        feature_matrix.obs['feature']='rna'
        atac_matrix.obs['feature']='atac'
        combined = ad.concat([feature_matrix, atac_matrixs])
        sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
        sc.tl.umap(combined)
        sc.pl.umap(combined)
        sc.tl.leiden(combined,resolution=0.8)
        combined_rna = combined[combined.obs.feature == 'rna']
        combined_rna.obs.loc[:,'new_leiden'] = np.nan
        glue_supercell(combined_rna,50)
        rna_leiden_clusters = combined_rna.obs['new_leiden']
        combined_atac = combined[combined.obs.feature == 'atac']
        distance_matrix = cdist(combined_atac.obsm['X_umap'], combined_rna.obsm['X_umap'], metric='euclidean')
        nearest_rna = np.argmin(distance_matrix, axis=1)
        atac_leiden_clusters = rna_leiden_clusters[nearest_rna]
        atac_leiden_clusters.index = combined_atac.obs.index
        feature_matrix.obs = combined_rna.obs
        cellgroup = pd.DataFrame(atac_leiden_clusters)
        if atac_file_type=='fragment':
            test_data.get_glue_cluster(rna_leiden_clusters)
            test_data.get_positive_marker_gene_parallel(cores=n_cores)
            get_supercell_fragment(cellgroup,'.',atac_path,chunksize = 10000000)
            process_tsv('./supercell_fragment/', 'hg38')
        else:
            dataframe_to_sparse_tsv(atac.to_df(), 'test.tsv')
            get_supercell_fragment(cellgroup,'.','./test.tsv',chunksize = 10000000)
            process_tsv('./supercell_fragment/', 'hg38')
    print("Calculating ISD")
    share_seq_data = SCRIPro_Multiome(8,'hg38',test_data)
    share_seq_data.cal_ISD_parallel('./bigwig/')
    share_seq_data.get_P_value_matrix()
    share_seq_data.get_chip_matrix()
    share_seq_data.get_tf()
    share_seq_data.tf_score.to_csv(filename+'.csv', index=False, sep=',')
    with open('scripro_result.pkl', 'wb') as file:
        pickle.dump(share_seq_data, file)
    return 


def get_target_score(args):
    pickle_file = args.scripro_result
    tf = args.tf_name
    project = args.project
    print("Data preprocessing...")
    with open(pickle_file, 'rb') as file:
        scripro_data = pickle.load(file)
    target = scripro_data.get_tf_target(tf)
    target.to_csv(project+'_target.csv', index=False, sep=',')
    return 



def main():
    warnings.filterwarnings("ignore")
    argparser = prepare_argparser()
    args = argparser.parse_args()
    subcommand  = args.subcommand
    if subcommand == "enrich_rna":
        try:
            run_enrich_only_rna(args)
        except MemoryError:
            sys.exit( "MemoryError occurred.")
    elif subcommand == "enrich_multiome":
        try:
            run_enrich_multiome(args)
        except MemoryError:
            sys.exit( "MemoryError occurred.")
    elif subcommand == "get_tf_target":
        try:
            get_target_score(args)
        except MemoryError:
            sys.exit( "MemoryError occurred.")

    return

def prepare_argparser():
    description = "%(prog)s"
    epilog = "For command line options of each command, type: %(prog)s COMMAND -h"
    argparser = argparse.ArgumentParser( description = description, epilog = epilog )
    argparser.add_argument( "--version", action="version", version="0.1.6")
    subparsers = argparser.add_subparsers( dest = 'subcommand' )
    subparsers.required = True
    add_enrich_parser(subparsers)
    add_enrich_parser_multiome(subparsers)
    add_target_parser(subparsers)
    return argparser


def add_enrich_parser( subparsers ):
    """Add main function 'enrich' argument parsers.
    """
    argparser_enrich_rna = subparsers.add_parser("enrich_rna", help="Calculate TF activation use scRNA-seq data")

    # group for input files
    group_input = argparser_enrich_rna.add_argument_group( "Input files arguments" )
    group_input.add_argument( "-i", "--input_feature_matrix", dest = "feature_matrix", type = str, required = True,
                              help = 'scRNA-seq data matrix . REQUIRED.' )
    group_input.add_argument( "-n", "--cell_number", dest = "cell_num", type = int, required = True,
                              help = 'Supercell Cell Number . REQUIRED.' )
    group_input.add_argument( "-s", "--species", dest = "species", choices= ['hg38', 'mm10'], required = True,
                              help = 'Species. "hg38"(human) or "mm10"(mouse). REQUIRED.' )
    # group for output files
    group_output = argparser_enrich_rna.add_argument_group( "Output arguments" )
    group_output.add_argument( "-p", "--project", dest = "project", type = str, default = "" ,required = True,
                               help = 'Project name, which will be used to generate output files.')
    
    # group for preprocessing
    group_preprocessing = argparser_enrich_rna.add_argument_group( "Preprocessing paramater arguments" )
    group_other = argparser_enrich_rna.add_argument_group( "Other options" )
    group_other.add_argument( "-t", '--thread', dest='n_cores', type = int, default = 8,
                              help="Number of cores use to run SCRIP. DEFAULT: 8.")

    
def add_enrich_parser_multiome( subparsers ):
    """Add main function 'enrich' argument parsers.
    """
    argparser_enrich_multiome = subparsers.add_parser("enrich_multiome", help="Calculate TF activation use scRNA-seq data and scATAC-seq data")

    # group for input files
    group_input = argparser_enrich_multiome.add_argument_group( "Input files arguments" )
    group_input.add_argument( "-i", "--input_feature_matrix", dest = "feature_matrix", type = str, required = True,
                              help = 'A cell by peak matrix . REQUIRED.' )
    group_input.add_argument( "-n", "--cell_number", dest = "cell_num", type = int, required = True,
                              help = 'Supercell Cell Number . REQUIRED.' )
    group_input.add_argument( "-s", "--species", dest = "species", choices= ['hg38', 'mm10'], required = True,
                              help = 'Species. "hg38"(human) or "mm10"(mouse). REQUIRED.' )
    group_input.add_argument( "-a", "--atac_file_type", dest = "atac_file_type", choices= ['fragment', 'matrix'], required = True,
                              help = 'atac_file_type,"fragment" or "matrix(h5,h5ad,mtx)"' )
    group_input.add_argument( "-b", "--barcode_corresponds", dest = "barcodes", choices= ['0', '1'], required = True,
                              help = 'Whether the rna barcode and atac barcode correspond to each other. "0"(Corresponding) or "1"(Not corresponding). REQUIRED.' )
    group_input.add_argument( "-f", "--atac_file", dest = "atac_path", type = str, required = True,
                              help = 'ATAC file,"fragment" or "matrix(h5,h5ad,mtx)". REQUIRED.' )
    group_input.add_argument( "-g", "--glue_annotation", dest = "glue_annotation", type = str,
                              help = 'glue_annotation,like \'gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz\'' )
    # group for output files
    group_output = argparser_enrich_multiome.add_argument_group( "Output arguments" )
    group_output.add_argument( "-p", "--project", dest = "project", type = str, default = "" ,required = True,
                               help = 'Project name, which will be used to generate output files. DEFAULT: Random generate.')
    # group for preprocessing
    group_preprocessing = argparser_enrich_multiome.add_argument_group( "Preprocessing paramater arguments" )
    group_other = argparser_enrich_multiome.add_argument_group( "Other options" )
    group_other.add_argument( "-t", '--thread', dest='n_cores', type = int, default = 8,
                              help="Number of cores use to run SCRIPros. DEFAULT: 8.")
    
def add_target_parser( subparsers ):
    """Add main function 'enrich' argument parsers.
    """
    argparser_target = subparsers.add_parser("get_tf_target", help="Calculate TF and target gene score")

    # group for input files
    group_input = argparser_target.add_argument_group( "Input files arguments" )
    group_input.add_argument( "-i", "--input_scripro_result", dest = "scripro_result", type = str, required = True,
                              help = 'scripro result pickle file. REQUIRED.' )
    group_input.add_argument( "-t", "--tf_name", dest = "tf_name", type = str, required = True,
                              help = 'Tf name to calculate the target . REQUIRED.' )
    
    # group for output files
    group_output = argparser_target.add_argument_group( "Output arguments" )
    group_output.add_argument( "-p", "--project", dest = "project", type = str, default = "" ,required = True,
                               help = 'Project name, which will be used to generate output file.')

if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    try:
        warnings.filterwarnings("ignore")
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted!\n")
        sys.exit(0)