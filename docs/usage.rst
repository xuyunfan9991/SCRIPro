Usage
============

SCRIPro includes three main functions.  

For enrichment of TF activity for single-cell or spatial RNA-seq data, you can use:

.. code:: 

    scripro enrich_rna [-h] -i FEATURE_MATRIX -n CELL_NUM -s {hg38,mm10} -p PROJECT [-t N_CORES]

    optional arguments:
    -h, --help            show this help message and exit

    Input files arguments:
    -i FEATURE_MATRIX, --input_feature_matrix FEATURE_MATRIX
                            scRNA-seq data matrix . REQUIRED.
    -n CELL_NUM, --cell_number CELL_NUM
                            Supercell Cell Number . REQUIRED.
    -s {hg38,mm10}, --species {hg38,mm10}
                            Species. "hg38"(human) or "mm10"(mouse). REQUIRED.

    Output arguments:
    -p PROJECT, --project PROJECT
                            Project name, which will be used to generate output files.

    Other options:
    -t N_CORES, --thread N_CORES
                            Number of cores use to run SCRIPro. DEFAULT: 8.

 For enrichment of TF activity for both RNA-seq and ATAC-seq for single-cell or spatial data, you can use:

.. code:: 

    scripro enrich_multiome [-h] -i FEATURE_MATRIX -n CELL_NUM -s {hg38,mm10} -a {fragment,matrix} -b {0,1} -f ATAC_PATH [-g GLUE_ANNOTATION] -p
                                PROJECT [-t N_CORES]

    optional arguments:
    -h, --help            show this help message and exit

    Input files arguments:
    -i FEATURE_MATRIX, --input_feature_matrix FEATURE_MATRIX
                            A cell by peak matrix . REQUIRED.
    -n CELL_NUM, --cell_number CELL_NUM
                            Supercell Cell Number . REQUIRED.
    -s {hg38,mm10}, --species {hg38,mm10}
                            Species. "hg38"(human) or "mm10"(mouse). REQUIRED.
    -a {fragment,matrix}, --atac_file_type {fragment,matrix}
                            atac_file_type,"fragment" or "matrix(h5,h5ad,mtx)"
    -b {0,1}, --barcode_corresponds {0,1}
                            Whether the scRNA-seq barcode matches the scATAC-seq barcode. "0"(Match) or "1"(Not match). REQUIRED.
    -f ATAC_PATH, --atac_file ATAC_PATH
                            ATAC file path.REQUIRED.
    -g GLUE_ANNOTATION, --glue_annotation GLUE_ANNOTATION
                            If the scRNA-seq barcodes \do not match the scATAC-seq barcodes, the glue_annotation file that will be used.,like
                            'gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz'

    Output arguments:
    -p PROJECT, --project PROJECT
                            Project name, which will be used to generate output files.

    Other options:
    -t N_CORES, --thread N_CORES
                            Number of cores use to run SCRIPros. DEFAULT: 8.

For getting the target of specific TR, you can use:

.. code:: 

    scripro get_tf_target [-h] -i SCRIPRO_RESULT -t TF_NAME -p PROJECT

    optional arguments:
    -h, --help            show this help message and exit

    Input files arguments:
    -i SCRIPRO_RESULT, --input_scripro_result SCRIPRO_RESULT
                            scripro result pickle file. REQUIRED.
    -t TF_NAME, --tf_name TF_NAME
                            Tf name to calculate the target . REQUIRED.

    Output arguments:
    -p PROJECT, --project PROJECT
                            Project name, which will be used to generate output file.
