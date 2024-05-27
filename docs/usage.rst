Usage
============

SCRIPro is dedicated for single-cell or spatial RNA-seq datasets and multiome datasets (chromatin accessibility and transcriptome) to explore the gene regulation relationships. The package includes four main functions.  


.. code:: 

    usage: scripro [-h] [--version] {enrich_rna,enrich_multiome,get_tf_target,install_reference} ...

    scripro

    positional arguments:
    {enrich_rna,enrich_multiome,get_tf_target,install_reference}
        enrich_rna          Calculate TF activation use scRNA-seq data
        enrich_multiome     Calculate TF activation use scRNA-seq data and scATAC-seq data
        get_tf_target       Calculate TF and target gene score
        install_reference   Install RP reference

    optional arguments:
    -h, --help            show this help message and exit
    --version             show program's version number and exit

    For command line options of each command, type: scripro COMMAND -h

Detailed usages are listed as follows:

scripro enrich_rna
~~~~~~~~~~~~~~~~~~

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

- In this function, you can input the feature count matrix in H5 or MTX format.   
- The ``-n`` parameter controls how many cells merge into a supercell. This parameter affects the resolution of the results.  
- This function will output a folder including these files:  
  + ``xxx``: xxxx.
  + ``xxx``: xxxx.


scripro enrich_multiome
~~~~~~~~~~~~~~~~~~

For enrichment of TF activity for both RNA-seq and ATAC-seq for single-cell or spatial data, you can use:

.. code:: 

    usage: scripro enrich_multiome [-h] -i FEATURE_MATRIX -n CELL_NUM -s {hg38,mm10} -a {fragment,matrix} -b {0,1} -f ATAC_PATH [-g GLUE_ANNOTATION] -p PROJECT [-t N_CORES]

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
                            If the scRNA-seq barcodes do not match the scATAC-seq barcodes, the glue_annotation file that will be used.,like 'gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz'

    Output arguments:
    -p PROJECT, --project PROJECT
                            Project name, which will be used to generate output files.

    Other options:
    -t N_CORES, --thread N_CORES
                            Number of cores use to run SCRIPros. DEFAULT: 8.

- In this function, you are allowed to input a transcriptome dataset and a chromatin accessibility dataset.  
- For transcriptome data, a feature count matrix is required. For chromatin accessibility, you can input a fragment file or feature count matrix either.   
- The ``-n`` parameter controls how many cells merge into a supercell. This parameter affects the resolution of the results.   
- For barcode matched multiome dataset, like SHARE-seq or 10X multiome dataset, the ``-b`` should be set to ``0``. Otherwise, this should be set as ``1``.  
- If ``-b`` is set as ``1``, a GTF annotation file need to provide.   

- This function will output a folder including these files:  
  + ``xxx``: xxxx.
  + ``xxx``: xxxx.


scripro get_tf_target
~~~~~~~~~~~~~~~~~~

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

- In this function, you can input the results of ``enrich_rna`` or ``enrich_multiome`` and a TF name and will output the target genes of the TF.  
- This function will output a folder including these files:  
  + ``xxx``: xxxx.
  + ``xxx``: xxxx.
