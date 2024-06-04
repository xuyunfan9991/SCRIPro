Usage
============

SCRIPro is dedicated for single-cell or spatial RNA-seq datasets and multiome datasets (chromatin accessibility and transcriptome) to explore the gene regulation relationships. The package includes four main functions.  


.. code:: 

    usage: scripro [-h] [--version] {enrich_rna,enrich_multiome,get_tf_target,install_reference,atac} ...

    scripro

    positional arguments:
    {enrich_rna,enrich_multiome,get_tf_target,install_reference,atac}
        enrich_rna          Calculate TF activation use scRNA-seq data
        enrich_multiome     Calculate TF activation use scRNA-seq data and scATAC-seq data
        get_tf_target       Calculate TF and target gene score
        install_reference   Install RP reference
        atac                Invoke SCRIP


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
- This function will output a pkl including the pvalue matrix and tf activity score matrix.


scripro enrich_multiome
~~~~~~~~~~~~~~~~~~~~~~~~

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

- This function will output a pkl including the pvalue matrix and tf activity score matrix and a dir named 'bigwig' will also be generated, containing the each supercell's corresponding chromatin landscape file


scripro get_tf_target
~~~~~~~~~~~~~~~~~~~~~~

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
- This function will output a csv containing the regulatory activity of tf downstream target genes within each supercell is generated


scripro atac
~~~~~~~~~~~~~~~~~~

For enrichment of TF activity for single-cell ATAC-seq data, you can use ``scripro atac`` command:
``scripro atac`` act same as ``SCRIP``, the function includs ``enrich``, ``impute``, ``target``, ``config``, and ``index``.

The reference files for SCRIP are different from SCRIPro, which you can download from `zenodo <https://zenodo.org/record/5840810>`_ and config with ``SCRIP config``.  

Using example:

.. code::

    scripro atac enrich [-h] -i FEATURE_MATRIX -s {hs,mm} [-p PROJECT] [--min_cells MIN_CELLS] [--min_peaks MIN_PEAKS] [--max_peaks MAX_PEAKS] [-t N_CORES] [-m {max,mean}] [-y] [--clean]

    optional arguments:
    -h, --help            show this help message and exit

    Input files arguments:
    -i FEATURE_MATRIX, --input_feature_matrix FEATURE_MATRIX
                            A cell by peak matrix . REQUIRED.
    -s {hs,mm}, --species {hs,mm}
                            Species. "hs"(human) or "mm"(mouse). REQUIRED.

    Output arguments:
    -p PROJECT, --project PROJECT
                            Project name, which will be used to generate output files folder. DEFAULT: Random generate.

    Preprocessing paramater arguments:
    --min_cells MIN_CELLS
                            Minimal cell cutoff for features. Auto will take 0.05% of total cell number.DEFAULT: "auto".
    --min_peaks MIN_PEAKS
                            Minimal peak cutoff for cells. Auto will take the mean-3*std of all feature number (if less than 500 is 500). DEFAULT: "auto".
    --max_peaks MAX_PEAKS
                            Max peak cutoff for cells. This will help you to remove the doublet cells. Auto will take the mean+5*std of all feature
                            number. DEFAULT: "auto".

    Other options:
    -t N_CORES, --thread N_CORES
                            Number of cores use to run SCRIP. DEFAULT: 16.
    -m {max,mean}, --mode {max,mean}
                            Deduplicate strategy. DEFAULT: max.
    -y, --yes             Whether ask for confirmation. DEFAULT: False.
    --clean               Whether delete tmp files(including bed and search results) generated by SCRIP. DEFAULT: False.

- In this function, you can input a peak count matrix in H5 or MTX format, with basic parameters of quality control. 
- This function will output a folder including these files:

    + beds: bed files of all cells
    + ChIP_result: txt files of Giggle search results
    + peaks_length.txt: peak total length of each cell
    + SCRIP_enrichment.txt: the result of the SCRIP score
    + dataset_overlap_df.pk: the raw number of overlaps of each cell to each dataset
    + dataset_cell_norm_df.pk: normalized scores
    + dataset_score_source_df.pk: matched reference datasets
    + tf_cell_score_df.pk: the same table to SCRIP_enrichment.txt but untransposed and in pickle format

- detail usage see `SCRIP documentation <https://scrip.readthedocs.io/en/latest/usage.html>`_
