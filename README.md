# SCRIPro

[![Documentation Status](https://readthedocs.org/projects/scripro/badge/?version=latest)](https://scripro.readthedocs.io/en/latest/?badge=latest)


The increasing availability of single-cell genomics data allows characterizing gene regulation networks (GRNs) at an unprecedented resolution. Previously we developed SCRIP, a computational method that integrates single-cell ATAC-seq data with a large-scale transcription regulator (TR) ChIP-seq data and motif reference for reconstructing single-cell TR activities and GRNs. Here, we present SCRIPro, an extended framework of SCRIP that suits both single-cell and spatial multi-ome data. SCRIPro first performed a density clustering based on expression and spatial similarity of the data to generate high coverage SuperCells. Next, SCRIPro performed in silico deletion analyses based on matched scATAC-seq or reconstructed chromatin landscapes from Cistrome accessibility data to evaluate the importance of TRs in regulating each SuperCell. Finally, SCRIPro combines the importance score of each TR with its gene expression, which generates the TR-centered GRNs at Supercell resolution. We applied SCRIPro on human PBMC and human B-cell lymphoma scMulti-ome data, as well as mouse developing embryo spatial transcriptomic data, and demonstrated that SCRIPro is able to identify cell-type-specific TR regulations and show superior performance compared to conventional motif-based methods such as SCENIC+. Taken together, SCRIPro is a convenient and fast method that could reconstruct TR activities and GRNs for both single-cell and spatial multi-omic data.

![avatar](docs/_static/img/workflow_new.png)

## Documentation

For the detailed usage and examples of SCRIPro, please refer to the [documentation](https://scripro.readthedocs.io/en/latest).  
For any problems encountered in using, feel free to open an [issue](https://github.com/xuyunfan9991/SCRIPro/issues).  


## Installation
### Use the following commands to install Minicoda3：
``` bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
It is recommended to install scripro to a conda virtual environment.:
~~~~
conda create -n scripro python=3.8
conda install -c liulab-dfci lisa2
lisa install hg38 oneshot ./hg38_scripro.h5 --force
lisa install mm10 oneshot ./mm10_scripro.h5 --force
~~~~~
hg38_scripro.h5 and mm10_scripro.h5 can be downloaded from https://zenodo.org/records/10554173.
### Install SCRIPro from Github (developing version)
```bash
git clone git@github.com:xuyunfan9991/SCRIPro.git
cd SCRIPro
python setup.py install
```

### Install SCRIPro from Pypi
~~~~
pip install scripro
~~~~~

## Usage
``` bash
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

```
``` bash
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

```
``` bash
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
```

