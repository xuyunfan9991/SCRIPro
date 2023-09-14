.. SCRIPro documentation master file, created by
   sphinx-quickstart on Tue Aug  3 19:13:29 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SCRIPro's documentation!
==================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

|Docs|

.. |Docs| image:: https://readthedocs.org/projects/scripro/badge/?version=latest
   :target: https://scripro.readthedocs.io

The increasing availability of single-cell genomics data allows characterizing gene regulation networks (GRNs) at an unprecedented resolution. Previously we developed `SCRIP <https://scrip.readthedocs.io/en/latest/>`_, a computational method that integrates single-cell ATAC-seq data with a large-scale transcription regulator (TR) ChIP-seq data and motif reference for reconstructing single-cell TR activities and GRNs. 

Here, we present SCRIPro, an extended framework of SCRIP that suits both single-cell and spatial multi-ome data. SCRIPro first performed a density clustering based on expression and spatial similarity of the data to generate high coverage SuperCells. Next, SCRIPro performed in silico deletion analyses based on matched scATAC-seq or reconstructed chromatin landscapes from Cistrome accessibility data to evaluate the importance of TRs in regulating each SuperCell. Finally, SCRIPro combines the importance score of each TR with its gene expression using GRNBOOST2, which generates the TR-centered GRNs at Supercell resolution. We applied SCRIPro on human PBMC and human B-cell lymphoma scMulti-ome data, as well as mouse developing embryo spatial transcriptomic data, and demonstrated that SCRIPro is able to identify cell-type-specific TR regulations and show superior performance compared to conventional motif-based methods such as SCENIC+. Taken together, SCRIPro is a convenient and fast method that could reconstruct TR activities and GRNs for both single-cell and spatial multi-omic data.

.. image:: _static/img/workflow_new.png
   :alt: Workflow of SCRIPro
   :width: 100%
   :align: center

.. toctree::
   :maxdepth: 1
   :hidden:

   installation
   usage
   examples
   release_notes
