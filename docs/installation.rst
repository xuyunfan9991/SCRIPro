Installation
==============


System requirements
~~~~~~~~~~~~~~~~~~~

* Linux/Unix
* Python >= 3.8

We recommend to create an independent conda environment for SCRIPro. If users do not have conda, please install Miniconda first:

.. code:: shell

   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh

Installation
~~~~~~~~~~~~~~~~~~~~~

.. code:: shell

   conda create -n scripro python=3.8
   conda activate scripro
   conda install -c liulab-dfci lisa2
   pip install scripro

Next, please download the reference files of SCRIPro from `Zenodo <https://zenodo.org/doi/10.5281/zenodo.10554172>`_ and configure it:  

.. code:: shell

   lisa install hg38 oneshot hg38_scripro.h5
   lisa install mm10 oneshot mm10_scripro.h5
   scripro install_reference -i TF_target_RP.h5

Install SCRIP
~~~~~~~~~~~~~~~~~~~~~

For scATAC-seq only datasets, we recommand to use `SCRIP <https://scrip.readthedocs.io/en/latest/>`_.  

SCRIP requires a separate conda environment, and we recommend that you create a new conda environment to install SCRIP. You can refer to this `webpage <https://scrip.readthedocs.io/en/latest/installation.html>`_ for the SCRIP installation instructions. Below is a simple workflow:  

.. code:: shell

   conda create -n scrip python=3.8
   conda activate scrip
   pip install scrip
   SCRIP install_giggle

The reference files for SCRIP are different from SCRIPro, which you can download from `zenodo <https://zenodo.org/record/5840810>`_ and config with ``SCRIP config``.  
