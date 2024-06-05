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



Install STAGATE
~~~~~~~~~~~~~~~~~~~~~


`STAGATE <https://github.com/QIFEIDKN/STAGATE>`_ is designed for spatial clustering and denoising expressions of spatial resolved transcriptomics (ST) data.

To install STAGATE, use:

.. code:: shell

   git clone git@github.com:QIFEIDKN/STAGATE.git
   cd STAGATE
   python setup.py build
   python setup.py install


