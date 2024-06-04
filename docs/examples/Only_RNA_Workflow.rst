RNA Only Workflow
-------------------------------


This is the most standard SCRIPro computation process, requiring only
the input of the corresponding scRNA-seq sequencing matrix.
To demonstrate SCRIP's ability to be applied to different tissue types and infer target genes for TRs, we applied SCRIP to 10X lymphoma sequencing data. Data are available on https://www.10xgenomics.com/datasets/fresh-frozen-lymph-node-with-b-cell-lymphoma-14-k-sorted-nuclei-1-standard-2-0


The resulting tf_score matrix can be obtained by using the following shell statement:

.. code:: ipython3

    scripro enrich -i ./data/rna/rna.h5ad -n 50 -s hs -p rna_workflow -t 32


 ========================
   
The resulting gata3_score matrix can be obtained through the following shell statement, where rna_workflow.pkl is the result of SCRIPro enrich:


.. code:: ipython3

    scripro get_tf_target -i rna_workflow.pkl -t GATA3 -p GATA3_target


 ========================


.. code:: ipython3

    import numpy as np
    import pandas as pd
    import scanpy as sc
    import h5py
    import scripro
    import anndata
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns
    import pandas as pd
    import scanpy as sc
    import warnings
    warnings.filterwarnings("ignore")

Load and preprocess data
========================

.. code:: ipython3

    rna = sc.read_h5ad('./data/rna/rna.h5ad')
    rna.var_names_make_unique()
    rna.raw = rna
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
    sc.pp.highly_variable_genes(rna, min_mean=0.0125, max_mean=3, min_disp=0.5)
    rna = rna[:, rna.var.highly_variable]
    sc.pp.scale(rna, max_value=10)
    sc.tl.pca(rna, svd_solver='arpack')
    sc.pp.neighbors(rna)
    sc.tl.umap(rna)
    sc.tl.leiden(rna)

Calculate the Supercell and the marker genes corresponding to Supercell

.. code:: ipython3

    test_data = scripro.Ori_Data(rna,Cell_num=50)

ad_all is the integrated counting matrix

.. code:: ipython3

    test_data.ad_all




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>MIR1302-2HG</th>
          <th>FAM138A</th>
          <th>OR4F5</th>
          <th>AL627309.1</th>
          <th>AL627309.3</th>
          <th>AL627309.2</th>
          <th>AL627309.5</th>
          <th>AL627309.4</th>
          <th>AP006222.2</th>
          <th>AL732372.1</th>
          <th>...</th>
          <th>AC133551.1</th>
          <th>AC136612.1</th>
          <th>AC136616.1</th>
          <th>AC136616.3</th>
          <th>AC136616.2</th>
          <th>AC141272.1</th>
          <th>AC023491.2</th>
          <th>AC007325.1</th>
          <th>AC007325.4</th>
          <th>AC007325.2</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>20_0</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>15_0</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>15_1</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>15_2</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>13_0</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>9_4</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>9_5</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>9_6</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>9_7</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>21_0</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
        </tr>
      </tbody>
    </table>
    <p>224 rows × 36621 columns</p>
    </div>



.. code:: ipython3

    test_data.get_positive_marker_gene_parallel()


.. parsed-literal::

    Compute marker gene


.. code:: ipython3

    rna_seq_data = scripro.SCRIPro_RNA(5,'hg38',test_data,assays=['Direct','DNase','H3K27ac'])

The computational process of In Silico Deletion
===============================================

.. code:: ipython3

  
    rna_seq_data.cal_ISD_cistrome()


.. parsed-literal::

    


The P-value matrix of each Supercell LISA is obtained according to the
calculation results

Get TF activity Score
=====================

.. code:: ipython3

    rna_seq_data.get_tf_score()

.. code:: ipython3

    rna_seq_data.P_value_matrix




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>ADNP</th>
          <th>AFF1</th>
          <th>AFF4</th>
          <th>AGO1</th>
          <th>AHR</th>
          <th>AIRE</th>
          <th>ALX1</th>
          <th>ALX3</th>
          <th>ALX4</th>
          <th>ANHX</th>
          <th>...</th>
          <th>ZSCAN22</th>
          <th>ZSCAN23</th>
          <th>ZSCAN29</th>
          <th>ZSCAN30</th>
          <th>ZSCAN31</th>
          <th>ZSCAN4</th>
          <th>ZSCAN5A</th>
          <th>ZSCAN5C</th>
          <th>ZXDB</th>
          <th>ZXDC</th>
        </tr>
        <tr>
          <th>row</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0_0</th>
          <td>1.982159e-05</td>
          <td>0.114342</td>
          <td>0.466165</td>
          <td>3.044442e-03</td>
          <td>0.065143</td>
          <td>0.116164</td>
          <td>0.261117</td>
          <td>0.090598</td>
          <td>0.043649</td>
          <td>0.070920</td>
          <td>...</td>
          <td>0.001946</td>
          <td>1.034024e-03</td>
          <td>0.000837</td>
          <td>0.023628</td>
          <td>0.187771</td>
          <td>0.130556</td>
          <td>0.000345</td>
          <td>0.072917</td>
          <td>9.929228e-07</td>
          <td>1.078112e-06</td>
        </tr>
        <tr>
          <th>0_1</th>
          <td>1.078489e-03</td>
          <td>0.045135</td>
          <td>0.541748</td>
          <td>4.741197e-02</td>
          <td>0.172083</td>
          <td>0.137448</td>
          <td>0.120097</td>
          <td>0.091863</td>
          <td>0.078125</td>
          <td>0.097334</td>
          <td>...</td>
          <td>0.027452</td>
          <td>6.524492e-02</td>
          <td>0.119130</td>
          <td>0.071906</td>
          <td>0.200513</td>
          <td>0.117636</td>
          <td>0.007210</td>
          <td>0.072906</td>
          <td>1.114402e-05</td>
          <td>3.193426e-03</td>
        </tr>
        <tr>
          <th>0_10</th>
          <td>1.945398e-04</td>
          <td>0.150389</td>
          <td>0.350183</td>
          <td>7.688059e-02</td>
          <td>0.089623</td>
          <td>0.316572</td>
          <td>0.277354</td>
          <td>0.399970</td>
          <td>0.437044</td>
          <td>0.195209</td>
          <td>...</td>
          <td>0.021498</td>
          <td>1.736244e-03</td>
          <td>0.091324</td>
          <td>0.003618</td>
          <td>0.320272</td>
          <td>0.071882</td>
          <td>0.000904</td>
          <td>0.098806</td>
          <td>2.213682e-06</td>
          <td>1.677967e-02</td>
        </tr>
        <tr>
          <th>0_11</th>
          <td>9.016532e-02</td>
          <td>0.124475</td>
          <td>0.635978</td>
          <td>2.211520e-02</td>
          <td>0.178290</td>
          <td>0.010232</td>
          <td>0.077026</td>
          <td>0.126848</td>
          <td>0.065793</td>
          <td>0.001066</td>
          <td>...</td>
          <td>0.211864</td>
          <td>4.717477e-02</td>
          <td>0.126473</td>
          <td>0.111667</td>
          <td>0.130438</td>
          <td>0.169036</td>
          <td>0.055158</td>
          <td>0.244485</td>
          <td>4.748398e-04</td>
          <td>1.358551e-02</td>
        </tr>
        <tr>
          <th>0_12</th>
          <td>1.508612e-01</td>
          <td>0.220131</td>
          <td>0.714978</td>
          <td>1.149924e-01</td>
          <td>0.166783</td>
          <td>0.000201</td>
          <td>0.019816</td>
          <td>0.003010</td>
          <td>0.003320</td>
          <td>0.003520</td>
          <td>...</td>
          <td>0.349635</td>
          <td>1.420289e-01</td>
          <td>0.171647</td>
          <td>0.123673</td>
          <td>0.080900</td>
          <td>0.042576</td>
          <td>0.047124</td>
          <td>0.017884</td>
          <td>1.611482e-01</td>
          <td>2.017362e-01</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>9_3</th>
          <td>1.481955e-05</td>
          <td>0.161472</td>
          <td>0.000004</td>
          <td>6.475927e-07</td>
          <td>0.004738</td>
          <td>0.093825</td>
          <td>0.145126</td>
          <td>0.158836</td>
          <td>0.204868</td>
          <td>0.006100</td>
          <td>...</td>
          <td>0.000030</td>
          <td>6.431066e-08</td>
          <td>0.041991</td>
          <td>0.001208</td>
          <td>0.000560</td>
          <td>0.012364</td>
          <td>0.000022</td>
          <td>0.036678</td>
          <td>5.952748e-08</td>
          <td>2.198499e-08</td>
        </tr>
        <tr>
          <th>9_4</th>
          <td>1.624109e-07</td>
          <td>0.304159</td>
          <td>0.185860</td>
          <td>1.608332e-02</td>
          <td>0.018612</td>
          <td>0.205191</td>
          <td>0.173053</td>
          <td>0.138393</td>
          <td>0.167866</td>
          <td>0.051846</td>
          <td>...</td>
          <td>0.006800</td>
          <td>1.012524e-04</td>
          <td>0.031388</td>
          <td>0.001566</td>
          <td>0.097648</td>
          <td>0.044065</td>
          <td>0.000073</td>
          <td>0.019923</td>
          <td>1.451613e-03</td>
          <td>7.308369e-03</td>
        </tr>
        <tr>
          <th>9_5</th>
          <td>1.541161e-06</td>
          <td>0.252129</td>
          <td>0.000368</td>
          <td>4.775720e-04</td>
          <td>0.036822</td>
          <td>0.136602</td>
          <td>0.147106</td>
          <td>0.204738</td>
          <td>0.165820</td>
          <td>0.031218</td>
          <td>...</td>
          <td>0.015975</td>
          <td>1.854799e-03</td>
          <td>0.069004</td>
          <td>0.008719</td>
          <td>0.092146</td>
          <td>0.088071</td>
          <td>0.000901</td>
          <td>0.005200</td>
          <td>1.631952e-04</td>
          <td>3.722424e-05</td>
        </tr>
        <tr>
          <th>9_6</th>
          <td>6.143819e-05</td>
          <td>0.349253</td>
          <td>0.150809</td>
          <td>3.164199e-02</td>
          <td>0.089277</td>
          <td>0.122468</td>
          <td>0.182552</td>
          <td>0.158537</td>
          <td>0.181882</td>
          <td>0.090961</td>
          <td>...</td>
          <td>0.012562</td>
          <td>5.747627e-03</td>
          <td>0.085607</td>
          <td>0.011577</td>
          <td>0.090943</td>
          <td>0.081455</td>
          <td>0.004634</td>
          <td>0.016923</td>
          <td>3.773492e-03</td>
          <td>5.942802e-02</td>
        </tr>
        <tr>
          <th>9_7</th>
          <td>6.450485e-04</td>
          <td>0.390047</td>
          <td>0.199128</td>
          <td>1.675784e-02</td>
          <td>0.132506</td>
          <td>0.096528</td>
          <td>0.102888</td>
          <td>0.107414</td>
          <td>0.135996</td>
          <td>0.100875</td>
          <td>...</td>
          <td>0.016645</td>
          <td>9.027264e-03</td>
          <td>0.067132</td>
          <td>0.021804</td>
          <td>0.122074</td>
          <td>0.053077</td>
          <td>0.000223</td>
          <td>0.008073</td>
          <td>8.117502e-03</td>
          <td>7.536773e-03</td>
        </tr>
      </tbody>
    </table>
    <p>224 rows × 1226 columns</p>
    </div>



The corresponding RP score and expression value are used to weight the
P-value obtained, and the final tf activity score is obtained

.. code:: ipython3

    rna_seq_data.tf_score




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>ADNP</th>
          <th>AFF1</th>
          <th>AFF4</th>
          <th>AGO1</th>
          <th>AHR</th>
          <th>AIRE</th>
          <th>ALX1</th>
          <th>ALX3</th>
          <th>ALX4</th>
          <th>ANHX</th>
          <th>...</th>
          <th>ZSCAN22</th>
          <th>ZSCAN23</th>
          <th>ZSCAN29</th>
          <th>ZSCAN30</th>
          <th>ZSCAN31</th>
          <th>ZSCAN4</th>
          <th>ZSCAN5A</th>
          <th>ZSCAN5C</th>
          <th>ZXDB</th>
          <th>ZXDC</th>
        </tr>
        <tr>
          <th>row</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0_0</th>
          <td>1.181346e-05</td>
          <td>0.060435</td>
          <td>0.307493</td>
          <td>1.462677e-04</td>
          <td>0.026594</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000038</td>
          <td>0.001076</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>8.489554e-06</td>
          <td>0.0</td>
          <td>4.385504e-08</td>
          <td>6.365249e-07</td>
        </tr>
        <tr>
          <th>0_1</th>
          <td>8.153228e-04</td>
          <td>0.028895</td>
          <td>0.455855</td>
          <td>2.507957e-03</td>
          <td>0.008484</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.006283</td>
          <td>0.005021</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>3.823604e-03</td>
          <td>0.0</td>
          <td>5.890852e-07</td>
          <td>1.917075e-03</td>
        </tr>
        <tr>
          <th>0_10</th>
          <td>1.138860e-04</td>
          <td>0.095834</td>
          <td>0.293383</td>
          <td>3.976904e-02</td>
          <td>0.037968</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.001849</td>
          <td>0.000498</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>6.564053e-05</td>
          <td>0.0</td>
          <td>1.429377e-07</td>
          <td>9.996831e-03</td>
        </tr>
        <tr>
          <th>0_11</th>
          <td>6.903511e-02</td>
          <td>0.076661</td>
          <td>0.422427</td>
          <td>1.190686e-03</td>
          <td>0.011600</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.006117</td>
          <td>0.008547</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>4.512259e-03</td>
          <td>0.0</td>
          <td>3.173963e-05</td>
          <td>8.298006e-03</td>
        </tr>
        <tr>
          <th>0_12</th>
          <td>8.898146e-02</td>
          <td>0.136908</td>
          <td>0.467959</td>
          <td>5.825133e-02</td>
          <td>0.009677</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.008715</td>
          <td>0.010824</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>1.991192e-03</td>
          <td>0.0</td>
          <td>9.687363e-03</td>
          <td>1.161664e-01</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>9_3</th>
          <td>8.815053e-06</td>
          <td>0.060870</td>
          <td>0.000002</td>
          <td>4.206756e-08</td>
          <td>0.001850</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000759</td>
          <td>0.000037</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>6.775830e-07</td>
          <td>0.0</td>
          <td>2.239746e-09</td>
          <td>1.277747e-08</td>
        </tr>
        <tr>
          <th>9_4</th>
          <td>1.216592e-07</td>
          <td>0.160054</td>
          <td>0.121961</td>
          <td>9.350271e-04</td>
          <td>0.007390</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.030882</td>
          <td>0.000199</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>2.069416e-06</td>
          <td>0.0</td>
          <td>6.039253e-05</td>
          <td>9.517900e-05</td>
        </tr>
        <tr>
          <th>9_5</th>
          <td>1.182939e-06</td>
          <td>0.095240</td>
          <td>0.000301</td>
          <td>2.557181e-05</td>
          <td>0.002016</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.002766</td>
          <td>0.000738</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>4.742670e-04</td>
          <td>0.0</td>
          <td>6.393928e-06</td>
          <td>2.261003e-05</td>
        </tr>
        <tr>
          <th>9_6</th>
          <td>3.662379e-05</td>
          <td>0.186604</td>
          <td>0.140389</td>
          <td>2.558394e-03</td>
          <td>0.038996</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.083027</td>
          <td>0.000672</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>3.029751e-04</td>
          <td>0.0</td>
          <td>1.696993e-04</td>
          <td>3.468848e-02</td>
        </tr>
        <tr>
          <th>9_7</th>
          <td>4.759191e-04</td>
          <td>0.254694</td>
          <td>0.129037</td>
          <td>9.390380e-04</td>
          <td>0.009117</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.063992</td>
          <td>0.002329</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>1.695396e-05</td>
          <td>0.0</td>
          <td>8.699077e-04</td>
          <td>4.172868e-03</td>
        </tr>
      </tbody>
    </table>
    <p>224 rows × 1226 columns</p>
    </div>



Calculate the downstream target gene of each TF in each Supercell
=================================================================

.. code:: ipython3

    gata3_score = rna_seq_data.get_tf_target('GATA3')

.. code:: ipython3

    gata3_score




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>SOS1</th>
          <th>ZNF487</th>
          <th>PPP1CA</th>
          <th>CFLAR</th>
          <th>WDR37</th>
          <th>CTLA4</th>
          <th>STK10</th>
          <th>NFKBIL1</th>
          <th>INO80B</th>
          <th>PPP2R5C</th>
          <th>...</th>
          <th>BCL2</th>
          <th>RPL18</th>
          <th>PRSS55</th>
          <th>UBL4B</th>
          <th>FAM13A</th>
          <th>WDR20</th>
          <th>SYTL3</th>
          <th>ASH1L</th>
          <th>APOC3</th>
          <th>CPNE8</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>3_10</th>
          <td>0.012644</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.096325</td>
          <td>0.000000</td>
          <td>0.026573</td>
          <td>0.000000</td>
          <td>0.067059</td>
          <td>0.021823</td>
          <td>...</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>4_1</th>
          <td>0.239298</td>
          <td>0.025236</td>
          <td>0.000000</td>
          <td>0.111141</td>
          <td>0.000000</td>
          <td>0.133851</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.077034</td>
          <td>...</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>1_0</th>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>...</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>12_4</th>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.120566</td>
          <td>0.209552</td>
          <td>0.093906</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>...</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>10_4</th>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>...</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>22_0</th>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>...</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>1_17</th>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>...</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>0_24</th>
          <td>0.095861</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.051342</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.012070</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>...</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>0_3</th>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.027181</td>
          <td>0.220751</td>
          <td>0.000000</td>
          <td>0.030307</td>
          <td>...</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>3_2</th>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.009417</td>
          <td>0.045787</td>
          <td>0.115108</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>...</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
        </tr>
      </tbody>
    </table>
    <p>224 rows × 3084 columns</p>
    </div>



