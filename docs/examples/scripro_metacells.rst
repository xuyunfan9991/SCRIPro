Metacell & Supercell Comparison
----------------------------------

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
    from sklearn.metrics import roc_curve, auc

Load and process data
=====================

Load data and raw data. Since the processed data downloaded does not
contain raw data, we download raw data for re-processing

.. code:: ipython3

    #rna = sc.read_h5ad('/fs/home/xuyunfan/data/ORF/select_tf.h5ad')
    rna = sc.read_h5ad('/fs/home/xuyunfan/Final/review/data/perturb_data.h5ad')
    rna_raw = sc.read_h5ad('/fs/home/xuyunfan/Final/review/data/raw.h5ad')
    rna_raw=rna_raw[rna.obs.index]
    sc.pp.normalize_total(rna_raw, target_sum=1e4)
    sc.pp.log1p(rna_raw)
    rna.raw = rna_raw
    test_data = scripro.Ori_Data(rna,Cell_num=50)
    rna.obs




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
          <th>gene</th>
          <th>n_genes</th>
          <th>n_genes_by_counts</th>
          <th>total_counts</th>
          <th>total_counts_mt</th>
          <th>pct_counts_mt</th>
          <th>leiden</th>
          <th>mixscape_class_p_ko</th>
          <th>mixscape_class</th>
          <th>mixscape_class_global</th>
          <th>pertclass</th>
          <th>hdbscan</th>
        </tr>
        <tr>
          <th>Cell_barcodes</th>
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
          <th>TAACCAGAGTAGAATC-8</th>
          <td>TRIM21</td>
          <td>3467</td>
          <td>3467</td>
          <td>10422.0</td>
          <td>755.0</td>
          <td>7.244291</td>
          <td>26</td>
          <td>1.0</td>
          <td>TRIM21 KO</td>
          <td>KO</td>
          <td>strong</td>
          <td>9</td>
        </tr>
        <tr>
          <th>CATAGACCAACACGAG-8</th>
          <td>CBY1</td>
          <td>2003</td>
          <td>2003</td>
          <td>4621.0</td>
          <td>392.0</td>
          <td>8.483012</td>
          <td>24</td>
          <td>1.0</td>
          <td>CBY1 KO</td>
          <td>KO</td>
          <td>strong</td>
          <td>10</td>
        </tr>
        <tr>
          <th>CTGTGAATCCGGTAAT-2</th>
          <td>LAT2</td>
          <td>4344</td>
          <td>4344</td>
          <td>16784.0</td>
          <td>1412.0</td>
          <td>8.412774</td>
          <td>9</td>
          <td>1.0</td>
          <td>LAT2 KO</td>
          <td>KO</td>
          <td>strong</td>
          <td>1</td>
        </tr>
        <tr>
          <th>GAGCTGCAGGTAGATT-8</th>
          <td>RELA</td>
          <td>2361</td>
          <td>2360</td>
          <td>6086.0</td>
          <td>380.0</td>
          <td>6.243838</td>
          <td>16</td>
          <td>1.0</td>
          <td>RELA KO</td>
          <td>KO</td>
          <td>strong</td>
          <td>4</td>
        </tr>
        <tr>
          <th>AAGTACCCAACTTCTT-3</th>
          <td>WT1</td>
          <td>2198</td>
          <td>2198</td>
          <td>5469.0</td>
          <td>545.0</td>
          <td>9.965259</td>
          <td>12</td>
          <td>1.0</td>
          <td>WT1 KO</td>
          <td>KO</td>
          <td>strong</td>
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
        </tr>
        <tr>
          <th>GGCTGTGAGGGCTAAC-5</th>
          <td>APOL2</td>
          <td>2503</td>
          <td>2503</td>
          <td>6126.0</td>
          <td>553.0</td>
          <td>9.027098</td>
          <td>30</td>
          <td>1.0</td>
          <td>APOL2 KO</td>
          <td>KO</td>
          <td>strong</td>
          <td>20</td>
        </tr>
        <tr>
          <th>ATATCCTCATCATTTC-8</th>
          <td>TNFRSF1B</td>
          <td>4380</td>
          <td>4379</td>
          <td>14271.0</td>
          <td>1052.0</td>
          <td>7.371593</td>
          <td>11</td>
          <td>1.0</td>
          <td>TNFRSF1B KO</td>
          <td>KO</td>
          <td>strong</td>
          <td>7</td>
        </tr>
        <tr>
          <th>CTAGGTAGTTGAGGAC-1</th>
          <td>CD27</td>
          <td>2385</td>
          <td>2385</td>
          <td>6731.0</td>
          <td>336.0</td>
          <td>4.991829</td>
          <td>19</td>
          <td>1.0</td>
          <td>CD27 KO</td>
          <td>KO</td>
          <td>strong</td>
          <td>17</td>
        </tr>
        <tr>
          <th>TGGGAAGGTGAGTTTC-6</th>
          <td>CTRL</td>
          <td>2988</td>
          <td>2988</td>
          <td>7500.0</td>
          <td>541.0</td>
          <td>7.213333</td>
          <td>3</td>
          <td>0.0</td>
          <td>CTRL</td>
          <td>CTRL</td>
          <td>CTRL</td>
          <td>30</td>
        </tr>
        <tr>
          <th>ATGATCGGTATCGTTG-7</th>
          <td>RELA</td>
          <td>2739</td>
          <td>2739</td>
          <td>7606.0</td>
          <td>311.0</td>
          <td>4.088877</td>
          <td>16</td>
          <td>1.0</td>
          <td>RELA KO</td>
          <td>KO</td>
          <td>strong</td>
          <td>4</td>
        </tr>
      </tbody>
    </table>
    <p>16707 rows × 12 columns</p>
    </div>



replace supercell data with metacell data
=========================================

Load the metacell data calculated by metacell, then replace supercell
data with metacell data (new_leiden column)

.. code:: ipython3

    metacell = pd.read_csv('./metacells.csv')
    metacell




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
          <th>Cell_barcodes</th>
          <th>gene</th>
          <th>excluded_umis</th>
          <th>metacell</th>
          <th>dissolved</th>
          <th>metacell_level</th>
          <th>cells_rare_gene_module</th>
          <th>rare_cell</th>
          <th>metacell_name</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>TAACCAGAGTAGAATC-8</td>
          <td>TRIM21</td>
          <td>849.0</td>
          <td>110</td>
          <td>False</td>
          <td>1</td>
          <td>-1</td>
          <td>False</td>
          <td>M110.11</td>
        </tr>
        <tr>
          <th>1</th>
          <td>CATAGACCAACACGAG-8</td>
          <td>CBY1</td>
          <td>459.0</td>
          <td>91</td>
          <td>False</td>
          <td>1</td>
          <td>-1</td>
          <td>False</td>
          <td>M91.30</td>
        </tr>
        <tr>
          <th>2</th>
          <td>CTGTGAATCCGGTAAT-2</td>
          <td>LAT2</td>
          <td>1571.0</td>
          <td>484</td>
          <td>False</td>
          <td>2</td>
          <td>-1</td>
          <td>False</td>
          <td>M484.07</td>
        </tr>
        <tr>
          <th>3</th>
          <td>GAGCTGCAGGTAGATT-8</td>
          <td>RELA</td>
          <td>408.0</td>
          <td>131</td>
          <td>False</td>
          <td>1</td>
          <td>-1</td>
          <td>False</td>
          <td>M131.58</td>
        </tr>
        <tr>
          <th>4</th>
          <td>AAGTACCCAACTTCTT-3</td>
          <td>WT1</td>
          <td>779.0</td>
          <td>412</td>
          <td>False</td>
          <td>1</td>
          <td>-1</td>
          <td>False</td>
          <td>M412.96</td>
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
        </tr>
        <tr>
          <th>16702</th>
          <td>GGCTGTGAGGGCTAAC-5</td>
          <td>APOL2</td>
          <td>712.0</td>
          <td>219</td>
          <td>False</td>
          <td>1</td>
          <td>-1</td>
          <td>False</td>
          <td>M219.04</td>
        </tr>
        <tr>
          <th>16703</th>
          <td>ATATCCTCATCATTTC-8</td>
          <td>TNFRSF1B</td>
          <td>1107.0</td>
          <td>5</td>
          <td>False</td>
          <td>1</td>
          <td>-1</td>
          <td>False</td>
          <td>M5.31</td>
        </tr>
        <tr>
          <th>16704</th>
          <td>CTAGGTAGTTGAGGAC-1</td>
          <td>CD27</td>
          <td>445.0</td>
          <td>235</td>
          <td>False</td>
          <td>1</td>
          <td>-1</td>
          <td>False</td>
          <td>M235.05</td>
        </tr>
        <tr>
          <th>16705</th>
          <td>TGGGAAGGTGAGTTTC-6</td>
          <td>CTRL</td>
          <td>650.0</td>
          <td>100</td>
          <td>False</td>
          <td>1</td>
          <td>-1</td>
          <td>False</td>
          <td>M100.73</td>
        </tr>
        <tr>
          <th>16706</th>
          <td>ATGATCGGTATCGTTG-7</td>
          <td>RELA</td>
          <td>332.0</td>
          <td>470</td>
          <td>False</td>
          <td>2</td>
          <td>-1</td>
          <td>False</td>
          <td>M470.82</td>
        </tr>
      </tbody>
    </table>
    <p>16707 rows × 9 columns</p>
    </div>



.. code:: ipython3

    test_data.adata.obs.new_leiden=list(metacell.metacell)

.. code:: ipython3

    test_data.adata.obs




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
          <th>gene</th>
          <th>n_genes</th>
          <th>n_genes_by_counts</th>
          <th>total_counts</th>
          <th>total_counts_mt</th>
          <th>pct_counts_mt</th>
          <th>leiden</th>
          <th>mixscape_class_p_ko</th>
          <th>mixscape_class</th>
          <th>mixscape_class_global</th>
          <th>pertclass</th>
          <th>hdbscan</th>
          <th>new_leiden</th>
        </tr>
        <tr>
          <th>Cell_barcodes</th>
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
          <th>TAACCAGAGTAGAATC-8</th>
          <td>TRIM21</td>
          <td>3467</td>
          <td>3467</td>
          <td>10422.0</td>
          <td>755.0</td>
          <td>7.244291</td>
          <td>26</td>
          <td>1.0</td>
          <td>TRIM21 KO</td>
          <td>KO</td>
          <td>strong</td>
          <td>9</td>
          <td>110</td>
        </tr>
        <tr>
          <th>CATAGACCAACACGAG-8</th>
          <td>CBY1</td>
          <td>2003</td>
          <td>2003</td>
          <td>4621.0</td>
          <td>392.0</td>
          <td>8.483012</td>
          <td>24</td>
          <td>1.0</td>
          <td>CBY1 KO</td>
          <td>KO</td>
          <td>strong</td>
          <td>10</td>
          <td>91</td>
        </tr>
        <tr>
          <th>CTGTGAATCCGGTAAT-2</th>
          <td>LAT2</td>
          <td>4344</td>
          <td>4344</td>
          <td>16784.0</td>
          <td>1412.0</td>
          <td>8.412774</td>
          <td>9</td>
          <td>1.0</td>
          <td>LAT2 KO</td>
          <td>KO</td>
          <td>strong</td>
          <td>1</td>
          <td>484</td>
        </tr>
        <tr>
          <th>GAGCTGCAGGTAGATT-8</th>
          <td>RELA</td>
          <td>2361</td>
          <td>2360</td>
          <td>6086.0</td>
          <td>380.0</td>
          <td>6.243838</td>
          <td>16</td>
          <td>1.0</td>
          <td>RELA KO</td>
          <td>KO</td>
          <td>strong</td>
          <td>4</td>
          <td>131</td>
        </tr>
        <tr>
          <th>AAGTACCCAACTTCTT-3</th>
          <td>WT1</td>
          <td>2198</td>
          <td>2198</td>
          <td>5469.0</td>
          <td>545.0</td>
          <td>9.965259</td>
          <td>12</td>
          <td>1.0</td>
          <td>WT1 KO</td>
          <td>KO</td>
          <td>strong</td>
          <td>0</td>
          <td>412</td>
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
        </tr>
        <tr>
          <th>GGCTGTGAGGGCTAAC-5</th>
          <td>APOL2</td>
          <td>2503</td>
          <td>2503</td>
          <td>6126.0</td>
          <td>553.0</td>
          <td>9.027098</td>
          <td>30</td>
          <td>1.0</td>
          <td>APOL2 KO</td>
          <td>KO</td>
          <td>strong</td>
          <td>20</td>
          <td>219</td>
        </tr>
        <tr>
          <th>ATATCCTCATCATTTC-8</th>
          <td>TNFRSF1B</td>
          <td>4380</td>
          <td>4379</td>
          <td>14271.0</td>
          <td>1052.0</td>
          <td>7.371593</td>
          <td>11</td>
          <td>1.0</td>
          <td>TNFRSF1B KO</td>
          <td>KO</td>
          <td>strong</td>
          <td>7</td>
          <td>5</td>
        </tr>
        <tr>
          <th>CTAGGTAGTTGAGGAC-1</th>
          <td>CD27</td>
          <td>2385</td>
          <td>2385</td>
          <td>6731.0</td>
          <td>336.0</td>
          <td>4.991829</td>
          <td>19</td>
          <td>1.0</td>
          <td>CD27 KO</td>
          <td>KO</td>
          <td>strong</td>
          <td>17</td>
          <td>235</td>
        </tr>
        <tr>
          <th>TGGGAAGGTGAGTTTC-6</th>
          <td>CTRL</td>
          <td>2988</td>
          <td>2988</td>
          <td>7500.0</td>
          <td>541.0</td>
          <td>7.213333</td>
          <td>3</td>
          <td>0.0</td>
          <td>CTRL</td>
          <td>CTRL</td>
          <td>CTRL</td>
          <td>30</td>
          <td>100</td>
        </tr>
        <tr>
          <th>ATGATCGGTATCGTTG-7</th>
          <td>RELA</td>
          <td>2739</td>
          <td>2739</td>
          <td>7606.0</td>
          <td>311.0</td>
          <td>4.088877</td>
          <td>16</td>
          <td>1.0</td>
          <td>RELA KO</td>
          <td>KO</td>
          <td>strong</td>
          <td>4</td>
          <td>470</td>
        </tr>
      </tbody>
    </table>
    <p>16707 rows × 13 columns</p>
    </div>



.. code:: ipython3

    test_data.adata.obs['new_leiden'] = test_data.adata.obs['new_leiden'].astype(str)
    test_data.get_positive_marker_gene_parallel()
    rna_seq_data = scripro.SCRIPro_RNA(12,'hg38',test_data,assays=['Direct','DNase','H3K27ac'])

Calculating ISD
===============

.. code:: ipython3

    rna_seq_data.cal_ISD_cistrome()
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
          <th>factor</th>
          <th>NELFA</th>
          <th>SUPT5H</th>
          <th>POLR2A</th>
          <th>TAF1</th>
          <th>E2F1</th>
          <th>MYC</th>
          <th>JMJD6</th>
          <th>TFDP1</th>
          <th>PHF8</th>
          <th>BRD4</th>
          <th>...</th>
          <th>ESCO2</th>
          <th>SOX8</th>
          <th>WWTR1</th>
          <th>ELF5</th>
          <th>ZIC3</th>
          <th>SOX6</th>
          <th>HOXA1</th>
          <th>TOP1</th>
          <th>FOXE3</th>
          <th>ETV2</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>54</th>
          <td>1.000000</td>
          <td>0.879108</td>
          <td>0.790975</td>
          <td>0.787378</td>
          <td>0.787086</td>
          <td>0.759230</td>
          <td>0.759144</td>
          <td>0.744869</td>
          <td>0.735694</td>
          <td>0.734767</td>
          <td>...</td>
          <td>3.128389e-11</td>
          <td>2.321808e-11</td>
          <td>1.636090e-11</td>
          <td>9.684857e-12</td>
          <td>6.359623e-12</td>
          <td>4.229628e-12</td>
          <td>3.758939e-12</td>
          <td>9.764622e-13</td>
          <td>7.729002e-13</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>568</th>
          <td>0.596106</td>
          <td>0.658273</td>
          <td>0.878588</td>
          <td>0.637187</td>
          <td>0.542015</td>
          <td>0.772881</td>
          <td>0.482539</td>
          <td>0.328786</td>
          <td>0.485381</td>
          <td>0.840747</td>
          <td>...</td>
          <td>5.529090e-02</td>
          <td>2.201403e-02</td>
          <td>6.034951e-02</td>
          <td>1.250646e-01</td>
          <td>9.094418e-02</td>
          <td>2.176096e-05</td>
          <td>1.033788e-01</td>
          <td>1.759315e-02</td>
          <td>8.293918e-02</td>
          <td>1.189246e-01</td>
        </tr>
        <tr>
          <th>171</th>
          <td>0.956295</td>
          <td>0.948642</td>
          <td>0.807200</td>
          <td>0.779516</td>
          <td>0.645101</td>
          <td>1.000000</td>
          <td>0.747362</td>
          <td>0.417036</td>
          <td>0.685183</td>
          <td>0.772588</td>
          <td>...</td>
          <td>0.000000e+00</td>
          <td>9.707026e-07</td>
          <td>4.524251e-08</td>
          <td>1.419998e-10</td>
          <td>2.228434e-07</td>
          <td>2.616552e-10</td>
          <td>5.366166e-07</td>
          <td>0.000000e+00</td>
          <td>1.553591e-07</td>
          <td>1.790557e-10</td>
        </tr>
        <tr>
          <th>106</th>
          <td>1.000000</td>
          <td>0.938330</td>
          <td>0.849668</td>
          <td>0.820421</td>
          <td>0.634343</td>
          <td>0.915763</td>
          <td>0.775252</td>
          <td>0.401171</td>
          <td>0.737218</td>
          <td>0.841678</td>
          <td>...</td>
          <td>0.000000e+00</td>
          <td>1.369131e-06</td>
          <td>2.459007e-08</td>
          <td>6.081526e-07</td>
          <td>2.266960e-05</td>
          <td>6.717899e-10</td>
          <td>2.948945e-04</td>
          <td>1.449065e-13</td>
          <td>1.059525e-07</td>
          <td>4.508395e-09</td>
        </tr>
        <tr>
          <th>79</th>
          <td>0.037224</td>
          <td>0.413154</td>
          <td>1.000000</td>
          <td>0.553018</td>
          <td>0.133872</td>
          <td>0.551101</td>
          <td>0.076116</td>
          <td>0.042632</td>
          <td>0.125038</td>
          <td>0.670272</td>
          <td>...</td>
          <td>6.291131e-04</td>
          <td>1.136219e-01</td>
          <td>7.055724e-02</td>
          <td>1.315009e-01</td>
          <td>5.778941e-02</td>
          <td>7.625599e-02</td>
          <td>1.515579e-01</td>
          <td>1.035683e-02</td>
          <td>1.189053e-01</td>
          <td>2.294130e-01</td>
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
          <th>337</th>
          <td>0.628819</td>
          <td>0.884519</td>
          <td>0.881672</td>
          <td>1.000000</td>
          <td>1.000000</td>
          <td>0.738973</td>
          <td>0.544195</td>
          <td>0.935860</td>
          <td>0.466437</td>
          <td>0.784920</td>
          <td>...</td>
          <td>2.697677e-06</td>
          <td>5.354145e-03</td>
          <td>1.452867e-02</td>
          <td>1.063817e-03</td>
          <td>2.141516e-05</td>
          <td>5.709512e-02</td>
          <td>2.659353e-02</td>
          <td>3.680871e-05</td>
          <td>5.829138e-02</td>
          <td>4.284606e-02</td>
        </tr>
        <tr>
          <th>222</th>
          <td>0.912849</td>
          <td>0.800940</td>
          <td>0.908895</td>
          <td>0.974336</td>
          <td>0.590447</td>
          <td>0.717522</td>
          <td>0.606849</td>
          <td>0.280193</td>
          <td>0.566149</td>
          <td>0.913928</td>
          <td>...</td>
          <td>2.674771e-03</td>
          <td>1.304035e-02</td>
          <td>1.189110e-01</td>
          <td>6.780617e-02</td>
          <td>3.949164e-03</td>
          <td>5.605106e-02</td>
          <td>3.067985e-02</td>
          <td>5.394981e-02</td>
          <td>1.673135e-01</td>
          <td>3.413787e-02</td>
        </tr>
        <tr>
          <th>304</th>
          <td>0.839100</td>
          <td>0.804796</td>
          <td>0.861108</td>
          <td>0.937081</td>
          <td>1.000000</td>
          <td>0.718710</td>
          <td>0.672186</td>
          <td>0.722735</td>
          <td>0.624386</td>
          <td>0.913212</td>
          <td>...</td>
          <td>2.928011e-02</td>
          <td>1.334770e-03</td>
          <td>2.150550e-02</td>
          <td>1.245458e-02</td>
          <td>1.675397e-04</td>
          <td>7.373146e-03</td>
          <td>6.442231e-03</td>
          <td>1.534435e-02</td>
          <td>5.141218e-02</td>
          <td>1.081799e-02</td>
        </tr>
        <tr>
          <th>375</th>
          <td>0.145227</td>
          <td>0.599807</td>
          <td>1.000000</td>
          <td>0.570559</td>
          <td>0.201262</td>
          <td>0.653092</td>
          <td>0.182292</td>
          <td>0.047675</td>
          <td>0.259976</td>
          <td>0.752675</td>
          <td>...</td>
          <td>5.605888e-04</td>
          <td>3.816154e-02</td>
          <td>3.498265e-02</td>
          <td>1.408408e-01</td>
          <td>2.847721e-02</td>
          <td>7.229822e-02</td>
          <td>9.771098e-02</td>
          <td>4.832158e-03</td>
          <td>5.258527e-02</td>
          <td>2.155378e-01</td>
        </tr>
        <tr>
          <th>220</th>
          <td>0.991618</td>
          <td>0.938332</td>
          <td>0.838142</td>
          <td>0.710450</td>
          <td>0.641620</td>
          <td>1.000000</td>
          <td>0.719913</td>
          <td>0.371703</td>
          <td>0.608496</td>
          <td>0.759212</td>
          <td>...</td>
          <td>1.202406e-14</td>
          <td>3.684495e-06</td>
          <td>4.245327e-07</td>
          <td>5.762894e-07</td>
          <td>4.323174e-08</td>
          <td>3.048350e-11</td>
          <td>7.881787e-09</td>
          <td>6.787328e-10</td>
          <td>6.295932e-06</td>
          <td>1.545003e-10</td>
        </tr>
      </tbody>
    </table>
    <p>592 rows × 1252 columns</p>
    </div>



.. code:: ipython3

    rna_seq_data.get_tf_score()
    tem_exp = rna_raw.to_df().merge(test_data.adata.obs.loc[:,'new_leiden'],left_index=True,right_index=True)
    grouped = tem_exp.groupby('new_leiden').mean()
    grouped




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
          <th>TNFRSF9-1</th>
          <th>TNFRSF9-2</th>
          <th>TRAF3IP2-1</th>
          <th>TRAF3IP2-2</th>
          <th>TRIM21-1</th>
          <th>TRIM21-2</th>
          <th>VAV1-1</th>
          <th>VAV1-2</th>
          <th>WT1-1</th>
          <th>WT1-2</th>
        </tr>
        <tr>
          <th>new_leiden</th>
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
          <th>0</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.022986</td>
          <td>0.020883</td>
          <td>0.141023</td>
          <td>0.019764</td>
          <td>0.070329</td>
          <td>0.045294</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>5.971172</td>
          <td>0.025028</td>
        </tr>
        <tr>
          <th>1</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.0</td>
          <td>0.028059</td>
          <td>0.000000</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.307047</td>
          <td>0.000000</td>
          <td>0.175104</td>
          <td>0.078949</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.031357</td>
          <td>3.475705</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>2</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.029139</td>
          <td>0.000000</td>
          <td>0.0</td>
          <td>0.061611</td>
          <td>0.029139</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.145290</td>
          <td>0.000000</td>
          <td>0.025556</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.082464</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>3</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.0</td>
          <td>0.030975</td>
          <td>0.000000</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.025993</td>
          <td>0.000000</td>
          <td>0.171809</td>
          <td>0.099170</td>
          <td>0.072312</td>
          <td>0.000000</td>
          <td>0.039627</td>
          <td>0.000000</td>
          <td>0.352380</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>4</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.0</td>
          <td>0.035995</td>
          <td>0.000000</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.032857</td>
          <td>0.000000</td>
          <td>0.288488</td>
          <td>0.000000</td>
          <td>0.043905</td>
          <td>0.044607</td>
          <td>0.000000</td>
          <td>0.035995</td>
          <td>0.072219</td>
          <td>0.000000</td>
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
          <th>587</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.037087</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.052179</td>
          <td>0.000000</td>
          <td>0.828584</td>
          <td>0.035307</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.094569</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>588</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.035962</td>
          <td>0.000000</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.039119</td>
          <td>0.000000</td>
          <td>1.642812</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.074946</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>589</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.141428</td>
          <td>0.029187</td>
          <td>0.287361</td>
          <td>0.202327</td>
          <td>0.059597</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.024091</td>
          <td>0.094148</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>590</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.031027</td>
          <td>0.000000</td>
          <td>0.155326</td>
          <td>0.051923</td>
          <td>0.000000</td>
          <td>0.105978</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.032125</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>-1</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.003845</td>
          <td>0.000000</td>
          <td>0.0</td>
          <td>0.009864</td>
          <td>0.000000</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.070795</td>
          <td>0.030611</td>
          <td>0.236199</td>
          <td>0.057267</td>
          <td>0.053324</td>
          <td>0.301895</td>
          <td>0.023147</td>
          <td>0.022444</td>
          <td>0.200147</td>
          <td>0.008683</td>
        </tr>
      </tbody>
    </table>
    <p>592 rows × 36755 columns</p>
    </div>



.. code:: ipython3

    rna_seq_data.Ori_Data.ad_all = grouped
    rna_seq_data.Ori_Data.super_gene_exp = grouped
    super_gene_exp = rna_seq_data.Ori_Data.super_gene_exp
    super_gene_mean = rna_seq_data.Ori_Data.super_gene_mean
    super_gene_std = rna_seq_data.Ori_Data.super_gene_std
    rna_seq_data.Ori_Data.super_gene_mean = rna_seq_data.Ori_Data.super_gene_exp.mean()
    rna_seq_data.Ori_Data.super_gene_std = rna_seq_data.Ori_Data.super_gene_exp.std()
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
          <th>-1</th>
          <td>8.563566e-04</td>
          <td>0.360750</td>
          <td>0.559371</td>
          <td>0.182508</td>
          <td>0.018380</td>
          <td>0.000006</td>
          <td>0.000029</td>
          <td>1.898788e-05</td>
          <td>0.000016</td>
          <td>0.000330</td>
          <td>...</td>
          <td>0.305053</td>
          <td>0.001561</td>
          <td>0.114077</td>
          <td>0.000631</td>
          <td>1.361497e-03</td>
          <td>3.585991e-06</td>
          <td>0.006103</td>
          <td>2.655798e-05</td>
          <td>0.261396</td>
          <td>1.215443e-01</td>
        </tr>
        <tr>
          <th>0</th>
          <td>3.446237e-12</td>
          <td>0.497215</td>
          <td>0.500529</td>
          <td>0.163227</td>
          <td>0.000100</td>
          <td>0.000003</td>
          <td>0.000129</td>
          <td>1.646497e-07</td>
          <td>0.000002</td>
          <td>0.000044</td>
          <td>...</td>
          <td>0.349755</td>
          <td>0.002214</td>
          <td>0.028432</td>
          <td>0.000275</td>
          <td>1.026750e-07</td>
          <td>3.986286e-06</td>
          <td>0.034384</td>
          <td>1.379630e-05</td>
          <td>0.219834</td>
          <td>1.525143e-12</td>
        </tr>
        <tr>
          <th>1</th>
          <td>7.096883e-02</td>
          <td>0.480607</td>
          <td>0.591523</td>
          <td>0.193952</td>
          <td>0.028921</td>
          <td>0.051201</td>
          <td>0.000002</td>
          <td>4.352112e-04</td>
          <td>0.000222</td>
          <td>0.000107</td>
          <td>...</td>
          <td>0.309516</td>
          <td>0.000656</td>
          <td>0.068821</td>
          <td>0.001493</td>
          <td>5.858211e-06</td>
          <td>8.907871e-07</td>
          <td>0.008626</td>
          <td>1.288200e-08</td>
          <td>0.212142</td>
          <td>2.165690e-01</td>
        </tr>
        <tr>
          <th>10</th>
          <td>5.333758e-02</td>
          <td>0.544322</td>
          <td>0.600351</td>
          <td>0.262902</td>
          <td>0.029763</td>
          <td>0.057495</td>
          <td>0.056367</td>
          <td>7.002899e-04</td>
          <td>0.000222</td>
          <td>0.029509</td>
          <td>...</td>
          <td>0.309205</td>
          <td>0.004910</td>
          <td>0.059540</td>
          <td>0.001900</td>
          <td>4.214787e-02</td>
          <td>3.509326e-07</td>
          <td>0.021020</td>
          <td>1.045370e-04</td>
          <td>0.176866</td>
          <td>2.342755e-01</td>
        </tr>
        <tr>
          <th>100</th>
          <td>7.107864e-02</td>
          <td>0.399857</td>
          <td>0.764365</td>
          <td>0.176755</td>
          <td>0.121364</td>
          <td>0.063864</td>
          <td>0.032153</td>
          <td>3.981877e-02</td>
          <td>0.026263</td>
          <td>0.010673</td>
          <td>...</td>
          <td>0.298497</td>
          <td>0.011864</td>
          <td>0.116632</td>
          <td>0.019865</td>
          <td>7.758212e-02</td>
          <td>9.608021e-03</td>
          <td>0.008497</td>
          <td>1.111107e-04</td>
          <td>0.092073</td>
          <td>1.953459e-01</td>
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
          <th>95</th>
          <td>1.291925e-01</td>
          <td>0.487603</td>
          <td>0.791613</td>
          <td>0.203145</td>
          <td>0.080655</td>
          <td>0.063231</td>
          <td>0.036989</td>
          <td>3.843538e-02</td>
          <td>0.045159</td>
          <td>0.022616</td>
          <td>...</td>
          <td>0.353437</td>
          <td>0.004247</td>
          <td>0.114041</td>
          <td>0.020601</td>
          <td>2.575347e-02</td>
          <td>5.994760e-05</td>
          <td>0.023438</td>
          <td>1.904101e-02</td>
          <td>0.140363</td>
          <td>9.592408e-02</td>
        </tr>
        <tr>
          <th>96</th>
          <td>1.031648e-01</td>
          <td>0.318737</td>
          <td>0.647392</td>
          <td>0.119693</td>
          <td>0.108568</td>
          <td>0.019088</td>
          <td>0.034342</td>
          <td>2.822314e-02</td>
          <td>0.049918</td>
          <td>0.056430</td>
          <td>...</td>
          <td>0.194852</td>
          <td>0.011646</td>
          <td>0.117273</td>
          <td>0.023761</td>
          <td>2.218386e-02</td>
          <td>3.685638e-03</td>
          <td>0.013456</td>
          <td>6.609630e-03</td>
          <td>0.204090</td>
          <td>1.425852e-01</td>
        </tr>
        <tr>
          <th>97</th>
          <td>1.062549e-01</td>
          <td>0.340830</td>
          <td>0.789539</td>
          <td>0.114368</td>
          <td>0.136252</td>
          <td>0.097897</td>
          <td>0.059248</td>
          <td>8.784248e-02</td>
          <td>0.156568</td>
          <td>0.042775</td>
          <td>...</td>
          <td>0.191674</td>
          <td>0.002716</td>
          <td>0.095484</td>
          <td>0.016832</td>
          <td>5.122700e-02</td>
          <td>8.143724e-02</td>
          <td>0.000837</td>
          <td>2.209640e-02</td>
          <td>0.089183</td>
          <td>1.023146e-01</td>
        </tr>
        <tr>
          <th>98</th>
          <td>6.693196e-02</td>
          <td>0.558301</td>
          <td>0.756792</td>
          <td>0.281248</td>
          <td>0.047730</td>
          <td>0.104862</td>
          <td>0.089017</td>
          <td>1.252955e-01</td>
          <td>0.083650</td>
          <td>0.035471</td>
          <td>...</td>
          <td>0.302286</td>
          <td>0.002517</td>
          <td>0.147804</td>
          <td>0.007596</td>
          <td>2.793163e-02</td>
          <td>4.735611e-03</td>
          <td>0.019666</td>
          <td>2.356419e-02</td>
          <td>0.289258</td>
          <td>1.428439e-01</td>
        </tr>
        <tr>
          <th>99</th>
          <td>6.862545e-02</td>
          <td>0.496669</td>
          <td>0.600869</td>
          <td>0.187158</td>
          <td>0.096891</td>
          <td>0.030484</td>
          <td>0.006576</td>
          <td>2.917082e-03</td>
          <td>0.001573</td>
          <td>0.015501</td>
          <td>...</td>
          <td>0.390762</td>
          <td>0.018814</td>
          <td>0.083292</td>
          <td>0.019330</td>
          <td>3.054598e-02</td>
          <td>4.523016e-06</td>
          <td>0.027676</td>
          <td>5.116583e-04</td>
          <td>0.212204</td>
          <td>2.528735e-01</td>
        </tr>
      </tbody>
    </table>
    <p>592 rows × 1226 columns</p>
    </div>



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
          <th>-1</th>
          <td>3.672389e-04</td>
          <td>0.178091</td>
          <td>0.239622</td>
          <td>0.070631</td>
          <td>0.007116</td>
          <td>1.357711e-06</td>
          <td>0.0</td>
          <td>2.845156e-06</td>
          <td>2.377675e-06</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.063907</td>
          <td>0.000153</td>
          <td>0.028755</td>
          <td>0.000139</td>
          <td>3.560237e-04</td>
          <td>0.0</td>
          <td>0.001316</td>
          <td>0.0</td>
          <td>0.052378</td>
          <td>4.111837e-02</td>
        </tr>
        <tr>
          <th>0</th>
          <td>2.414882e-12</td>
          <td>0.261679</td>
          <td>0.205148</td>
          <td>0.078429</td>
          <td>0.000027</td>
          <td>5.387451e-07</td>
          <td>0.0</td>
          <td>2.461573e-08</td>
          <td>2.656785e-07</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.126825</td>
          <td>0.000271</td>
          <td>0.014292</td>
          <td>0.000077</td>
          <td>1.620604e-08</td>
          <td>0.0</td>
          <td>0.003788</td>
          <td>0.0</td>
          <td>0.063796</td>
          <td>1.009363e-12</td>
        </tr>
        <tr>
          <th>1</th>
          <td>5.021123e-02</td>
          <td>0.086393</td>
          <td>0.228981</td>
          <td>0.160712</td>
          <td>0.007152</td>
          <td>8.023564e-03</td>
          <td>0.0</td>
          <td>7.339054e-05</td>
          <td>3.607554e-05</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.040119</td>
          <td>0.000085</td>
          <td>0.005075</td>
          <td>0.000094</td>
          <td>9.613452e-07</td>
          <td>0.0</td>
          <td>0.000946</td>
          <td>0.0</td>
          <td>0.085753</td>
          <td>6.809414e-02</td>
        </tr>
        <tr>
          <th>10</th>
          <td>3.808361e-02</td>
          <td>0.175952</td>
          <td>0.306938</td>
          <td>0.168117</td>
          <td>0.009945</td>
          <td>1.062668e-02</td>
          <td>0.0</td>
          <td>1.275641e-04</td>
          <td>4.347075e-05</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.037240</td>
          <td>0.004910</td>
          <td>0.033682</td>
          <td>0.000281</td>
          <td>6.428773e-03</td>
          <td>0.0</td>
          <td>0.002387</td>
          <td>0.0</td>
          <td>0.081476</td>
          <td>8.746815e-02</td>
        </tr>
        <tr>
          <th>100</th>
          <td>3.513059e-02</td>
          <td>0.293045</td>
          <td>0.449179</td>
          <td>0.129605</td>
          <td>0.061245</td>
          <td>8.614222e-03</td>
          <td>0.0</td>
          <td>5.562580e-03</td>
          <td>3.426428e-03</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.019991</td>
          <td>0.001214</td>
          <td>0.053916</td>
          <td>0.005434</td>
          <td>9.358260e-03</td>
          <td>0.0</td>
          <td>0.000642</td>
          <td>0.0</td>
          <td>0.008968</td>
          <td>1.053210e-01</td>
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
          <th>95</th>
          <td>1.177092e-01</td>
          <td>0.209708</td>
          <td>0.538176</td>
          <td>0.134788</td>
          <td>0.020214</td>
          <td>1.010271e-02</td>
          <td>0.0</td>
          <td>6.525869e-03</td>
          <td>7.493207e-03</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.028318</td>
          <td>0.000345</td>
          <td>0.100466</td>
          <td>0.007886</td>
          <td>3.624425e-03</td>
          <td>0.0</td>
          <td>0.001972</td>
          <td>0.0</td>
          <td>0.081238</td>
          <td>4.317646e-02</td>
        </tr>
        <tr>
          <th>96</th>
          <td>2.874359e-02</td>
          <td>0.148875</td>
          <td>0.243910</td>
          <td>0.046832</td>
          <td>0.029400</td>
          <td>2.040194e-03</td>
          <td>0.0</td>
          <td>2.741001e-03</td>
          <td>4.072542e-03</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.058285</td>
          <td>0.000712</td>
          <td>0.037024</td>
          <td>0.010316</td>
          <td>2.618379e-03</td>
          <td>0.0</td>
          <td>0.000991</td>
          <td>0.0</td>
          <td>0.020238</td>
          <td>5.766247e-02</td>
        </tr>
        <tr>
          <th>97</th>
          <td>1.906883e-02</td>
          <td>0.174282</td>
          <td>0.399717</td>
          <td>0.026488</td>
          <td>0.035957</td>
          <td>1.174490e-02</td>
          <td>0.0</td>
          <td>1.196944e-02</td>
          <td>1.987612e-02</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.010206</td>
          <td>0.000062</td>
          <td>0.016045</td>
          <td>0.006168</td>
          <td>6.485161e-03</td>
          <td>0.0</td>
          <td>0.000341</td>
          <td>0.0</td>
          <td>0.004761</td>
          <td>2.665132e-03</td>
        </tr>
        <tr>
          <th>98</th>
          <td>1.754907e-02</td>
          <td>0.471482</td>
          <td>0.416241</td>
          <td>0.115289</td>
          <td>0.022126</td>
          <td>1.411043e-02</td>
          <td>0.0</td>
          <td>1.252955e-01</td>
          <td>1.405428e-02</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.021562</td>
          <td>0.000198</td>
          <td>0.069066</td>
          <td>0.005385</td>
          <td>4.197360e-03</td>
          <td>0.0</td>
          <td>0.010347</td>
          <td>0.0</td>
          <td>0.024765</td>
          <td>8.091507e-02</td>
        </tr>
        <tr>
          <th>99</th>
          <td>2.922981e-02</td>
          <td>0.146265</td>
          <td>0.232052</td>
          <td>0.074721</td>
          <td>0.056472</td>
          <td>3.580166e-03</td>
          <td>0.0</td>
          <td>4.096966e-04</td>
          <td>2.107131e-04</td>
          <td>0.0</td>
          <td>...</td>
          <td>0.029510</td>
          <td>0.001427</td>
          <td>0.026408</td>
          <td>0.005165</td>
          <td>3.396306e-03</td>
          <td>0.0</td>
          <td>0.001820</td>
          <td>0.0</td>
          <td>0.108117</td>
          <td>1.003711e-01</td>
        </tr>
      </tbody>
    </table>
    <p>592 rows × 1226 columns</p>
    </div>



Calculate the AUPRC and AUROC
=============================

.. code:: ipython3

    scripro_score = test_data.adata.obs.merge(rna_seq_data.tf_score,left_on='new_leiden',right_index=True).iloc[:,13:]
    commontf = set(test_data.adata.obs['gene']).intersection(set(scripro_score.columns))
    scripro_auroc_dic = {}
    for k in commontf:
        y_true = []
        for i in scripro_score.index:
            if test_data.adata.obs.loc[i,'gene'] == k:
                y_true.append(1)
            else: 
                y_true.append(0)
        y_scores = list(scripro_score.loc[:,k])
        fpr, tpr, thresholds = roc_curve(y_true, y_scores)
        roc_auc = auc(fpr, tpr)
        scripro_auroc_dic[k]=roc_auc
    
    scripro_auroc_score = pd.DataFrame([scripro_auroc_dic]).T.sort_values(ascending = False,by = 0)
    scripro_auroc_score.columns = ['auroc']
    scripro_auroc_score




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
          <th>auroc</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>EOMES</th>
          <td>0.993938</td>
        </tr>
        <tr>
          <th>GATA3</th>
          <td>0.951226</td>
        </tr>
        <tr>
          <th>RELA</th>
          <td>0.941175</td>
        </tr>
        <tr>
          <th>FOXD2</th>
          <td>0.916245</td>
        </tr>
        <tr>
          <th>PRDM1</th>
          <td>0.915120</td>
        </tr>
        <tr>
          <th>TBX21</th>
          <td>0.863753</td>
        </tr>
        <tr>
          <th>LHX4</th>
          <td>0.818389</td>
        </tr>
        <tr>
          <th>FOXQ1</th>
          <td>0.743893</td>
        </tr>
        <tr>
          <th>LHX6</th>
          <td>0.735591</td>
        </tr>
        <tr>
          <th>WT1</th>
          <td>0.729839</td>
        </tr>
        <tr>
          <th>JMJD1C</th>
          <td>0.691111</td>
        </tr>
        <tr>
          <th>ALX4</th>
          <td>0.662308</td>
        </tr>
        <tr>
          <th>NOTCH1</th>
          <td>0.583613</td>
        </tr>
        <tr>
          <th>IKZF3</th>
          <td>0.564741</td>
        </tr>
        <tr>
          <th>FOSB</th>
          <td>0.413174</td>
        </tr>
      </tbody>
    </table>
    </div>



.. code:: ipython3

    import pandas as pd
    from sklearn.metrics import precision_recall_curve, auc
    
    scripro_auprc_dic = {}
    for k in commontf:
        y_true = []
        for i in scripro_score.index:
            if test_data.adata.obs.loc[i, 'gene'] == k:
                y_true.append(1)
            else: 
                y_true.append(0)
        y_scores = list(scripro_score.loc[:,k])
        precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
        auprc = auc(recall, precision)
        scripro_auprc_dic[k] = auprc
    scripro_auprc_score = pd.DataFrame([scripro_auprc_dic]).T.sort_values(ascending=False, by=0)
    scripro_auprc_score.columns = ['auprc']
    scripro_auprc_score




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
          <th>auprc</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>EOMES</th>
          <td>0.871598</td>
        </tr>
        <tr>
          <th>RELA</th>
          <td>0.760414</td>
        </tr>
        <tr>
          <th>GATA3</th>
          <td>0.731483</td>
        </tr>
        <tr>
          <th>WT1</th>
          <td>0.598684</td>
        </tr>
        <tr>
          <th>LHX6</th>
          <td>0.495087</td>
        </tr>
        <tr>
          <th>FOXD2</th>
          <td>0.280925</td>
        </tr>
        <tr>
          <th>TBX21</th>
          <td>0.097302</td>
        </tr>
        <tr>
          <th>ALX4</th>
          <td>0.061300</td>
        </tr>
        <tr>
          <th>FOXQ1</th>
          <td>0.003539</td>
        </tr>
        <tr>
          <th>JMJD1C</th>
          <td>0.000962</td>
        </tr>
        <tr>
          <th>PRDM1</th>
          <td>0.000350</td>
        </tr>
        <tr>
          <th>FOSB</th>
          <td>0.000324</td>
        </tr>
        <tr>
          <th>NOTCH1</th>
          <td>0.000202</td>
        </tr>
        <tr>
          <th>LHX4</th>
          <td>0.000163</td>
        </tr>
        <tr>
          <th>IKZF3</th>
          <td>0.000128</td>
        </tr>
      </tbody>
    </table>
    </div>



