.. LIANA x Tensor-cell2cell documentation master file, created by
   sphinx-quickstart on Tue Feb 28 14:38:49 2023.

LIANA x Tensor-cell2cell tutorials
----------------------------

Documentation for the combined use of LIANA and Tensor-cell2cell to analyze cell-cell communication events from omics data.


Background
==========

In recent years, the data-driven inference of cell-cell communication (CCC), 
specifically when using single-cell transcriptomics data, has enabled the study of coordinated biological processes across cell types. 
Yet, as the capabilities to generate large single-cell and spatial transcriptomics datasets continue to increase,
together with the interest in studying intercellular programmes, the need to easily and robustly decipher CCC is essential.
Here, we integrate our tools, LIANA and Tensor-cell2cell, to enable identification of intercellular programmes across multiple samples and contexts. 
We show how our unified framework facilitates the choice of method to infer cell-cell communication 
and the application of factor analysis to obtain and summarize biological insights. 


.. image:: _static/intro.png
  :width: 550


Tool Availability
=================

Tensor-cell2cell is available at:
https://github.com/earmingol/cell2cell


LIANA is available in:
Python: https://github.com/saezlab/liana-py
R: https://github.com/saezlab/liana


References
==========

To cite this work: 
Baghdassarian, H. M., Dimitrov, D., Armingol, E., Saez-Rodriguez, J. & Lewis, N. E. Combining LIANA and Tensor-cell2cell to decipher cell-cell communication across multiple samples. BioRxiv (2023) doi:10.1101/2023.04.28.538731.

To cite LIANA:
Dimitrov, D., Türei, D., Garrido-Rodriguez, M., Burmedi, P.L., Nagai, J.S., Boys, C., Ramirez Flores, R.O., Kim, H., Szalai, B., Costa, I.G. and Valdeolivas, A., 2022. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nature Communications, 13(1), p.3224.

To cite Tensor-cell2cell:
Armingol, E., Baghdassarian, H.M., Martino, C., Perez-Lopez, A., Aamodt, C., Knight, R. and Lewis, N.E., 2022. Context-aware deconvolution of cell–cell communication with Tensor-cell2cell. Nature communications, 13(1), p.3665.



.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: COVID-19 - QuickStart Tutorials (BALF)

   notebooks/ccc_python/QuickStart.ipynb

   notebooks/ccc_R/QuickStart.ipynb

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Lupus - QuickStart Tutorials (PBMCs)

   notebooks/ccc_python/QuickStart_pbmc.ipynb

   notebooks/ccc_R/QuickStart_pbmc.ipynb


.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Python tutorials

   notebooks/ccc_python/01-Preprocess-Expression.ipynb

   notebooks/ccc_python/02-Infer-Communication-Scores.ipynb

   notebooks/ccc_python/03-Generate-Tensor.ipynb

   notebooks/ccc_python/04-Perform-Tensor-Factorization.ipynb

   notebooks/ccc_python/05-Downstream-Visualizations.ipynb

   notebooks/ccc_python/06-Functional-Footprinting.ipynb


.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: R tutorials

   notebooks/ccc_R/01-Preprocess-Expression.ipynb

   notebooks/ccc_R/02-Infer-Communication-Scores.ipynb

   notebooks/ccc_R/03-Generate-Tensor.ipynb

   notebooks/ccc_R/04-Perform-Tensor-Factorization.ipynb

   notebooks/ccc_R/05-Downstream-Visualizations.ipynb

   notebooks/ccc_R/06-Functional-Footprinting.ipynb


.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Further Information

   ./faq.rst
   
.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Python Supplements
   
   notebooks/ccc_python/S1_Batch_Correction.ipynb

   notebooks/ccc_python/S2_Interoperability.ipynb
   
   notebooks/ccc_python/S3A_Score_Consistency.ipynb
   
   notebooks/ccc_python/S3B_Score_Consistency.ipynb

   notebooks/ccc_python/S4_Spatial-Decomposition.ipynb 
   

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: R Supplements
   
   notebooks/ccc_R/S1_Batch_Correction.ipynb

   notebooks/ccc_R/S2_Interoperability.ipynb
   
   notebooks/ccc_R/S3_Score_Consistency.ipynb
