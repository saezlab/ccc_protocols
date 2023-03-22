.. LIANA x Tensor-c2c documentation master file, created by
   sphinx-quickstart on Tue Feb 28 14:38:49 2023.

LIANA x Tensor-c2c tutorials
----------------------------

Documentation for the combined use of LIANA and Tensor-cell2cell to analyze cell-cell communication events from omics data.


Tool Availability
=================

Tensor-cell2cell is available at:
https://github.com/earmingol/cell2cell


LIANA is available in:
- R: https://github.com/saezlab/liana
- Python: https://github.com/saezlab/liana-py


References
==========

To cite this work: 
TODO Baghdassarian, ...

To cite LIANA:
Dimitrov, D., Türei, D., Garrido-Rodriguez, M., Burmedi, P.L., Nagai, J.S., Boys, C., Ramirez Flores, R.O., Kim, H., Szalai, B., Costa, I.G. and Valdeolivas, A., 2022. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nature Communications, 13(1), p.3224.

To cite Tensor-c2c:
Armingol, E., Baghdassarian, H.M., Martino, C., Perez-Lopez, A., Aamodt, C., Knight, R. and Lewis, N.E., 2022. Context-aware deconvolution of cell–cell communication with Tensor-cell2cell. Nature communications, 13(1), p.3665.



.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: QuickStart Tutorials

   notebooks/ccc_python/QuickStart.ipynb

   notebooks/ccc_R/QuickStart.ipynb


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

   notebooks/ccc_python/S1_Interoperability.ipynb

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

   notebooks/ccc_R/S1_Interoperability.ipynb


.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Further Information

   ./faq.rst