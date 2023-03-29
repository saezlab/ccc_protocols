Frequently Asked Questions
--------------------------

1. When to use LIANA?
==================================

LIANA is a ligand-receptor inference tool that can be used on any type of single-cell data following preprocessing and the annotation of the cells into biologically-meaningful cell types.
It enables the use of different ligand-receptor methods & resources, and also provides a consensus for both.

2. When to use Tensor-cell2cell?
==================================
Tensor-cell2cell, put simply, is a factorization approach that can be used to reduce the dimensionality of inferred ligand-receptor interactions across samples.
This is done by decomposing the inferred ligand-receptor interactions into a set of latent factors that can be used to reconstruct the original tensor.
This serves as a way to improve the interpretability of the inferred interactions, and to obtain coordinates intercellular programmes that can distinguish patients samples, or other contexts.

3. When to use LIANA & Tensor-cell2cell?
==================================

LIANA & Tensor-cell2cell can be used together to infer ligand-receptor interactions across samples, and then to reduce the dimensionality of the inferred interactions.
They are thus fit for use when working with single-cell atlases that contain multiple samples (coming from e.g. multiple patients or time points).


4. What kind of data can be used with LIANA & Tensor-cell2cell
========================================================

LIANA typically works with log1p normalized data and requires annotated cell groups. Tensor-cell2cell then takes LIANA's output as input.


5. What kind of experimental design can be used with Tensor-cell2cell
========================================================================

In theory, Tensor-cell2cell would work with any experimental design. 
It will then identify hidden factors associated with the experimental conditions in an untargeted manner, given the condition at hand drives differences between contexts/samples.