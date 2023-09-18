# 3D Array Colocalization
Code for quantifying colocalization of diffraction-limited arrays in 3D.

From Watson, J.L. & Krüger, L.K. et al. 2023

Synthetic Par polarity induces cytoskeleton asymmetry in unpolarized mammalian cells

### Overview ###

This code takes multichannel xyz images and assesses colocalization between the different channels, of diffraction-limited spots in 3D.

The main function for Gaussian fitting and detecting colocalization is `fit_func.m`. This takes image and region of interest (ROI) inputs (the latter can be generated in imageJ), along with analysis- and microscope-specific inputs.

`Colloc_3D_batch.m` can we used to execute `fit_func.m` on a whole directory of images.

Cite: Watson, J.L. & Krüger, L.K. et al. Synthetic Par polarity induces cytoskeleton asymmetry in unpolarized mammalian cells. Cell, 2023
