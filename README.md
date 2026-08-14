# SpatialSurvTrees

A two-stage nonparametric method for recovering geographic variation in survival
risk after adjustment for clinical covariates. A survival tree fit to the clinical
covariates alone supplies Nelson-Aalen cumulative hazard residuals, and a second
tree fit to those residuals on the spatial coordinates returns a discrete risk map
with sharp, possibly curved boundaries, imposing neither a functional form on the
clinical hazard nor smoothness on the spatial effect.

Code for **"Separating Spatial and Clinical Risk with Node-Splitting SVM Survival
Trees"**, Drew Lazar and Aye Aye Maung. [arXiv:XXXX.XXXXX](https://arxiv.org/abs/XXXX.XXXXX)
