# method/

Implementation of the node-splitting survival tree and the two-stage
procedure. Everything in 'application/' and 'simulations/' sources from here.

| File | Contents |
|---|---|
| 'dipole_tree.R' | The tree class. Dipole classification, kernel splitting criterion, orientation algorithm, growing, prediction. |
| 'bootstrap_pruning.R' | Split-complexity pruning with bootstrap selection of the complexity parameter. |
| 'computeresiduals_NA.R' | Leaf Nelson-Aalen cumulative hazard residuals, with censoring carried over. |
| 'metrics.R' | Zone-recovery scoring: accuracy, far-boundary accuracy, ARI, standardized MSE, span ratio, AUC. |

## The tree class

'DipolarSurvivalTree_MixedGaussPolyLinearKernel_recpom' differs from the
criterion of Maung et al. in two respects, both described in Section 3.4 of
the paper. The pure-or-mixed dipole classification is recomputed at every
non-terminal node from the observations in that node, rather than once at
the root. The margin is likewise recomputed locally from feature-space
distances, floored at a proportion of the root value, controlled by
'adaptive', 'probs' and 'epsilon_floor_alpha'.

The spatial kernel is a mixture of Gaussian and polynomial components,
weighted by 'gaussweight' and 'polyweight'. Setting 'linearweight = 1' with
the other two at zero recovers oblique splits, which is what the clinical
stage uses.

## Usage

```r
source("method/dipole_tree.R")
source("method/bootstrap_pruning.R")
source("method/computeresiduals_NA.R")
source("method/metrics.R")

Dip <- DipolarSurvivalTree_MixedGaussPolyLinearKernel_recpom$new(
  data, time, censor, covariates, quantiles, tolerance,
  epsilon = epsilon, kappa = exp(-2), nsize = nsize,
  pureweight = 1, mixedweight = 1, metric = "kernel",
  Ksigma = ksigma, Kconstant = 1, Kpoly_order = 2,
  gaussweight = 0.60, polyweight = 0.40, linearweight = 0,
  gausscovariates_index = 1:2, polycovariates_index = 1:2,
  linearcovariates_index = 1:2, ncovariatestosearch = 2,
  adaptive = TRUE, epsilon_floor_alpha = 0.2, probs = 0.23)

tree <- Dip$createtree(1:n)
tree <- bootstrapPruning(tree, Dip, 3.7, time = time, censor = censor)[[4]]
pred <- Dip$predicttime(data, tree)
```

The two-stage procedure is this twice: a linear-kernel tree on the clinical
covariates, then 'computeresiduals_NA' on its leaves, then a mixed-kernel
tree fit to those residuals on the coordinates.

## Dependencies

'R6', 'data.tree', 'quadprog' (or 'osqp'), 'fields', survival'.
