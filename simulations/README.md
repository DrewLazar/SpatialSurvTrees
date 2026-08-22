# simulations/

Four designs with a known spatial partition, plus the robustness sweep.
Table 2 and Figures 2 to 5 of the paper.

Each design is a pair. The generate script builds the data and leaves the
truth objects in the environment; the analysis script fits all three
methods, scores them, and draws the recovery figure. 

| Design | Scripts | Condition |
|---|---|---|
| Degenerate | `sim_degenerate_*.R` | Covariates independent of location, `S ⊥ X` |
| Two-zone | `sim_zperpx_*.R` | Exogeneity holds, `tpi` spatially structured |
| Three-zone | `sim_threezone_*_v2.R` | Exogeneity, three nested zones |
| Four-zone | `sim_fourzone_*_v2.R` | Exogeneity, four nested zones |
| Sweep | `sweep/` | Graded violation of exogeneity, 8 levels x 30 replicates |

## Common settings

All designs use the real LeukSurv coordinates, so the spatial support is
irregular. Hyperparameters are held fixed across designs rather than tuned.
The spatial-only and two-stage analyses are grown and pruned identically, so 
they differ only in the clinical adjustment.

The clinical hazard carries a crossing age effect whose log hazard ratio
changes sign over follow-up. This violates proportional hazards, which the
`spBayesSurv` benchmark assumes and Theorem 1 does not require.

## Metrics

Scoring comes from `../method/metrics.R` so that every design is scored
identically. Far-boundary performance is reported as accuracy among the
half of observations farthest from a true boundary, matching Table 2, where
higher is better in every column except standardized MSE.

## Running

```r
source("sim_degenerate_generate.R")
source("sim_degenerate_analysis.R")
```
