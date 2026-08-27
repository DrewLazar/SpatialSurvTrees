simulations/sweep/
The robustness sweep of Section 6.3. Starting from the three-zone design,
the clinical covariate `tpi` is progressively coupled to membership in the
high-frailty zone, so that the exogeneity condition fails by degrees.
Eight coupling levels, 30 replicate datasets each, three analyses per
replicate: two-stage, spatial-only, and the `spBayesSurv` PH+GRF
benchmark. The two arms are fit separately, giving 480 fits in total.
File	Role
`sim\_sweep\_generate.R`	Builds one replicate at a given coupling level. `make\_sweep\_data()`.
`sim\_sweep\_analysis.R`	Fits and scores one replicate. `collect\_sweep()`, `summarise\_sweep()`.
`run\_sweep\_node.R`	Cluster worker, tree arm. `Rscript run\_sweep\_node.R <1..20>`
`run\_sweep\_node\_sb.R`	Cluster worker, `spBayesSurv` arm.
`fig\_sweep\_result.R`	Collates and draws `fig\_sweep\_result.pdf`.
`fig\_contamination\_sidebyside.R`	The confounding mechanism at the low and high ends. Not a figure in the paper, but it shows what the sweep varies.
`sweep\_all.rds`	All 480 replicates, one row each.
