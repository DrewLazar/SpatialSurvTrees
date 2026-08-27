application/
The LeukSurv analysis of Section 5. Two scripts, run in order.
File	Role
`leuksurv_realdata.R`	Fits everything: the `spBayesSurv` PH+GRF benchmark, the same model without the spatial term, the two-stage method, and an unadjusted spatial-only tree. Produces Table 1 and the numbers quoted in Section 5.
`leuksurv_figure_4panel.R`	Draws `leuksurv_comparison_4panel.pdf` from the fitted objects.
```r
source("leuksurv_realdata.R")
source("leuksurv_figure_4panel.R")
```
The first script ends by writing two artifacts the second reads:
`alldata_all_models.rds` — LeukSurv plus every fitted column, so the
figure script runs from that file alone
`twostage_session.RData` — `surf` and `spatial_risk`, needed only for
the smooth benchmark panel. Without them the figure falls back to
three panels rather than four.
Both are committed, so the figure can be rebuilt without refitting. The
MCMC fit is the slow step, and `leuksurv_realdata.R` caches it: if
`fit_ph.rds` is present it is reused rather than refit.
Boundaries
The district outlines come from `nwengland.bnd`, read directly out of the
installed `spBayesSurv` package with `read.bnd()` from `BayesX`. A copy
also sits in `../data/` for reference. The outlines are drawn on the map
only; both spatial fits use the exact patient coordinates, not the
districts.
Other output
`leuksurv_realdata.R` also writes, into `outputs/`, an overall
Kaplan-Meier curve, the two-stage risk map and the spatial-only risk map
as standalone figures, and a three-panel comparison. Only the four-panel
version appears in the paper.
