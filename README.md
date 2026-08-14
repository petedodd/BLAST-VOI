# BLAST-VOI #

BLAST VOI analyses

# Analyses

Every `output/Figure{n}.png` is saved alongside a matching `output/Figure{n}.pdf`.

## 0 = background figures ##

### "0_Figure1.R" ###

As the name suggests, this code creates `Figure 1` (map and background data).
Outputs: `output/Figure1.png`, `output/Figure1.pdf`

### "0_Figure2.R" ###

As the name suggests, this code creates `Figure 2` (conceptual VOI figure).
Outputs: `output/Figure2.png`, `output/Figure2.pdf`

### "0_flux_approx.R" ###

Explores and justifies the flux approximations used elsewhere in the pipeline.
Outputs: `output/x_fluxapp.png`, `output/x_fluxapp_model.png`, `output/x_fluxapp_model_opt.png`, `output/x_fluxapp_note_v_FOI.png`

## 1 = Blantyre analyses ##

### "1a_hiv_art_checks.R" ###

Sources `utils/common_fitting_setup.R` and checks HIV/ART/demographic model
dynamics against (scaled) national survey data.
Outputs: `output/x_hiv_art_dynamics.png`, `output/x_hiv_zone_prevalence.png`, `output/x_hiv_in_tb.png`

### "1b_fit_pmcmc.R" ###

PMCMC fit of the Blantyre TB model to zone-level notification rates
(2015-2023) plus a single cross-sectional TBI-prevalence point.
Outputs: `tmpdata/map_fit.Rdata`, `tmpdata/pmcmc_posterior.Rdata`, `output/x_pmcmc_trace.png`, `output/x_pmcmc_corplot.png`

### "1c_posterior_diagnostics.R" ###

Posterior-based fit-quality check and ACF-vs-No-ACF counterfactual
projection, using the fitted PMCMC output from `1b_fit_pmcmc.R`; creates
`Figure 3` and `Figure 4`.
Outputs: `tmpdata/fit_diagnostics.Rdata`, `tmpdata/simulation_summaries.Rdata`, `tmpdata/counterfactual_projection.Rdata`, `output/Figure3.png`, `output/Figure3.pdf`, `output/Figure4.png`, `output/Figure4.pdf`

## 2 = PSA of VOI ##

### "2a_PSA_sample.R" ###

Creates the PSA parameter sample.
Outputs: `tmpdata/PARGS.Rdata`

### "2b_run_PSAs.sh" ###

Runs `utils/PSA_single_N.R` in parallel over different numbers of zones
(5, 10, 15, ..., 50).
Outputs: `tmpdata/R{n}.Rdata`, `tmpdata/P{n}.Rdata` for each n

### "2c_PSA_analysis.R" ###

Reads in PSA run results and calculates VOI metrics.
Outputs: `tmpdata/PB.Rdata`, `output/VOI_psa_prccN.csv`, `output/VOI_psa_prccNR.csv`

### "2d_VOI_figure.R" ###

Creates PSA `Figure 5`.
Outputs: `output/Figure5.png`, `output/Figure5.pdf`

# Suggested run order #

```
Rscript R/0_Figure1.R
Rscript R/0_Figure2.R
Rscript R/0_flux_approx.R          # optional, exploratory

Rscript R/1a_hiv_art_checks.R      # optional, diagnostic
Rscript R/1b_fit_pmcmc.R
Rscript R/1c_posterior_diagnostics.R

Rscript R/2a_PSA_sample.R
bash R/2b_run_PSAs.sh
Rscript R/2c_PSA_analysis.R
Rscript R/2d_VOI_figure.R
```

# Dependencies #

R version 4.4+ recommended. The following packages were used (versions as
currently installed in the development environment):

| Package | Version |
|---------|---------|
| bayesplot | 1.15.0 |
| BLASTtbmod | 0.0.1 |
| coda | 0.19-4.1 |
| corplot | 0.0.0.9000 |
| data.table | 1.18.2.1 |
| epiR | 2.0.96 |
| expm | 1.0-0 |
| ggplot2 | 4.0.3 |
| ggpubr | 0.6.3 |
| ggrepel | 0.9.8 |
| ggthemes | 5.2.0 |
| glue | 1.8.1 |
| here | 1.0.2 |
| lhs | 1.3.0 |
| logitnorm | 0.8.39 |
| lubridate | 1.9.5 |
| mcstate | 0.9.22 |
| paletteer | 1.7.0 |
| patchwork | 1.3.2 |
| scales | 1.4.0 |
| sf | 1.1-0 |
| tictoc | 1.2.1 |

`ggthemes` is not loaded directly but is required by `paletteer` for the
`"ggthemes::calc"` palette used in `utils/flux_utils.R`.

The packages `BLASTtbmod` and `corplot` are not on CRAN and are available from:
- https://github.com/petedodd/BLASTtbmod
- https://github.com/petedodd/corplot
