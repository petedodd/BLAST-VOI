# BLAST-VOI #

BLAST VOI analyses

# Analyses

## 0 = background figures ##


### "0_Figure1.R" ###

As the name suggests, this code creates `Figure 1` (map and background data)

### "0_Figure2.R ###

As the name suggests, this code creates `Figure 2` (conceptual VOI figure)


## 1 = Blantyre analyses ##

### "1a_fitting.R" ###

This does the PMCMC fitting of the model, calculates the 7-zone return slopes, and saves out data

### "1b_fit_analysis.R" ###

Reads in output from `fitting.R` and creates `Figure 3`
Computes return slopes from fit data, and uses `utils/benefit.R` on these to calculate VOI; creates `Figure 4`

## 2 = PSA of VOI ##

### "2a_PSA_sample.R" ###

Creates PSA parameter sample

### "2b_run_PSAs.sh" ###

Runs PSAs over different numbers of zones using `utils/PSA_single_N.R`

### "2c_PSA_analysis.R" ###

Reads in PSA run results and calculates VOI metrics

### "2d_VOI_figure.R" ###

Creates PSA `Figure 5`

# Dependencies #

TODO
