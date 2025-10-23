# BLAST-VOI
BLAST VOI analyses

# R files


## "0_Figure1.R"

As the name suggests, this code creates `Figure 1` (map and background data).

## "0_Figure2.R

As the name suggests, this code creates `Figure 2` (conceptual VOI figure).


## "1a_fitting.R" ##

This does the PMCMC fitting of the model, calculates the 7-zone return slopes, and saves out data.

## "1b_fit_analysis.R" ##

Reads in output from `fitting.R` and creates `Figure 3`.
Computes return slopes from fit data, and uses `utils/benefit.R` on these to calculate VOI; creates `Figure 4`.


## "2a_VOI_PSA.R" ##

Runs VOI PSA?

## "2b_VOI_analysis.R" ##

Creates PSA `Figure 5`

