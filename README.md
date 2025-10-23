# BLAST-VOI
BLAST VOI analyses

# R files

Looking at contents the other way around.

## "0_Figure1.R"

As the name suggests, this code creates `Figure 1` (map and background data).

## "0_Figure2.R

As the name suggests, this code creates `Figure 2` (conceptual VOI figure).


## "1a_fitting.R" ##

This does the PMCMC fitting of the model, calculates the 7-zone return slopes, and saves out data.

### "1b_fit_analysis.R" ###

Reads in output from `fitting.R` and creates `Figure 3`.
Computes return slopes from fit data, and uses `utils/benefit.R` on these to calculate VOI; creates `Figure 4`.


## "2a_VOI_PSA.R" ##

Runs VOI PSA?

## "2b_VOI_analysis.R" ##

Creates PSA `Figure 5`



# Figure notes (old)

## Figure 1
*data & setting*

"Figure1.png"

Done by:
"Figure1.R"

## Figure 2
*VOI concept*

"VOI_concept.png"

Done by:
"VOIfigure.R"

## Figure 3
*fit*

"fit_fig_simple.png"

Done by:
"fitting.R"

## Figure 4
*Blantyre VOI*

"fig4.png"

Done by:
"matchedfitplots.R"

## Figure 5
*VOI PSA*

"VOI2_psa_fig.png"

Done by (I think):
"VOI_N_aggTask3.R"
"VOIfigure.R"

