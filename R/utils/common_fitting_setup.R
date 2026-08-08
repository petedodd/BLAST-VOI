## ============================================================
## Shared setup for the TB notification-fitting pipeline
##
## Callers must set `years` (the total simulated horizon, in years)
## before sourcing this file:
## 1b_fit_pmcmc.R only needs to reach the TBI survey point (18+1/12 years);
## 1c_posterior_diagnostics.R needs 21 years for the ACF projection window
## ============================================================

library(BLASTtbmod)
library(data.table)
library(here)
library(ggplot2)

data(in_args)
data(pf_data7)

mwi_vs_blantyre <-
  mean(BLASTtbmod::blantyre$hivpre) /
    hivp_mwi[variable == "HIVpc" & Period == 2015, value]

get_ITL <- function(pars) {
  ne <- pars$sim_length
  list(
    (ne - 1 * 12):ne, (ne - 2 * 12):(ne - 1 * 12), (ne - 3 * 12):(ne - 2 * 12),
    (ne - 4 * 12):(ne - 3 * 12), (ne - 5 * 12):(ne - 4 * 12),
    (ne - 6 * 12):(ne - 5 * 12), (ne - 7 * 12):(ne - 6 * 12)
  )
}

start_year <- 2010
args00 <- args <- get.parms(
  start_year = start_year, years = years, ari0 = 1e-2,
  Dinit = cbind(rep(0, 7), rep(150e-5, 7), rep(150e-5, 7)),
  hivfac = mwi_vs_blantyre, hivdecline = 0, hiv_init_override = 0.21,
  ART_haz = 0.18, ART_init_override = 1e-1, hiv_checking = FALSE, debug = FALSE
)
hirr <- 25
args$Hirr <- c(1, hirr, hirr * 0.43)
args$beta <- 1.25
args$cdr <- 0.8
ITL <- get_ITL(args)
for (i in 1:7) args$ACFhaz0[i, ITL[[i]]] <- args$ACFhaz1[i, ITL[[i]]] <- 0.2
args0 <- args
out <- run.model(args, args$tt, n.particles = 200)

before <- 6
args <- restart_parms(args0, 60 - before, out, 72 + before)

restart_step <- 60 - before
tbi_survey_step_args00 <- 162
new_sim_length <- tbi_survey_step_args00 - restart_step + 1
extend_fit_window <- function(p, sim_length_new) {
  keep <- restart_step:(restart_step + sim_length_new - 1)
  newparms <- p
  newparms$tt <- args00$tt[keep] - args00$tt[keep][1]
  newparms$sim_length <- length(keep)
  newparms$births_int <- args00$births_int[keep]
  newparms$HIV_int <- args00$HIV_int[keep]
  newparms$ART_int <- args00$ART_int[keep]
  newparms$mu_noHIV_int <- args00$mu_noHIV_int[, keep]
  newparms$mu_HIV_int <- args00$mu_HIV_int[, keep]
  newparms$mu_ART_int <- args00$mu_ART_int[, keep]
  newparms$m_in_int <- args00$m_in_int[, keep]
  newparms$ACFhaz0 <- args00$ACFhaz0[, keep]
  newparms$ACFhaz1 <- args00$ACFhaz1[, keep]
  newparms
}

note_names <- paste0("notifrate_", 1:7)
ndat7_orig <- as.data.table(pf_data7)[, ..note_names]
ndat7_orig <- rbind(ndat7_orig[rep(1, before)], ndat7_orig)
n_time_data <- nrow(ndat7_orig) # months with real notification data
start_yeare <- 2015 - before / 12

## data for fitting, including MAP-derived
## zone-specific progression IRR (incidence rate ratio)
X00 <- args$popinit
BASE_IRR_VEC <- rep(1, 7)
BASE_IRR_VEC[c(3, 4, 5, 6)] <- c(
  0.650190565, 3.033807548, 1.150418365, 1.809569580
)

in_argsrealA <- args
in_argsrealA$beta <- NULL
in_argsrealA$cdr <- NULL
in_argsrealA$popinit <- NULL
in_argsrealA$pDf <- NULL
in_argsrealA$pDs <- NULL
in_argsrealA$IRR <- NULL
in_argsrealA_ext <- extend_fit_window(in_argsrealA, new_sim_length)

build_pars_flat <- function(beta, cdr, pDf, pDs, popinit) {
  c(in_argsrealA_ext, list(
    beta = unname(beta), cdr = unname(cdr), popinit = popinit,
    pDf = unname(pDf), pDs = unname(pDs), IRR = BASE_IRR_VEC
  ))
}

