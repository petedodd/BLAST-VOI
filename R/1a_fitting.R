## generating data for plot
library(BLASTtbmod)
library(data.table)
library(ggplot2)
library(here)

data(in_args) # example
## NOTE this data is for 72 months from 2015, so 2015-2020
data(pf_data7)

## define other parameters
n_particles <- 10
start_year <- 2010
years <- 18 + 1 / 12 #duration of runs
## this is to accommodate data and ACF experiment phases

## ========= HIV background & set-up

## how much higher is Blantyre HIV prevalence than national (data)?
mwi_vs_blantyre_PR <-
  mean(BLASTtbmod::blantyre$hivpre) / # 2015
    hivp_mwi[variable == "HIVpc" & Period == 2015, value]
## how much higher is Blantyre HIV inc than national: use prev data
mwi_vs_blantyre <- mwi_vs_blantyre_PR


## change parms to match HIV data
args <- get.parms(
  start_year = start_year, years = years,
  hivfac = mwi_vs_blantyre, # taken from data
  hivdecline = 0, hiv_init_override = 0.21,
  ART_haz = 0.18, ART_init_override = 1e-1,
  hiv_checking = TRUE
)
hirr <- 40
args$Hirr <- c(1, hirr, hirr * 0.43)
args$ari0 <- 9e-2
args$initD[, 2:3] <- 500e-5
args$beta <- 2
args$cdr <- 0.8
## ACF: doing this makes B
## this sets when ACF (measurements) are on
ne <- args$sim_length
ITL <- list(
  (ne - 1 * 12):ne,
  (ne - 2 * 12):(ne - 1 * 12),
  (ne - 3 * 12):(ne - 2 * 12),
  (ne - 4 * 12):(ne - 3 * 12),
  (ne - 5 * 12):(ne - 4 * 12),
  (ne - 6 * 12):(ne - 5 * 12),
  (ne - 7 * 12):(ne - 6 * 12)
)
for (i in 1:7) args$ACFhaz0[i, ITL[[i]]] <- args$ACFhaz1[i, ITL[[i]]] <- 0.2
## break points in ITL as data for plots
brks <- sort(c(unlist(lapply(ITL, min)), unlist(lapply(ITL, max))))
brk_yrs <- data.table(t = brks, yr = start_year + brks / 12)

## run fwd simulation & test (un-calibrated)
test <- run.model(args, args$tt, n.particles = 200)

## check un-calibrated notifications & ACF timing
gp <- plot_compare_noterate_agrgt(test,
  realdata = TRUE,
  start_year = start_year
)
gp + geom_vline(
  data = brk_yrs,
  aes(xintercept = yr),
  linetype = "dashed", col = "grey"
)

## ----------- other checks
## --- inspect demographic outputs
plot_compare_demog(test, start_year = 2015, by_comp = "age")
## NOTE we don't expect perfect agreement here because
## we data are scaled national demographic change
## and there is much higher HIV prevalence in Blantyre


## --- HIV comparisons
gp <- plot_HIV_dynamic(test,
  start_year = start_year,
  by_patch = FALSE
)
## comparison data: need to scale national
data(hivp_mwi)
hivp_mwi[, step := (Period - start_year) * 12 + 1]
xdta <- hivp_mwi[Period >= 2015]
## NOTE scales to match difference seen in data
xdta[
  variable == "HIVpc",
  c("value", "lo", "hi") := .(
    value * mwi_vs_blantyre_PR,
    lo * mwi_vs_blantyre_PR,
    hi * mwi_vs_blantyre_PR
  )
]
gp <- gp +
  geom_pointrange(data = xdta, aes(ymin = lo, ymax = hi), shape = 1) ## +
  ## xlim(c(2015, 2025))
gp

ggsave(gp, filename = here("output/x_hivart.png"), w = 7, h = 5)

## --- zone-wise comparison
hivpd <- BLASTtbmod::blantyre$hivpre
hivpd <- data.table(
  zone = paste0("Zone ", 1:7),
  step = hivp_mwi[Period == 2015 & variable == "HIVpc", step],
  variable = "HIVpc", value = hivpd
)
gp <- plot_HIV_dynamic(test,
  start_year = start_year,
  show_ART = FALSE
)
gp <- gp +
  geom_point(data = hivpd, pch = 1, size = 2, stroke = 2) +
  xlim(c(2015, 2025)) + ylim(c(0, NA))
gp

ggsave(gp, filename = here("output/x_hivpatch.png"), w = 7, h = 7)

## HIV in TB
gp <- plot_HIV_in_TB(test, start_year = start_year) + xlim(c(2015, 2021))
gp

ggsave(gp, filename = here("output/x_hivintb.png"), w = 7, h = 5)

## ========= INFERENCE

## update parameters for restart in 2015
args0 <- args #for safe-keeping
args <- restart_parms(args0, 60, test)

## === reincorporating old version from 2015
start_year <- 2015
years <- 13 + 1 / 6
## CHECK
length(as.double(seq(0, 12 * years)))
args$sim_length


## update parameters for restart in 2014
args0 <- args # for safe-keeping
args <- restart_parms(args0, 60 - 6, test)


## === reincorporating old version from 2015
start_year <- 2015 - 1 / 2
years <- 13 + 1 / 6 + 1 / 2
## CHECK
length(as.double(seq(0, 12 * years)))
args$sim_length


## Redo ACF
## ACF: doing this makes B
ne <- args$sim_length
ITL <- list(
  (ne - 1 * 12):ne,
  (ne - 2 * 12):(ne - 1 * 12),
  (ne - 3 * 12):(ne - 2 * 12),
  (ne - 4 * 12):(ne - 3 * 12),
  (ne - 5 * 12):(ne - 4 * 12),
  (ne - 6 * 12):(ne - 5 * 12),
  (ne - 7 * 12):(ne - 6 * 12)
)
for (i in 1:7) args$ACFhaz0[i, ITL[[i]]] <- args$ACFhaz1[i, ITL[[i]]] <- 0.2

## ## fwd simulation & re-test
test2 <- run.model(args, args$tt, n.particles = 200)
## plot_compare_noterate_agrgt(test2, realdata = TRUE)



## create common compare function and filter
S <- 10
case_compare7v <- function(state, observed, pars = NULL) {
  ## NOTE internal: state is not like test
  ## state is: n states x n particles
  ans <- rep(0, dim(state)[2]) # number of particles
  ## loop over zones
  for (i in 1:7) {
    totnotes <- colSums(
      state[BLASTtbmod::ln7[[i]], , drop = TRUE] *
        state[BLASTtbmod::bn7[[i]], , drop = TRUE]
    )
    totpops <- colSums(state[BLASTtbmod::bn7[[i]], , drop = TRUE])
    notes_modelled <- totnotes / totpops
    notes_observed <- observed[[paste0("notifrate_", i)]]
    if (is.null(observed$SV)) {
      SV <- S # pull from top environment
    } else {
      SV <- observed$SV # NOTE
    }
    ## measurement likelihood contribution
    ans <- ans + dnorm(
      x = notes_modelled, mean = notes_observed,
      sd = SV, log = TRUE
    )
  }
  ans
}

## create PF
filter <- create.particlefilter(
  pf_data7,
  case_compare7v,
  n_particles = 100,
  n_threads = 4
)
## NOTE or see below



## --- d0 inference
## args$initD
curve(dlnorm(x, log(150 / 1e5), 0.99), n = 1e3, from = 0, to = 0.01)
D0pr <- function(x) dlnorm(x, log(1e0 / 1e5), 1)

## scale prior
D00 <- args$initD
curve(dlnorm(x, -0.5, 1), n = 1e3, from = 0, to = 5)
SDpr <- function(x) dlnorm(x, -1 / 2, 1)


## transform
## making an option outside for greater flexibility
make_transform <- function(ARGS) {
  function(theta) {
    c(
      ARGS,
      list(
        beta = unname(theta[1]),
        ari0 = unname(theta[2]),
        initD = D00 * theta[3:9]
      )
    )
  }
}


## common inference priors
## beta
betamn <- 1.678
betasg <- 0.371
beta <- mcstate::pmcmc_parameter("beta",
  initial = qlnorm(0.5, betamn, betasg),
  min = 1e-6, max = 15,
  prior = function(x) dlnorm(x, betamn, betasg)
)

## ari0
ari0 <- mcstate::pmcmc_parameter("ari0",
  initial = qlnorm(0.5, log(0.05), 0.75),
  min = 1e-6, max = 0.5,
  prior = function(x) dlnorm(x, log(0.05), 0.75)
)

proposal_matrix <- diag(c(
  0.05, # beta
  1e-5, # ari0
  rep(0.05, 7) # SDpr
))



## A: 'TAU' inference
for (i in 1:7) args$ACFhaz0[i, ITL[[i]]] <- args$ACFhaz1[i, ITL[[i]]] <- 0 # zero again
in_argsrealA <- args
in_argsrealA$beta <- in_argsrealA$ari0 <- in_argsrealA$initD <- NULL

## proposal_matrix <- A$proposal_matrix
## save(proposal_matrix, file = here("tmpdata/proposal_matrix.Rdata"))

init_scale <- 0.5
prior_list <- list(
  beta = beta, ari0 = ari0,
  d1 = mcstate::pmcmc_parameter("d1",
    initial = init_scale,
    min = 1e-6, max = 10, prior = SDpr
  ),
  d2 = mcstate::pmcmc_parameter("d2",
    initial = init_scale,
    min = 1e-6, max = 10, prior = SDpr
  ),
  d3 = mcstate::pmcmc_parameter("d3",
    initial = init_scale,
    min = 1e-6, max = 10, prior = SDpr
  ),
  d4 = mcstate::pmcmc_parameter("d4",
    initial = init_scale,
    min = 1e-6, max = 10, prior = SDpr
  ),
  d5 = mcstate::pmcmc_parameter("d5",
    initial = init_scale,
    min = 1e-6, max = 10, prior = SDpr
  ),
  d6 = mcstate::pmcmc_parameter("d6",
    initial = init_scale,
    min = 1e-6, max = 10, prior = SDpr
  ),
  d7 = mcstate::pmcmc_parameter("d7",
    initial = init_scale,
    min = 1e-6, max = 10, prior = SDpr
  )
)

## as list
mcmc_pars <- mcstate::pmcmc_parameters$new(
  prior_list,
  proposal_matrix,
  transform = make_transform(in_argsrealA)
)

## check
mcmc_pars$initial()
mcmc_pars$model(mcmc_pars$initial()) #looks OK


## A: run inference
pmcmc_out <- run.pmcmc(
  particle.filter = filter,
  parms = in_argsrealA,
  n.steps = 500, n.burnin = 250, n.chains = 1,
  n.threads = 4, n.epochs = 1,
  mcmc_pars = mcmc_pars,
  save_restart = 72, returnall = TRUE
)


mcmc1 <- coda::as.mcmc(cbind(
  pmcmc_out$processed_chains$probabilities,
  pmcmc_out$processed_chains$pars
))

summary(mcmc1)
coda::effectiveSize(mcmc1)
coda::rejectionRate(mcmc1)
## plot(mcmc1)

## check
## plot_compare_noterate_agrgt( pmcmc_out$processed_chains$trajectories$state )


## approach to running counterfactual
transform_counterfactual <- function(p) {
  pars <- pmcmc_out$pmcmc_run$predict$transform(p)
  for (i in 1:7) pars$ACFhaz0[i, ITL[[i]]] <- pars$ACFhaz1[i, ITL[[i]]] <- 0.2
  pars
}
transform_basecase <- function(p) {
  pars <- pmcmc_out$pmcmc_run$predict$transform(p)
  for (i in 1:7) pars$ACFhaz0[i, ITL[[i]]] <- pars$ACFhaz1[i, ITL[[i]]] <- 0.0
  pars
}


## loop over parameters for basecase
basecase_pars <- lapply(
  seq_len(nrow(pmcmc_out$pmcmc_run$pars)),
  function(i) transform_basecase(pmcmc_out$pmcmc_run$pars[i, ])
)

## loop over parameters for counterfactual
counterfactual_pars <- lapply(
  seq_len(nrow(pmcmc_out$pmcmc_run$pars)),
  function(i) transform_counterfactual(pmcmc_out$pmcmc_run$pars[i, ])
)

## basecase model
bcmod <- BLASTtbmod:::stocm$new(basecase_pars, 0, ## initial_time,
  n_particles = 1, seed = 1,
  deterministic = TRUE, pars_multi = TRUE
)

## counterfactual model
cfmod <- BLASTtbmod:::stocm$new(counterfactual_pars, 0, ## initial_time,
  n_particles = 1, seed = 1,
  deterministic = TRUE, pars_multi = TRUE
)


## simulate models
(maxt <- dim(args$ACFhaz0)[2])
output_steps <- seq(0, maxt)

## res will have dimensions (states x particles x samples x time)
res0 <- bcmod$simulate(output_steps)
res1 <- cfmod$simulate(output_steps)
str(res0)
ress0 <- res0[, 1, , ]
str(ress0)
ress1 <- res1[, 1, , ]
str(ress1)

fn <- here("tmpdata")
if (!file.exists(fn)) dir.create(fn)

save(ress0, file = here("tmpdata/ress0.Rdata"))
save(ress1, file = here("tmpdata/ress1.Rdata"))

## ===========================================


## === YET ANOTHER GO

## adjusting filter data to match simulation times
## years within sim when there are actualy data
time_cols <- c("month_start", "month_end", "time_start", "time_end")
note_cols <- setdiff(colnames(pf_data7), time_cols)

## trying a version that is long
mn_notes <- colMeans(pf_data7[, note_cols])
pf_data7_adj <- pf_data7[rep(1, nrow(pf_data7) + 6), ] # template
pf_data7_adj[, note_cols] <- matrix(
  mn_notes,
  nrow = nrow(pf_data7) + 6, ncol = 7,
  byrow = TRUE
)

seen <- (2015 - start_year) * 12 # offset
seen <- (seen + 1):(seen + 72)
pf_data7_adj[seen, note_cols] <- pf_data7[, note_cols]
pf_data7_adj <- cbind(pf_data7_adj, SV = 10)
pf_data7_adj[seen, "SV"] <- 10 #true observations
pf_data7_adj <- pf_data7_adj[, c(note_cols, "SV")]
pf_data7_adj <- cbind(pf_data7_adj, month = seq_len(nrow(pf_data7_adj)))

pf_data7_adj <- mcstate::particle_filter_data(
  pf_data7_adj,
  time = "month",
  rate = 1,
  initial_time = 0
)


head(pf_data7_adj)
tail(pf_data7_adj)

## create PF
filter <- create.particlefilter(
  pf_data7_adj,
  case_compare7v,
  n_particles = 100,
  n_threads = 4
)
