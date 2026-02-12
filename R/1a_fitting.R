## generating data for plot
library(BLASTtbmod)
library(data.table)
library(ggplot2)
library(here)

data(in_args) #example
data(pf_data7)

## define other parameters
n_particles <- 10
time <- 1 #NOTE not used TODO check
start_year <- 2015
years <- 6

## create common compare function and filter
S <- 10
case_compare7v <- function(state, observed, pars = NULL) {
  ans <- rep(0, dim(state)[2])
  for (i in 1:7) {
    totnotes <- colSums(state[BLASTtbmod::ln7[[i]], , drop = TRUE] *
      state[BLASTtbmod::bn7[[i]], , drop = TRUE])
    totpops <- colSums(state[BLASTtbmod::bn7[[i]], , drop = TRUE])
    notes_modelled <- totnotes / totpops
    notes_observed <- observed[[paste0("notifrate_", i)]]
    ans <- ans + dnorm(
      x = notes_modelled, mean = notes_observed,
      sd = S, log = TRUE
    )
  }
  ans
}


filter <- create.particlefilter(
  pf_data7,
  case_compare7v,
  n_particles = 100,
  n_threads = 4
)


## change parms
args <- get.parms(start_year = start_year, years = years + 7, hivfac = 2)
args$ari0 <- 5e-2
args$initD[, 2:3] <- 150e-5
args$beta <- 5
args$cdr <- 0.8
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
## fwd simulation & test
test <- run.model(args, args$tt, n.particles = 200)


## check
plot_compare_noterate_agrgt(test, realdata = FALSE)


## ----------- other checks
## --- inspect demographic outputs
plot_compare_demog(test, by_comp = "age")
## NOTE we don't expect perfect agreement here because
## we data are scaled national demographic change
## and there is much higher HIV prevalence in Blantyre


## --- HIV comparisons
gp <- plot_HIV_dynamic(test, start_year = 2015, by_patch = FALSE)
## comparison data: need to scale national
data(hivp_mwi)
hivp_mwi[, step := (Period - 2015) * 12 + 1]
xdta <- hivp_mwi[Period >= 2015]
## NOTE scales to match initial:
fac <- xdta[step == min(step) & variable == "HIVpc", value]
fac <- gp@data[step == min(step) & variable == "HIVpc", value] / fac
xdta[
  variable == "HIVpc",
  c("value", "lo", "hi") := .(value * fac, lo * fac, hi * fac)
]
gp <- gp +
  geom_pointrange(data = xdta, aes(ymin = lo, ymax = hi), shape = 1) +
  xlim(c(2015, 2025))
gp

ggsave(gp, filename = here("output/x_hivart.png"), w = 7, h = 5)

## --- zone-wise comparison
hivpd <- BLASTtbmod::blantyre$hivpre
hivpd <- data.table(
  patch = paste0("Patch ", 1:7), step = 1, variable = "HIVpc", value = hivpd
)

gp <- plot_HIV_dynamic(test, start_year = 2015, show_ART = FALSE)
gp <- gp + geom_point(data = hivpd, pch = 1, size = 2, stroke = 2)
gp

ggsave(gp, filename = here("output/x_hivpatch.png"), w = 7, h = 7)

## HIV in TB
gp <- plot_HIV_in_TB(test, start_year = 2015) + xlim(c(2015, 2021))
gp

ggsave(gp, filename = here("output/x_hivintb.png"), w = 7, h = 5)

## --- D0 inference
## args$initD
curve(dlnorm(x, log(150 / 1e5), 0.99), n = 1e3, from = 0, to = 0.01)
D0pr <- function(x) dlnorm(x, log(1e0 / 1e5), 1)


## transform
## making an option outside for greater flexibility
make_transform <- function(ARGS) {
  function(theta) {
    c(
      ARGS,
      list(
        beta = unname(theta[1]),
        ari0 = unname(theta[2]),
        initD = cbind(
          rep(0, 7),
          unname(theta[3:9]),
          unname(theta[3:9])
        )
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
  rep(1e-9, 7) # TODO D0
))


## A: 'TAU' inference
for (i in 1:7) args$ACFhaz0[i, ITL[[i]]] <- args$ACFhaz1[i, ITL[[i]]] <- 0 # zero again
in_argsrealA <- args
in_argsrealA$beta <- in_argsrealA$ari0 <- in_argsrealA$initD <- NULL

## proposal_matrix <- A$proposal_matrix
## save(proposal_matrix, file = here("tmpdata/proposal_matrix.Rdata"))

prior_list <- list(
  beta = beta, ari0 = ari0,
  d1 = mcstate::pmcmc_parameter("d1",
    initial = 2 * 150e-5,
    min = 1e-6, max = 2e-2, prior = D0pr
  ),
  d2 = mcstate::pmcmc_parameter("d2",
    initial = 150e-5,
    min = 1e-6, max = 2e-2, prior = D0pr
  ),
  d3 = mcstate::pmcmc_parameter("d3",
    initial = 4 * 150e-5,
    min = 1e-6, max = 2e-2, prior = D0pr
  ),
  d4 = mcstate::pmcmc_parameter("d4",
    initial = 4 * 150e-5,
    min = 1e-6, max = 2e-2, prior = D0pr
  ),
  d5 = mcstate::pmcmc_parameter("d5",
    initial = 150e-5,
    min = 1e-6, max = 2e-2, prior = D0pr
  ),
  d6 = mcstate::pmcmc_parameter("d6",
    initial = 150e-5,
    min = 1e-6, max = 2e-2, prior = D0pr
  ),
  d7 = mcstate::pmcmc_parameter("d7",
    initial = 150e-5,
    min = 1e-6, max = 2e-2, prior = D0pr
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


## A <- run.pmcmc(
##   particle.filter = filter,
##   parms = in_argsrealA,
##   prior.list = list(
##     beta = beta,
##     ari0 = ari0
##   ),
##   initial.proposal.matrix = proposal_matrix,
##   n.steps = 500, n.burnin = 250, n.chains = 1,
##   n.threads = 4, n.epochs = 1,
##   save_restart = 72, returnall = TRUE
## )

## processed_chains <- mcstate::pmcmc_thin(A)
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
maxt <- dim(args$ACFhaz0)[2]
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
