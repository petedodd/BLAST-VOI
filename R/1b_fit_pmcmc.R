## ============================================================
## PMCMC fit of the Blantyre TB model to zone-level notification rates
## (2015-2023) plus a single cross-sectional TBI-prevalence point.
##
## theta = (beta, cdr, pDf, pDs); fixed zone-specific progression IRR
##
## Output:
## tmpdata/map_fit.Rdata
## tmpdata/pmcmc_posterior.Rdata
## output/x_pmcmc_trace.png
## ============================================================

library(coda)
library(logitnorm)

years <- 18 + 1 / 12
source(here::here("R/utils/common_fitting_setup.R"))

## ========= fit-specific data setup =========
ndat7_ext <- data.table(month = seq_len(new_sim_length))
for (nm in note_names) ndat7_ext[[nm]] <- NA_real_
ndat7_ext[["tbi_adult"]] <- NA_real_
orig_n <- args$sim_length
ndat7_ext[seq_len(orig_n), (note_names) := ndat7_orig]
ndat7_ext[new_sim_length, tbi_adult := 0.21]

pf_data7e_ext <- mcstate::particle_filter_data(
  data = ndat7_ext, time = "month", rate = 1, initial_time = 0
)

## ========= compare function (log-likelihood) =========
F_zone <- ndat7_orig[, vapply(.SD, max, numeric(1)), .SDcols = note_names]
SD_FRAC <- 0.30
TBI_SD <- 0.10

pop_zone <- BLASTtbmod::blantyre$population
pop_scale <- sqrt(mean(pop_zone) / pop_zone)

gc <- BLASTtbmod::get_cols
stvrs <- grep("U\\[|LR\\[|LL\\[|D\\[|SC\\[|Tr\\[|R\\[", gc, value = TRUE)
stvrs_kid <- grep(",1,", stvrs, value = TRUE)
stvrs_adult <- setdiff(stvrs, stvrs_kid)
stvrsn_adult <- seq_along(gc)[gc %in% stvrs_adult]
stvrsu <- grep("U\\[", gc, value = TRUE)
stvrsu_adult <- setdiff(stvrsu, grep(",1,", stvrsu, value = TRUE))
stvrsun_adult <- seq_along(gc)[gc %in% stvrsu_adult]

case_compare <- function(state, observed, pars = NULL) {
  ans <- rep(0, dim(state)[2])
  ## notifications:
  for (i in 1:7) {
    obs_i <- observed[[paste0("notifrate_", i)]]
    if (!is.na(obs_i)) {
      totnotes <- colSums(state[BLASTtbmod::ln7[[i]], , drop = FALSE])
      totpops <- colSums(state[BLASTtbmod::bn7[[i]], , drop = FALSE])
      notes_modelled <- 1e5 * totnotes / totpops
      notes_observed <- rep(obs_i, dim(state)[2])
      ans <- ans + dnorm(
        x = notes_modelled, mean = notes_observed,
        sd = SD_FRAC * F_zone[i] * pop_scale[i], log = TRUE
      )
    }
  }
  ## TBI prevalence:
  if (!is.na(observed$tbi_adult)) {
    prop_inf <- 1 - colSums(state[stvrsun_adult, , drop = FALSE]) /
      colSums(state[stvrsn_adult, , drop = FALSE])
    ans <- ans + dnorm(prop_inf, observed$tbi_adult, TBI_SD, log = TRUE)
  }
  ans
}

## ========= fixed element: initial TB prevalence =========
TBD_PREV_FIXED <- 0.0004245566

draw_popinit <- function(tbd) {
  X0 <- X00
  stat_nos <- apply(X0, c(2, 3, 4), sum)
  adult_nos <- stat_nos[, 2:3, ]
  mean_count <- tbd * adult_nos
  floor_count <- floor(mean_count)
  frac <- mean_count - floor_count
  total_diseased <- array(
    floor_count + rbinom(length(mean_count), 1, as.vector(frac)),
    dim = dim(adult_nos), dimnames = dimnames(adult_nos)
  )
  n_patch <- dim(adult_nos)[1]
  for (p in seq_len(n_patch)) {
    cell_pop <- adult_nos[p, , ]
    idx <- arrayInd(which.max(cell_pop), dim(cell_pop))
    total_diseased[p, idx[1], idx[2]] <- total_diseased[p, idx[1], idx[2]] + 1
  }
  D_count <- array(
    rbinom(
      length(total_diseased),
      size = as.vector(total_diseased),
      prob = 0.5
    ),
    dim = dim(total_diseased), dimnames = dimnames(total_diseased)
  )
  SC_count <- total_diseased - D_count
  X0["D", , 2:3, ] <- D_count
  X0["SC", , 2:3, ] <- SC_count
  X0
}

index_fn <- function(info) {
  idx <- seq_along(gc)
  names(idx) <- gc
  list(run = idx, state = idx)
}

## deterministic filter -- MAP optimization + Hessian only
filter_det <- mcstate::particle_deterministic$new(
  data = pf_data7e_ext, model = BLASTtbmod:::stocm,
  compare = case_compare, index = index_fn, n_threads = 10
)

## ========= priors =========
prior_beta <- function(x) dlnorm(x, log(0.5), .25, log = TRUE)
prior_cdr <- function(x) dbeta(x, 80, 20, log = TRUE)
pDf_logit_mu <- -2.773863
pDf_logit_sigma <- 0.1428601
prior_pDf <- function(x) {
  dlogitnorm(
    x, pDf_logit_mu, pDf_logit_sigma,
    log = TRUE
  )
}
prior_pDs <- function(x) dlnorm(x, -6.89, 0.58, log = TRUE)

MAP_SEED <- 20240718
bounds <- rbind(
  beta = c(min = 1e-6, max = 15),
  cdr = c(min = 0.2, max = 0.99),
  pDf = c(min = 1e-8, max = 1 - 1e-8),
  pDs = c(min = 1e-6, max = 0.05)
)
logit <- function(p) log(p / (1 - p))
ilogit <- function(x) 1 / (1 + exp(-x))
box_forward <- function(theta) {
  p <- (theta - bounds[, "min"]) / (bounds[, "max"] - bounds[, "min"])
  logit(pmin(pmax(p, 1e-8), 1 - 1e-8))
}
box_inverse <- function(y) {
  v <- bounds[, "min"] + ilogit(y) * (bounds[, "max"] - bounds[, "min"])
  names(v) <- rownames(bounds)
  v
}
neg_log_posterior <- function(y) {
  theta <- box_inverse(y)
  set.seed(MAP_SEED)
  X0 <- draw_popinit(TBD_PREV_FIXED)
  pars <- build_pars_flat(
    theta[["beta"]], theta[["cdr"]], theta[["pDf"]], theta[["pDs"]], X0
  )
  ll <- filter_det$run(pars = pars)
  lp <- 0 +
    prior_beta(theta[["beta"]]) +
    prior_cdr(theta[["cdr"]]) +
    prior_pDf(theta[["pDf"]]) +
    prior_pDs(theta[["pDs"]])
  val <- ll + lp
  if (!is.finite(val)) return(1e10)
  -val
}
optimize_from <- function(theta0, n_restarts = 2, maxit = 1500) {
  y <- box_forward(theta0)
  for (i in seq_len(n_restarts + 1)) {
    fit <- optim(
      y, neg_log_posterior,
      method = "Nelder-Mead", control = list(maxit = maxit, reltol = 1e-12)
    )
    y <- fit$par
  }
  list(theta = box_inverse(fit$par), value = -fit$value)
}

## ========= STEP 1: MAP estimate (multi-start) =========
cat("\n=== STEP 1: MAP estimation (theta = beta, cdr, pDf, pDs) ===\n")
starts <- list(
  start_default = c(
    beta = 0.213345546, cdr = 0.766849843, pDf = 0.056921031, pDs = 0.002309459
  ),
  start_alt1 = c(
    beta = 0.3, cdr = 0.85, pDf = 0.06, pDs = 0.005
  ),
  start_alt2 = c(
    beta = 0.150567146, cdr = 0.923109085, pDf = 0.051103972, pDs = 0.007208723
  )
)
results <- lapply(names(starts), function(nm) {
  cat("optimizing from:", nm, "...\n")
  r <- optimize_from(starts[[nm]])
  cat(
    "  log-posterior:", r$value,
    " theta:", paste(signif(r$theta, 5), collapse = ", "), "\n"
  )
  r
})
names(results) <- names(starts)
best <- results[[which.max(vapply(results, `[[`, numeric(1), "value"))]]
cat("\nBEST MAP:\n")
print(best$theta)
cat("log-posterior:", best$value, "\n")

y_opt <- box_forward(best$theta)
H <- optimHess(y_opt, neg_log_posterior)
## in case singular:
vcov_unconstrained <- tryCatch(solve(H), error = function(e) NULL)
save(
  results, best, vcov_unconstrained, F_zone, SD_FRAC, TBI_SD,
  file = here("tmpdata/map_fit.Rdata")
)

cat("Saved MAP to tmpdata/map_fit.Rdata\n")

## ========= STEP 2: Hessian-derived correlated proposal =========
map_theta <- best$theta[rownames(bounds)]
J <- (
  map_theta - bounds[, "min"]) *
  (bounds[, "max"] - map_theta) /
  (bounds[, "max"] - bounds[, "min"]
  )
Cov_theta <- diag(J) %*% vcov_unconstrained %*% diag(J)
dimnames(Cov_theta) <- list(names(map_theta), names(map_theta))
cat("\nimplied sd at MAP:\n")
print(sqrt(diag(Cov_theta)))

eig_decomp <- eigen(Cov_theta, symmetric = TRUE)
floor_ratio <- 1e-6
floor_val <- floor_ratio * max(eig_decomp$values)
n_floored <- sum(eig_decomp$values < floor_val)
if (n_floored > 0) {
  cat("flooring", n_floored, "eigenvalue(s) below", floor_val, "\n")
  adj_values <- pmax(eig_decomp$values, floor_val)
  Cov_theta <- eig_decomp$vectors %*% diag(adj_values) %*% t(eig_decomp$vectors)
  dimnames(Cov_theta) <- list(names(map_theta), names(map_theta))
}


## empirical correction for pseudo-marginal (particle-filter) noise
## inflating the effective proposal variance beyond the idealized
## 2.38^2/d asymptotic scaling
PROPOSAL_SHRINK <- 0.5
proposal_matrix <- Cov_theta * (2.38^2 / length(map_theta)) * PROPOSAL_SHRINK
dimnames(proposal_matrix) <- dimnames(Cov_theta)

## ========= STEP 3: PMCMC transform + per-particle stochastic initial state =========
make_transform_pmcmc <- function() {
  function(theta) {
    build_pars_flat(
      theta[["beta"]], theta[["cdr"]], theta[["pDf"]], theta[["pDs"]], X00
    )
  }
}
make_initial_stoch <- function(tbd_prev_fixed) {
  force(tbd_prev_fixed)
  function(info, n_particles, pars) {
    pars_list <- vector("list", n_particles)
    for (i in seq_len(n_particles)) {
      p_i <- pars
      p_i$popinit <- draw_popinit(tbd_prev_fixed)
      pars_list[[i]] <- p_i
    }
    mod <- BLASTtbmod:::stocm$new(
      pars_list,
      time = 0,
      n_particles = 1,
      seed = sample.int(1e8, 1),
      pars_multi = TRUE
    )
    s <- mod$state()
    matrix(s, nrow = dim(s)[1])
  }
}

prior_list <- list(
  beta = mcstate::pmcmc_parameter(
    "beta",
    initial = map_theta[["beta"]],
    min = bounds["beta", "min"],
    max = bounds["beta", "max"],
    prior = prior_beta
  ),
  cdr = mcstate::pmcmc_parameter(
    "cdr",
    initial = map_theta[["cdr"]],
    min = bounds["cdr", "min"],
    max = bounds["cdr", "max"],
    prior = prior_cdr
  ),
  pDf = mcstate::pmcmc_parameter(
    "pDf",
    initial = map_theta[["pDf"]],
    min = bounds["pDf", "min"],
    max = bounds["pDf", "max"],
    prior = prior_pDf
  ),
  pDs = mcstate::pmcmc_parameter(
    "pDs",
    initial = map_theta[["pDs"]],
    min = bounds["pDs", "min"],
    max = bounds["pDs", "max"],
    prior = prior_pDs
  )
)

mcmc_pars <- mcstate::pmcmc_parameters$new(
  prior_list, proposal_matrix,
  transform = make_transform_pmcmc()
)

filter <- create.particlefilter(
  pf_data7e_ext, case_compare,
  n_particles = 1000, n_threads = 10,
  index = index_fn, initial = make_initial_stoch(TBD_PREV_FIXED)
)


report <- function(label, pc) {
  mcmc1 <- coda::as.mcmc(cbind(pc$probabilities, pc$pars))
  cat("\n===", label, "===\n")
  cat("n retained:", nrow(pc$pars), "\n")
  cat("acceptance rate:", mean(1 - coda::rejectionRate(mcmc1)), "\n")
  cat("ESS:\n")
  print(round(coda::effectiveSize(mcmc1), 1))
  cat("posterior mean:\n")
  print(colMeans(pc$pars))
  cat("correlation matrix:\n")
  print(round(cor(pc$pars), 3))
  mcmc1
}


## ========= STEP 4: PMCMC run =========
N_STEPS <- 2500
N_BURNIN <- 500
cat("\n=== STEP 4: PMCMC (", N_STEPS, "steps, n_particles=1000) ===\n")
t0 <- Sys.time()
pmcmc_out <- run.pmcmc(
  particle.filter = filter,
  parms = in_argsrealA_ext,
  n.steps = N_STEPS,
  n.burnin = N_BURNIN,
  n.chains = 1,
  n.threads = 10,
  n.epochs = 1,
  n.workers = 1,
  mcmc_pars = mcmc_pars,
  save_restart = in_argsrealA_ext$sim_length,
  returnall = TRUE
)

cat("elapsed:", format(Sys.time() - t0), "\n")
pc <- pmcmc_out$processed_chains
mcmc1 <- report("PMCMC posterior", pc)

## MCMC traceplot
p <- bayesplot::mcmc_trace(mcmc1[, -(1:3)])
ggsave(p, file = here("output/x_pmcmc_trace.png"), w = 10, h = 8)

## MCMC correlation plot
png(here("output/x_pmcmc_corplot.png"), w = 800, h = 800)
corplot::corplot( as.matrix(mcmc1[, -(1:3)]), points = TRUE)
dev.off()

save(
  pmcmc_out, proposal_matrix, Cov_theta, map_theta, F_zone, SD_FRAC, TBI_SD,
  file = here("tmpdata/pmcmc_posterior.Rdata")
)

cat("Saved PMCMC posterior to tmpdata/pmcmc_posterior.Rdata\n")
cat("=== DONE ===\n")
