## ============================================================
## Posterior-based fit-quality check + ACF-vs-No-ACF counterfactual
## projection, using the fitted PMCMC output from fit_pmcmc.R.
##
## The historical fit-quality/hump/trend check uses the REAL filtered
## trajectory from the production PMCMC
## (pmcmc_out$processed_chains$trajectories$state, all retained
## posterior samples) -- not a fresh re-simulation at a point estimate.
## The projection continues each posterior sample from its own filtered
## state at the end of the fit window via a stochastic ensemble, and
## prepends that same real filtered historical trajectory rather than
## re-simulating history from scratch. This matters: particle-filter
## resampling can track data patterns (e.g. a rise-then-fall in
## notifications) that a single deterministic trajectory at a point
## estimate may not show as clearly.
##
## Every zone-rate below is 1e5*sum(notes)/sum(N), a genuine per-capita
## rate -- see fit_pmcmc.R's compare-function comment for why this
## matters over a multi-year projection horizon.
##
## NOTE on closures: pmcmc_out$processed_chains$predict$transform is
## the transform closure saved inside the PMCMC output. It calls
## build_pars_flat(), a free variable resolved in whatever R global
## environment is current when the closure is invoked (it is not
## serialized with the object, since its enclosing environment is
## .GlobalEnv). common_setup.R's build_pars_flat(), in_argsrealA_ext,
## X00 and BASE_IRR_VEC must therefore match fit_pmcmc.R's exactly
## (they do, since both source the same file) for predict$transform()
## to work here.
##
## Output: fit_diagnostics.Rdata, counterfactual_projection.Rdata,
## output/counterfactual_projection.png
## ============================================================

years <- 21 # extended horizon for ACF window
source(here::here("R/utils/common_fitting_setup.R"))

cat("Loading fitted PMCMC posterior...\n")
load(here("tmpdata/pmcmc_posterior.Rdata")) # pmcmc_out, ...
cat("Loaded. n posterior samples:", nrow(pmcmc_out$processed_chains$pars), "\n")

## per-capita notification rate for zone i, summed over all age x
## HIV strata, at time index tix, across every sample/particle in Y
zone_notification_rate <- function(Y, i, tix) {
  n <- dim(Y)[2]
  notes_tot <- colSums(matrix(Y[BLASTtbmod::ln7[[i]], , tix], ncol = n))
  pop_tot <- colSums(matrix(Y[BLASTtbmod::bn7[[i]], , tix], ncol = n))
  1e5 * notes_tot / pop_tot
}

## ========= PART 1: historical fit-quality / trend / hump, full posterior =========
traj <- pmcmc_out$processed_chains$trajectories$state
cat(
  "historical trajectory dim:", paste(dim(traj), collapse = " x "),
  "(comparison on first", n_time_data, "months with real data)\n"
)


hist_check <- rbindlist(lapply(1:7, function(i) {
  rate_mat <- vapply(
    seq_len(n_time_data),
    function(k) mean(zone_notification_rate(traj, i, k)),
    numeric(1)
  )
  rate_lo <- vapply(
    seq_len(n_time_data),
    function(k) quantile(zone_notification_rate(traj, i, k), 0.05),
    numeric(1)
  )
  rate_hi <- vapply(
    seq_len(n_time_data),
    function(k) quantile(zone_notification_rate(traj, i, k), 0.95),
    numeric(1)
  )
  obs <- ndat7_orig[[note_names[i]]]
  data.table(
    zone = i, month = seq_len(n_time_data), obs = obs,
    mod_mean = rate_mat, mod_lo = rate_lo, mod_hi = rate_hi
  )
}))


bias_tab <- hist_check[
  , .(obs_mean = mean(obs), mod_mean = mean(mod_mean)),
  by = zone
]
bias_tab[, pct_bias := 100 * (mod_mean - obs_mean) / obs_mean]
cat("\n=== historical fit-quality (n=", n_time_data, "months) ===\n")
print(bias_tab)
cat("overall mean % bias:", bias_tab[, mean(pct_bias)], "\n")

zoom_months <- 24
cutoff <- n_time_data - zoom_months
trend_dt <- rbindlist(lapply(1:7, function(i) {
  d <- hist_check[zone == i & month >= cutoff]
  data.table(zone = i, data_slope = coef(lm(obs ~ month, data = d))[["month"]],
             model_slope = coef(lm(mod_mean ~ month, data = d))[["month"]])
}))
trend_dt[, match := sign(data_slope) == sign(model_slope)]
cat("\nzones matching recent-trend direction:", trend_dt[, sum(match)], "/7\n")
print(trend_dt)

shape_check <- hist_check[, {
  mid <- floor(.N / 2)
  fh <- coef(lm(mod_mean ~ month, data = .SD[1:mid]))[["month"]]
  sh <- coef(lm(mod_mean ~ month, data = .SD[(mid + 1):.N]))[["month"]]
  list(hump = fh > 0 && sh < 0)
}, by = zone]
cat("zones with hump (rise then fall):", shape_check[, sum(hump)], "/7\n")
print(shape_check)

save(
  bias_tab, trend_dt, shape_check,
  file = here("tmpdata/fit_diagnostics.Rdata")
)
cat("Saved fit-quality diagnostics to tmpdata/fit_diagnostics.Rdata\n")

## ========= PART 2: ACF-vs-No-ACF counterfactual projection =========
transform_basecase <- function(p) {
  pars <- pmcmc_out$processed_chains$predict$transform(p)
  pars <- extend_fit_window(pars, args00$sim_length - restart_step + 1)
  ITLf <- get_ITL(pars)
  for (i in 1:7) pars$ACFhaz0[i, ITLf[[i]]] <- pars$ACFhaz1[i, ITLf[[i]]] <- 0.0
  pars
}
transform_counterfactual <- function(p) {
  pars <- pmcmc_out$processed_chains$predict$transform(p)
  pars <- extend_fit_window(pars, args00$sim_length - restart_step + 1)
  ITLf <- get_ITL(pars)
  for (i in 1:7) pars$ACFhaz1[i, ITLf[[i]]] <- 0.2
  for (i in 1:7) pars$ACFhaz0[i, ITLf[[i]]] <- 0.2
  pars
}


cat(
  "\nBuilding basecase (No ACF) and counterfactual (ACF) pars for",
  nrow(pmcmc_out$processed_chains$pars),
  "posterior samples...\n"
)
basecase_pars <- lapply(
  seq_len(nrow(pmcmc_out$processed_chains$pars)),
  function(i) transform_basecase(pmcmc_out$processed_chains$pars[i, ])
)
counterfactual_pars <- lapply(
  seq_len(nrow(pmcmc_out$processed_chains$pars)),
  function(i) transform_counterfactual(pmcmc_out$processed_chains$pars[i, ])
)

predict_time <- pmcmc_out$processed_chains$predict$time
filtered_state <- pmcmc_out$processed_chains$state
output_steps <- predict_time:length(basecase_pars[[1]]$tt)
cat(
  "predict_time:",
  predict_time,
  " projecting to month",
  max(output_steps),
  "(", (max(output_steps) - predict_time) / 12,
  "years past predict_time)\n"
)

N_PROJ_PARTICLES <- 8
summarize_array <- function(Y, months, label) {
  rbindlist(lapply(1:7, function(i) {
    rate_mat <- vapply(
      seq_along(months),
      function(k) zone_notification_rate(Y, i, k),
      numeric(dim(Y)[2])
    )
    data.table(
      zone = i, month = months, scenario = label,
      mod_mean = colMeans(rate_mat),
      mod_lo = apply(rate_mat, 2, quantile, 0.025),
      mod_hi = apply(rate_mat, 2, quantile, 0.975)
    )
  }))
}


T1 <- dim(traj)[3]
cat("Summarizing shared historical trajectory (", T1, "months)...\n")
hist_summary <- summarize_array(traj, seq_len(T1), "historical")

simulate_and_summarize <- function(pars_list, label) {
  cat("Simulating", label, "forward (n_particles=", N_PROJ_PARTICLES, ")...\n")
  mod <- BLASTtbmod:::stocm$new(
    pars_list,
    predict_time,
    n_particles = N_PROJ_PARTICLES,
    seed = 1,
    deterministic = FALSE,
    pars_multi = TRUE
  )
  mod$update_state(state = filtered_state)
  res <- mod$simulate(output_steps)
  dim(res) <- c(dim(res)[1], dim(res)[2] * dim(res)[3], dim(res)[4])
  T2 <- dim(res)[3]
  out <- summarize_array(
    res[, , 2:T2, drop = FALSE],
    (T1 + 1):(T1 + T2 - 1),
    label
  )
  rm(res, mod)
  gc()
  out
}
summ0 <- simulate_and_summarize(basecase_pars, "No ACF")
summ1 <- simulate_and_summarize(counterfactual_pars, "ACF")

hist_summary_bc <- copy(hist_summary)[, scenario := "No ACF"]
hist_summary_cf <- copy(hist_summary)[, scenario := "ACF"]
cf_summary <- rbind(hist_summary_bc, hist_summary_cf, summ0, summ1)
cf_summary[, yr := start_yeare + (month - 1) / 12]

ITLf <- get_ITL(basecase_pars[[1]])
window_bounds <- rbindlist(lapply(
  seq_along(ITLf),
  function(i) {
    data.table(zone = i, lo = min(ITLf[[i]]), hi = max(ITLf[[i]]))
  }
))
window_bounds[, yr_lo := start_yeare + (lo - 1) / 12]
window_bounds[, yr_hi := start_yeare + (hi - 1) / 12]

real_dat <- data.table(
  zone = rep(1:7, each = n_time_data),
  month = rep(seq_len(n_time_data), 7),
  obs = unlist(lapply(note_names, function(nm) ndat7_orig[[nm]]))
)
real_dat[, yr := start_yeare + (month - 1) / 12]

cat("\n=== long-run (No ACF) projection slope per zone ===\n")
long_dt <- rbindlist(lapply(
  1:7,
  function(i) {
    d <- cf_summary[zone == i & scenario == "No ACF" & month > T1]
    data.table(
      zone = i,
      longrun_slope = coef(lm(mod_mean ~ month, data = d))[["month"]]
    )
  }
))
print(long_dt)


## create plot
plt <- "Accent"
p_cf <- ggplot(
  cf_summary,
  aes(
    x = yr,
    y = mod_mean, color = scenario, fill = scenario
  )
) +
  geom_rect(
    data = window_bounds,
    aes(xmin = yr_lo, xmax = yr_hi, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = NA_real_,
    lty = 2,
    col = "darkgrey",
    alpha = 0.95
  ) +
  geom_ribbon(
    aes(
      ymin = mod_lo,
      ymax = mod_hi,
      group = scenario
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(lwd = 1) +
  geom_point(
    data = real_dat,
    aes(x = yr, y = obs),
    inherit.aes = FALSE,
    col = 2, shape = 1
  ) +
  scale_fill_brewer(palette = plt) +
  scale_color_brewer(palette = plt) +
  facet_wrap(~zone) +
  xlab("Year") +
  ylab("Tuberculosis notifications per 100,000 per month") +
  theme_linedraw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.6, 0.15),
    legend.title = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 1), col = guide_legend(nrow = 1)) +
  xlim(2015, NA)


ggsave(file = here("output/Figure3.png"), w = 8, h = 7)


cat(
  "\nSaved counterfactual projection plot ",
  "to output/Figure3.png\n"
)

save(
  cf_summary, window_bounds, long_dt,
  file = here("tmpdata/counterfactual_projection.Rdata")
)
cat("Saved to tmpdata/counterfactual_projection.Rdata\n")
cat("=== DONE ===\n")


## TODO this needs revising to match the new workflow
## ============ looking at figure 4 =========
source(here("R/utils/benefit.R"))
load(here("data/pops.Rdata"))

## EB ddf are the slopes
EB[qty == "ddf" & !is.na(mid)]


## calculate slopes
slpd <- EB[qty == "ddf" & !is.na(mid),
  {
    list(slp = coef(lm(data = .SD, mid ~ t))[2])
  },
  by = patch
] # slopes
true_slopes <- slpd$slp / 1e4

## tests
halfcov <- DbenefitFromSlps(true_slopes,
  pops, # use pops as ranking
  pops, 0.5 * sum(pops),
  verbose = TRUE,
  separate = TRUE
)
halfcov
## loop over coverage
covz <- seq(from = 5e-2, to = 95e-2, by = 5e-2)
iterz <- 1:10e3
B <- list()
k <- 1
for (i in 1:length(covz)) {
  for (j in 1:length(iterz)) {
    B[[k]] <- DbenefitFromSlps(true_slopes,
      sample(7), ## pops, #use pops as ranking
      pops, covz[i] * sum(pops),
      verbose = FALSE,
      separate = TRUE
    )
    B[[k]][["coverage"]] <- covz[i]
    B[[k]][["iter"]] <- iterz[j]
    k <- k + 1
  }
}


B <- rbindlist(B)
B[, db := bnft_opt - bnft_sub]
eps <- 0.25
BS <- B[, .(
  db = mean(db),
  db.lo = mean(db) - sd(db),
  db.hi = mean(db) + sd(db)
  ## db.lo=quantile(db,eps),
  ## db.hi=quantile(db,1-eps)
), by = coverage]
BS <- rbind(
  BS,
  data.table(
    coverage = c(0.0, 1.0), db = rep(0.0, 2), db.lo = rep(0.0, 2),
    db.hi = rep(0.0, 2)
  )
)

## see utils/benefit.R for daly
BS[, DALY := db * (1 - exp(-daly$r * daly$LE)) / daly$r]
BS[, DALY.lo := db.lo * (1 - exp(-daly$r * daly$LE)) / daly$r]
BS[, DALY.hi := db.hi * (1 - exp(-daly$r * daly$LE)) / daly$r]


fig4a <- ggplot(
  BS,
  aes(
    x = coverage, y = DALY,
    ymin = DALY.lo, ymax = DALY.hi
  )
) +
  geom_ribbon(alpha = 0.3, col = NA) +
  geom_line() +
  geom_point() +
  theme_classic() +
  ggpubr::grids() +
  scale_x_continuous(label = scales::percent) +
  ylab("VOI in DALYs averted")

fig4a

## USD version
BT <- data.table(CET = seq(from = 0, to = 1e3, by = 1))
BT[, c("voi.mid", "voi.lo", "voi.hi") := BS[
  coverage == 0.5,
  .(
    CET * DALY, CET * DALY.lo,
    CET * DALY.hi
  )
]]

fig4b <- ggplot(BT, aes(x = CET, y = voi.mid, ymin = voi.lo, ymax = voi.hi)) +
  geom_ribbon(alpha = 0.3, col = NA) +
  geom_line() +
  theme_classic() +
  ggpubr::grids() +
  xlab("Cost-effectiveness threshold (USD)") +
  scale_y_continuous(label = scales::comma) +
  ylab("VOI in USD")

fig4b

fig4 <- ggarrange(fig4a, fig4b, nrow = 1, ncol = 2, labels = LETTERS[1:2])
fig4

ggsave(fig4, file = here("output/Figure4.png"), w = 10, h = 5)

