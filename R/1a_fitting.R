## generating data for plot
library(BLASTtbmod)
library(data.table)
library(ggplot2)
library(here)

data(in_args) # example
## NOTE this data is for 72 months from 2015, so 2015-2020
data(pf_data7)

## ========= HIV background & set-up

## how much higher is Blantyre HIV prevalence than national (data)?
mwi_vs_blantyre_PR <-
  mean(BLASTtbmod::blantyre$hivpre) / # 2015
    hivp_mwi[variable == "HIVpc" & Period == 2015, value]

## how much higher is Blantyre HIV inc than national: use prev data
mwi_vs_blantyre <- mwi_vs_blantyre_PR

## ACF helper
get_ITL <- function(pars) {
  ne <- pars$sim_length
  list(
    (ne - 1 * 12):ne,
    (ne - 2 * 12):(ne - 1 * 12),
    (ne - 3 * 12):(ne - 2 * 12),
    (ne - 4 * 12):(ne - 3 * 12),
    (ne - 5 * 12):(ne - 4 * 12),
    (ne - 6 * 12):(ne - 5 * 12),
    (ne - 7 * 12):(ne - 6 * 12)
  )
}

## define other parameters
n_particles <- 10
start_year <- 2010
years <- 18 + 1 / 12 #duration of runs
## this is to accommodate data and ACF experiment phases
## change parms to match HIV data
args00 <- args <- get.parms(
  start_year = start_year, years = years,
  ari0 = 1e-2,
  Dinit = cbind(rep(0, 7), rep(150e-5, 7), rep(150e-5, 7)),
  hivfac = mwi_vs_blantyre, # taken from data
  hivdecline = 0, hiv_init_override = 0.21,
  ART_haz = 0.18, ART_init_override = 1e-1,
  hiv_checking = FALSE, debug = FALSE
)
hirr <- 25
args$Hirr <- c(1, hirr, hirr * 0.43)
args$beta <- 1.25
args$cdr <- 0.8
## ACF: doing this makes B
## this sets when ACF (measurements) are on
ITL <- get_ITL(args)
for (i in 1:7) args$ACFhaz0[i, ITL[[i]]] <- args$ACFhaz1[i, ITL[[i]]] <- 0.2
## break points in ITL as data for plots
brks <- sort(c(unlist(lapply(ITL, min)), unlist(lapply(ITL, max))))
brk_yrs <- data.table(t = brks, yr = start_year + brks / 12)
args0 <- args

## run fwd simulation & test (un-calibrated)
out <- run.model(args, args$tt, n.particles = 200)
## plot_compare_noterate_agrgt(out, realdata = TRUE)
## TBI 21% in adults in 23/24
(tbi23 <- get_TBI_prev(out, 162, grp = "adult"))
fwrite(data.table(tbi23), file = here("output/tbi23.csv"))
## TBI over time
tbiz <- 1:160
for(i in seq_along(tbiz)) tbiz[i] <- get_TBI_prev(out, i, grp = "adult")[1]
plot(tbiz, type = "l")
## plot_HIV_in_TB(out, start_year = start_year) + xlim(c(2015, 2021))

## check un-calibrated notifications & ACF timing
gp <- plot_compare_noterate_agrgt(out,
  realdata = TRUE,
  start_year = start_year
)
gp <- gp + geom_vline(
  data = brk_yrs,
  aes(xintercept = yr),
  linetype = "dashed", col = "grey"
)
gp <- gp + labs(subtitle = "Un-calibrated model: notifications & ACF timing")

ggsave(gp, file = here("output/fit_prefit0.png"), w = 12, h = 10)


## ----------- other checks
## --- inspect demographic outputs
plot_compare_demog(out, start_year = 2015, by_comp = "age")
## NOTE we don't expect perfect agreement here because
## we data are scaled national demographic change
## and there is much higher HIV prevalence in Blantyre

## --- HIV comparisons
gp <- plot_HIV_dynamic(out,
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
gp <- plot_HIV_dynamic(out,
  start_year = start_year,
  show_ART = FALSE
)
gp <- gp +
  geom_point(data = hivpd, pch = 1, size = 2, stroke = 2) +
  xlim(c(2015, 2025)) + ylim(c(0, NA))
gp <- gp + theme(axis.text.x = element_text(angle = 55, hjust = 1))

ggsave(gp, filename = here("output/x_hivpatch.png"), w = 12, h = 7)

## HIV in TB
gp <- plot_HIV_in_TB(out, start_year = start_year) + xlim(c(2015, 2021))
gp

ggsave(gp, filename = here("output/x_hivintb.png"), w = 7, h = 5)

## === comparing trends with HIV vs notifications
real_dat <- as.data.table(BLASTtbmod::TB_notes_HIV_patch)
real_dat[, `:=`(c("year", "month"), tstrsplit(monthyr, split = "-"))]
htb <- real_dat[, .(total = sum(total)),
  by = .(year,
    hiv = ifelse(hiv == "HIV-", "HIV-", "HIV+")
  )
]
htb[, `:=`(n, sum(total)), by = year]
htb <- htb[hiv == "HIV+"]
htb[, `:=`(p, total / n)]
htb[, `:=`(s, sqrt(p * (1 - p) / n))]
htb[, `:=`(c("lo", "hi"), .(p - 1.96 * s, p + 1.96 * s))]
htb[, `:=`(year, as.numeric(year))]
lmo <- lm(p ~ year, data = htb)
lmo <- coef(lmo)[2]

real_dat <- BLASTtbmod::md7
tmp <- real_dat[round(yr) == 2017, .(cf = mean(mid)), by = comid]
real_dat <- merge(real_dat, tmp, by = "comid", all.x = TRUE)
tmp <- real_dat[yr >= 2017, .(b = coef(lm(mid ~ yr))[2]), by = comid]
real_dat <- merge(real_dat, tmp, by = "comid", all.x = TRUE)
real_dat[, y := cf + (yr - 2017) * b]
real_dat[, cf := cf + (yr - 2015) * lmo]


ggplot(real_dat, aes(yr, mid, group = comid)) +
  geom_point(shape = 21) +
  facet_wrap(~comid) +
  geom_line(aes(y = cf), col = 2) +
  geom_line(aes(y = y), col = 4) +
  geom_vline(xintercept = 2017, lty = 2, col = 1) +
  theme_linedraw() +
  expand_limits(y = 0) +
  labs(x = "Year", y = "TB notifications per month")

ggsave(file = here("output/TrendCF.png"), w = 8, h = 7)

## ========= INFERENCE

## filter extend to 6 months before 2015

## (2021 - 2015) * 12

before <- 6
args1 <- args <- restart_parms(args0, 60 - before, out, 72 + before)
str(args)

## ## ## fwd simulation & re-test
## out2 <- run.model(args, args$tt, n.particles = 200)
## ## plot_compare_noterate_agrgt(out2, realdata = TRUE)
## get_TBI_prev(out2, 1, grp = "adult") # 2015
## get_TBI_prev(out2, 71, grp = "adult") # 2021


ndat7 <- as.data.table(pf_data7)
note_names <- paste0("notifrate_", 1:7)
ndat7 <- ndat7[,..note_names]
ndat7 <- rbind(ndat7[rep(1, before)], ndat7) # FOCB
ndat7[, month := seq(1, by = 1, length.out = .N)]


pf_data7e <- mcstate::particle_filter_data(
  data = ndat7, # now using real data
  time = "month",
  rate = 1,
  initial_time = 0
)


## === reincorporating old version from 2015
start_yeare <- 2015 - before / 12
yearse <- 13 + 1 / 6 + before / 12

## ==== bunched together for easier experimentation
S <- 5 # 10*2*2
case_compare7v <- function(state, observed, pars = NULL) {
  ans <- rep(0, dim(state)[2])
  for (i in 1:7) {
    totnotes <- colSums(                           # total per month x 100K
      state[BLASTtbmod::ln7[[i]], , drop = TRUE] * # monthly per 100K
        state[BLASTtbmod::bn7[[i]], , drop = TRUE] # total population
    )
    ## print(str(totnotes))
    totpops <- colSums(state[BLASTtbmod::bn7[[i]], , drop = TRUE])
    ## print(str(totpops))
    notes_modelled <- totnotes / 1e5 # totpops
    notes_observed <- observed[[paste0("notifrate_", i)]]
    notes_observed <- rep(notes_observed, dim(state)[2]) # to match particles
    ans <- ans + dnorm(
      x = notes_modelled,
      mean = notes_observed,
      sd = (1 / S) * (notes_observed+1),
      ## sd = (1 / S) * (notes_observed + 1),
      ## sd = S,
      log = TRUE
    )
  }
  ans
}

filter <- create.particlefilter(
  pf_data7e,
  case_compare7v,
  n_particles = 50,
  n_threads = 10
)

filter$run(save_history = TRUE, pars = args)
fout <- filter$history()
plot_compare_noterate_agrgt(
  fout,
  realdata = TRUE,
  start_year = start_yeare
)


in_argsrealA <- args
in_argsrealA$beta <- NULL
in_argsrealA$cdr <- NULL
X00 <- args$popinit
in_argsrealA$popinit <- NULL
in_argsrealA$pDf <- NULL
in_argsrealA$pDs <- NULL

## common inference priors
make_transform <- function(ARGS) {
  function(theta) {
    ## remake the initial state
    X0 <- X00
    ## === TBD done as simple prevalence
    tbd <- unname(theta[3])
    ## ## === tbdi
    ## tbd <- 50e-5
    stat_nos <- apply(X0, c(2, 3, 4), sum) # total in each state
    X0["D", , 2:3, ] <- round(tbd * stat_nos[, 2:3, ])
    X0["SC", , 2:3, ] <- round(tbd * stat_nos[, 2:3, ])
    c(
      ARGS,
      list(
        beta = unname(theta[1]),
        cdr = unname(theta[2]),
        popinit = X0,
        pDf = unname(theta[4]),
        pDs = unname(theta[5])
      )
    )
  }
}

## beta
betamn <- 1.678
betasg <- 0.371
initd <- rep(50e-5,7) #1/(1+exp(-lans$par[2:8])) #NOTE
initp <- qlnorm(0.5, -2.837, 0.32) # 1/(1+exp(-lans$par[1])) #NOTE qlnorm(0.5, -2.837, 0.32)
ldm <- log(40e-5)
lds <- 0.05
sc <- 1 # 1.75
proposal_scale <- 1e5
prior_list <- list(
  beta = mcstate::pmcmc_parameter("beta",
    ## initial = results$par[1],
    initial = qlnorm(0.5, log(0.5), .25),
    min = 1e-6, max = 15,
    prior = function(x) dlnorm(x, log(0.5), .25, log = TRUE) # log(x)##
    ## prior = function(x) dlnorm(x, betamn, betasg, log = TRUE) # log(x)##
  ),
  cdr = mcstate::pmcmc_parameter("cdr",
    ## initial = results$par[2],
    initial = qbeta(0.5, 80, 20),#.75,
    min = 0.2, max = .9,
    prior = function(x) dbeta(x, 80, 20, log = TRUE)
  ),
  tbd_prev = mcstate::pmcmc_parameter("tbd_prev",
    ## initial = results$par[3],
    initial = qlnorm(0.5, ldm, lds),#20e-5,
    min = 1e-5, max = 2e-2,
    prior = function(x) dlnorm(x, ldm, lds, log = TRUE)
  ),
  pDf = mcstate::pmcmc_parameter("pDf",
    ## initial = results$par[4],
    initial = initp,
    min = 1e-6, max = 0.1,
    prior = function(x) dlnorm(x, -2.837, 0.32 / sc, log = TRUE)
  ),
  pDs = mcstate::pmcmc_parameter("pDs",
    ## initial = results$par[5],
    initial = qlnorm(0.5, -6.89 - .5 * 0, 0.58),
    min = 1e-6, max = 1e-2,#,2e-3,
    prior = function(x) dlnorm(x, -6.89 - .5 * 0, 0.58 / sc, log = TRUE)
  )
)
## proposal_matrix <- result$cov
inits <- rep(1, length(prior_list))
proposal_matrix <- diag(inits)
names(inits) <- names(prior_list)
colnames(proposal_matrix) <- rownames(proposal_matrix) <- names(prior_list)
for(nm in names(prior_list)) {
  inits[nm] <- prior_list[[nm]]$initial
  proposal_matrix[nm, nm] <- inits[nm] / proposal_scale
}


## as list
mcmc_pars <- mcstate::pmcmc_parameters$new(
  prior_list,
  proposal_matrix,
  transform = make_transform(in_argsrealA)
)
## check
tsta <- mcmc_pars$model(mcmc_pars$initial()) #looks OK
mcmc_pars$initial()
mcmc_pars$prior(mcmc_pars$initial())
sts <- apply(tsta$popinit, c(1), sum)
1e2 * sts / sum(sts)

## a: run inference
pmcmc_out <- run.pmcmc(
  particle.filter = filter,
  parms = in_argsrealA,
  n.steps = 250, n.burnin = 100, n.chains = 1,
  ## n.steps = 500, n.burnin = 250, n.chains = 1,
  ## n.steps = 1500, n.burnin = 750, n.chains = 1,
  n.threads = 10, n.epochs = 1,
  n.workers = 1,
  mcmc_pars = mcmc_pars,
  save_restart = in_argsrealA$sim_length, returnall = TRUE
)
mcmc1 <- coda::as.mcmc(cbind(
  pmcmc_out$processed_chains$probabilities,
  pmcmc_out$processed_chains$pars
))
cat("ESS:\n")
coda::effectiveSize(mcmc1)
cat("Rejection rate:\n")
coda::rejectionRate(mcmc1)
cat("================ Average acceptance rate ====================:\n")
mean(1 - coda::rejectionRate(mcmc1))
beepr::beep("coin")
(smy <- summary(mcmc1))
plot_compare_noterate_agrgt(
  pmcmc_out$pmcmc_run$trajectories$state,
  realdata = TRUE,
  start_year = start_yeare
)

pmcmc_out$proposal_matrix

pairs(pmcmc_out$processed_chains$pars)

## bayesplot::mcmc_trace(mcmc1)
p <- bayesplot::mcmc_trace(mcmc1[, c(2, 4:ncol(mcmc1))])
p

## ggsave(p, file = here("tmpdata/mcmc_diagnostics.pdf"), w = 12, h = 10)

## before/after data comparison
par_before <- inits
par_after <- smy$statistics[-c(1:3), "Mean"]
par_after / par_before


## approach to running continuation & counterfactual continuation
## NOTE this extends the *background* time-varying parameter arrays
## layers in (re)initial state from MCMC
extend_times <- function(p) {
  newparms <- p
  restart_step <- 60 - before # start of the fitting window, in args00 indexing
  sim_end <- args00$sim_length # end of args00
  keep <- restart_step:sim_end
  newparms$tt <- args00$tt[keep]
  newparms$tt <- newparms$tt - newparms$tt[1]
  newparms$sim_length <- length(newparms$tt)
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
transform_counterfactual <- function(p) {
  pars <- pmcmc_out$processed_chains$predict$transform(p)
  pars <- extend_times(pars)
  ITL <- get_ITL(pars)
  for (i in 1:7) pars$ACFhaz1[i, ITL[[i]]] <- 0.2
  for (i in 1:7) pars$ACFhaz0[i, ITL[[i]]] <- 0.2
  pars
}
transform_basecase <- function(p) {
  pars <- pmcmc_out$processed_chains$predict$transform(p)
  pars <- extend_times(pars)
  ITL <- get_ITL(pars)
  for (i in 1:7) pars$ACFhaz0[i, ITL[[i]]] <- pars$ACFhaz1[i, ITL[[i]]] <- 0.0
  pars
}
## loop over post-burnin posterior samples for basecase
basecase_pars <- lapply(
  seq_len(nrow(pmcmc_out$processed_chains$pars)),
  function(i) transform_basecase(pmcmc_out$processed_chains$pars[i, ])
)
## loop over post-burnin posterior samples for counterfactual
counterfactual_pars <- lapply(
  seq_len(nrow(pmcmc_out$processed_chains$pars)),
  function(i) transform_counterfactual(pmcmc_out$processed_chains$pars[i, ])
)

## Continuing from the actual particle-filtered state that PMCMC ended on
predict_time <- pmcmc_out$processed_chains$predict$time
filtered_state <- pmcmc_out$processed_chains$state

## basecase model: continue from the filtered state
bcmod <- BLASTtbmod:::stocm$new(
  basecase_pars,
  predict_time,
  n_particles = 1, seed = 1,
  deterministic = TRUE, pars_multi = TRUE
)
bcmod$update_state(state = filtered_state)
## counterfactual model: continue from the same filtered state
cfmod <- BLASTtbmod:::stocm$new(
  counterfactual_pars,
  predict_time,
  n_particles = 1, seed = 1,
  deterministic = TRUE, pars_multi = TRUE
)
cfmod$update_state(state = filtered_state)

## simulate forward from predict_time to the end of the extended horizon
## dim(res) = (states x particles x samples x time)
output_steps <- predict_time:length(basecase_pars[[1]]$tt)
res0 <- bcmod$simulate(output_steps)
res1 <- cfmod$simulate(output_steps)

## prepend the PMCMC's own fitted (filtered) trajectory
combine_traj <- function(traj, cont) {
  T1 <- dim(traj)[3]
  T2 <- dim(cont)[3]
  out <- array(NA_real_, dim = c(dim(traj)[1], dim(traj)[2], T1 + T2 - 1))
  out[, , 1:T1] <- traj
  out[, , (T1 + 1):(T1 + T2 - 1)] <- cont[, , 2:T2]
  out
}
traj <- pmcmc_out$processed_chains$trajectories$state
ress0 <- combine_traj(traj, res0[, 1, , ])
ress1 <- combine_traj(traj, res1[, 1, , ])
get_TBI_prev(ress0, 1, grp = "adult")


## fn <- here("tmpdata")
## if (!file.exists(fn)) dir.create(fn)

## save(ress0, file = here("tmpdata/ress0.Rdata"))
## save(ress1, file = here("tmpdata/ress1.Rdata"))

## ===========================================
## CHECKS
## NOTE uses the function from 1b_fit_analysis.R

## extract data (NOTE memory hungry & time consuming)
E0 <- formplotdata(ress0, eps = 0.05, use_median = FALSE)
E1 <- formplotdata(ress1, eps = 0.05, use_median = FALSE)
E0[, acf := "No ACF"]
E1[, acf := "ACF"]
EB <- rbind(E0, E1)

## notifications
real_dat <- BLASTtbmod::md7
real_dat[["patch"]] <- paste("Patch", md7$comid)
real_dat$qty <- "noterate"
real_dat$acf <- "No ACF"
real_dat[, yr := 2015 + t / 12]
## veritcal lines data
VL <- data.table(
  patch = rep(EB[, unique(patch)], each = 2),
  t = EB[, max(t)] - c(0, rep(1:6, each = 2), 7) * 12
)
VL[, item := rep(2:1, 7)]
VLW <- dcast(VL, patch ~ item, value.var = "t")
names(VLW)[2:3] <- c("bot", "top")
VL[, yr := start_yeare + t / 12]
## difference
ED <- dcast(EB[qty == "cummort", .(t, patch, mid, acf)],
  t + patch ~ acf,
  value.var = "mid"
)
ED[, ddf := `No ACF` - ACF]
ED <- merge(ED, VLW, by = "patch")
ED[!(t <= top & t >= bot), ddf := NA_real_]
ED[, mnddf := min(ddf, na.rm = TRUE), by = patch]
ED[!is.na(ddf)]
ED[, ddf := ddf - mnddf]
ED[
  ,
  c("qty", "mid", "lo", "hi", "acf") := .(
    "ddf", ddf, NA_real_, NA_real_, "No ACF - ACF"
  )
]
EB <- rbind(EB, ED[, .(t, patch, qty, hi, lo, mid, acf)])
EB[, yr := start_yeare + t / 12]
## renaming
EB[, zone := gsub("Patch", "Zone", patch)]
VL[, zone := gsub("Patch", "Zone", patch)]
real_dat[, zone := gsub("Patch", "Zone", patch)]

EB[, nqty := fcase(
  qty == "noterate", "Notifications per 100,000 per month",
  qty == "mortrate", "Deaths per 100,000 per month",
  qty == "cummort", "Cumulative deaths",
  qty == "ddf", "Difference in cumulative deaths"
)]
real_dat[, nqty := fcase(
  qty == "noterate", "Notifications per 100,000 per month",
  qty == "mortrate", "Deaths per 100,000 per month",
  qty == "cummort", "Cumulative deaths",
  qty == "ddf", "Difference in cumulative deaths"
)]
## simplified version
EBR <- EB[qty == "noterate"]

## plot
plt <- "Accent"
tst <- ggplot(EBR, aes(yr,
  y = mid, ymin = lo, ymax = hi, col = acf, fill = acf,
  group = paste(patch, qty, acf)
)) +
  geom_ribbon(alpha = 0.3, col = NA) +
  geom_line(lwd = 1) +
  scale_fill_brewer(palette = plt) +
  scale_color_brewer(palette = plt) +
  geom_vline(
    data = VL,
    aes(xintercept = yr),
    lty = 2, col = "darkgrey"
  ) +
  facet_wrap(~zone) +
  xlab("Time") +
  ylab("TB notifications per month") +
  geom_point(data = real_dat, col = 2, shape = 1) +
  theme_linedraw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.6, 0.15),
    legend.title = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 1), col = guide_legend(nrow = 1)) +
  xlim(2015, NA)
  ## xlim(start_yeare, NA)
tst
beepr::beep("coin")

ggsave(tst, file = here("tmpdata/fitng3.png"), w = 12, h = 10)
