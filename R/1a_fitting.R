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
out <- run.model(args, args$tt, n.particles = 200)

## check un-calibrated notifications & ACF timing
gp <- plot_compare_noterate_agrgt(out,
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
gp

ggsave(gp, filename = here("output/x_hivpatch.png"), w = 7, h = 7)

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

## update parameters for restart in 2015
args0 <- args #for safe-keeping
args1 <- args <- restart_parms(args0, 60, out)

## === reincorporating old version from 2015
start_year <- 2015
years <- 13 + 1 / 6
## CHECK
length(as.double(seq(0, 12 * years)))
args$sim_length


## D00 <- args$initD
## args$initD <- D00 / 4
## out2 <- run.model(args, args$tt, n.particles = 100)
## plot_compare_noterate_agrgt(out2, realdata = TRUE)

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
out2 <- run.model(args, args$tt, n.particles = 200)
## plot_compare_noterate_agrgt(out2, realdata = TRUE)

## proposal_matrix <- A$proposal_matrix
## save(proposal_matrix, file = here("tmpdata/proposal_matrix.Rdata"))


## ## A: 'TAU' inference
## for (i in 1:7) args$ACFhaz0[i, ITL[[i]]] <- args$ACFhaz1[i, ITL[[i]]] <- 0 # zero again
## in_argsrealA <- args
## in_argsrealA$beta <- in_argsrealA$initD <- NULL # in_argsrealA$ari0 <-
## ## in_argsrealA$ari0 <- in_argsrealA$initD <- NULL
## ## in_argsrealA$pDf <- NULL

## ======== NOTES
## this is working OK ish
## not replicating trend down
## args$HIV_dur_ratio <- 6 # BUG checking
## args$tfr <- 0.02661832
## args$tfr / 0.5
## args$cfr / args$dur

## ====================================== working
##   ## pDs = mcstate::pmcmc_parameter("pDs",
##   ##   initial = qlnorm(0.5, -6.89, 0.58),
##   ##   min = 1e-6, max = 1, prior = function(x) dlnorm(x, -6.89, 0.58)
##   ##   ),
##   pDf = mcstate::pmcmc_parameter("pDf",
##     initial = qlnorm(0.5, -2.837, 0.32),
##     min = 1e-6, max = 1, prior = function(x) dlnorm(x, -2.837, 0.32)
##     ),
## ==== bunched together for easier experimentation
S <- 5 # 10*2*2
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
      sd = (1 / 5) * notes_observed, log = TRUE
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
for (i in 1:7) args$ACFhaz0[i, ITL[[i]]] <- args$ACFhaz1[i, ITL[[i]]] <- 0 # zero again
in_argsrealA <- args
in_argsrealA$beta <- NULL
in_argsrealA$pDf <- NULL
## in_argsrealA$pDs <- NULL
in_argsrealA$initD <- NULL # in_argsrealA$ari0 <-
## common inference priors
make_transform <- function(ARGS) {
  function(theta) {
    c(
      ARGS,
      list(
        beta = unname(theta[1]),
        pDf = unname(theta[2]),
        ## pDs = unname(theta[3]),
        ## ari0 = unname(theta[2]),
        initD = cbind(
          rep(0, 7),
          unname(theta[3:9]),
          unname(theta[3:9])
        )
      )
    )
  }
}
## curve(dlnorm(x, log(150 / 1e5), 0.99), n = 1e3, from = 0, to = 0.01)
D0pr <- function(x) dlnorm(x, log(1e1 / 1e5), 1)
## beta
betamn <- 1.678
betasg <- 0.371
initd <- 50e-5
prior_list <- list(
  beta = mcstate::pmcmc_parameter("beta",
    initial = qlnorm(0.25, betamn, betasg),
    min = 1e-6, max = 15,
    prior = function(x) dlnorm(x, betamn, betasg)
  ),
  pDf = mcstate::pmcmc_parameter("pDf",
    initial = qlnorm(0.5, -2.837, 0.32),
    min = 1e-6, max = 1, prior = function(x) dlnorm(x, -2.837, 0.32)
  ),
  ## pDs = mcstate::pmcmc_parameter("pDs",
  ##   initial = qlnorm(0.5, -6.89, 0.58),
  ##   min = 1e-6, max = 1, prior = function(x) dlnorm(x, -6.89, 0.58)
  ##   ),
  d1 = mcstate::pmcmc_parameter("d1",
    initial = initd,
    min = 1e-6, max = 2e-2, prior = D0pr
  ),
  d2 = mcstate::pmcmc_parameter("d2",
    initial = initd,
    min = 1e-6, max = 2e-2, prior = D0pr
  ),
  d3 = mcstate::pmcmc_parameter("d3",
    initial = initd,
    min = 1e-6, max = 2e-2, prior = D0pr
  ),
  d4 = mcstate::pmcmc_parameter("d4",
    initial = initd,
    min = 1e-6, max = 2e-2, prior = D0pr
  ),
  d5 = mcstate::pmcmc_parameter("d5",
    initial = initd,
    min = 1e-6, max = 2e-2, prior = D0pr
  ),
  d6 = mcstate::pmcmc_parameter("d6",
    initial = initd,
    min = 1e-6, max = 2e-2, prior = D0pr
  ),
  d7 = mcstate::pmcmc_parameter("d7",
    initial = initd,
    min = 1e-6, max = 2e-2, prior = D0pr
  )
)
proposal_matrix <- diag(c(
  0.05, # beta
  1e-3, # pDf
   ## 1e-5, # pDs
  ## 1e-5, # ari0
  rep(1e-9, 7)
))


## as list
mcmc_pars <- mcstate::pmcmc_parameters$new(
  prior_list,
  proposal_matrix,
  transform = make_transform(in_argsrealA)
)
## check
mcmc_pars$model(mcmc_pars$initial()) #looks OK
mcmc_pars$initial()

## A: run inference
pmcmc_out <- run.pmcmc(
  particle.filter = filter,
  parms = in_argsrealA,
  n.steps = 250, n.burnin = 100, n.chains = 1,
  ## n.steps = 500, n.burnin = 250, n.chains = 1,
  ## n.steps = 1500, n.burnin = 750, n.chains = 1,
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
  for (i in 1:7) pars$ACFhaz1[i, ITL[[i]]] <- 0.2
  for (i in 1:7) pars$ACFhaz0[i, ITL[[i]]] <- 0.2
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
## CHECKS
## NOTE uses the function from 1b_fit_analysis.R

## extract data (NOTE memory hungry & time consuming)
E0 <- formplotdata(ress0, eps = 0.05)
E1 <- formplotdata(ress1, eps = 0.05)
E0[, acf := "No ACF"]
E1[, acf := "ACF"]
EB <- rbind(E0, E1)
start_year <- 2015 # start year used in simulation

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
VL[, yr := start_year + t / 12]
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
EB[, yr := start_year + t / 12]
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
ggplot(EBR, aes(yr,
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
  ## xlim(2015, NA)
  xlim(start_year, NA)



