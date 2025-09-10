## generating data for plot
library(BLASTtbmod)
library(data.table)
library(ggplot2)

data(in_args) #example
data(pf_data7)

## define other parameters
n_particles <- 10
time <- 1 #NOTE not used TODO check
start_year <- 2015
years <- 6

## create common compare function and filter
S <- 10
case_compare7v <- function (state, observed, pars = NULL) {
  ans <- rep(0, dim(state)[2])
  for (i in 1:7) {
    totnotes <- colSums(state[BLASTtbmod::ln7[[i]], , drop = TRUE] * 
                        state[BLASTtbmod::bn7[[i]], , drop = TRUE])
    totpops <- colSums(state[BLASTtbmod::bn7[[i]], , drop = TRUE])
    notes_modelled <- totnotes/totpops
    notes_observed <- observed[[paste0("notifrate_", i)]]
    ans <- ans + dnorm(x = notes_modelled, mean = notes_observed, 
                       sd = S, log = TRUE)
  }
  ans
}
filter <- create.particlefilter( pf_data7, case_compare7v, n_particles = 100, n_threads = 4)

## change parms
args <- get.parms(start_year=start_year,years=years+7)
args$ari0 <- 5e-2
args$initD[,2:3] <- 150e-5
args$beta <- 5
args$cdr <- 0.8
str(args)

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
for(i in 1:7) args$ACFhaz0[i,ITL[[i]]] <- args$ACFhaz1[i,ITL[[i]]] <- 0.2

## fwd simulation & test
test <- run.model(args,args$tt,n.particles=200)
plot_compare_noterate_agrgt(test,realdata=FALSE)


## --- D0 inference?
## args$initD

curve(dlnorm(x,log(150/1e5),0.99),n=1e3,from=0,to=0.01)
D0pr <- function(x) dlnorm(x, log(1e0/1e5), 1)

## transform?
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
for(i in 1:7) args$ACFhaz0[i,ITL[[i]]] <- args$ACFhaz1[i,ITL[[i]]] <- 0 #zero again
in_argsrealA <- args
in_argsrealA$beta <- in_argsrealA$ari0 <- in_argsrealA$initD <- NULL

## proposal_matrix <- A$proposal_matrix
## save(proposal_matrix, file = here("tmpdata/proposal_matrix.Rdata"))


prior_list <- list(
  beta = beta, ari0 = ari0,
  d1 = mcstate::pmcmc_parameter("d1", initial = 2*150e-5, min = 1e-6, max = 2e-2, prior = D0pr),
  d2 = mcstate::pmcmc_parameter("d2", initial = 150e-5, min = 1e-6, max = 2e-2, prior = D0pr),
  d3 = mcstate::pmcmc_parameter("d3", initial = 4*150e-5, min = 1e-6, max = 2e-2, prior = D0pr),
  d4 = mcstate::pmcmc_parameter("d4", initial = 4*150e-5, min = 1e-6, max = 2e-2, prior = D0pr),
  d5 = mcstate::pmcmc_parameter("d5", initial = 150e-5, min = 1e-6, max = 2e-2, prior = D0pr),
  d6 = mcstate::pmcmc_parameter("d6", initial = 150e-5, min = 1e-6, max = 2e-2, prior = D0pr),
  d7 = mcstate::pmcmc_parameter("d7", initial = 150e-5, min = 1e-6, max = 2e-2, prior = D0pr)
)


## as list
mcmc_pars <- mcstate::pmcmc_parameters$new( prior_list, proposal_matrix, transform = make_transform(in_argsrealA) )
## CHECK
mcmc_pars$initial()
mcmc_pars$model(mcmc_pars$initial()) #looks OK

## want to be able extract transform, prior list, and proposal matric
## str(mcmc_pars)
## mcmc_pars$.__enclos_env__$private$transform
## mcmc_pars$.__enclos_env__$private$proposal_kernel
## mcmc_pars$.__enclos_env__$private$parameters


## A: run inference
A <- run.pmcmc(
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


## str(A)
## pmcmc_run=pmcmc_run,
## processed_chains=processed_chains,
## proposal_matrix=proposal.matrix
## TODO inference results
## processed_chains <- mcstate::pmcmc_thin(A)
mcmc1 <- coda::as.mcmc(cbind(A$processed_chains$probabilities, A$processed_chains$pars))
summary(mcmc1)
coda::effectiveSize(mcmc1)
coda::rejectionRate(mcmc1)
## plot(mcmc1)

## check
## plot_compare_noterate_agrgt( A$processed_chains$trajectories$state )


## approach to running counterfactual
transform_counterfactual <- function(p) {
  pars <- A$pmcmc_run$predict$transform(p)
  for (i in 1:7) pars$ACFhaz0[i, ITL[[i]]] <- pars$ACFhaz1[i, ITL[[i]]] <- 0.2
  pars
}
transform_basecase <- function(p) {
  pars <- A$pmcmc_run$predict$transform(p)
  for (i in 1:7) pars$ACFhaz0[i, ITL[[i]]] <- pars$ACFhaz1[i, ITL[[i]]] <- 0.0
  pars
}


## loop over parameters for basecase
basecase_pars <- lapply(
  seq_len(nrow(A$pmcmc_run$pars)),
  function(i) transform_basecase(A$pmcmc_run$pars[i, ])
)

## loop over parameters for counterfactual
counterfactual_pars <- lapply(
  seq_len(nrow(A$pmcmc_run$pars)),
  function(i) transform_counterfactual(A$pmcmc_run$pars[i, ])
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
output_steps <- seq(0,maxt)

## res will have dimensions (states x particles x samples x time)
res0 <- bcmod$simulate(output_steps)
res1 <- cfmod$simulate(output_steps)
str(res0)

## NOTE new safeties seem to have helped

ress0 <- res0[,1,,]
str(ress0)
ress1 <- res1[,1,,]
str(ress1)

## save(ress0,file="~/Downloads/ress0.Rdata")
## save(ress1,file="~/Downloads/ress1.Rdata")


## ===== look!
## NOTE this is horibly memory hungry
## see matched fitplots
formplotdata <- function(Y, eps = 0.25) {
  n_chain_steps <- dim(Y)[2]
  D1 <- extract.pops.multi(Y, n_chain_steps, out_type = "N")
  D1 <- D1[, .(N = sum(N)), by = .(chain_step, t, patch)]
  D2 <- extract.pops.multi(Y, n_chain_steps, out_type = "notes")
  D2 <- D2[, .(notes = sum(notes)), by = .(chain_step, t, patch)]
  D3 <- extract.pops.multi(Y, n_chain_steps, out_type = "TB_deaths")
  D3 <- D3[, .(TB_deaths = sum(TB_deaths)), by = .(chain_step, t, patch)]
  D <- merge(D1, D2, by = c("chain_step", "t", "patch"))
  D <- merge(D, D3, by = c("chain_step", "t", "patch"))
  D[, cummort := cumsum(TB_deaths), by = .(chain_step, patch)]
  aggD <- D[, .(
    noterate.mid = 1e5 * median(notes / N),
    noterate.lo = 1e5 * quantile(notes / N, eps),
    noterate.hi = 1e5 * quantile(notes / N, 1 - eps),
    mortrate.mid = 1e5 * median(TB_deaths / N),
    mortrate.lo = 1e5 * quantile(TB_deaths / N, eps),
    mortrate.hi = 1e5 * quantile(TB_deaths / N, 1 - eps),
    cummort.mid = median(cummort),
    cummort.lo = quantile(cummort, eps),
    cummort.hi = quantile(cummort, 1 - eps)
  ), by = .(t, patch)]
  aggD <- aggD[t > 1]
  aggD <- melt(aggD, id = c("t", "patch"))
  aggD[, c("qty", "type") := tstrsplit(variable, "\\.")]
  aggD <- dcast(aggD, t + patch + qty ~ type, value.var = "value")
  aggD
}
E0 <- formplotdata(ress0,eps=0.05)
E1 <- formplotdata(ress1,eps=0.05)
E0[,acf:="No ACF"]
E1[,acf:="ACF"]
EB <- rbind(E0,E1)
real_dat <- BLASTtbmod::md7
real_dat[["patch"]] <- paste("Patch", md7$comid)
real_dat$qty <- "noterate"
real_dat$acf <- "No ACF"
VL <- data.table(
  patch = rep(EB[, unique(patch)], each = 2),
  t = EB[, max(t)] - c(0, rep(1:6, each = 2), 7) * 12
)
VL[, item := rep(2:1, 7)]
VLW <- dcast(VL, patch ~ item, value.var = "t")
names(VLW)[2:3] <- c("bot", "top")
ED <- dcast(EB[qty == "cummort", .(t, patch, mid, acf)], t + patch ~ acf, value.var = "mid")
ED[, ddf := `No ACF` - ACF]
ED <- merge(ED, VLW, by = "patch")
ED[!(t <= top & t >= bot), ddf := NA_real_]
ED[, mnddf := min(ddf, na.rm = TRUE), by = patch]
ED[!is.na(ddf)]
ED[, ddf := ddf - mnddf]
ED[, c("qty", "mid", "lo", "hi", "acf") := .("ddf", ddf, NA_real_, NA_real_, "No ACF - ACF")]
EB <- rbind(EB, ED[, .(t, patch, qty, hi, lo, mid, acf)])
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
EBR <- EB[qty == "noterate"]


## plot
plt <- "Accent"
ggplot(EBR,
       aes(2015 + t / 12,
           y = mid, ymin = lo, ymax = hi,
           col = acf, fill = acf,
           group = paste(patch, qty, acf))) +
  geom_ribbon(alpha = 0.3, col = NA) +
  geom_line(lwd = 1) +
  scale_fill_brewer(palette = plt) +
  scale_color_brewer(palette = plt) +
  geom_vline(data = VL, aes(xintercept = 2015 + t / 12),
             lty = 2, col = "darkgrey") +
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
  guides(fill = guide_legend(nrow = 1), col = guide_legend(nrow = 1))

ggsave(file = "~/Downloads/BLASTnewfit.png", w = 7, h = 5)
