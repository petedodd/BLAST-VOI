log_posterior <- function(newpars, basepars, pf, compare, prior, verbose = FALSE) {
  ## update parameters
  pars <- basepars
  for(p in names(newpars)) pars[p] <- newpars[p]
  ## run particle filter
  ll <- pf$run(pars)
  if(verbose) cat("ll =", ll, "\n")
  ## extract log-prior
  lp <- prior(newpars)
  if(verbose) cat("lp =", lp, "\n")
  ## return log-posterior
  ll + lp
}

## testing
log_posterior(c(pDf = 0.1, ari0 = 0.01),
              args,
              filter,
              case_compare7,
              mcmc_pars$prior,
              verbose = TRUE)

## test in
testin <- c(pDf = 0.1,
            d1 = 50e-5, d2 = 50e-5, d3 = 50e-5, d4 = 50e-5,
            d5 = 50e-5, d6 = 50e-5, d7 = 50e-5)

log_posterior(testin,
              args,
              filter,
              case_compare7,
              mcmc_pars$prior,
              verbose = TRUE)


## test in
testin <- c(beta = 1, cdr = 0.5,
            d1 = 50e-5, d2 = 50e-5, d3 = 50e-5, d4 = 50e-5,
            d5 = 50e-5, d6 = 50e-5, d7 = 50e-5)

log_posterior(testin,
              args,
              filter,
              case_compare7,
              mcmc_pars$prior,
              verbose = TRUE)


## ----
lpf <- function(x){
  x <- 1/(1+exp(-x))
  names(x) <- names(testin)
  print(x)
  log_posterior(x, args, filter, case_compare7, mcmc_pars$prior, verbose = FALSE)
}


lpf(log(testin))


lans <- optim(
  par = log(testin),
  fn = lpf,
  method = "Nelder-Mead",
  control = list(maxit = 5000, trace = TRUE, REPORT = 1, fnscale = -1)
)
beepr::beep('coin')


lans
1/(1+exp(-lans$par))
1e5/(1+exp(-lans$par[-1]))

## ============ working with beta

## testing
log_posterior(c(pDf = 0.1, ari0 = 0.01),
              args,
              filter,
              case_compare7,
              mcmc_pars$prior,
              verbose = TRUE)

## test in
testin <- c(beta = 2,
            d1 = 50e-5, d2 = 50e-5, d3 = 50e-5, d4 = 50e-5,
            d5 = 50e-5, d6 = 50e-5, d7 = 50e-5)

log_posterior(testin,
              args,
              filter,
              case_compare7,
              mcmc_pars$prior,
              verbose = TRUE)



## ----
lpf <- function(x){
  x <- 1/(1+exp(-x))
  x[1] <- x[1] / (1-x[1]) #convert to exp
  names(x) <- names(testin)
  print(x)
  log_posterior(x, args, filter, case_compare7, mcmc_pars$prior, verbose = FALSE)
}


lpf(log(testin))


lans <- optim(
  par = log(testin),
  fn = lpf,
  method = "Nelder-Mead",
  control = list(maxit = 5000, trace = TRUE, REPORT = 1, fnscale = -1)
)
beepr::beep('coin')


lans <- optim(
  par = log(testin),
  fn = lpf,
  method = "SANN",
  control = list(maxit = 5000, trace = TRUE, REPORT = 1, fnscale = -1)
)
beepr::beep('coin')


lpf(lans$par)

## TODO visualize

argso <- copy(args)
argso$beta <- exp(lans$par[1])
argso$cdr <- 1/(1 + exp(-lans$par[2]))
argso$initD[,2] <- argso$initD[,3] <- 1/(1+exp(-lans$par[3:9]))

argso$initD[,2] * 1e5

outo <- run.model(argso, args$tt, n.particles = 200)
gp <- plot_compare_noterate_agrgt(outo, realdata = TRUE)
gp <- gp + scale_y_continuous(limits = c(0, 50))
ggsave(gp,file = here("tmpdata/optout.png"), w = 12, h = 10)

## very odd in  couple of communities at start



lans
1/(1+exp(-lans$par))
1e5/(1+exp(-lans$par[-1]))




## ----
## test in
testin <- c(beta = 2, pDf = 0.06,
            d1 = 50e-5, d2 = 50e-5, d3 = 50e-5, d4 = 50e-5,
            d5 = 50e-5, d6 = 50e-5, d7 = 50e-5)


lpf <- function(x){
  x <- 1/(1+exp(-x))
  x[1] <- x[1] / (1-x[1]) #convert to exp
  names(x) <- names(testin)
  print(x)
  log_posterior(x, args, filter, case_compare7, mcmc_pars$prior, verbose = FALSE)
}


lpf(log(testin))


## lans <- optim(
##   par = log(testin),
##   fn = lpf,
##   method = "Nelder-Mead",
##   control = list(maxit = 5000, trace = TRUE, REPORT = 1, fnscale = -1)
## )
## beepr::beep('coin')


lans <- optim(
  par = log(testin),
  fn = lpf,
  method = "SANN",
  control = list(maxit = 5000, trace = TRUE, REPORT = 1, fnscale = -1)
)
beepr::beep('coin')



lans
1/(1+exp(-lans$par))
1e5/(1+exp(-lans$par[-1]))


args3 <- args
args3$initD[,2] <- args3$initD[,3] <- 1/(1+exp(-lans$par[3:9]))
args3$pDf <- 1/(1+exp(-lans$par[2]))
args3$beta <- exp(lans$par[1])

out3 <- run.model(args3, args3$tt, n.particles = 200)
gp <- plot_compare_noterate_agrgt(out3, realdata = TRUE)
ggsave(gp,filename = here("tmpdata/optimization.png"), width = 8, height = 6)


BLASTtbmod::TBN[,summary(notifrate), by = zone]


## =========================================
## extracting TBI prevalence
BLASTtbmod::get_cols

recent <- c(241:303)
nonrecent <- c(304:366)
pop <- c(619:681)

## np x na x nh
BLASTtbmod::get_cols[recent]
BLASTtbmod::get_cols[nonrecent]
BLASTtbmod::get_cols[pop]

## 15-24
BLASTtbmod::get_cols[recent]
BLASTtbmod::get_cols[nonrecent]
BLASTtbmod::get_cols[pop]


recentA <- grep(",2,",BLASTtbmod::get_cols[recent], value  = TRUE)
recentA <- which(BLASTtbmod::get_cols %in% recentA)
nonrecentA <- grep(",2,",BLASTtbmod::get_cols[nonrecent], value  = TRUE)
nonrecentA <- which(BLASTtbmod::get_cols %in% nonrecentA)
popA <- grep(",2,",BLASTtbmod::get_cols[pop], value = TRUE)
popA <- which(BLASTtbmod::get_cols %in% popA)

## CHECK
BLASTtbmod::get_cols[recentA]
BLASTtbmod::get_cols[nonrecentA]
BLASTtbmod::get_cols[popA]

get_TBI <- function(X, tindex = 1){
  totpop <- sum(X[popA,,tindex])
  recpop <- sum(X[recentA,,tindex])
  nonpop <- sum(X[nonrecentA,,tindex])
  c(recent_tbi = recpop / totpop,
    nonrecent_tbi = nonpop / totpop,
    tot_tbi = (recpop + nonpop) / totpop)
}

get_TBI(outo, tindex = 1)
get_TBI(out, tindex = 1) #90+%
get_TBI(out, tindex = dim(out)[3])      #70%
get_TBI(out2, tindex = 1) #90+% NOTE should not be
get_TBI(out2, tindex = dim(out2)[3])      #70%

## TODO check Peter MacP data 21% in adults in 2023/2024
get_TBI(out, tindex = 1) #90+%
get_TBI(out, tindex = dim(out)[3])      #70%


dim(out2)

## 2023.5 is mid 2023, 2015 is start of simulation
(2023.5 - 2015)*12 #102
get_TBI(out2, tindex = 102)  #want 21%


get_TBI(out, tindex = 1)
get_TBI(out, tindex = dim(out)[3])
get_TBI(out2, tindex = 1)
get_TBI(out2, tindex = 102)


restart_parms2 <- function (parms, restart_step, end_state) {
    cat("Creating new parameters to restart at given step...\n")
    cat("Summarising end state...\n")
    denom <- extract.pops.multi(end_state, dim(end_state)[2], 
                                out_type = "N")
    numer <- extract.pops.multi(end_state, dim(end_state)[2], 
        out_type = "D")
    hnumr <- denom[t == restart_step, .(N = sum(N)), by = .(patch, 
        age, hiv, particle = chain_step)]
    denom <- denom[t == restart_step, .(N = sum(N)), by = .(patch, 
        age, particle = chain_step)]
    numer <- numer[t == restart_step, .(D = sum(D)), by = .(patch, 
        age, particle = chain_step)]
    both <- data.table::merge.data.table(denom, numer, by = c("patch", 
        "age", "particle"))
    boths <- both[, .(prev = mean(D/N)), by = .(patch, age)]
    D00 <- data.table::dcast(boths, patch ~ age, value.var = "prev")
    D00 <- as.matrix(D00[, -1])
    hnumr[, `:=`(tot, sum(N)), by = .(patch, age, particle)]
    hnumr[, `:=`(p, N/tot)]
    hnumr <- hnumr[, .(p = mean(p)), by = .(patch, age, hiv)]
    H00 <- array(hnumr[order(hiv, age, patch)]$p, c(7, 3, 3), 
        dimnames = list(patch = unique(hnumr$patch), age = unique(hnumr$age), 
                        hiv = unique(hnumr$hiv)))
    ## get ari0 to match TBI rates
    tbiz <- get_TBI(end_state, tindex = restart_step)
    ari <- - log(1 - tbiz["tot_tbi"]) / 45
    cat("Estimated ARI at restart step:", ari, "\n")
    cat("Updating parameter object...\n")
    newparms <- parms
    newparms$initD <- D00
    newparms$propinit_hiv <- H00
    keep <- restart_step:length(parms$tt)
    newparms$tt <- newparms$tt[keep]
    newparms$tt <- newparms$tt - newparms$tt[1]
    newparms$sim_length <- length(newparms$tt)
    newparms$births_int <- newparms$births_int[keep]
    newparms$HIV_int <- newparms$HIV_int[keep]
    newparms$ART_int <- newparms$ART_int[keep]
    newparms$mu_noHIV_int <- newparms$mu_noHIV_int[, keep]
    newparms$mu_HIV_int <- newparms$mu_HIV_int[, keep]
    newparms$mu_ART_int <- newparms$mu_ART_int[, keep]
    newparms$m_in_int <- newparms$m_in_int[, keep]
    newparms$ACFhaz0 <- newparms$ACFhaz0[, keep]
    newparms$ACFhaz1 <- newparms$ACFhaz1[, keep]
    newparms$ari0 <- ari
    return(newparms)
}

