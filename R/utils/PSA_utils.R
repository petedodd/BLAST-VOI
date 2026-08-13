## various helper functions for the N-patches VOI

## ====== CHANGING PARMS
## variables with 7 in:
## popinit
## patch_dims
## initD
## IRR
## MM
## ACFhaz0
## ACFhaz1

## more realistic (average, non-zone-specific) age x HIV population structure
## typical adult HIV prevalence is ~13-21% (BLASTtbmod::blantyre$hivprev)
## TB progression hazards are multiplied by Hirr = c(1, 25, 10.75)
## *average* used as a one-off default: from get.parms()'s own real 7-zone
## popinit (Malawi age pyramid +  Blantyre adult HIV prevalence), by summing
realisticAgeHIVfrac <- function() {
  refargs <- get.parms(start_year = 2015, years = 6)
  ## age x HIV, summed over state & zone
  f <- apply(refargs$popinit, c(3, 4), sum)
  f / sum(f)
}


## function to update args to new default
updateARGS <- function(args, numpatches = 7, simlength = 73) {
  ## defn things with 7 in:
  args$patch_dims <- numpatches

  ## defn things with 73 in:
  args$sim_length <- simlength
  args$tt <- seq(0, by = 1, length.out = simlength)

  ## adjust args depending on number of patches
  ## popinit must be array (compartment, patch, age, HIV)
  age_dims <- 3
  HIV_dims <- 3
  patchpop <- 500e3 / numpatches
  ageHIVfrac <- realisticAgeHIVfrac() # age x HIV, sums to 1 -- see above
  diseasefrac <- 0.5e-2 # matches the old (now-unused) initD default
  popinit <- array(0, dim = c(7, numpatches, age_dims, HIV_dims))
  for (a in 1:age_dims) {
    for (h in 1:HIV_dims) {
      cellpop <- patchpop * ageHIVfrac[a, h]
      if (a == 1) { # age class 1 (children): fully uninfected (U)
        popinit[1, , a, h] <- cellpop
      } else {
        popinit[1, , a, h] <- cellpop * (1 - diseasefrac) # adults: mostly U
        popinit[4, , a, h] <- cellpop * diseasefrac # ...with a small D seed
      }
    }
  }

  ## construct ensure integer NOTE some are later overwritten
  args$popinit <- round(popinit)
  args$initD <- matrix(0, nrow = numpatches, ncol = 3)
  args$initD[, 2:3] <- 0.5e-2
  args$IRR <- rep(1, numpatches)
  args$MM <- diag(rep(1, numpatches)) +
    matrix(0.1, nrow = numpatches, ncol = numpatches)
  args$ACFhaz0 <- args$ACFhaz1 <- matrix(0, nrow = numpatches, ncol = simlength)

  ## simpler demographics
  args$HIV_int <- args$ART_int <- rep(0, simlength) # HIV/ART incidence
  args$births_int <- rep(40, simlength) # births
  args$mu_HIV_int <-
    args$mu_ART_int <- matrix(0, nrow = 3, ncol = simlength) # HIV/ART mortality
  args$mu_noHIV_int <- matrix(1 / 50, nrow = 3, ncol = simlength) # background
  args$m_in_int <- matrix(0, nrow = 3, ncol = simlength) # migration

  ## return:
  args
}


## overwrite the disease (D) vs. uninfected (U) split
## split *across HIV strata* within each patch/age cell
setPrevalence <- function(args, prevz) {
  popinit <- args$popinit
  poptot <- apply(popinit, c(2, 3), sum) # patch x age, total pop
  for (p in seq_along(prevz)) {
    for (a in 2:3) { # adult age classes only, matching updateARGS()
      byhiv <- colSums(popinit[c(1, 4), p, a, ]) # current pop by HIV stratum
      hivfrac <- byhiv / sum(byhiv) # preserve this cell's realistic HIV split
      cellpop <- poptot[p, a] * hivfrac
      popinit[1, p, a, ] <- cellpop * (1 - prevz[p])
      popinit[4, p, a, ] <- cellpop * prevz[p]
    }
  }
  args$popinit <- round(popinit) # ensure integers
  args
}



## make a schedule for ACF
makeITZ <- function(pops0, screenrate, burnin = 0) {
  nzones <- length(pops0) # number of zones
  TZ <- 0.5 * pops0 / screenrate # time in each zone (in years)
  ## round up times in months/steps,
  TZ <- ceiling(12 * TZ - 1e-3)
  ## cap window length at 10 months
  ## but avoid population-rounding noise
  TZ[TZ > 10] <- 10
  sum(TZ) # ~60 months
  CTZ <- c(0, cumsum(TZ)) # cumulative
  etz <- stz <- CTZ + 2
  stz <- rev(rev(stz)[-1])
  etz <- etz[-1]
  etz <- etz - 1
  ## cbind(stz,etz) #check
  ## set timings base on above TZ
  itz <- list()
  for (i in 1:nzones) itz[[i]] <- (burnin + stz[i]):(burnin + etz[i])
  itz
}



## assign ACF hazard schedule
setACF <- function(args, screenrate, itz) {
  npatch <- args$patch_dims
  ## popinit is (compartment, patch, age, HIV)
  patchpop <- apply(args$popinit, 2, sum)
  for (i in 1:npatch) {
    args$ACFhaz0[i, ] <- args$ACFhaz1[i, ] <- 0
    if (length(itz[[i]]) > 0) {
      args$ACFhaz0[i, itz[[i]]] <-
        args$ACFhaz1[i, itz[[i]]] <-
        screenrate / patchpop[i]
    }
  }
  args
}




## ====== DATA EXTRACTION

## a function to get the new column names
getnewcolnames <- function(args) {
  tmpmod <- BLASTtbmod:::stocm$new(
    pars = args,
    time = min(args$tt) + 1,
    n_particles = 1
  )
  return(tmpmod$.__enclos_env__$private$info_) # column names
}


## functions to get deaths & differences quickly from model output
getMeanDeaths <- function(X, indices) {
  nz <- length(indices)
  D <- matrix(NA, nrow = dim(X)[3], ncol = nz)
  for (j in 1:nz) {
    tt <- apply(X[indices[[j]], , ], 2:3, sum) # sum over strata in each zone
    D[, j] <- colMeans(tt) # mean across simulations at each t
  }
  D
}


## reformat death differences
getMeanDD <- function(X, X0, itz, indices) {
  nz <- length(itz)
  ## apply
  D <- getMeanDeaths(X, indices) # month x zone, per-timestep flow
  D0 <- getMeanDeaths(X0, indices)
  ## cumulative deaths averted by ACF
  DD <- apply(D0, 2, cumsum) - apply(D, 2, cumsum) # benefit
  ## differences
  DD <- as.data.table(DD)
  names(DD) <- paste("Zone", 1:nz)
  DD[, t := seq_len(nrow(DD))]
  DD <- melt(DD, id = "t")
  DD[, intervention := "no"]
  for (i in 1:nz) {
    if (length(itz[[i]]) > 0) {
      DD[
        variable == paste("Zone", i) & t %in% itz[[i]],
        intervention := "yes"
      ]
    }
  }
  DD
}


## make a list of indices for TB deaths in each zone
getIndexList <- function(args) {
  nz <- args$patch_dims
  cnmz <- getnewcolnames(args)
  idz <- cnmz$index$TB_deaths # column ids of deaths
  zdz <- rep(1:args$patch_dims, 3 * 3) # zones x (HIVxage)
  wcz <- list() # zone key
  for (i in 1:nz) wcz[[i]] <- idz[which(zdz == i)] # list key to zone death cols
  wcz
}


getIndexListArray <- function(args, var33 = "TB_deaths") {
  nz <- args$patch_dims
  cnmz <- getnewcolnames(args)
  idz <- cnmz$index[[var33]] # column ids of deaths
  zdz <- rep(1:args$patch_dims, 3 * 3) # zones x (HIVxage)
  wcz <- list() # zone key
  for (i in 1:nz) wcz[[i]] <- idz[which(zdz == i)] # list key to zone death cols
  wcz
}


getIndexListVector <- function(args, varnp = "PrevByPatch") {
  nz <- length(args$popinit)
  ## get updated column data: see also cnmz$dim$TB_deaths
  cnmz <- getnewcolnames(args)
  idz <- cnmz$index[[varnp]] # column ids of deaths
  idz
}


## calculates weighted means
aggregator <- function(N, k, x, w = 1) {
  span <- N / k
  if (N %% k) stop("k must divide N!")
  id <- rep(1:k, each = span)
  D <- data.table(x, w, id = id)
  D <- D[, .(M = weighted.mean(x, w)), by = id]
  D$M
}

