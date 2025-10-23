## implementing an aggregation over patches

## get args
## R --slave --vanilla --args < 04.calculations.R npatches
args <- commandArgs(trailingOnly = TRUE) #get argsm
print(args)

npatches <- 35## 50
totN <- 7 * 5^2 * 3^2 * 2^3

npatches <- as.integer(args[1])

cat("Using npatches = ", npatches, "...\n")


## NOTE this version is looking at varying N also
## ====================================================
## ** task 1 **: Set up model and check if slopes apply

## trying to create something to test interventions:
library(BLASTtbmod)
library(ggplot2)
library(data.table)
library(tictoc)
library(here)
library(glue)

gh <- function(x) glue(here(x))

## utility functions
source(here("R/utils/PSA_utils.R"))

## PSA data for 50 comms
load(here("tmpdata/PARGS.Rdata"))

## baseparms
args <- get.parms(start_year = 2015, years = 6)


## --- change patch no
num_patches <- npatches
simlength <- 73

## check
args <- updateARGS(args, numpatches = num_patches, simlength = 74)
test0 <- run.model(args, args$tt, n.particles = 5)

## strength NOTE needed in analysis also
screenrate0 <- 1e4 / args$dt # 10K per months as people per year
## tone down ACF to avoid saturating small patches
screenrate <- 0.5 * screenrate0 / (num_patches/10)

## check
itz <- makeITZ(args$popinit,
  screenrate = screenrate,
  burnin = 12
) # each one is >=5 long
args <- updateARGS(args,
  numpatches = num_patches,
  simlength = 1 + max(itz[[num_patches]])
)
indexlist <- getIndexList(args) #TB deaths
indexlistPrev <- getIndexListVector(args, varnp = "PrevByPatch")
indexlistCFOI <- getIndexListVector(args, varnp = "cum_inf_ByPatch")


## ====================================================
## ** task 2 **: Generate data & emulate model
## see other work:

## ================== PSA
## --- sample loop
Nreps <- 5
R <- P <- list() #results containers
NK <- length(PARGS)
tic()
for (k in 1:NK) {
  if (!k %% 10) cat("k = ", k, "\n")

  ## prevalence
  prevz <- aggregator(totN, num_patches, PARGS[[k]]$prevz)

  ## initial state
  ari0 <- PARGS[[k]]$ari0

  ## transmission
  bet <- PARGS[[k]]$bet
  alph <- PARGS[[k]]$alph
  RRbeta <- aggregator(totN, num_patches, PARGS[[k]]$RRbeta)

  ## test
  MM <- diag(1 / num_patches, num_patches) * alph +
    (1 - alph) * matrix(1 / num_patches, num_patches, num_patches)
  MM <- diag(RRbeta) %*% MM # heterogeneity in beta

  ## progression
  pf <- aggregator(totN, num_patches, PARGS[[k]]$pf)

  ## run
  args$initD[, 2:3] <- prevz
  args$ari0 <- ari0
  args$pf <- pf
  args$beta <- bet
  args$MM <- MM

  args0 <- args <- setACF(args, screenrate, itz) # ACF
  args0$ACFhaz0[, ] <- args0$ACFhaz1[, ] <- 0 # no ACF
  YZ <- run.CF(args0, args, args$tt, n.particles = Nreps)
  Y0 <- YZ[[1]]
  Y <- YZ[[2]]

  ## apply:
  DD <- getMeanDD(Y, Y0, itz, indexlist)

  ## extract other metrics
  PrevEnds <- rowMeans(Y0[indexlistPrev, , dim(Y0)[3]])
  FoiEnds <- rowMeans(Y0[indexlistCFOI, , dim(Y0)[3]]) /
    (args$popinit * dim(Y0)[3] / 12)

  ## record
  R[[k]] <- DD[intervention == "yes",
    {
      list(slp = coef(lm(data = .SD, value ~ t))[2])
    },
    by = variable
  ] # slopes
  R[[k]][, c("prev", "iter", "numpatches") := .(prevz, k, num_patches)]
  P[[k]] <- data.table(
    iter = k,
    real.prev.mean = mean(prevz),
    real.prev.CV = sd(prevz) / mean(prevz),
    real.Fprev.mean = mean(PrevEnds),
    real.Fprev.CV = sd(PrevEnds) / mean(PrevEnds),
    real.foi.mean = mean(FoiEnds),
    real.foi.CV = sd(FoiEnds) / mean(FoiEnds),
    ari0 = ari0,
    pf = pf,
    alph = alph,
    bet = bet,
    numpatches = num_patches,
    RRbetaSD = sd(RRbeta),
    RRbeta.CV = sd(RRbeta) / mean(RRbeta)
  )
} # end k-loop
toc()

R <- rbindlist(R)
P <- rbindlist(P)

## save out
save(R, file = gh("tmpdata/R{num_patches}.Rdata"))
save(P, file = gh("tmpdata/P{num_patches}.Rdata"))
