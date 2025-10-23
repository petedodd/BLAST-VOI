## this will create a set of 50 populations for each parameter
library(data.table)
library(here)
library(lhs)

## tmpdata/

(tot <- 7 * 5^2 * 3^2 * 2^3) #smallest integer divisible by all of below
exx <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
tot / exx

num_patches <- tot #fixed

## --- top level stuff
## prevalence
prev.mean <- 200e-5
prev.CV <- sqrt(num_patches / 10) # 1 for npatch=10, growing with sqrt npatches
prev.mu <- log(prev.mean / sqrt(1 + prev.CV^2))
prev.sg <- sqrt(log(1 + prev.CV^2))
## curve(dlnorm(x, prev.mu, prev.sg), from = 0, to = 1e-2, n = 1e3)

## initial state
ari0.mu <- log(0.05)
## sqrt(log(2)) for npatch=10, CV growing with sqrt npatches
ari0.sg <- sqrt(log(1 + num_patches / 10))
## curve(dlnorm(x, ari0.mu, ari0.sg), from = 0, to = 0.1, n = 1e3)

## transmission
beta.mu <- log(10)
## sqrt(log(2)) for npatch=10, CV growing with sqrt npatches
beta.sg <- sqrt(log(1+num_patches/10))
##curve(dlnorm(x, beta.mu, beta.sg), from = 0, to = 20, n = 1e3)
alph.a <- 2
alph.b <- 2

## progression
mn <- 0.05860121
pf.mu <- log(mn / (sqrt(1 + 1.5)))
## sqrt(log(1.5)) for npatch=10, CV growing with sqrt npatches
pf.sg <- sqrt(log(1 + 0.5 * num_patches / 10))
## curve(dlnorm(x, pf.mu, pf.sg), from = 0, to = 0.1, n = 1e3)

## ================== PSA
## --- sample loop
set.seed(123)
plen <- 40
prvmns <- seq(from = 10, to = 700, len = plen) * 1e-5
prvCVs <- seq(from = 0.2, 5, len = plen)
G <- expand.grid(prvmns, prvCVs)
NK <- nrow(G)
LH <- randomLHS(NK, 4)
PARGS <- list()
for (k in 1:NK) {
  PARGS[[k]] <- list() # list of lists

  ## prevalence
  prev.mean <- G[k, 1]
  prev.CV0 <- G[k, 2]
  prev.CV <- prev.CV0 * sqrt(num_patches / 10)
  prev.mu <- log(prev.mean / sqrt(1 + prev.CV^2))
  prev.sg <- sqrt(log(1 + prev.CV^2))
  prevz <- rlnorm(num_patches, meanlog = prev.mu, sdlog = prev.sg)
  PARGS[[k]][["prevz"]] <- prevz

  ## initial state
  ari0.mean <- 0.05
  ari0.CV <- 1
  ari0.mu <- log(ari0.mean / sqrt(1 + ari0.CV^2))
  ari0.sg <- sqrt(log(1 + ari0.CV^2))
  ari0 <- qlnorm(LH[k, 1], ari0.mu, ari0.sg)
  PARGS[[k]][["ari0"]] <- ari0

  ## transmission
  beta.mu <- log(10)
  beta.sg <- sqrt(log(2))
  bet <- qlnorm(LH[k, 2], beta.mu, beta.sg)
  RRbeta <- rlnorm(num_patches, meanlog = 0, sdlog = prev.sg)
  alph <- qbeta(LH[k, 3], alph.a, alph.b)
  PARGS[[k]][["bet"]] <- bet
  PARGS[[k]][["alph"]] <- alph
  PARGS[[k]][["RRbeta"]] <- RRbeta

  ## progression
  mn <- 0.05860121
  pf.mu <- log(mn / (sqrt(1 + 1.5)))
  pf.sg <- sqrt(log(1.5)) # sqrt(log(1.5))
  pf <- qlnorm(LH[k, 4], pf.mu, pf.sg)
  PARGS[[k]][["pf"]] <- RRbeta
} # end k-loop


save(PARGS, file = here("tmpdata/PARGS.Rdata"))


