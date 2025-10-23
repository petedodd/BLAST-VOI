## Adapting to look at N-patch simulations
library(here)
library(data.table)
library(ggplot2)
library(epiR)
library(glue)

## === trying different data for
fn <- here("data/pops.Rdata")
if (!file.exists(fn)) {
  zones <- sf::read_sf(
    here("data/blantyre_7zone_update/blantyre_7zone_update.shp")
  )
  zones <- zones[order(zones$zone), ] # order by zone number
  pops <- zones$population # populations
  save(pops, file = fn)
} else {
  load(fn)
}


## === utils
source(here("R/utils/benefit.R")) #for benefit()

## for tables etc
gr <- function(m, l, h) {
  m <- formatC(m, format = "f", digits = 2)
  l <- formatC(l, format = "f", digits = 2)
  h <- formatC(h, format = "f", digits = 2)
  glue("{m} ({l} to {h})")
}


## ====================================================
## ** task 3 **: Analyse strategies & VOI with emulator

## TODO
## these are from task 2

## patch veresion
np <- c(
  5,
  10,
  15,
  20,
  25,
  30,
  35,
  40,
  45,
  50
)

RL <- PL <- list()
for (n in np) {
  load(glue("~/Dropbox/Holocron/tmp/RP2/R{n}.Rdata"))
  load(glue("~/Dropbox/Holocron/tmp/RP2/P{n}.Rdata"))
  RL[[n]] <- R
  PL[[n]] <- P
}
R <- rbindlist(RL)
P <- rbindlist(PL)


## check alignment
R[, .N, by = .(iter, numpatches)]
P[, .N, by = .(iter, numpatches)]
all(R[, .N, by = .(iter, numpatches)] == P[, .N, by = .(iter, numpatches)])


RP <- R

## NOTE new conversion of slopes to have harmonized screenrate
## screenrate = 0.5 * screenrate0 / (num_patches/10)
## to get slope under screenrate0 -> x (num_patches/10) / 0.5 = num_patches/5
RP[, slp := slp * numpatches / 5]


## difference in benefits (DbenefitFromSlps in utils/benefit.R)

## tests
true_slopes <- RP[iter == 1 & numpatches == 10, slp] / 1e4
pops <- rep(500e3 / 10, 10)
DbenefitFromSlps(true_slopes, pops, # use pops as ranking
  pops, 0.5 * sum(pops),
  verbose = TRUE,
  separate = TRUE
)


## loop at half coverage:
B <- RP[,
  {
    DbenefitFromSlps(pmax(0, slp) / 1e4,
      ## rep(1,numpatches),
      1:numpatches,
      rep(500e3 / numpatches, numpatches),
      250e3,
      separate = TRUE
    ) # using com id as rank
  },
  by = .(iter, numpatches)
]


PB <- merge(P, B, by = c("iter", "numpatches"))
save(PB, file = here("tmpdata/PB.Rdata"))

## still issue with 35?! TODO
ggplot(PB, aes(numpatches, bnft_opt, group = numpatches)) +
  geom_boxplot()


BRO <- PB[, .(real.Fprev.mean, real.Fprev.CV, real.foi.mean, real.foi.CV, alph,
  numpatches,
  DB = bnft_opt - bnft_sub
)]


tab <- as.data.table(epi.prcc(BRO))
tab[, txt := gr(est, lower, upper)]
tab

tabr <- tab[rev(order(abs(est)))][1:5][, .(variable = var, PRCC = txt)]
tabr
fwrite(tabr,file=here("output/VOI_psa_prccN.csv"))

BRO2 <- BRO[numpatches == 50]
BRO2[, numpatches := NULL]

tab <- as.data.table(epi.prcc(BRO2))
tab[, txt := gr(est, lower, upper)]
tab

tabr <- tab[rev(order(abs(est)))][1:5][, .(variable = var, PRCC = txt)]
tabr
fwrite(tabr, file = here("output/VOI_psa_prccNR.csv"))
