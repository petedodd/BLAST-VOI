## exploring and justifying flux approximations
library(here)
library(data.table)
library(expm)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(readxl)
library(sf)
library(BLASTtbmod)


## ========= PART 1: analytical approximations to dynamics

## example parameters
sigma <- 0.62 # stabilisation
alpha <- 0.05 # fast progression
mu <- 1 / 60 # mortality
epsilon <- 1e-3 # slow progn
gamma <- 0.5 # symptom progn
omega <- 1 / 3 # tb cessation
delta <- 2 / 3 # detection rate
rho <- 0.1 # relapse rate
tau <- 2 # 1/treatment dur


## transition matrix
M <- matrix(
  c(
    -(sigma + alpha + mu), 0, 0, 0, 0, 0, # E
    sigma, -(epsilon + mu), 0, 0, 0, 0, # L
    alpha, epsilon, -(gamma + mu), 0, 0, 0, # SC
    0, 0, gamma, -(omega + delta + mu), 0, rho, # CD
    0, 0, 0, delta, -(tau + mu), 0, # T
    0, 0, 0, 0, tau, -(rho + mu) # R
  ),
  nrow = 6, ncol = 6, byrow = TRUE
)
M

## dynamics
nmz <- c("E", "L", "SC", "CD", "T", "R")
tz <- seq(from = 0, to = 20, by = 0.1)
Y <- matrix(nrow = length(tz), ncol = 6)
for (i in 1:nrow(Y)) {
  Y[i, ] <- expm(M * tz[i])[, 1]
}
colnames(Y) <- nmz

## inspect:
matplot(tz, Y, type = "l")
legend("topright", lty = 1:6, col = 1:6, legend = nmz)


matplot(tz, Y[, 3:6], type = "l", lty = 3:6, col = 3:6)
legend("topright", lty = 3:6, col = 3:6, legend = nmz[3:6])


Eapp <- exp(-(sigma + alpha + mu) * tz)
Lapp <- (sigma / (epsilon - sigma - alpha)) *
  (exp(-(sigma + alpha + mu) * tz) - exp(-(mu + epsilon) * tz))

tm <- tau + mu
rr <- rho + mu
om <- omega + delta + mu
sam <- sigma + alpha + mu
gm <- gamma + mu
em <- epsilon + mu
B <- epsilon * sigma / (epsilon - sigma - alpha)
A <- B + alpha
SCapp <- (exp(-sam * tz) - exp(-gm * tz)) * A / (gamma - sigma - alpha) -
  (exp(-em * tz) - exp(-gm * tz)) * B / (gamma - epsilon)

CDapp0 <- gamma * A * (exp(-sam * tz) - exp(-om * tz)) /
  ((gamma - sigma - alpha) * (omega + delta - sigma - alpha)) -
  gamma * B * (exp(-em * tz) - exp(-om * tz)) /
    ((gamma - epsilon) * (omega + delta - epsilon)) +
  gamma * (B / (gamma - epsilon) - A / (gamma - sigma - alpha)) *
    (exp(-gm * tz) - exp(-om * tz)) /
    (omega + delta - gamma)

F <- delta * gamma * A /
  ((gamma - sigma - alpha) * (omega + delta - sigma - alpha))
G <- -delta * gamma * B /
  ((gamma - epsilon) * (omega + delta - epsilon))
H <- delta * gamma *
  (B / (gamma - epsilon) -
   A / (gamma - sigma - alpha)) / (omega + delta - gamma)
J <- delta * gamma * (A * (1 / (omega + delta - gamma) -
  1 / (omega + delta - sigma -
    alpha)) / (gamma - sigma - alpha) +
  B * (1 / (omega + delta - epsilon) -
    1 / (omega + delta - gamma)) / (gamma - epsilon))
K <- -F / (tau - sigma - alpha) -
  G / (tau - epsilon) -
  H / (tau - gamma) - J / (tau - omega - delta)


Tapp0 <- F * exp(-sam * tz) +
  G * exp(-em * tz) +
  H * exp(-gm * tz) +
  J * exp(-om * tz)


plot(Tapp0)
## Tapp0 ~ CDapp0
plot(CDapp0, Tapp0)
abline(a = 0, b = 1, col = 2)

## presume because tau large enough to make it like delta convolution
Rapp0 <- tau * F * (exp(-sam * tz) - exp(-rr * tz)) /
  ((rho - sigma - alpha) * (tau - sigma - alpha)) +
  tau * G * (exp(-em * tz) - exp(-rr * tz)) /
    ((rho - epsilon) * (tau - epsilon)) +
  tau * H * (exp(-gm * tz) - exp(-rr * tz)) /
    ((rho - gamma) * (tau - gamma)) +
  tau * J * (exp(-om * tz) - exp(-rr * tz)) /
    ((rho - omega - delta) * (tau - omega - delta)) +
  tau * K * (exp(-tm * tz) - exp(-rr * tz)) /
    (rho - tau)


Yapp <- cbind(Eapp, Lapp, SCapp, CDapp0, Tapp0, Rapp0)



## inspect
matplot(tz, Y[, 3:6], type = "p", pch = 1, col = 3:6, xlab = "Years", ylab = "")
matplot(tz, Yapp[, 3:6], type = "l", add = TRUE, lty = 3:6, col = 3:6)
legend("topright", fill = 3:6, legend = nmz[3:6])

## 1 loop corx
corx <- rho * tau * (
  F * ((exp(-sam * tz) - exp(-om * tz)) / (omega + delta - sigma - alpha) -
    (exp(-rr * tz) - exp(-om * tz)) / (omega + delta - rho)) /
    ((rho - sigma - alpha) * (tau - sigma - alpha)) +
    G * ((exp(-em * tz) - exp(-om * tz)) / (omega + delta - epsilon) -
      (exp(-rr * tz) - exp(-om * tz)) / (omega + delta - rho)) /
      ((rho - epsilon) * (tau - epsilon)) +
    H * ((exp(-gm * tz) - exp(-om * tz)) / (omega + delta - gamma) -
      (exp(-rr * tz) - exp(-om * tz)) / (omega + delta - rho)) /
      ((rho - gamma) * (tau - gamma)) +
    J * ((tz * exp(-om * tz)) - (exp(-rr * tz) -
      exp(-om * tz)) / (omega + delta - rho)) /
      ((rho - omega - delta) * (tau - omega - delta)) +
    K * ((exp(-tm * tz) - exp(-om * tz)) / (omega + delta - tau) -
      (exp(-rr * tz) - exp(-om * tz)) / (omega + delta - rho)) /
      (rho - tau)
)

Yapp1 <- cbind(Eapp, Lapp, SCapp, CDapp0 + corx, Tapp0, Rapp0)


matplot(tz, Y[, 3:6], type = "p", pch = 1, col = 3:6, xlab = "Years", ylab = "")
matplot(tz, Yapp1[, 3:6], type = "l", add = TRUE, lty = 3:6, col = 3:6)
legend("topright", fill = 3:6, legend = nmz[3:6])

## neater plot:
dta <- data.table(
  time = rep(tz, 4),
  value = c(Y[, 3:6]),
  quantity = rep(nmz[3:6], each = nrow(Y)),
  type = c("truth")
)

app <- data.table(
  time = rep(tz, 8),
  value = c(Yapp[, 3:6], Yapp1[, 3:6]),
  quantity = rep(rep(nmz[3:6], each = nrow(Y), 2)),
  type = c(rep("1st order", 4 * length(tz)), rep("2nd order", 4 * length(tz)))
)

ggplot(dta, aes(time, value,
  col = quantity, shape = quantity, lty = type
)) +
  geom_point() +
  geom_line(data = app) +
  scale_shape(solid = FALSE) +
  theme_linedraw()

## save out plot for relevant quantity
gp <- ggplot(dta[quantity == "CD"], aes(time, value,
  col = quantity, shape = quantity, lty = type
)) +
  geom_point() +
  geom_line(data = app[quantity == "CD"]) +
  scale_shape(solid = FALSE) +
  theme_linedraw() +
  theme(legend.position = "none") +
  xlab("Years") +
  ylab("Value")
gp

ggsave(gp, filename = here("output/x_fluxapp.png"), w = 7, h = 5)


## ========= PART 2: fluxes as inputs
## === trying different data for
zones <- read_sf(here("data/shp/blantyre_7zone_update.shp"))
zones <- zones[order(zones$zone), ] #order by zone number
pops <- zones$population #populations
prev <- zones$rate_15
plot(zones$tb_pre, zones$rate_15)
abline(a = 0, b = 1, col = 2)

## === phyloflows output
## proportion of all transmission that is:
## from row/patch
## to col/patch
FM <- read_excel(here("data/PhyloFlows_7Zone.xlsx"))
FM <- FM[, 2:ncol(FM)]
FM <- as.matrix(FM)
FM <- FM / sum(FM) # correct rounding


## initial prevalences
dinit <- cbind(rep(0, 7), prev, prev) * 1e-5
## mixing from flux
## (transpose because rows are from, cols are to)
MM <- t(FM)
for (i in 1:nrow(MM)) { #loops for transparency
  for (j in 1:nrow(MM)) {
    MM[i, j] <- MM[i, j] / (pops[i] * prev[j])
  }
}
MM <- MM / max(Re(eigen(MM)$values)) #renormalize

## create arguments
n_particles <- 100
start_year <- 2015
years <- 6
end <- years * 12
args <- get.parms(
  start_year = start_year,
  years = years, Dinit = dinit, ari0 = 5e-3
)
args$beta <- 20
args$cdr <- 0.7
args$MM <- MM
args$popinit <- pops #CHECK correct order in data/function

## run model
set.seed(1234)
output <- run.model(
  args,
  0:end,
  n.particles = n_particles,
  convert = FALSE
)

## get notification flux
note_flux_idx <- grep("cum_note_flux", BLASTtbmod::get_cols)
(nmz <- grep("cum_note_flux", BLASTtbmod::get_cols, value = TRUE)) # check
tmp <- t(apply(
  output[note_flux_idx, , ],
  MARGIN = c(1, 3), FUN = mean
)) # matrix over time, mean over particles
tmp_sd <- t(apply(
  output[note_flux_idx, , ],
  MARGIN = c(1, 3), FUN = sd
)) # matrix over time, SD over particles
matplot(tmp, type = "l")           #plot
note_flux <- matrix(
  tail(tmp, n = 1), # last for cumulative
  ncol = args$patch_dims,
  nrow = args$patch_dims
)
note_flux_sd <- matrix(
  tail(tmp_sd, n = 1), # last for cumulative
  ncol = args$patch_dims,
  nrow = args$patch_dims
)
note_flux_sd <- note_flux_sd / max(note_flux) #normalize
note_flux <- note_flux / sum(note_flux) #normalize

## get infection flux
inf_flux_idx <- grep("cum_inf_flux", BLASTtbmod::get_cols)
(nmz <- grep("cum_inf_flux", BLASTtbmod::get_cols, value = TRUE)) # check
tmp <- t(apply(
  output[inf_flux_idx, , ],
  MARGIN = c(1, 3), FUN = mean
)) # matrix over time, mean over particles
tmp_sd <- t(apply(
  output[inf_flux_idx, , ],
  MARGIN = c(1, 3), FUN = sd
)) # matrix over time, SD over particles
matplot(tmp, type = "l")          #plot
inf_flux <- matrix(
  tail(tmp, n = 1),
  ncol = args$patch_dims, nrow = args$patch_dims
)
inf_flux_sd <- matrix(
  tail(tmp_sd, n = 1),
  ncol = args$patch_dims, nrow = args$patch_dims
)
inf_flux_sd <- inf_flux_sd / max(inf_flux) # normalize
inf_flux <- inf_flux / sum(inf_flux) #normalize

## comparison
CF <- data.table(
  genomic.flux = c(t(FM) / max(FM)), #transpose due to different convention
  output.FOI.flux = c(inf_flux),
  output.FOI.flux.sd = c(inf_flux_sd),
  output.note.flux = c(note_flux),
  output.note.flux.sd = c(note_flux_sd)
)
CF[, variable := gsub("cum_inf_flux", "", nmz)]
CF[, `recipient zone` := gsub("\\[", "", gsub(",.\\]", "", variable))]
CF[, `source zone` := gsub("\\]", "", gsub("\\[.,", "", variable))]

## FOI vs note
GP1 <- ggplot(
  CF,
  aes(
    x = output.FOI.flux,
    xmin = output.FOI.flux - 1.96 * output.FOI.flux.sd,
    xmax = output.FOI.flux + 1.96 * output.FOI.flux.sd,
    y = output.note.flux,
    ymin = output.note.flux - 1.96 * output.note.flux.sd,
    ymax = output.note.flux + 1.96 * output.note.flux.sd,
    label = variable, shape = `source zone`, col = `source zone`
  )
) +
  geom_abline(intercept = 0, slope = 1, col = 2, lty = 2) +
  scale_shape_manual(values = 0:6) +
  geom_errorbar(width = 1e-2, alpha = 0.4) +
  geom_errorbarh(width = 1e-2, alpha = 0.4) +
  geom_point(size = 2) +
  geom_text_repel(show.legend = FALSE) +
  theme_classic() +
  ggpubr::grids() +
  xlab("Model infection flux") +
  ylab("Model notification flux") +
  theme(legend.position = "top") +
  guides(shape = guide_legend(nrow = 1, byrow = TRUE))



## data vs note
GP2 <- ggplot(
  CF,
  aes(
    x = genomic.flux,
    y = output.note.flux,
    ymin = output.note.flux - 1.96 * output.note.flux.sd,
    ymax = output.note.flux + 1.96 * output.note.flux.sd,
    label = variable, shape = `source zone`, col = `source zone`
  )
) +
  geom_abline(intercept = 0, slope = 1, col = 2, lty = 2) +
  scale_shape_manual(values = 0:6) +
  geom_errorbar(width = 1e-2, alpha = 0.4) +
  geom_point(size = 2) +
  geom_text_repel(show.legend = FALSE) +
  theme_classic() +
  ggpubr::grids() +
  xlab("Genomic-derived flux") +
  ylab("Model notification flux") +
  theme(legend.position = "top") +
  guides(shape = guide_legend(nrow = 1, byrow = TRUE))
GP2

## combine plot
GP <- ggarrange(
  GP1, GP2,
  nrow = 1, ncol = 2,
  common.legend = TRUE,
  labels = c("A", "B")
)

ggsave(GP,filename = here("output/x_fluxapp_model.png"), w = 10, h = 5)


## ==== exploring over different mixing patterns
get_single_results <- function(n) {
  ## create random mixing matrix
  MM <- FM <- matrix(runif(49), nrow = 7)
  ## sample as assortativity + random
  ## FM <- runif(1) * FM / 7
  ## diag(FM) <- runif(7)
  FM <- FM / sum(FM) # normalize
  for (i in 1:nrow(MM)) { # loops for transparency
    for (j in 1:nrow(MM)) {
      MM[i, j] <- FM[i, j] / (pops[i] * prev[j])
    }
  }
  MM <- MM / max(Re(eigen(MM)$values)) # renormalize
  ## modify arguments
  args$MM <- MM

  ## run model
  output <- run.model(
    args,
    0:end,
    n.particles = n_particles,
    convert = FALSE
  )

  ## get infection flux
  tmp <- t(apply(
    output[inf_flux_idx, , ],
    MARGIN = c(1, 3), FUN = mean
  )) # matrix over time, mean over particles
  inf_flux <- matrix(
    tail(tmp, n = 1),
    ncol = args$patch_dims, nrow = args$patch_dims
  )
  inf_flux <- inf_flux / sum(inf_flux) # normalize

  ## get notification flux
  tmp <- t(apply(
    output[note_flux_idx, , ],
    MARGIN = c(1, 3), FUN = mean
  )) # matrix over time, mean over particles
  note_flux <- matrix(
    tail(tmp, n = 1), # last for cumulative
    ncol = args$patch_dims, nrow = args$patch_dims
  )
  note_flux <- note_flux / sum(note_flux) # normalize

  ## comparison
  CF <- data.table(
    iter = n,
    genomic.flux = c(t(FM) / sum(FM)), # transpose
    output.FOI.flux = c(inf_flux),
    output.note.flux = c(note_flux)
  )
  CF[, variable := gsub("cum_inf_flux", "", nmz)]
  CF[, `recipient zone` := gsub("\\[", "", gsub(",.\\]", "", variable))]
  CF[, `source zone` := gsub("\\]", "", gsub("\\[.,", "", variable))]

  return(CF)
}


## loop this as simple PSA
n_sims <- 100
set.seed(1234)
results <- list()
## loop over simulations
for (n in 1:n_sims) {
  cat("Running simulation", n, "\n")
  results[[n]] <- get_single_results(n)
}
results <- rbindlist(results)

## summarize over iterations
smy <- results[,
  .(
    truth = mean(genomic.flux),
    infection = mean(output.FOI.flux),
    notification = mean(output.note.flux),
    truth.lo = quantile(genomic.flux, 0.025),
    infection.lo = quantile(output.FOI.flux, 0.025),
    notification.lo = quantile(output.note.flux, 0.025),
    truth.hi = quantile(genomic.flux, 1 - 0.025),
    infection.hi = quantile(output.FOI.flux, 1 - 0.025),
    notification.hi = quantile(output.note.flux, 1 - 0.025)
  ),
  by = .(variable, `source zone`, `recipient zone`)
]


## plot data vs note
GP <- ggplot(
  smy,
  aes(
    x = truth, xmin = truth.lo, xmax = truth.hi,
    y = notification, ymin = notification.lo, ymax = notification.hi,
    label = variable, shape = `source zone`, col = `source zone`
  )
) +
  geom_abline(intercept = 0, slope = 1, col = 2, lty = 2) +
  scale_shape_manual(values = 0:6) +
  geom_errorbar(width = 1e-2, alpha = 0.4) +
  geom_errorbar(width = 1e-2, alpha = 0.4, orientation = "y") +
  geom_point(size = 2) +
  geom_text_repel(show.legend = FALSE) +
  theme_classic() +
  ggpubr::grids() +
  ylab("Mixing flux") +
  xlab("Model notification flux") +
  theme(legend.position = "top") +
  guides(shape = guide_legend(nrow = 1, byrow = TRUE)) +
  coord_fixed()
GP
