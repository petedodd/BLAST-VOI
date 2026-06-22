## exploring and justifying flux approximations
library(here)
library(data.table)
library(expm)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(sf)
library(BLASTtbmod)
library(paletteer)

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
for (i in seq_len(nrow(Y))) {
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
source(here("R/utils/flux_utils.R"))

## === trying different data for
zones <- BLASTtbmod::blantyre
pops <- zones$population #populations
prev <- zones$rate_15
plot(zones$population, zones$rate_15, pch = as.character(zones$zone))


## === phyloflows output
## proportion of all transmission that is:
## from row/patch
## to col/patch
FM <- read.csv(here("data/blantyre_flow_april2026.csv"))
FMse <- read.csv(here("data/blantyre_SE_april2026.csv"))
FM <- FM[, 2:ncol(FM)]
FM <- as.matrix(FM)
FMse <- FMse[, 2:ncol(FMse)]
FMse <- as.matrix(FMse)
FMse <- FMse / sum(FM) # NOTE ignores correlation?
FM <- FM / sum(FM) # correct rounding

## initial prevalences
dinit <- cbind(rep(0, 7), prev, prev) * 1e-5
## mixing from flux
## (transpose because rows are from, cols are to)
MM <- t(FM)
for (i in seq_len(nrow(MM))) { #loops for transparency
  for (j in seq_len(ncol(MM))) {
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


## run model
set.seed(1234)
output <- run.model(
  args,
  0:end,
  n.particles = n_particles,
  convert = FALSE
)

## comparison
flx <- get_fluxes(output)

CF <- data.table(
  genomic.flux = c(t(FM) / sum(FM)), #transpose due to different convention
  genomic.flux.sd = c(t(FMse) / sum(FM)), #transpose due to different convention
  output.FOI.flux = c(flx$inf_flux),
  output.FOI.flux.sd = c(flx$inf_flux_sd),
  output.note.flux = c(flx$note_flux),
  output.note.flux.sd = c(flx$note_flux_sd)
)
CF[, variable := gsub("cum_inf_flux", "", nmz)]
CF[, `recipient zone` := gsub("\\[", "", gsub(",.\\]", "", variable))]
CF[, `source zone` := gsub("\\]", "", gsub("\\[.,", "", variable))]


GP <- make_flux_compare_graph(CF)

ggsave(GP,filename = here("output/x_fluxapp_model.png"), w = 10, h = 5)


## === optimization over fluxes
xx <- rep(0.5, 49)
maxit <- 1e3
tol <- 1e-3
k <- 1
err <- 1
while (k < maxit && err > tol) {
  if(k %% 10 == 0) {
    cat("Iteration", k, "Error", err, "\n")
  }
  ratios <- get_ratio_from_x(xx)
  xx <- xx / ratios
  err <- sum((ratios - 1)^2)
  k <- k + 1
}
xx
ratios

## corresponding mixing
MMo <- make_MM_from_x(xx)

plot(MMo / sum(MMo), MM / sum(MM))
abline(a = 0, b = 1, col = 2)


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

argso <- copy(args)
argso$MM <- MMo


## run model
set.seed(1234)
outputo <- run.model(
  argso,
  0:end,
  n.particles = n_particles,
  convert = FALSE
)

## comparison
flxo <- get_fluxes(outputo)

CF <- data.table(
  genomic.flux = c(t(FM) / sum(FM)), # transpose due to different convention
  genomic.flux.sd = c(t(FMse) / sum(FM)),
  output.FOI.flux = c(flxo$inf_flux),
  output.FOI.flux.sd = c(flxo$inf_flux_sd),
  output.note.flux = c(flxo$note_flux),
  output.note.flux.sd = c(flxo$note_flux_sd)
)
CF[, variable := gsub("cum_inf_flux", "", nmz)]
CF[, `recipient zone` := gsub("\\[", "", gsub(",.\\]", "", variable))]
CF[, `source zone` := gsub("\\]", "", gsub("\\[.,", "", variable))]


GP <- make_flux_compare_graph(CF)

ggsave(GP, filename = here("output/x_fluxapp_model_opt.png"), w = 10, h = 5)


## ==== exploring over different mixing patterns impact on fluxes
## inf vs note

## loop this as simple PSA
n_sims <- 1000
set.seed(1234)
results <- list()
## loop over simulations
for (n in 1:n_sims) {
  if (!n %% 10) cat("Running simulation", n, "\n")
  results[[n]] <- get_single_results(n, n_particles = 1) #random MM
}
results <- rbindlist(results)
results
shft <- mean(results$output.note.flux)

smy <- results[,
  .(
    ratio = mean((shft + output.note.flux) / (shft + output.FOI.flux))
  ),
  by = .(variable, `source zone`, `recipient zone`)
]
smy[, summary(ratio)]

smy[, source_pop := pops[as.integer(`source zone`)]]
smy[, receipt_pop := pops[as.integer(`recipient zone`)]]
smy[, source_prev := prev[as.integer(`source zone`)]]
smy[, receipt_prev := prev[as.integer(`recipient zone`)]]


## ----
A <- ggplot(smy, aes(receipt_pop, ratio, group = receipt_pop)) +
  geom_boxplot() +
  theme_classic() +
  ggpubr::grids() +
  scale_x_continuous(label = scales::comma) +
  geom_hline(yintercept = 1, col = 2, lty = 2) +
  labs(
    x = "Recipient population",
    y = "Offset ratio of notification flux to FOI flux"
  )

B <- ggplot(smy, aes(source_pop, ratio, group = source_pop)) +
  geom_boxplot() +
  theme_classic() +
  ggpubr::grids() +
  scale_x_continuous(label = scales::comma) +
  geom_hline(yintercept = 1, col = 2, lty = 2) +
  labs(
    x = "Source population",
    y = "Offset ratio of notification flux to FOI flux"
  )


GP <- ggarrange(
  A, B,
  nrow = 1, ncol = 2,
  common.legend = TRUE,
  labels = c("A", "B")
)
GP

ggsave(GP, filename = here("output/x_fluxapp_note_v_FOI.png"), w = 10, h = 5)
