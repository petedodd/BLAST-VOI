## exploring and justifying flux approximations
library(expm)
library(ggplot2)
library(data.table)
library(here)

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


## ========= PART 2:


## TODO simulation of this and comparison with mixing
rm(list = ls())

library(here)
library(BLASTtbmod)
library(data.table)
library(readxl)
library(ggplot2)
library(ggrepel)
library(sf)

## === trying different data for

zones <- read_sf(here("data/shp/blantyre_7zone_update.shp"))
zones <- zones[order(zones$zone), ] #order by zone number
pops <- zones$population #populations
prev <- zones$rate_15
plot(zones$tb_pre, zones$rate_15)
abline(a = 0, b = 1, col = 2)

## === phyloflows output
FM <- read_excel(here("data/PhyloFlows_7Zone.xlsx"))
FM <- FM[, 2:ncol(FM)]
FM <- as.matrix(FM)

## === running the model

## define other parameters
n_particles <- 10
time <- 1 # NOTE not used TODO check
start_year <- 2015
years <- 6

## NOTE on prevalence for initial state
## prev is split across SC and D in 1:1 ratio - so need to double expectation (see init_D/SC)
## is only applied to (1-exp(-ari0*ageMids)) of pop, so need to divide by this
## PPP <- fread(here("data/PPP.csv"))
## PPP$ab <- PPP$posterior.a + PPP$posterior.b
## PPP$prev.mid <- PPP$posterior.a / PPP$ab
## PPP$patch <- paste("Patch", PPP$zone)
ageMids <- c(7.5, 32.5, 60) #get from get.parms CHECK used?
ariparm <- 5e-2 # high?

## plot(prev)


## ## TESTING
## prev <- rep(mean(prev),7) #uniform
## prev[2:7] <- 0

## initial prevalences
dinit <- cbind(rep(0, 7), prev, prev) * 1e-5
## mixing from flux
MM <- t(FM) ## use rescaled flux
## (transpose because rows are from, cols are to)
for (i in 1:nrow(MM)) { #loops for transparency
  for (j in 1:nrow(MM)) {
    ## MM[i, j] <- MM[i, j] / (pops[i] * prev[j])
    MM[i, j] <- pops[j] #this is random mixing
    ## MM[i, j] <- 1 #should give fixed lam_ij
  }
}
MM <- MM / max(Re(eigen(MM)$values)) #renormalize

## create arguments
args <- get.parms(start_year = start_year, years = years, Dinit = dinit, ari0 = ariparm/10)
args$beta <- 20
args$cdr <- 0.7
args$MM <- MM
args$popinit <- pops #CHECK correct order in data/function


 ## ## ## mixing matrix from FM
## ## mean(FM) / mean(args$MM)
## ## max(Re(eigen(FM)$values)) / max(Re(eigen(args$MM)$values))
## ## args$MM <- (max(Re(eigen(args$MM)$values)) / max(Re(eigen(FM)$values))) * t(FM) ## use rescaled flux
## ## (transpose because rows are from, cols are to)


## === new comparisons ===
end <- 60
test <- run.model(
  args,
  0:end,
  n.particles = 5,
  convert = FALSE
)


## check infection flux:
note_flux_idx <- grep("cum_note_flux", BLASTtbmod::get_cols)
(nmz <- grep("cum_note_flux", BLASTtbmod::get_cols, value = TRUE)) # check
tmp <- t(test[note_flux_idx, 1, ]) # matrix over time for particle 1 TODO
matplot(tmp, type = "l")           #plot
note_flux <- matrix(
  tail(tmp, n = 1), # last for cumulative
  ncol = args$patch_dims,
  nrow = args$patch_dims
)
note_flux <- note_flux / max(note_flux) #normalize

inf_flux_idx <- grep("cum_inf_flux", BLASTtbmod::get_cols)
(nmz <- grep("cum_inf_flux", BLASTtbmod::get_cols, value = TRUE)) # check
tmp <- t(test[inf_flux_idx, 1, ]) # matrix over time for particle 1 TODO
matplot(tmp, type = "l")          #plot
## empirical infection flux
inf_flux <- matrix(
  tail(tmp, n = 1),
  ncol = args$patch_dims, nrow = args$patch_dims
)
inf_flux <- inf_flux / max(inf_flux) #normalize

## comparison
CF <- data.table(
  genomic.flux = c(FM / max(FM)),
  output.FOI.flux = c(inf_flux),
  output.note.flux = c(note_flux)
)
CF[, variable := gsub("cum_inf_flux", "", nmz)]
CF[, fromid := gsub("\\[", "", gsub(",.\\]", "", variable))]


## data vs FOI
ggplot(
  CF,
  aes(
    x = genomic.flux, y = output.FOI.flux,
    label = variable, shape = fromid, col = fromid
  )
) +
  geom_abline(intercept = 0, slope = 1, col = 2) +
  scale_shape_manual(values = 0:6) +
  geom_point(size = 2) +
  geom_text_repel(show.legend = FALSE) +
  theme_classic() +
  ggpubr::grids()

## ggsave(here("plots/flux_cfFOI.png"), w = 7, h = 6)


## data vs note
ggplot(
  CF,
  aes(
    x = genomic.flux, y = output.note.flux,
    label = variable, shape = fromid, col = fromid
  )
) +
  geom_abline(intercept = 0, slope = 1, col = 2) +
  scale_shape_manual(values = 0:6) +
  geom_point(size = 2) +
  geom_text_repel(show.legend = FALSE) +
  theme_classic() +
  ggpubr::grids()

## ggsave(here("plots/flux_cfnote.png"), w = 7, h = 6)

## FOI vs note
ggplot(
  CF,
  aes(
    x = output.FOI.flux, y = output.note.flux,
    label = variable, shape = fromid, col = fromid
  )
) +
  geom_abline(intercept = 0, slope = 1, col = 2) +
  scale_shape_manual(values = 0:6) +
  geom_point(size = 2) +
  geom_text_repel(show.legend = FALSE) +
  theme_classic() +
  ggpubr::grids()

## ggsave(here("plots/flux_notesVfoi.png"), w = 7, h = 6)
1


## (nmz <- grep("Tijk", BLASTtbmod::get_cols, value = TRUE)) # check
## (nmz <- grep("Sij", BLASTtbmod::get_cols, value = TRUE)) # check
## (nmz <- grep("cum_note_flux", BLASTtbmod::get_cols, value = TRUE)) # check
## (nmz <- grep("Sij", BLASTtbmod::get_cols, value = TRUE)) # check
## idx <- grep("Sij", BLASTtbmod::get_cols)
## tmp <- t(test[idx, 1, ]) # matrix over time for particle 1
## matplot(tmp, type = "l")
## str(tmp)


## (nmz <- grep("Tijk", BLASTtbmod::get_cols, value = TRUE)) # check
## idx <- grep("Tijk", BLASTtbmod::get_cols)
## tmp <- t(test[idx, 1, ]) # matrix over time for particle 1
## matplot(tmp, type = "l")
## str(tmp)

## TIJK <- array(tmp[end, ], dim = c(7, 7, 6)) # last time
## TIJ <- matrix(0,7,7)
## for(k in 1:6) TIJ <- TIJ + TIJK[,,k]

## TIJ


## (nmz2 <- grep("cum_inf_flux", BLASTtbmod::get_cols, value = TRUE)) # check
## idx2 <- grep("cum_inf_flux", BLASTtbmod::get_cols)
## tmp2 <- t(test[idx2, 1, ]) # matrix over time for particle 1
## matplot(tmp2, type = "l")
## str(tmp2)

## FF <- matrix(tmp2[end,],7,7) #last time

## plot(
##   TIJ / sum(TIJ),
##   FF / sum(FF)
## )

## abline(a=0,b=1,col=2)

## ## :-)

## ## NOTE see bottom for parameter computations
## ## what about ellij?
## ## OK: calculated with version with cum_note_flux acumulates ellij
## ## OK: noisier with sampling against foitemp but OK
## (nmz3 <- grep("cum_note_flux", BLASTtbmod::get_cols, value = TRUE)) # check
## idx3 <- grep("cum_note_flux", BLASTtbmod::get_cols)
## tmp3 <- t(test[idx3, 1, ]) # matrix over time for particle 1
## matplot(tmp3, type = "l")

## GG <- matrix(tmp3[end,],7,7) #last time

## plot(
##   GG/sum(GG),
##   FF/sum(FF)
## )
## abline(a=0,b=1,col=2)

## ## :-)

## ========= upper checks


## ## run dust model not as particle filter
## test <- run.model(args,
##   0:12, # CHECK start time, times etc
##   n.particles = 50,
##   convert = FALSE
## )

## ## checking ordering (prevalence)
## testo <- extract.pops(test, 1:dim(test)[2], out_type = c("notes", "N", "D", "SC"))
## testo[varname == "SC", varname := "D"] # merge SC and D
## testo <- testo[age != "0-14",
##   .(value = sum(value)),
##   by = .(particle, step, varname, patch)
## ] # aggr/adult
## testo <- dcast(testo, particle+step + patch ~ varname, value.var = "value")
## ## ggplot(testo[], aes(step, 1e5 * D / N,group=particle)) +
## ##   geom_line() +
## ##   facet_wrap(~patch) +
## ##   expand_limits(y = c(0, NA))

## L <- testo[step == 1]
## L[, prevmod := 1e5 * D / N]
## L <- merge(L, data.table(prevdat = prev, patch = paste0("Patch ", 1:7)), by = "patch")
## LM <- L[,.(prevmod=mean(prevmod),prevdat=mean(prevdat)),by=patch]
## ggplot(L, aes(prevdat, prevmod, group = patch)) +
##   geom_boxplot() +
##   geom_point() +
##   geom_text(data=LM,aes(label=patch),col="magenta")+
##   geom_abline(intercept = 0, slope = 1, col = 2)

## ms <- max(testo$step)
## L <- testo[step %in% c(1,ms)]
## L[, prevmod := 1e5 * D / N]
## L <- dcast(L, particle + patch ~ step, value.var = "prevmod")
## names(L)[3:4] <- c("start","end")

## ## ggplot(L, aes(start, end, group = patch,col=patch)) +
## ##   geom_point() +
## ##   stat_ellipse()+
## ##   geom_abline(intercept = 0, slope = 1, col = 2)+
## ##   expand_limits(x=c(0,NA),y=c(0,NA))

## L <- merge(L,
##   data.table(prevdat = prev, patch = paste0("Patch ", 1:7)),
##   by = "patch"
## )

## ## ggplot(L, aes(prevdat, end, group = patch, col = patch, label = patch)) +
## ##   geom_point() +
## ##   geom_text()

## ## check initial states: pops and prev
## testp <- extract.pops(test, 1:dim(test)[2], out_type = c( "N", "D", "SC"))
## testp[varname=="SC",varname:="D"] #merge SC and D
## testp <- testp[step==1, .(value = sum(value)), by = .(particle, age, varname, patch)] # aggr/adult
## testp <- testp[, .(value = mean(value)), by = .(age, varname, patch)] # aggr/adult
## testp[varname=='N',sum(value),by=patch]


## idxp <- grep("N\\[1", BLASTtbmod::get_cols)               #pop 1
## nmzp <- grep("N\\[1", BLASTtbmod::get_cols, value = TRUE) # check

## sum(test[idxp,1,1])
## pops[1]


## ## names
## ## unique(gsub("[[:punct:]]+|[[:digit:]]+", "", BLASTtbmod::get_cols))
## ## str(test) # 1262, particles, times if convert = FALSE
## ## dim(test)
## ## length(BLASTtbmod::get_cols) #need to update

## ## NOTE
## ## CHECK get_cols made loc see dev which refers to test/clean
## ## get_cols <- colnames(y_stocmmod1)[2:ncol(y_stocmmod1)]

## ## cuminf patch
## idx0 <- grep("cum_inf_ByPatch", BLASTtbmod::get_cols)
## nmz0 <- grep("cum_inf_ByPatch", BLASTtbmod::get_cols, value = TRUE) # check
## ## foitemp matrix
## idx1 <- grep("cum_note_flux", BLASTtbmod::get_cols)
## nmz1 <- grep("cum_note_flux", BLASTtbmod::get_cols, value = TRUE) # check
## ## prev
## idx2 <- grep("PrevByPatch", BLASTtbmod::get_cols)
## nmz2 <- grep("PrevByPatch", BLASTtbmod::get_cols, value = TRUE) # check
## ## check infection flux:
## idx <- grep("cum_inf_flux", BLASTtbmod::get_cols)
## nmz <- grep("cum_inf_flux", BLASTtbmod::get_cols,value=TRUE) #check
## tmp <- t(test[idx, 1, ]) # matrix over time for particle 1

## ## matplot(tmp, type = "l")

## ## ## CHECK
## ## p <- 2
## ## t(test[idx0, p, 12]) #cumulative infections by patch
## ## rowSums(matrix(test[idx, p, 12], 7, 7)) #cumulative infections by patch summed
## ## ## test[idx0, p, 12] / sum(test[idx0, p, 12])
## ## ## pops/sum(pops)

## ## matrix(nmz,7,7)
## ## matrix(test[idx, p, 2], 7, 7)
## ## rowSums(matrix(test[idx, p, 2], 7, 7))
## ## test[idx0, p, 2]
## ## ## sum(test[idx0, p, 2]) * pops / sum(pops)

## ## ## FOI:
## ## matrix(test[idx1, p, 12], 7, 7)
## ## rowSums(matrix(test[idx1, p, 12], 7, 7)) # OK cst
## ## chk <- matrix(test[idx1, p, 12], 7, 7)[1, ]
## ## chk / sum(chk) #from
## ## pops / sum(pops)

## ## plot(chk / sum(chk), pops * prev / sum(pops * prev),pch=paste0(1:7))
## ## abline(a = 0, b = 1)

## ## plot(chk / sum(chk), prev / sum(prev), pch = paste0(1:7))

## ## (mldp <- test[idx2, p, 2]) #prev
## ## plot(1e6 * mldp / prev)
## ## ## dinit
## ## prev

## ## plot(prev, 1e5 * mldp, pch = paste0(1:7))
## ## abline(a = 0, b = 1)

## ## TH = test here
## ## ## for rerunning
## test <- run.model(args,
##   0:12,
##   n.particles = 5,
##   convert = FALSE
## )


## ## empirical infection flux
## tmp <- test[idx, 1:dim(test)[2], dim(test)[3]]   #foi flux matrix
## tmp0 <- test[idx0, 1:dim(test)[2], dim(test)[3]] #foi flux totals (vec)
## keep <- abs(colSums(tmp)) < 1e10 # NOTE ditch buggy runs
## tmp <- rowMeans(test[idx, keep, dim(test)[3]]) # mean over particles at last time
## tmp0 <- rowMeans(test[idx0, keep, dim(test)[3]]) # mean over particles at last time
## ## matrix(nmz, ncol = args$patch_dims, nrow = args$patch_dims) # CHECK
## EIF <- matrix(tmp, ncol = args$patch_dims, nrow = args$patch_dims)
## nEIF <- EIF / max(EIF)
## ## FM / max(FM)
## ## sums
## ## sum(EIF)
## ## rowSums(EIF)
## ## tmp0
## 1e2 * sum(EIF) / sum(pops)
## inZ <- rowSums(EIF)
## ## rowSums(FM) * sum(EIF)
## fromZ <- colSums(EIF)
## ## colSums(FM) * sum(EIF)

## ## random mix, even prev expt:
## pops / sum(pops)
## inZ / sum(inZ) #flux into == tmp0/sum(tmp0)
## fromZ / sum(fromZ)

## ## random-mixing expectation: proportional to pop: OK
## plot(sum(EIF) * pops / sum(pops), inZ, pch = paste0(1:7), main = "in/to", xlab = "model", ylab = "data")
## abline(a = 0, b = 1, col = 2)

## ## from follows pop with random mixing
## plot(sum(EIF) * pops / sum(pops), fromZ, pch = paste0(1:7), main = "from", xlab = "model", ylab = "data")
## abline(a = 0, b = 1, col = 2)


## ## random-mixing expectation: proportional to tot TB: lo for 4+5
## ## NOTE pbm even when prevalence is uniform
## plot(sum(EIF) * pops * prev / sum(pops * prev), fromZ,
##   pch = paste0(1:7),
##   main = "from", xlab = "model", ylab = "data"
## )
## abline(a = 0, b = 1, col = 2)

## ## ## OK BUG: chk missing
## ## plot(chk / sum(chk), fromZ / sum(fromZ))
## ## abline(a = 0, b = 1)


## ## NOTE CHECK DONE:
## ## cst MM, cst prev:     I/F/
## ## random mix, cst prev: I/F/
## ## random mix, var prev: I/FX - works for chk vs from, ie FOI vs from
## ## => issue in FOI computation


## ## test[idx0, p, 12]
## ## test[idx0, p, 2]

## ## where are infections happening in/to?
## plot(colSums(FM) * sum(EIF), inZ, pch = paste0(1:7), main = "in/to", xlab = "model", ylab = "data")
## abline(a = 0, b = 1, col = 2)

## ## where are infections happening from?
## plot(rowSums(FM) * sum(EIF), fromZ, pch = paste0(1:7), main = "from", xlab = "model", ylab = "data")
## abline(a = 0, b = 1, col = 2)

## ## comparison
## CF <- data.table(genomic.flux = c(FM / max(FM)), output.FOI.flux = c(nEIF))
## CF[, variable := gsub("cum_inf_flux", "", nmz)]
## CF[, fromid := gsub("\\[", "", gsub(",.\\]", "", variable))]

## ## plot:
## ggplot(
##   CF,
##   aes(
##     x = genomic.flux, y = output.FOI.flux,
##     label = variable, shape = fromid, col = fromid
##   )
## ) +
##   geom_abline(intercept = 0, slope = 1, col = 2) +
##   ## scale_shape(solid=FALSE)+
##   scale_shape_manual(values = 0:6) +
##   geom_point(size = 2) +
##   geom_text_repel(show.legend = FALSE) +
##   theme_classic() +
##   ggpubr::grids()

## ## ggsave(here("plots/flux_cfFOI_new.png"), w = 7, h = 6)


