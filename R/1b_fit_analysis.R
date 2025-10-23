## working towards a publication figure
library(here)
library(ggplot2)
library(ggpubr)
library(data.table)
library(BLASTtbmod)

## load data in this case
load(ress0, file = here("tmpdata/ress0.Rdata"))
load(ress1, file = here("tmpdata/ress1.Rdata"))


## ## plots
## plot_compare_noterate_agrgt(ress0)
## plot_compare_noterate_agrgt(ress1)


## str(A$trajectories$state)

## unique(gsub("\\[.+\\]","",BLASTtbmod::get_cols))

## test <- extract.pops.multi(A$trajectories$state, 250, out_type = "notes")

## function to get out relevant data
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


## extract data (NOTE memory hungry)
E0 <- formplotdata(ress0, eps = 0.05)
E1 <- formplotdata(ress1, eps = 0.05)
E0[, acf := "No ACF"]
E1[, acf := "ACF"]
EB <- rbind(E0, E1)


fn <- here("tmpdata/EB.Rdata")
save(EB, file = fn)

load(fn)

## notifications
real_dat <- BLASTtbmod::md7
real_dat[["patch"]] <- paste("Patch", md7$comid)
real_dat$qty <- "noterate"
real_dat$acf <- "No ACF"


## veritcal lines data
VL <- data.table(
  patch = rep(EB[, unique(patch)], each = 2),
  t = EB[, max(t)] - c(0, rep(1:6, each = 2), 7) * 12
)
VL[, item := rep(2:1, 7)]
VLW <- dcast(VL, patch ~ item, value.var = "t")
names(VLW)[2:3] <- c("bot", "top")


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

## plot
ggplot(EB, aes(2015 + t / 12,
  y = mid, ymin = lo, ymax = hi, col = acf, fill = acf,
  group = paste(patch, qty, acf)
)) +
  geom_ribbon(alpha = 0.3, col = NA) +
  geom_line() +
  geom_vline(
    data = VL,
    aes(xintercept = 2015 + t / 12),
    lty = 2, col = "darkgrey"
  ) +
  facet_grid(
    factor(nqty,
      levels = c(
        "Notifications per 100,000 per month",
        "Deaths per 100,000 per month",
        "Cumulative deaths",
        "Difference in cumulative deaths"
      )
    ) ~ zone,
    scales = "free_y", switch = "y"
  ) +
  xlab("Time") +
  ylab("Value") +
  geom_point(data = real_dat, col = 2, shape = 1) +
  theme_linedraw() +
  theme(legend.position = "top", legend.title = element_blank())

ggsave(file = here("output/fit_fig.png"), w = 12, h = 10)


## simplified version
EBR <- EB[qty == "noterate"]


## plot
plt <- "Accent"
ggplot(EBR, aes(2015 + t / 12,
  y = mid, ymin = lo, ymax = hi, col = acf, fill = acf,
  group = paste(patch, qty, acf)
)) +
  geom_ribbon(alpha = 0.3, col = NA) +
  geom_line(lwd = 1) +
  scale_fill_brewer(palette = plt) +
  scale_color_brewer(palette = plt) +
  geom_vline(
    data = VL,
    aes(xintercept = 2015 + t / 12),
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
  guides(fill = guide_legend(nrow = 1), col = guide_legend(nrow = 1))

ggsave(file = here("output/Figure3.png"), w = 8, h = 7)


## ============ looking at figure 4 =========
source(here("R/utils/benefit.R"))
load(here("data/pops.Rdata"))

## EB ddf are the slopes
EB[qty == "ddf" & !is.na(mid)]


## calculate slopes
slpd <- EB[qty == "ddf" & !is.na(mid),
  {
    list(slp = coef(lm(data = .SD, mid ~ t))[2])
  },
  by = patch
] # slopes
true_slopes <- slpd$slp / 1e4

## tests
halfcov <- DbenefitFromSlps(true_slopes,
  pops, # use pops as ranking
  pops, 0.5 * sum(pops),
  verbose = TRUE,
  separate = TRUE
)
halfcov
## loop over coverage
covz <- seq(from = 5e-2, to = 95e-2, by = 5e-2)
iterz <- 1:10e3
B <- list()
k <- 1
for (i in 1:length(covz)) {
  for (j in 1:length(iterz)) {
    B[[k]] <- DbenefitFromSlps(true_slopes,
      sample(7), ## pops, #use pops as ranking
      pops, covz[i] * sum(pops),
      verbose = FALSE,
      separate = TRUE
    )
    B[[k]][["coverage"]] <- covz[i]
    B[[k]][["iter"]] <- iterz[j]
    k <- k + 1
  }
}

B <- rbindlist(B)
B[, db := bnft_opt - bnft_sub]
eps <- 0.25
BS <- B[, .(
  db = mean(db),
  db.lo = mean(db) - sd(db),
  db.hi = mean(db) + sd(db)
  ## db.lo=quantile(db,eps),
  ## db.hi=quantile(db,1-eps)
), by = coverage]
BS <- rbind(
  BS,
  data.table(
    coverage = c(0.0, 1.0), db = rep(0.0, 2), db.lo = rep(0.0, 2),
    db.hi = rep(0.0, 2)
  )
)


r <- 3e-2
LE <- 35
BS[, DALY := db * (1 - exp(-r * LE)) / r]
BS[, DALY.lo := db.lo * (1 - exp(-r * LE)) / r]
BS[, DALY.hi := db.hi * (1 - exp(-r * LE)) / r]


fig4a <- ggplot(
  BS,
  aes(
    x = coverage, y = DALY,
    ymin = DALY.lo, ymax = DALY.hi
  )
) +
  geom_ribbon(alpha = 0.3, col = NA) +
  geom_line() +
  geom_point() +
  theme_classic() +
  ggpubr::grids() +
  scale_x_continuous(label = scales::percent) +
  ylab("VOI in DALYs averted")

fig4a

## USD version
BT <- data.table(CET = seq(from = 0, to = 1e3, by = 1))
BT[, c("voi.mid", "voi.lo", "voi.hi") := BS[
  coverage == 0.5,
  .(
    CET * DALY, CET * DALY.lo,
    CET * DALY.hi
  )
]]

fig4b <- ggplot(BT, aes(x = CET, y = voi.mid, ymin = voi.lo, ymax = voi.hi)) +
  geom_ribbon(alpha = 0.3, col = NA) +
  geom_line() +
  theme_classic() +
  ggpubr::grids() +
  xlab("Cost-effectiveness threshold (USD)") +
  scale_y_continuous(label = scales::comma) +
  ylab("VOI in USD")

fig4b

fig4 <- ggarrange(fig4a, fig4b, nrow = 1, ncol = 2, labels = LETTERS[1:2])
fig4

ggsave(fig4, file = here("output/Figure4.png"), w = 10, h = 5)

