## libraries (also need lubridate)
library(here)
library(ggplot2)
library(ggpubr)
library(data.table)
library(sf)

## === read data
## prevalence & map
prev_dat <- readRDS(here("data/dat_scale.rds"))
zones <- BLASTtbmod::blantyre

## notifications
load(file = here("data/TBN.Rdata"))
load(file = here("data/TB_notes_HIV_patch.Rdata"))
TB_notes_HIV_patch <- as.data.table(TB_notes_HIV_patch)
TB_notes_HIV_patch[, range(monthyr)]
TBN[, monthyr := lubridate::ymd("2015-01-01") + lubridate::dmonths(month)]

## HIV
load(here("data/HIV_inc_1990_2021.Rdata"))

## === phyloflows output
FM <- read.csv(here("data/blantyre_flow_april2026.csv"))
FM <- FM[, 2:ncol(FM)]
FM <- as.matrix(FM)
FMD <- as.data.table(FM)
names(FMD) <- paste0(1:7)
FMD[, fromZone := 1:7]
FMDM <- melt(FMD, id = "fromZone")
names(FMDM)[2] <- "toZone"
FMDM[, c("From zone", "To zone") := .(fromZone, toZone)]


## === panel A: map++
## construct map 10⁻³ 10³ 10⁵
mz <- zones
mz$txt <- paste0(
  "zone: ", zones$zone, "\n",
  "pop: ", round(zones$population / 1e3), "K\n",
  "HIV: ", round(zones$hivpre * 1e2), "%\n",
  "TB: ", round(zones$rate_15), " /10⁵y\n"
)
XY <- matrix(
  unlist(st_geometry(st_centroid(mz$geometry))),
  nrow = 7, ncol = 2, byrow = TRUE
)

## plot
yvc <- 1e-2 * c(1, 1, 1 / 2, 0, 1 / 2, 1, 0)
xvc <- rep(1e-2, 7)
xvc[c(1, 5)] <- xvc[c(1, 5)] + 0.5e-2
xvc[5] <- xvc[5] + 0.5e-2
fig1a <- ggplot() +
  geom_sf(data = mz, aes(fill = as.factor(zone))) +
  scale_fill_brewer(palette = "Set3") +
  geom_text(
    aes(x = XY[, 1] - xvc, y = XY[, 2] - yvc, label = mz$txt),
    hjust = 0
  ) +
  xlab("") +
  ylab("") +
  theme_minimal() +
  theme(legend.position = "none")

fig1a

## === panel B: per capita notifications
TBNR <- TBN[, .(
  pop = sum(population), tbn = 12 * sum(population * notifrate)
),
by = monthyr
]
TBNR[, notes := tbn / pop]

fig1b <- ggplot(TBNR, aes(monthyr, notes)) +
  geom_line() +
  geom_point(shape = 1) +
  theme_classic() +
  grids() +
  xlab("Year") +
  ylab("TB notification\n rate per 10⁵y") +
  expand_limits(y = c(0, NA))

fig1b

## === panel C: HIV/ART
fig1c <- ggplot(
  HIV_inc_1990_2021,
  aes(Year, Adult_incidence_1000_uninfected_est)
) +
  geom_line() +
  geom_point(shape = 1) +
  expand_limits(y = c(0, NA)) +
  theme_classic() +
  grids() +
  xlab("Year") +
  ylab("HIV incidence\n  rate per 10³y")

fig1c

## === panel D: something genomic/fluxy
FMDM[, value_text := paste0(round(100 * value, 1), "%")]

fig1d <- ggplot(
  FMDM,
  aes(`From zone`, `To zone`, fill = value, label = value_text)
) +
  geom_tile() +
  geom_text(col = 2) +
  theme_minimal() +
  theme(legend.position = "none")


fig1d

## ============= COMBINED FIGURE
sidecol <- ggarrange(
  plotlist = list(fig1b, fig1c, fig1d),
  labels = LETTERS[2:4], ncol = 1, nrow = 3
)

mainfig <- ggarrange(
  plotlist = list(fig1a),
  labels = LETTERS[1], ncol = 1, nrow = 1
)


Fig1 <- ggarrange(
  plotlist = list(mainfig, sidecol),
  nrow = 1, ncol = 2,
  widths = c(2, 1)
)

Fig1

## save
ggsave(here("output/Figure1.png"),
  Fig1,
  width = 10, height = 6, dpi = 300
)

## save
ggsave(here("output/Figure1.pdf"),
  Fig1,
  width = 10, height = 6,
  device = cairo_pdf
)
