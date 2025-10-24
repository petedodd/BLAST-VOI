library(here)
library(data.table)
library(ggplot2)
library(ggpubr)
library(epiR)
library(patchwork)

## load data
load(file = here("tmpdata/PB.Rdata"))


## ============ individual figures

## ------------ fig 1a

## data:
BRO <- PB[, .(real.Fprev.mean, real.Fprev.CV, real.foi.mean, real.foi.CV, alph,
  numpatches,
  DB = bnft_opt - bnft_sub
)]



tab <- as.data.table(epi.prcc(BRO))
tabr <- tab[rev(order(abs(est)))][1:5]
tabr$var <- factor(tabr$var, levels = rev(tabr$var), ordered = TRUE)
tabr[, var2 := fcase(
  var == "numpatches", "Number of patches",
  var == "real.Fprev.mean", "TB prevalence mean",
  var == "real.Fprev.CV", "TB prevalence CV",
  var == "real.foi.CV", "FOI CV",
  var == "real.foi.mean", "FOI mean"
)]

tabr$var2 <- factor(tabr$var2, levels = rev(tabr$var2), ordered = TRUE)


## plot:
fig1a <- ggplot(
  tabr,
  aes(var2,
    y = est, ymin = ifelse(est > 0, est, lower),
    ymax = ifelse(est > 0, upper, est)
  )
) +
  geom_bar(stat = "identity") +
  geom_errorbar(width = 0, col = 2, linewidth = 1.5) +
  xlab("") +
  ylab("Partial rank correlation coefficient with VOI") +
  coord_flip() +
  theme_classic() +
  grids()

## fig1a  + theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 2))
## fig1a <- fig1a  + theme(axis.text.y = element_blank()) #maybe?
fig1a

## ------------ fig 1b
BN <- PB[,
  .(
    DB = mean(bnft_opt - bnft_sub) * (1 - exp(-daly$r * daly$LE)) / daly$r,
    DB.sd = sd(bnft_opt - bnft_sub) * (1 - exp(-daly$r * daly$LE)) / daly$r
  ),
  by = numpatches
]


fig1b <- ggplot(BN, aes(numpatches, DB)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  expand_limits(y = c(0, NA)) +
  grids() +
  xlab("Number of zones") +
  ylab("Value of information\n(VOI) in DALYs")

fig1b


## ------------ fig 1c
## smooth
NB <- 10
PBR <- PB[numpatches > 5 & real.Fprev.mean < 0.02]
brks <- seq(from = 0, to = 2000, by = 250)
PBR[, v1c := cut(1e5 * real.Fprev.mean,
  include.lowest = TRUE, ordered_result = TRUE, breaks = brks,
  dig.lab = 4
)]
PBR[, unique(v1c)]

XY <- PBR[,
  .(
    `VOI in DALYs` =
      mean(bnft_opt - bnft_sub) * (1 - exp(-daly$r * daly$LE)) / daly$r
  ),
  by = .(v1c, numpatches)
]

## plot
fig1c <- ggplot(XY, aes(v1c, numpatches, fill = `VOI in DALYs`)) +
  geom_tile() +
  theme_minimal() +
  xlab("TB prevalence per 100,000") +
  ylab("Number of patches") +
  scale_fill_viridis_c(option = "plasma") + # plas,mag
  theme(
    legend.position = "top",
    ## legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

fig1c


## ============ combined figure
((fig1a / fig1b) | fig1c) +
  plot_layout(ncol = 2, widths = c(1, 1.5)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

ggsave(here("output/Figure5.png"), w = 10, h = 6)

