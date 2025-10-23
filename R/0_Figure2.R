library(here)
library(ggplot2)
library(ggpubr)
library(data.table)
library(patchwork)


## ---- timeseries bits

## create fake data
tz <- seq(from = 0, to = 5, len = 1e3)
ton <- 3
toff <- 3.5
D1 <- data.table(
  tz = tz,
  value = 150 * exp(-0.05 * tz),
  qty = "notifications",
  scenario = "basecase"
)
D2 <- data.table(
  tz = tz,
  value = 15 * exp(-0.05 * tz),
  qty = "mortality",
  scenario = "basecase"
)
D3 <- data.table(
  tz = tz,
  value = 150 * fcase(tz < ton, exp(-0.05 * tz),
    tz > toff, 0.9 * exp(-(0.05 / 1.5) * tz),
    default = exp(-0.05 * tz * 2) * 1.3
  ),
  qty = "notifications",
  scenario = "ACF"
)
D4 <- data.table(
  tz = tz,
  value = 15 * fcase(tz < ton, exp(-0.05 * tz),
    tz > toff & tz <= (toff + 0.5), 2 * (tz - toff) * 0.9 * exp(-0.05 * tz) +
      2 * (0.5 + toff - tz) * exp(-0.05 * tz) / 1.3,
    tz > (toff + 0.5), 0.9 * exp(-0.05 * tz),
    default = exp(-0.05 * tz) / 1.3
  ),
  qty = "mortality",
  scenario = "ACF"
)
DALL <- rbindlist(list(D1, D2, D3, D4))
DALL$qty <- factor(DALL$qty)
DALL$senario <- factor(DALL$scenario)
DALL$quantity <- DALL$qty
sft <- 2
DALL <- DALL[tz >= sft]
DALL[, tz := tz - sft]
top <- 150

F1 <- ggplot(
  DALL,
  aes(tz, value, col = quantity, lty = scenario)
) +
  geom_line() +
  annotate("rect",
    xmin = ton - sft, xmax = toff - sft,
    ymin = 0, ymax = top,
    alpha = .2
  ) +
  scale_y_continuous(limits = c(0, top)) +
  geom_vline(xintercept = ton - sft, col = "blue", lty = 3) +
  geom_vline(xintercept = toff - sft, col = "blue", lty = 3) +
  theme_classic() +
  grids() +
  theme(
    axis.text = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  guides(col = guide_legend(ncol = 1), lty = guide_legend(ncol = 1)) +
  xlab("Time") +
  ylab("Notifications or deaths")

F1

DALL2 <- dcast(DALL, tz + qty ~ scenario, value.var = "value")
DALL2[, df := basecase - ACF]
DALL2[, cumdf := cumsum(df), by = qty]
DALL2[tz < ton - sft | tz > toff - sft, cumdf := NA_real_]

top2 <- 300
F2 <- ggplot(
  DALL2[qty == "mortality"],
  aes(tz, cumdf, col = qty)
) +
  geom_line() +
  annotate("rect",
    xmin = ton - sft, xmax = toff - sft,
    ymin = 0, ymax = top2,
    alpha = .2
  ) +
  scale_y_continuous(limits = c(0, top2)) +
  geom_vline(xintercept = ton - sft, col = "blue", lty = 3) +
  geom_vline(xintercept = toff - sft, col = "blue", lty = 3) +
  theme_classic() +
  grids() +
  theme(legend.position = "none") +
  theme(axis.text = element_blank()) +
  xlab("Time") +
  ylab("Cumulative deaths averted")
F2


## ------- VOI
## multi zone: 49 aggregating to 7
set.seed(12)
slps49 <- rlnorm(7^2, sdlog = 3)
pops49 <- rlnorm(7^2, sdlog = 0.5)
pops <- slps <- rep(0, 7)
for (i in 1:7) {
  this <- (7 * (i - 1) + 1):(7 * i)
  pops[i] <- sum(pops49[this])
  slps[i] <- weighted.mean(slps49[this], pops49[this])
}
slps

der <- rev(order(slps))
slps[der]
der49 <- rev(order(slps49))
slps49[der49]
nmn <- mean(pops)

DR <- D7 <- list()
yr <- xr <- y <- x <- 0
for (i in 1:length(slps)) {
  ## ordered
  pop <- pops[der[i]]
  ben <- slps[der[i]] * pops[der[i]] / nmn
  D7[[i]] <- data.table(xs = x, xe = x + pop, ys = y, ye = y + ben)
  x <- x + pop
  y <- y + ben
  ## random
  popr <- pops[i]
  benr <- slps[i] * pops[i] / nmn
  DR[[i]] <- data.table(xs = xr, xe = xr + popr, ys = yr, ye = yr + benr)
  xr <- xr + popr
  yr <- yr + benr
}
D7 <- rbindlist(D7)
DR <- rbindlist(DR)

D49 <- list()
y <- x <- 0
for (i in 1:length(slps49)) {
  ## ordered
  pop <- pops49[der49[i]]
  ben <- slps49[der49[i]] * pops49[der49[i]] / nmn
  D49[[i]] <- data.table(xs = x, xe = x + pop, ys = y, ye = y + ben)
  x <- x + pop
  y <- y + ben
}
D49 <- rbindlist(D49)


F4 <- ggplot(
  D7,
  aes(x = xs, xend = xe, y = ys, yend = ye)
) +
  geom_point(shape = 1, col = "blue") +
  geom_segment(col = "blue") +
  geom_segment(data = D49, col = "blue", lty = 2) +
  geom_point(data = D49, shape = 20, col = "blue") +
  theme_classic() +
  grids() +
  theme(axis.text = element_blank()) +
  xlab("Population screened = budget used") +
  ylab("Deaths averted")
F4

## 7 zone using above
X <- 30
F3 <- ggplot(
  D7,
  aes(x = xs, xend = xe, y = ys, yend = ye)
) +
  geom_point(shape = 1, col = "blue") +
  geom_segment(col = "blue") +
  geom_segment(data = DR, col = "blue", lty = 2) +
  geom_point(data = DR, shape = 1, col = "blue") +
  geom_segment(
    x = X, y = 41, xend = X, yend = 127,
    arrow = arrow(length = unit(0.03, "npc"), ends = "both"), col = 2
  ) +
  geom_vline(xintercept = X, col = 2, lty = 2) +
  annotate("text", x = X - 12, y = 0, label = "Available budget", col = 2) +
  annotate("text", x = X - 5, y = 75, label = "VOI", col = 2) +
  theme_classic() +
  grids() +
  theme(axis.text = element_blank()) +
  xlab("Population screened = budget used") +
  ylab("Deaths averted")
F3

## ======== combined figure
FC <- (
  F1 + F2 + F3 + F4 +
    plot_layout(ncol = 2, byrow = FALSE)
) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))
FC

ggsave(here("output/Figure2.png"), w = 8, h = 8)
