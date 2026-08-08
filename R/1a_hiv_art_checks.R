## ============================================================
## Output:
## output/hiv_art_dynamics.png
## output/hiv_zone_prevalence.png
## output/hiv_in_tb.png
## ============================================================

years <- 18 + 1 / 12
source(here::here("R/utils/common_fitting_setup.R"))

## --- demographic outputs (age structure) vs scaled national data ---
## NOTE: not expected to match exactly -- comparison data is scaled
## national demographic change, while Blantyre has much higher HIV
## prevalence than the national average.
gp <- plot_compare_demog(out, start_year = 2015, by_comp = "age")
gp


## --- overall HIV/ART dynamics vs (scaled) national survey data ---
gp <- plot_HIV_dynamic(out, start_year = start_year, by_patch = FALSE)
data(hivp_mwi) # fresh copy -- the next line mutates it by reference
hivp_mwi[, step := (Period - start_year) * 12 + 1]
xdta <- hivp_mwi[Period >= 2015]
xdta[
  variable == "HIVpc",
  c("value", "lo", "hi") := .(
    value * mwi_vs_blantyre, lo * mwi_vs_blantyre, hi * mwi_vs_blantyre
  )
]
gp <- gp + geom_pointrange(data = xdta, aes(ymin = lo, ymax = hi), shape = 1)
ggsave(gp, filename = here("output/hiv_art_dynamics.png"), w = 7, h = 5)

## --- zone-wise HIV prevalence vs Blantyre survey data ---
hivpd <- BLASTtbmod::blantyre$hivpre
hivpd <- data.table(
  zone = paste0("Zone ", 1:7),
  step = hivp_mwi[Period == 2015 & variable == "HIVpc", step],
  variable = "HIVpc", value = hivpd
)
gp <- plot_HIV_dynamic(out, start_year = start_year, show_ART = FALSE)
gp <- gp +
  geom_point(data = hivpd, pch = 1, size = 2, stroke = 2) +
  xlim(c(2015, 2025)) + ylim(c(0, NA))
gp <- gp + theme(axis.text.x = element_text(angle = 55, hjust = 1))
ggsave(gp, filename = here("output/hiv_zone_prevalence.png"), w = 12, h = 7)

## --- proportion HIV+ among TB notifications ---
gp <- plot_HIV_in_TB(out, start_year = start_year) + xlim(c(2015, 2021))
ggsave(gp, filename = here("output/hiv_in_tb.png"), w = 7, h = 5)

cat("Saved HIV/ART diagnostic plots to output/\n")
cat("=== DONE ===\n")
