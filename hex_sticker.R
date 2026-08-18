# hex_sticker.R
# This script creates the heatwave3 hex sticker.
# Adapted from heatwaveR's hex_sticker.R


# Load libraries ----------------------------------------------------------

library(heatwaveR)
library(ggplot2)
library(dplyr)
library(cowplot)
library(hexSticker)
library(showtext)


# The category flame --------------------------------------------------------
# Identical to heatwaveR's own hex_sticker.R (same data, same geom_flame,
# same category colours) so the flame is the same shape as the original
# heatwaveR logo.

ts_dat <- ts2clm(sst_WA, climatologyPeriod = c("1983-01-01", "2012-12-31"))
res <- detect_event(ts_dat)

mhw <- res$clim[10627:10690, ]

clim_cat <- mhw %>%
  dplyr::mutate(diff = thresh - seas,
                thresh_2x = thresh + diff,
                thresh_3x = thresh_2x + diff,
                thresh_4x = thresh_3x + diff)

fillColCat <- c(
  "Moderate" = "#ffc866",
  "Strong"   = "#ff6900",
  "Severe"   = "#9e0000",
  "Extreme"  = "#2d0000"
)

gf3 <- ggplot(data = clim_cat, aes(x = t, y = temp)) +
  geom_flame(aes(y2 = thresh, fill = "Moderate")) +
  geom_flame(aes(y2 = thresh_2x, fill = "Strong")) +
  geom_flame(aes(y2 = thresh_3x, fill = "Severe")) +
  geom_flame(aes(y2 = thresh_4x, fill = "Extreme")) +
  scale_fill_manual(name = NULL, values = fillColCat, guide = "none") +
  theme_void()
gf3


# The pixel grid ------------------------------------------------------------
# A nod to heatwave3's gridded, per-pixel detection: a pyramid of "pixels"
# narrowing towards the most extreme category.


pixel_counts <- c(12, 9, 6, 3)
pixels <- dplyr::bind_rows(lapply(seq_along(fillColCat), function(i) {
  n <- pixel_counts[i]
  data.frame(row = i, x = seq_len(n) - (n + 1) / 2, cat = names(fillColCat)[i])
}))
pixels$cat <- factor(pixels$cat, levels = names(fillColCat))

pixel_grid <- ggplot(data = pixels, aes(x = x, y = -row, fill = cat)) +
  geom_tile(width = 0.8, height = 0.8, colour = NA) +
  scale_fill_manual(name = NULL, values = fillColCat, guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  theme_void() +
  # nudge the grid a little to the right relative to the flame above it
  theme(plot.margin = margin(l = -10, r = 00, unit = "pt"))
pixel_grid

# Flame above, pixel pyramid below, as a single combined image
gf3_grid <- cowplot::plot_grid(gf3, pixel_grid, ncol = 1, rel_heights = c(1, 0.4))
gf3_grid


# Place it on a sticker ---------------------------------------------------

# OperatorMono download location
# https://ifonts.xyz/operator-font-family.html
# Then put the file in the ~/.fonts folder and unzip
# font_add(family = "OperatorMono-MediumItalic", regular = "~/.fonts/OperatorMono-MediumItalic.otf")

# Create hex sticker
sticker(gf3_grid, package = "heatwave3", p_size = 66, p_y = 0.70,
        s_x = 1.06, s_y = 1.02, s_width = 1.6, s_height = 1.8,
        h_fill = "royalblue4", h_color = "#9e0000", p_family = "OperatorMono-MediumItalic",
        filename = "~/heatwave3/logo.png", dpi = 900)
# devtools::check() creates notes when there are things in the root that it doesn't expect
# This includes any .png files that are not labeled "logo.png"
