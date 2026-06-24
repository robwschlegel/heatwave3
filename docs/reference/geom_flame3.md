# Flame geom for MHW visualisation

A ggplot2 geom that fills the area between two lines (temperature and
threshold) where the first exceeds the second, creating "flame"
polygons. For cold-spells, set `reverse = TRUE` to fill where `y < y2`.

## Usage

``` r
geom_flame3(
  mapping = NULL,
  data = NULL,
  stat = "identity",
  position = "identity",
  ...,
  n = 0,
  n_gap = 0,
  reverse = FALSE,
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE
)
```

## Arguments

- mapping:

  Set of aesthetic mappings. Requires `x`, `y`, and `y2`.

- data:

  The data to display.

- stat:

  The statistical transformation. Default `"identity"`.

- position:

  Position adjustment. Default `"identity"`.

- ...:

  Additional arguments passed to the layer.

- n:

  Minimum number of consecutive exceeding points. Default `0`.

- n_gap:

  Maximum gap to bridge between exceeding segments. Default `0`.

- reverse:

  Logical. If `TRUE`, fill where `y < y2` (for cold-spells). Default
  `FALSE` (fill where `y > y2`).

- na.rm:

  Remove NAs? Default `FALSE`.

- show.legend:

  Show legend? Default `NA`.

- inherit.aes:

  Inherit aesthetics from the plot? Default `TRUE`.

## Value

A ggplot2 layer.

## Examples

``` r
library(ggplot2)
# Simple example with synthetic data
df <- data.frame(
  x = 1:100,
  y = sin(seq(0, 4 * pi, length.out = 100)) + 2,
  y2 = 1.5
)
ggplot(df, aes(x = x, y = y, y2 = y2)) +
  geom_flame3() +
  geom_line()

```
