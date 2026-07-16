<div id="main" class="col-md-9" role="main">

# Lollipop geom for event metric visualisation

<div class="ref-description section level2">

A ggplot2 geom that draws lollipop (stem + point) plots for event
metrics.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
geom_lolli3(
  mapping = NULL,
  data = NULL,
  stat = "identity",
  position = "identity",
  ...,
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   mapping:

    Aesthetic mappings. Requires `x` and `y`.

-   data:

    The data to display.

-   stat:

    Statistical transformation. Default `"identity"`.

-   position:

    Position adjustment. Default `"identity"`.

-   ...:

    Additional arguments.

-   na.rm:

    Remove NAs? Default `FALSE`.

-   show.legend:

    Show legend? Default `NA`.

-   inherit.aes:

    Inherit aesthetics? Default `TRUE`.

</div>

<div class="section level2">

## Value

A ggplot2 layer.

</div>

<div class="section level2">

## Details

Each lollipop is drawn independently from its own row – unlike
`geom_flame3`, nothing here depends on row order or grouping within a
series. For depth-resolved data, map `depth` to `colour` (or `shape`) to
distinguish levels, or use `ggplot2::facet_wrap(~ depth)`; both work
with the standard ggplot2 aesthetics and need no special handling from
this geom.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
library(ggplot2)
df <- data.frame(x = as.Date("2020-01-01") + seq(0, 300, by = 30),
                 y = c(1.2, 2.1, 0.8, 3.5, 1.9, 2.7, 0.5, 1.1, 2.3, 1.6, 0.9))
ggplot(df, aes(x = x, y = y)) + geom_lolli3()

```

</div>

</div>

</div>
