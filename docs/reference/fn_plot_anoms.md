# Plot weight-at-age anomalies

Plot weight-at-age anomalies

## Usage

``` r
fn_plot_anoms(dfin, maxage = 10, firstyr = 1982, minage = 3, txtsize = 3)
```

## Arguments

- dfin:

  Data frame output from
  [`fn_get_pred()`](https://afsc-assessments.github.io/WtAgeRe/reference/fn_get_pred.md).

- maxage:

  Maximum age to plot.

- firstyr:

  First year to include.

- minage:

  Minimum age to plot.

- txtsize:

  Text label size.

## Value

A ggplot object.
