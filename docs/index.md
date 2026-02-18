# Fish body mass at age estimation

A random effects model for estimating variance in weight-at-age.

## Package usage

``` r

# install.packages("devtools")
# devtools::load_all("WtAgeRe")

library(WtAgeRe)

# Example: read predictions and plot anomalies
setwd("examples/ebspollock")
pred <- fn_get_pred(file = "wt", dat = NULL, source = "model")
fn_plot_anoms(pred)
```

## ADMB model

Examples are in `examples/`. Each example directory contains a
`Makefile` and `wt.tpl` for the ADMB model. Run `make` inside the
example directory to build and run the model.
