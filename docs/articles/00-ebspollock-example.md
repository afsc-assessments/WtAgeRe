# WtAgeRe ebspollock example

## Example: EBS pollock

This vignette shows how to run the weight-at-age model using the
`examples/ebspollock` dataset and then visualize anomalies.

``` r

library(WtAgeRe)

# Run ADMB model in the example directory (requires ADMB toolchain)
setwd(here::here("examples", "ebspollock"))
# system("make")

# Read predictions and plot anomalies (requires wt.rep)
if (file.exists("wt.rep")) {
  pred <- fn_get_pred(file = "wt", dat = NULL, source = "model")
  fn_plot_anoms(pred)
}
```
