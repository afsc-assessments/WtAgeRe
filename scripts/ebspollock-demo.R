# WtAgeRe ebspollock demo
# Run from repo root or source in RStudio

library(WtAgeRe)

setwd(here::here("examples", "ebspollock"))
# system("make")

pred <- fn_get_pred(file = "wt", dat = NULL, source = "model")
fn_plot_anoms(pred)
