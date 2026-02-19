#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
example_dir <- if (length(args) >= 1) args[1] else file.path("examples", "ebspollock")
base <- if (length(args) >= 2) args[2] else "wt"

datfile <- file.path(example_dir, paste0(base, ".dat"))
repfile <- file.path(example_dir, paste0(base, ".rep"))

if (!file.exists(datfile)) stop("Missing dat file: ", datfile)
if (!file.exists(repfile)) stop("Missing rep file: ", repfile)

source(file.path("R", "readadmb.R"))
source(file.path("R", "rtmb.R"))

dat_raw <- read_dat(datfile)
admb <- read_rep(repfile)
fit <- fit_wt_rtmb(datfile, compute_sdrep = FALSE, control = list(iter.max = 500, eval.max = 1000))
rtmb <- fit$report

ad_year <- as.integer(as.vector(admb$yr))
rt_year <- as.integer(as.vector(rtmb$yr))
years <- intersect(ad_year, rt_year)

if (length(years) == 0) stop("No overlapping years between ADMB and RTMB outputs.")

ad_i <- match(years, ad_year)
rt_i <- match(years, rt_year)

ad_wt <- as.matrix(admb$wt_pre[ad_i, , drop = FALSE])
rt_wt <- as.matrix(rtmb$wt_pre[rt_i, , drop = FALSE])
wt_diff <- rt_wt - ad_wt

mnwt_diff <- as.numeric(rtmb$mnwt) - as.numeric(admb$mnwt)
sigma_yr_diff <- as.numeric(rtmb$sigma_yr) - as.numeric(admb$sigma_yr)
sigma_coh_diff <- as.numeric(rtmb$sigma_coh) - as.numeric(admb$sigma_coh)

res_metrics <- data.frame(
  dataset = integer(0),
  max_abs = numeric(0),
  rmse = numeric(0)
)

for (h in seq_len(fit$dat$ndat)) {
  nm <- paste0("residuals_", h)
  if (!is.null(admb[[nm]]) && !is.null(rtmb[[nm]])) {
    ad_r <- as.matrix(admb[[nm]])
    rt_r <- as.matrix(rtmb[[nm]])
    nr <- min(nrow(ad_r), nrow(rt_r))
    nc <- min(ncol(ad_r), ncol(rt_r))
    d <- rt_r[seq_len(nr), seq_len(nc), drop = FALSE] - ad_r[seq_len(nr), seq_len(nc), drop = FALSE]
    res_metrics <- rbind(
      res_metrics,
      data.frame(
        dataset = h,
        max_abs = max(abs(d), na.rm = TRUE),
        rmse = sqrt(mean(d^2, na.rm = TRUE))
      )
    )
  }
}

cat("RTMB parity check\n")
cat("  example_dir: ", example_dir, "\n", sep = "")
cat("  base:        ", base, "\n", sep = "")
cat(
  "  dat cur/end: ",
  as.integer(dat_raw$cur_yr[1]), "/", as.integer(dat_raw$endyr[1]),
  "\n",
  sep = ""
)
cat(
  "  rep cur/end: ",
  as.integer(admb$cur_yr[1]), "/", as.integer(admb$endyr[1]),
  "\n",
  sep = ""
)
cat("  convergence: ", fit$convergence, "\n", sep = "")
cat("  objective:   ", signif(fit$objective, 7), "\n", sep = "")
cat("  years:       ", min(years), "-", max(years), " (n=", length(years), ")\n", sep = "")
cat("  wt_pre max|diff|: ", signif(max(abs(wt_diff), na.rm = TRUE), 6), "\n", sep = "")
cat("  wt_pre RMSE:      ", signif(sqrt(mean(wt_diff^2, na.rm = TRUE)), 6), "\n", sep = "")
cat("  mnwt max|diff|:   ", signif(max(abs(mnwt_diff), na.rm = TRUE), 6), "\n", sep = "")
cat("  sigma_yr diff:    ", signif(sigma_yr_diff, 6), "\n", sep = "")
cat("  sigma_coh diff:   ", signif(sigma_coh_diff, 6), "\n", sep = "")

if (nrow(res_metrics) > 0) {
  cat("  residual metrics:\n")
  for (i in seq_len(nrow(res_metrics))) {
    cat(
      "    dataset ", res_metrics$dataset[i],
      ": max|diff|=", signif(res_metrics$max_abs[i], 6),
      ", rmse=", signif(res_metrics$rmse[i], 6),
      "\n",
      sep = ""
    )
  }
}

if (
  as.integer(dat_raw$cur_yr[1]) != as.integer(admb$cur_yr[1]) ||
    as.integer(dat_raw$endyr[1]) != as.integer(admb$endyr[1])
) {
  cat("  NOTE: dat/rep year metadata differ; this usually indicates an unmatched pair.\n")
}
