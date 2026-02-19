#' Build RTMB-ready data for the weight-at-age model
#'
#' Converts the package's ADMB-style data list into a standardized structure
#' used by `fit_wt_rtmb()`.
#'
#' @param datfile Either a path to an ADMB `.dat` file or a list returned by
#'   [read_dat()].
#' @param trim_retro Logical. If `TRUE`, mimic the ADMB `retro` trimming rule.
#'
#' @return A standardized list for RTMB fitting.
#' @export
make_wt_rtmb_data <- function(datfile, trim_retro = TRUE) {
  dat <- if (is.character(datfile)) read_dat(datfile) else datfile
  if (!is.list(dat)) stop("`datfile` must be a file path or list.")

  required <- c(
    "cur_yr", "styr", "endyr", "N_data_sets", "N_yrs_data_sets",
    "stage", "endage", "fishery", "fishery_std", "fshry_yrs"
  )
  missing_required <- setdiff(required, names(dat))
  if (length(missing_required) > 0) {
    stop("Missing required ADMB data fields: ", paste(missing_required, collapse = ", "))
  }

  ndat <- as.integer(dat$N_data_sets[1])
  if (ndat < 1L || ndat > 2L) {
    stop("Current RTMB scaffold supports 1-2 datasets; found ndat=", ndat, ".")
  }

  nyrs_decl <- as.integer(as.vector(dat$N_yrs_data_sets))
  if (length(nyrs_decl) < ndat) {
    stop("`N_yrs_data_sets` has fewer entries than `N_data_sets`.")
  }
  nyrs_decl <- nyrs_decl[seq_len(ndat)]

  yrs_data <- list(as.integer(dat$fshry_yrs))
  wt_obs <- list(as.matrix(dat$fishery))
  sd_obs <- list(as.matrix(dat$fishery_std))

  if (ndat > 1L) {
    needed <- c("survey", "survey_std", "survey_yrs")
    missing_survey <- setdiff(needed, names(dat))
    if (length(missing_survey) > 0) {
      stop("Second dataset requested but missing: ", paste(missing_survey, collapse = ", "))
    }
    yrs_data[[2]] <- as.integer(dat$survey_yrs)
    wt_obs[[2]] <- as.matrix(dat$survey)
    sd_obs[[2]] <- as.matrix(dat$survey_std)
  }

  nyrs_data <- integer(ndat)
  for (h in seq_len(ndat)) {
    nyrs_data[h] <- min(
      nyrs_decl[h],
      length(yrs_data[[h]]),
      nrow(wt_obs[[h]]),
      nrow(sd_obs[[h]])
    )
    if (nyrs_data[h] < 1L) {
      stop("Dataset ", h, " has no usable rows after alignment.")
    }
    keep <- seq_len(nyrs_data[h])
    yrs_data[[h]] <- yrs_data[[h]][keep]
    wt_obs[[h]] <- wt_obs[[h]][keep, , drop = FALSE]
    sd_obs[[h]] <- sd_obs[[h]][keep, , drop = FALSE]
  }

  cur_yr <- as.integer(dat$cur_yr[1])
  styr <- as.integer(dat$styr[1])
  endyr <- as.integer(dat$endyr[1])
  age_st <- as.integer(dat$stage[1])
  age_end <- as.integer(dat$endage[1])
  retro <- as.integer(endyr - cur_yr - 2L)

  if (trim_retro && retro > 0L) {
    for (h in seq_len(ndat)) {
      nuse <- max(0L, nyrs_data[h] - retro)
      if (nuse < 1L) stop("Retro trimming removed all rows for dataset ", h, ".")
      keep <- seq_len(nuse)
      nyrs_data[h] <- nuse
      yrs_data[[h]] <- yrs_data[[h]][keep]
      wt_obs[[h]] <- wt_obs[[h]][keep, , drop = FALSE]
      sd_obs[[h]] <- sd_obs[[h]][keep, , drop = FALSE]
    }
  }

  nages <- age_end - age_st + 1L
  for (h in seq_len(ndat)) {
    if (ncol(wt_obs[[h]]) != nages || ncol(sd_obs[[h]]) != nages) {
      stop("Dataset ", h, " does not match expected age span (nages=", nages, ").")
    }
    if (any(yrs_data[[h]] < styr) || any(yrs_data[[h]] > endyr)) {
      stop("Dataset ", h, " has years outside model range [styr, endyr].")
    }
  }

  list(
    cur_yr = cur_yr,
    retro = retro,
    styr = styr,
    endyr = endyr,
    ndat = ndat,
    nyrs_data = as.integer(nyrs_data),
    yrs_data = yrs_data,
    age_st = age_st,
    age_end = age_end,
    wt_obs = wt_obs,
    sd_obs = sd_obs
  )
}

.wt_build_states <- function(dat, par, bias_correct_year = TRUE) {
  nages <- dat$age_end - dat$age_st + 1L
  nyears <- dat$endyr - dat$styr + 1L
  ages <- seq.int(dat$age_st, dat$age_end)

  alpha <- exp(par$log_alpha)
  K <- exp(par$log_K)
  sigma_coh <- exp(par$log_sd_coh)
  sigma_yr <- exp(par$log_sd_yr)

  inv_logit <- function(x) 1 / (1 + exp(-x))
  L1 <- 5 + (50 - 5) * inv_logit(par$L1_u)
  L2 <- 10 + (110 - 10) * inv_logit(par$L2_u)

  age_idx <- seq.int(0L, nages - 1L)
  denom <- 1 - K^(nages - 1L)
  mnwt <- alpha * (L1 + (L2 - L1) * (1 - K^age_idx) / denom)^3
  wt_inc <- mnwt[2:nages] - mnwt[1:(nages - 1L)]

  wt_pre <- matrix(NA, nrow = nyears, ncol = nages)
  wt_pre[1, ] <- mnwt

  yr_bias <- if (isTRUE(bias_correct_year)) 0.5 * sigma_yr^2 else 0
  for (iy in 2:nyears) {
    wt_pre[iy, 1] <- mnwt[1] * exp(0.5 * sigma_coh^2 + sigma_coh * par$coh_eff[iy])
    wt_pre[iy, 2:nages] <-
      wt_pre[iy - 1, 1:(nages - 1L)] + wt_inc * exp(yr_bias + sigma_yr * par$yr_eff[iy])
  }

  wt_hat <- vector("list", dat$ndat)
  residuals <- vector("list", dat$ndat)
  nll_data <- 0

  for (h in seq_len(dat$ndat)) {
    nh <- dat$nyrs_data[h]
    pred_h <- matrix(NA, nrow = nh, ncol = nages)
    for (i in seq_len(nh)) {
      yidx <- dat$yrs_data[[h]][i] - dat$styr + 1L
      pred_i <- wt_pre[yidx, ]
      if (h > 1L) {
        idx <- ((h - 2L) * nages + 1L):((h - 1L) * nages)
        pred_i <- pred_i * par$d_scale[idx]
      }
      pred_h[i, ] <- pred_i
    }
    wt_hat[[h]] <- pred_h
    residuals[[h]] <- dat$wt_obs[[h]] - pred_h

    sd_h <- dat$sd_obs[[h]]
    obs_h <- dat$wt_obs[[h]]
    mn_h <- matrix(mnwt, nrow = nh, ncol = nages, byrow = TRUE)
    nll_data <- nll_data + sum((obs_h - mn_h)^2 / (2 * sd_h^2))
    nll_data <- nll_data + sum((obs_h - pred_h)^2 / (2 * sd_h^2))
  }

  list(
    ages = ages,
    alpha = alpha,
    K = K,
    sigma_coh = sigma_coh,
    sigma_yr = sigma_yr,
    L1 = L1,
    L2 = L2,
    mnwt = mnwt,
    wt_pre = wt_pre,
    wt_hat = wt_hat,
    residuals = residuals,
    nll_data = nll_data
  )
}

#' Convert RTMB fit output to ADMB-style report fields
#'
#' @param fit A fitted object returned by [fit_wt_rtmb()].
#'
#' @return A named list mirroring key fields from ADMB `.rep` output.
#' @export
wt_rtmb_report <- function(fit) {
  if (!inherits(fit, "wt_rtmb_fit")) stop("`fit` must be an object from `fit_wt_rtmb()`.")
  dat <- fit$dat
  par <- fit$par
  states <- .wt_build_states(dat, par, bias_correct_year = fit$bias_correct_year)

  out <- list(
    lof1 = 0,
    lof2 = 0,
    lof3 = 0,
    endyr = dat$endyr,
    retro = dat$retro,
    cur_yr = dat$cur_yr,
    data = do.call(rbind, dat$wt_obs),
    yr = matrix(seq.int(dat$styr, dat$endyr), ncol = 1),
    wt_pre = states$wt_pre
  )

  for (h in seq_len(dat$ndat)) {
    out[[paste0("residuals_", h)]] <- states$residuals[[h]]
  }

  out$sigma_yr <- states$sigma_yr
  out$yr_eff <- par$yr_eff
  out$sigma_coh <- states$sigma_coh
  out$cohort <- ((dat$styr + 1L):dat$endyr) - dat$age_st
  out$coh_eff <- states$sigma_coh * par$coh_eff
  out$ages <- states$ages
  out$mnwt <- states$mnwt
  out$K <- states$K
  out$L1 <- states$L1
  out$L2 <- states$L2
  out
}

#' Fit the weight-at-age model with RTMB
#'
#' Mirrors the ADMB `wt.tpl` likelihood using RTMB random effects.
#'
#' @param datfile Either a path to an ADMB `.dat` file or a list returned by
#'   [read_dat()].
#' @param trim_retro Logical. If `TRUE`, apply ADMB-style retrospective trimming.
#' @param bias_correct_year Logical. If `TRUE`, use the final-phase ADMB form with
#'   `0.5 * sigma_yr^2` in the yearly increment multiplier.
#' @param start Optional named list overriding default starting values.
#' @param control Optional `nlminb` control list.
#' @param compute_sdrep Logical. If `TRUE`, compute `RTMB::sdreport()`.
#'
#' @return A fitted object with optimization output and ADMB-style report fields.
#' @export
fit_wt_rtmb <- function(
    datfile,
    trim_retro = TRUE,
    bias_correct_year = TRUE,
    start = NULL,
    control = NULL,
    compute_sdrep = FALSE) {
  if (!requireNamespace("RTMB", quietly = TRUE)) {
    stop("Package `RTMB` is required for `fit_wt_rtmb()`.")
  }

  dat <- make_wt_rtmb_data(datfile = datfile, trim_retro = trim_retro)
  nages <- dat$age_end - dat$age_st + 1L
  nyears <- dat$endyr - dat$styr + 1L
  nscale <- max(0L, dat$ndat - 1L)

  par <- list(
    L1_u = qlogis((27 - 5) / (50 - 5)),
    L2_u = qlogis((46 - 10) / (110 - 10)),
    log_alpha = -11,
    log_K = -0.13,
    d_scale = rep(1, nscale * nages),
    log_sd_coh = -1.8,
    log_sd_yr = -0.82,
    coh_eff = rep(0, nyears),
    yr_eff = rep(0, nyears)
  )

  if (!is.null(start)) {
    bad_names <- setdiff(names(start), names(par))
    if (length(bad_names) > 0) {
      stop("Unknown start values: ", paste(bad_names, collapse = ", "))
    }
    for (nm in names(start)) par[[nm]] <- start[[nm]]
  }
  if (nscale > 0L) {
    if (is.matrix(par$d_scale)) {
      par$d_scale <- as.vector(t(par$d_scale))
    }
    if (length(par$d_scale) != nscale * nages) {
      stop("`d_scale` must have length ", nscale * nages, " (or matrix ", nscale, " x ", nages, ").")
    }
  }

  objective <- function(parms) {
    RTMB::getAll(parms, dat)

    nages_obj <- dat$age_end - dat$age_st + 1L
    nyears_obj <- dat$endyr - dat$styr + 1L
    inv_logit <- function(x) 1 / (1 + exp(-x))

    alpha <- exp(log_alpha)
    K <- exp(log_K)
    sigma_coh <- exp(log_sd_coh)
    sigma_yr <- exp(log_sd_yr)
    L1 <- 5 + (50 - 5) * inv_logit(L1_u)
    L2 <- 10 + (110 - 10) * inv_logit(L2_u)

    age_idx <- seq.int(0L, nages_obj - 1L)
    mnwt <- alpha * (L1 + (L2 - L1) * (1 - K^age_idx) / (1 - K^(nages_obj - 1L)))^3
    wt_inc <- mnwt[2:nages_obj] - mnwt[1:(nages_obj - 1L)]

    wt_pre <- vector("list", nyears_obj)
    wt_pre[[1]] <- mnwt
    yr_bias <- if (isTRUE(bias_correct_year)) 0.5 * sigma_yr^2 else 0
    for (iy in 2:nyears_obj) {
      wt_pre[[iy]] <- c(
        mnwt[1] * exp(0.5 * sigma_coh^2 + sigma_coh * coh_eff[iy]),
        wt_pre[[iy - 1]][1:(nages_obj - 1L)] + wt_inc * exp(yr_bias + sigma_yr * yr_eff[iy])
      )
    }

    nll <- 0
    for (h in seq_len(dat$ndat)) {
      for (i in seq_len(dat$nyrs_data[h])) {
        yidx <- dat$yrs_data[[h]][i] - dat$styr + 1L
        pred <- wt_pre[[yidx]]
        if (h > 1L) {
          idx <- ((h - 2L) * nages_obj + 1L):((h - 1L) * nages_obj)
          pred <- pred * d_scale[idx]
        }
        obs <- dat$wt_obs[[h]][i, ]
        sdo <- dat$sd_obs[[h]][i, ]
        nll <- nll + sum((obs - mnwt)^2 / (2 * sdo^2))
        nll <- nll + sum((obs - pred)^2 / (2 * sdo^2))
      }
    }

    nll <- nll + 0.5 * sum(coh_eff^2) + 0.5 * sum(yr_eff^2)

    RTMB::REPORT(mnwt)
    RTMB::REPORT(sigma_coh)
    RTMB::REPORT(sigma_yr)
    RTMB::REPORT(K)
    RTMB::REPORT(L1)
    RTMB::REPORT(L2)
    RTMB::ADREPORT(sigma_coh)
    RTMB::ADREPORT(sigma_yr)
    RTMB::ADREPORT(K)
    nll
  }

  obj <- RTMB::MakeADFun(
    objective,
    par = par,
    random = c("coh_eff", "yr_eff"),
    silent = TRUE
  )

  opt_control <- modifyList(list(iter.max = 500, eval.max = 1000), control %||% list())
  opt <- withCallingHandlers(
    stats::nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, control = opt_control),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("number of items to replace is not a multiple", msg, fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )

  par_full <- obj$env$parList()
  par_fixed <- par_full
  par_fixed$coh_eff <- NULL
  par_fixed$yr_eff <- NULL

  fit <- list(
    dat = dat,
    obj = obj,
    opt = opt,
    par_fixed = par_fixed,
    par = par_full,
    bias_correct_year = bias_correct_year,
    convergence = opt$convergence,
    objective = opt$objective
  )
  class(fit) <- "wt_rtmb_fit"

  fit$report <- wt_rtmb_report(fit)
  fit$sdrep <- NULL
  if (isTRUE(compute_sdrep)) {
    fit$sdrep <- tryCatch(RTMB::sdreport(obj), error = function(e) e)
  }

  fit
}

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}
