#' Get predicted weight-at-age
#'
#' Read ADMB output or RTMB fit objects and return predicted weight-at-age by year.
#'
#' @param file Character vector of base filenames (no extension). Default is "wt".
#'   If `fit` is `NULL`, this can also be a `wt_rtmb_fit` object or a list of
#'   `wt_rtmb_fit` objects.
#' @param fit Optional `wt_rtmb_fit` object or list of `wt_rtmb_fit` objects.
#' @param dat Optional data list used to infer age range (indices 8 and 9).
#' @param source Character vector of model labels for each file.
#' @param age_range Optional numeric vector of ages to use as column names.
#'
#' @return A data frame with year, age columns, and a source column.
#' @export
fn_get_pred <- function(file = c("wt"), fit = NULL, dat = NULL, source = c("model"), age_range = NULL) {
  if (!is.null(fit) && missing(file)) {
    file <- character(0)
  }

  if (is.null(fit)) {
    if (inherits(file, "wt_rtmb_fit")) {
      fit <- list(file)
      file <- character(0)
    } else if (
      is.list(file) &&
      length(file) > 0 &&
      all(vapply(file, inherits, logical(1), what = "wt_rtmb_fit"))
    ) {
      fit <- file
      file <- character(0)
    }
  }

  if (!is.null(fit) && !(
    inherits(fit, "wt_rtmb_fit") ||
      (is.list(fit) && length(fit) > 0 && all(vapply(fit, inherits, logical(1), what = "wt_rtmb_fit")))
  )) {
    stop("`fit` must be a `wt_rtmb_fit` object or list of `wt_rtmb_fit` objects.")
  }

  if (inherits(fit, "wt_rtmb_fit")) fit <- list(fit)

  n_models <- length(file) + if (is.null(fit)) 0 else length(fit)
  if (n_models < 1) stop("Provide at least one ADMB file base name or RTMB fit object.")

  if (length(source) == 1 && n_models > 1) {
    source <- rep(source, n_models)
  }
  if (length(source) != n_models) {
    stop("Length of `source` (", length(source), ") must be 1 or match number of models (", n_models, ").")
  }

  build_pred_df <- function(year, wt_pre, src) {
    wt_pre <- as.matrix(wt_pre)
    colnames(wt_pre) <- paste0("age_", seq_len(ncol(wt_pre)))
    data.frame(year = as.integer(as.vector(year)), wt_pre, source = src, check.names = FALSE)
  }

  df_tmp <- data.frame()
  imodel <- 0

  for (i in seq_along(file)) {
    imodel <- imodel + 1
    dd <- read_rep(paste0(file[i], ".rep"))
    df_tmp <- rbind(df_tmp, build_pred_df(dd$yr, dd$wt_pre, source[imodel]))
  }

  if (!is.null(fit)) {
    for (i in seq_along(fit)) {
      imodel <- imodel + 1
      dd <- if (!is.null(fit[[i]]$report)) fit[[i]]$report else wt_rtmb_report(fit[[i]])
      df_tmp <- rbind(df_tmp, build_pred_df(dd$yr, dd$wt_pre, source[imodel]))
    }
  }

  if (is.null(age_range)) {
    if (!is.null(dat) && !is.null(dat$stage) && !is.null(dat$endage)) {
      age_range <- seq(as.integer(dat$stage[1]), as.integer(dat$endage[1]))
    } else if (!is.null(dat) && !is.null(dat$age_st) && !is.null(dat$age_end)) {
      age_range <- seq(as.integer(dat$age_st[1]), as.integer(dat$age_end[1]))
    } else if (!is.null(dat) && length(dat) >= 9) {
      age_range <- seq(as.integer(unlist(dat[8])), as.integer(unlist(dat[9])))
    } else if (!is.null(df_tmp)) {
      age_range <- seq_len(ncol(df_tmp) - 2)
    }
  }

  if (length(age_range) != (ncol(df_tmp) - 2)) {
    stop("`age_range` length does not match number of age columns in prediction output.")
  }

  names(df_tmp) <- c("year", as.character(age_range), "source")
  df_tmp
}

#' Plot weight-at-age anomalies
#'
#' @param dfin Data frame output from `fn_get_pred()`.
#' @param maxage Maximum age to plot.
#' @param firstyr First year to include.
#' @param minage Minimum age to plot.
#' @param txtsize Text label size.
#'
#' @return A ggplot object.
#' @export
fn_plot_anoms <- function(dfin, maxage = 10, firstyr = 1982, minage = 3, txtsize = 3) {
  p1 <- dplyr::mutate(
    tidyr::pivot_longer(dfin, cols = 2:(ncol(dfin) - 1), names_to = "age", values_to = "wt"),
    age = as.numeric(age)
  ) |>
    dplyr::group_by(age, source) |>
    dplyr::mutate(mnwt = mean(wt, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::filter(year >= firstyr, age >= minage, age <= maxage) |>
    dplyr::mutate(anom = wt / mnwt - 1, Anomaly = ifelse(abs(anom) > 0.5, NA, anom)) |>
    ggplot2::ggplot(ggplot2::aes(y = year, x = age, fill = Anomaly, label = round(wt, 2))) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"), na.value = "white") +
    ggplot2::geom_text(size = txtsize) +
    ggplot2::ylab("Year") +
    ggplot2::xlab("Age") +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_minimal(base_size = 18)

  if (length(unique(dfin$source)) > 1) {
    p1 <- p1 + ggplot2::facet_grid(. ~ source)
  }
  p1
}

#' Replace zeros with near-two for SD
#'
#' @param column Numeric vector.
#' @return Numeric vector with zeros replaced by 1.9999.
#' @export
replace_zeros_with_two <- function(column) {
  column[column == 0] <- 1.9999
  column
}

#' Replace outliers with mean for mean weight-at-age
#'
#' @param column Numeric vector.
#' @return Numeric vector with outliers replaced by the mean of non-zero values.
#' @export
replace_zeros_with_mean <- function(column) {
  mean_without_zeros <- mean(column[column != 0], na.rm = TRUE)
  column[column >= 1.8 * mean_without_zeros | column <= 0.5 * mean_without_zeros] <- mean_without_zeros
  column
}
