#' Get predicted weight-at-age
#'
#' Read ADMB output and return predicted weight-at-age by year.
#'
#' @param file Character vector of base filenames (no extension). Default is "wt".
#' @param dat Optional data list used to infer age range (indices 8 and 9).
#' @param source Character vector of model labels for each file.
#' @param age_range Optional numeric vector of ages to use as column names.
#'
#' @return A data frame with year, age columns, and a source column.
#' @export
fn_get_pred <- function(file = c("wt"), dat = NULL, source = c("model"), age_range = NULL) {
  df_tmp <- data.frame()
  for (i in seq_along(file)) {
    dd <- read_rep(paste0(file[i], ".rep"))
    df_tmp <- rbind(df_tmp, data.frame(year = dd$yr, dd$wt_pre, source = source[i]))
  }

  if (is.null(age_range)) {
    if (!is.null(dat) && length(dat) >= 9) {
      age_range <- seq(as.integer(unlist(dat[8])), as.integer(unlist(dat[9])))
    } else if (!is.null(df_tmp)) {
      age_range <- seq_len(ncol(df_tmp) - 2)
    }
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
