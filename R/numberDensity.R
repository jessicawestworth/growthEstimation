#' von Mises Probability Density Function
#'
#' A lightweight implementation used to model circular seasonality
#' (e.g., spawning day within a year) without adding dependencies.
#' @param x Angle in radians.
#' @param mu Mean direction in radians.
#' @param kappa Concentration parameter (higher means more concentrated).
#' @return The probability density evaluated at `x`.
#' @export
#' @examples
#' # Density at angle pi when centered at pi with moderate concentration
#' von_mises_pdf(pi, mu = pi, kappa = 2)
von_mises_pdf <- function(x, mu, kappa) {
    denominator <- 2 * pi * besselI(kappa, 0)
    numerator <- exp(kappa * cos(x - mu))
    return(numerator / denominator)
}

#' Spawning Density Function S(d)
#'
#' Calculates relative spawning intensity for numeric dates using a
#' von Mises density on the unit circle of the year.
#' @param numeric_dates A vector of numeric dates (e.g., 2023.45). Only the
#'   fractional part is used to determine day-of-year.
#' @param mu The mean spawning date as a fraction of a year in \[0, 1). For
#'   example, 0.5 is mid-year.
#' @param kappa The concentration parameter for the spawning distribution.
#' @return A numeric vector of relative spawning intensities with the same
#'   length as `numeric_dates`.
#' @export
#' @examples
#' # Peak at mid-year with moderate spread
#' spawning_density(numeric_dates = c(2020.45, 2020.50, 2020.55),
#'                  mu = 0.5, kappa = 4)
spawning_density <- function(numeric_dates, mu, kappa) {
    day_fraction <- numeric_dates %% 1
    day_rad <- day_fraction * 2 * pi
    mu_rad <- mu * 2 * pi
    density <- von_mises_pdf(day_rad, mu = mu_rad, kappa = kappa)
    return(density)
}

#' Get number density as a function of length and time
#'
#' Convolves the Green's function obtained with `getGreens()` with the
#' spawning density from `spawning_density()` to produce the number density
#' as a function of length and time.
#' @inheritParams getGreens
#' @return A matrix holding the number density u(t,l). Rows correspond to time and columns to length.
getNumberDensity <- function(pars, l_max, Delta_l = 1,
                             t_max = 10, Delta_t = 0.05) {
    G <- getGreens(pars, l_max = l_max, Delta_l = Delta_l,
                   t_max = t_max, Delta_t = Delta_t)
    # G has rows = time (age) 0..N_t and cols = size classes

    # Build age/time grid
    N_t <- nrow(G) - 1L
    a_grid <- (0:N_t) * Delta_t

    # Output matrix u with same shape as G (rows time, cols size)
    u <- matrix(0, nrow = N_t + 1L, ncol = ncol(G))

    # Spawning parameters
    mu <- pars$spawning_mu
    kappa <- pars$spawning_kappa

    # For each time t_n, convolve S(t_n - a) with G(a, l) over a in [0, t_n]
    for (n in 0:N_t) {
        # ages considered up to current time
        idx <- 0:n
        # birth dates as numeric years (fractional part encodes day-of-year)
        birth_dates <- a_grid[n + 1L] - a_grid[idx + 1L]
        w <- spawning_density(birth_dates, mu = mu, kappa = kappa)
        # Weighted sum over ages (rows) to get u(t_n, l)
        # Using crossprod to get 1 x N_l row
        u[n + 1L, ] <- as.numeric(crossprod(w, G[idx + 1L, , drop = FALSE])) * Delta_t
    }

    return(u)
}
