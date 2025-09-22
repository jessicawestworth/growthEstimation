#' Fit growth parameters by minimizing the negative log likelihood with TMB
#'
#' Optimizes k, L_inf, d, m, and annuli_min_age using nlminb().
#' Spawning parameters and annuli_date are treated as data.
#'
#' @param surveys Data frame with columns survey_date (numeric), Length, K, count.
#' @param Delta_l Numeric size step (cm), default 1.
#' @param Delta_t Numeric time step (years), default 0.05.
#' @param spawning_mu Spawning mean date in [0,1).
#' @param spawning_kappa Spawning concentration parameter.
#' @param annuli_date Annual ring formation date in [0,1).
#' @param start Named numeric vector for parameters c(k, L_inf, d, m, annuli_min_age).
#' @param lower Named numeric vector of lower bounds (same names).
#' @param upper Named numeric vector of upper bounds (same names).
#' @param compile Logical, compile the TMB template if needed.
#' @return A list with optimizer result, sdreport, obj, and data used.
#' @export
fit_tmb_nll <- function(
    surveys,
    Delta_l = 1,
    Delta_t = 0.05,
    spawning_mu,
    spawning_kappa,
    annuli_date,
    start,
    lower = c(k = 1e-6, L_inf = 1e-3, d = 1e-6, m = 1e-6, annuli_min_age = 0.0),
    upper = c(k = Inf, L_inf = Inf, d = Inf, m = Inf, annuli_min_age = 5.0),
    compile = TRUE
) {
    stopifnot(all(c("survey_date", "Length", "K", "count") %in% names(surveys)))
    surveys <- as.data.frame(surveys)

    # Build grids similarly to getLogLik()
    l_max <- ceiling(max(surveys$Length) * 1.1)
    t_max <- max(surveys$K) + 2
    N_l <- ceiling(l_max / Delta_l)
    N_t <- ceiling(t_max / Delta_t)
    l_grid <- (1:N_l - 0.5) * Delta_l
    a_grid <- (0:N_t) * Delta_t

    # Map observed Length to nearest grid cell index
    length_index <- pmax(1L, pmin(N_l, as.integer(floor(surveys$Length / Delta_l + 0.5))))

    # Unique surveys and indexing
    survey_levels <- sort(unique(surveys$survey_date))
    survey_index <- match(surveys$survey_date, survey_levels)

    # Observations as vectors, aggregated
    obs_df <- data.frame(
        s = survey_index,
        j = length_index,
        K = as.integer(surveys$K),
        count = as.numeric(surveys$count)
    )
    obs_df <- stats::aggregate(count ~ s + j + K, data = obs_df, sum)
    obs_survey_index <- as.integer(obs_df$s)
    obs_length_index <- as.integer(obs_df$j)
    obs_K <- as.integer(obs_df$K)
    obs_count <- as.numeric(obs_df$count)

    # Prepare TMB data and parameters
    tmb_data <- list(
        obs_survey_index = obs_survey_index,
        obs_length_index = obs_length_index,
        obs_K = obs_K,
        obs_count = obs_count,
        survey_dates = as.numeric(survey_levels),
        l_grid = as.numeric(l_grid),
        a_grid = as.numeric(a_grid),
        spawning_mu = as.numeric(spawning_mu),
        spawning_kappa = as.numeric(spawning_kappa),
        annuli_date = as.numeric(annuli_date),
        Delta_l = as.numeric(Delta_l),
        Delta_t = as.numeric(Delta_t),
        log_eps = log(1e-9)
    )

    tmb_parameters <- list(
        k = as.numeric(start[["k"]]),
        L_inf = as.numeric(start[["L_inf"]]),
        d = as.numeric(start[["d"]]),
        m = as.numeric(start[["m"]]),
        annuli_min_age = as.numeric(start[["annuli_min_age"]])
    )

    obj <- TMB::MakeADFun(
        data = tmb_data,
        parameters = tmb_parameters,
        DLL = "growthEstimation",
        silent = TRUE
    )

    opt <- nlminb(
        start = obj$par,
        objective = obj$fn,
        gradient = obj$gr,
        lower = as.numeric(lower[names(tmb_parameters)]),
        upper = as.numeric(upper[names(tmb_parameters)])
    )

    sdr <- try(TMB::sdreport(obj), silent = TRUE)

    list(
        opt = opt,
        sdreport = sdr,
        obj = obj,
        data = tmb_data
    )
}



