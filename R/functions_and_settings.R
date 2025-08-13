# functions  ----------------------------------------------

#' Apply transformation to runoff data
#'
#' Supports inverse, square root, power, and Box-Cox transformations of runoff data.
#' Log transformation is included but generally not recommended for KGE calculations
#' (see airGR or Santos, 2018). Use Use \code{\link{KGEtang}} instead.
#'
#' @param Q A numeric vector, matrix, or data frame of runoff values.
#' @param q_transfo_type A string specifying the transformation type. Options are:
#'   "none", "sqrt", "inv", "log", "boxcox", "boxcoxsantos", and "power".
#' @param lambda A numeric value used for Box-Cox and power transformations. Default is 0.25.
#' @param ... Additional arguments passed to mean(), e.g., na.rm = TRUE, used in boxcoxsantos.
#' @note Consider how to handle infinite values resulting from transformations.
#' @return Transformed runoff data in the same format as input.
#' @export
#' @examples
#' transfo_q(array(0:10, c(2, 5)), "sqrt")
transfo_q <- function(Q, q_transfo_type = "none", lambda = 0.25, ...) {
  if (is.null(lambda)) lambda <- 0.25
  if (is.na(lambda)) lambda <- 0.25

  switch(q_transfo_type,
    none = Q_transfo <- Q,
    sqrt = Q_transfo <- sqrt(Q),
    inv = Q_transfo <- 1 / Q,
    log = Q_transfo <- log(Q),
    boxcox =
      if (lambda == 0) {
        Q_transfo <- log(Q)
      } else {
        Q_transfo <- (Q^lambda - 1) / lambda
      },
    boxcoxsantos =
      if (lambda == 0) {
        Q_transfo <- log(Q)
      } else {
        Q_transfo <- (Q^lambda - (.01 * mean(Q, ...))^lambda) / lambda
      },
    power =
      Q_transfo <- Q^lambda,
    stop("Invalid transformation choice: \"",
      q_transfo_type, "\". Choose one of the following: none, sqrt, inv, log, boxcox, power.",
      call. = FALSE
    )
  )

  return(Q_transfo)
}


#' Kling-Gupta Efficiency (KGE'') after Tang et al. (2021)
#'
#' Computes the modified Kling-Gupta Efficiency (KGE) as proposed by Tang et al. (2021),
#' https://doi.org/10.1175/jcli-d-21-0067.1, Equation (5)
#'
#' @param Qsim A numeric vector of simulated runoff values.
#' @param Qobs A numeric vector of observed runoff values.
#' @param single_output A boolean if a single combined or
#' all components should be output
KGEtang <- function(Qsim, Qobs, single_output = TRUE, ...) {
  or_NA <- is.na(Qsim) | is.na(Qobs)
  Qsim_valid <- Qsim[!or_NA]
  Qobs_valid <- Qobs[!or_NA]

  meanSim <- mean(Qsim_valid)
  meanObs <- mean(Qobs_valid)
  varSim <- var(Qsim_valid)
  varObs <- var(Qobs_valid)
  rProd <- cor(Qsim_valid, Qobs_valid)

  xBeta <- meanSim / meanObs
  yBeta <- (meanObs - meanSim) / sqrt(varObs)
  alpha <- sqrt(varSim) / sqrt(varObs)

  KGEtang <- 1 - sqrt(yBeta^2 + (alpha - 1)^2 +
    (rProd - 1)^2)

  if (single_output) {
    out <- KGEtang
  } else {
    out <- list()
    out$NSE <- 2 * alpha * rProd - yBeta^2 - alpha^2

    out$KGE <- 1 - sqrt((xBeta - 1)^2 + (alpha - 1)^2 +
      (rProd - 1)^2)

    out$KGEtang <- KGEtang
    out$meanSim <- meanSim
    out$meanObs <- meanObs
    out$varSim <- varSim
    out$varObs <- varObs
    out$rProd <- rProd
    out$xBeta <- xBeta
    out$yBeta <- yBeta
    out$alpha <- alpha
  }
  return(out)
}


#' Calculate Performance Criterion
#'
#' Computes a performance criterion for comparing simulated and observed runoff.
#' If the criterion is `"KGEtang"`, it uses the \code{\link{KGEtang}} function.
#' Otherwise, it delegates to \code{calc_hydroGOF()} for other supported criteria.
#'
#' @param error_crit A string specifying the error criterion to compute.
#'   Supported values include `"KGEtang"` and any criterion accepted by \code{calc_hydroGOF()}.
#' @param Qsim A numeric vector of simulated runoff values.
#' @param Qobs A numeric vector of observed runoff values.
#'
#' @return A numeric value representing the selected performance criterion.
#'
#' @seealso \code{\link{KGEtang}}, \code{calc_hydroGOF}
#'
#' @export
#' @examples
#' Qsim <- c(1, 2, 3, 4, 5)
#' Qobs <- c(1.1, 2.1, 2.9, 4.2, 5.1)
#' calc_crit("KGEtang", Qsim, Qobs)

calc_crit <- function(error_crit, Qsim, Qobs) {
  if (error_crit == "KGEtang") {
    crit <- KGEtang(Qsim, Qobs)
  } else {
    crit <- calc_hydroGOF(error_crit, Qsim, Qobs)
  }
  return(crit)
}


#' Wrapper function around \code{\link{hydroGOF}} functions
#'
#' Calculates Goodness-of-Fit functions for two runoff series
#'
#' @param GOF_fun a function, or (todo consider only a functional) a string with
#'   function name, of the format \code{GOF_fun(Qsim, Qobs, na.rm = "TRUE")},
#'   typically from the \code{\link{hydroGOF}} package
#' @param Qsim vector, matrix, data.frame etc of simulated runoff values
#' @param Qobs vector, matrix, data.frame etc of observed runoff values
#' @param na.rm a logical value indicating if NA should be removed
#'
#' @seealso \code{\link{hydroGOF}}
#' @note This function requires the `hydroGOF` package to be installed. Not
#'   imported as this package depends on deprecated packages as sp etc.
#'
#' @return transformed runoff in same format as input
#' @export
#' @examples
#' # hydroGOF must be loaded
#' library(hydroGOF)
#' calc_hydroGOF(KGE, 1:10, seq(0, 9))
#' # this is NA
#' calc_hydroGOF("KGE", 1:10, rep(0, 10))
#' # this is also NA
#' calc_hydroGOF(KGE, 1:10, as.numeric(rep(NA, 10)))
calc_hydroGOF <- function(GOF_fun, Qsim, Qobs, na.rm = TRUE) {
  tryCatch(
    error = function(cnd) {
      stop("Goodness of fit function not applicable. \"", ifelse(is.function(GOF_fun), "", GOF_fun),
        "\". Choose for example one of the hydroGOF package functions.",
        call. = FALSE
      )
    },
    {
      GOF_fun <- ifelse(is.function(GOF_fun), GOF_fun, get(GOF_fun))
    }
  )
  # infinite values to NA
  Qobs[is.infinite(Qobs)] <- NA
  Qsim[is.infinite(Qsim)] <- NA

  out <- GOF_fun(Qsim, Qobs, na.rm = na.rm)
  return(out)
}

#' Get the Family or Package of a Hydrological Model
#'
#' Performs a reverse lookup to identify the package or model family associated
#' with a given hydrological model name.
#'
#' @param model_str A string specifying the name of a hydrological model (e.g., "GR4J").
#' @return A string indicating the model's package or family name. Returns `"none"` if the model is not found.
#' @export
#' @examples
#' get_family("GR4J")
#' get_family("UnknownModel")
get_family <- function(model_str) {
  family <- hydro_family$family[hydro_family$model == model_str]
  if (identical(family, character(0))) family <- "none"
  return(family)
}


#' Load Meteorological Data
#'
#' Loads a data frame stored in an `.rds` file containing meteorological data
#' formatted similar to the airGR model family. The input must include specific
#' columns as described below.
#'
#' @details
#' The input data frame must contain the following columns:
#' \itemize{
#'   \item \code{DatesR}: Dates in \code{Date} or \code{POSIXt} format
#'   \item \code{P}: Average precipitation \[mm/day\]
#'   \item \code{T}: Catchment average air temperature \[°C\]
#'   \item \code{E}: Catchment average potential evapotranspiration \[mm/day\]
#'   \item \code{Qmm}: Outlet discharge \[mm/day\]
#' }
#'
#' For more information on the required format, see the airGR documentation:
#' \url{https://hydrogr.github.io/airGR/page_1_get_started.html}
#'
#' @param file A string specifying the path to the `.rds` file containing the meteorological data.
#' @param tzone A string specifying the time zone to assign to the \code{DatesR} column.
#'   If \code{"UTC"} (default), a \code{Date} vector is not converted to another time zone.
#'
#' @return A data frame named \code{BasinObs} containing the loaded and time-adjusted data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' load_meteo_data("D:/input/airGR/HSU_2044.rds")
#' }
load_meteo_data <- function(file, tzone = "UTC") {
  BasinObs <- readr::read_rds(file)

  # convert Date to POSIXct
  BasinObs$DatesR <- as.POSIXct(BasinObs$DatesR)
  attr(BasinObs$DatesR, "tzone") <- tzone

  return(BasinObs)
}


#' Create Model Input Structure
#'
#' Creates an input structure tailored to the specified hydrological model.
#' The function supports models from the \pkg{airGR}, \pkg{TUWmodel}, \pkg{hydromad},
#' and \pkg{topmodel} packages. The structure and content of the input depend on the
#' model's requirements.
#'
#' @param model A string specifying the hydrological model (e.g., \code{"GR4J"}, \code{"TUW"}, \code{"topmodel"}).
#'   For a complete list see table in `vignette("model_overview")`.
#' @param BasinObs A data frame containing time series of meteorological and hydrological data,
#'   typically created using \code{\link{load_meteo_data}}.
#' @param BasinInfo A list containing spatial and catchment-specific information such as
#'   \code{HypsoData}, \code{delay}, or \code{topidx}, depending on the model.
#'
#' @return A model-specific input object, either a data frame, list, or S3 class, depending on the model type.
#'
#' @details
#' \itemize{
#'   \item For \code{airGR} models, the function returns an object created by \code{airGR::CreateInputsModel()}.
#'   \item For \code{TUW}, \code{TUWsnow}, \code{snowsnow}, and \code{hydromad} models, a renamed data frame is returned.
#'   \item For \code{topmodel}, a list is returned including additional elements like \code{delay} and \code{topidx}.
#' }
#'
#' @seealso \code{\link{load_meteo_data}}, \code{\link{get_family}}, \code{airGR::CreateInputsModel()}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' BasinObs <- load_meteo_data("D:/input/airGR/HSU_2044.rds")
#' BasinInfo <- list(HypsoData = c(500, 600, 700), delay = 2, topidx = runif(10))
#' input <- create_input("GR4J", BasinObs, BasinInfo)
#' }
create_input <- function(model, BasinObs, BasinInfo) {
  if (model %in% c("TUW", "TUWsnow", "snowsnow", "topmodel") | get_family(model) == "hydromad") {
    lookup <- c(Q = "Qmm", date = "DatesR")

    input <- dplyr::select(BasinObs, dplyr::any_of(c("DatesR", "P", "E", "T", "Qmm"))) %>%
      dplyr::rename(dplyr::any_of(lookup))

    # convert to a list for topmodel to also include delay and topidx
    if (model == "topmodel") {
      input <- input %>% as.list()
      input$delay <- BasinInfo$delay
      input$topidx <- BasinInfo$topidx
    }
  } else if (get_family(model) == "airGR") {
    input <- airGR::CreateInputsModel(
      FUN_MOD = paste0("RunModel_", model), DatesR = BasinObs$DatesR,
      Precip = BasinObs$P, PotEvap = BasinObs$E,
      TempMean = BasinObs$T,
      ZInputs = median(BasinInfo$HypsoData),
      HypsoData = BasinInfo$HypsoData, NLayers = 1
    )
  } else {
    stop("Invalid model choice: \"", model,
      "\". Currently implemented are those from the airGR, TUWmodel, hydromad and topmodel packages.",
      call. = FALSE
    )
  }
  return(input)
}


#' Split a Date Vector into Warm-up, Calibration, and Validation Periods
#'
#' Splits a date vector or a data frame with a date column into three time periods:
#' warm-up, calibration, and validation. The function allows for optional adjustment
#' of the calibration and validation periods based on the presence of missing data
#' at the beginning of the time series.
#'
#' @param df A vector of dates (e.g., `Date`, `POSIXt`) or a data frame containing a
#'   `DatesR` column and a `Qmm` column (used to detect the first non-NA value).
#' @param start_end_date_vec A character vector of length six, specifying the start
#'   and end dates for the warm-up, calibration, and validation periods, in that order.
#' @param ensure_warm_up Logical. If `TRUE`, adjusts the warm-up period to start at
#'   the first non-NA value in `Qmm`, if applicable. Default is `TRUE`.
#' @param adjust_cal_end Logical. If `TRUE`, the end date of the calibration period is adjusted
#'   proportionally to the shift in the warm-up period, preserving the original calibration-to-validation
#'   duration ratio. This ensures that the calibration period remains representative even if the warm-up
#'   period is shifted due to missing data.
#' @param adjust_val_start Logical. If `TRUE`, the start date of the validation period is adjusted
#'   to immediately follow the (potentially shifted) calibration period. This ensures continuity between
#'   calibration and validation periods when the calibration end date has been modified.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{ind_warm}{Indices corresponding to the warm-up period.}
#'   \item{ind_cal}{Indices corresponding to the calibration period.}
#'   \item{ind_val}{Indices corresponding to the validation period.}
#' }
#'
#' @examples
#' \dontrun{
#' dates <- seq(as.Date("2000-01-01"), as.Date("2010-12-31"), by = "month")
#' df <- data.frame(DatesR = dates, Qmm = c(rep(NA, 12), runif(length(dates) - 12)))
#' periods <- split_data_set(df, c("2000-01-01", "2002-12-31", "2003-01-01", "2006-12-31", "2007-01-01", "2010-12-31"))
#' }
#'
#' @export
split_data_set <- function(df, start_end_date_vec, ensure_warm_up = TRUE,
                           adjust_cal_end = FALSE,
                           adjust_val_start = FALSE) {
  # Convert the dates to Date objects
  date_vec <- as.Date(start_end_date_vec)

  # determine relation between dates
  diffs <- diff(date_vec)[c(1, 3, 5)] %>% as.numeric()

  # Calculate ratios
  ratios <- diffs / sum(diffs)

  # backwards compatibility, do nothing with adjusting dates otherwise
  if (is.data.frame(df)) {
    # Get the date of the first non-NA value
    first_non_na_date <- df$DatesR[which(!is.na(df$Qmm))[1]]

    # only shift things if first date is NA
    if (first_non_na_date != df$DatesR[1] & ensure_warm_up) {
      # determine warm-up length in years
      warm_up_length <- date_vec[2] - date_vec[1]

      # Get the start year, increment if the month is not January
      start_year <- as.integer(format(first_non_na_date, "%Y"))
      if (format(first_non_na_date, "%m") != "01") {
        start_year <- start_year + 1
      }

      # re-define warm up period
      warm_up_start_date_orig <- date_vec[1]
      date_vec[1] <- as.Date(paste0(start_year, "-01-01"))
      shift_dates <- date_vec[1] - warm_up_start_date_orig
      date_vec[2] <- date_vec[1] + warm_up_length
      # Shift the start of calibration period back, starting one day later as warmup
      date_vec[3] <- date_vec[2] + 1

      # Shift the end of calibration period back, ~ the relation calibration to total in input
      if (adjust_cal_end) {
        date_vec[4] <- date_vec[4] + shift_dates * ratios[3]

        # Get the start year, increment if the month is not January
        start_year_val <- as.integer(format(date_vec[4], "%Y"))
        if (format(date_vec[4], "%m") != "01") {
          start_year_val <- start_year_val + 1
        }
        date_vec[4] <- as.Date(paste0(start_year_val, "-01-01")) - 1
      }

      # Shift the start of validation period back, starting one day later as calibration
      if (adjust_val_start) {
        # start of validation, one day later
        date_vec[5] <- date_vec[4] + 1
      }
    }
    # backwards compatibility
  } else {
    # df will be actually a date vector
    df$DatesR <- df
  }

  ind_warm <- seq(
    which(df$DatesR %>% as.Date() == date_vec[1]),
    which(df$DatesR %>% as.Date() == date_vec[2])
  )

  ind_cal <- seq(
    which(df$DatesR %>% as.Date() == date_vec[3]),
    which(df$DatesR %>% as.Date() == date_vec[4])
  )

  ind_val <- seq(
    which(df$DatesR %>% as.Date() == date_vec[5]),
    which(df$DatesR %>% as.Date() == date_vec[6])
  )

  return(list(ind_warm = ind_warm, ind_cal = ind_cal, ind_val = ind_val))
}


#' Find Monthly Indices in a Date Vector
#'
#' Returns the indices of dates that fall within specified months. Useful for
#' subsetting time series data by month.
#'
#' @param date A vector of dates (e.g., `Date`, `POSIXt`).
#' @param months A character vector of two-digit month strings (e.g., \code{c("02", "03")}).
#' @param ind An optional vector of indices used to subset the input date vector before filtering.
#'   Defaults to all indices in the date vector.
#'
#' @return An integer vector of indices corresponding to
find_monthly_indices <- function(date, months, ind = seq_along(date)) {
  dates_df <- data.frame(date = date[ind]) %>%
    dplyr::mutate(ind = ind) %>%
    dplyr::mutate(month = format(as.Date(date), "%m")) %>%
    dplyr::filter(month %in% months)

  return(dates_df$ind)
}


#' Min-Max Normalization and Re-Transformation
#'
#' Scales numeric data to the \[0, 1\] range using min-max normalization, or
#' re-transforms normalized data back to its original scale. The behavior depends
#' on the specified direction.
#'
#' @param x A numeric vector, matrix, or array to be scaled or re-transformed.
#' @param min The minimum value used for scaling. Required for both directions,
#'   but typically set to \code{min(x)} when \code{direction = "RT"}.
#' @param max The maximum value used for scaling. Required for both directions,
#'   but typically set to \code{max(x)} when \code{direction = "RT"}.
#' @param direction A character string indicating the direction of transformation:
#'   \code{"RT"} (real to transformed) for normalization, or \code{"TR"} (transformed to real)
#'   for re-scaling. Default is \code{"RT"}.
#'
#' @return A numeric object of the same shape as \code{x}, either normalized or
norm_minmax <- function(x, min, max, direction = "RT") {
  if (direction == "RT") {
    xt <- (x - min) / (max - min)
  } else if (direction == "TR") {
    xt <- x * (max - min) + min
  } else {
    stop("Invalid direction choice: \"", direction,
      "\". Choose one of the following: RT or TR.",
      call. = FALSE
    )
  }
  return(xt)
}


#' Parameter Transformation to and from Hypercube
#'
#' Transforms model parameters to the unit hypercube \[0, 1\] and back, depending on the
#' specified direction. For `airGR` models, the corresponding transformation functions
#' from the \pkg{airGR} package are used. For other models, a min-max normalization
#' approach is applied using \code{\link{norm_minmax}}.
#'
#' Model combinations such as \code{"CemaNeigeGR4J"} are supported. Parameters for snow
#' modules can also be included and transformed independently.
#'
#' @param param A numeric vector of model parameters.
#' @param direction A character string indicating the direction of transformation:
#'   \code{"RT"} (real to transformed) or \code{"TR"} (transformed to real).
#' @param model A character string specifying the hydrological model (e.g., \code{"GR4J"},
#'   \code{"TUW"}, \code{"CemaNeigeGR4J"}).
#'   For a complete list see table in `vignette("model_overview")`.
#' @param snow_module Optional. A character string specifying the snow module to be included
#'   in the transformation (currently \code{"CemaNeige"} and \code{"TUWsnow"}).
#' @param add_snow_par Logical. If \code{TRUE}, snow module parameters are included in the
#'   transformation. Default is \code{FALSE}.
#' @param cal_parameter A list containing calibration parameter bounds (lower and upper)
#'   for each model and snow module. Defaults to \code{default_cal_par}.
#'
#' @return A numeric vector of transformed or re-transformed parameters.
#'
#' @seealso \code{\link[airGR]{TransfoParam}}, \code{\link{norm_minmax}}
#'
#' @note
#' \enumerate{
#'   \item \code{CemaNeigeHyst} is not yet implemented.
#'   \item The \code{airGR} transformation functions require the \pkg{airGR} package to be installed.
#'   \item Future versions may include an option to bypass \pkg{airGR} transformations entirely.
#' }
#'
#' @export
#'
#' @examples
#' # Scale a parameter set for model "TUW" to [0,1] and back
#' param <- c(1, 2, 3, -1, 1, 1, 200, 10, 1, 15, 100, 50, 2, 15, 50)
#' scaled <- transfo_param(param, "RT", "TUW")
#' rescaled <- transfo_param(scaled, "TR", "TUW")
transfo_param <- function(param, direction, model, snow_module = NULL, add_snow_par = FALSE,
                          cal_parameter = default_cal_par) {
  # airGR models
  if (get_family(model) == "airGR") {
    # find the transformation function based on the model, remove CemaNeige
    if (stringr::str_detect(model, "CemaNeige")) {
      with_CemaNeige <- TRUE
      nof_param_hydro <- readr::parse_number(model)

      # find the right tranformation function based on the hydro model
      transfo_fn_hdyro <- get(paste0(
        "TransfoParam_",
        stringr::str_remove(model, "CemaNeige")
      ))
      param_trans_hdyro <- transfo_fn_hdyro(
        ParamIn = param[1:nof_param_hydro],
        Direction = direction
      )

      # transform the last two parameters assumed to be for CemaNeige
      # todo: include CemaNeigeHyst
      param_trans_snow <- TransfoParam_CemaNeige(
        ParamIn = param[c(nof_param_hydro + 1, nof_param_hydro + 2)],
        Direction = direction
      )

      # combine hydro an snow parameters
      param_trans <- c(param_trans_hdyro, param_trans_snow)

      # without CemaNeige
    } else {
      transfo_fn <- get(paste0("TransfoParam_", model))
      param_trans <- transfo_fn(ParamIn = param, Direction = direction)
    }

    # no airGR model to be scaled to [0,1]
  } else {
    # snow module is also added and transformed to [0,1], independently if there is  a method for airGR transformation
    cal_par <- list()
    if (!is.null(snow_module) & add_snow_par) {
      cal_par$lower <- c(cal_parameter[[snow_module]][["lower"]], cal_parameter[[model]][["lower"]])
      cal_par$upper <- c(cal_parameter[[snow_module]][["upper"]], cal_parameter[[model]][["upper"]])
    } else {
      cal_par$lower <- cal_parameter[[model]][["lower"]]
      cal_par$upper <- cal_parameter[[model]][["upper"]]
    }

    param_trans <- purrr::pmap_dbl(list(
      param, cal_par[["lower"]], cal_par[["upper"]], direction
    ), norm_minmax)
  }

  return(param_trans)
}


#' Call Calibration Function for Hydrological Model
#'
#' Executes a specified optimization algorithm to calibrate a hydrological model.
#' Supports multiple calibration methods and Monte Carlo sampling strategies.
#'
#' @param cal_fn A string specifying the calibration function or method to use.
#'   Supported options include \code{"DEoptim"}, \code{"hydroPSO"}, and \code{"malschains"}.
#'   Monte Carlo variants can be specified using the format \code{"method__sampling__nruns"}.
#'   For a complete list see table in `vignette("calibration_methods_overview")`.
#' @param hydro_data A list or data frame containing observed runoff data, typically
#'   loaded using \code{\link{load_meteo_data}}.
#' @param split_indices A list of index vectors (e.g., from \code{\link{split_data_set}})
#'   indicating warm-up and calibration periods.
#' @param model A string specifying the hydrological model to calibrate.
#'    For a complete list see table in `vignette("model_overview")`.
#' @param input A list of model input data, typically created using \code{\link{create_input}}.
#' @param snow_module Optional. A string specifying the snow module used (currently \code{"CemaNeige"} and \code{"TUWsnow"}).
#' @param snow_input Optional. Input data for the snow module.
#' @param snow_parameters Optional. Initial or fixed parameters for the snow module.
#' @param error_crit A string naming the error criterion function (e.g., \code{"KGE"}).
#'   Must be compatible with \code{\link{calc_hydroGOF}} or from the \pkg{hydroGOF} package.
#' @param cal_maximize Logical. If \code{TRUE}, the calibration maximizes the objective function.
#' @param cal_q_transfo A string indicating how runoff should be transformed
#'   (see \code{\link{transfo_q}}).
#' @param lambda Optional. A numeric value or vector used for regularization or multi-objective weighting.
#' @param do_transfo_param Logical. If \code{TRUE}, parameters are transformed to a unit hypercube
#'   before calibration.
#' @param cal_par A list of calibration settings specific to the chosen calibration function.
#'   Defaults to \code{default_cal_par}, but can be customized by the user.
#'
#' @return A list or object containing the results of the calibration, including optimized
#'   parameters and performance metrics.
#'
#' @note This function requires the `hydroPSO` package to be installed. Not
#'   imported as this package depends on deprecated packages as sp etc.
#'
#' @import airGR
#' @import DEoptim
#' @import Rmalschains
#' @importFrom lhs randomLHS
#' @importFrom randtoolbox sobol
#' @importFrom stats nlminb optim
#' @export
#'
#' @examples
#' \dontrun{
#' cal_output <- call_cal_fn(
#'   cal_fn = "DEoptim",
#'   hydro_data = hydro_data,
#'   split_indices = split_data_set(...),
#'   model = "GR4J",
#'   input = create_input(...),
#'   error_crit = "KGE",
#'   cal_maximize = TRUE,
#'   cal_q_transfo = "none",
#'   do_transfo_param = TRUE,
#'   cal_par = default_cal_par
#' )
#' }
call_cal_fn <- function(cal_fn, hydro_data, split_indices, model, input,
                        snow_module = NULL, snow_input = NULL, snow_parameters = NULL,
                        error_crit, cal_maximize, cal_q_transfo, lambda, do_transfo_param,
                        cal_par = default_cal_par) {
  # functions --------------
  # determine mc_method
  get_mc_parameters <- function(cal_fn) {
    dummy <- stringr::str_split(cal_fn, pattern = "__", simplify = TRUE)

    if (!length(dummy) %in% c(1, 2, 3)) stop("cal_fn monte carlo not correctly specified, needs to be max 2 strings separated by _.")

    cal_fn <- dummy[1]
    nruns <- NULL
    mc_method <- NULL

    if (length(dummy) >= 2) {
      mc_method <- dummy[2]
    }
    if (length(dummy) == 3) {
      nruns <- dummy[3] %>% as.numeric()
    }

    return(list(cal_fn = cal_fn, mc_method = mc_method, nruns = nruns))
  }

  # run nruns with random parameters, returns a df with single column
  create_random_sample <- function(nruns = 100, method = "random") {
    if (is.null(nruns)) {
      nruns <- 100
    }

    if (is.null(method)) {
      method <- "random"
    }

    # latin hypercube
    switch(method,
      lhs = par_samples <- lhs::randomLHS(nruns, k = length(lower)),
      sobol = par_samples <- randtoolbox::sobol(nruns, dim = length(lower), scrambling = 3),
      random = par_samples <- purrr::map_vec(1:nruns, ~ t(runif(n = length(lower)))), # runif works with min and max, but I need to blow it up for the other methods anyway
      stop("Invalid monte carlo sampling choice: \"",
        method, "\". Choose one of the following: random, lhs, sobol.",
        call. = FALSE
      )
    )

    # to from 0, 1 to lower and upper, t(kronecker(matrix(1, 1, nruns), lower)) allows element wise multiplication, as it replicates the vector to the shape of par_samples
    par_samples <- par_samples %>%
      norm_minmax(min = t(kronecker(matrix(1, 1, nruns), lower)), max = t(kronecker(matrix(1, 1, nruns), upper)), direction = "TR")

    return(par_samples)
  }


  run_random_optim <- function(nruns = 100, method = "random") {
    par_samples <- create_random_sample(nruns, method)

    start <- Sys.time()

    par_comb <- par_samples %>%
      tidyr::as_tibble(.name_repair = "unique") %>%
      dplyr::rowwise() %>%
      dplyr::mutate(crit = optim_fn(
        dplyr::c_across(dplyr::everything() %>% c()),
        hydro_data, split_indices, model, input,
        snow_module, snow_input, snow_parameters,
        error_crit, cal_maximize, cal_q_transfo, lambda, do_transfo_param,
        airGR_RunOptions, airGR_RunOptions_snow_module
      ))

    end <- Sys.time()

    # todo: check if arrange must be desc or not depending of cal_maximize
    par_crit <- par_comb %>%
      dplyr::arrange(crit) %>%
      dplyr::first()

    cat(sprintf(
      "Random sampling with method %s finished in %s %s with best value %s...\n",
      cal_fn, round(end - start, 3), attributes(end - start)$units, ifelse(cal_maximize, -signif(par_crit$crit, 2), signif(par_crit$crit, 2))
    ))

    return(par_crit)
  }

  # check if always minimize is appropriate
  run_random_optim_threshold <- function(threshold_crit, maxnruns, method = "random", nr_par_sets = 10) {
    par_samples <- create_random_sample(maxnruns, method)

    start <- Sys.time()

    # always minimize at this point as if maximize crit is negative
    # todo: check this
    par_crit <- NULL
    for (i in seq_along(par_samples[, 1])) {
      crit_run <- optim_fn(
        par_samples[i, ],
        hydro_data, split_indices, model, input,
        snow_module, snow_input, snow_parameters,
        error_crit, cal_maximize, cal_q_transfo, lambda, do_transfo_param,
        airGR_RunOptions, airGR_RunOptions_snow_module
      )

      if (crit_run < threshold_crit) {
        # add row to df, t() is needed as par_sample[i,] is not a matrix but a numeric atomic vector
        par_crit <- par_crit %>% dplyr::bind_rows(par_samples[i, ] %>% t() %>% tidyr::as_tibble(.name_repair = "unique") %>% dplyr::mutate(crit = crit_run))

        # stop for loop if enough parameters are found
        if (nrow(par_crit) == nr_par_sets) {
          break
        }
      }
    }
    end <- Sys.time()

    if (is.null(par_crit)) {
      cat(sprintf(
        "No parameters found with method %s finished after %s runs in %s %s ...\n",
        cal_fn, i, signif(end - start, 3), attributes(end - start)$units
      ))
    } else {
      cat(sprintf(
        "Random sampling with method %s finished after %s runs, found %s parameter sets in %s %s ...\n",
        cal_fn, i, nrow(par_crit), round(end - start, 3), attributes(end - start)$units
      ))
    }
    print(par_crit)

    return(par_crit)
  }

  steepest_descent <- function(ParamStartR, ParamStartT, CritStart, verbose = TRUE) {
    NParam <- length(ParamStartT)
    OptimParam <- rep(TRUE, NParam)

    NRuns <- 0
    CritName <- NULL
    CritBestValue <- NULL

    HistParamR <- matrix(NA, nrow = 500 * NParam, ncol = NParam)
    HistParamT <- matrix(NA, nrow = 500 * NParam, ncol = NParam)
    HistCrit <- matrix(NA, nrow = 500 * NParam, ncol = 1)

    HistParamR[1, ] <- ParamStartR
    HistParamT[1, ] <- ParamStartT
    HistCrit[1, ] <- CritStart
    ProposeCandidatesLoc <- function(NewParamOptimT, OldParamOptimT,
                                     RangesT, OptimParam, Pace) {
      if (nrow(NewParamOptimT) != 1 | nrow(OldParamOptimT) !=
        1) {
        stop("each input set must be a matrix of one single line")
      }
      if (ncol(NewParamOptimT) != ncol(OldParamOptimT) | ncol(NewParamOptimT) !=
        length(OptimParam)) {
        stop("each input set must have the same number of values")
      }
      NParam <- ncol(NewParamOptimT)
      VECT <- NULL
      for (I in 1:NParam) {
        if (OptimParam[I]) {
          for (J in 1:2) {
            Sign <- 2 * J - 3
            Add <- TRUE
            PotentialCandidateT <- NewParamOptimT
            PotentialCandidateT[1, I] <- NewParamOptimT[I] +
              Sign * Pace
            if (PotentialCandidateT[1, I] < RangesT[
              1,
              I
            ]) {
              PotentialCandidateT[1, I] <- RangesT[
                1,
                I
              ]
            }
            if (PotentialCandidateT[1, I] > RangesT[
              2,
              I
            ]) {
              PotentialCandidateT[1, I] <- RangesT[
                2,
                I
              ]
            }
            if (NewParamOptimT[I] == RangesT[1, I] & Sign <
              0) {
              Add <- FALSE
            }
            if (NewParamOptimT[I] == RangesT[2, I] & Sign >
              0) {
              Add <- FALSE
            }
            if (identical(PotentialCandidateT, OldParamOptimT)) {
              Add <- FALSE
            }
            if (Add) {
              VECT <- c(VECT, PotentialCandidateT)
            }
          }
        }
      }
      Output <- NULL
      Output$NewCandidatesT <- matrix(VECT,
        ncol = NParam,
        byrow = TRUE
      )
      return(Output)
    }
    if (verbose) {
      message("Steepest-descent local search in progress")
    }
    Pace <- 0.64
    PaceDiag <- rep(0, NParam)
    CLG <- 0.7^(1 / NParam)
    Compt <- 0
    CritOptim <- CritStart
    RangesT <- rbind(lower, upper)
    RangesR <- RangesT
    NewParamOptimT <- ParamStartT
    OldParamOptimT <- ParamStartT
    for (ITER in 1:(100 * NParam)) {
      if (Pace < 0.01) {
        break
      }
      CandidatesParamT <- ProposeCandidatesLoc(
        NewParamOptimT,
        OldParamOptimT, RangesT, OptimParam, Pace
      )$NewCandidatesT


      CandidatesParamR <- CandidatesParamT
      # CandidatesParamR <- FUN_TRANSFO(CandidatesParamT, "TR")
      # CandidatesParamR <- apply(CandidatesParamR, 1, function(x) {
      #   x[!OptimParam] <- CalibOptions$FixedParam[!OptimParam]
      #   return(x)
      # })
      if (NParam > 1) {
        # CandidatesParamR <- t(CandidatesParamR)
      } else {
        CandidatesParamR <- cbind(CandidatesParamR)
      }
      iNewOptim <- 0
      for (iNew in 1:nrow(CandidatesParamT)) {
        Param <- CandidatesParamR[iNew, ]


        crit <- optim_fn(
          Param,
          hydro_data, split_indices, model, input,
          snow_module, snow_input, snow_parameters,
          error_crit, cal_maximize, cal_q_transfo, lambda, do_transfo_param,
          airGR_RunOptions, airGR_RunOptions_snow_module
        )

        # OutputsModel <- RunModel(InputsModel, RunOptions,
        #   Param,
        #   FUN_MOD = FUN_MOD, ...
        # )
        # OutputsCrit <- ErrorCrit(InputsCrit, OutputsModel,
        #   verbose = FALSE
        # )


        if (!is.na(crit)) {
          # if (!is.na(OutputsCrit$CritValue)) {
          if (crit <
            # if (OutputsCrit$CritValue * OutputsCrit$Multiplier <
            CritOptim) {
            CritOptim <- crit
            # CritOptim <- OutputsCrit$CritValue * OutputsCrit$Multiplier
            iNewOptim <- iNew
          }
        }
      }
      NRuns <- NRuns + nrow(CandidatesParamR)
      if (iNewOptim != 0) {
        OldParamOptimT <- NewParamOptimT
        NewParamOptimT <- matrix(CandidatesParamT[
          iNewOptim,
          1:NParam
        ], nrow = 1)
        Compt <- Compt + 1
        if (Compt > 2 * NParam) {
          Pace <- Pace * 2
          Compt <- 0
        }
        VectPace <- NewParamOptimT - OldParamOptimT
        for (iC in 1:NParam) {
          if (OptimParam[iC]) {
            PaceDiag[iC] <- CLG * PaceDiag[iC] + (1 -
              CLG) * VectPace[iC]
          }
        }
      } else {
        Pace <- Pace / 2
        Compt <- 0
      }
      if (ITER > 4 * NParam) {
        NRuns <- NRuns + 1
        iNewOptim <- 0
        iNew <- 1
        CandidatesParamT <- NewParamOptimT + PaceDiag
        if (!is.matrix(CandidatesParamT)) {
          CandidatesParamT <- matrix(CandidatesParamT,
            nrow = 1
          )
        }
        for (iC in 1:NParam) {
          if (OptimParam[iC]) {
            if (CandidatesParamT[iNew, iC] < RangesT[
              1,
              iC
            ]) {
              CandidatesParamT[iNew, iC] <- RangesT[
                1,
                iC
              ]
            }
            if (CandidatesParamT[iNew, iC] > RangesT[
              2,
              iC
            ]) {
              CandidatesParamT[iNew, iC] <- RangesT[
                2,
                iC
              ]
            }
          }
        }
        CandidatesParamR <- CandidatesParamT
        # CandidatesParamR <- FUN_TRANSFO(
        #   CandidatesParamT,
        #   "TR"
        # )
        Param <- CandidatesParamR[iNew, ]

        crit <- optim_fn(
          Param,
          hydro_data, split_indices, model, input,
          snow_module, snow_input, snow_parameters,
          error_crit, cal_maximize, cal_q_transfo, lambda, do_transfo_param,
          airGR_RunOptions, airGR_RunOptions_snow_module
        )

        if (!is.na(crit)) {
          if (crit <
            CritOptim) {
            CritOptim <- crit
            iNewOptim <- iNew
          }
        }

        # OutputsModel <- RunModel(InputsModel, RunOptions,
        #   Param,
        #   FUN_MOD = FUN_MOD, ...
        # )
        # OutputsCrit <- ErrorCrit(InputsCrit, OutputsModel,
        #   verbose = FALSE
        # )
        # if (OutputsCrit$CritValue * OutputsCrit$Multiplier <
        #   CritOptim) {
        #   CritOptim <- OutputsCrit$CritValue * OutputsCrit$Multiplier
        #   iNewOptim <- iNew
        # }


        if (iNewOptim != 0) {
          OldParamOptimT <- NewParamOptimT
          NewParamOptimT <- matrix(CandidatesParamT[
            iNewOptim,
            1:NParam
          ], nrow = 1)
        }
      }
      # NewParamOptimR <- FUN_TRANSFO(NewParamOptimT, "TR")
      NewParamOptimR <- NewParamOptimT
      HistParamR[ITER + 1, ] <- NewParamOptimR
      HistParamT[ITER + 1, ] <- NewParamOptimT
      HistCrit[ITER + 1, ] <- CritOptim
    }
    ITER <- ITER - 1
    if (CritOptim == CritStart & verbose) {
      message("\t No progress achieved")
    }
    ParamFinalR <- NewParamOptimR
    ParamFinalT <- NewParamOptimT
    CritFinal <- CritOptim
    NIter <- 1 + ITER

    if (!is.null(lambda)) {
      CritName <- paste0(error_crit, "_", cal_q_transfo, "_", lambda)
    } else {
      CritName <- paste0(error_crit, "_", cal_q_transfo)
    }
    if (verbose) {
      message(sprintf(
        "\t Calibration completed (%s iterations, %s runs)",
        NIter, NRuns
      ))
      message("\t     Param = ", paste(sprintf("%8.3f", ParamFinalR),
        collapse = ", "
      ))
      message(sprintf(
        "\t     Crit. %-12s = %.4f", CritName,
        CritFinal
      ))
    }
    HistParamR <- cbind(HistParamR[1:NIter, ])
    colnames(HistParamR) <- paste0("Param", 1:NParam)
    HistParamT <- cbind(HistParamT[1:NIter, ])
    colnames(HistParamT) <- paste0("Param", 1:NParam)
    HistCrit <- cbind(HistCrit[1:NIter, ])
    # BoolCrit_Actual <- InputsCrit$BoolCrit
    # BoolCrit_Actual[OutputsCrit$Ind_notcomputed] <- FALSE
    # MatBoolCrit <- cbind(InputsCrit$BoolCrit, BoolCrit_Actual)
    # colnames(MatBoolCrit) <- c("BoolCrit_Requested", "BoolCrit_Actual")
    OutputsCalib <- list(
      ParamFinalR = as.double(ParamFinalR),
      CritFinal = CritFinal, NIter = NIter, NRuns = NRuns,
      HistParamR = HistParamR, HistCrit = HistCrit,
      # MatBoolCrit = MatBoolCrit,
      CritName = CritName, CritBestValue = CritBestValue
    )
    # class(OutputsCalib) <- c("OutputsCalib", "HBAN")
    return(OutputsCalib)
  }



  # main --------------
  # determine upper and lower boundaries for parameter set
  lower <- cal_par[[model]][["lower"]]
  upper <- cal_par[[model]][["upper"]]

  # add snow parameters
  add_snow_par <- FALSE
  airGR_RunOptions_snow_module <- NULL
  if (!is.null(snow_module) & is.null(snow_parameters)) {
    lower <- c(cal_par[[snow_module]][["lower"]], lower)
    upper <- c(cal_par[[snow_module]][["upper"]], upper)
    add_snow_par <- TRUE

    # create RunOptions for airGR models to speed up runtime
    if (get_family(snow_module) == "airGR") {
      ind_warm_cal <- c(split_indices$ind_warm, split_indices$ind_cal)

      airGR_RunOptions_snow_module <- airGR::CreateRunOptions(
        FUN_MOD = paste0("RunModel_", snow_module),
        InputsModel = snow_input,
        IndPeriod_Run = ind_warm_cal,
        # output only PliqAndMelt during calibration, much faster
        Outputs_Sim = "PliqAndMelt"
      )
    }
  }

  # create RunOptions for airGR models to speed up runtime
  if (get_family(model) == "airGR") {
    ind_warm_cal <- c(split_indices$ind_warm, split_indices$ind_cal)

    airGR_RunOptions <- airGR::CreateRunOptions(
      FUN_MOD = paste0("RunModel_", model),
      InputsModel = input,
      IndPeriod_Run = ind_warm_cal,
      # output only Qsim during calibration, much faster
      Outputs_Sim = "Qsim",
    )
  }

  # par_start <- model_parameters[[main_basin]][[model]]$Calibration_Michel$KGE$log

  # model_parameters[[main_basin]][[model]]$Calibration_Michel$KGE$log %>% transfo_param("RT", model)
  # model_parameters[[main_basin]][[model]]$DEoptim$KGE$log %>% transfo_param("RT", model)

  # CalibOptions <- CreateCalibOptions(FUN_MOD = paste0("RunModel_", model), FUN_CALIB = Calibration_Michel)
  # par_start_dist <- CalibOptions$StartParamDistrib


  # if (do_transfo_param) {
  #   for (i in 1:nrow(par_start_dist)) {
  #     par_start_dist[i,] <- transfo_param(par_start_dist[i,], "RT", model)
  #   }
  # }

  # transform to hypercube if chosen
  if (do_transfo_param) {
    lower <- transfo_param(lower, "RT", model, snow_module, add_snow_par, cal_parameter = cal_par)
    upper <- transfo_param(upper, "RT", model, snow_module, add_snow_par, cal_parameter = cal_par)

    # par_start <- transfo_param(par_start, "RT", model, snow_module, add_snow_par)
  }
  if (cal_fn == "steepest_descent") {
    # monte carlo simulations
    initial_par_crit <- run_random_optim()

    initial_crit <- initial_par_crit %>%
      dplyr::select(crit) %>%
      as.numeric()
    par_start <- initial_par_crit %>%
      dplyr::select(-crit) %>%
      as.numeric()

    cal_results <- steepest_descent(par_start %>% as.matrix() %>% t(), par_start %>% as.matrix() %>% t(), initial_crit)

    # return best parameters
    model_param <- cal_results$ParamFinalR

    # return best error criterion
    error_crit_val <- cal_results$CritFinal
  } else if (stringr::str_starts(cal_fn, "montecarlo")) {
    # get mc_pars
    mc_pars_list <- get_mc_parameters(cal_fn)

    # monte carlo simulations
    par_start <- run_random_optim(mc_pars_list$nruns, method = mc_pars_list$mc_method)

    # # test threshold option
    # par_start <- run_random_optim_threshold(-.85, maxnruns = mc_pars_list$nruns, method = mc_pars_list$mc_method)
    # par_start <- par_start %>% dplyr::arrange(crit) %>% dplyr::first()

    # return best parameters
    model_param <- par_start %>%
      dplyr::select(-crit) %>%
      as.numeric()

    # return best error criterion
    error_crit_val <- par_start %>% dplyr::pull(crit)

    print(error_crit_val)

    cal_results <- mc_pars_list
  } else if (stringr::str_starts(cal_fn, "nlminb")) {
    # get mc_pars
    mc_pars_list <- get_mc_parameters(cal_fn)

    if (is.null(mc_pars_list$mc_method)) {
      # start from mean
      par_start <- (upper - lower) / 2 + lower
      # add criterium of initial parameters and store it as data frame
      initial_par_crit <- par_start %>%
        as.matrix() %>%
        t() %>%
        tidyr::as_tibble(.name_repair = "unique") %>%
        dplyr::mutate(crit = optim_fn(
          dplyr::c_across(dplyr::everything() %>% c()),
          hydro_data, split_indices, model, input,
          snow_module, snow_input, snow_parameters,
          error_crit, cal_maximize, cal_q_transfo, lambda, do_transfo_param,
          airGR_RunOptions, airGR_RunOptions_snow_module
        ))
    } else {
      # determine best parameter set on random samples
      initial_par_crit <- run_random_optim(mc_pars_list$nruns, method = mc_pars_list$mc_method)
      par_start <- initial_par_crit %>%
        dplyr::select(-crit) %>%
        as.numeric()
    }

    # # mode of all catchments
    # par_start_mode <- readRDS(paste0("D:/gitlabext/hydroens_droughtch/analysis/data/median_mode_parameters_KGE_log_", model, ".rds")) %>% dplyr::pull(mode_par)

    # PORT
    cal_results <- stats::nlminb(
      start = par_start,
      # start = par_start_mode,
      objective = optim_fn,
      hydro_data = hydro_data, split_indices = split_indices,
      model = model, input = input,
      snow_module = snow_module, snow_input = snow_input, snow_parameters = snow_parameters,
      error_crit = error_crit, cal_maximize = cal_maximize,
      cal_q_transfo = cal_q_transfo, lambda = lambda, do_transfo_param = do_transfo_param,
      lower = lower, upper = upper,
      airGR_RunOptions = airGR_RunOptions, airGR_RunOptions_snow_module = airGR_RunOptions_snow_module,
      control = list(trace = 10)
    )

    # return best parameters
    model_param <- cal_results$par

    # return best error criterion
    error_crit_val <- cal_results$objective

    cal_results$initial_par_crit <- initial_par_crit
  } else if (stringr::str_starts(cal_fn, "optim")) {
    # get mc_pars
    mc_pars_list <- get_mc_parameters(cal_fn)

    if (is.null(mc_pars_list$mc_method)) {
      # start from mean
      par_start <- (upper - lower) / 2 + lower
      # add criterium of initial parameters and store it as data frame
      initial_par_crit <- par_start %>%
        as.matrix() %>%
        t() %>%
        tidyr::as_tibble(.name_repair = "unique") %>%
        dplyr::mutate(crit = optim_fn(
          dplyr::c_across(dplyr::everything() %>% c()),
          hydro_data, split_indices, model, input,
          snow_module, snow_input, snow_parameters,
          error_crit, cal_maximize, cal_q_transfo, lambda, do_transfo_param,
          airGR_RunOptions, airGR_RunOptions_snow_module
        ))
    } else {
      # determine best parameter set on random samples
      initial_par_crit <- run_random_optim(mc_pars_list$nruns, method = mc_pars_list$mc_method)
      par_start <- initial_par_crit %>%
        dplyr::select(-crit) %>%
        as.numeric()
    }

    # # mode of all catchments
    # par_start_mode <- readRDS(paste0("D:/gitlabext/hydroens_droughtch/analysis/data/median_mode_parameters_KGE_log_", model, ".rds")) %>% dplyr::pull(mode_par)

    # PORT
    cal_results <- stats::optim(
      par = par_start,
      # par = par_start_mode,
      fn = optim_fn,
      hydro_data = hydro_data, split_indices = split_indices,
      model = model, input = input,
      snow_module = snow_module, snow_input = snow_input, snow_parameters = snow_parameters,
      error_crit = error_crit, cal_maximize = cal_maximize,
      cal_q_transfo = cal_q_transfo, lambda = lambda, do_transfo_param = do_transfo_param,
      lower = lower, upper = upper,
      airGR_RunOptions = airGR_RunOptions, airGR_RunOptions_snow_module = airGR_RunOptions_snow_module,
      method = "L-BFGS-B",
      control = list(trace = 1)
    )

    # return best parameters
    model_param <- cal_results$par

    # return best error criterion
    error_crit_val <- cal_results$value

    # add par_start_df
    cal_results$initial_par_crit <- initial_par_crit
  } else if (cal_fn == "DEoptim") {
    # differential evolution
    cal_results <- DEoptim::DEoptim(
      fn = optim_fn, hydro_data = hydro_data, split_indices = split_indices,
      model = model, input = input,
      snow_module = snow_module, snow_input = snow_input, snow_parameters = snow_parameters,
      error_crit = error_crit, cal_maximize = cal_maximize,
      cal_q_transfo = cal_q_transfo, lambda = lambda, do_transfo_param = do_transfo_param,
      lower = lower, upper = upper,
      airGR_RunOptions = airGR_RunOptions, airGR_RunOptions_snow_module = airGR_RunOptions_snow_module,
      control = DEoptim::DEoptim.control(
        itermax = cal_par[[model]][["DEoptim"]][["itermax"]],
        trace = 10
      )
    )

    # return best parameters
    model_param <- cal_results$optim$bestmem

    # return best error criterion
    error_crit_val <- cal_results$optim$bestval
  } else if (cal_fn == "hydroPSO") {
    # particle swarm, this one needs hydro_data = hydro_data, DEoptim  not
    # todo: hydroPSO has a normalize option, maybe this can be used here
    #   instead of own transformations?
    cal_results <- hydroPSO::hydroPSO(
      fn = optim_fn, hydro_data = hydro_data, split_indices = split_indices,
      model = model, input = input,
      snow_module = snow_module, snow_input = snow_input, snow_parameters = snow_parameters,
      error_crit = error_crit, cal_maximize = cal_maximize,
      cal_q_transfo = cal_q_transfo, lambda = lambda, do_transfo_param = do_transfo_param,
      airGR_RunOptions = airGR_RunOptions, airGR_RunOptions_snow_module = airGR_RunOptions_snow_module,
      lower = lower, upper = upper,
      control = cal_par[[model]][["hydroPSO"]][["control"]]
    )

    # return best parameters
    model_param <- cal_results$par

    # return best error criterion
    error_crit_val <- cal_results$value
  } else if (cal_fn == "malschains") {
    # one needs to have the optimization function in the environment of this very function
    # (i.e.call_cal_fn) to work, at least this is how I got it to work
    # also, it needs only one argument, this is why this wrapper is here
    # todo: any better solutions to both above mentioned points?
    optim_fn_single <- function(ParamOptim) {
      return(optim_fn(
        ParamOptim, hydro_data, split_indices, model, input,
        snow_module, snow_input,
        snow_parameters = snow_parameters,
        error_crit, cal_maximize, cal_q_transfo, lambda, do_transfo_param,
        airGR_RunOptions, airGR_RunOptions_snow_module
      ))
    }

    # MA-LS chains
    env_single <- environment(fun = optim_fn_single)

    cal_results <- Rmalschains::malschains(
      fn = optim_fn_single,
      env = env_single,
      lower = lower, upper = upper,
      maxEvals = cal_par[[model]][["malschains"]][["maxEvals"]]
    )

    # return best parameters
    model_param <- cal_results$sol

    # return best error criterion
    error_crit_val <- cal_results$fitness
  } else {
    stop("Invalid choice of calibration function: \"", cal_fn,
      "\". Choose one of the following: Calibration_Michel, DEoptim, hydroPSO, malschains, nlminb, optim, montecarlo_XXX.",
      call. = FALSE
    )
  }

  # rescale parameters back to real if needed
  if (do_transfo_param) {
    model_param <- transfo_param(model_param, "TR", model, snow_module, add_snow_par, cal_parameter = cal_par)
  }

  #  create output list
  return(list(
    model_param = model_param,
    error_crit_val = error_crit_val,
    cal_results = cal_results
  ))
}


#' Objective Function for Hydrological Model Calibration
#'
#' This function is used as the objective function during model calibration.
#' It simulates a hydrological model over the warm-up and calibration periods
#' and evaluates an error criterion (e.g., KGE) on the calibration period only.
#'
#' The function handles parameter transformation, snow module simulation,
#' and error handling for cases where the simulated runoff is invalid (e.g., all NAs or zeros).
#'
#' @param ParamOptim A numeric vector of model parameters to be optimized.
#' @param hydro_data A list or data frame containing observed runoff, typically loaded
#'   using \code{\link{load_meteo_data}}.
#' @param split_indices A list of index vectors indicating warm-up and calibration periods,
#'   usually from \code{\link{split_data_set}}.
#' @param model A string specifying the hydrological model.
#'   For a complete list see table in `vignette("model_overview")`.
#' @param input A list of model input data, typically created using \code{\link{create_input}}.
#' @param snow_module Optional. A string specifying the snow module (currently \code{"CemaNeige"} and \code{"TUWsnow"}).
#' @param snow_input Optional. Input data for the snow module.
#' @param snow_parameters Optional. A vector of fixed snow parameters. If \code{NULL},
#'   snow parameters are assumed to be part of \code{ParamOptim}.
#' @param error_crit A string naming the error criterion function (e.g., \code{"KGE"}).
#'   Must be compatible with \code{\link{calc_hydroGOF}} or from the \pkg{hydroGOF} package.
#' @param cal_maximize Logical. If \code{TRUE}, the calibration maximizes the objective function.
#' @param cal_q_transfo A string indicating how runoff should be transformed
#'   (see \code{\link{transfo_q}}).
#' @param lambda Optional. A numeric value or vector used for regularization or transformation.
#' @param do_transfo_param Logical. If \code{TRUE}, parameters are transformed from the unit hypercube
#'   to real-world values before simulation.
#' @param airGR_RunOptions Optional. Run options for \pkg{airGR} models.
#' @param airGR_RunOptions_snow_module Optional. Run options for the snow module if using \pkg{airGR}.
#'
#' @return A single numeric value representing the error criterion to be minimized or maximized.
#'
#' @note
#' \itemize{
#'   \item If \code{Qsim} is entirely \code{NA}, a large penalty value (\code{+/- 1e10}) is returned.
#'   \item If \code{Qsim} is all zeros and the criterion is \code{KGE}, the asymptotic value \code{1 - sqrt(3)} is returned.
#'   \item \code{Qsim} is coerced to numeric to ensure compatibility with \code{Qobs}, especially for models like TUW.
#'   \item Future improvements may include spatially explicit versions and refined NA handling.
#' }
#'
#' @seealso \code{\link{calibrate_model}}, \code{\link{create_input}}, \code{\link{load_meteo_data}},
#'   \code{\link{calc_hydroGOF}}, \code{\link{transfo_param}}, \code{\link{simulate_model}}, \code{\link{simulate_snow}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cal_results <- DEoptim::DEoptim(
#'   fn = optim_fn,
#'   lower = lower,
#'   upper = upper,
#'   control = DEoptim::DEoptim.control(NP = 50, itermax = 100),
#'   hydro_data = hydro_data,
#'   split_indices = split_indices,
#'   model = "GR4J",
#'   input = input,
#'   error_crit = "KGE",
#'   cal_maximize = TRUE,
#'   cal_q_transfo = "none",
#'   do_transfo_param = TRUE
#' )
#' }
optim_fn <- function(ParamOptim, hydro_data, split_indices, model, input,
                     snow_module = NULL, snow_input, snow_parameters = NULL,
                     error_crit, cal_maximize, cal_q_transfo, lambda, do_transfo_param,
                     airGR_RunOptions = NULL, airGR_RunOptions_snow_module = NULL) {
  tryCatch(
    error = function(cnd) {
      # todo: see if this is an error below or just complicated
      # error_crit_val <- ifelse(cal_maximize, -1e10, 1e10)
      # ifelse(cal_maximize, return(-error_crit_val), return(error_crit_val))

      return(1e10)
    },
    {
      # transform Qobs if needed and subset only calibration period
      # todo: replace hydro_data with Qobs or even with input$Q. attention with special airGR calibration in calibrate_model
      Qobs <- transfo_q(hydro_data$BasinObs$Qmm[split_indices$ind_cal], cal_q_transfo, lambda)

      # re transform parameters to real values if needed
      if (do_transfo_param) {
        add_snow_par <- FALSE
        if (is.null(snow_parameters)) {
          add_snow_par <- TRUE
        }
        ParamOptim <- transfo_param(ParamOptim, "TR", model, snow_module, add_snow_par, cal_parameter = cal_par)
      }

      # simulate for warm up and calibration
      ind_warm_cal <- c(split_indices$ind_warm, split_indices$ind_cal)

      # model snow
      if (!is.null(snow_module)) {
        # adapt precipitation
        # ensure that P is present in input
        if (!"P" %in% names(input)) stop("P is not an entry of a list like input")
        nof_param_snow <- cal_par[[snow_module]][["nof_param"]]

        # override input$P
        input$P <- rep(NA, length(input$P))
        # either simulate with fixed or with free snow parameters
        if (is.null(snow_parameters)) {
          input$P[ind_warm_cal] <- simulate_snow(snow_module, ParamOptim[1:nof_param_snow], snow_input, ind_warm_cal,
            airGR_RunOptions = airGR_RunOptions_snow_module
          )$surface_water_runoff
          # remaining runoff model parameters
          ParamOptim <- ParamOptim[(nof_param_snow + 1):length(ParamOptim)]
        } else {
          input$P[ind_warm_cal] <- simulate_snow(snow_module, snow_parameters, snow_input, ind_warm_cal,
            airGR_RunOptions = airGR_RunOptions_snow_module
          )$surface_water_runoff
        }
      }


      # transform Qsim and take only values from calibration period
      Qsim <- transfo_q(suppressMessages(simulate_model(
        model, ParamOptim, input, ind_warm_cal,
        airGR_RunOptions = airGR_RunOptions
      )$Qsim[-seq_along(split_indices$ind_warm)]), cal_q_transfo, lambda) %>%
        # Qsim  converted better to numeric
        as.numeric()


      # deal with solely NA from Qsim
      if (all(is.na(Qsim))) {
        error_crit_val <- ifelse(cal_maximize, -1e10, 1e10)
      } else {
        error_crit_val <- calc_crit(error_crit, Qsim, Qobs)

        # set to asymptotic value, i.e. 1 - sqrt(3)
        if (is.na(error_crit_val)) {
          # todo: figure out for KGEtang the asymptotic value
          if (!stringr::str_starts(error_crit, "KGE")) stop("optim_fn: NA in error criterion")
          error_crit_val <- 1 - sqrt(3)
          if (any(Qsim > 0)) {
            warning("optim_fn: something else causes NAs as well in error_crit_val")
          }
        }
      }

      ifelse(cal_maximize, return(-error_crit_val), return(error_crit_val))
    }
  )
}

#' Calibrate a Hydrological Model
#'
#' Performs calibration of a hydrological model using a specified optimization algorithm.
#' Supports models from the \pkg{airGR}, \pkg{TUWmodel}, \pkg{hydromad}, and \pkg{topmodel}
#' packages, with optional snow module integration. The function supports both native
#' calibration routines (e.g., \code{Calibration_Michel} for \pkg{airGR}) and general-purpose
#' optimizers (e.g., \code{DEoptim}, \code{hydroPSO}, \code{malschains}).
#'
#' @param hydro_data A list or data frame containing observed runoff and meteorological data,
#'   typically loaded using \code{\link{load_meteo_data}}.
#' @param split_indices A list of index vectors (e.g., from \code{\link{split_data_set}})
#'   indicating warm-up and calibration periods.
#' @param model A string specifying the hydrological model to calibrate.
#'   For a complete list see table in `vignette("model_overview")`.
#' @param input A list of model input data, typically created using \code{\link{create_input}}.
#' @param snow_module Optional. A string specifying the snow module (currently \code{"CemaNeige"} and \code{"TUWsnow"}).
#' @param snow_input Optional. Input data for the snow module.
#' @param snow_parameters Optional. A vector of fixed snow parameters. If \code{NULL},
#'   snow parameters are assumed to be part of the calibration.
#' @param error_crit_transfo A string combining the error criterion and runoff transformation,
#'   separated by \code{"__"} (e.g., \code{"KGE__none"}). Optionally, a third value (e.g., lambda)
#'   can be included for transformations like Box-Cox.
#' @param cal_maximize Logical. If \code{TRUE}, the calibration maximizes the objective function.
#' @param cal_fn A string specifying the calibration function. Supported options include
#'   \code{"Calibration_Michel"} (for \pkg{airGR} models), \code{"DEoptim"}, \code{"hydroPSO"},
#'   \code{"malschains"}, and other supported optimizers.
#'   For a complete list see table in `vignette("calibration_methods_overview")`.
#' @param do_transfo_param Logical. If \code{TRUE}, parameters are transformed to a unit hypercube
#'   before calibration.
#' @param cal_par A list of calibration settings specific to the chosen calibration function.
#'   Defaults to \code{default_cal_par}, but can be customized by the user.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{model_param}: Calibrated model parameters.
#'   \item \code{error_crit_transfo}: The error criterion and transformation used.
#'   \item \code{error_crit_val}: The final value of the error criterion.
#'   \item \code{cal_fn}: The calibration function used.
#'   \item \code{do_transfo_param}: Whether parameter transformation was applied.
#'   \item \code{duration}: Duration of the calibration process.
#'   \item \code{cal_par}: Calibration settings used.
#'   \item \code{more_info}: Additional model- or method-specific output.
#' }
#'
#' @note
#' \itemize{
#'   \item \code{Calibration_Michel} is only available for \pkg{airGR} models.
#'   \item If \code{Calibration_Michel} is used, parameters are assumed to be transformed.
#'   \item The function supports power and Box-Cox runoff transformations, with lambda optionally specified.
#'   \item Future improvements may simplify access to \code{cal_par} for end users.
#' }
#'
#' @seealso \code{\link{call_cal_fn}}, \code{\link{optim_fn}}, \code{\link{create_input}},
#'   \code{\link{load_meteo_data}}, \code{\link{split_data_set}}, \code{\link{transfo_q}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' calibration_results <- calibrate_model(
#'   hydro_data = hydro_data,
#'   split_indices = split_data_set(...),
#'   model = "GR4J",
#'   input = create_input(...),
#'   error_crit_transfo = "KGE__none",
#'   cal_maximize = TRUE,
#'   cal_fn = "DEoptim",
#'   do_transfo_param = TRUE,
#'   cal_par = default_cal_par
#' )
#' }
calibrate_model <- function(hydro_data, split_indices, model, input,
                            snow_module = NULL, snow_input = NULL, snow_parameters = NULL,
                            error_crit_transfo = "KGE__none", cal_maximize = TRUE,
                            cal_fn = "DEoptim", do_transfo_param = FALSE,
                            cal_par = default_cal_par) {
  start <- Sys.time()

  # split error_crit entry
  error_crit_transfo_split <- stringr::str_split(error_crit_transfo, "__", simplify = TRUE)

  if (!length(error_crit_transfo_split) %in% c(2, 3)) {
    stop("error_crit_transfo not correctly specified")
  }


  # error crit needs to be global for the optim_fn
  error_crit <- error_crit_transfo_split[1]
  cal_q_transfo <- error_crit_transfo_split[2]
  if (length(error_crit_transfo_split) == 3) {
    lambda <- error_crit_transfo_split[3] %>% as.numeric()
  } else {
    lambda <- NULL
  }

  if (cal_fn %in% c("DEoptim", "hydroPSO", "malschains", "steepest_descent") | any(stringr::str_starts(cal_fn, c("montecarlo", "nlminb", "optim")))) {
    # call an R package based calibration function
    cal_output <- call_cal_fn(
      cal_fn, hydro_data, split_indices, model, input, snow_module, snow_input, snow_parameters,
      error_crit, cal_maximize, cal_q_transfo, lambda, do_transfo_param, cal_par
    )

    # best value of chosen error criterion
    error_crit_val <- ifelse(cal_maximize, -cal_output$error_crit_val, cal_output$error_crit_val)

    model_param <- cal_output$model_param

    # additional output
    more_info <- cal_output$cal_results

    # Calibration Michel just for airGR models
  } else if (cal_fn == "Calibration_Michel") {
    if (get_family(model) == "airGR") {
      # this is applied to transformed parameter set by default as long a FuN_MOD is native, throw a warning
      if (!do_transfo_param) {
        warning("Calibration_Michel is automatically applied to transformed model parameters")
      }

      # adapt cal_q_transfo power transformations
      if (cal_q_transfo == "power") {
        cal_q_transfo <- lambda
        # only boxcox_santos is implemented, if user supplies boxcox, throw warning
      } else if (cal_q_transfo == "boxcoxsantos") {
        cal_q_transfo <- "boxcox"
        # only boxcox_santos is implemented, if user supplies boxcox, throw warning
      } else if (cal_q_transfo == "boxcox") {
        warning("only boxcoxsantos is implemented for Calibration_Michel")
      } else {
        # if "none" in cal_q_transfo, it needs to be "" here
        cal_q_transfo <- stringr::str_replace(cal_q_transfo, "none", "")
      }

      RunOptions <- airGR::CreateRunOptions(
        FUN_MOD = paste0("RunModel_", model),
        InputsModel = input,
        IndPeriod_Run = split_indices$ind_cal,
        IndPeriod_WarmUp = split_indices$ind_warm,
        Outputs_Sim = "Qsim"
      )

      InputsCrit <- airGR::CreateInputsCrit(
        FUN_CRIT = paste0("ErrorCrit_", error_crit),
        InputsModel = input,
        RunOptions = RunOptions,
        VarObs = "Q",
        Obs = hydro_data$BasinObs$Qmm[split_indices$ind_cal],
        transfo = cal_q_transfo
      )

      CalibOptions <- airGR::CreateCalibOptions(
        FUN_MOD = paste0("RunModel_", model),
        FUN_CALIB = cal_fn
      )

      cal_fn <- get(cal_fn)
      OutputsCalib <- cal_fn(
        InputsModel = input,
        RunOptions = RunOptions,
        InputsCrit = InputsCrit, CalibOptions = CalibOptions,
        FUN_MOD = paste0("RunModel_", model)
      )

      model_param <- OutputsCalib$ParamFinalR

      error_crit_val <- OutputsCalib$CritFinal

      more_info <- list(calibration = OutputsCalib, InputsModel = input)

      # for all other models return NULL and throw a warning
    } else {
      # skip as this is not available
      calibration_results <- NULL

      warning("Calibration_Michel is not available for model \"", model, "\".", call. = FALSE)

      return(calibration_results)
    }
  } else {
    stop("Invalid choice of calibration function: \"", cal_fn,
      "\". Choose one of the following: Calibration_Michel, DEoptim, hydroPSO, malschains, nlminb, optim, montecarlo_XXX.",
      call. = FALSE
    )
  }

  end <- Sys.time()

  # define what to output
  calibration_results <- list(
    model = model,
    snow_module = snow_module,
    model_param = model_param,
    preset_snow_parameters = !is.null(snow_parameters),
    cal_fn = cal_fn,
    error_crit_transfo = error_crit_transfo,
    error_crit_val = error_crit_val,
    cal_maximize = cal_maximize,
    do_transfo_param = do_transfo_param,
    duration = round(end - start, 3),
    cal_par = cal_par[[model]],
    more_info = more_info
  )

  return(calibration_results)
}


#' Simulate Snow Module
#'
#' Simulates snow accumulation and melt processes using a specified snow module.
#' Currently supports \code{"CemaNeige"} (from \pkg{airGR}) and \code{"TUWsnow"} (from \pkg{TUWmodel}).
#' Returns surface water runoff and other snow-related variables such as SWE, melt, and precipitation components.
#'
#' @param snow_module A string specifying the snow module to use. Supported options are
#'   \code{"CemaNeige"} and \code{"TUWsnow"}. \code{"snowsnow"} is not yet implemented.
#' @param model_param A numeric vector of snow module parameters.
#' @param input A list containing input data for the snow module, including precipitation (\code{P}),
#'   air temperature (\code{T}), and optionally evapotranspiration (\code{E}).
#' @param ind A vector of indices specifying the time steps to simulate. Defaults to the full time range.
#' @param airGR_RunOptions Optional. A pre-created \code{RunOptions} object for \pkg{airGR} models.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{surface_water_runoff}: Simulated liquid water output (rain + melt).
#'   \item \code{SWE}: Snow water equivalent.
#'   \item \code{psolid}: Solid precipitation.
#'   \item \code{pliquid}: Liquid precipitation.
#'   \item \code{melt}: Meltwater.
#'   \item \code{more_info}: A list with additional model-specific output.
#' }
#'
#' @note
#' \itemize{
#'   \item For \code{"CemaNeige"}, the function uses \pkg{airGR}'s \code{RunModel_CemaNeige}.
#'   \item For \code{"TUWsnow"}, the function uses \pkg{TUWmodel} and appends constant runoff parameters.
#'   \item The \code{"snowsnow"} module is not yet implemented due to missing runoff output.
#'   \item Multilayer or spatially distributed snowpack outputs are not yet supported.
#' }
#'
#' @seealso \code{\link{calibrate_model}}, \code{\link{simulate_model}}, \code{\link[airGR]{RunModel_CemaNeige}}, \code{\link[TUWmodel]{TUWmodel}}
#'
#' @export
simulate_snow <- function(snow_module, model_param, input, ind = seq_along(input[[1]]), airGR_RunOptions = NULL) {
  if (snow_module == "CemaNeige") {
    # speed up simulation runs if RunOptions were not always newly created
    if (is.null(airGR_RunOptions)) {
      ## preparation of the RunOptions object for the validation period
      RunOptions <- airGR::CreateRunOptions(
        FUN_MOD = paste0("RunModel_", snow_module),
        InputsModel = input, IndPeriod_Run = ind
      )
    } else {
      RunOptions <- airGR_RunOptions
    }

    model_fn <- get(paste0("RunModel_", snow_module))
    output_model <- model_fn(InputsModel = input, RunOptions = RunOptions, Param = model_param)

    # get simulated snowpack runoff in case of one layer
    # todo: think about multilayer output (i.e. distributed in space, not a multilayer snowpack)
    surface_water_runoff <- output_model$CemaNeigeLayers$Layer01$PliqAndMelt

    # create output list
    snow_module_results <- list(
      surface_water_runoff = surface_water_runoff,
      # todo: think about multilayer output (i.e. distributed in space, not a multilayer snowpack)
      SWE = output_model$CemaNeigeLayers$Layer01$SnowPack,
      psolid = output_model$CemaNeigeLayers$Layer01$Psol,
      pliquid = output_model$CemaNeigeLayers$Layer01$Pliq,
      melt = output_model$CemaNeigeLayers$Layer01$Melt,
      more_info = list(output_model = output_model)
    )
  } else if (snow_module == "TUWsnow") {
    # take lower for constant runoff parameters
    model_param_TUW <- cal_par[["TUW"]][["lower"]]
    # override with variable snow parameters
    model_param_TUW <- c(model_param, model_param_TUW[(length(model_param) + 1):length(model_param_TUW)])
    output_model <- TUWmodel::TUWmodel(
      prec = input$P[ind], airt = input$T[ind], ep = input$E[ind],
      area = 1, param = model_param_TUW
    )

    # get simulated runoff
    surface_water_runoff <- output_model$rain + output_model$melt

    # create output list
    snow_module_results <- list(
      surface_water_runoff = surface_water_runoff %>% as.numeric(),
      SWE = output_model$swe %>% as.numeric(),
      psolid = output_model$snow %>% as.numeric(),
      pliquid = output_model$rain %>% as.numeric(),
      melt = output_model$melt %>% as.numeric(),
      more_info = list(output_model = output_model)
    )
  } else if (snow_module == "snowsnow") {
    stop("not implemented yet, as the surface_water_runoff is not returned")

    output_model <- create_hydromad_model("snow", cal_par[["snow"]])

    model_param_hydromad <- cal_par$snow$lower
    # override with variable snow parameters
    names(model_param) <- cal_par$snowsnow$lower %>% names()
    model_param_hydromad <- c(model_param, model_param_hydromad[(length(model_param) + 1):length(model_param_hydromad)])

    output_model <- hydromad::update(output_model,
      newdata = zoo::read.zoo(input[ind, ]),
      newpars = model_param_hydromad,
      warmup = 0
    )

    # get simulated states
    output_model <- hydromad::predict(output_model, return_state = TRUE)

    # create output list
    snow_module_results <- list(
      # todo: find out what the surface water runoff is
      # surface_water_runoff = output_model$TF %>% as.numeric(),
      SWE = output_model$SWE %>% as.numeric(),
      more_info = list(output_model = output_model)
    )
  } else {
    stop(sprintf("the model choice %s is not implemented yet, only TUWsnow and CemaNeige are implemented", model))
  }

  return(snow_module_results)
}


#' Simulate a Hydrological Model
#'
#' Simulates discharge using a specified hydrological model over a given time period.
#' Supports models from the \pkg{airGR}, \pkg{TUWmodel}, \pkg{hydromad}, and \pkg{topmodel} packages.
#' Optionally includes snow-related outputs if the model supports it.
#'
#' @param model A string specifying the hydrological model to use. Supported models include
#'   \code{"TUW"}, \code{"topmodel"}, \code{"hydromad"} models (e.g., \code{"sacramento"}, \code{"cwi"}),
#'   and \pkg{airGR} models (e.g., \code{"GR4J"}, \code{"CemaNeigeGR4J"}).
#'   For a complete list see table in `vignette("model_overview")`.
#' @param model_param A numeric vector of model parameters specific to the chosen model.
#' @param input A list of model input data, typically created using \code{\link{create_input}}.
#'   Must include time series of precipitation (\code{P}), temperature (\code{T}), and optionally
#'   evapotranspiration (\code{E}), as well as spatial data like catchment area or topography.
#' @param ind A vector of indices specifying the time steps to simulate. Defaults to the full range.
#' @param Qobs Optional. A vector of observed discharge values to include in the output for comparison.
#' @param airGR_RunOptions Optional. A pre-created \code{RunOptions} object for \pkg{airGR} models.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{date}: Vector of simulation dates.
#'   \item \code{Qsim}: Simulated discharge.
#'   \item \code{Qobs}: Observed discharge (if provided).
#'   \item \code{SWE}: Snow water equivalent (if available).
#'   \item \code{psolid}: Solid precipitation (if available).
#'   \item \code{pliquid}: Liquid precipitation (if available).
#'   \item \code{melt}: Meltwater (if available).
#'   \item \code{more_info}: A list with model-specific output.
#' }
#'
#' @note
#' \itemize{
#'   \item For \pkg{airGR} models, the appropriate \code{RunModel_*} function is called.
#'   \item For \pkg{hydromad} models, missing dates are filled to ensure compatibility.
#'   \item Snow-related outputs are only available for models that simulate snow processes.
#'   \item If \code{Qobs} is provided, its length must match the length of \code{Qsim}.
#' }
#'
#' @seealso \code{\link{simulate_snow}}, \code{\link{create_input}}, \code{\link{calibrate_model}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' simulation_results <- simulate_model(
#'   model = "TUW",
#'   model_param = calibration_results$model_param,
#'   input = input,
#'   ind = split_indices$ind_cal
#' )
#' }
simulate_model <- function(model, model_param, input, ind = seq_along(input[[1]]), Qobs = NULL, airGR_RunOptions = NULL) {
  if (model == "TUW") {
    output_model <- TUWmodel::TUWmodel(
      prec = input$P[ind], airt = input$T[ind], ep = input$E[ind],
      area = 1, param = model_param
    )

    # get simulated runoff
    Qsim <- as.numeric(output_model$q)
    date <- input$date[ind] %>% as.Date()

    # snow information
    SWE <- as.numeric(output_model$swe)
    psolid <- as.numeric(output_model$snow)
    pliquid <- as.numeric(output_model$rain)
    melt <- as.numeric(output_model$melt)
  } else if (model == "topmodel") {
    model_param[11] <- 24 # time step expressed in hours
    output_model <- topmodel::topmodel(
      parameters = model_param,
      topidx = input$topidx,
      delay = input$delay,
      rain = input$P[ind] * 0.001, # in meters per time step
      ETp = input$E[ind] * 0.001, # in meters per time step
      verbose = FALSE
    ) # verbose = TRUE provides a list with entires Q, q0, qs, S, fex, Ea

    # get simulated runoff in mm/d
    # Qsim <- output_model$Q*1000
    Qsim <- output_model * 1000
    date <- input$date[ind] %>% as.Date()

    # snow information
    SWE <- NA
    psolid <- NA
    pliquid <- NA
    melt <- NA
  } else if (get_family(model) == "airGR") {
    # speed up simulation runs if RunOptions were not always newly created
    if (is.null(airGR_RunOptions)) {
      ## preparation of the RunOptions object for the validation period
      RunOptions <- airGR::CreateRunOptions(
        FUN_MOD = paste0("RunModel_", model),
        InputsModel = input, IndPeriod_Run = ind
      )
    } else {
      RunOptions <- airGR_RunOptions
    }

    model_fn <- get(paste0("RunModel_", model))
    output_model <- model_fn(InputsModel = input, RunOptions = RunOptions, Param = model_param)

    # get simulated runoff
    Qsim <- output_model$Qsim
    date <- input$DatesR[ind] %>% as.Date()

    # snow information
    if (stringr::str_detect(model, "CemaNeige")) {
      SWE <- output_model$CemaNeigeLayers$Layer01$SnowPack
      psolid <- output_model$CemaNeigeLayers$Layer01$Psol
      pliquid <- output_model$CemaNeigeLayers$Layer01$Pliq
      melt <- output_model$CemaNeigeLayers$Layer01$Melt
    } else {
      SWE <- NA
      psolid <- NA
      pliquid <- NA
      melt <- NA
    }
  } else if (get_family(model) == "hydromad") {
    init_model <- create_hydromad_model(model, cal_par[[model]], routing = cal_par[[model]]$routing)

    free_param_names <- init_model %>%
      hydromad::getFreeParsRanges() %>%
      names()

    # test if scale parameter need to be added for model cwi
    if (model == "cwi" & length(model_param) == length(free_param_names) + 1) {
      free_param_names <- c(free_param_names, "scale")
    }

    names(model_param) <- free_param_names

    # the date vector for output, put on top as input may change later
    date <- input$date[ind] %>% as.Date()

    # all models except of sacramento do not like missing dates (i.e. a real gap in the dates),
    # mimic behaviour of other models which do not care about dates at all and join dates and simulated values at the end
    if (model != "sacramento") {
      # fill NAs at the end of the data frame
      dates_with_gaps <- input$date

      # Create a sequence of dates from min to max with no gaps, and shorter P and T
      input_list <- list(
        date = seq(min(input$date), max(input$date), by = "day"),
        P = input$P, T = input$T, E = input$E
      )

      # fill NA for shorter columns
      fill_na <- function(data_list) {
        # the column with maximum length
        max_length <- max(purrr::map_int(data_list, length))
        # looping over all list entries (i.e. vectors and fill NA at the end)
        df <- data_list %>%
          map_if(~ length(.x) < max_length, ~ c(.x, rep(NA, max_length - length(.x)))) %>%
          bind_cols()
      }

      # data frame with complete dates, values at incorrect dates, NAs at end
      input <- input_list %>% fill_na()
    }

    fitted_model <- hydromad::update(init_model,
      newdata = zoo::read.zoo(input[ind, ]),
      newpars = model_param,
      warmup = 0
    )

    # get simulated runoff
    Qsim <- hydromad::fitted(fitted_model) %>% as.numeric()

    # todo: get snow information from "snow"
    # snow information
    SWE <- NA
    psolid <- NA
    pliquid <- NA
    melt <- NA

    output_model <- list()
    output_model$init_model <- init_model
    output_model$fitted_model <- fitted_model
  } else {
    stop(sprintf("the model choice %s is not implemented yet", model))
  }

  if (!is.null(Qobs)) {
    if (length(Qobs) != length(Qsim)) {
      stop(sprintf("Qobs and Qsim do not have the same length."))
    }
  }

  # create output list
  simulation_results <- list(
    date = date,
    Qsim = Qsim,
    Qobs = Qobs,
    SWE = SWE,
    psolid = psolid,
    pliquid = pliquid,
    melt = melt,
    more_info = list(output_model = output_model)
  )

  return(simulation_results)
}


#' Merge Snow and Runoff Simulation Results
#'
#' Combines the output of a snow module simulation with the results of a hydrological
#' runoff simulation. Adds snow-related variables and snow model metadata to the
#' runoff simulation results.
#'
#' @param simulation_results A list returned by \code{\link{simulate_model}}, containing
#'   simulated runoff and related metadata.
#' @param snow_module_results A list returned by \code{\link{simulate_snow}}, containing
#'   snow-related outputs such as SWE, melt, and precipitation components.
#'
#' @return A list identical to \code{simulation_results}, but with additional elements:
#' \itemize{
#'   \item \code{SWE}: Snow water equivalent.
#'   \item \code{psolid}: Solid precipitation.
#'   \item \code{pliquid}: Liquid precipitation.
#'   \item \code{melt}: Meltwater.
#'   \item \code{surface_water_runoff}: Combined snow runoff.
#'   \item \code{more_info$snow_module}: Snow model output metadata.
#' }
#'
#' @seealso \code{\link{simulate_model}}, \code{\link{simulate_snow}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' merged_results <- merge_snow_runoff_sim(simulation_results, snow_module_results)
#' }
merge_snow_runoff_sim <- function(simulation_results, snow_module_results) {
  # add to simulation_results
  simulation_results$SWE <- snow_module_results$SWE
  simulation_results$psolid <- snow_module_results$psolid
  simulation_results$pliquid <- snow_module_results$pliquid
  simulation_results$melt <- snow_module_results$melt
  simulation_results$surface_water_runoff <- snow_module_results$surface_water_runoff
  simulation_results$more_info$snow_module <- snow_module_results$more_info$output_model

  return(simulation_results)
}


#' Subset Simulation Results by Index
#'
#' This function subsets the elements of a simulation results list to the specified indices.
#' It excludes the `"more_info"` field by default but can optionally retain it.
#'
#' @param ind A vector of indices to subset the simulation results.
#' @param simulation_results A list containing simulation results, where each element (except `"more_info"`) is assumed to be indexable.
#' @param keep_more_info Logical; if `TRUE`, the `"more_info"` field is retained in the output. Default is `FALSE`.
#'
#' @return A list of simulation results subset to the specified indices. If `keep_more_info` is `TRUE`, the `"more_info"` field is included unchanged.
#'
#' @details
#' This function assumes that all elements of `simulation_results`, except `"more_info"`, can be subset using the provided indices.
#' A future enhancement could include a check to ensure that each field is indeed subsettable.
#'
#' @examples
#' \dontrun{
#' results <- list(
#'   data = list(1:10, 11:20, 21:30),
#'   metrics = list(a = 1:3, b = 4:6),
#'   more_info = list(description = "Full simulation run")
#' )
#' subset_simulations(1:2, results, keep_more_info = TRUE)
#' }
#'
#' @export
subset_simulations <- function(ind, simulation_results, keep_more_info = FALSE) {
  # subset all elements except "more info"
  elements <- names(simulation_results)
  simulation_results_subset <- purrr::map(simulation_results[elements[elements != "more_info"]], ~ (.x[ind]))

  # add more info again
  if (keep_more_info) {
    simulation_results_subset$more_info <- simulation_results$more_info
  }

  return(simulation_results_subset)
}


#' Validate model
#'
#' Calculates validation measures for different transformation types
#'
#' @param Qsim vector with simulated runoff
#' @param Qobs vector with observed runoff
#' @param val_crit_transfo a vector of strings specifying validation criteria and a
#'   runoff transformation separated by a \code{"__"}. Supported are validation
#'   criteria from the \code{\link{hydroGOF}} package usable by the
#'   \code{\link{calc_hydroGOF}} function, for supported runoff transformations
#' please refer to \code{\link{transfo_q}}
#'
#' @return a long data frame with columns \code{crit} indicating the used
#'   validation criterion, \code{transfo} for the used runoff transformation and
#'   \code{value}.
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' validate_model(
#'   1:10, seq(2, 11),
#'   c("KGE__log", "NSE__inv", "VE__none", "pbias__none")
#' )
validate_model <- function(Qsim, Qobs, val_crit_transfo = "KGE__none") {
  val_crit_transfo_split <- stringr::str_split(val_crit_transfo, "__", simplify = TRUE)

  if (!ncol(val_crit_transfo_split) %in% c(2, 3)) stop("val_crit_transfo not correctly specified")

  val_crit <- val_crit_transfo_split[, 1]
  val_q_transfo <- val_crit_transfo_split[, 2]
  if (ncol(val_crit_transfo_split) == 3) {
    val_lambda <- val_crit_transfo_split[, 3] %>% as.numeric()
  } else {
    # here it is important to have NA and not NULL
    val_lambda <- NA
  }
  names(val_crit) <- val_crit_transfo

  # loop over validation criteria and q transformations
  # it is required to delete looping variables
  rm(crit, q_transfo, lambda) %>% suppressWarnings()
  validation_results <- purrr::pmap_df(
    list(crit = val_crit, q_transfo = val_q_transfo, lambda = val_lambda),
    \(crit, q_transfo, lambda) calc_crit(crit, transfo_q(Qsim, q_transfo, lambda), transfo_q(Qobs, q_transfo, lambda))
  ) %>%
    # a long data_frame
    tidyr::pivot_longer(cols = everything(), values_to = "value") %>%
    # separate the name columns in crit and transfo, too few align start because sometimes lambda is not present and  NA is inserted
    tidyr::separate_wider_delim(name, names = c("crit", "transfo", "lambda"), delim = "__", too_few = "align_start")

  # # delete lambda column if completely NA,
  # # no preferred as save_cal_val_plot has hardcoded text
  # if (all(is.na(val_lambda))) {
  #   validation_results <- validation_results %>% dplyr::select(-lambda)
  # }

  return(validation_results)
}


#' Validate Hydrological Model Performance
#'
#' Calculates validation metrics for simulated versus observed runoff using various transformation types.
#'
#' @param Qsim A numeric vector of simulated runoff values.
#' @param Qobs A numeric vector of observed runoff values.
#' @param val_crit_transfo A character vector specifying validation criteria and runoff transformations,
#'   separated by `"__"`. Optionally, a third part can specify a lambda parameter. Supported criteria
#'   are those from the \code{\link[hydroGOF]{hydroGOF}} package and must be compatible with
#'   \code{\link{calc_hydroGOF}}. Supported transformations are described in \code{\link{transfo_q}}.
#'
#' @return A long-format data frame with columns:
#'   \itemize{
#'     \item \code{crit} – the validation criterion used (e.g., "KGE", "NSE").
#'     \item \code{transfo} – the runoff transformation applied (e.g., "log", "inv", "none").
#'     \item \code{lambda} – the lambda parameter used for transformation, if applicable.
#'     \item \code{value} – the resulting validation metric value.
#'   }
#'
#' @details
#' The function splits each entry in `val_crit_transfo` into components: validation criterion,
#' transformation type, and optionally a lambda value. It applies the transformation to both
#' `Qsim` and `Qobs`, computes the specified metric, and returns the results in a tidy format.
#'
#' @importFrom magrittr %>%
#' @importFrom stringr str_split
#' @importFrom tidyr pivot_longer separate_wider_delim
#' @importFrom purrr pmap_df
#'
#' @examples
#' validate_model(
#'   Qsim = 1:10,
#'   Qobs = seq(2, 11),
#'   val_crit_transfo = c("KGE__log", "NSE__inv", "VE__none", "pbias__none")
#' )
#'
#' @export
calc_validation_results <- function(ind, period_name, col_name = "period",
                                    Qsim, Qobs, val_crit_transfo = "KGE__none") {
  validation_results <- validate_model(Qsim[ind], Qobs[ind], val_crit_transfo) %>%
    dplyr::mutate("{col_name}" := period_name)

  return(validation_results)
}


#' Calculate Subseasonal Validation Metrics
#'
#' Computes validation metrics for specified subseasonal periods within a hydrological dataset.
#' For each named period in `val_subseason`, the function subsets the data using the provided indices
#' and calculates performance metrics using \code{\link{calc_validation_results}}.
#'
#' @param val_subseason A named list where each element is a character vector of two-digit month codes
#'   (e.g., `"06"`, `"07"`) defining the months for each subseasonal period.
#' @param dates A vector of dates (e.g., of class `Date` or `POSIXct`) corresponding to the time series.
#' @param ind A vector of indices used to subset the hydrological data.
#' @param period_name A string indicating the name of the period (e.g., `"calibration"` or `"validation"`),
#'   which will be added as a column in the output.
#' @param col_name A string specifying the name of the additional column to be added to the output
#'   (default is `"period"`).
#' @param Qsim A numeric vector of simulated runoff values.
#' @param Qobs A numeric vector of observed runoff values.
#' @param val_crit_transfo A character vector specifying validation criteria and runoff transformations,
#'   separated by `"__"`. Optionally, a third part can specify a lambda parameter. Supported criteria
#'   are those from the \code{\link[hydroGOF]{hydroGOF}} package and must be compatible with
#'   \code{\link{calc_hydroGOF}}. Supported transformations are described in \code{\link{transfo_q}}.
#'
#' @return A data frame similar to the output of \code{\link{validate_model}}, but with two additional
#'   columns: one for the subseasonal period (e.g., `"spring"`, `"summer"`) and one for the period name
#'   (e.g., `"calibration"`), as specified by `col_name` and `period_name`.
#'
#' @seealso \code{\link{calc_validation_results}}, \code{\link{validate_model}},
#'   \code{\link{calc_hydroGOF}}, \code{\link{transfo_q}}
#'
#' @examples
#' perf_cal <- calc_subseasonal_validation_results(
#'   val_subseason = list(
#'     spring = c("02", "03", "04", "05"),
#'     summer = c("06", "07", "08", "09")
#'   ),
#'   dates = hydro_data$BasinObs$DatesR,
#'   ind = split_indices$ind_cal,
#'   period_name = "calibration",
#'   col_name = "period",
#'   Qsim = simulation_results$Qsim,
#'   Qobs = Qobs,
#'   val_crit_transfo = c(
#'     "KGE__none", "NSE__none", "VE__none", "pbias__none",
#'     "KGE__inv", "NSE__inv",
#'     "KGE__sqrt", "NSE__sqrt"
#'   )
#' )
#'
#' @export
calc_subseasonal_validation_results <- function(val_subseason, dates, ind, period_name, col_name = "period",
                                                Qsim, Qobs, val_crit_transfo = "KGE__none") {
  ind_list <- purrr::map(
    val_subseason,
    ~ find_monthly_indices(dates, .x, ind)
  )
  ind_list$all <- ind

  validation_results <- purrr::imap_dfr(
    ind_list, calc_validation_results, "season",
    Qsim, Qobs,
    val_crit_transfo
  ) %>%
    dplyr::mutate("{col_name}" := period_name)

  return(validation_results)
}


#' Write ASCII Summary of Calibration and Validation Results
#'
#' Writes a plain text (ASCII) file summarizing model calibration parameters and validation results.
#' The output can be written in either tab-separated or fixed-width format.
#'
#' @param file A string specifying the path to the output file.
#' @param calibration_results A list containing the calibration results from \code{\link{calibrate_model}}.
#'   Only the vector of calibrated model parameters (\code{model_param}) is written.
#' @param validation_results A data frame containing validation results, typically from
#'   \code{\link{validate_model}}.
#' @param equally_spaced Logical; if `TRUE` (default), attempts to write a fixed-width formatted file.
#'   If `FALSE` or if fixed-width writing fails, a tab-separated file is written instead.
#'
#' @return A logical value indicating whether the file was successfully written.
#'
#' @details
#' If `equally_spaced = TRUE`, the function attempts to write a fixed-width formatted file using
#' \code{\link[gdata]{write.fwf}}. If an error occurs during this process, it falls back to writing
#' a tab-separated file using \code{\link[readr]{write_tsv}}.
#'
#' @importFrom readr write_tsv
#' @importFrom dplyr mutate_if arrange
#' @importFrom gdata write.fwf
#'
#' @examples
#' write_ascii(
#'   file = "results.txt",
#'   calibration_results = calibration_results,
#'   validation_results = validation_results
#' )
#'
#' @export
write_ascii <- function(file, calibration_results, validation_results, equally_spaced = TRUE) {
  if (!equally_spaced) {
    readr::write_tsv(dplyr::mutate_if(validation_results, is.numeric, round, digits = 3), file)
  } else {
    # try to write a equally spaced file, and write the non-equally spaced file if something goes wrong
    # todo is there a way to not duplicate the !equally spaced code?
    tryCatch(
      error = function(cnd) {
        readr::write_tsv(dplyr::mutate_if(validation_results, is.numeric, round, digits = 3), file)
      },
      {
        library(gdata)

        validation_results <- dplyr::arrange(validation_results, season)

        gdata::write.fwf(data.frame(t(names(validation_results))),
          file = file, width = rep(15, ncol(validation_results)), colnames = FALSE
        )
        gdata::write.fwf(as.data.frame(dplyr::mutate_if(validation_results, is.numeric, round, digits = 3)),
          file = file, width = rep(15, ncol(validation_results)), colnames = FALSE, append = TRUE
        )
        write(" ", file = file, append = TRUE)
        gdata::write.fwf(matrix(round(calibration_results$model_param, digits = 3), nrow = 1),
          file = file, append = TRUE
        )
      }
    )
  }
  return(TRUE)
}


#' Save Calibration and Validation Plot
#'
#' Creates and saves a plot comparing observed and simulated runoff over the calibration and validation periods.
#' The plot includes performance metrics such as KGE, NSE, and percent bias, and displays them for each period.
#'
#' @param file A string specifying the filename for the saved plot (e.g., `"cal_val.pdf"`).
#' @param BasinObs A data frame containing observed runoff and corresponding dates, typically from
#'   \code{\link{load_meteo_data}}. Must include columns \code{Qmm} and \code{DatesR}.
#' @param Qsim A numeric vector of simulated runoff values.
#' @param split_indices A list of indices from \code{\link{split_data_set}}, containing elements
#'   \code{ind_cal} and \code{ind_val} for calibration and validation periods, respectively.
#'
#' @return A logical value indicating whether the plot was successfully saved.
#'
#' @details
#' The function generates a two-panel plot showing observed and simulated runoff for both calibration
#' and validation periods. It also computes and displays selected validation metrics using
#' \code{\link{calc_validation_results}}. The metrics are shown as annotations on the plot.
#'
#' @note Future improvements could include adding a seasonal rolling mean to the plot.
#'
#' @importFrom purrr imap map_df iwalk
#' @importFrom dplyr mutate filter bind_rows mutate_if
#' @importFrom graphics plot lines mtext par
#' @importFrom grDevices pdf dev.off
#'
#' @examples
#' save_cal_val_plot(
#'   file = "cal_val.pdf",
#'   BasinObs = BasinObs,
#'   Qsim = simulation_results$Qsim,
#'   split_indices = split_indices
#' )
#'
#' @export
save_cal_val_plot <- function(file, BasinObs, Qsim, split_indices) {
  Qobs <- BasinObs$Qmm
  plot_df <- data.frame(dates = BasinObs$DatesR, Qsim = Qsim, Qobs = Qobs)

  ind_list <- list()
  ind_list$Calibration <- split_indices$ind_cal
  ind_list$Validation <- split_indices$ind_val

  val_crit_transfo <- c("pbias__none", "KGE__none", "NSE__none", "KGEtang__log", "KGE__power__-0.5")
  val_crit_label <- c("%", "", "", "", "")

  validation_results <- purrr::imap(
    ind_list, calc_validation_results,
    "period", Qsim, Qobs, val_crit_transfo
  ) %>%
    purrr::map_df(dplyr::bind_rows) %>%
    # rename none and pbias and lambda = NA and power and KGEtang
    dplyr::mutate(transfo = replace(transfo, transfo == "none", "")) %>%
    dplyr::mutate(crit = replace(crit, crit == "pbias", "dV")) %>%
    dplyr::mutate(crit = replace(crit, crit == "KGEtang", "KGE2")) %>%
    dplyr::mutate(lambda = replace(lambda, is.na(lambda), "")) %>%
    dplyr::mutate(transfo = replace(transfo, transfo == "power", ""))

  adj_val <- seq(0, 1, 1 / (length(val_crit_transfo) - 1))

  plot_single_Qsim_Qobs <- function(plot_df, validation_results, adj_val) {
    plot(plot_df$dates, plot_df$Qobs,
      type = "l", xlab = "", ylab = "Discharge [mm/d]",
      main = validation_results[1, 5], col = "blue"
    )
    lines(plot_df$dates, plot_df$Qsim, col = "red")
    for (i in 1:nrow(validation_results)) {
      mtext(paste0(
        validation_results[i, 1], validation_results[i, 2], validation_results[i, 3], " = ",
        validation_results[i, 4]$value %>% signif(digits = 2), val_crit_label[i]
      ), adj = adj_val[i], cex = .8)
    }
  }

  # map solution
  filter_val_plot <- function(ind, season_str, validation_results) {
    validation_results <- dplyr::filter(validation_results, period == season_str) %>%
      dplyr::mutate_if(is.numeric, round, digits = 2)
    plot_single_Qsim_Qobs(plot_df[ind, ], validation_results, adj_val)
  }

  pdf(file = file)
  par(mfrow = c(2, 1), mar = c(3, 4, 3, 2))
  purrr::iwalk(ind_list, filter_val_plot, validation_results)
  dev.off()

  return(TRUE)
}


#' Save airGR Diagnostic Plots
#'
#' Generates and saves diagnostic plots for hydrological model simulations using the \pkg{airGR} plotting functions.
#' This function supports both airGR and non-airGR models by converting simulation results into a compatible format.
#'
#' @param file A string specifying the filename for the saved plot (e.g., `"airGR_plot.pdf"`).
#' @param model A string indicating the model name (e.g., `"CemaNeigeGR4J"`), used to determine if the model is from the airGR family.
#'   For a complete list see table in `vignette("model_overview")`.
#' @param simulation_results A list containing simulation outputs, typically from \code{\link{simulate_model}}.
#'   Must include at least \code{Qsim}, \code{Qobs}, and optionally snow-related variables like \code{SWE}, \code{psolid}, and \code{pliquid}.
#' @param ind A vector of indices specifying the time period to plot. Defaults to the full time range.
#' @param hydro_data Optional list containing \code{BasinObs} and \code{BasinInfo}, required if the model is not from the airGR family.
#'
#' @return A logical value indicating whether the plot was successfully saved.
#'
#' @details
#' If the model is not from the airGR family, the function reconstructs an \code{airGR::OutputsModel} object
#' using \code{\link{create_input}} and \code{\link{simulate_model}}, and overrides it with the provided simulation results.
#' If snow-related outputs are present, additional snow plots are included.
#'
#' @note
#' \itemize{
#'   \item Currently supports only one CemaNeige snow layer.
#'   \item A future version may support plotting without requiring \code{Qobs}.
#' }
#'
#' @importFrom grDevices pdf dev.off
#'
#' @seealso \code{\link{simulate_model}}, \code{\link[airGR]{plot}}, \code{\link{create_input}}
#'
#' @examples
#' save_airGR_plot(
#'   file = "airGR_plot.pdf",
#'   model = "CemaNeigeGR4J",
#'   simulation_results = simulation_results,
#'   hydro_data = hydro_data
#' )
#'
#' @export
save_airGR_plot <- function(file, model, simulation_results, ind = seq_along(simulation_results$date), hydro_data = NULL) {
  if (get_family(model) != "airGR") {
    # create an airGR::OutputsModel
    airGR_input <- create_input("CemaNeigeGR4J", hydro_data$BasinObs, hydro_data$BasinInfo) %>% suppressWarnings()
    airGR_OutputsModel <- simulate_model("CemaNeigeGR4J", default_cal_par$CemaNeigeGR4J$lower, airGR_input, Qobs = simulation_results$Qobs)$more_info$output_model %>%
      suppressWarnings()
    # override values
    airGR_OutputsModel$Qsim <- simulation_results$Qsim
    airGR_OutputsModel$Qobs <- simulation_results$Qobs
    which_plots <- c("Flows", "Regime", "CumFreq", "CorQQ")

    # add and plot also snow results, if no snow results, simulate_model adds a single NA
    if (!is.na(simulation_results$SWE[1])) {
      # override values for plotting
      airGR_OutputsModel$CemaNeigeLayers$Layer01$SnowPack <- simulation_results$SWE
      airGR_OutputsModel$CemaNeigeLayers$Layer01$Psol <- simulation_results$psolid
      airGR_OutputsModel$CemaNeigeLayers$Layer01$Pliq <- simulation_results$pliquid
      which_plots <- "synth"
      which_plots <- c("Precip", "SnowPack", "Flows", "Error", "Regime", "CumFreq", "CorQQ")
    }
    # for airGR models just take the already created output model
  } else {
    airGR_OutputsModel <- simulation_results$more_info$output_model
    which_plots <- c("Precip", "SnowPack", "Flows", "Error", "Regime", "CumFreq", "CorQQ")
  }

  # plot pdf
  pdf(file = file)
  plot(airGR_OutputsModel, Qobs = simulation_results$Qobs, IndPeriod_Plot = ind, which = which_plots)
  dev.off()
}


#' Set Calibration Parameter Value
#'
#' Updates a specific entry in a nested calibration parameter list for a given model.
#'
#' @param model A string specifying the model name (e.g., `"TUW"`), used to access the corresponding sublist in `cal_par`.
#'   For a complete list see table in `vignette("model_overview")`.
#' @param setting_name_value A string representing the setting to update, in the format `"sublist$parameter = value"`.
#'   This string is parsed and evaluated to modify the calibration settings.
#' @param cal_par A list containing calibration settings for one or more models.
#'
#' @return The updated calibration settings list.
#'
#' @details
#' This function uses non-standard evaluation to dynamically update a nested element in the `cal_par` list.
#' It is useful for programmatically modifying calibration settings without manually navigating the list structure.
#'
#' @examples
#' cal_par_updated <- set_cal_par("TUW", "DEoptim$itermax = 5", cal_par)
#'
#' @export
set_cal_par <- function(model, setting_name_value, cal_par) {
  subset_list_expr <- rlang::parse_expr(paste0("cal_par$", model, "$", setting_name_value))

  eval(subset_list_expr)

  return(cal_par)
}


#' Create Default Calibration Settings for airGR Models
#'
#' Generates a list of default calibration settings for a specified model from the airGR family.
#' This includes parameter bounds and optimizer configurations for DEoptim, malschains, and hydroPSO.
#'
#' @param model A string specifying the airGR model (e.g., `"GR4J"`, `"CemaNeigeGR4J"`).
#'   For a complete list see table in `vignette("model_overview")`.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{lower}, \code{upper} – numeric vectors of parameter bounds.
#'   \item \code{nof_param} – number of parameters.
#'   \item \code{DEoptim}, \code{malschains}, \code{hydroPSO} – optimizer settings.
#'   \item \code{has_snow_module} – logical indicating whether the model includes a snow module.
#' }
#'
#' @details
#' The function uses \code{\link[airGR]{CreateCalibOptions}} to extract parameter bounds and sets
#' default configurations for several optimizers. It also detects whether the model includes a snow
#' module based on the model name.
#'
#' @note
#' This function only supports models from the airGR family. An error is thrown otherwise.
#'
#' @importFrom airGR CreateCalibOptions
#' @importFrom stringr str_detect
#'
#' @examples
#' set_airGR_par("GR4J")
#' set_airGR_par("CemaNeigeGR4J")
#'
#' @export
set_airGR_par <- function(model) {
  # todo check if this is the correct way to include airGR, otherwise the lines will not work, @import airGR does not work
  # model_fn <- get(paste0("RunModel_", model))
  # CalibOptions <- airGR::CreateCalibOptions(FUN_MOD = model_fn)
  if (get_family(model) == "airGR") {
    # call CreateCalibOptions from airGR package to retrieve upper and lower boundaries
    model_fn <- get(paste0("RunModel_", model))
    CalibOptions <- airGR::CreateCalibOptions(FUN_MOD = model_fn)

    out <- list()
    out$lower <- CalibOptions$SearchRanges[1, ]
    out$upper <- CalibOptions$SearchRanges[2, ]
    out$nof_param <- length(out$lower)

    # DEoptim
    out$DEoptim$NP <- out$nof_param * 10
    out$DEoptim$itermax <- 200

    # malschains
    out$malschains$maxEvals <- 2000

    # hydroPSO
    out$hydroPSO$control$write2disk <- FALSE
    out$hydroPSO$control$verbose <- FALSE
    # lower number of variables, thus leave it to standard
    out$hydroPSO$control$npart <- 80
    out$hydroPSO$control$maxit <- 50
    out$hydroPSO$control$reltol <- 1E-10

    # define if snow module is included
    if (stringr::str_detect(model, "CemaNeige")) {
      out$has_snow_module <- TRUE
    } else {
      out$has_snow_module <- FALSE
    }


    return(out)
  } else {
    stop("not an airGR model")
  }
}


#' Create Default Calibration Settings for hydromad Models
#'
#' Generates a list of default calibration settings for a given hydromad model, including parameter bounds
#' and optimizer configurations for DEoptim, malschains, and hydroPSO. Supports both SMA-only and SMA-routing model combinations.
#'
#' @param model A string specifying the soil moisture accounting (SMA) model (e.g., `"gr4j"`, `"sacramento"`, `"snow"`),
#'   or a pre-defined \code{hydromad} model object. For a complete list see table in `vignette("model_overview")`.
#' @param routing A string specifying the routing model (e.g., `"expuh"`, `"lambda"`). Only used if \code{model} is a string.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{lower}, \code{upper} – numeric vectors of parameter bounds.
#'   \item \code{nof_param} – number of parameters.
#'   \item \code{routing} – the routing model used.
#'   \item \code{DEoptim}, \code{malschains}, \code{hydroPSO} – optimizer settings.
#'   \item \code{has_snow_module} – logical indicating whether the model includes a snow module.
#' }
#'
#' @details
#' If a model string is provided, a temporary hydromad model is created using the specified SMA and routing.
#' Parameter bounds are extracted using \code{\link[hydromad]{getFreeParsRanges}}.
#' Routing parameters are only added for models that support routing.
#'
#' @importFrom hydromad hydromad update getFreeParsRanges
#' @importFrom purrr map_dbl
#' @importFrom methods is
#'
#' @examples
#' set_hydromad_par("gr4j")
#' set_hydromad_par("sacramento")
#' set_hydromad_par(hydromad(DATA = NULL, sma = "snow", routing = "expuh"))
#'
#' @export
set_hydromad_par <- function(model, routing = "expuh") {
  # if model is not a hydromad object, it is a string defining an sma
  if (!is(model, "hydromad")) {
    # no routing  for these smas
    if (model == "sacramento") {
      hydromad_model <- hydromad::hydromad(
        DATA = NULL,
        sma = model
      )
    } else {
      hydromad_model <- hydromad::hydromad(
        DATA = NULL,
        sma = model,
        routing = routing
      )
    }
  } else {
    hydromad_model <- model
  }

  # no routing parameters for these smas
  if (hydromad_model$sma != "sacramento") {
    if (hydromad_model$routing == "expuh") {
      hydromad_model <- hydromad::update(hydromad_model,
        tau_s = c(5, 100),
        tau_q = c(0, 5),
        v_s = c(0, 1)
      )
    } else if (hydromad_model$routing == "lambda") {
      hydromad_model <- hydromad::update(hydromad_model,
        tau_s = c(5, 100),
        tau_q = c(0, 5),
        v_s = c(0, 1),
        lambda = c(-1, 0)
      )
    }
  }

  out <- list()
  out$routing <- hydromad_model$routing
  out$lower <- hydromad_model %>%
    hydromad::getFreeParsRanges() %>%
    purrr::map_dbl(\(x) x[1])
  out$upper <- hydromad_model %>%
    hydromad::getFreeParsRanges() %>%
    purrr::map_dbl(\(x) x[2])

  out$nof_param <- length(out$lower)

  # DEoptim settings
  out$DEoptim$NP <- out$nof_param * 10
  if (hydromad_model$sma == "sacramento") {
    out$DEoptim$itermax <- 200
  } else {
    out$DEoptim$itermax <- 50
  }

  # malschains settings
  out$malschains$maxEvals <- 2000

  # hydroPSO settings
  out$hydroPSO$control$write2disk <- FALSE
  out$hydroPSO$control$verbose <- FALSE
  out$hydroPSO$control$npart <- 80
  out$hydroPSO$control$maxit <- 50
  out$hydroPSO$control$reltol <- 1E-10

  # define if snow module is included
  if (hydromad_model$sma == "snow") {
    out$has_snow_module <- TRUE
  } else {
    out$has_snow_module <- FALSE
  }

  return(out)
}


#' Create a hydromad Model Object
#'
#' Constructs a \code{hydromad} model using a specified soil moisture accounting (SMA) routine
#' and a set of default calibration parameters.
#'
#' @param sma A string specifying the soil moisture accounting model (e.g., `"cmd"`, `"gr4j"`, `"sacramento"`).
#' @param cal_par A list of calibration parameters, typically created using \code{\link{set_hydromad_par}}.
#'   Must contain \code{lower} and \code{upper} bounds for each parameter.
#' @param routing A string specifying the routing model (default is `"expuh"`).
#'
#' @return A \code{hydromad} model object with parameter ranges set.
#'
#' @details
#' The function initializes a \code{hydromad} model with no data and updates it with the parameter
#' bounds provided in \code{cal_par}. This is useful for preparing a model structure before calibration.
#'
#' @importFrom hydromad hydromad update
#' @importFrom purrr map2
#'
#' @examples
#' create_hydromad_model("cmd", default_cal_par[["cmd"]])
#'
#' @export
create_hydromad_model <- function(sma, cal_par, routing = "expuh") {
  pars <- purrr::map2(cal_par$lower, cal_par$upper, ~ c(.x, .y))

  output_model <- hydromad::hydromad(
    DATA = NULL,
    sma = sma,
    routing = routing
  )
  output_model <- hydromad::update(output_model, newpars = pars)
}


# Settings ---------------------------------
library(airGR)

#' Hydrofamily or Package
#'
#' @export hydro_family
hydro_family <- data.frame(
  model = c(
    "GR4J", "GR5J", "GR6J", "CemaNeigeGR4J", "CemaNeigeGR5J", "CemaNeigeGR6J", "CemaNeige",
    "TUW", "TUWsnow", "topmodel",
    "cmd", "cwi", "awbm", "bucket", "sacramento", "snow"
  ),
  family = c(rep("airGR", 7), "TUWmodel", "TUWmodel", "topmodel", rep("hydromad", 6))
)

# airGR models
airGR_models <- dplyr::filter(hydro_family, family == "airGR")$model
names(airGR_models) <- airGR_models

# hydromad_models
hydromad_models <- dplyr::filter(hydro_family, family == "hydromad")$model
names(hydromad_models) <- hydromad_models

#' Default calibration parameters
#'
#' @export default_cal_par
default_cal_par <- purrr::map(airGR_models, set_airGR_par)

default_cal_par <- c(default_cal_par, purrr::map(hydromad_models, set_hydromad_par))

## set snow snow module (without f, e and routing parameters)
default_cal_par$snowsnow$lower <- default_cal_par$snow$lower[1:7]
default_cal_par$snowsnow$upper <- default_cal_par$snow$upper[1:7]
default_cal_par$snowsnow$nof_param <- length(default_cal_par$snowsnow$lower)

## set TUW model
default_cal_par$TUW$lower <- c(SCF = 0.9, DDF = 0.0, Tr = 1.0, Ts = -3.0, Tm = -2.0, LPrat = 0.0, FC = 0.0, BETA = 0.0, k0 = 0.0, k1 = 2.0, k2 = 30.0, lsuz = 1.0, cperc = 0.0, bmax = 0.0, croute = 0.0)
default_cal_par$TUW$upper <- c(SCF = 1.5, DDF = 5.0, Tr = 3.0, Ts = 1.0, Tm = 2.0, LPrat = 1.0, FC = 600.0, BETA = 20.0, k0 = 2.0, k1 = 30.0, k2 = 250.0, lsuz = 100.0, cperc = 8.0, bmax = 30.0, croute = 50.0)
default_cal_par$TUW$nof_param <- length(default_cal_par$TUW$lower)
default_cal_par$TUW$has_snow_module <- TRUE

## set TUW model snow module
default_cal_par$TUWsnow$lower <- default_cal_par$TUW$lower[1:5]
default_cal_par$TUWsnow$upper <- default_cal_par$TUW$upper[1:5]
default_cal_par$TUWsnow$nof_param <- length(default_cal_par$TUWsnow$lower)

# DEoptim settings
default_cal_par$TUW$DEoptim$NP <- 300
default_cal_par$TUW$DEoptim$itermax <- 200

# malchains settings
default_cal_par$TUW$malschains$maxEvals <- 2000

# hydroPSO settings
default_cal_par$TUW$hydroPSO$control$write2disk <- FALSE
default_cal_par$TUW$hydroPSO$control$verbose <- FALSE
default_cal_par$TUW$hydroPSO$control$npart <- 80
default_cal_par$TUW$hydroPSO$control$maxit <- 50
default_cal_par$TUW$hydroPSO$control$reltol <- 1E-10

# from hyroPSO vignette
# npart=80,
# maxit=50,
# reltol=1E-10
# Using npart=40 should be good enough for most model applications. However, if the number of parameters of
# your model is ~10 or more, I suggest to explore the use of a larger number of particles in combination
# with a lower number of model runs (e.g., npart=80 and maxit=50).


## set topmodel
default_cal_par$topmodel$lower <- c(qs0 = 1e-04, lnTe = 2, m = 0, Sr0 = 0, Srmax = 0, td = 0, vch = 100, vr = 100, k0 = 0, CD = 0)
default_cal_par$topmodel$upper <- c(qs0 = 2.5e-04, lnTe = 3, m = 0.2, Sr0 = 0.2, Srmax = 2, td = 3, vch = 2500, vr = 2500, k0 = 10, CD = 5)
default_cal_par$topmodel$nof_param <- length(default_cal_par$topmodel$lower)
default_cal_par$topmodel$has_snow_module <- FALSE

# DEoptim settings
default_cal_par$topmodel$DEoptim$NP <- 100
default_cal_par$topmodel$DEoptim$itermax <- 50

# malchains settings
default_cal_par$topmodel$malschains$maxEvals <- 2000

# hydroPSO settings
default_cal_par$topmodel$hydroPSO$control$write2disk <- FALSE
default_cal_par$topmodel$hydroPSO$control$verbose <- FALSE
default_cal_par$topmodel$hydroPSO$control$npart <- 80
default_cal_par$topmodel$hydroPSO$control$maxit <- 50
default_cal_par$topmodel$hydroPSO$control$reltol <- 1E-10
