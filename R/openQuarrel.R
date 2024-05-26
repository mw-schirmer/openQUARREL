# functions  ----------------------------------------------

#' Transforms runoff
#'
#' Inverse, sqrt and no transform of runoff, log is not advised for KGE, see
#' airGR or Santos 2018
#'
#' @param Q vector, matrix, data.frame etc of runoff values
#' @param q_transfo_type a string indicating how the runoff should be
#'   transformed. Currently \code{"none"}, \code{"sqrt"}, \code{"inv"} and
#'   \code{"log"} are supported.
#' @note todo: discuss what to do with Infinity
#' @return transformed runoff in same format as input
#' @export
#' @examples
#' transfo_q(array(0:10, c(2, 5)), "log")
transfo_q <- function(Q, q_transfo_type = "none") {
  switch(q_transfo_type,
    none = Q_transfo <- Q,
    sqrt = Q_transfo <- sqrt(Q),
    inv = Q_transfo <- 1 / Q,
    log = Q_transfo <- log(Q),
    stop("Invalid transformation choice: \"", q_transfo_type,
      "\". Choose one of the following: none, sqrt, inv, log.",
      call. = FALSE
    )
  )
  return(Q_transfo)
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
#' @import hydroGOF
#'
#' @return transformed runoff in same format as input
#' @export
#' @examples
#' # for the first example the function KGE must be loaded,
#' #  for example with library(hydroGOF)
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
  out <- GOF_fun(Qsim, Qobs, na.rm = na.rm)
  return(out)
}

#' Get the family of a model
#'
#' Reverse lookup table to get the package/family of a hydrological model
#'
#' @param model_str a string indicating a hydrological model
#' @return the package or family name
#' @export
#' @examples get_family("GR4J")
get_family <- function(model_str) {
  family <- dplyr::filter(hydro_family, model == model_str)$family
  return(family)
}




#' Load meteo data
#'
#' Loads a data frame stored in an .rds file containing meteo data as required
#' by the airGR family (see
#' \url{https://odelaigue.github.io/airGR/page_1_get_started.html)}. Required
#' are columns listed under Details
#'
#' Required columns:
#' * DatesR: dates in the POSIXt or Date format
#' * P: average precipitation \[mm/day\]
#' * T: catchment average air temperature \[<U+2103>\]
#' * E: catchment average potential evapotranspiration \[mm/day\]
#' * Qmm: outlet discharge \[mm/day\]
#'
#' @param file file
#' @param tzone transfers the \code{DatesR} column to POSiXct to time zone
#'   \code{tzone} if \code{"UTC"} (default) is chosen, a Date vector is not
#'   transferred to another time zone
#' @return a data frame BasinObs
#' @export
#' @examples
#' load_meteo_data("D:/input/airGR/HSU_2044.rds")
load_meteo_data <- function(file, tzone = "UTC") {
  BasinObs <- readr::read_rds(file)

  # convert Date to POSIXct
  BasinObs$DatesR <- as.POSIXct(BasinObs$DatesR)
  attr(BasinObs$DatesR, "tzone") <- tzone

  return(BasinObs)
}


#' Create model input
#'
#' Create input structure dependent on hydrological model choice
#'
#' @param model a string indicating a hydrological model
#' @param BasinObs data frame with time series of input data from
#'   \code{\link{load_meteo_data}}
#' @param BasinInfo a list with spatial basin information
#' @return model a string specifying the hydrological model
#' @export
#' @examples
#' create_input("TUW", BasinObs, BasinInfo)
create_input <- function(model, BasinObs, BasinInfo) {
  if (model == "TUW") {
    input <- data.frame(P = BasinObs$P, T = BasinObs$T, E = BasinObs$E, BasinArea = BasinInfo$BasinArea)
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
      "\". Currently implemented are those from the airGR and TUWmodel package.",
      call. = FALSE
    )
  }
}

#' split a date vector
#'
#' Splits a date vector hardcoded into three time periods, i.e. warm up,
#' calibration and validation period
#'
#' @param date a vector of dates (e.g. Dates, Posix)
#' @param start_end_date_vec vector of length six, with entries start and end
#'   date for warm-up, start and end date for calibration and start and end date
#'   for validation, in this order.
#' @return a list with three vectors containing indices \code{ind_warm},
#'   \code{ind_cal}, \code{ind_val}
#' @export
#' @examples
#' split_data_set(
#'   seq(from = as.Date("1981-01-01"), to = as.Date("2020-12-31"), by = "days"),
#'   c("1981-01-01", "1982-12-31", "1983-01-01", "2000-12-31", "2001-01-01", "2020-12-31")
#' )
split_data_set <- function(date, start_end_date_vec) {
  ind_warm <- seq(
    which(format(date, format = "%Y-%m-%d") == start_end_date_vec[1]),
    which(format(date, format = "%Y-%m-%d") == start_end_date_vec[2])
  )

  ind_cal <- seq(
    which(format(date, format = "%Y-%m-%d") == start_end_date_vec[3]),
    which(format(date, format = "%Y-%m-%d") == start_end_date_vec[4])
  )

  ind_val <- seq(
    which(format(date, format = "%Y-%m-%d") == start_end_date_vec[5]),
    which(format(date, format = "%Y-%m-%d") == start_end_date_vec[6])
  )

  return(list(ind_warm = ind_warm, ind_cal = ind_cal, ind_val = ind_val))
}

#' Monthly indices
#'
#' Returns monthly indices for specified months within a date vector
#'
#' @param date a vector of dates (e.g. Dates, Posix)
#' @param months a vector of months in double digits strings, e.g.
#' \code{c("02",  "03")}
#' @param ind a vector of indices which can be used to subset the time dependent
#'   elements in input.
#' @return a vector with indices
#'
#' @export
#' @examples
#' find_monthly_indices(seq(from = as.Date("1981-01-01"), to = as.Date("2020-12-31"), by = "days"), c("02", "03"))
find_monthly_indices <- function(date, months, ind = seq_along(date)) {
  dates_df <- data.frame(date = date[ind]) %>%
    dplyr::mutate(ind = ind) %>%
    dplyr::mutate(month = format(as.Date(date), "%m")) %>%
    dplyr::filter(month %in% months)

  return(dates_df$ind)
}



#' Min-max normalization and vice versa
#'
#' Scales data to \[0,1\] with min and max of applied data or with specified
#' bounds
#'
#' @note todo: \code{function(x, min = min(x), max = max(x), direction = "RT")}
#'   does not work, other solutions?
#'   default values for min and max only valid if
#'   direction is also RT, how can this be required
#' @param x data to be scaled, can be a vector, array or matrix
#' @param min lower bound
#' @param max upper bound
#' @param direction string indicating from real to transformed \code{"RT"} or
#'   vice-versa \code{"TR"}. For direction \code{"TR"} the lower and upper
#'   bounds are required inputs.
#'
#'
#' @return transformed or re-transformed data of same format as \code{x}
#' @export
#'
#' @examples
#' vec <- runif(10, -5, 10)
#' norm_minmax(vec, min(vec), max(vec))
#' vecT <- norm_minmax(vec, -5, 10)
#' # re-scale to vec
#' norm_minmax(vecT, -5, 10, "TR")
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


#' Parameter transformation
#'
#' Transforms model parameters to hypercube and vice versa
#'
#' For airGR models the functions \code{\link[airGR]{TransfoParam}} of package
#' \code{\link{airGR}}  are applied.
#' Model combinations as \code{"CemaNeigeGR4J"} are allowed.
#'
#' For other models the parameter space is transformed to \[0,1\] with
#' \code{\link{norm_minmax}}
#'
#' @note todo: 1) CemaNeigeHyst is not implemented yet
#'   2) how do the links work without airGR installed?
#'   3) implement functionality that airGR functions are not used at all
#' @param param vector of model parameters
#' @param direction string indicating from real to transformed \code{"RT"} or
#'   vice-versa \code{"TR"}.
#' @param model a string specifying the hydrological model
#' @seealso \code{\link[airGR]{TransfoParam}}
#'
#'
#' @return a vector with transformed parameters
#' @export
#'
#' @examples
#' # scale a parameter set for model "TUW" to [0,1] and back
#' param <- c(1, 2, 3, -1, 1, 1, 200, 10, 1, 15, 100, 50, 2, 15, 50)
#' scaled <- transfo_param(param, "RT", "TUW")
#' rescaled <- transfo_param(scaled, "TR", "TUW")
#' # scale a parameter set for "CemaNeigeGR4J"  to [-9.99,9.99] and back
#' param <- c(1000, 2, 250, 7, .2, 109.0365)
#' scaled <- transfo_param(param, "RT", "CemaNeigeGR4J")
#' rescaled <- transfo_param(scaled, "TR", "CemaNeigeGR4J")
transfo_param <- function(param, direction, model, cal_parameter = default_cal_par) {

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
    param_trans <- purrr::pmap_dbl(list(
      param, cal_parameter[[model]][["lower"]],
      cal_parameter[[model]][["upper"]], direction
    ), norm_minmax)
  }

  return(param_trans)
}




#' Call calibration function
#'
#' Call several optimization functions for model calibration which are available
#' as R packages
#'
#' @param cal_fn a string specifying optimization algorithm to be used for the
#'   model calibration Currently \code{"DEoptim"}, \code{"hydroPSO"} and
#'   \code{"malschains"} are supported.
#' @param hydro_data a list or data frame containing observed runoff loaded
#'   with \code{\link{load_meteo_data}}
#' @param split_indices indices indicating the
#'   warm-up and calibration period, i.e. the output of
#'   \code{\link{split_data_set}}
#' @param model a string indicating the hydrological model
#' @param input model input returned by \code{\link{create_input}}
#' @param error_crit a string of a function indicating an error criterion, e.g.
#'   \code{\link[hydroGOF]{KGE}}, or others from the \code{\link{hydroGOF}}
#'   which are usable by \code{\link{calc_hydroGOF}}
#'   functions
#' @param cal_maximize  a logical indicating if function should be
#'   maximized
#' @param cal_q_transfo a string indicating how runoff should be transformed
#'   with \code{\link{transfo_q}}
#' @param do_transfo_param logical indicating
#'   whether calibration is done with real or transformed parameters
#' @param cal_par a list of calibration settings dependent on the calibration
#'   function. Package default values \code{default_cal_par} can be overwritten
#'   by users.
#' @export
#' @examples
#' cal_output <- call_cal_fn(
#'   cal_fn, hydro_data, split_indices, model, input,
#'   error_crit, cal_maximize, cal_q_transfo, do_transfo_param, cal_par
#' )
call_cal_fn <- function(cal_fn, hydro_data, split_indices, model, input,
                        error_crit, cal_maximize, cal_q_transfo, do_transfo_param,
                        cal_par = default_cal_par) {

  # determine upper and lower boundaries for parameter set
  lower <- cal_par[[model]][["lower"]]
  upper <- cal_par[[model]][["upper"]]

  # transform to hypercube if chosen
  if (do_transfo_param) {
    lower <- transfo_param(lower, "RT", model)
    upper <- transfo_param(upper, "RT", model)
  }

  if (cal_fn == "DEoptim") {

    # differential evolution
    cal_results <- DEoptim::DEoptim(
      fn = optim_fn, hydro_data = hydro_data, split_indices = split_indices,
      model = model, input = input, error_crit = error_crit, cal_maximize = cal_maximize,
      cal_q_transfo = cal_q_transfo, do_transfo_param = do_transfo_param,
      lower = lower, upper = upper,
      control = DEoptim::DEoptim.control(
        NP = cal_par[[model]][["DEoptim"]][["NP"]],
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
      model = model, input = input, error_crit = error_crit, cal_maximize = cal_maximize,
      cal_q_transfo = cal_q_transfo, do_transfo_param = do_transfo_param,
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
        error_crit, cal_maximize, cal_q_transfo, do_transfo_param
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
      "\". Choose one of the following: DEoptim, hydroPSO and malschains.",
      call. = FALSE
    )
  }

  # rescale parameters back to real if needed
  if (do_transfo_param) {
    model_param <- transfo_param(model_param, "TR", model)
  }

  #  create output list
  return(list(
    model_param = model_param,
    error_crit_val = error_crit_val,
    cal_results = cal_results
  ))
}


#' Function to be optimized during calibration
#'
#' Function that takes global parameters (see Details) and performs
#' optimization to find best parameter set
#'
#' Simulates a hydrological model during warm-up and calibration
#' period, and calculates an error criterion \code{error_crit} on the
#' calibration period only.
#'
#' Replace NA in error criterion, which is generated when simulated runoff is
#' * solely NAs: set a particularly bad value, i.e. +/- 1e10
#'   (dependent on minimization or maximization)
#' * solely 0s for \code{\link[hydroGOF]{KGE}}: set to asymptotic value,
#'   i.e. 1 - sqrt(3)
#' * for all other cases it throws an error
#'
#' @param ParamOptim vector of model parameters to be optimized
#' @param hydro_data a list or data frame containing observed runoff loaded
#'   with \code{\link{load_meteo_data}}
#' @param split_indices indices indicating the
#'   warm-up and calibration period, i.e. the output of
#'   \code{\link{split_data_set}}
#' @param model a string indicating the hydrological model
#' @param input model input returned by \code{\link{create_input}}
#' @param error_crit a string of a function indicating an error criterion, e.g.
#'   \code{\link[hydroGOF]{KGE}}, or others from the \code{\link{hydroGOF}}
#'   which are usable by \code{\link{calc_hydroGOF}}
#'   functions
#' @param cal_maximize  a logical indicating if function should be
#'   maximized
#' @param cal_q_transfo a string indicating how runoff should be transformed
#'   with \code{\link{transfo_q}}
#' @param do_transfo_param logical indicating
#'   whether calibration is done with real or transformed parameters
#' @note: Qsim is converted to numeric to similar to Qobs, this is important for
#' TUW for example, which returns a matrix 1 x time, for the lumped case.
#' todo:
#' 1) include a spatial explicit version
#' 2) exclude warning for other NA cases
#'   in error criterion as described above
#' 3) Consider to link KGE etc with cal_maximize
#' @seealso \code{\link{calibrate_model}}, \code{\link{calibrate_model}},
#'   \code{\link{create_input}}, \code{\link{load_meteo_data}},
#'   \code{\link{calc_hydroGOF}}
#'
#' @return error criterion to be optimized
#' @export
#'
#' @examples
#' # differential evolution
#' cal_results <- DEoptim::DEoptim(
#'   fn = optim_fn, hydro_data = hydro_data, split_indices = split_indices,
#'   model = model, input = input, error_crit = error_crit, cal_maximize = cal_maximize,
#'   cal_q_transfo = cal_q_transfo, do_transfo_param = do_transfo_param,
#'   lower = lower, upper = upper,
#'   control = DEoptim::DEoptim.control(
#'     NP = cal_par[[model]][["DEoptim"]][["NP"]],
#'     itermax = cal_par[[model]][["DEoptim"]][["itermax"]],
#'     trace = 10
#'   )
#' )
optim_fn <- function(ParamOptim, hydro_data, split_indices, model, input,
                     error_crit, cal_maximize, cal_q_transfo, do_transfo_param) {

  # transform Qobs if needed and subset only calibration period
  Qobs <- transfo_q(hydro_data$BasinObs$Qmm[split_indices$ind_cal], cal_q_transfo)

  # re transform parameters to real values if needed
  if (do_transfo_param) {
    ParamOptim <- transfo_param(ParamOptim, "TR", model)
  }

  # simulate for warm up and calibration
  ind_warm_cal <- c(split_indices$ind_warm, split_indices$ind_cal)

  # transform Qsim and take only values from calibration period
  Qsim <- transfo_q(suppressMessages(simulate_model(
    model, ParamOptim, input, ind_warm_cal
  )$Qsim[-seq_along(split_indices$ind_warm)]), cal_q_transfo)

  # deal with solely NA from Qsim
  if (all(is.na(Qsim))) {
    error_crit_val <- ifelse(cal_maximize, -1e10, 1e10)
  } else {

    # Qsim  converted better to numeric
    error_crit_val <- calc_hydroGOF(error_crit, as.numeric(Qsim), Qobs)

    # set to asymptotic value, i.e. 1 - sqrt(3)
    if (is.na(error_crit_val)) {
      if (error_crit != "KGE") stop("optim_fn: NA in error criterion")
      error_crit_val <- 1 - sqrt(3)
      if (any(Qsim > 0)) {
        warning("optim_fn: something else causes NAs as well in error_crit_val")
      }
    }
  }

  ifelse(cal_maximize, return(-error_crit_val), return(error_crit_val))
}


#' Model calibration
#'
#' Calibrates a hydrological model
#'
#' @param model a string specifying the hydrological model, currently
#'   implemented are airGR and TUWmodel package models
#' @param input the output of the \code{create_input} function, dependent on the
#'   model choice, in general containing information about date, precipitation
#'   air temperature and potential evapotranspiration and spatial information as
#'   area of the catchment or hypsometric curves
#' @param BasinObs data frame with time series of input data from
#'   \code{\link{load_meteo_data}}
#' @param split_indices a list of indices from \code{\link{split_data_set}}
#'   containing the elements \code{ind_cal} and \code{ind_warm}
#' @param error_crit_transfo string specifying the error criterion for
#'   calibration and the runoff transformation separated by a \code{"_"}
#'   supported are validation criteria from the \code{\link{hydroGOF}} package
#'   usable by the \code{\link{calc_hydroGOF}} function, for supported runoff
#'   transformations please refer to \code{\link{transfo_q}}
#' @param cal_maximize a logical indicating if the calibration
#'   error criterion should be maximized (or minimized)
#' @param cal_fn a string specifying optimization algorithm to be used
#'   for the model calibration Currently \code{\link[airGR]{Calibration_Michel}}
#'   is implemented for \code{\link{airGR}} models, and \code{\link{DEoptim}} ,
#'   \code{\link{hydroPSO}} and \code{\link{malschains}} are supported for all
#'   models.
#' @param do_transfo_param logical indicating if parameter transformation to a
#'   hypercube should be applied during calibration
#' @param cal_par a list of calibration settings dependent on the calibration
#'   function. Package default values \code{default_cal_par} can be overwritten
#'   by users.
#' @note Calibration_Michel does not call \code{\link{optim_fn}}, this is why
#'  the large list of otherwise global variables is still required as input
#'  todo: 1) check if some input is really needed for Calibration_Michel as
#'    e.g. BasinObs
#'  2) think about how to easily access and change cal_par values for
#'   end user
#' @import airGR
#' @seealso \code{\link{call_cal_fn}}, \code{\link{optim_fn}}
#'
#' @return a list of calibration results with elements called
#' \itemize{
#'   \item \code{model_param} a vector of calibrated model parameters,
#'   \item \code{error_crit_transfo} the error criterion and transformation used,
#'   \item \code{error_crit_val} the value of the error criterion,
#'   \item \code{more info} calibration function and model specific output
#'      information.
#' }
#'
#' @export
#'
#' @examples
#' calibration_results <- calibrate_model(
#'   hydro_data, split_indices, model, input, error_crit_transfo, cal_maximize,
#'   cal_fn, do_transfo_param, cal_par
#' )
calibrate_model <- function(hydro_data, split_indices, model, input,
                            error_crit_transfo = "KGE_none", cal_maximize = TRUE,
                            cal_fn = "DEoptim", do_transfo_param = FALSE,
                            cal_par = default_cal_par) {

  # split error_crit entry
  # todo can the dummy variable be avoided?
  error_crit_transfo_split <- stringr::str_split(error_crit_transfo, "_", simplify = TRUE)

  if (length(error_crit_transfo_split) != 2) stop("error_crit_transfo not correctly specified")

  # error crit needs to be global for the optim_fn
  error_crit <- error_crit_transfo_split[1]
  cal_q_transfo <- error_crit_transfo_split[2]


  if (cal_fn %in% c("DEoptim", "hydroPSO", "malschains")) {

    # call an R package based calibration function
    cal_output <- call_cal_fn(
      cal_fn, hydro_data, split_indices, model, input,
      error_crit, cal_maximize, cal_q_transfo, do_transfo_param, cal_par
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

      RunOptions <- airGR::CreateRunOptions(
        FUN_MOD = paste0("RunModel_", model),
        InputsModel = input,
        IndPeriod_Run = split_indices$ind_cal,
        IndPeriod_WarmUp = split_indices$ind_warm,
        Outputs_Sim = "Qsim"
      )

      # if "none" in cal_q_transfo, it needs to be "" here
      InputsCrit <- airGR::CreateInputsCrit(
        FUN_CRIT = paste0("ErrorCrit_", error_crit),
        InputsModel = input,
        RunOptions = RunOptions,
        VarObs = "Q",
        Obs = hydro_data$BasinObs$Qmm[split_indices$ind_cal],
        transfo = stringr::str_replace(cal_q_transfo, "none", "")
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
      "\". Choose one of the following: Calibration_Michel, DEoptim, hydroPSO and malschains.",
      call. = FALSE
    )
  }

  # define what to output
  calibration_results <- list(
    model_param = model_param,
    error_crit_transfo = error_crit_transfo,
    error_crit_val = error_crit_val,
    more_info = more_info
  )

  return(calibration_results)
}


#' Simulate a hydrological model
#'
#' Takes a model choice, model parameters, model input and indices indicating a
#' time subset of the model input which should be used and simulates discharge
#' values and model specific output
#'
#' @param model a string specifying the hydrological model, currently
#'   implemented are airGR and TUWmodel package models
#' @param model_param a vector of model parameters specific for each model
#'   choice
#' @param input the output of the \code{create_input} function, dependent on the
#'   model choice, in general containing information about date, precipitation
#'   air temperature and potential evapotranspiration and spatial information as
#'   area of the catchment or hypsometric curves
#' @param ind a vector of indices which can be used to subset the time dependent
#'   elements in input. The default is taking all indices from the first element
#'   of input
#' @return a list with the elements \code{Qsim}, the simulated runoff, and
#'   \code{more_info} with model specific output
#' @export
#' @examples
#' simulation_results <- simulate_model("TUWmodel", calibration_results$model_param, input, ind = split_indices$ind_cal)
simulate_model <- function(model, model_param, input, ind = seq_along(input[[1]])) {
  if (model == "TUW") {

    # todo consider not rounding
    output_model <- TUWmodel::TUWmodel(
      prec = input$P[ind], airt = input$T[ind], ep = input$E[ind],
      area = input$BasinArea, param = signif(model_param, 2)
    )

    # get simulated runoff
    Qsim <- as.numeric(output_model$q)
  } else if (get_family(model) == "airGR") {

    ## preparation of the RunOptions object for the validation period
    RunOptions <- airGR::CreateRunOptions(
      FUN_MOD = paste0("RunModel_", model),
      InputsModel = input, IndPeriod_Run = ind
    )

    model_fn <- get(paste0("RunModel_", model))
    output_model <- model_fn(InputsModel = input, RunOptions = RunOptions, Param = model_param)

    # get simulated runoff
    Qsim <- output_model$Qsim
  } else {
    stop(sprintf("the model choice %s is not implemented yet", model))
  }

  # create output list
  simulation_results <- list(
    Qsim = Qsim,
    more_info = list(output_model = output_model)
  )
}

#' Validate model
#'
#' Calculates validation measures for different transformation types
#'
#' @param Qsim vector with simulated runoff
#' @param Qobs vector with observed runoff
#' @param val_crit_transfo a vector of strings specifying validation criteria and a
#'   runoff transformation separated by a \code{"_"}. Supported are validation
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
#'   c("KGE_log", "NSE_inv", "VE_none", "pbias_none")
#' )
validate_model <- function(Qsim, Qobs, val_crit_transfo = "KGE_none") {
  val_crit_transfo_split <- stringr::str_split(val_crit_transfo, "_", simplify = TRUE)

  if (ncol(val_crit_transfo_split) != 2) stop("val_crit_transfo not correctly specified")

  val_crit <- val_crit_transfo_split[, 1]
  val_q_transfo <- val_crit_transfo_split[, 2]
  names(val_crit) <- val_crit_transfo

  # loop over validation criteria and q transformations
  validation_results <- purrr::map2_df(
    val_crit, val_q_transfo,
    ~ calc_hydroGOF(.x, transfo_q(Qsim, .y), transfo_q(Qobs, .y))
  ) %>%
    # a long data_frame
    tidyr::pivot_longer(cols = everything(), values_to = "value") %>%
    # separate the name columns in crit and transfo
    tidyr::separate(name, c("crit", "transfo"), "_")

  return(validation_results)
}


#' Calculates validation results for a subset period
#'
#' Helper function which calls \code{\link{validate_model}} and subsets \code{Qsim}
#' and \code{Qobs} with \code{ind}, and returns a new data frame with an
#' additional column called \code{col_name} with entries \code{period_name}
#'
#' @param ind indices used for subsetting \code{Qsim} and \code{Qobs}
#' @param col_name additional column name in returned data frame (default is
#'   period)
#' @param period_name entries in column \code{col_name} naming the subset period
#' @param Qsim vector with simulated runoff
#' @param Qobs vector with observed runoff
#' @param val_crit_transfo a vector of strings specifying validation criteria and a
#'   runoff transformation separated by a \code{"_"}. supported are validation
#'   criteria from the \code{\link{hydroGOF}} package usable by the
#'   \code{\link{calc_hydroGOF}} function, for supported runoff transformations
#'   please refer to \code{\link{transfo_q}}
#'
#' @seealso \code{\link{validate_model}}, \code{\link{calc_hydroGOF}},
#'   \code{\link{transfo_q}}
#'
#' @return a data frame as \code{\link{validate_model}} but with an additional
#'   columns naming the period
#' @export
#'
#' @examples
#' calc_validation_results(1:5, "summer", "season", 1:10, seq(0, 9))
calc_validation_results <- function(ind, period_name, col_name = "period",
                                    Qsim, Qobs, val_crit_transfo = "KGE_none") {
  validation_results <- validate_model(Qsim[ind], Qobs[ind], val_crit_transfo) %>%
    dplyr::mutate("{col_name}" := period_name)

  return(validation_results)
}

#' Calculates subseasonal validation results
#'
#' Within a subset of hydrological input data \code{hydro_data} subset with
#' \code{ind} it calculates performances metrics for monthly defined periods by
#' calling \code{\link{calc_validation_results}} with \code{col_name = "season"}
#' and the period names from the names in \code{val_subseason}.
#' Returns a new data frame with an
#' additional column called \code{col_name} with entries \code{period_name}
#'
#' @param val_subseason a list with named arrays of two digits describing
#' months used to calculate subseasonal validation metrics
#' @param date date a vector of dates (e.g. Dates, Posix)
#'   with \code{\link{load_meteo_data}}
#' @param ind indices used for subsetting \code{hydro_data}
#' @param col_name additional column name in returned data frame (default is
#'   period)
#' @param period_name entries in column \code{col_name} naming the subset period
#' @param Qsim vector with simulated runoff
#' @param Qobs vector with observed runoff
#' @param val_crit_transfo a vector of strings specifying validation criteria and a
#'   runoff transformation separated by a \code{"_"}. supported are validation
#'   criteria from the \code{\link{hydroGOF}} package usable by the
#'   \code{\link{calc_hydroGOF}} function, for supported runoff transformations
#'   please refer to \code{\link{transfo_q}}
#' @seealso \code{\link{calc_validation_results}}, \code{\link{validate_model}},
#'  \code{\link{calc_hydroGOF}}, \code{\link{transfo_q}}
#'
#' @return a data frame as \code{\link{validate_model}} but with two additional
#'   columns naming the period and the season
#' @export
#'
#' @examples
#' perf_cal <- calc_subseasonal_validation_results(
#'   val_subseason = list(
#'     spring = c("02", "03", "04", "05"),
#'     summer = c("06", "07", "08", "09"),
#'     hydro_data$BasinObs$DatesR,
#'     split_indices$ind_cal, "calibration",
#'     col_name = "period",
#'     simulation_results$Qsim, Qobs,
#'     val_crit_transfo =
#'       c(
#'         "KGE_none", "NSE_none", "VE_none", "pbias_none",
#'         "KGE_inv", "NSE_inv",
#'         "KGE_sqrt", "NSE_sqrt"
#'       )
#'   )
#' )
calc_subseasonal_validation_results <- function(val_subseason, dates, ind, period_name, col_name = "period",
                                                Qsim, Qobs, val_crit_transfo = "KGE_none") {
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


#' Writes ascii results
#'
#' Writes an ascii overview of parameters and validation results
#'
#' @param file filename
#' @param calibration_results a list containing the calibration results from
#'   \code{\link{calibrate_model}}. Only a vector with calibrated model
#'   parameters are written out
#' @param validation_results a data frame containing the validation results from
#'   \code{\link{validate_model}}.
#' @param equally_spaced (optional) a logical indicating if a equally spaced
#'   output should be written out
#'
#' @return a logical if the file was successfully written
#' @export
#'
#' @examples
#' write_ascii("results.txt", calibration_results, validation_results)
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

# save cal_val_plot
save_cal_val_plot_old <- function(file, data, Qsim, split_indices) {
  pdf(file = file)

  par(mfrow = c(2, 1))
  plot(data$DatesR[split_indices$ind_cal], data$Qmm[split_indices$ind_cal],
    type = "l", xlab = "", ylab = "Discharge [mm/d]",
    main = "Calibration", col = "blue"
  )
  lines(data$DatesR[split_indices$ind_cal], Qsim[split_indices$ind_cal], col = "red")
  mtext(paste("NSE =", signif(NSE(Qsim[split_indices$ind_cal], data$Qmm[split_indices$ind_cal]), 3)), 3, 0, adj = 0.95, cex = 1)
  mtext(paste("KGE =", signif(KGE(Qsim[split_indices$ind_cal], data$Qmm[split_indices$ind_cal]), 3)), 3, 0, adj = 0.5, cex = 1)
  mtext(paste("dV =", signif(pbias(Qsim[split_indices$ind_cal], data$Qmm[split_indices$ind_cal]), 2), "%"), 3, 0, adj = 0.05, cex = 1)

  plot(data$DatesR[split_indices$ind_val], data$Qmm[split_indices$ind_val],
    type = "l", xlab = "", ylab = "Discharge [mm/d]",
    main = "Validation", col = "blue"
  )
  lines(data$DatesR[split_indices$ind_val], Qsim[split_indices$ind_val], col = "red")
  mtext(paste("NSE =", signif(NSE(Qsim[split_indices$ind_val], data$Qmm[split_indices$ind_val]), 3)), 3, 0, adj = 0.95, cex = 1)
  mtext(paste("KGE =", signif(KGE(Qsim[split_indices$ind_val], data$Qmm[split_indices$ind_val]), 3)), 3, 0, adj = 0.5, cex = 1)
  mtext(paste("dV =", signif(pbias(Qsim[split_indices$ind_val], data$Qmm[split_indices$ind_val]), 2), "%"), 3, 0, adj = 0.05, cex = 1)

  dev.off()
}


#' Save a calibration validation plot
#'
#' Create and save a plot containing over the calibration and validation period
#'
#' Plots observed and simulated runoff and displays certain error metrics
#'
#' @param file the filename of the saved plot
#' @param BasinObs data frame with time series of input data from
#'   \code{\link{load_meteo_data}}
#' @param Qsim vector with simulated runoff
#' @param split_indices a list of indices from \code{\link{split_data_set}}
#'   containing the elements \code{ind_cal} and \code{ind_val}
#' @note todo: include a seasonal rolling mean
#'
#' @return a logical if the file was successfully written
#' @export
#'
#' @examples
#' save_cal_val_plot("cal_val.pdf", BasinObs, split_indices)
save_cal_val_plot <- function(file, BasinObs, Qsim, split_indices) {
  Qobs <- BasinObs$Qmm
  plot_df <- data.frame(dates = BasinObs$DatesR, Qsim = Qsim, Qobs = Qobs)

  ind_list <- list()
  ind_list$Calibration <- split_indices$ind_cal
  ind_list$Validation <- split_indices$ind_val

  val_crit_transfo <- c("pbias_none", "KGE_none", "NSE_none", "KGE_inv", "NSE_inv")
  val_crit_label <- c("%", "", "", "", "")

  validation_results <- purrr::imap(
    ind_list, calc_validation_results,
    "period", Qsim, Qobs, val_crit_transfo
  ) %>%
    purrr::map_df(dplyr::bind_rows) %>%
    # rename none and pbias
    dplyr::mutate(transfo = replace(transfo, transfo == "none", "")) %>%
    dplyr::mutate(crit = replace(crit, crit == "pbias", "dV"))

  adj_val <- seq(0, 1, 1 / (length(val_crit_transfo) - 1))

  plot_single_Qsim_Qobs <- function(plot_df, validation_results, adj_val) {
    plot(plot_df$dates, plot_df$Qobs,
      type = "l", xlab = "", ylab = "Discharge [mm/d]",
      main = validation_results[1, 4], col = "blue"
    )
    lines(plot_df$dates, plot_df$Qsim, col = "red")
    for (i in 1:nrow(validation_results)) {
      mtext(paste0(
        validation_results[i, 1], validation_results[i, 2], " = ",
        validation_results[i, 3], val_crit_label[i]
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


#' save an airGR overview plot to pdf
#'
#' use the \code{\link{airGR}} \code{\link[airGR]{plot}} function
#'
#' @param file the filename of the saved plot
#' @param simulation_results output from \code{\link{simulate_model}}, i.e. the
#'   output object of the \code{\link{airGR}} function
#'   \code{\link[airGR]{RunModel}}
#' @param Qobs a vector of observed runoff
#'
#' @return a logical if file was written successfully
#' @export
#'
#' @examples
#' save_airGR_plot("airGR_plot.pdf", simulation_results, Qobs)
save_airGR_plot <- function(file, simulation_results, Qobs) {
  if (get_family(model) == "airGR") {
    pdf(file = file)
    plot(simulation_results$more_info$output_model, Qobs = Qobs)
    dev.off()
  }
}


# Settings ---------------------------------

#' Hydrofamily or Package
#'
#' @export hydro_family
hydro_family <- data.frame(
  model = c("GR4J", "GR5J", "GR6J", "CemaNeigeGR4J", "CemaNeigeGR5J", "CemaNeigeGR6J", "TUW"),
  family = c(rep("airGR", 6), "TUW")
)


#' Create default calibration settings for the airGR family
#'
#' @param model a model string
#' @return a list of calibration settings
#' @export
#'
#' @examples
#' set_airGR_par("GR4J")
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
    nof_param <- length(out$lower)

    # DEoptim
    out$DEoptim$NP <- nof_param * 10
    out$DEoptim$itermax <- 200

    # malschains
    out$malschains$maxEvals <- 2000

    # hydroPSO
    out$hydroPSO$control$write2disk <- FALSE
    out$hydroPSO$control$verbose <- FALSE

    return(out)
  } else {
    stop("not an airGR model")
  }
}

airGR_models <- dplyr::filter(hydro_family, family == "airGR")$model
names(airGR_models) <- airGR_models


#' Default calibration parameters
#'
#' @export default_cal_par
default_cal_par <- purrr::map(airGR_models, set_airGR_par)

## set TUW model
default_cal_par$TUW$lower <- c(0.9, 0.0, 1.0, -3.0, -2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 30.0, 1.0, 0.0, 0.0, 0.0)
default_cal_par$TUW$upper <- c(1.5, 5.0, 3.0, 1.0, 2.0, 1.0, 600.0, 20.0, 2.0, 30.0, 250.0, 100.0, 8.0, 30.0, 50.0)

# DEoptim settings
default_cal_par$TUW$DEoptim$NP <- 300
default_cal_par$TUW$DEoptim$itermax <- 1000

# malchanins settings
default_cal_par$TUW$malschains$maxEvals <- 2000

# hydroPSO settings
default_cal_par$TUW$hydroPSO$control$write2disk <- FALSE
default_cal_par$TUW$hydroPSO$control$verbose <- FALSE
