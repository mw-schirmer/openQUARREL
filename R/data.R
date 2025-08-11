#' Input data for the hydrological models in vignette
#'
#' Precipitation, temperature, potential evaporation and observed discharge
#'
#'
#' @format ## `input_data`
#' A data frame with 148,530 rows and 7 columns:
#' \describe{
#'   \item{HSU_ID}{The ID of the catchment used by the federal office for the environment (FOEN)}
#'   \item{DatesR}{Dates}
#'   \item{P}{Precipitation in mm/d}
#'   \item{T}{Temperature in deg Celsius}
#'   \item{E}{Potential evaporation in mm/d}
#'   \item{Qm3s}{observed discharge in m3/s}
#'   \item{Qmm}{observed discharge in mm/d}
#' }
#' @source <https://www.meteoswiss.admin.ch/climate/the-climate-of-switzerland/spatial-climate-analyses.html>
#' @source <https://www.bafu.admin.ch/bafu/en/home/topics/water/data-and-maps/water-monitoring-data/hydrological-data-service-for-watercourses-and-lakes.html>
"input_data"
