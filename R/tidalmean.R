######
#' Create a tidalmean class object
#' 
#' Prepare water quality data for weighted regression for the mean response by creating a tidalmean class object
#' 
#' @param dat_in Input data frame for a water quality time series with four columns for date (Y-m-d format), response variable, salinity/flow, and detection limit for left-censored data
#' @param ind four element numeric vector indicating column positions of date, response variable, salinity/flow, and detection limit of input data frame
#' @param reslab character string or expression for labelling the response variable in plots, defaults to log-chlorophyll in ug/L
#' @param flolab character string or expression for labelling the flow variable in plots, defaults to Salinity
#' @param reslog logical indicating if input response variable is already in log-space, default \code{TRUE}
#' @param rm_miss logical indicating if missing observations in the input data are removed
#' @param ... arguments passed from other methods
#' 
#' @return A tidalmean object as a data frame and attributes.  The data frame has columns ordered as date, response variable, salinity/flow (rescaled to 0, 1 range), detection limit, logical for detection limit, day number, month, year, and decimal time.  The attributes are as follows:
#' \describe{
#'  \item{\code{names}}{Column names of the data frame}
#'  \item{\code{row.names}}{Row names of the data frame}
#'  \item{\code{class}}{Class of the object}
#'  \item{\code{half_wins}}{List of numeric values used for half-window widths for model fitting, in the same order as the wt_vars argument passed to \code{\link{getwts}}. Initially will be \code{NULL} if \code{\link{wrtds}} has not been used.}
#'  \item{\code{fits}}{List with a single element with fits for the WRTDS mean interpolation grid.  Initially will be NULL if \code{wrtds} has not been used.}
#'  \item{\code{predonobs}}{A \code{data.frame} of predictions using the observed data that were used to fit the model.  This is required for \code{wrtdsperf} if a novel dataset is used for predictions after fitting the model. Initially will be NULL if \code{respred} has not been used.}
#'  \item{\code{bt_fits}}{List with a single element with back-transformed fits for the WRTDS mean interpolation grid.  Initially will be NULL if \code{wrtds} has not been used.}
#'  \item{\code{flo_grd}}{Numeric vector of salinity/flow values that was used for the interpolation grids}
#'  \item{\code{floobs_rng}}{Two element vector indicating the salinity/flow range of the observed data}
#'  \item{\code{nobs}}{List with one matrix showing the number of weights greater than zero for each date and salinity/flow combination used to create the fit matrices in \code{fits}.  Initially will be \code{NULL} if \code{\link{wrtds}} has not been used.}
#'  \item{\code{reslab}}{expression or character string for response variable label in plots}
#'  \item{\code{flolab}}{expression or character string for flow variable label in plots}
#' }
#'
#' @details
#' This function is a simple wrapper to \code{\link[base]{structure}} that is used to create a tidalmean object for use with weighted regression in tidal waters, specifically to model the mean response as compared to a conditional quantile. Input data should be a four-column \code{\link[base]{data.frame}} with date, response variable, salinity/flow data, and detection limit for each observation of the response.  The response data are assumed to be log-transformed, otherwise use \code{reslog = FALSE}.  Salinity data can be provided as fraction of freshwater or as parts per thousand.  The limit column can be entered as a sufficiently small number if all values are above the detection limit or no limit exists.  The current implementation of weighted regression for tidal waters only handles left-censored data.  Missing observations are also removed.  
#'  
#' The tidalmean object structure is almost identical to the tidal object, with the exception of an additional attribute for the back-transformed interpolation grid.  This is included to account for retransformation bias of log-transformed variables associated with mean models.
#' 
#' @export
#' 
#' @examples
#' ## raw data
#' 
#' data(chldat)
#' 
#' ## format
#' chldat <- tidalmean(chldat)
#' 
tidalmean <- function(dat_in, ind = c(1, 2, 3, 4), reslab = NULL, flolab = NULL, reslog = TRUE, rm_miss = FALSE, ...){
  
  # sanity checks
  if(!any(c('Date', 'POSIXct', 'POSIXlt') %in% class(dat_in[, ind[1]])))
    stop('Class for time column must be Date, POSIXct, or POSIXlt')
  if(ncol(dat_in) > 4)
    stop('Only four input columns are allowed')
  
  # get relevant columns and set names
  dat_in <- dat_in[, ind, drop = F]
  names(dat_in) <- c('date', 'res', 'flo', 'lim')

  # columns as numeric
  dat_in[, 2:4] <- apply(dat_in[, 2:4], 2, function(x) as.numeric(as.character(x)))
  
  # retain only complete data
  if(nrow(na.omit(dat_in)) != nrow(dat_in) & rm_miss){
    warning('Missing observations removed from original dataset')
    dat_in <- na.omit(dat_in)
  }

  # log transform res if T
  if(!reslog) dat_in$res <- log(dat_in$res)
  
  # rescale salinity/flow to 0 - 1, save rng
  floobs_rng <- range(dat_in$flo, na.rm = TRUE)
  dat_in$flo <- with(dat_in, (flo - floobs_rng[1])/diff(floobs_rng))
  
  # TF column of limits for surv regression
  dat_in$not_cens <- with(dat_in, res > lim)
  
  # get decimal time
  dect <- dec_time(dat_in$date)
  dat_in$day_num <- dect$day_num
  dat_in$month <- dect$month
  dat_in$year <- dect$year
  dat_in$dec_time <- dect$dec_time
  
  # organize by date
  dat_in <- dat_in[order(dat_in$date), ]
  
  # plot labels
  if(is.null(reslab))
    reslab <- chllab(reslog)
  if(is.null(flolab))
    flolab <- 'Salinity'
  
  # create class, with multiple attributes
  tidalmean <- structure(
    .Data = dat_in, 
    class = c('tidalmean', 'data.frame'),
    half_wins = NULL, 
    fits = NULL, 
    predonobs = NULL,
    bt_fits = NULL,
    flo_grd = NULL,
    floobs_rng = floobs_rng,
    nobs = NULL, 
    reslab = reslab,
    flolab = flolab
    )
  
  return(tidalmean)
  
}
