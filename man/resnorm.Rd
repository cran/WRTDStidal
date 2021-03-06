% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resnorm.R
\name{resnorm}
\alias{resnorm}
\alias{resnorm.tidal}
\alias{resnorm.tidalmean}
\title{Get salinity/flow normalized WRTDS predictions from interpolation grids}
\usage{
resnorm(dat_in, ...)

\method{resnorm}{tidal}(dat_in, trace = TRUE, ...)

\method{resnorm}{tidalmean}(dat_in, trace = TRUE, ...)
}
\arguments{
\item{dat_in}{input tidal or tidalmean object}

\item{...}{arguments passed to or from other methods}

\item{trace}{logical indicating if progress is shown in the console}
}
\value{
Appends columns to the data.frame for normalized values.  For, tidal objects, columns are named starting with the prefix `norm', e.g., `norm0.5' are the normalized values for the fit through the median.  For tidalmean objects, columns are appended for the log-transformed and back-transformed normalized values, named `norm' and `bt_norm'.
}
\description{
Get normalized model predictions from WRTDS to remove the effect of salinity/flow on the response variable.  Predicted values in the interpolation grids are averaged across dates.
}
\details{
This function is used after \code{wrtds} to normalize predicted values of the response variable from the interpolation grid for each model.  The normalized values are based on the average of all predicted estimates across the range of salinity/flow values that have occurred on the same date throughout each year.  For example, normalized values for July 2000 are the mean predicted response at that date using the observed salinity/flow values that occur in July of all years.  The normalized values allow an interpretation of trends in the response variable that are independent of changes in salinity or freshwater inputs.
}
\examples{
\dontrun{
##

# load a tidal object
data(tidobj)

# get flow-normalized values for each quantile
res <- resnorm(tidobj)

# load a tidalmean object
data(tidobjmean)

# get flow-normalized values
res <- resnorm(tidobjmean)
}
}
