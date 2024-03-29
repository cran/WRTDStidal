% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kendallSeasonalTrendTest.R
\name{kendallSeasonalTrendTest}
\alias{kendallSeasonalTrendTest}
\alias{kendallSeasonalTrendTest.default}
\alias{kendallSeasonalTrendTest.data.frame}
\alias{kendallSeasonalTrendTest.formula}
\alias{kendallSeasonalTrendTest.matrix}
\title{Kendall seasonal trend test}
\usage{
kendallSeasonalTrendTest(y, ...)

\method{kendallSeasonalTrendTest}{default}(
  y,
  season,
  year,
  alternative = "two.sided",
  correct = TRUE,
  ci.slope = TRUE,
  conf.level = 0.95,
  independent.obs = TRUE,
  data.name = NULL,
  season.name = NULL,
  year.name = NULL,
  parent.of.data = NULL,
  subset.expression = NULL,
  ...
)

\method{kendallSeasonalTrendTest}{data.frame}(y, ...)

\method{kendallSeasonalTrendTest}{formula}(y, data = NULL, subset, na.action = na.pass, ...)

\method{kendallSeasonalTrendTest}{matrix}(y, ...)
}
\arguments{
\item{y}{an object containing data for the trend test. In the default method, the argument \code{y} must be numeric vector of observations. When \code{y} is a \code{data frame}, all columns must be numeric. When \code{y} is a \code{matrix}, it must be a numeric matrix. In the formula method, \code{y} must be a formula of the form \code{y ~ season + year}, where \code{y}, \code{season}, and \code{year} specify what variables to use for the these arguments in the call to \code{kendallSeasonalTrendTest.default}. Missing (\code{NA}), undefined (\code{NaN}), and infinite (\code{Inf, -Inf}) values are allowed but will be removed.}

\item{...}{methods passed to or from other methods}

\item{season}{numeric or character vector or a factor indicating the seasons in which the observations in y were taken. The length of \code{season} must equal the length of \code{y}.}

\item{year}{numeric vector indicating the years in which the observations in \code{y} were taken. The length of \code{year} must equal the length of \code{y}.}

\item{alternative}{character string indicating the kind of alternative hypothesis. The possible values are \code{"two.sided"} (tau not equal to 0; the default), \code{"less"} (tau less than 0), and \code{"greater"} (tau greater than 0).}

\item{correct}{logical scalar indicating whether to use the correction for continuity in computing the z-statistic that is based on the test statistic S'. The default value is \code{TRUE}.}

\item{ci.slope}{logical scalar indicating whether to compute a confidence interval for the slope. The default value is \code{TRUE}.}

\item{conf.level}{numeric scalar between 0 and 1 indicating the confidence level associated with the confidence interval for the slope. The default value is \code{0.95}.}

\item{independent.obs}{logical scalar indicating whether to assume the observations in y are seially independent. The default value is \code{TRUE}.}

\item{data.name}{character string indicating the name of the data used for the trend test. The default value is \code{deparse(substitute(y))}.}

\item{season.name}{character string indicating the name of the data used for the season. The default value is \code{deparse(substitute(season))}.}

\item{year.name}{character string indicating the name of the data used for the year. The default value is \code{deparse(substitute(year))}.}

\item{parent.of.data}{character string indicating the source of the data used for the trend test.}

\item{subset.expression}{character string indicating the expression used to subset the data.}

\item{data}{specifies an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{kendallTrendTest} is called.}

\item{subset}{specifies an optional vector specifying a subset of observations to be used.}

\item{na.action}{specifies a function which indicates what should happen when the data contain \code{NAs}. The default is \code{\link[stats]{na.pass}}.}
}
\value{
A list object with elements for results of the test
}
\description{
Nonparametric test for monotonic trend Within each season based on Kendall's Tau statistic
}
\details{
Perform a nonparametric test for a monotonic trend within each season based on Kendall's tau statistic, and optionally compute a confidence interval for the slope across all seasons.
}
\examples{
kendallSeasonalTrendTest(res ~ month + year, tidfitmean)
}
\references{
Hirsch, R.M., Slack, J.R., Smith, R.A. 1982. Techniques of trend analysis for monthly water quality data. Water Resources Research, 18:107-121.
 
Millard, S. P. 2013. EnvStats: An R Package for Environmental Statistics. Springer, New York.
}
