## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = T, message = F, warning = F---------------------------------
library(WRTDStidal)

## ------------------------------------------------------------------------
# import data
data(chldat)

# data format
str(chldat)

## ----eval = F------------------------------------------------------------
#  # load a fitted model, quantiles
#  data(tidfit)
#  
#  # load a fitted model, mean
#  data(tidfitmean)
#  
#  # or recreate the quanitile models from chldat
#  tidfit <- modfit(chldat, tau = c(0.1, 0.5, 0.9))
#  
#  # or recreate the mean model from chldat
#  tidfitmean <- modfit(chldat, resp_type = 'mean')

## ----fig.height = 6, fig.width = 8---------------------------------------
# create a tidal object from a data frame, or use tidalmean function
tidobj <- tidal(chldat)

# plot the raw data
obsplot(tidobj)

## ------------------------------------------------------------------------
# data
head(tidobj)

# names of the attributes
names(attributes(tidobj))

# load a fitted tidal object, or use tidfitmean
data(tidfit)

# fitted data
head(tidfit)

# fitted attributes
names(attributes(tidfit))


## ----message = FALSE, cache = TRUE---------------------------------------
# get wrtds results, quantile model
mod <- modfit(chldat)

# get wrtds mean model
mod <- modfit(chldat, resp_type = 'mean')

## ----message = FALSE, eval = FALSE---------------------------------------
#  # this is equivalent to running modfit
#  # modfit is a wrapper for tidal, wrtds, respred, and resnorm functions
#  
#  # pipes from the dplyr (magrittr) package are used for simplicity
#  library(dplyr)
#  
#  # quantile model
#  mod <- tidal(chldat) %>%  # creates a tidal object
#    wrtds %>% # creates wrtds interpolation grids
#    respred %>% # get predictions from grids
#    resnorm # get normalized predictions from grids
#  
#  # mean model
#  mod <- tidalmean(chldat) %>%  # creates a tidal object
#    wrtds %>% # creates wrtds interpolation grids
#    respred %>% # get predictions from grids
#    resnorm # get normalized predictions from grids

## ----eval = FALSE--------------------------------------------------------
#  ## fit the model and get predicted/normalized chlorophyll data
#  # default median fit, quantile model
#  # grids predicted across salinity range with ten values
#  mod <- modfit(chldat)
#  
#  ## fit different quantiles and smaller interpolation grid
#  mod <- modfit(chldat, tau = c(0.2, 0.8), flo_div = 5)
#  
#  ## fit with different window widths
#  # half-window widths of one day, five years, and 0.3 salinity
#  mod <- modfit(chldat, wins = list(1, 5, 0.3))
#  
#  ## suppress console output
#  mod <- modfit(chldat, trace = FALSE)

## ----fig.height=4, fig.width=8, message=FALSE, warning=FALSE-------------
# load data from the package for the example
data(tidfit)

# plot using fitplot function
fitplot(tidfit)

# plot non-aggregated results
fitplot(tidfit, annuals = FALSE)

## ----fig.height=4, fig.width=8, message=FALSE, warning=FALSE-------------
# plot january, july as defaults
sliceplot(tidfit)

## ----fig.height=6, fig.width=8, message=FALSE, warning=FALSE-------------
# fits by month, normalized
fitmoplot(tidfit, predicted = F)

## ----fig.height=4, fig.width=8, message=FALSE, warning=FALSE-------------
seasplot(tidfit)

## ----fig.height=4, fig.width=8, message=FALSE, warning=FALSE-------------
seasyrplot(tidfitmean, predicted = F)

## ----fig.height=4, fig.width=8, message=FALSE, warning=FALSE-------------
# plot predicted, normalized results for each quantile
prdnrmplot(tidfit)

# plot as monthly values
prdnrmplot(tidfit, annuals = FALSE)

## ----fig.height=6, fig.width=8, message=FALSE, warning=FALSE-------------
# plot using defaults
# defaults to the fiftieth quantile
dynaplot(tidfit)

## ----fig.height=6, fig.width=8, message=FALSE, warning=FALSE-------------
# create a gridded plot
# defaults to the fiftieth quantile
gridplot(tidfit)
gridplot(tidfit, month = 'all')

## ----fig.height=6, fig.width=8, message=FALSE, warning=FALSE-------------
library(dplyr)
library(plotly)

dat <- attr(tidfitmean, 'fits') %>% 
  .[[1]] %>% 
  select(-date, -year, -month, -day) %>% 
  as.matrix

scene <- list(
  aspectmode = 'manual', 
  aspectratio = list(x = 0.5, y = 1, z = 0.3), 
  xaxis = list(title = 'Salinity'), 
  yaxis = list(title = 'Time'), 
  zaxis = list(title = 'log-Chl')
  )

p <- plot_ly(z = ~dat) %>% 
  add_surface(colors = rev(RColorBrewer::brewer.pal(11, 'Spectral'))) %>% 
  layout(scene = scene)
p

## ----fig.height=7, fig.width=7, message=FALSE, warning=FALSE-------------
# wt plot
wtsplot(tidfit, ref = '1995-07-01')

## ----fig.height=3, fig.width=5, message=FALSE, warning=FALSE-------------
# create a nobsplot
nobsplot(tidfit)

## ------------------------------------------------------------------------
wrtdsperf(tidfit)
# setup month, year categories for trend summaries
mobrks <- list(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9), c(10, 11, 12))
yrbrks <- c(-Inf, 1985, 1994, 2003, Inf)
molabs <- c('JFM', 'AMJ', 'JAS', 'OND')
yrlabs <- c('1974-1985', '1986-1994', '1995-2003', '2004-2012')
wrtdstrnd(tidfit, mobrks, yrbrks, molabs, yrlabs)

## ------------------------------------------------------------------------
wrtdstrnd_sk(tidfit, mobrks, yrbrks, molabs, yrlabs)

