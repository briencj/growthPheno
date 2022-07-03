\name{as.smooths.frame}
\alias{as.smooths.frame}
\title{Forms a \code{\link{smooths.frame}} from a \code{data.frame}, ensuring that the correct 
       columns are present.}
\description{Creates a \code{\link{smooths.frame}} from a \code{\link{data.frame}} by adding the 
             \code{class} \code{\link{smooths.frame}} and a set of \code{\link{attributes}} to it.}

\usage{
as.smooths.frame(data, individuals = NULL, times = NULL)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the results of smoothing the 
            data on a set of \code{individuals} over time, the data being arranged
            in long format both with respect to the times and the 
            smoothing-parameter values. It must contain the columns \code{Type}, 
            \code{TunePar}, \code{TuneVal}, \code{Tuning} and \code{Method} that 
            give the smoothing-parameter values that were used to produce each 
            smooth of the data, as well as the columns identifying the 
            \code{individuals}, the observation times of the responses and 
            the unsmoothed and smoothed responses. 
            Each response occupies a single column.}
 \item{individuals}{A \code{\link{character}} giving the name of the 
            \code{\link{factor}} that defines the subsets of the \code{data} 
            for which each subset corresponds to the \code{response} values for 
            an individual (plant/cart/plot/unit).}
 \item{times}{A \code{\link{character}} giving the name of the 
            \code{\link{numeric}}, or \code{\link{factor}} with numeric levels, that 
            contains the values of the predictor variable to be supplied to 
            \code{\link{smooth.spline}} and to be plotted on the x-axis.}
}
\value{A \code{\link{smooths.frame}}}


\author{Chris Brien}
\seealso{\code{\link{validSmoothsFrame}}, \code{\link{as.smooths.frame}}}
\examples{
dat <- read.table(header = TRUE, text = "
Type TunePar TuneVal Tuning Method       ID  DAP   PSA      sPSA
NCSS      df       4   df-4 direct 045451-C   28 57.446 51.18456
NCSS      df       4   df-4 direct 045451-C   30 89.306 87.67343
NCSS      df       7   df-7 direct 045451-C   28 57.446 57.01589
NCSS      df       7   df-7 direct 045451-C   30 89.306 87.01316
")
dat[1:7] <- lapply(dat[1:6], factor)
dat <- as.smooths.frame(dat, individuals = "ID", times = "DAP")
is.smooths.frame(dat)
validSmoothsFrame(dat)
}

\keyword{asreml}
\keyword{htest}