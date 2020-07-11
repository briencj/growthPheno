\name{GrowthRates}
\alias{GrowthRates}
\alias{AGRdiff}
\alias{PGR}
\alias{RGRdiff}
\title{Calculates growth rates (AGR, PGR, RGRdiff) between pairs of values in a vector}
\description{Calculates either the Absolute Growth Rate (AGR), Proportionate Growth 
             Rate (PGR) or Relative Growth Rate (RGR) between pairs of time points, 
             the second of which is \code{lag} positions before the first. 
             in \code{x}.}
\usage{AGRdiff(x, time.diffs, lag=1)
PGR(x, time.diffs, lag=1)
RGRdiff(x, time.diffs, lag=1)}
\arguments{
 \item{x}{A \code{\link{numeric}} from which the growth rates are to 
          be calculated.}
 \item{time.diffs}{a \code{numeric} giving the time differences between successive 
       values in \code{x}.}
 \item{lag}{A \code{integer} specifying, for the second value in the pair to be operated on, 
            the number positions it is ahead of the current value.}
}
\value{A \code{\link{numeric}} containing the growth rates which is the same length as \code{x} 
       and in which the first \code{lag} values \code{NA}.}
\details{
The AGRdiff is calculated as the difference between a pair of values divided by the \code{time.diffs}. 
The PGR is calculated as the ratio of a value to a second value which is \code{lag} values 
ahead of the first in \code{x} and the ratio raised to the 
power of the reciprocal of \code{time.diffs}. 
The RGRdiff is calculated as the \code{log} of the PGR and so is equal to the difference between 
the logarithms of a pair of values divided by the \code{time.diffs}.
The differences and ratios are obtained using \code{calcLagged} with \code{lag = 1}.
}
\author{Chris Brien}
\seealso{\code{\link{intervalGRaverage}}, \code{\link{intervalGRdiff}}, \code{\link{splitContGRdiff}}, \code{\link{splitSplines}}, \code{calcLagged}}
\examples{
data(exampleData)
longi.dat$Area.AGR <- with(longi.dat, AGRdiff(Area, time.diffs = Days.diffs))
}
\keyword{data}
\keyword{manip}