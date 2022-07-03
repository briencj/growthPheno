\name{smooths.frame}
\alias{smooths.frame}
\alias{smooths.frame-class}
\title{Description of a smooths.frame object}
\description{A \code{data.frame} of S3-class \code{smooths.frame} that stores the 
             smooths of one or more responses for several sets of smoothing parameters.

             \code{\link{as.smooths.frame}} is function that converts a 
             \code{\link{data.frame}} to an object of this class.

             \code{\link{is.smooths.frame}} is the membership function for this class; it tests 
             that an object has class \code{smooths.frame}.

             \code{\link{validSmoothsFrame}} can be used to test the validity of a 
             \code{smooths.frame}.
}
\value{A \code{\link{data.frame}} that is also inherits the S3-class 
       \code{\link{smooths.frame}}. It contains the results of smoothing a response 
       over time from a set of \code{individuals}, the data being arranged in long 
       format both with respect to the times and the smoothing-parameter values used in 
       the smoothing. That is, each response occupies a single column. The 
       \code{\link{smooths.frame}} must include the columns \code{Type}, \code{TunePar}, 
       \code{TuneVal}, \code{Tuning} (the combination of \code{TunePar} and 
       \code{TuneVal}) and \code{Method}, and the columns that would be nominated using 
       the \code{\link{probeSmooths}} arguments \code{individuals}, the \code{plots} and 
       \code{facet} arguments, \code{times}, \code{response}, \code{response.smoothed}, and, 
       if requested, the AGR and the RGR of the \code{response} and \code{response.smoothed}. 
       The names of the growth rates should be formed from \code{response} and 
       \code{response.smoothed} by adding \code{.AGR} and \code{.RGR} to both of them. 
       The function \code{\link{probeSmooths}} produces a \code{\link{smooths.frame}} 
       for a response. 

       A \code{\link{smooths.frame}} has the following attributes: 
       \enumerate{
       \item \code{individuals}, the \code{\link{character}} giving the name of the 
            \code{\link{factor}} that define the subsets of the \code{data} 
            for which each subset corresponds to the \code{response} values for 
            an individual; 
       \item \code{n}, the number of unique \code{individuals}; 
       \item \code{times}, the \code{\link{character}} giving the name of the 
            \code{\link{numeric}}, or \code{\link{factor}} with numeric levels, that 
            contains the values of the predictor variable plotted on the x-axis;
       \item \code{t}, the number of unique values in the \code{times}; 
       \item \code{nschemes}, the number of unique combinations of the 
       smoothing-parameter values in the \code{smoothsframe}. }
}
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

data(exampleData)
vline <- list(ggplot2::geom_vline(xintercept=29, linetype="longdash", size=1))
smths <- probeSmooths(data = longi.dat, 
                      response = "PSA", response.smoothed = "sPSA", 
                      times = "DAP", df = c(4,7), 
                      facet.x.pf = "Tuning", facet.y.pf = "Treatment.1", 
                      ggplotFuncsProfile=vline)
is.smooths.frame(smths)
validSmoothsFrame(smths)
}

\author{Chris Brien}
\seealso{\code{\link{probeSmooths}}, \code{\link{is.smooths.frame}}, 
          \code{\link{as.smooths.frame}}, \code{\link{validSmoothsFrame}}}
\keyword{asreml}
\keyword{htest}