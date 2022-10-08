\name{args4chosen_smooth}
\alias{args4chosen_smooth}
\title{Creates a \code{list} of the values for the smoothing parameters for 
             which a smooth is to be extracted}
\description{Creates a \code{\link{list}} of the values for the smoothing parameters for 
             which a smooth is to be extracted. 
             The default settings for these are optimized for \code{\link{traitSmooth}}.}
\usage{
args4chosen_smooth(smoothing.methods = "logarithmic", 
                   spline.types = "PS", 
                   df = NULL, 
                   lambdas = NULL,
                   combinations = "single")}
\arguments{
 \item{smoothing.methods}{A \code{\link{character}} giving the smoothing method 
      for the chosen smooth. The two possibilites are (i) \code{"direct"}, for directly 
      smoothing the observed \code{response}, and (ii) \code{"logarithmic"}, for 
      smoothing the \code{log}-transformed \code{response}.}
 \item{spline.types}{A \code{\link{character}} giving the type of spline  
      for the chosen smooth. Currently, the possibilites are (i) \code{"NCSS"}, for natural 
      cubic smoothing splines, and (ii) \code{"PS"}, for P-splines.}
 \item{df}{A \code{\link{numeric}} with single value that specifies, for 
      natural cubic smoothing splines (\code{NCSS}), the desired equivalent 
      numbers of degrees of freedom of the chosen smooth (trace of the smoother matrix). 
      Lower values result in more smoothing. }
 \item{lambdas}{A named \code{\link{list}} or a \code{\link{numeric}} specifying 
      the positive penalties for which the chosen smooth is required.}
 \item{combinations}{Generally, this argument should be set to \code{single} so that 
      ony one value should be supplied to the functions arguments. Also, only one of 
      \code{df} or \code{lambdas} should be set.}
} 
\value{A named \code{\link{list}}.}
\author{Chris Brien}
\seealso{\code{\link{traitSmooth}} and \code{\link{probeSmooths}}.}
\examples{

args4chosen_smooth(smoothing.methods = "direct", 
                   spline.types = "NCSS", df = 4)
args4chosen_smooth(smoothing.methods = "log", 
                   spline.types = "PS", lambdas = 0.36)
}
\keyword{data}
\keyword{manip}