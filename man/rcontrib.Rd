\name{rcontrib}
\alias{rcontrib}
\title{Computes a measure of how correlated each variable in a set is with the other variable, 
       conditional on a nominated subset of them}
\description{A measure of how correlated a variable is with those in a set is given by the 
             square root of the sum of squares of the correlation coefficients between the 
             variables and the other variables in the set (Cummings, 2007). Here, the partial 
             correlation  between the subset of the variables listed in \code{response} that 
             are not listed in \code{include} is calculated from the partial correlation matrix 
             for the subset, adjusting for those variables in \code{include}. This is useful 
             for manually deciding which of the variables not in \code{include} should next be 
             added to it.}
\usage{rcontrib(responses, data, include = NULL)}
\arguments{
 \item{responses}{A \code{\link{character}} giving the names of the columns in 
             \code{data} from which the correlation measure is to be calculated.}
 \item{data}{A \code{\link{data.frame}} containing the columns of variables from which the 
             correlation measure is to be calculated.}
 \item{include}{A \code{\link{character}} giving the names of the columns in 
             \code{data} for the variables for which other variables are to be adjusted.}
}
\value{A \code{\link{numeric}} giving the correlation measures.
}
\author{Chris Brien}
\seealso{\code{\link{PVA}}, \code{\link{intervalPVA}}}
\references{
Cumming, J. A. and D. A. Wood (2007) Dimension reduction via principal variables. \emph{Computational Statistics
and Data Analysis}, \bold{52}, 550--565.
}
\examples{
data(exampleData)
responses <- c("Area","Area.SV","Area.TV", "Image.Biomass", "Max.Height","Centre.Mass",
               "Density", "Compactness.TV", "Compactness.SV")
h <-  rcontrib(responses, longi.dat, include = "Area")
}
\keyword{data}
\keyword{manip}