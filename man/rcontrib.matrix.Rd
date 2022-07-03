\name{rcontrib.matrix}
\alias{rcontrib.matrix}
\title{Computes a measure of how correlated each variable in a set is with the other variable, 
       conditional on a nominated subset of them}
\description{A measure of how correlated a variable is with those in a set is given by the 
             square root of the sum of squares of the correlation coefficients between the 
             variables and the other variables in the set (Cumming and Wooff, 2007). Here, the partial 
             correlation  between the subset of the variables listed in \code{response} that 
             are not listed in \code{include} is calculated from the partial correlation matrix 
             for the subset, adjusting for those variables in \code{include}. This is useful 
             for manually deciding which of the variables not in \code{include} should next be 
             added to it.}
\usage{
\method{rcontrib}{matrix}(obj, responses, include = NULL, ...)}
\arguments{
 \item{obj}{A \code{\link{matrix}} containing the correlations of the variables from which the 
             correlation measure is to be calculated.}
 \item{responses}{A \code{\link{character}} giving the names of the columns in 
             \code{data} from which the correlation measure is to be calculated.}
 \item{include}{A \code{\link{character}} giving the names of the columns in 
             \code{data} for the variables for which other variables are to be adjusted.}
 \item{...}{allows passing of arguments to other functions.}
}
\value{A \code{\link{numeric}} giving the correlation measures.
}
\author{Chris Brien}
\seealso{\code{\link{rcontrib}}, \code{\link{rcontrib.data.frame}}, \code{\link{PVA}}, \code{\link{intervalPVA.data.frame}}}
\references{
Cumming, J. A. and D. A. Wooff (2007) Dimension reduction via principal variables. \emph{Computational Statistics
and Data Analysis}, \bold{52}, 550--565.
}
\examples{
data(exampleData)
longi.dat <- prepImageData(data=raw.dat, smarthouse.lev=1)
longi.dat <- within(longi.dat, 
                    {
                      Max.Height <- pmax(Max.Dist.Above.Horizon.Line.SV1,  
                                         Max.Dist.Above.Horizon.Line.SV2)
                      Density <- PSA/Max.Height
                      PSA.SV = (PSA.SV1 + PSA.SV2) / 2
                      Image.Biomass = PSA.SV * (PSA.TV^0.5)
                      Centre.Mass <- (Center.Of.Mass.Y.SV1 + Center.Of.Mass.Y.SV2) / 2
                      Compactness.SV = (Compactness.SV1 + Compactness.SV2) / 2
                    })
responses <- c("PSA","PSA.SV","PSA.TV", "Image.Biomass", "Max.Height","Centre.Mass",
               "Density", "Compactness.TV", "Compactness.SV")
R <- Hmisc::rcorr(as.matrix(longi.dat[responses]))$r
h <- rcontrib(R, responses, include = "PSA")
}
\keyword{data}
\keyword{manip}