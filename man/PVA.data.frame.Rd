\name{PVA.data.frame}
\alias{PVA.data.frame}
\title{Selects a subset of variables stored in a data.frame using Principal Variable Analysis (PVA)}
\description{Principal Variable Analysis (PVA) (Cumming and Wooff, 2007) selects a subset 
  from a set of the variables such that the variables in the subset are as 
  uncorrelated as possible, in an effort to ensure that all aspects of the 
  variation in the data are covered. }
\usage{
\method{PVA}{data.frame}(obj, responses, nvarselect = NULL, p.variance = 1, include = NULL, 
    plot = TRUE, ...)}
\arguments{
 \item{obj}{A \code{\link{data.frame}} containing the columns of variables from which the 
             selection is to be made.}
 \item{responses}{A \code{\link{character}} giving the names of the columns in 
             \code{data} from which the variables are to be selected.}
 \item{nvarselect}{A \code{\link{numeric}} specifying the number of variables to be selected, 
                   which includes those listed in \code{include}. If \code{nvarselect = 1}, as 
                   many variables are selected as is need to satisfy \code{p.variance}.}
 \item{p.variance}{A \code{\link{numeric}} specifying the minimum proportion of the variance 
                   that the selected variables must account for,}
 \item{include}{A \code{\link{character}} giving the names of the columns in 
             \code{data} for the variables whose selection is mandatory.}
 \item{plot}{A \code{\link{logical}} indicating whether a plot of the cumulative proportion 
            of the variance explained is to be produced.}
 \item{...}{allows passing of arguments to other functions}
}
\value{A \code{\link{data.frame}} giving the results of the variable selection. 
       It will contain the columns \code{Variable}, \code{Selected}, 
       \code{h.partial}, \code{Added.Propn} and \code{Cumulative.Propn}.
}
\details{
The variable that is most correlated with the other variables is selected first for inclusion. The partial correlation for each of the remaining variables, given the first selected variable, is calculated and the most correlated of these variables is selects for inclusion next. Then the partial correlations are adjust for the second included variables. This process is repeated until the specified criteria have been satisfied. The possibilities are:
\enumerate{
  \item the default (\code{nvarselect = NULL} and \code{p.variance = 1}), which selects all 
        variables in increasing order of amount of information they provide;
  \item to select exactly \code{nvarselect} variables;
  \item to select just enough variables, up to a maximum of \code{nvarselect} variables, to explain 
        at least \code{p.variance}*100 per cent of the total variance.
}
}
\author{Chris Brien}
\seealso{\code{\link{PVA}}, \code{\link{PVA.matrix}}, \code{\link{intervalPVA.data.frame}}, \code{\link{rcontrib}}}
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
results <-  PVA(longi.dat, responses, p.variance=0.9, plot = FALSE)
}
\keyword{data}
\keyword{manip}