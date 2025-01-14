\name{plotCorrmatrix}
\alias{plotCorrmatrix}
\title{Calculates and plots correlation matrices for a set of responses}
\description{Having calculated the correlations a heat map indicating the magnitude of the 
       correlations is produced using \code{ggplot}. In this heat map, the darker the red in 
       a cell then the closer the correlation is to -1, while the deeper the blue in the cell,
       then the closer the correlation is to 1. Matrix plots of all pairwise 
       combinations of the variables can be produced that includes the values of the 
       the correlation coefficients. If \code{pairs.sets} is set, a matrix 
       plot, along with the values of the correlation coefficients, is produced  for 
       each of the \code{pair.sets}. That is, the argument \code{pairs.sets} can be used to 
       restrict the pairs in a matrix plot to those combinations within each set. 
}
\usage{
plotCorrmatrix(data, responses, which.plots = c("heatmap","matrixplots"), 
               title = NULL, labels = NULL, labelSize = 4, pairs.sets = NULL, 
               show.sig = TRUE, axis.text.size = 20, ggplotFuncs = NULL, 
               printPlot = TRUE, ...)
}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the columns of variables to be correlated.}
 \item{responses}{A \code{\link{character}} giving the names of the columns in 
             \code{data} containing the variables to be correlated.}
 \item{which.plots}{A \code{\link{character}} specifying the plots of the correlations to be 
                    produced. The possibilities are one or both of \code{heatmap} and 
                    \code{matrixplots}.}
 \item{title}{Title for the plots.}
 \item{labels}{A \code{\link{character}} specifying the labels to be used in the plots. 
               If \code{labels} is \code{NULL}, \code{responses} is used for the \code{labels}.}
 \item{labelSize}{A \code{\link{numeric}} giving the size of the \code{labels} in the 
                  \code{matrixplots}.}
\item{pairs.sets}{A \code{\link{list}} each of whose components is a \code{\link{numeric}} giving 
                the position of the variable names in \code{responses} that are to be included 
                in the set. All pairs of variables in this pairs.set will be included in  
                \code{matrixplots}.}
 \item{show.sig}{A \code{\link{logical}} indicating whetherto give asterisks on the 
                 \code{heatmap} and \code{matrixplots} that indicate that the correlations are 
                 significantly different from zero.}
 \item{axis.text.size}{A \code{\link{numeric}} giving the size of the \code{labels} on the axes   
                  of the \code{heatmap}.}
 \item{ggplotFuncs}{A \code{list}, each element of which contains the 
                results of evaluating a \code{ggplot} function. 
                It is created by calling the \code{list} function with 
                a \code{ggplot} function call for each element.
                These functions are applied in creating the \code{ggplot} 
                object.}
 \item{printPlot}{A \code{\link{logical}} indicating whether or not to print the 
                  plots.}
 \item{...}{allows passing of arguments to other functions; not used at present.}
}
\value{A \code{\link{list}} object that has components \code{heatmap} and \code{matrixplots}.  The component \code{heatmap} will contain the heatmap plot, if produced, as an object of class "\code{ggplot}", which can be plotted using \code{print}; otherwise \code{NULL} is returned. 
Similarly, if not \code{NULL}, the component \code{matrixplots} will contain a list with one or more components, depending on the setting of \code{pair.sets}, each of which is a scatterplot matrix stored as an object of class "\code{ggmatrix}".}
\details{The correlations and their p-values are producced using \code{rcorr} 
         from the \code{Hmisc} package. The \code{heatmap} is produced using 
         \code{ggplot} and the matrixplots are produced using \code{GGally}.
}
\author{Chris Brien}
\seealso{\code{rcorr}, \code{GGally}, \code{ggplot}.}
\examples{\donttest{
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
plotCorrmatrix(longi.dat, responses, pairs.sets=list(c(1:4),c(5:7)))

}}
\keyword{hplot}

