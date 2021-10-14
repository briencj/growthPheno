\name{plotCorrmatrix}
\alias{plotCorrmatrix}
\title{Calculates and plots correlation matrices for a set of responses}
\description{Having calculated the correlations a heat map indicating the magnitude of the 
       correlations is produced using \code{ggplot}. In this heat map, the darker the red in 
       a cell then the closer the correlation is to -1, while the deeper the blue in the cell,
       then the closer the correlation is to 1. A matrix plot of all pairwise 
       combinations of the variables can be produced. The matrix plot contains a scatter diagram 
       for each pair, as well as the value of the correlation coefficient. The argument 
       \code{pairs.sets} can be used to restrict the pairs in the matrix plot to those 
       combinations within each set. 
}
\usage{
plotCorrmatrix(data, responses, which.plots = c("heatmap","matrixplot"), 
               title = NULL, labels = NULL, labelSize = 4, pairs.sets = NULL, 
               show.sig = FALSE, axis.text.size = 20, ggplotFuncs = NULL, 
               printPlot = TRUE, ...)
}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the columns of variables to be correlated.}
 \item{responses}{A \code{\link{character}} giving the names of the columns in 
             \code{data} containing the variables to be correlated.}
 \item{which.plots}{A \code{\link{character}} specifying the plots of the correlations to be 
                    produced. The possibilities are one or both of \code{heatmap} and 
                    \code{matrixplot}.}
 \item{title}{Title for the plots.}
 \item{labels}{A \code{\link{character}} specifying the labels to be used in the plots. 
               If \code{labels} is \code{NULL}, \code{responses} is used for the \code{labels}.}
 \item{labelSize}{A \code{\link{numeric}} giving the size of the \code{labels} in the 
                  \code{matrixplot}.}
\item{pairs.sets}{A \code{\link{list}} each of whose components is a \code{\link{numeric}} giving 
                the position of the variable names in \code{responses} that are to be included 
                in the set. All pairs of variables in this pairs.set will be included in a 
                \code{matrixplot}.}
 \item{show.sig}{A \code{\link{logical}} indicating whether or not to give asterisks on the 
                 \code{heatmap} indicating the correlations are significantly different from zero.}
 \item{axis.text.size}{A \code{\link{numeric}} giving the size of the \code{labels} on the axes   
                  of the \code{heatmap}.}
 \item{ggplotFuncs}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \code{\link{ggplot}} function. 
                It is created by calling the \code{\link{list}} function with 
                a \code{\link{ggplot}} function call for each element.
                These functions are applied in creating the \code{ggplot} object.}
 \item{printPlot}{A \code{\link{logical}} indicating whether or not to print the 
                  plot.}
 \item{...}{allows passing of arguments to other functions; not used at present.}
}
\value{The \code{heatmap} plot, if produced, as an object of class "\code{ggplot}", which 
       can be plotted using \code{print}; otherwise \code{NULL} is returned.}
\details{The correlations and their p-values are producced using \code{\link{rcorr}} 
         from the \code{Hmisc} package. The \code{heatmap} is produced using 
         \code{\link{ggplot}} and the matrixplot is produced using \code{\link{GGally}}.
}
\author{Chris Brien}
\seealso{\code{\link{rcorr}}, \code{\link{GGally}}, \code{\link{ggplot}}.}
\examples{\donttest{
data(exampleData)
responses <- c("Area","Area.SV","Area.TV", "Image.Biomass", "Max.Height","Centre.Mass",
               "Density", "Compactness.TV", "Compactness.SV")
plotCorrmatrix(longi.dat, responses, pairs.sets=list(c(1:4),c(5:7)))

}}
\keyword{hplot}

