\name{plotDeviationsBoxes}
\alias{plotDeviationsBoxes}
\title{Produces boxplots of the deviations of the observed values from the smoothed values over values of x.}
\description{Produces boxplots of the deviations of the observed values from the smoothed values over values of x.}
\usage{
plotDeviationsBoxes(data, observed, smoothed, x.factor, 
                    x.title = NULL, y.titles = NULL,
                    facet.x = ".", facet.y = ".", labeller = NULL, 
                    df, deviations.plots = "absolute", 
                    ggplotFuncs = NULL, printPlot = TRUE, ...)
 }
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the observed and smoothed values 
             from which the deviations are to be computed.}
 \item{observed}{A \code{character} specifying the response variable for which the 
                 observed values are supplied.}
 \item{smoothed}{A \code{character} specifying the smoothed response variable, 
                 corresponding to \code{observed}, for which values are supplied.}
 \item{x.factor}{A \code{\link{character}} giving the \code{\link{factor}} to be 
                plotted on the x-axis.}
 \item{x.title}{Title for the x-axis. If \code{NULL} then set to \code{x}.}
 \item{y.titles}{A \code{\link{character}} giving the titles for the y-axis, 
             one for each plot specified \code{deviations.plots}.}
 \item{facet.x}{A \code{\link{data.frame}} giving the variable to be used to 
                form subsets to be plotted in separate columns of plots. 
                Use \code{"."} if a split into columns is not wanted. For 
                \code{which.plots} set to \code{methodcompare} or \code{dfcompare} 
                \code{facet.x.pf} is ignored.}
 \item{facet.y}{A \code{\link{data.frame}} giving the variable to be used to 
                form subsets to be plotted in separate rows of plots.
                Use \code{"."} if a split into columns is not wanted.}
 \item{labeller}{A \code{\link{ggplot}} \code{\link{function}} for labelling the 
                 facets of a plot produced using the \code{\link{ggplot}} function. 
                 For more information see \code{\link{ggplot}}.}
 \item{df}{A \code{\link{numeric}} specifying the smoothing degrees of freedom used in 
            producing the \code{response.smoothed} and which is to be used in labelling 
            the plot.}
 \item{deviations.plots}{A \code{\link{character}} specifying whether 
              \code{absolute} and/or \code{relative} deviations are to be plotted.}
 \item{ggplotFuncs}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \code{\link{ggplot}} function. 
                It is created by calling the \code{\link{list}} function with 
                a \code{\link{ggplot}} function call for each element.
                These functions are applied in creating the \code{ggplot} object for plotting.}
 \item{printPlot}{A \code{\link{logical}} indicating whether or not to print the 
                  plots.}
\item{...}{allows passing of arguments to \code{\link{ggplot}}.}
}
\value{A list whose components are named \code{absolute} and \code{relative}; a 
       component will contain an object of class "\code{ggplot}" when the plot has been 
       requested using the \code{deviations.plots} argument and a \code{NULL} otherwise. 
       The objects can be plotted using \code{print}.}
\author{Chris Brien}
\seealso{\code{\link{plotMedianDeviations}}, \code{\link{probeSmoothing}}, \code{\link{ggplot}}.}
\examples{
data(exampleData)

plotDeviationsBoxes(longi.dat, observed = "PSA", smoothed = "sPSA",
                    x.factor="DAP", facet.x.pf = ".", facet.y= ".", df =5)
}
\keyword{hplot}
\keyword{manip}