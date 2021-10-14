\name{plotLongitudinal}
\alias{plotLongitudinal}
\title{Plots longitudinal data for a set of indiividuals}
\description{Produce profile or longitudinal plots of a response using \code{ggplot}. 
             A line is drawn for the data for each \code{individual} and the plot 
             can be facetted so that a grid of plots is produced. For each facet a line for 
             the medians over time can be added, along with the vaue of the outer whiskers 
             (median +/- 1.5 * IQR).
}
\usage{
plotLongitudinal(data, x = "xDays+44.5", response = "Area", 
                 individuals = "Snapshot.ID.Tag", title = NULL, 
                 x.title = "Days", y.title = "Area (kpixels)", 
                 facet.x = "Treatment.1", facet.y = "Smarthouse", 
                 labeller = NULL, colour = "black", 
                 colour.column = NULL, colour.values = NULL, 
                 alpha = 0.1, addMediansWhiskers = FALSE, 
                 xname = "xDays", ggplotFuncs = NULL, 
                 printPlot = TRUE)
}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the data to be plotted.}
 \item{x}{A \code{\link{character}} giving the variable to be plotted on the x-axis.}
 \item{response}{A \code{character} specifying the response variable that 
                 is to be plotted on the y-axis.}
 \item{individuals}{A \code{\link{character}} giving the name(s) of the 
             \code{\link{factor}}(s) that define the subsets of the 
             \code{data} for which each subset corresponds to the 
             \code{response} values for an individual.}
 \item{x.title}{Title for the x-axis.}
 \item{y.title}{Title for the y-axis.}
 \item{title}{Title for the plot.}
 \item{facet.x}{A \code{\link{data.frame}} giving the variable to be used to 
                form subsets to be plotted in separate columns of plots. 
                Use \code{"."} if a split into columns is not wanted.}
 \item{facet.y}{A \code{\link{data.frame}} giving the variable to be used to 
                form subsets to be plotted in separate rows of plots.
                Use \code{"."} if a split into columns is not wanted.}
 \item{labeller}{A \code{\link{ggplot}} \code{\link{function}} for labelling the 
                 facets of a plot produced using the \code{\link{ggplot}} function. 
                 For more information see \code{\link{ggplot}}.}
 \item{colour}{A \code{\link{character}} specifying a single colour to use in 
               drawing the lines for the profiles. If colouring according to the 
               values of a variable is required then use \code{colour.column}.}
 \item{colour.column}{A \code{\link{character}} giving the name of a column 
             in \code{data} over whose values the colours of the lines are to be 
             varied. The colours can be specified using \code{colour.values} .}
 \item{colour.values}{A \code{\link{character}} vector specifying the values of 
               the colours to use in drawing the lines for the profiles.  
               If this is a named vector, then the values will be matched based 
               on the names. If unnamed, values will be matched in order 
               (usually alphabetical) with the limits of the scale.}
 \item{alpha}{A \code{\link{numeric}} specifying the degrees of transparency to 
              be used in plotting. It is a ratio in which the denominator 
              specifies the number of points (or lines) that must be overplotted 
              to give a solid cover.}
 \item{addMediansWhiskers}{A \code{\link{logical}} indicating whether plots over 
             time of the medians and outer whiskers are to be added to the plot. 
             The outer whiskers are related to the whiskers on a box-and-whisker 
             and are defined as the median plus (and minus) 1.5 times the 
             interquartile range (IQR). Points lying outside the whiskers are 
             considered to be potential outliers.}
 \item{xname}{A \code{\link{character}} giving the name of the
             \code{\link{numeric}} that contains the values of the predictor 
              variable from which \code{x} is derived, it being that \code{x} 
             may incorporate an expression.}
 \item{ggplotFuncs}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \code{\link{ggplot}} function. 
                It is created by calling the \code{\link{list}} function with 
                a \code{\link{ggplot}} function call for each element.
                These functions are applied in creating the \code{ggplot} object.}
 \item{printPlot}{A \code{\link{logical}} indicating whether or not to print the 
                  plot.}
}
\value{An object of class "\code{\link{ggplot}}", which can be plotted using 
       \code{print}.}
\author{Chris Brien}
\seealso{\code{\link{ggplot}}, \code{\link{labeller}}.}
\examples{
data(exampleData)
plotLongitudinal(data = longi.dat, response = "Area.smooth")

plt <- plotLongitudinal(data = longi.dat, response = "Area.smooth", x.title = "DAP",  
                        y.title = "Area.smooth", x="xDays+35.42857143", printPlot=FALSE)
plt <- plt + ggplot2::geom_vline(xintercept=29, linetype="longdash", size=1) +
             ggplot2::scale_x_continuous(breaks=seq(28, 42, by=2)) + 
             ggplot2::scale_y_continuous(limits=c(0,750))
print(plt)

plotLongitudinal(data = longi.dat, response = "Area.smooth", x.title = "DAP",  
                 y.title = "Area.smooth", x="xDays+35.42857143", 
                 ggplotFuncs = list(ggplot2::geom_vline(xintercept=29, linetype="longdash", 
                                                        size=1), 
                                    ggplot2::scale_x_continuous(breaks=seq(28, 42, by=2)), 
                                    ggplot2::scale_y_continuous(limits=c(0,750))))
}
\keyword{hplot}

