\name{plotProfiles}
\alias{plotProfiles}
\title{Produces profile plots of longitudinal data for a set of individuals}
\description{Produce profile plots of longitudinal data for a response using \code{ggplot}. 
             A line is drawn for the data for each \code{individual} and the plot 
             can be faceted so that a grid of plots is produced. For each facet a line for 
             the medians over time can be added, along with the vaue of the outer whiskers 
             (median +/- 1.5 * IQR).
}
\usage{
plotProfiles(data, response = "PSA", 
             individuals = "Snapshot.ID.Tag", times = "DAP", 
             x = NULL, title = NULL, 
             x.title = "DAP", y.title = "PSA (kpixels)", 
             facet.x = ".", facet.y = ".", 
             labeller = NULL, scales = "fixed", 
             breaks.spacing.x = -2, angle.x = 0, 
             colour = "black", 
             colour.column = NULL, colour.values = NULL, 
             alpha = 0.1, addMediansWhiskers = FALSE, 
             ggplotFuncs = NULL, 
             printPlot = TRUE)
}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the data to be plotted.}
 \item{response}{A \code{character} specifying the response variable that 
                 is to be plotted on the y-axis.}
 \item{individuals}{A \code{\link{character}} giving the name of the 
            \code{\link{factor}} that defines the subsets of the \code{data} 
            for which each subset corresponds to the \code{response} values for 
            an individual (e.g. plant, pot, cart, plot or unit).}
 \item{times}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the times at which the data was collected, 
             either as a \code{\link{numeric}}, \code{\link{factor}}, or 
             \code{\link{character}}. If not a \code{\link{numeric}}, it will be 
             converted to a \code{\link{numeric}} and used to provide the values 
             to be plotted on the x-axis. If a \code{\link{factor}} or 
             \code{\link{character}}, the values should be numerics stored as 
             characters.}
 \item{x}{A \code{\link{character}} specifying a variable, or a function of 
             variables, to be plotted on the x-axis. If \code{NULL}, it will be 
             set to the value of \code{times}, which it can be assumed will be 
             converted to a \code{\link{numeric}}.}
 \item{x.title}{Title for the x-axis.}
 \item{y.title}{Title for the y-axis.}
 \item{title}{Title for the plot.}
 \item{facet.x}{A \code{\link{data.frame}} giving the variable to be used to 
                form subsets to be plotted in separate columns of plots. 
                Use \code{"."} if a split into columns is not wanted.}
 \item{facet.y}{A \code{\link{data.frame}} giving the variable to be used to 
                form subsets to be plotted in separate rows of plots.
                Use \code{"."} if a split into rows is not wanted.}
 \item{labeller}{A \code{\link{ggplot}} \code{\link{function}} for labelling the 
                 facets of a plot produced using the \code{\link{ggplot}} function. 
                 For more information see \code{\link{ggplot}}.}
 \item{scales}{A \code{\link{character}} specifying whether the scales are shared 
              across all facets of a plot (the default, "fixed"), or do they vary across 
              rows ("free_x"), columns ("free_y"), or both rows and columns ("free")?}
 \item{breaks.spacing.x}{A \code{\link{numeric}} whose absolute values specifies the 
             distance between major breaks for the x-axis in a sequence beginning with the 
             minimum x value and continuing up to the maximum x value. If it is negative, 
             the breaks that do not have x values in \code{data} will be omitted. 
             Minor breaks will be at half major break value or, if these do not correspond to 
             x-values in \code{data} when \code{breaks.spacing.x} is negative, have 
             a spacing of one. Thus, when \code{breaks.spacing.x} is negative,  
             grid lines will only be included for x-values that occur in \code{data}.  
             These settings can be overwritten by supplying, in \code{ggplotFuncs}, 
             a \code{scale_x_continuous} function from \code{ggplot2}.}
 \item{angle.x}{A \code{\link{numeric}} between 0 and 360 that gives the angle of the 
             x-axis text to the x-axis. It can also be set by supplying, in 
             \code{ggplotFuncs}, a \code{theme} function from \code{ggplot2}.}
 \item{colour}{A \code{\link{character}} specifying a single colour to use in 
               drawing the lines for the profiles. If colouring according to the 
               values of a variable is required then use \code{colour.column}.}
 \item{colour.column}{A \code{\link{character}} giving the name of a column 
             in \code{data} over whose values the colours of the lines are to be 
             varied. The colours can be specified using \code{colour.values}.}
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
plotProfiles(data = longi.dat, response = "sPSA", times = "DAP")

plt <- plotProfiles(data = longi.dat, response = "sPSA", 
                    y.title = "sPSA (kpixels)", 
                    facet.x = "Treatment.1", facet.y = "Smarthouse", 
                    breaks.spacing.x = 2,
                    printPlot=FALSE)
plt <- plt + ggplot2::geom_vline(xintercept=29, linetype="longdash", linewidth=1) +
             ggplot2::scale_y_continuous(limits=c(0,750))
print(plt)

plotProfiles(data = longi.dat, response = "sPSA", times = "DAP", 
             x.title = "DAP", y.title = "sPSA (kpixels)", 
             facet.x = "Treatment.1", facet.y = "Smarthouse", 
                 ggplotFuncs = list(ggplot2::geom_vline(xintercept=29, 
                                                        linetype="longdash", 
                                                        size=1), 
                                    ggplot2::scale_x_continuous(breaks=seq(28, 42, 
                                                                           by=2)), 
                                    ggplot2::scale_y_continuous(limits=c(0,750))))
}
\keyword{hplot}

