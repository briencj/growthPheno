\name{plotSmoothsDevnBoxplots}
\alias{plotSmoothsDevnBoxplots}
\title{Produces boxplots for several sets of deviations of the smoothed values from a response, possibly along with growth rates.}
\description{ Calculates and produces, using \code{\link{plotDeviationsBoxes}}, boxplots 
             of the deviations of the supplied smoothed values from the observed response values for the 
             traits and for combinations of the different smoothing parameters and for subsets of 
             non-smoothing-\code{\link{factor}} combinations. Which traits are plotted is controlled by  
             \code{trait.types} and may include the (\code{response}and the computed traits of the 
             Absolute Growth Rates (AGR) and/or the Relative Growth Rates (RGR). The observed and smoothed 
             values are supplied in long format i.e. with the values for each set of smoothing parameters  
             stacked one under the other in the supplied \code{\link{smooths.frame}}. Such data can be 
             generated using \code{\link{probeSmooths}}.

             Multiple plots, possibly each having multiple facets, are produced using \code{ggplot2}. 
             The layout of these plots is controlled via the arguments \code{plots.by}, 
             \code{facet.x} and \code{facet.y}. The basic principle is that the number of levels 
             combinations of the smoothing-parameter \code{\link{factor}}s \code{Type}, \code{TunePar}, 
             \code{TuneVal}, \code{Tuning} (the combination of (\code{TunePar} and \code{TuneVal}), and 
             \code{Method} that are included  in \code{plots.by}, \code{facet.x} and \code{facet.y} 
             must be the same as those covered by the combinations of the values incorporated 
             into the \code{\link{smooths.frame}} input to \code{plotSmoothsDevnBoxplots} via the 
             \code{data} argument. This ensures that smooths from different parameter sets are not 
             pooled into the same plot. The \code{\link{factor}}s other than the smoothing-parameter 
             \code{\link{factor}}s can be supplied to the \code{plots.by} and \code{facet} arguments.
}
\usage{
plotSmoothsDevnBoxplots(data, response, response.smoothed = NULL, 
                        individuals = "Snapshot.ID.Tag", times = "DAP", 
                        trait.types = c("response", "AGR", "RGR"), 
                        which.plots =  "absolute.boxplots", 
                        x.title = NULL, y.titles = NULL, 
                        devnboxes.plot.args = 
                          args4devnboxes_plot(plots.by = NULL, 
                                              facet.x = ".", facet.y = "."),
                        printPlot = TRUE, ...)
}
\arguments{
 \item{data}{A \code{\link{smooths.frame}}, such as is produced by 
             \code{\link{probeSmooths}} and that contains the data resulting from 
             smoothing a response over time for a set of \code{individuals}, the data 
             being arranged in long format both with respect to the 
             times and the smoothing-parameter values used in the smoothing. That is, 
             each response occupies a single column. The unsmoothed \code{response} and 
             the \code{response.smoothed} are to be plotted for different sets of values 
             for the smoothing parameters. The \code{\link{smooths.frame}} must include 
             the columns \code{Type}, \code{TunePar}, \code{TuneVal}, \code{Tuning} and 
             \code{Method}, and the columns nominated using the arguments 
             \code{individuals}, \code{times}, \code{plots.by}, \code{facet.x}, \code{facet.y}, 
             \code{response}, \code{response.smoothed}, and, if requested, 
             the AGR and the RGR of the \code{response} and \code{response.smoothed}. 
             The names of the growth rates should be formed from \code{response} and 
             \code{response.smoothed} by adding \code{.AGR} and \code{.RGR} to both of them.}
 \item{response}{A \code{\link{character}} specifying the response variable for which the 
             observed values are supplied.}
 \item{response.smoothed}{A \code{\link{character}} specifying the name of the column 
             containing the values of the smoothed response variable, corresponding 
             to \code{response} and obtained for the combinations of 
             \code{smoothing.methods} and \code{df}, usually using smoothing splines. 
             If \code{response.smoothed} is \code{NULL}, then 
             \code{response.smoothed} is set to the \code{response} to which is added 
             the prefix \code{s}. }
 \item{times}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the times at which the data was 
             collected, either as a \code{\link{numeric}}, \code{\link{factor}}, or 
             \code{\link{character}}. It will be used to provide the values to be plotted 
             on the x-axis. If a \code{\link{factor}} or \code{\link{character}}, 
             the values should be numerics stored as characters.}
 \item{individuals}{A \code{\link{character}} giving the name of the 
            \code{\link{factor}} that defines the subsets of the \code{data} 
            for which each subset corresponds to the \code{response} values for 
            an individual (e.g. plant, pot, cart, plot or unit).}
 \item{trait.types}{A \code{\link{character}} giving the \code{trait.types} that 
             are to be plotted. If \code{all}, each of \code{response}, \code{AGR} and 
             \code{RGR} is plotted.}
\item{which.plots}{A \code{\link{logical}} indicating which plots are to be produced. 
             The options are either \code{none} or \code{absolute.deviations} and/or 
             \code{relative.deviations}. Boxplots of the absolute deviations 
             are specified by \code{absolute.boxplots}, the absolute deviations being 
             the values of a trait minus their smoothed values (observed - smoothed). 
             Boxplots of the relative deviations are specified by \code{relative.boxplots}, 
             the relative deviations being the absolute deviations divided by the 
             smoothed values of the trait.}
 \item{x.title}{Title for the x-axis, used for all plots. If \code{NULL} then set to 
             \code{times}.}
 \item{y.titles}{A \code{\link{character}} giving the titles for the y-axis, 
             one for each trait specified by \code{trait.types} and used for all plots. 
             If \code{NULL}, then set to the traits derived for \code{response} 
             from \code{trait.types}.}
 \item{devnboxes.plot.args}{A named \code{\link{list}} that is most easily 
             generated using \code{\link{args4devnboxes_plot}}, it documenting the 
             options available for varying the boxplots. \emph{Note 
             that if \code{\link{args4devnboxes_plot}} is to be called to change 
             from the default settings given in the default \code{probeSmooths} 
             call and some of those settings are to be retained, then the arguments 
             whose settings are to be retained must also be included in the call 
             to \code{\link{args4devnboxes_plot}}; be aware that if you call 
             \code{\link{args4devnboxes_plot}}, then the defaults for this call 
             are those for \code{\link{args4devnboxes_plot}}, \bold{NOT} the call 
             to \code{\link{args4devnboxes_plot}} shown as the default for 
             \code{probeSmooths}.}}
 \item{printPlot}{A \code{\link{logical}} indicating whether or not to print any 
             plots.}
 \item{...}{allows passing of arguments to \code{\link{plotProfiles}}.}
}
\value{A multilevel \code{\link{list}} that contains the \code{\link{ggplot}} 
        objects for the plots produced. The first-level \code{\link{list}} 
        has a component for each \code{trait.types} and each of these is a 
        second-level \code{\link{list}} with contains the deviations boxplots 
        for a \code{response}. Each plot is in an object of class 
        \code{\link{ggplot}}, which can be plotted using \code{print}. 
}
\author{Chris Brien}
\seealso{\code{\link{traitSmooth}}, \code{\link{probeSmooths}}, \code{\link{args4profile_plot}}, \code{\link{plotDeviationsBoxes}}, \code{\link{plotSmoothsMedianDevns}}, \code{\link{ggplot}}.}
\examples{
data(exampleData)
traits <- probeSmooths(data = longi.dat, 
                       response = "PSA", response.smoothed = "sPSA", 
                       times = "DAP", 
                       #only df is changed from the probeSmooth default
                       smoothing.args = 
                         args4smoothing(smoothing.methods = "direct", 
                                        spline.types = "NCSS", 
                                        df = c(4,7), lambdas = NULL), 
                       which.plots = "none")
plotSmoothsDevnBoxplots(data = traits, 
                        response = "PSA", response.smoothed = "sPSA", 
                        times = "DAP", x.title = "DAP",
                        #only facet.x is changed from the probeSmooth default
                        devnboxes.plot.args = 
                          args4devnboxes_plot(plots.by = NULL, 
                                              facet.x = "Tuning", facet.y = "."))
}
\keyword{hplot}
\keyword{manip}