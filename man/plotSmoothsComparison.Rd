\name{plotSmoothsComparison}
\alias{plotSmoothsComparison}
\title{Plots  several sets of smoothed values for a response, possibly along with growth rates and optionally including the unsmoothed values, as well as deviations boxplots.}
\description{Plots the smoothed values for an observed \code{response} and, optionally, the unsmoothed 
             observed \code{response} using \code{\link{plotProfiles}}. Depending on the setting of 
             \code{trait.types} (\code{response}, \code{AGR} or \code{RGR}), the computed traits of the 
             Absolute Growth Rates (AGR) and/or the Relative Growth Rates (RGR) are plotted.  This 
             function will also calculate and produce, using \code{\link{plotDeviationsBoxes}}, boxplots 
             of the deviations of the supplied smoothed values from the observed response values for the 
             traits and for combinations of the different smoothing parameters and for subsets of 
             non-smoothing-\code{\link{factor}} combinations. The observed and smoothed values are 
             supplied in long format i.e. with the values for each set of smoothing parameters stacked 
             one under the other in the supplied \code{\link{smooths.frame}}. Such data can be generated 
             using \code{\link{probeSmooths}};  to prevent \code{\link{probeSmooths}} producing the 
             plots, which it is does using \code{plotSmoothsComparison}, \code{\link{plotDeviationsBoxes}} 
             and \code{\link{plotSmoothsMedianDevns}}, set \code{which.plots} to \code{none}. 
             The smoothing parameters include \code{spline.types}, \code{df}, \code{lambdas} and 
             \code{smoothing.methods} (see \code{\link{probeSmooths}}). 

             Multiple plots, possibly each having multiple facets, are produced using \code{ggplot2}. 
             The layout of these plots is controlled via the arguments \code{plots.by}, 
             \code{facet.x} and \code{facet.y}. The basic principle is that the number of levels 
             combinations of the smoothing-parameter \code{\link{factor}}s \code{Type}, \code{TunePar}, 
             \code{TuneVal}, \code{Tuning} (the combination of (\code{TunePar} and \code{TuneVal}), and 
             \code{Method} that are included  in \code{plots.by}, \code{facet.x} and 
             \code{facet.y} must be the same as those covered by the combinations of the values 
             supplied to \code{spline.types}, \code{df}, \code{lambdas} and \code{Method} and incorporated 
             into the \code{\link{smooths.frame}} input to \code{plotSmoothsComparison} via the 
             \code{data} argument. This ensures that smooths from different parameter sets are not 
             pooled into the same plot. The \code{\link{factor}}s other than the smoothing-parameter 
             \code{\link{factor}}s can be supplied to the \code{plots.by} and \code{facet} arguments.


             The following profiles plots can be produced: (i) separate plots of the 
             smoothed traits for each combination of the smoothing parameters 
             (include \code{Type}, \code{Tuning} and \code{Method} in \code{plots.by}); 
             (ii) as for (i), with the corresponding plot for the unsmoothed trait 
             preceeding the plots for the smoothed trait (also set \code{include.raw} to 
             \code{alone}); (iii) profiles plots that compare a smoothed trait for all 
             combinations of the values of the smoothing parameters, arranging the plots 
             side-by-side or one above the other (include \code{Type}, \code{Tuning} and 
             \code{Method} in \code{facet.x} and/or \code{facet.y} - to include the 
             unsmoothed trait set \code{include.raw} to one of \code{facet.x} or 
             \code{facet.y}; (iv) as for (iii), except that separate plots are 
             produced for each combination of the levels of the \code{\link{factor}}s 
             in \code{plot.by} and each plot compares the smoothed traits for the 
             smoothing-parameter \code{\link{factor}}s included in \code{facet.x} 
             and/or \code{facet.y} (set both \code{plots.by} and one or more of 
             \code{facet.x} and \code{facet.y}). 
}
\usage{
plotSmoothsComparison(data, response, response.smoothed = NULL, 
                      individuals = "Snapshot.ID.Tag", times = "DAP", 
                      trait.types = c("response", "AGR", "RGR"), 
                      x.title = NULL, y.titles = NULL, 
                      profile.plot.args = 
                        args4profile_plot(plots.by = NULL, 
                                          facet.x = ".", facet.y = ".", 
                                          include.raw = "no"),
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
             are to be plotted when \code{which.plots} is \code{profiles}. 
             Irrespective of the setting of \code{get.rates}, the nominated traits 
             are plotted. If \code{all}, each of \code{response}, \code{AGR} and 
             \code{RGR} is plotted.}
 \item{x.title}{Title for the x-axis, used for all plots. If \code{NULL} then set to 
             \code{times}.}
 \item{y.titles}{A \code{\link{character}} giving the titles for the y-axis, 
             one for each trait specified by \code{trait.types} and used for all plots. 
             If \code{NULL}, then set to the traits derived for \code{response} 
             from \code{trait.types}.}
 \item{profile.plot.args}{A named \code{\link{list}} that is most easily 
             generated using \code{\link{args4profile_plot}}, it documenting the 
             options available for varying profile plots and boxplots. \emph{Note 
             that if \code{\link{args4profile_plot}} is to be called to change 
             from the default settings given in the default 
             \code{plotSmoothsComparison} call 
             and some of those settings are to be retained, then the arguments 
             whose settings are to be retained must also be included in the call 
             to \code{\link{args4profile_plot}}; be aware that if you call 
             \code{\link{args4profile_plot}}, then the defaults for this call are 
             those for \code{\link{args4profile_plot}}, \bold{NOT} the call to 
             \code{\link{args4profile_plot}} shown as the default for 
             \code{plotSmoothsComparison}.}}
 \item{printPlot}{A \code{\link{logical}} indicating whether or not to print any 
             plots.}
 \item{...}{allows passing of arguments to \code{\link{plotProfiles}}.}
}
\value{A multilevel \code{\link{list}} that contains the \code{\link{ggplot}} 
        objects for the plots produced. The first-level \code{\link{list}} 
        has a component for each \code{trait.types} and each of these is a 
        second-level \code{\link{list}} that contains the trait 
        profile plots and for a \code{trait}. It may contain components labelled 
        \code{Unsmoothed}, \code{all} or for one of the levels of the 
        \code{\link{factor}}s in \code{plots.by}; each of these third-level 
        \code{\link{list}}s contains a \code{\link{ggplot}} object that can 
        be plotted using \code{print}. 
}
\author{Chris Brien}
\seealso{\code{\link{traitSmooth}}, \code{\link{probeSmooths}}, \code{\link{args4profile_plot}}, \code{\link{plotDeviationsBoxes}}, \code{\link{plotSmoothsMedianDevns}}, \code{\link{ggplot}}.}
\examples{
data(exampleData)
vline <- list(ggplot2::geom_vline(xintercept=29, linetype="longdash", size=1))
traits <- probeSmooths(data = longi.dat, 
                       response = "PSA", response.smoothed = "sPSA", 
                       times = "DAP", 
                       #only df is changed from the probeSmooth default
                       smoothing.args = 
                         args4smoothing(smoothing.methods = "direct", 
                                        spline.types = "NCSS", 
                                        df = c(4,7), lambdas = NULL), 
                       which.plots = "none")
plotSmoothsComparison(data = traits, 
                      response = "PSA", response.smoothed = "sPSA", 
                      times = "DAP", x.title = "DAP",
                      #only facet.x is changed from the probeSmooth default
                      profile.plot.args = 
                        args4profile_plot(plots.by = NULL, 
                                          facet.x = "Tuning", facet.y = ".", 
                                          include.raw = "no", 
                                          ggplotFuncs = vline))
}
\keyword{hplot}
\keyword{manip}