\name{plotSmoothsMedianDevns}
\alias{plotSmoothsMedianDevns}
\title{Calculates and plots the medians of the deviations from the observed values for 
       several sets of smoothed values  stored in a \code{data.frame} in long format.}
\description{Calculates and plots the medians of the deviations of the supplied smoothed values from the 
             supplied observed values for traits and combinations of different smoothing parameters, 
             possibly for subsets of non-smoothing-\code{\link{factor}} combinations. The observed and smoothed values are supplied in 
             long format i.e. with the values for each set of smoothing parameters stacked one under the other 
             in the supplied \code{\link{data.frame}}. Such data can be generated using \code{\link{probeSmooths}};  
             to prevent \code{\link{probeSmooths}} producing the plots, which it is does using 
             \code{\link{plotSmoothsComparison}}, \code{\link{plotDeviationsBoxes}} and 
             \code{plotSmoothsMedianDevns}, set \code{which.plots} to \code{none}. 
             The smoothing parameters include \code{spline.types}, \code{df}, \code{lambdas} and 
             \code{smoothing.methods} (see \code{\link{probeSmooths}}). 

             Multiple plots, possibly each having multiple facets, are produced using \code{ggplot2}. 
             The layout of these plots is controlled via the smoothing-parameter \code{\link{factor}}s 
             \code{Type}, \code{Tuning} (the combination of \code{TunePar} and \code{TuneVal}) and 
             \code{Method} that can be supplied to the arguments \code{plots.by.med}, 
             \code{plots.group.med}, \code{facet.x.med} and \code{facet.y.med}. These \code{plots} and 
             \code{facet} arguments can also include \code{\link{factor}}s other than the 
             smoothing-parameter \code{\link{factor}}s, that are also associated with the data. 
             The basic principle is that the number of levels combinations of the smoothing-parameter 
             \code{\link{factor}}s included in the \code{plots} and \code{facet} 
             arguments must be the same as those covered by the combinations of the values supplied to 
             \code{spline.types}, \code{df}, \code{lambdas} and \code{Method} and incorporated into the 
             \code{\link{smooths.frame}} input to \code{plotSmoothsMedianDevns} via the \code{data} argument. 
             This ensures that smooths from different parameter sets are not pooled in a single plot. 
             Envelopes of the median value of a trait for each \code{\link{factor}} combination can be added.}
\usage{
plotSmoothsMedianDevns(data, response, response.smoothed = NULL, 
                       individuals = "Snapshot.ID.Tag",  times = "DAP", 
                       trait.types = c("response", "AGR", "RGR"), 
                       x.title = NULL, y.titles = NULL, labeller = NULL, 
                       breaks.spacing.x = 4, 
                       plots.by.med = NULL, plots.group.med = NULL, 
                       facet.x.med = ".", facet.y.med = ".", 
                       colour.values.med = NULL, shape.values.med = NULL, 
                       alpha.med = 0.5, 
                       propn.note.med = TRUE, 
                       propn.types.med = c(0.1, 0.5, 0.75), 
                       ggplotFuncsMedDevn = NULL, printPlot = TRUE, ...)
}
\arguments{
 \item{data}{A \code{\link{smooths.frame}}, such as is produced by 
             \code{\link{probeSmooths}} and that contains the data resulting from 
             smoothing a response over time for a set of \code{individuals}, the data 
             being arranged in long format both with respect to the 
             times and the smoothing-parameter values used in the smoothing. That is, 
             each response occupies a single column. The \code{\link{smooths.frame}} 
             must include the columns \code{Type}, \code{TunePar}, \code{TuneVal}, 
             \code{Tuning} and \code{Method}, and the columns nominated using the 
             arguments \code{individuals}, \code{times}, \code{plots.by.med}, 
             \code{facet.x.med}, \code{facet.y.med}, \code{plots.group.med}, \code{response}, 
             \code{response.smoothed}, and, if requested, the AGR and the RGR of the 
             \code{response} and \code{response.smoothed}. The names of the growth 
             rates should be formed from \code{response} and \code{response.smoothed} 
             by adding \code{.AGR} and \code{.RGR} to both of them.}
 \item{response}{A \code{\link{character}} specifying the response variable for which the 
             observed values are supplied. Depending on the setting of 
             \code{trait.types}, the observed values of related \code{trait.types} 
             may also need to be be supplied.}
 \item{response.smoothed}{A \code{\link{character}} specifying the name of the column 
             containing the values of the smoothed response variable, corresponding 
             to \code{response} and obtained for the combinations of 
             \code{smoothing.methods} and \code{df}, usually using smoothing splines. 
             If \code{response.smoothed} is \code{NULL}, then 
             \code{response.smoothed} is set to the \code{response} to which is added 
             the prefix \code{s}. 
             Depending on the setting of \code{trait.types}, the smoothed 
             values of related \code{trait.types} may also need to be be supplied.}
 \item{individuals}{A \code{\link{character}} giving the name of the 
            \code{\link{factor}} that defines the subsets of the \code{data} 
            for which each subset corresponds to the \code{response} values for 
            an individual (e.g. plant, pot, cart, plot or unit).}
 \item{times}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the times at which the data was 
             collected, either as a \code{\link{numeric}}, \code{\link{factor}}, or 
             \code{\link{character}}. It will be used to provide the values to be plotted 
             on the x-axis. If a \code{\link{factor}} or \code{\link{character}}, 
             the values should be numerics stored as characters.}
 \item{trait.types}{A \code{\link{character}} giving the traits types that are to be 
             plotted. While \code{AGR} and \code{RGR} are commonly used, the names 
             can be arbitrary, except that \code{response} is a special case that 
             indicates that the original response is to be plotted. If \code{all}, 
             each of \code{response}, \code{AGR} and \code{RGR} is plotted.}
 \item{x.title}{Title for the x-axis. If \code{NULL} then set to \code{times}.}
 \item{y.titles}{A \code{\link{character}} giving the titles for the y-axis, 
             one for each trait specified by \code{trait.types}. 
             If \code{NULL}, then set to the traits derived for \code{response} 
             from \code{trait.types}.}
\item{labeller}{A \code{\link{ggplot}} \code{\link{function}} for labelling the 
             facets of a plot produced using the \code{\link{ggplot}} function. 
             For more information see \code{\link{ggplot}}.}
\item{breaks.spacing.x}{A \code{\link{numeric}} specifying the distance between breaks for 
             the x-axis in a sequnce beginning with the minimum x value and continuing 
             up to the maximum x value.}
\item{plots.by.med}{A \code{\link{character}} that give the names of the set of 
             \code{\link{factor}}s by which medians deviations data is to be 
             grouped and a separate plot produced for each group. If \code{NULL}, 
             no groups are formed. If a set of \code{\link{factor}}s, such as  
             \code{Type}, \code{Tuning} and \code{Method}, that uniquely index the 
             combinations of the smoothing-parameter values is specified, then 
             groups are formed for each combination of the levels of the these 
             \code{\link{factor}}s, and and a separate plot is produced for each combination.}
 \item{plots.group.med}{A \code{\link{character}} that gives the names of the set 
             of \code{\link{factor}}s by which the subset of medians deviations data 
             within a single facet in a single plot is to be grouped for plotting  
             as separate lines.}
  \item{facet.x.med}{A \code{\link{character}} giving the \code{\link{factor}}s 
             to be used to form subsets to be plotted in separate columns of the 
             medians deviations plots. The default of \code{"."} results 
             in no split into columns.}
 \item{facet.y.med}{A \code{\link{character}} giving the \code{\link{factor}}s 
             to be used to form subsets to be plotted in separate rows of the 
             medians deviations plots. The default of \code{"."} results 
             in no split into rows.}
 \item{colour.values.med}{A \code{\link{character}} vector specifying the values of 
              the colours to use in drawing the lines for the medians deviations
              within a facet. If this is a named vector, then the values will be 
              matched based on the names. If unnamed, values will be matched in order 
              (usually alphabetical) within the limits of the scale.}
 \item{shape.values.med}{A \code{\link{numeric}} vector specifying the values of 
              the shapes to use in drawing the points for the medians deviations
              within a facet. If this is a named vector, then the values will be 
              matched based on the names. If unnamed, values will be matched in 
              order. If \code{NULL}, as many of the values 
              \code{c(21:24,7,9,10,11,3,4)} as are needed will be used. If more 
              than 10 values are needed, \code{shape.values.med} will need to be 
              set to to a \code{\link{numeric}} of the required length.}
 \item{alpha.med}{A \code{\link{numeric}} specifying the degrees of 
            transparency to be used in plotting a median deviations plot. 
            It is a ratio in which the denominator specifies the number of 
            points (or lines) that must be overplotted to give a solid cover.}
 \item{propn.note.med}{A \code{\link{logical}} indicating whether a note giving the 
             proportion of the median value of the \code{response} for each time is 
             to be included in the \code{medians.deviations} plots.}
 \item{propn.types.med}{A \code{\link{numeric}} giving, for each of the 
             \code{trait.types}, the proportion of the median value of the 
             \code{response} for each time to be used to plot envelopes in the 
             median deviations plots . If set to \code{NULL}, the plots of the 
             proprotion envelopes are omitted.}
 \item{ggplotFuncsMedDevn}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \code{\link{ggplot}} function. 
                It is created by calling the \code{\link{list}} function with 
                a \code{\link{ggplot}} function call for each element.
                These functions are applied in creating the \code{ggplot} object.}
 \item{printPlot}{A \code{\link{logical}} indicating whether or not to print the 
                  plot.}
 \item{...}{allows passing of arguments to \code{\link{plotProfiles}}.}
}
\value{A \code{\link{list}} that consists of two components:  (i) a componenent named  
       \code{plots} that stores a two-level \code{\link{list}} of the median deviations 
       plots; the first-level \code{\link{list}} has a component for each 
       \code{trait.types} and each of these \code{\link{list}}(s) is a second-level 
       \code{\link{list}} that contains the set of plots specified by \code{plots.by.med} 
       (if \code{plots.by.med} is \code{NULL}, a single plot is stored); (ii) a component named 
       \code{med.dev.dat} that stores the \code{\link{data.frame}} containing the 
       median deviations that have been plotted.  Each plot in the \code{plots} 
       \code{\link{list}} is in an object of class \code{\link{ggplot}}, which can be 
       plotted using \code{print}.}
\author{Chris Brien}
\seealso{\code{\link{probeSmooths}}, \code{\link{plotSmoothsComparison}}, \code{\link{plotDeviationsBoxes}}, \code{\link{ggplot}}.}
\examples{
data(exampleData)
vline <- list(ggplot2::geom_vline(xintercept=29, linetype="longdash", size=1))
traits <- probeSmooths(data = longi.dat, 
                       response = "PSA", response.smoothed = "sPSA", 
                       times = "DAP", df = c(4:7), 
                       which.plots = "none", 
                       propn.types.med = NULL)
med <- plotSmoothsMedianDevns(data = traits, 
                              response = "PSA", response.smoothed = "sPSA", 
                              times = "DAP", 
                              plots.group.med = "Tuning",  
                              trait.types = "response", propn.types.med = 0.02,
                              ggplotFuncsMedDevn = vline)
}
\keyword{hplot}
\keyword{manip}