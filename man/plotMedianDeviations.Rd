\name{plotMedianDeviations}
\alias{plotMedianDeviations}
\title{Calculates and plots  the median of the deviations of the smoothed values from the observed values.}
\description{Calculates and plots the median of the deviations of the supplied smoothed values from the 
             supplied observed values for traits and combinations of different smoothing methods and 
             smoothing degrees of freedom, possibly for subsets of \code{\link{factor}} combinations. 
             The requisite values can be generated using \code{\link{probeSmoothing}} with \code{which.plots} 
             set to \code{none}. The results of smoothing methods applied externally to 
             \code{growthPheno} can be included via the \code{extra.smooths} argument. Envelopes of the 
             median value of a trait for each \code{\link{factor}} combination can be added.

             \strong{Note: this function is soft deprecated and may be removed in 
             future versions. \cr Use \code{\link{plotSmoothsMedianDevns}}.}
}
\usage{
plotMedianDeviations(data, response, response.smoothed, 
                     x = NULL, xname="xDays", 
                     individuals = "Snapshot.ID.Tag",  
                     x.title = NULL, y.titles = NULL,
                     facet.x = "Treatment.1", facet.y = "Smarthouse", 
                     labeller = NULL, 
                     trait.types = c("response", "AGR", "RGR"), 
                     propn.types = c(0.1, 0.5, 0.75), propn.note = TRUE, 
                     alpha.med.devn = 0.5, 
                     smoothing.methods = "direct", df, extra.smooths = NULL, 
                     ggplotFuncsMedDevn = NULL, printPlot = TRUE, ...)
}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the observed and smoothed values 
             from which the deviations are to be computed. There should be a column 
             of smoothed values for each combination of \code{smoothing.methods}, 
             \code{df} and the types specified by \code{trait.types}. In addition, 
             there should be a column of values for each element of 
             \code{extra.smooths} in combination with the elements of 
             \code{trait.types}. Also, there should be a column of observed values 
             for each of the types specified by \code{trait.types}}. The naming of the 
             columns for smoothed traits should follow the convention that a name is 
             made up, in the order sepcified, of (i) a \code{response.smoothed}, 
             (ii) the \code{trait.type} if not just a \code{response} trait type, a 
             \code{smoothing.method} or an \code{extra.smooths} and, 
             (iii) if a \code{smoothing.method}, a df.  
             Each component should be separated by a period (.).
 \item{response}{A \code{\link{character}} specifying the response variable for which the 
             observed values are supplied. Depending on the setting of 
             \code{trait.types}, the observed values of related \code{trait.types} 
             may also need to be be supplied.}
 \item{response.smoothed}{A \code{\link{character}} specifying the name of the column 
             containing the values of the smoothed response variable, corresponding 
             to \code{response} and obtained for the combinations of 
             \code{smoothing.methods} and \code{df}, usually using smoothing splines. 
             If \code{response.smoothed} is \code{NULL}, then 
             \code{response.smoothed} is set to the \code{response} to which 
             \code{.smooth} is added. 
             Depending on the setting of \code{trait.types}, the smoothed 
             values of related \code{trait.types} may also need to be be supplied.}
 \item{x}{A \code{\link{character}} giving the variable to be plotted on the 
          x-axis; it may incorporate an expression. If \code{x} is \code{NULL} then 
          \code{xname} is used.}
 \item{xname}{A \code{\link{character}} giving the name of the
             \code{\link{numeric}} that contains the values from which
             \code{x} is derived, it being that \code{x} may incorporate an 
             expression.}
 \item{individuals}{A \code{\link{character}} giving the name of the 
            \code{\link{factor}} that defines the subsets of the \code{data} 
            for which each subset corresponds to the \code{response} values for 
            an individual (plant/cart/plot/unit).}
 \item{x.title}{A \code{\link{character}} giving the title for the x-axis. 
             If \code{NULL} then set to \code{xname}.}
 \item{y.titles}{A \code{\link{character}} giving the titles for the y-axis, 
             one for each trait specified by \code{trait.types}. 
             If \code{NULL} then set to the traits derived for \code{response} 
             from \code{trait.types}.}
 \item{facet.x}{A \code{\link{data.frame}} giving the variable to be used to 
                form subsets to be plotted in separate columns of plots. 
                Use \code{"."} if a split into columns is not wanted. For 
                \code{which.plots} set to \code{methodcompare} or \code{dfcompare} 
                \code{facet.x} is ignored.}
 \item{facet.y}{A \code{\link{data.frame}} giving the variable to be used to 
                form subsets to be plotted in separate rows of plots.
                Use \code{"."} if a split into columns is not wanted.}
 \item{labeller}{A \code{\link{ggplot}} \code{\link{function}} for labelling the 
                 facets of a plot produced using the \code{\link{ggplot}} function. 
                 For more information see \code{\link{ggplot}}.}
 \item{trait.types}{A \code{\link{character}} giving the traits types that are to be 
             plotted. While \code{AGR} and \code{RGR} are commonly used, the names 
             can be arbitrary, except that \code{response} is a special case that 
             indicates that the original response is to be plotted.}
 \item{propn.types}{A \code{\link{numeric}} giving the proportion of the medians 
             values of each of the \code{trait.types} that are to be plotted 
             in the median deviations plots. If set to \code{NULL}, the plots of 
             the proprotions are omitted.}
 \item{propn.note}{A \code{\link{logical}} indicating whether a note giving the 
             proportion of the median values plotted in the \code{compare.medians} 
             plots.}
 \item{alpha.med.devn}{A \code{\link{numeric}} specifying the degrees of 
              transparency to be used in plotting a median deviations plot. 
              It is a ratio in which the denominator specifies the number of 
              points (or lines) that must be overplotted to give a solid cover.}
 \item{smoothing.methods}{A \code{\link{character}} giving the smoothing method used in 
            producing the \code{response.smoothed} and which is to be used in labelling 
            the plot.}
 \item{df}{A \code{\link{numeric}} specifying the smoothing degrees of freedom used in 
            producing the \code{response.smoothed} and which is to be used in labelling 
            the plot.}
 \item{extra.smooths}{A \code{\link{character}} specifying one or more 
            \code{smoothing.method} labels that have been used in naming of columns of 
            smooths of the \code{response} obtained by methods other than the smoothing 
            spline methods provided by \code{growthPheno}. Depending on the setting of 
            \code{trait.types}, the smoothed values of related trait types must also be 
            supplied, with names constructed according to the convention described under 
            \code{data}.}
 \item{ggplotFuncsMedDevn}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \code{\link{ggplot}} function. 
                It is created by calling the \code{\link{list}} function with 
                a \code{\link{ggplot}} function call for each element.
                These functions are applied in creating the \code{ggplot} object.}
 \item{printPlot}{A \code{\link{logical}} indicating whether or not to print any 
                  plots.}
 \item{...}{allows passing of arguments to \code{\link{plotLongitudinal}}.}
}
\value{A \code{\link{list}} that consists of two components:  (i) a componenent named  
       \code{plots} that stores a \code{\link{list}} of the median deviations plots, 
       one for each \code{trait.types}; (ii) a component named \code{med.dev.dat} that stores the 
       \code{\link{data.frame}} containing the median deviations that have been plotted.  
       Each plot in the \code{plots} \code{\link{list}} is in an object of class 
       "\code{\link{ggplot}}", which can be plotted using \code{print}.}
\author{Chris Brien}
\seealso{\code{\link{plotDeviationsBoxes}}, \code{\link{probeSmoothing}}, \code{\link{ggplot}}.}
\examples{
data(exampleData)
vline <- list(ggplot2::geom_vline(xintercept=29, linetype="longdash", size=1),
              ggplot2::scale_x_continuous(breaks=seq(28, 42, by=2)))
traits <- probeSmoothing(data = longi.dat, xname = "xDAP", times.factor = "DAP", 
                         response = "PSA", response.smoothed = "sPSA", 
                         df = c(4:7), 
                         facet.x = ".", facet.y = ".",
                         which.plots = "none",
                         propn.types = NULL)
med <- plotMedianDeviations(data = traits, 
                            response = "PSA", response.smoothed = "sPSA", 
                            x="xDAP", xname = "xDAP", 
                            df = c(4,7), x.title = "DAP", 
                            facet.x = ".", facet.y = ".",
                            trait.types = "response", propn.types = 0.05,
                            ggplotFuncsMedDevn = vline)
}
\keyword{hplot}
\keyword{manip}