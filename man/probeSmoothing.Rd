\name{probeSmoothing}
\alias{probeSmoothing}
\title{Compares, for a set of specified values of df and different smoothing methods, 
       a response and the smooths of it, possibly along with growth rates calculated 
       from the smooths}
\description{Takes a \code{response} and, for each \code{individual}, uses 
             \code{\link{splitSplines}} to smooth its values for each individual 
             using  the degrees of freedom values in \code{df}. 
             Provided \code{get.rates} is \code{TRUE},
             both the Absolute Growth Rates (AGR) and the Relative Growth Rates (RGR) 
             are calculated for each smooth, either using differences or first 
             derivatives. A combination of the unsmoothed and smoothed values, 
             as well as the AGR and RGR, can be plotted for each value in 
             \code{smoothing methods} in combination with \code{df}. Note that the 
             arguments that modify the plots apply to all 
             plots that are produced. The handling of missing values is controlled 
             via \code{na.x.action} and \code{na.y.action}}
\usage{
probeSmoothing(data, response = "Area", x = NULL, xname="xDays", 
               times.factor = "Days", individuals="Snapshot.ID.Tag", 
               na.x.action="exclude", na.y.action = "exclude", 
               df, smoothing.methods = "direct", correctBoundaries = FALSE, 
               get.rates = TRUE, rates.method="differences", 
               facet.x = "Treatment.1", facet.y =   "Smarthouse", 
               labeller = NULL, x.title = NULL, 
               colour = "black", colour.column=NULL, 
               colour.values=NULL, alpha = 0.1, 
               trait.types = c("response", "AGR", "RGR"), 
               propn.types = c(0.1, 0.5, 0.75), propn.note = TRUE, 
               which.plots = "smoothedonly", 
               deviations.plots = "none", alpha.med.devn = 0.5,
               ggplotFuncs = NULL, ggplotFuncsMedDevn = NULL, 
               ...)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the data.}
 \item{response}{A \code{character} specifying the response variable to be 
            supplied to \code{\link{smooth.spline}} and that is to be plotted 
            on the y-axis.}
  \item{x}{A \code{\link{character}} giving the variable to be plotted on the 
            x-axis; it may incorporate an expression, it being that \code{x} may 
            incorporate an expression. If \code{x} is \code{NULL} then 
            \code{xname} is used.}
 \item{xname}{A \code{\link{character}} giving the name of the
            \code{\link{numeric}} that contains the values of the predictor 
            variable to be supplied to \code{\link{smooth.spline}} and 
            from which \code{x} is derived.}
 \item{times.factor}{A \code{\link{character}} giving the name of the column in 
            \code{data} containing the factor for times at which the data was 
            collected. Its levels will be used in calculating growth rates and 
            should be numeric values stored as characters.}
 \item{individuals}{A \code{\link{character}} giving the name(s) of the 
            \code{\link{factor}}(s) that define the subsets of the \code{data} 
            for which each subset corresponds to the \code{response} values for 
            an individual.}
 \item{na.x.action}{A \code{character} string that specifies the action to 
            be taken when values of \code{x} are \code{NA}. The possible 
            values are \code{fail}, \code{exclude} or \code{omit}. 
            For \code{exclude} and \code{omit}, predictions and derivatives 
            will only be obtained for nonmissing values of \code{x}. 
            The difference between these two codes is that for \code{exclude} the returned 
            \code{data.frame} will have as many rows as \code{data}, the 
            missing values have been incorporated. }
 \item{na.y.action}{A \code{character} string that specifies the action to 
            be taken when values of \code{y}, or the \code{response}, are 
            \code{NA}.  The possible values are \code{fail}, \code{exclude}, 
            \code{omit}, \code{allx}, \code{trimx}, \code{ltrimx} or 
            \code{rtrimx}. For  all options, except \code{fail}, missing 
            values in \code{y} will be removed before smoothing. 
            For \code{exclude} and \code{omit}, predictions 
            and derivatives will be obtained only for nonmissing values of 
            \code{x} that do not have missing \code{y} values. Again, the 
            difference between these two is that, only for \code{exclude} 
            will the missing values be incorporated into the 
            returned \code{data.frame}. For \code{allx}, predictions and 
            derivatives will be obtained for all nonmissing \code{x}. 
            For \code{trimx}, they will be obtained for all nonmissing 
            \code{x} between the first and last nonmissing \code{y} values 
            that have been ordered for \code{x}; for \code{ltrimx} and 
            \code{utrimx} either the lower or upper missing \code{y} 
            values, respectively, are trimmed.}
 \item{df}{A \code{\link{numeric}} specifying the set of degrees of freedom to 
              be probed.}
 \item{smoothing.methods}{A \code{\link{character}} giving one or more methods 
            to use for smoothing. Currently, the two possibilites 
            are (i) \code{"direct"}, for directly smoothing the observed 
            \code{response}, and (ii) \code{"logarithmic"}, for smoothing the 
            \code{log}-transformed \code{response}, followed by taking 
            exponentials of the fitted values to back-trasnform them.}
 \item{correctBoundaries}{A \code{\link{logical}} indicating whether the fitted 
            spline values are to have the method of Huang (2001) applied 
            to them to correct for estimation bias at the end-points. Note that 
            if \code{rates.method} is set to \code{"derivatives"} then it is not 
            possible to have \code{correctBoundaries} set to \code{TRUE}.}
 \item{get.rates}{A \code{\link{logical}} specifying whether or not the growth 
              rates (AGR and RGR) are to be computed and stored.}
 \item{rates.method}{A \code{\link{character}} specifying the method to use in 
               calculating the growth rates. The two possibilities are 
               \code{"differences"} and \code{"derivates"}.}
 \item{facet.x}{A \code{\link{data.frame}} giving the variable to be used to 
                form subsets to be plotted in separate columns of plots. 
                Use \code{"."} if a split into columns is not wanted. For 
                \code{which.plots} set to \code{methodscompare} or 
                \code{dfcompare} , \code{facet.x} is ignored.}
 \item{facet.y}{A \code{\link{data.frame}} giving the variable to be used to 
                form subsets to be plotted in separate rows of plots.
                Use \code{"."} if a split into columns is not wanted.}
 \item{labeller}{A \code{\link{ggplot}} \code{\link{function}} for labelling the 
                 facets of a plot produced using the \code{\link{ggplot}} function. 
                 For more information see \code{\link{ggplot}}.}
 \item{x.title}{Title for the x-axis. If \code{NULL} then set to \code{times.factor}.}
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
 \item{trait.types}{A \code{\link{character}} giving the \code{trait.types} that 
             are to be produced. One of more of \code{response}, \code{AGR} and 
             \code{RGR}. If \code{all}, all three traits are produced. 
             The unsmoothed growth rates are not calculated if only smoothed plots 
            are requested.}
 \item{propn.types}{A \code{\link{numeric}} giving the proportion of the median 
             values of each of the \code{trait.types} that are to be plotted in 
             the \code{compare.medians} plots of the deviations of the observed 
             values from the smoothed values. If set 
             to \code{NULL}, the plots of the proprotions of the median values of
             the traits are omitted.}
 \item{propn.note}{A \code{\link{logical}} indicating whether a note giving the 
             proportion of the median values plotted in the \code{compare.medians} 
             plots.}
 \item{which.plots}{A \code{\link{character}} giving the plots that are to be 
             produced. If \code{none}, no plots are produced. If \code{smoothedonly}, 
             plots of the smoothed traits are plotted. If \code{bothseparately}, 
             plots of the unsmoothed trait followed by the smoothed traits are 
             produced for each trait. If \code{methodscompare}, a combined plot of 
             the smoothed traits for each \code{smoothing.methods} is produced, 
             for each value of \code{df}. If \code{methods+rawcompare}, the unsmoothed 
             trait is added to the combined plot. if \code{dfcompare}, a combined 
             plot of the smoothed trait for each \code{df} is produced, for each
             \code{smoothing.methods}. If \code{df+rawcompare}, the unsmoothed 
             trait is added to the combined plot.}
 \item{deviations.plots}{A \code{\link{character}} is either \code{none} or any 
             combination of \code{absolute.boxplots}, \code{relative.boxplots} and 
             \code{compare.medians}. If \code{none}, no plots are produced. 
             Boxplots of the absolute and relative deviations are specified by 
             \code{absolute.boxplots} and  \code{relative.boxplots}. The absolute 
             deviations are the values of a trait minus their smoothed values 
             (observed - smoothed); the relative deviations are the absolute 
             deviations divided by the smoothed values of the trait. The option 
             \code{compare.medians} results in a plot that compares the medians 
             of the deviations over the \code{times.factor} for each combination 
             of the \code{smoothing.methods} and the \code{df}. The argument 
             \code{trait.types} controls the traits for which boxplots are produced.}
 \item{alpha.med.devn}{A \code{\link{numeric}} specifying the degrees of 
              transparency to be used in plotting a median deviations plot. 
              It is a ratio in which the denominator specifies the number of 
              points (or lines) that must be overplotted to give a solid cover.}
 \item{ggplotFuncs}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \code{\link{ggplot}} function. 
                It is created by calling the \code{\link{list}} function with 
                a \code{\link{ggplot}} function call for each element. Note that 
                these functions are applied to all three plots produced.}
 \item{ggplotFuncsMedDevn}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \code{\link{ggplot}} function. 
                It is created by calling the \code{\link{list}} function with 
                a \code{\link{ggplot}} function call for each element. Note that 
                these functions are applied to the \code{compare.median} deviations plots only.}
\item{...}{allows passing of arguments to \code{\link{plotLongitudinal}}.}
}
\value{A \code{\link{data.frame}} containing \code{individuals}, 
       \code{times.factor}, \code{facet.x}, \code{facet.y}, \code{xname}, 
       \code{response}, and, for each \code{df}, the smoothed 
       response, the AGR and the RGR. It is returned invisibly. The names of 
       the new data are constructed by joining elements separated by full 
       stops (\code{.}). In all cases, the last element is the value of 
       \code{df}. For the smoothed response, the other elements are 
       \code{response} and \code{"smooth"}; for AGR and RGR, the other elements 
       are the name of the smoothed response and either \code{"AGR"} or 
       \code{"RGR"}.}
\author{Chris Brien}
\seealso{\code{\link{splitSplines}}, \code{\link{splitContGRdiff}}, \code{\link{smooth.spline}}, \code{\link{ggplot}}.}
\examples{
data(exampleData)
vline <- list(ggplot2::geom_vline(xintercept=20, linetype="longdash", size=1),
              ggplot2::scale_x_continuous(breaks=seq(12, 36, by=2)))
probeSmoothing(data = longi.dat, response = "Area", df = c(4,7), x="xDays+24.16666667", 
               ggplotFuncs=vline)
}\keyword{hplot}
\keyword{manip}