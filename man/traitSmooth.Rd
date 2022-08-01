\name{traitSmooth}
\alias{traitSmooth}
\title{Obtain smooths for a trait by fitting spline functions and, having compared 
       several smooths, allows one of them to be chosen and returned in a \code{data.frame}.}
\description{Takes a \code{response} that has been observed for a set of 
             \code{individuals} over a number \code{times} and produces 
             \code{response.smoothed}, using \code{\link{probeSmooths}}, for a default 
             set of smoothing parameter settings (see Details below for the default set of 
             parameters). The settings can be varied from the defaults by 
             specifying alternate values for the smoothing parameters, the parameters 
             being the type of spline (\code{spline.types}), the degrees of freedom 
             (\code{df}) or smoothing penalty (\code{lambdas}) and \code{smoothing.methods} 
             (for details see \code{\link{probeSmooths}}). The secondary traits of 
             the absolute growth rate (\code{AGR}) and relative growth rate 
             (\code{RGR}) are calculated from the two primary traits, the 
             \code{response} and \code{response.smoothed}.

             This function is a wrapper function for \code{\link{probeSmooths}} that 
             (i) uses preset values for many of the \code{\link{probeSmooths}} arguments 
             and  (ii) extends \code{\link{probeSmooths}} to return a single smooth and 
             asociated profile plots. In order to adapt the call to \code{traitSmooth} to 
             the requirements of particular data sets, it can include any of 
             the \code{\link{probeSmooths}} arguments. The \code{\link{probeSmooths}} 
             arguments can be used to modify the Data description, Smoothing parameters, 
             General plot features, Profile plots (pf) and Median-deviations (med) plots.  
             Arguments for \code{\link{plotProfiles}} can also be included in the call, 
             although the \code{facet} and \code{colour} arguments for the chosen smooth 
             must be supplied using the \code{traitSmooth} arguments that end in 
             \code{chosen}.

             Three sets of plots are produced for the \code{response.smoothed} and its 
             \code{AGR} and \code{RGR}: (i) a set of profile plots for each trait for 
             each of the combinations of the smoothing parameters, a profile plot 
             having a single line for each individual (e.g. plant, pot, cart, plot 
             or unit); (ii) a set of medians deviations plots for each trait covering 
             the combinations of the smoothing parameters; (ii) optionally, 
             profile plots of the three smoothed traits for a single combination of the 
             smoothing parameters.

             By default, the single smooth for an arbitrarily chosen combination of the 
             smoothing parameters is returned by the function. The smooth for a single 
             combination other than default combination can be nominated for return 
             using the \code{chosen.smooth} argument. This combination must involve 
             only the supplied values of the smoothing parameters. The values for 
             \code{response}, the \code{response.smoothed} and their \code{AGR}s 
             and \code{RGR}s are are added to \code{data}, after any pre-existing  
             columns of these have been removed from \code{data}. Profile plots of 
             the three smoothed traits are produced using \code{\link{plotProfiles}}. 
             However, if \code{chosen.smooth} is \code{NULL}, all of the smooths will be 
             returned, and plots for the single combination of the smoothing parameters 
             will not be produced. 
}
\usage{
traitSmooth(data, response, response.smoothed, individuals, times, 
            get.rates = TRUE, keep.columns = NULL, 
            x.title = NULL, y.titles = NULL, 
            chosen.smooth = list(spline.type = "PS", 
                                 df = NULL, 
                                 lambda = NULL,
                                 smoothing.method = "logarithmic"), 
             facet.x.chosen = ".", facet.y.chosen = ".", 
             labeller.chosen = NULL, breaks.spacing.x.chosen = 2, 
             colour.chosen = "black", colour.column.chosen = NULL, 
             colour.values.chosen = NULL, alpha.chosen = 0.3, 
             addMediansWhiskers.chosen = TRUE, 
             ggplotFuncsChosen = NULL, 
             ...)
}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the data.}
 \item{response}{A \code{\link{character}} specifying the response variable to be 
              smoothed.}
 \item{response.smoothed}{A \code{\link{character}} specifying the name of the column 
              to contain the values of the smoothed response variable, corresponding 
              to \code{response}.}
 \item{individuals}{A \code{\link{character}} giving the name of the 
              \code{\link{factor}} that defines the subsets of the \code{data} 
              for which each subset corresponds to the \code{response} values for 
              an individual (e.g. plant, pot, cart, plot or unit).}
 \item{times}{A \code{\link{character}} giving the name of the \code{\link{numeric}}, 
              or \code{\link{factor}} with numeric levels, that contains the values 
              of the predictor variable to be supplied to \code{\link{smooth.spline}} 
              and to be plotted on the x-axis.}
 \item{get.rates}{A \code{\link{logical}} specifying whether or not the growth 
             rates (AGR and RGR) of the \code{response} and the 
             \code{response.smoothed} are to be computed and stored. By default, for
             \code{get.rates} set to \code{TRUE}, both the \code{AGR} and \code{RGR} 
             of the \code{response} amd \code{response.smoothed} are computed by 
             differencing. This can be changed using the arguments 
             \code{traits.types} and \code{rates.method} from 
             \code{\link{probeSmooths}}.}
 \item{keep.columns}{A \code{\link{character}} vector giving the names of columns 
              from \code{data} that are to be included in the 
              \code{\link{smooths.frame}} that will be returned.}
 \item{x.title}{Title for the x-axis, used for all plots. If \code{NULL} then set 
              to \code{times}.}
 \item{y.titles}{A \code{\link{character}} giving the titles for the y-axis, one 
              for each the response, the AGE and the RGR. They are used for all plots. 
              If \code{NULL} then they are set to the \code{response} and the  
              \code{response} with \code{.AGR} and \code{.RGR} appended.}
 \item{chosen.smooth}{A named \code{\link{list}} with four components named 
              \code{spline.type}, \code{df},  \code{lambda} and \code{smoothing.method} 
              and each with a single value that specifies the smooth. amongst those 
              investigated, that is to be returned. Set \code{chosen.method} to 
              \code{NULL} to have all smooths returned. Otherwise, the value for 
              \code{spline.type} must be either \code{NCSS} (Natural Cubic Smoothing 
              Spline) or \code{PS} (P=spline). The value of \code{smoothing.method} 
              must be either \code{direct} or \code{logarithmic}. If both \code{df} 
              and \code{lambda} in \code{chosen.smooth} are \code{NULL}, then, 
              depending on the settings for \code{spline.type} and 
              \code{smoothng.method}, the value of either \code{df} or \code{lambdas} 
              that is the median value or the observed value immediatly below the 
              median value will be added to \code{chosen.smooth}. Otherwise, one of 
              \code{df} and \code{lambda} should be \code{NULL} and the other should 
              be a single \code{numeric} value. The four values in \code{chosen.smooth} 
              must be amongst those for which the smooths have been generated. 
              If a value in \code{chosen.smooth} is not amongst those investigated, 
              a value that was investigated will be substituted.}
 \item{facet.x.chosen}{A \code{\link{character}} giving the names of the \code{\link{factor}}s 
              to be used to form subsets to be plotted in separate columns of the 
              profile plots for the chosen smooth. The default of \code{"."} results 
              in no split into columns.}
 \item{facet.y.chosen}{A \code{\link{character}} giving the \code{\link{factor}}s 
              to be used to form subsets to be plotted in separate rows of the 
              profile plots for the chosen smooth. The default of \code{"."} results 
              in no split into rows.}
 \item{labeller.chosen}{A \code{\link{ggplot}} \code{\link{function}} for labelling the 
              facets of plots produced using the \code{\link{ggplot}} function, 
              including profile plots for the chosen smooth and plots produced using 
              \code{\link{probeSmooths}}. For more information see \code{\link{ggplot}}.}
 \item{breaks.spacing.x.chosen}{A \code{\link{numeric}} specifying the distance between 
             breaks for the x-axis in a sequnce beginning with the minimum x value and 
             continuing up to the maximum x value.}
 \item{colour.chosen}{A \code{\link{character}} specifying a single colour to use in 
              drawing the lines for the profiles. If colouring according to the 
              values of a variable is required then use \code{colour.column.chosen}.}
 \item{colour.column.chosen}{A \code{\link{character}} giving the name of a column 
              in \code{data} over whose values the colours of the lines for the profiles 
              are to be varied. The colours can be specified using 
              \code{colour.values.chosen}.}
 \item{colour.values.chosen}{A \code{\link{character}} vector specifying the values of 
              the colours to use in drawing the lines for the profiles.  
              If this is a named vector, then the values will be matched based 
              on the names. If unnamed, values will be matched in order 
              (usually alphabetical) within the limits of the scale.}
 \item{alpha.chosen}{A \code{\link{numeric}} specifying the degrees of transparency to 
              be used in plotting. It is a ratio in which the denominator 
              specifies the number of points (or lines) that must be overplotted 
              to give a solid cover.}
 \item{addMediansWhiskers.chosen}{A \code{\link{logical}} indicating whether plots over 
             time of the medians and outer whiskers are to be added to the profile plot 
             for the chosen smooth. The outer whiskers are related to the whiskers on 
             a box-and-whisker and are defined as the median plus (and minus) 1.5 times 
             the interquartile range (IQR). Points lying outside the whiskers are 
             considered to be potential outliers.}
 \item{ggplotFuncsChosen}{A \code{\link{list}}, each component of which contains the 
              results of evaluating a \code{\link{ggplot}} function. 
              It is created by calling the \code{\link{list}} function with 
              a \code{\link{ggplot}} function call for each component. Note that 
              these functions are applied to the profile plots for the chosen smooth.}
 \item{...}{allows arguments to be passed to \code{\link{probeSmooths}} and 
       to \code{\link{plotProfiles}}.}
}
\details{
             This function is a wrapper function for \code{\link{probeSmooths}}. 
             The default values for the arguments differ between the two functions, 
             the default values for the \code{traitSmooth} arguments that are different being: 

             \describe{
             \item{smoothing.methods:}{"logarithmic"}
             \item{spline.types:}{c("NCSS","PS") }
             \item{df:}{5:7}
             \item{lambdas:}{list(PS = round(10^c(-0.5, 0, 0.5, 1), digits = 3))}
             \item{which.plots:}{c("traits", "medians.deviations")}
             \item{plots.by.pf:}{"Type"}
             \item{facet.x.pf:}{c("Method", "Tuning")}
             \item{include.raw.pf:}{"facet.x"}
             \item{facet.x.med:}{c("Method","Type")} 
             \item{plots.group.med:}{"Tuning"}
             \item{propn.types.med:}{NULL}
             \item{printPlot:}{TRUE}
}
}
\value{A \code{\link{smooths.frame}} or a \code{\link{data.frame}} that contains 
       the unsmoothed and smoothed data in long format. That is, all the values for 
       either an unsmoothed or a smoothed trait are in a single column. 
       
       A \code{\link{smooths.frame}} will be returned when \code{chosen.smooth} is 
       \code{NULL}. It will contains the smooths for a trait for the different 
       combinatons of the smoothing parameters, the values for the different smooths 
       being placed in rows one below the other. The columns that are included in the 
       \code{\link{smooths.frame}} are \code{Type}, \code{TunePar}, \code{TuneVal}, 
       \code{Tuning} and \code{Method}, as well as those specified by 
       \code{individuals}, \code{times}, \code{response}, and \code{response.smoothed}, 
       and any included in the \code{keep.columns}, \code{plots} and \code{facet} 
       arguments. The \code{AGR} or \code{RGR} for the \code{response} and 
       \code{response.smoothed} will also be included. A \code{\link{smooths.frame}} 
       has the attributes described in \code{\link{smooths.frame}}. 

       A \code{\link{data.frame}} will be returned when \code{chosen.smooth} is 
       not \code{NULL}. It will be the same as \code{data}, but with \code{response}, 
       \code{response.smoothed} for the chosen smooth and their AGRs and RGRs 
       added, the added columns replacing any pre-existing columns of the same name 
       in \code{data}.

       The \code{\link{smooths.frame}}/\code{\link{data.frame}} is returned invisibly.}
\author{Chris Brien}
\seealso{\code{\link{probeSmooths}}.}
\examples{
data(exampleData)
vline <- list(ggplot2::geom_vline(xintercept=29, linetype="longdash", size=1))
yfacets <- c("Smarthouse", "Treatment.1")
smth.dat <- traitSmooth(data = longi.dat, 
                        response = "PSA", response.smoothed = "sPSA",
                        individuals = "Snapshot.ID.Tag",times = "DAP", 
                        keep.columns = yfacets, 
                        facet.y.pf = yfacets, 
                        facet.y.chosen = yfacets, 
                        ggplotFuncsProfile=vline, ggplotFuncsChosen = vline)

}\keyword{hplot}
\keyword{manip}