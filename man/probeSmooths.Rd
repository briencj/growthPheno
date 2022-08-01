\name{probeSmooths}
\alias{probeSmooths}
\title{For a response in a \code{data.frame} in long format, computes and compares, for 
       sets of smoothing parameters, smooths of the response, possibly along with 
       growth rates calculated from the smooths.}
\description{Takes an observed \code{response} and, for each \code{individual}, uses 
             \code{byIndv4Times_SplinesGRs} to smooth its values employing the smoothing 
             parameters specified by (i) \code{spline.types}, (ii) the tuning 
             parameters, being the degrees of freedom values in \code{df} or the 
             smoothing penalties in \code{lambdas}, and (iii) the 
             \code{smoothing.methods}. An alternative to specifying the values of the 
             smoothing parameters is to supply a customized set of the 
             smoothing-parameter values in a \code{\link{data.frame}} that has a column 
             for each smoothing parameter (see the \code{smoothing.schemes} argument). 
             All combinations of the smoothing parameters that make up the rows of the 
             \code{smoothing.schemes} \code{\link{data.frame}} will be compared.

             Provided \code{get.rates} is \code{TRUE} and depending on the setting 
             of \code{trait.types}, the Absolute Growth Rates (AGR) and/or the 
             Relative Growth Rates (RGR) are calculated  for each individual from the 
             unsmoothed, observed \code{response} using \code{differences} and from 
             the smooths of the \code{response}, using either differences 
             or first derivatives, as specified by \code{rates.method}. 

             Generally, profile plots for the traits (a response, an AGR or an RGR) 
             specified in \code{traits.types} are produced if \code{which.plots} is 
             \code{profiles}; if \code{which.plots} specifies one or more deviations 
             plots, then those deviations plots, based on the unsmoothed data from 
             which the smoothed data has been subtracted, will also be produced. 
             The layout of the plots is controlled via combinations of one or more 
             of the smoothing-parameter \code{\link{factor}}s \code{Type}, 
             \code{TunePar}, \code{TuneVal}, \code{Tuning} (the combination 
             of \code{TunePar} and \code{TuneVal}) and \code{Method}, as well as 
             other \code{\link{factor}}s associated with the \code{data}. 
             The \code{\link{factor}}s that are to be used are supplied to the arguments 
             \code{plots.by.pf}, \code{facet.x.pf}, and \code{facet.y.pf} for the 
             profile plots and deviations boxplots; for the plots of the medians of the 
             deviations, the \code{\link{factor}}s are supplied using the  arguments 
             \code{plots.by.med}, \code{facet.x.med}, \code{facet.y.med} and 
             \code{plots.group.med}. Here, the basic principle is that the number of 
             levels combinations of the smoothing-parameter \code{\link{factor}}s 
             included in the set of \code{plots} and \code{facets} arguments with the 
             same suffix (\code{pf} for profile plots and \code{med} for medians 
             deviations plots) must be the same as those covered by the combinations of 
             the values supplied to \code{spline.types}, \code{df}, \code{lambdas} and 
             \code{smoothing.method} and incorporated into the \code{\link{smooths.frame}}, 
             such as is returned by \code{probeSmooths}. This ensures that smooths 
             from different parameter sets are not pooled together in a single plot. 
             It is also possible to include \code{\link{factor}}s that are not 
             smoothing-parameter \code{\link{factor}}s in the \code{plots} amd 
              \code{facets} arguments. 

             The following profiles plots can be produced: (i) separate plots of the 
             smoothed traits for each combination of the smoothing parameters 
             (include \code{Type}, \code{Tuning} and \code{Method} in \code{plots.by.pf}); 
             (ii) as for (i), with the corresponding plot for the unsmoothed trait 
             preceeding the plots for the smoothed trait (also set \code{include.raw.pf} to 
             \code{alone}); (iii) profiles plots that compare a smoothed trait for all 
             combinations of the values of the smoothing parameters, arranging the plots 
             side-by-side or one above the other (include \code{Type}, \code{Tuning} and 
             \code{Method} in \code{facet.x.pf} and/or \code{facet.y.pf} - to include the 
             unsmoothed trait set \code{include.raw.pf} to one of \code{facet.x} or 
             \code{facet.y}; (iv) as for (iii), except that separate plots are 
             produced for each combination of the levels of the \code{\link{factor}}s 
             in \code{plot.by} and each plot compares the smoothed traits for the 
             smoothing-parameter \code{\link{factor}}s included in \code{facet.x.pf} 
             and/or \code{facet.y.pf} (set both \code{plots.by.pf} and one or more of 
             \code{facet.x.pf} and \code{facet.y.pf}). 

             Deviation plots that can be produced are the absolute and relative deviations 
             boxplots and plots of medians deviations (see \code{which.plots}). 

             The handling of missing values is controlled 
             via \code{na.x.action} and \code{na.y.action}.

             The \code{probeSmooths} arguments are grouped according to function in the following order:
             \enumerate{
                \item{\bold{Data description arguments:} data, response, response.smoothed, 
                      individuals, times, keep.columns, trait.types, get.rates, 
                      rates.method, ntimes2span.}
                \item{\bold{Smoothing parameters:} smoothing.methods, smoothing.segments, 
                      spline.types, df, lambdas, npspline.segments, smoothing.schemes, 
                      na.x.action, na.y.action, external.smooths.}
                \item{\bold{General plot control:} x.title, y.titles, labeller, which.plots.} 
                \item{\bold{Profile plots (pf) features:} plots.by.pf, facet.x.pf, facet.y.pf, 
                      collapse.facets.x.pf, collapse.facets.y.pf, include.raw.pf, colour.pf, 
                      colour.column.pf, colour.values.pf, alpha.pf, addMediansWhiskers.pf, 
                      ggplotFuncsProfile.} 
                \item{\bold{Median-deviations (med) plots features:} plots.by.med, 
                      plots.group.med, facet.x.med, facet.y.med,
                      colour.values.med, shape.values.med, alpha.med, 
                      propn.note.med, propn.types.med, ggplotFuncsMedDevn.}
             }
}
\usage{
probeSmooths(data, response = "PSA", response.smoothed = NULL, 
             individuals="Snapshot.ID.Tag", times = "DAP", 
             keep.columns = NULL, 
             trait.types = c("response", "AGR", "RGR"), 
             get.rates = TRUE, rates.method="differences", 
             ntimes2span = NULL, 
             smoothing.methods = "direct", smoothing.segments = NULL, 
             spline.types = "NCSS", df = NULL, lambdas = NULL, 
             npspline.segments = NULL, smoothing.schemes = NULL, 
             na.x.action = "exclude", na.y.action = "trimx", 
             external.smooths = NULL, 
             correctBoundaries = FALSE, 
             x.title = NULL, y.titles = NULL, labeller = NULL, 
             which.plots = "profiles", 
             breaks.spacing.x = 4, 
             plots.by.pf = NULL, facet.x.pf = ".", facet.y.pf = ".", 
             collapse.facets.x.pf = TRUE, collapse.facets.y.pf = FALSE, 
             include.raw.pf = "no", 
             colour.pf = "black", colour.column.pf = NULL, 
             colour.values.pf = NULL, alpha.pf = 0.3, 
             addMediansWhiskers.pf = TRUE,
             ggplotFuncsProfile = NULL, 
             plots.by.med = NULL, plots.group.med = NULL, 
             facet.x.med = ".", facet.y.med = ".",
             colour.values.med = NULL, shape.values.med = NULL, 
             alpha.med = 0.5, 
             propn.note.med = TRUE, propn.types.med = c(0.1, 0.5, 0.75), 
             ggplotFuncsMedDevn = NULL, 
             ggplotFuncsDevnBoxes = NULL, ...)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the data.}
 \item{response}{A \code{\link{character}} specifying the response variable to be 
            supplied to \code{\link{smoothSpline}} and that is to be plotted 
            on the y-axis.}
 \item{response.smoothed}{A \code{\link{character}} specifying the name of the column 
             containing the values of the smoothed response variable, corresponding 
             to \code{response}. If \code{response.smoothed} is \code{NULL}, then 
             \code{response.smoothed} is set to the \code{response} to which is added 
             the prefix \code{s}.}
 \item{individuals}{A \code{\link{character}} giving the name of the 
            \code{\link{factor}} that defines the subsets of the \code{data} 
            for which each subset corresponds to the \code{response} values for 
            an individual (e.g. plant, pot, cart, plot or unit).}
 \item{times}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the times at which the data was 
             collected, either as a \code{\link{numeric}}, \code{\link{factor}}, or 
             \code{\link{character}}. It will be used as the values of the predictor 
             variable to be supplied to \code{\link{smooth.spline}} and to be plotted 
             on the x-axis. If a \code{\link{factor}} or \code{\link{character}}, 
             the values should be numerics stored as characters.}
\item{keep.columns}{A \code{\link{character}} vector giving the names of columns 
            from \code{data} that are to be included in the 
            \code{\link{smooths.frame}} that will be returned. Its main use is 
            when no plots are being produced by \code{probeSmooths}, but there 
            are columns in the supplied \code{\link{data.frame}} that are likely 
            to be needed for the \code{plots} and \code{facets} arguments when 
            producing plots subsequently.}
 \item{trait.types}{A \code{\link{character}} giving the \code{trait.types} that 
             are to be computed when \code{get.rates} is \code{TRUE}. Irrespective 
             of the setting of \code{get.rates}, the nominated profiles are plotted. 
             If \code{all}, the \code{AGR} and the \code{RGR} of the of a 
             \code{response} and the \code{response.smoothed} are computed and, 
             for each, a \code{response} and its \code{AGR} and \code{RGR} are 
             plotted.}
 \item{get.rates}{A \code{\link{logical}} specifying whether or not the growth 
             rates (AGR and RGR) of the \code{response}, which must be computed 
             by differencing, and the \code{response.smoothed} are to be 
             computed and stored.}
 \item{rates.method}{A \code{\link{character}} specifying the method to use in 
             calculating the growth rates for \code{response.smoothed}. The two 
             possibilities are \code{"differences"} and \code{"derivatives"}.}
 \item{ntimes2span}{A \code{\link{numeric}} giving the number of values in 
             \code{times} to span in calculating growth rates by differencing. 
             For \code{ntimes2span} set to \code{NULL}, if \code{rates.method} 
             is set to \code{differences} then \code{ntimes2span} is set to 2; 
             if \code{rates.method} is set to \code{derivatives} then 
             \code{ntimes2span} is set to 3. 
             Note that when \code{get.rates} is \code{TRUE}, the growth rates 
             for the unsmoothed response must be calculated by differencing, 
             even if the growth rates for the smoothed response are computed 
             using \code{derivatives}. When differencing, 
             each growth rate is calculated as the difference in the values of 
             one of the \code{responses} for pairs of \code{times} values that 
             are spanned by \code{ntimes2span} \code{times} values divided by 
             the difference between this pair of \code{times} values. For 
             \code{ntimes2span} set to 2, a growth rate is the  
             difference between consecutive pairs of values of one of the 
             \code{responses} divided by the difference between consecutive 
             pairs of \code{times} values.}
 \item{smoothing.methods}{A \code{\link{character}} giving the smoothing method 
      to use. The two possibilites are (i) \code{"direct"}, for directly 
      smoothing the observed \code{response}, and (ii) \code{"logarithmic"}, for 
      smoothing the \code{log}-transformed \code{response} and then 
      back-transforming by taking the exponentional of the fitted values.}
 \item{smoothing.segments}{A named \code{\link{list}}, each of whose components 
      is a numeric pair specifying the first and last values of an 
      \code{times}-interval whose data is to be subjected as an entity to smoothing 
      using splines. The separate smooths will be combined to form a whole 
      smooth for each individual. If \code{get.rates} is \code{TRUE}, 
      \code{rates.method} is \code{differences} and \code{ntimes2span} is 2, 
      the smoothed growth rates will be computed over the set of segments; 
      otherwise, they will be computed within segments. 
      If \code{smoothing.segments} is \code{NULL}, the data is not 
      segmented for smoothing.}
 \item{spline.types}{A \code{\link{character}} giving the type of spline  
      to use. Currently, the possibilites are (i) \code{"NCSS"}, for natural 
      cubic smoothing splines, and (ii) \code{"PS"}, for P-splines.}
 \item{df}{A \code{\link{numeric}} with at least one value that specifies, for 
      natural cubic smoothing splines (\code{NCSS}), the desired equivalent 
      numbers of degrees of freedom of the smooths (trace of the smoother matrix). 
      Lower values result in more smoothing. If \code{df = NULL}, the amount of 
      smoothing can be controlled by including a component named \code{NCSS} 
      in the \code{\link{list}} for \code{lambdas}. If \code{df} is \code{NULL} 
      and \code{lambda} does not include a component named \code{NCSS}, 
      then an error is issued.}
 \item{lambdas}{A named \code{\link{list}} or a \code{\link{numeric}} specifying 
      the positive penalties to apply in order to control the amount of smoothing. 
      The amount of smoothing decreases as \code{lamda} decreases. 
      If \code{lambdas} is a \code{\link{list}}, then include a 
      components with \code{lambdas} values and  named for each of the specified 
      values of \code{spline.types} for which \code{lambdas} are to be used. 
      If \code{spline.types} includes \code{PS}, then a component named 
      \code{PS} with at least one numeric value must be present. If a 
      \code{\link{numeric}}, then it will be converted to a \code{\link{list}} 
      with the single component named \code{PS}.}
 \item{npspline.segments}{A \code{\link{numeric}} specifying, for P-splines (\code{PS}), 
      the number of equally spaced segments between \code{min(x)} and \code{max(x)}, 
      excluding missing values, to use in constructing the B-spline basis for the 
      spline fitting. If \code{npspline.segments} is NULL, \code{npspline.segments} 
      is set to the maximum of 10 and \code{ceiling((nrow(data)-1)/2)} i.e. there will 
      be at least 10 segments and, for more than 22 \code{times} values, there will be 
      half as many segments as there are \code{times} values. The amount of smoothing 
      decreases as \code{npspline.segments} increases. When the data has been 
      segmented for smoothing (\code{smoothing.segments} is not \code{NULL}), 
      an \code{npspline.segments} value can be supplied for each segment.}
 \item{smoothing.schemes}{A \code{\link{data.frame}} whose rows each specify smoothing 
      schemes to be compared. The columns of the \code{\link{data.frame}} must be named 
      \code{Type}, \code{TunePar}, \code{TuneVals} and \code{Method}. The values in the 
      \code{Type} and \code{Method} columns must contain at least enough of the legal 
      options for the \code{spline.types} and \code{smoothing.methods} arguments of 
      \code{\link{smoothSpline}} so that they are unique amongst the possible options. 
      The column \code{TunePar} should contain, for \code{Method} set to \code{NCSS}, 
      either \code{"df"} or \code{"lambda"} and, for \code{Method} set to \code{PS}, 
      \code{"lambda"}, or truncated versions of these. The column \code{TuneVals} 
      should contain values appropriate to the corresponding settings of \code{TunePar}. 
      The column \code{Tuning} will be formed from \code{TunePar} and \code{TuneVal} 
      and added to the \code{\link{data.frame}}.}
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
 \item{external.smooths}{A \code{\link{data.frame}} containing the one or more smooths of 
      a \code{response} in the column specified by \code{smoothed.response}. 
      Multiple smoooths should be supplied in long.format with the same columns 
      as the \code{\link{smooths.frame}} \code{data}, except for the 
      smoothing-parameter columns \code{Type}, \code{TunePar}, \code{TuneVal}, 
      \code{Tuning} and \code{Method}. Only those smoothing-parameter columns that 
      are to be used in any of \code{plots.by.pf}, \code{plots.group.med}, \code{facet.x.pf} 
      and \code{facet.x.med} should be included with labels appropriate to the 
      \code{external.smooths}. Those smoothing-parameter columns not included in 
      \code{external.smooths} will have columns of \code{"Other"} added to 
      \code{external.smooths}.}
 \item{correctBoundaries}{A \code{\link{logical}} indicating whether the fitted spline 
             values are to have the method of Huang (2001) applied 
             to them to correct for estimation bias at the end-points. Note that 
             \code{spline.type} must be \code{NCSS} and \code{lambda} and \code{deriv} 
             must be \code{NULL} for \code{correctBoundaries} to be set to \code{TRUE}.}
 \item{x.title}{Title for the x-axis, used for all plots. If \code{NULL} then set to 
             \code{times}.}
 \item{y.titles}{A \code{\link{character}} giving the titles for the y-axis, 
             one for each trait specified by \code{trait.types} and used for all plots. 
             If \code{NULL} then set to the traits derived for \code{response} 
             from \code{trait.types}.}
\item{labeller}{A \code{\link{ggplot}} \code{\link{function}} for labelling the 
             facets of a plot produced using the \code{\link{ggplot}} function. 
             For more information see \code{\link{ggplot}}.}
 \item{which.plots}{A \code{\link{logical}} indicating which plots are to be produced. 
             The options are either \code{none} or some combination of \code{profiles}, 
             \code{absolute.boxplots}, \code{relative.boxplots} and 
             \code{medians.deviations}. The various \code{profiles} plots that can be 
             poduced are described in the introduction to this function. 

             Boxplots of the absolute deviations are specified by 
             \code{absolute.boxplots}, the absolute deviations being the values of a 
             trait minus their smoothed values (observed - smoothed). Boxplots of the 
             relative deviations are specified by \code{relative.boxplots}, the 
             relative deviations being the absolute deviations divided by the smoothed 
             values of the trait. 

             The option \code{medians.deviations} results 
             in a plot that compares the medians of the absolute deviations over the 
             values of \code{times} for each combination of the smoothing-parameter values. 
             The arguments to \code{probeSmooths} that apply to \code{medians.deviations} 
             plots have the suffix \code{med}.}
 \item{breaks.spacing.x}{A \code{\link{numeric}} specifying the distance between breaks for 
             the x-axis in a sequnce beginning with the minimum x value and continuing 
             up to the maximum x value.}
\item{plots.by.pf}{A \code{\link{character}} that gives the names of the set of 
             \code{\link{factor}}s by which the data is to be grouped and a separate 
             plot produced for each group. If \code{NULL}, no groups are formed. If a 
             set of \code{\link{factor}}s, such as  \code{Type}, \code{Tuning} and 
             \code{Method}, that uniquely index the combinations of the 
              smoothing-parameter values is specified, then groups are formed for each 
             combination of the levels of the these \code{\link{factor}}s, and a separate plot is 
             produced for each combination.}
 \item{facet.x.pf}{A \code{\link{character}} giving the names of the \code{\link{factor}}s 
             to be used to form subsets to be plotted in separate columns of the 
             profiles plots and deviations boxplots. The default of \code{"."} results 
             in no split into columns.}
 \item{facet.y.pf}{A \code{\link{character}} giving the \code{\link{factor}}s 
             to be used to form subsets to be plotted in separate rows of the 
             profiles plots and deviations boxplots. The default of \code{"."} results 
             in no split into rows.}
 \item{collapse.facets.x.pf}{A \code{\link{logical}} to indicate whether all variables 
             specified by \code{facets.x.pf} are to be collapsed to a single variable. 
             Note that the smoothing-parameters factors, if present, are always collapsed.}
 \item{collapse.facets.y.pf}{A \code{\link{logical}} to indicate whether all variables 
             specified by \code{facets.y.pf} are to be collapsed to a single variable. 
             Note that the smoothing-parameters factors, if present, are always collapsed.}
 \item{include.raw.pf}{A \code{\link{character}} indicating whether plots of the 
             raw (unsmoothed) trait, corresponding to the plots of the smoothed traits, 
             are to be included in profile plots. The options are \code{no}, \code{alone}, 
             \code{facet.x}, or \code{facet.y}. That is, the plots of the raw traits are 
             plotted separately or as part of either \code{facet.x.pf} or 
             \code{facet.y.pf}.}
 \item{colour.pf}{A \code{\link{character}} specifying a single colour to use in 
              drawing the lines for the profiles. If colouring according to the 
              values of a variable is required then use \code{colour.column.pf}.}
 \item{colour.column.pf}{A \code{\link{character}} giving the name of a column 
              in \code{data} over whose values the colours of the lines are to be 
              varied. The colours can be specified using \code{colour.values.pf}.}
 \item{colour.values.pf}{A \code{\link{character}} vector specifying the values of 
              the colours to use in drawing the lines for the profiles.  
              If this is a named vector, then the values will be matched based 
              on the names. If unnamed, values will be matched in order 
              (usually alphabetical) within the limits of the scale.}
 \item{alpha.pf }{A \code{\link{numeric}} specifying the degrees of transparency to 
              be used in plotting. It is a ratio in which the denominator 
              specifies the number of points (or lines) that must be overplotted 
              to give a solid cover.}
 \item{addMediansWhiskers.pf}{A \code{\link{logical}} indicating whether plots over 
             time of the medians and outer whiskers are to be added to the profile plot. 
             The outer whiskers are related to the whiskers on a box-and-whisker 
             and are defined as the median plus (and minus) 1.5 times the 
             interquartile range (IQR). Points lying outside the whiskers are 
             considered to be potential outliers.}
 \item{ggplotFuncsProfile}{A \code{\link{list}}, each element of which contains the 
              results of evaluating a \code{\link{ggplot}} function. 
              It is created by calling the \code{\link{list}} function with 
              a \code{\link{ggplot}} function call for each element. Note that 
              these functions are applied to the profiles plots only.}
 \item{plots.by.med}{A \code{\link{character}} that give the names of the set of 
             \code{\link{factor}}s by which medians deviations data is to be 
             grouped and a separate plot produced for each group. If \code{NULL}, 
             no groups are formed. If a set of \code{\link{factor}}s, such as  
             \code{Type}, \code{Tuning} and \code{Method}, that uniquely index the 
             combinations of the smoothing-parameter values is specified, then 
             groups are formed for each combination of the levels of the these 
             \code{\link{factor}}s, and a separate plot is produced for each combination.}
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
              order.}
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
              a \code{\link{ggplot}} function call for each element. Note that 
              these functions are applied to the \code{compare.median} deviations plots only.}
 \item{ggplotFuncsDevnBoxes}{A \code{\link{list}}, each element of which contains 
              the results of evaluating a \code{\link{ggplot}} function. 
              It is created by calling the \code{\link{list}} function with 
              a \code{\link{ggplot}} function call for each element.
              These functions are applied in creating the \code{ggplot} 
              object for deviations boxplot only.}
 \item{...}{allows passing of arguments to \code{\link{plotProfiles}}.}
}
\value{A \code{\link{smooths.frame}} that contains the unsmoothed and smoothed data 
       in long format. That is, all the values for either an unsmoothed or a 
       smoothed trait are in a single column. The smooths for a trait for the 
       different combinatons of the smoothing parameters are placed in rows one 
       below the other. The columns that are included in the \code{\link{smooths.frame}} 
       are \code{Type}, \code{TunePar}, \code{TuneVal}, \code{Tuning} and \code{Method}, 
       as well as those specified by \code{individuals}, \code{times}, \code{response}, and 
       \code{response.smoothed}. and any included in the \code{keep.columns}, \code{plots} 
       and \code{facet} arguments. 
       If \code{traits.types} includes \code{AGR} or \code{RGR}, then the included growth 
       rate(s) of the \code{response} and \code{response.smoothed} must be present, unless 
       \code{get.rates} is TRUE when they and the differences between the \code{times} 
       used in calculating the rates will be computed and added. Then, the names of the 
       growth rates are formed from \code{response} and \code{response.smoothed} by 
       appending \code{.AGR} and \code{.RGR} as appropriate; the name of the column with 
       the \code{times} differences will be formed by appending \code{.diffs} to the 
       value of \code{times}. The \code{external.smooths} will also be included. 
       A \code{\link{smooths.frame}} has the attributes described in 
       \code{\link{smooths.frame}}. 

       Columns in the supplied \code{\link{data.frame}} that have not been used in 
       \code{probeSmooths} will not be included in the returned \code{\link{smooths.frame}}. 
       If they might be needed subsequently, such as when extra plots are produced, they 
       can be included in the \code{\link{smooths.frame}} by listing them in a 
       \code{\link{character}} vector for the \code{keep.columns} argument.

       The \code{\link{smooths.frame}} is returned invisibly.}
\author{Chris Brien}
\seealso{\code{\link{traitSmooth}}, \code{\link{smoothSpline}}, \code{\link{byIndv4Times_SplinesGRs}}, 
         \code{\link{byIndv4Times_GRsDiff}}, \code{\link{smooth.spline}}, \code{\link{psNormal}}, \cr 
         \code{\link{plotSmoothsComparison}}, \code{\link{plotSmoothsMedianDevns}},  
         \code{\link{ggplot}}.}
\examples{
data(exampleData)
vline <- list(ggplot2::geom_vline(xintercept=29, linetype="longdash", size=1))
probeSmooths(data = longi.dat, response = "PSA", response.smoothed = "sPSA", 
             times = "DAP", df = c(4,7), 
             facet.x.pf = "Tuning", facet.y.pf = c("Smarthouse", "Treatment.1"), 
             alpha.pf  = 0.4, colour.column.pf = "Method", 
             colour.values.pf = c("orange", "olivedrab"), 
             ggplotFuncsProfile=vline)

#An example that includes comparisons with logistic model fits

if (requireNamespace("nlme", quietly = TRUE))
{
  extra.dat <- longi.dat[, -grep("sPSA", names(longi.dat), fixed = TRUE)]  
  logist.grp <- nlme::groupedData(PSA ~ xDAP | Snapshot.ID.Tag, data = longi.dat)
  #Fit the simple logistic model  
  logist.lis <- nlme::nlsList(SSlogis, logist.grp, na.action = na.pass)
  logist.dat <- within(extra.dat, sPSA <- fitted(logist.lis))
  logist.dat <- byIndv4Times_GRsDiff(data = logist.dat, 
                                     response = "sPSA", 
                                     individuals="Snapshot.ID.Tag",
                                     times="DAP", 
                                     which.rates = c("AGR", "RGR"))
  logist.dat <- cbind(Tuning = factor("Logistic"), logist.dat)
  #Fit the four-parameter logistic model - generates warnings
  logis4.lis <- suppressWarnings(nlme::nlsList(SSfpl, logist.grp, na.action = na.pass))
  logis4.dat <- within(extra.dat, sPSA <- fitted(logis4.lis))
  logis4.dat <- byIndv4Times_GRsDiff(data = logis4.dat, 
                                     response = "sPSA", 
                                     individuals="Snapshot.ID.Tag",
                                     times="DAP", 
                                     which.rates = c("AGR", "RGR"))
  logis4.dat <- cbind(Tuning = factor("Logis-4par"), logis4.dat)
  #Combine the logistic fits
  extra.dat <- rbind(logist.dat,logis4.dat)
  extra.dat <- cbind(Type = factor("NonLinear"), extra.dat)

  #Compare spline and logistic fits
  probeSmooths(data = longi.dat, response = "PSA", response.smoothed = "sPSA", 
               times = "DAP", 
               df = c(4,7), external.smooths = extra.dat, 
               which.plots = "median", 
               plots.by.med = "Type", plots.group.med = "Tuning",  
               propn.types.med = c(0.02,0.1, 0.2),
               ggplotFuncsMedDevn = vline)
}

#An example that supplies three smoothing schemes to be compared
data(tomato.dat)
spar.schemes <- data.frame(Type = c("N", "NCS", "P"), 
                           TunePar = c("df", "df", "lam"),
                           TuneVal = c(4, 6, 1),
                           Method = c("dir", "log", "log"))
probeSmooths(data = tomato.dat, 
             response = "PSA", response.smoothed = "sPSA", 
             times = "DAP", 
             smoothing.schemes = spar.schemes,
             which.plots = "medians.deviations",
             plots.group.med = c("Type", "Tuning", "Method"),              
             propn.types.med = NULL)

}
\keyword{hplot}
\keyword{manip}