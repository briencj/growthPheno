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
             \code{smoothing.methods}. The values of these, and other, smoothing 
             arguments are set using the helper function 
             \code{\link{args4smoothing}}.

             Provided \code{get.rates} is \code{TRUE} or includes \code{raw} and/or 
             \code{smoothed} and depending on the setting 
             of \code{trait.types}, the Absolute Growth Rates (AGR) and/or the 
             Relative Growth Rates (RGR) are calculated for each individual from the 
             unsmoothed, observed \code{response} and from the smooths of the 
             \code{response}, using either differences or first derivatives, as 
             specified by \code{rates.method}. 

             Generally, profile plots for the traits (a response, an AGR or an RGR) 
             specified in \code{traits.types} are produced if \code{which.plots} is 
             \code{profiles}; if \code{which.plots} specifies one or more deviations 
             plots, then those deviations plots will also be produced, these being 
             based on the unsmoothed data from which the smoothed data has been 
             subtracted. The layout of the plots is controlled via combinations of 
             one or more of the smoothing-parameter \code{\link{factor}}s \code{Type}, 
             \code{TunePar}, \code{TuneVal}, \code{Tuning} (the combination 
             of \code{TunePar} and \code{TuneVal}) and \code{Method}, as well as 
             other \code{\link{factor}}s associated with the \code{data}. 
             The \code{\link{factor}}s that are to be used for the profile plots and 
             deviations boxplots are supplied via the argument \code{profile.plot.args} 
             using the helper function \code{\link{args4profile.plot}} to set 
             \code{plots.by}, \code{facet.x}, and \code{facet.y}; for the plots of 
             the medians of the deviations, the \code{\link{factor}}s are supplied via 
             the argument \code{meddevn.plot.args} using the helper function 
             \code{\link{args4meddevn.plot}} to set \code{plots.by}, \code{facet.x}, 
             \code{facet.y} and \code{plots.group}. Here, the basic principle is that 
             the number of levels combinations of the smoothing-parameter 
             \code{\link{factor}}s included in the set of \code{plots} and 
             \code{facets} arguments to one of these helper functions must be the same 
             as those covered by the combinations of the values supplied to 
             \code{spline.types}, \code{df}, \code{lambdas} and \code{smoothing.methods} 
             and incorporated into the \code{\link{smooths.frame}}, such as is 
              returned by \code{probeSmooths}. This ensures that smooths 
             from different parameter sets are not pooled together in a single plot. 
             It is also possible to include \code{\link{factor}}s that are not 
             smoothing-parameter \code{\link{factor}}s in the \code{plots} amd 
             \code{facets} arguments. 

             The following profiles plots can be produced using 
             \code{\link{args4profile.plot}}: (i) separate plots of the 
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

             Deviation plots that can be produced are the absolute and relative deviations 
             boxplots and plots of medians deviations (see \code{which.plots}). 

             The handling of missing values is controlled via \code{na.x.action} and 
             \code{na.y.action} supplied to the helper function \code{\link{args4smoothing}}.

             The \code{probeSmooths} arguments are grouped according to function in the following order:
             \enumerate{
                \item{\bold{Data description arguments:} data, response, response.smoothed, 
                      individuals, times, keep.columns, trait.types, get.rates, 
                      rates.method, ntimes2span.}
                \item{\bold{Smoothing arguments:} smoothing.args (see \code{\link{args4smoothing}}).}
                \item{\bold{General plot control:} x.title, y.titles, facet.labeller, which.plots.} 
                \item{\bold{Profile plots (pf) features:} profile.plot.args (see \code{\link{args4profile.plot}})} 
                \item{\bold{Median-deviations (med) plots features:} meddevn.plot.args (see \code{\link{args4meddevn.plot}})} 
             }
}
\usage{
probeSmooths(data, response = "PSA", response.smoothed = NULL, 
             individuals="Snapshot.ID.Tag", times = "DAP", 
             keep.columns = NULL, 
             get.rates = TRUE, 
             rates.method="differences", ntimes2span = NULL, 
             trait.types = c("response", "AGR", "RGR"), 
             smoothing.args = 
               args4smoothing(smoothing.methods = "direct", 
                              spline.types = "NCSS", 
                              df = NULL, lambdas = NULL), 
             x.title = NULL, y.titles = NULL, 
             which.plots = "profiles", 
             profile.plot.args = 
               args4profile.plot(plots.by = NULL, 
                                 facet.x = ".", facet.y = ".", 
                                 include.raw = "no"), 
             meddevn.plot.args = 
               args4meddevn.plot(plots.by = NULL, plots.group = NULL, 
                                 facet.x = ".", facet.y = ".",
                                 propn.note = TRUE, 
                                 propn.types = c(0.1, 0.5, 0.75)), 
             ggplotFuncsDevnBoxes = NULL, ...)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the data or a \code{\link{smooths.frame}} 
             as is produced by \code{probeSmooths}. if \code{data} is not a 
             \code{\link{smooths.frame}}, then smoothing will be performed. If \code{data} 
             is a \code{\link{smooths.frame}}, then the plotting and selection of smooths 
             will be performed as specified by \code{smoothing.args} and \code{which.plots}.}
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
 \item{get.rates}{A \code{\link{logical}} or a \code{\link{character}} specifying 
             which of the \code{response} and the \code{response.smoothed} are to 
             have growth rates (AGR and/or RGR) computed and stored. If set to 
             \code{TRUE} or \code{c("raw", "smoothed")}, growth rates will be 
             obtained for both. Setting to only one of \code{raw} or 
             \code{smoothed}, results in the growth rates for either the 
             \code{response} or the \code{response.smoothed} being computed. 
             If set to \code{none} or \code{FALSE}, no growth rates ar computed. 
             Which growth.rates are computed can be changed using the arguments 
             \code{traits.types} and the method used for computing them for the 
             \code{response.smooth} by \code{rates.method}. The growth rates for 
             the \code{response} can only be computed by \code{differencing}.}
 \item{rates.method}{A \code{\link{character}} specifying the method to use in 
             calculating the growth rates for \code{response.smoothed}. The two 
             possibilities are \code{"differences"} and \code{"derivatives"}.}
 \item{ntimes2span}{A \code{\link{numeric}} giving the number of values in 
             \code{times} to span in calculating growth rates by differencing. 
             For \code{ntimes2span} set to \code{NULL}, if \code{rates.method} 
             is set to \code{differences} then \code{ntimes2span} is set to 2; 
             if \code{rates.method} is set to \code{derivatives} then 
             \code{ntimes2span} is set to 3. 
             Note that when \code{get.rates} is includes \code{raw} or is 
             \code{TRUE}, the growth rates for the unsmoothed response must be 
             calculated by differencing, even if the growth rates for the 
             smoothed response are computed using \code{derivatives}. When 
             differencing, each growth rate is calculated as the difference in 
             the values of one of the \code{responses} for pairs of \code{times} 
             values that are spanned by \code{ntimes2span} \code{times} values 
             divided by the difference between this pair of \code{times} values. 
             For \code{ntimes2span} set to 2, a growth rate is the  
             difference between consecutive pairs of values of one of the 
             \code{responses} divided by the difference between consecutive 
             pairs of \code{times} values.}
 \item{trait.types}{A \code{\link{character}} giving the \code{trait.types} that 
             are to be plotted. If growth rates are included in \code{trait.types}, 
             then they will be computed for either the \code{response} and/or the 
             \code{response.smoothed}, depending on the setting of \code{get.rates}. 
             Any growth rates included in \code{trait.types} for the \code{response} 
             that are available in \code{data}, but have not been specified for 
             computation in \code{get.rates}, will be retained in the returned 
             \code{smooths.frame}.
             If \code{all}, the \code{response.smoothed}, its \code{AGR} and 
             \code{RGR}, will be plotted. The \code{response}, and its \code{AGR} 
             and \code{RGR}, will be plotted as the plotting options require it.}
 \item{smoothing.args}{A \code{\link{list}} that is most easily generated 
             using \code{\link{args4smoothing}}, it documenting the 
             options available for smoothing the data. It gives the settings of 
             \code{smoothing.methods}, \code{spline.types}, \code{df}, \code{lambdas}, 
             \code{smoothing.segments}, \code{npspline.segments}, \code{na.x.action}, 
             \code{na.y.action}, \code{external.smooths}, and \code{correctBoundaries}, 
             to be used in smoothing the \code{response} or in selecting a subset 
             of the smooths in \code{data}, depending on whether \code{data} is a 
             \code{\link{data.frame}} or a \code{\link{smooths.frame}}, respectively. 
             Set \code{smoothing.args} to \code{NULL} if \code{data} is a 
             \code{\link{smooths.frame}} and only plotting or extraction of a chosen 
             smooth is required.}
 \item{x.title}{Title for the x-axis, used for all plots. If \code{NULL} then set to 
             \code{times}.}
 \item{y.titles}{A \code{\link{character}} giving the titles for the y-axis, 
             one for each trait specified by \code{trait.types} and used for all plots. 
             If \code{NULL} then set to the traits derived for \code{response} 
             from \code{trait.types}.}
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
 \item{profile.plot.args}{A named \code{\link{list}} that is most easily 
             generated using \code{\link{args4profile.plot}}, it documenting the 
             options available for varying profile plots and boxplots. \emph{Note 
             that if \code{\link{args4profile.plot}} is to be called to change 
             from the default settings given in the default \code{probeSmooths} call 
             and some of those settings are to be retained, then the arguments 
             whose settings are to be retained must also be included in the call 
             to \code{\link{args4profile.plot}}; be aware that if you call 
             \code{\link{args4profile.plot}}, then the defaults for this call are 
             those for \code{\link{args4profile.plot}}, \bold{NOT} the call to 
             \code{\link{args4profile.plot}} shown as the default for 
             \code{probeSmooths}.}}
 \item{meddevn.plot.args}{A named \code{\link{list}} that is most easily 
             generated using \code{\link{args4meddevn.plot}}, it documenting the 
             options available for varying median deviations plots. \emph{Note 
             that if \code{\link{args4meddevn.plot}} is to be called to change 
             from the default settings given in the default \code{probeSmooths} call 
             and some of those settings are to be retained, then the arguments 
             whose settings are to be retained must also be included in the call 
             to \code{\link{args4meddevn.plot}}; be aware that if you call 
             \code{\link{args4meddevn.plot}}, then the defaults for this call are 
             those for \code{\link{args4meddevn.plot}}, \bold{NOT} the call to 
             \code{\link{args4meddevn.plot}} shown as the default for 
             \code{probeSmooths}.}}
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
       If \code{trait.types} includes \code{AGR} or \code{RGR}, then the included growth 
       rate(s) of the \code{response} and \code{response.smoothed} must be present, unless 
       \code{get.rates} is TRUE or includes \code{raw} and/or \code{smoothed}. In this case, 
       the growth rates specified by \code{trait.types} will be calculated for the responses 
       nominated by \code{get.rates} and the differences between the \code{times} 
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
\seealso{\code{\link{args4smoothing}}, , \code{\link{args4meddevn.plot}}, \code{\link{args4profile.plot}}, 
         \code{\link{traitSmooth}}, \cr 
         \code{\link{smoothSpline}}, \code{\link{byIndv4Times_SplinesGRs}}, \code{\link{byIndv4Times_GRsDiff}},  
         \code{\link{smooth.spline}}, \code{\link{psNormal}}, \code{\link{plotSmoothsComparison}}, \cr
         \code{\link{plotSmoothsMedianDevns}}, \code{\link{ggplot}}.}
\examples{
data(exampleData)
vline <- list(ggplot2::geom_vline(xintercept=29, linetype="longdash", size=1))
probeSmooths(data = longi.dat, response = "PSA", response.smoothed = "sPSA", 
             times = "DAP", 
             smoothing.args = args4smoothing(df = c(4,7)), 
             profile.plot.args = 
               args4profile.plot(plots.by = NULL, 
                                 facet.x = "Tuning", 
                                 facet.y = c("Smarthouse", "Treatment.1"), 
                                 include.raw = "no", 
                                 alpha  = 0.4, 
                                 colour.column = "Method", 
                                 colour.values = c("orange", "olivedrab"), 
                                 ggplotFuncs = vline))

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
  probeSmooths(data = longi.dat, 
               response = "PSA", response.smoothed = "sPSA", 
               times = "DAP", 
               smoothing.args = 
                 args4smoothing(smoothing.methods = "direct", 
                                spline.types = "NCSS", 
                                df = c(4,7), lambdas = NULL, 
                                external.smooths = extra.dat), 
               which.plots = "median", 
               meddevn.plot.args = 
                 args4meddevn.plot(plots.by = "Type", 
                                   plots.group = "Tuning",  
                                   facet.x = ".", facet.y = ".", 
                                   propn.types = c(0.02,0.1, 0.2),
                                   ggplotFuncs = vline))
}

#An example that supplies three smoothing schemes to be compared
data(tomato.dat)
probeSmooths(data = tomato.dat, 
             response = "PSA", response.smoothed = "sPSA", 
             times = "DAP", 
             smoothing.args = 
               args4smoothing(spline.types      = c(  "N", "NCS",   "P"), 
                              df                = c(   4,     6,    NA), 
                              lambdas           = c(  NA,    NA,     1), 
                              smoothing.methods = c("dir", "log", "log"),
                              combinations      = "parallel"),
             which.plots = "medians.deviations",
             meddevn.plot.args = 
               args4meddevn.plot(plots.by = NULL, 
                                 plots.group = c("Type", "Tuning", "Method"),              
                                 facet.x = ".", facet.y = ".", 
                                 propn.note = FALSE, propn.types = NULL))

}
\keyword{hplot}
\keyword{manip}