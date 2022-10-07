\name{traitSmooth}
\alias{traitSmooth}
\title{Obtain smooths for a trait by fitting spline functions and, having compared 
       several smooths, allows one of them to be chosen and returned in a \code{data.frame}.}
\description{Takes a \code{response} that has been observed for a set of 
             \code{individuals} over a number \code{times} and carries out the following steps:
             \describe{
                \item{Smooth:}{Produces \code{response.smoothed} using splines for a set of 
                      smoothing parameter settings and, optionally, computes growth 
                      rates either as differences or derivatives. 
                      (see \code{smoothing.args} below and \code{\link{args4smoothing}})} 
                      This step is bypassed if a \code{\link{data.frame}} that is also of 
                      class \code{\link{smooths.frame}} is supplied to \code{data}.
                \item{Profile plots:}{Produces profile plots of \code{response.smoothed} and its 
                      growth rates that compare the smooths; also, boxplots of the 
                      deviations of the observed from smoothed data can be obtained. 
                      (see \code{profile.plot.args} below and \code{\link{args4profile.plot}})} 
                      Whether these plots are produced is controlled via \code{which.plots}. 
                \item{Median deviations plots:}{Produces plots of the medians of the deviations of the observed 
                      \code{response}, and its growth rates, from \code{response.smoothed},  
                      and its growth rates. These aid in the assessment of the different smooths. 
                      (see \code{meddevn.plot.args} below and \code{\link{args4meddevn.plot}})}
                      Whether these plots are produced is controlled via \code{which.plots}. 
                \item{Choose a smooth:}{Extract a single, favoured \code{response.smoothed}, and its growth rates, 
                      for a chosen set of smoothing parameter settings. 
                      (see \code{chosen.smooth.args} below and \code{\link{args4chosen.smooth}})} 
                      This step will be omitted if \code{chosen.smooth.args} is \code{NULL}. 
                \item{Chosen smooth plot:}{Produces profile plots of the chosen smooth and its growth rates. 
                      (see \code{chosen.plot.args} below and \code{\link{args4chosen.plot}})} 
                      The \code{\link{logical}} argument \code{chosen.plot} determines whether 
                      or not this plot is produced,
}
             Each of the `\code{args4}' functions has a set of defaults that will be used if 
             the corresponding argument, ending in `\code{.args}', is omitted. The defaults 
             have been optimized for \code{traitSmooth}. 

             Input to the function can be either a \code{\link{data.frame}}, that contains data to 
             be smoothed, or a \code{\link{smooths.frame}}, that contains data that has been smoothed. 
             The function can be run (i) without saving any output, (ii) saving the complete set 
             of smooths in a \code{\link{data.frame}} that is also of class \code{\link{smooths.frame}}, 
             (iii) saving a subset of the smooths in a supplied \code{\link{smooths.frame}}, or 
             (iv) saving a single smooth in a \code{\link{data.frame}}, which can be merged 
             with a pre-existing \code{\link{data.frame}} such as the \code{\link{data.frame}} 
             that contains the unsmoothed data.  
}
\usage{
traitSmooth(data, response, response.smoothed, individuals, times, 
            keep.columns = NULL, 
            get.rates = TRUE, 
            rates.method="differences", ntimes2span = NULL, 
            trait.types = c("response", "AGR", "RGR"), 
            smoothing.args = args4smoothing(), 
            x.title = NULL, y.titles = NULL, 
            which.plots = c("profiles", "medians.deviations"), 
            profile.plot.args = args4profile.plot(), 
            meddevn.plot.args = args4meddevn.plot(), 
            chosen.smooth.args = args4chosen.smooth(), 
            chosen.plot.args = args4chosen.plot(), 
            ggplotFuncsDevnBoxes = NULL,
            mergedata = NULL, 
            ...)
}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the data or a \code{\link{smooths.frame}} 
             as is produced by \code{probeSmooths}. if \code{data} is not a 
             \code{\link{smooths.frame}}, then smoothing will be performed. If \code{data} 
             is a \code{\link{smooths.frame}}, then the plotting and selection of smooths 
             will be performed as specified by \code{smoothing.args}, \code{which.plots} 
             and \code{chosen.smooth.args}.}
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
 \item{keep.columns}{A \code{\link{character}} vector giving the names of columns 
             from \code{data} that are to be included in the 
             \code{\link{smooths.frame}} that will be returned.}
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
 \item{x.title}{Title for the x-axis, used for all plots. If \code{NULL} then set 
             to \code{times}.}
 \item{y.titles}{A \code{\link{character}} giving the titles for the y-axis, one 
             for each the response, the AGE and the RGR. They are used for all plots. 
             If \code{NULL} then they are set to the \code{response} and the  
             \code{response} with \code{.AGR} and \code{.RGR} appended.}
 \item{which.plots}{A \code{\link{logical}} indicating which plots of the smooths 
             specified by \code{smoothing.args} are to be produced. 
             The options are either \code{none} or some combination of \code{profiles}, 
             \code{absolute.boxplots}, \code{relative.boxplots} and 
             \code{medians.deviations}. The various \code{profiles} plots that can be 
             poduced are described in the introduction to this function.} 
 \item{profile.plot.args}{A named \code{\link{list}} that is most easily 
             generated using \code{\link{args4profile.plot}}, it documenting the 
             options available for varying profile plots and boxplots. \emph{Note 
             that if \code{\link{args4profile.plot}} is to be called to change 
             from the default settings given in the default \code{traitSmooth} call 
             and some of those settings are to be retained, then the arguments 
             whose settings are to be retained must also be included in the call 
             to \code{\link{args4profile.plot}}; be aware that if you call 
             \code{\link{args4profile.plot}}, then the defaults for this call are 
             those for \code{\link{args4profile.plot}}, \bold{NOT} the call to 
             \code{\link{args4profile.plot}} shown as the default for 
             \code{traitSmooth}.}}
 \item{meddevn.plot.args}{A named \code{\link{list}} that is most easily 
             generated using \code{\link{args4meddevn.plot}}, it documenting the 
             options available for varying median deviations plots. \emph{Note 
             that if \code{\link{args4meddevn.plot}} is to be called to change 
             from the default settings given in the default \code{traitSmooth} 
             call and some of those settings are to be retained, then the arguments 
             whose settings are to be retained must also be included in the call 
             to \code{\link{args4meddevn.plot}}; be aware that if you call 
             \code{\link{args4meddevn.plot}}, then the defaults for this call are 
             those for \code{\link{args4meddevn.plot}}, \bold{NOT} the call to 
             \code{\link{args4meddevn.plot}} shown as the default for 
             \code{traitSmooth}.}}
\item{chosen.smooth.args}{A named \code{\link{list}} with just one element or 
             \code{NULL} for each component. It is most easily generated using 
             \code{\link{args4chosen.smooth}} with \code{combinations} set to 
             \code{single}. The call to \code{\link{args4smoothing}} should give 
             the settings of \code{smoothing.methods}, \code{spline.types}, 
             \code{df} and \code{lambdas} for a single smooth that is to be 
             extracted and that is amongst the smooths that have been produced 
             for the settings specified in \code{smoothing.methods}. If both 
             \code{df} and \code{lambda} in \code{chosen.smooth.args} are \code{NULL}, 
             then, depending on the settings for \code{spline.type} and 
             \code{smoothng.method}, the value of either \code{df} or \code{lambdas} 
             that is the median value or the observed value immediatly below the 
             median value will be added to \code{chosen.smooth.args}. Otherwise, one of 
             \code{df} and \code{lambda} should be \code{NULL} and the other should 
             be a single \code{numeric} value. If a value in \code{chosen.smooth.args} 
             is not amongst those investigated, a value that was investigated will 
             be substituted.}
 \item{chosen.plot.args}{A named \code{\link{list}} that is most easily 
             generated using \code{\link{args4chosen.plot}}, it documenting the 
             options available for varying profile plots. Because this plot 
             includes only a single smooth, the \code{chosen.smooth.args}, the 
             smoothing-parameter \code{\link{factor}}s are unnecessary and an 
             error will be given if any are included. \emph{Note that if 
             \code{\link{args4profile.plot}} is to be called to change from the 
             default settings given in the default \code{traitSmooth} call, 
             then you will need to set all of the following areguments: 
             \code{plots.by}, \code{facet.x}, \code{facet.y} and 
             \code{include.raw}.} If set to \code{NULL}, then no chosen-smooth plot 
             will be produced.}
 \item{ggplotFuncsDevnBoxes}{A \code{\link{list}}, each element of which contains 
              the results of evaluating a \code{\link{ggplot}} function. 
              It is created by calling the \code{\link{list}} function with 
              a \code{\link{ggplot}} function call for each element.
              These functions are applied in creating the \code{ggplot} 
              object for deviations boxplots only.}
 \item{mergedata}{A \code{\link{data.frame}} that is to have the values for the 
             \code{trait.types} for the smooth specified by \code{chosen.smooth.args} 
             \code{\link{merge}}d with it. It must contain columns with the names given 
             in \code{individuals} and \code{times}, and for which there is only one 
             row for each combination of unique values in these columns. In general, 
             it will be that the number of rows in \code{mergedata} is equal to the 
             number of unique combinations of the values in the columns of the 
             \code{chosen.smooth.args} whose names are given by \code{individuals} and 
             \code{times}, but this is not mandatory. If only one smooth has been 
             produced, then it will be merged with \code{data} provided {mergedata} is 
             \code{NULL} and \code{data} is not  a \code{\link{smooths.frame}}. Othewrwise, 
             a single smooth will be be merged with \code{mergedata}.}
 \item{...}{allows arguments to be passed to \code{\link{plotProfiles}}.}
}
\value{A \code{\link{smooths.frame}} or a \code{\link{data.frame}} that contains 
       the unsmoothed and smoothed data in long format. That is, all the values for 
       either an unsmoothed or a smoothed trait are in a single column. 
       
       A \code{\link{smooths.frame}} will be returned when \code{chosen.smooth.args} is 
       \code{NULL} and there is more than one smooth specified by the smoothing 
       parameter arguments. It will contain the smooths for a trait for the different 
       combinatons of the smoothing parameters, the values for the different smooths 
       being placed in rows one below the other. The columns that are included in the 
       \code{\link{smooths.frame}} are \code{Type}, \code{TunePar}, \code{TuneVal}, 
       \code{Tuning} and \code{Method}, as well as those specified by 
       \code{individuals}, \code{times}, \code{response}, and \code{response.smoothed}, 
       and any included in the \code{keep.columns}, \code{plots} and \code{facet} 
       arguments when the smooths were produced. The \code{AGR} or \code{RGR} for the 
       \code{response} and \code{response.smoothed}, if obtained, will also be 
       included. A \code{\link{smooths.frame}} has the attributes described in 
       \code{\link{smooths.frame}}. 

       A \code{\link{data.frame}} will be returned when \code{chosen.smooth.args} is 
       not \code{NULL} or it and \code{mergedata} are \code{NULL}, \code{data} is not a 
       \code{smooths.frame} and there is only one smooth specified by 
       the smoothing parameter arguments. It will include the columns \code{individuals}, 
       \code{times}, \code{response}, and \code{response.smoothed}, and any included in 
       the \code{keep.columns}, \code{plots} and \code{facet} arguments, as well as any 
       growth rates calculated as a result of \code{get.rates} and \code{trait.type}. If 
       \code{mergedata} is not \code{NULL}, the \code{\link{data.frame}} containing the 
       \code{chosen.smooth} or a single smooth, when that is all that is produced, 
       will be merged with the \code{\link{data.frame}} specified by \code{mergedata}.

       The \code{\link{smooths.frame}}/\code{\link{data.frame}} is returned invisibly.}
\details{This function is a wrapper function for \code{\link{probeSmooths}}, 
         \code{\link{plotSmoothsComparison}}, \code{\link{plotSmoothsComparison}} 
         and \code{\link{plotDeviationsBoxes}}. It uses the helper functions 
         \code{\link{args4smoothing}}, \code{\link{args4profile.plot}} and 
         \code{\link{args4meddevn.plot}} to se arguments that control the smoothing 
         and plotting.

             It takes a \code{response} that has been observed for a set of 
             \code{individuals} over a number \code{times} and produces 
             \code{response.smoothed}, using \code{\link{probeSmooths}}, for a default 
             set of smoothing parameter settings (see \code{\link{args4smoothing}} 
             for the defaults). The settings can be varied from the defaults by 
             specifying alternate values for the smoothing parameters, the parameters 
             being the type of spline (\code{spline.types}), the degrees of freedom 
             (\code{df}) or smoothing penalty (\code{lambdas}) and \code{smoothing.methods}. 
             There are also several other smoothing arguments that can be manipulated to 
             affect the smooth (for details see \code{\link{args4smoothing}}). 
             The secondary traits of the absolute growth rate (\code{AGR}) and relative 
             growth rate (\code{RGR}) are calculated from the two primary traits, the 
             \code{response} and \code{response.smoothed}.

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

             By default, the single smooth for an arbitrarily chosen combination of the 
             smoothing parameters is returned by the function. The smooth for a single 
             combination other than default combination can be nominated for return 
             using the \code{chosen.smooth.args} argument. This combination must involve 
             only the supplied values of the smoothing parameters. The values for 
             \code{response}, the \code{response.smoothed} and their \code{AGR}s 
             and \code{RGR}s are are added to \code{data}, after any pre-existing  
             columns of these have been removed from \code{data}. Profile plots of 
             the three smoothed traits are produced using \code{\link{plotProfiles}}. 
             However, if \code{chosen.smooth.args} is \code{NULL}, all of the smooths will be 
             returned in a \code{\link{smooths.frame}}, and plots for the single 
             combination of the smoothing parameters will not be produced. 
}
\author{Chris Brien}
\seealso{\code{\link{args4smoothing}}, \code{\link{args4meddevn.plot}}, \code{\link{args4profile.plot}}, 
        \code{\link{args4chosen.smooth}}, \cr 
        \code{\link{args4chosen.plot}}, \code{\link{probeSmooths}} 
        \code{\link{plotSmoothsComparison}} and \cr 
        \code{\link{plotSmoothsMedianDevns}}, \code{\link{ggplot}}.}
\examples{
data(exampleData)
vline <- list(ggplot2::geom_vline(xintercept=29, linetype="longdash", size=1))
yfacets <- c("Smarthouse", "Treatment.1")
smth.dat <- traitSmooth(data = longi.dat, 
                        response = "PSA", response.smoothed = "sPSA",
                        individuals = "Snapshot.ID.Tag",times = "DAP", 
                        keep.columns = yfacets, 
                        profile.plot.args = 
                          args4profile.plot(facet.y = yfacets, 
                                            ggplotFuncs = vline),
                        chosen.plot.args = 
                          args4chosen.plot(facet.y = yfacets,  
                                           ggplotFuncs = vline))
}
\keyword{hplot}
\keyword{manip}