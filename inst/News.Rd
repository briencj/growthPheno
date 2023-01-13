\name{news}
\title{News for Package growthPheno}
\section{growthPheno Version 2.1-17 (2023-01-13)}{
  \itemize{
  \item{Fix bug in plotting of Method in a propn.type envelope.}
  \item{Fix bug in plotComparison when \code{include.raw} is set to include the raw data plot in a facet of the plot and the facet is set to \code{"."}.}
  \item{Replace \code{size} with \code{linewidth} in ggplot line function calls to comply with changes made to ggplot2.}
  }}
\section{growthPheno Version 2.1-16 (2022-11-13)}{
  \itemize{
  \item{Build for R 4.2.2 and for CRAN submission.}
  \item{Add check that the factors names included in \code{chosen.plots.args} are available in \code{data}.}
  }}
\section{growthPheno Version 2.1-15 (2022-11-6)}{
  \itemize{
  \item{Revise documentation of traitSmooth to clarify how to specify \code{smoothing.args} when \code{data} is set to a \code{smooths.frame}.}
  }}
\section{growthPheno Version 2.1-14 (2022-10-31)}{
  \itemize{
  \item{Fix bug in the return of the chosen smooth. If there is more than one smooth specified by \code{smoothing.args}, then the chosen.smooth will only be returned if \code{mergedata} is not \code{NULL}.}
  }}
\section{growthPheno Version 2.1-13 (2022-10-16)}{
  \itemize{
  \item{Reorganized the plotting of deviation boxplots in probeSmooths and TraitSmooth, in the process creating new functions args4devnboxed_plot and plotSmoothsDevnBoxes. 
  \enumerate{
     \item{In plotDeviationsBoxes, renamed \code{labeller} to \code{facet.labeller}, added \code{facet.scales} and \code{angle.x}, and removed \code{df}.}
     \item{Transferred deviations plotting from plotSmoothsComparison to the new function plotSmoothsDevnBoxplots. Removed \code{which.plot} and \code{ggplotFuncsDevnBoxes} arguments from plotSmoothComparisons}
     \item{Added \code{devnboxes.plot.args} to probeSmooths and traitSmooth and removed \code{ggplotFuncsDevnBoxes} argument.}
     \item{Added function args4devnboxes_plot to set \code{devnboxes.plot.args}.}
  }}
  \item{Rebuild the Tomato vignette.}
  }}
\section{growthPheno Version 2.1-12 (2022-10-15)}{
  \itemize{
  \item{Fix bug in args4chosen_plot where \code{angle.x} was being ignored.}
  \item{Fix bug in traitSmooth where \code{trait.types} was being ignored by chosen.cmooth.args.}
  \item{Fix bug in traitSmooth where \code{combinations} was not being passed correctly.}
  \item{Fix bug in probeSmmoths and TraitSmooth in which a smooth was not produced when the \code{combinations} was set to  \code{single} for the \code{smoothing.args} argument.}
  \item{Rebuild the Tomato vignette.}
  }}
\section{growthPheno Version 2.1-11 (2022-10-08)}{
  \itemize{
  \item{Revise the names of args4chosen.plot, args4chosen.smooth and args4meddevn.plot, args4profile.plot to args4chosen_plot, args4chosen_smooth, args4meddevn_plot and args4profile_plot.}
  }}
\section{growthPheno Version 2.1-10 (2022-10-07)}{
  \itemize{
  \item{Revamp the traitSmooth, ProbeSmooths, plotSmoothsComparisons, PlotSmoothsMedians arguments so that groups of them are specified by the new functions args4chosen.plot, args4chosen.smooth, args4meddevn.plot, args4profile.plot and args4smoothing.}
  \item{Revise the Tomato vignette for this version.}
  }}
\section{growthPheno Version 2.0-20 (2022-09-26)}{
  \itemize{
  \item{Allow a \code{smooths.frame} as input to \code{probeSmooths} and \code{traitSmooth} via the \code{data} argument. Modify traitSmooth so that if a \code{smooths.frame} is supplied as input, then only plotting and the selection and plotting of a \code{chosen.smooth} is executed.}
  \item{Add argument \code{mergedata} to traitSmooth and add allow the chosen smooth to be plotted and returned, possibly merged with the \code{data.frame} nominated by \code{mergedata}.}
  \item{Modify \code{get.rates} argument to allow the specification of the options \code{raw} and \code{smoothed} so that the type of responses for which the rates that are to be computed can be specified.}
  \item{Add \code{scales} argument to plotProfiles.}
  \item{Add arguments \code{scales.pf}, \code{scales.med} and \code{scales.chosen} to traitSmooth,  probeSmooths, plotSmoothsComparison and plotMedianDEviatons.}
  }}
\section{growthPheno Version 2.0-19 (2022-09-14)}{
  \itemize{
  \item{Changed the return value for traitSmooth in the case that there is only one smooth specified by the smoothing parameters and the argument \code{chosen.smooth.args} \code{NULL}. Now the smooth will be merged with the data and returned as a \code{data.frame}.}
  \item{Modify \code{breaks.spacing.x} in plotSmoothsComparison, plotMedianDeviations, probeSmooths and traitSmooth to allow for negative values that signify breaks or minor breaks are to be omitted for values that do not occur in the data. Change defaults for \code{breaks.spacing.x} to be negative.}
  \item{Fix bug in the error message in plotImagetimes for columns not in \code{data}.}
  }}
\section{growthPheno Version 2.0-18 (2022-09-06)}{
  \itemize{
  \item{Move argument \code{responses4singletimes} and add argument \code{times.single} in traitExtractFeatures.}
  \item{Rename arguments \code{times.overall} to \code{intvl.overall} in traitExtractFestures.}
  }}
\section{growthPheno Version 2.0-17 (2022-09-04)}{
  \itemize{
  \item{Fix a bug in traitExtractFeatures associated with responses.water and allow parallel looping over multiple values for both water.use and responses.water.}
  \item{Add new arguments \code{responses4overall.rates}, \code{water.use4overall.water} and \code{responses4overall.water} to traitExtractFeatures and make changes to byIndv4Intvl_GRsDiff to allow suffix.interval to be NULL.}
  \item{Change arguments \code{responses.singletimes}. \code{response.rates}, \code{rates.method}, \code{water.use}, \code{responses.water}, \code{responses.total}, \code{responses.max}, \code{times.whole} and \code{suffix.whole} to \code{responses4singletimes}, \code{responses4intvl.rates}, \code{growth.rates.method}, \code{water.use4intvl.traits}, \code{responses4water}, \code{responses4overall.totals}, \code{responses4overall.max}, \code{times.overall} and \code{suffix.overall} in traitExtractFeatures.}
  \item{Add new argument \code{sep.rates} to byIndv4Intvl_GRsDiff, byIndv4Intvl_GRsAvg, byIndv4Times_GRsDiff and byIndv4Times_SplinesGRs, and \code{sep.growth.rates} to traitExtractFeatures.}
  \item{Add new argument \code{sep.water.traits} to byIndv4Intvl_WaterUse and traitExtractFeatures.}
  \item{Add new argument \code{sep.suffix.interval} to byIndv4Intvl_GRsDiff, byIndv4Intvl_GRsAvg, byIndv4Intvl_ValueCalc and byIndv4Intvl_WaterUse, as well as \code{sep.suffix.interval} togetTimesSubset and  traitExtractFeatures.}
  \item{Add new argument \code{sep.FUNname} to byIndv_ValueCalc and byIndv4Intvl_ValueCalc.}
  }}
\section{growthPheno Version 2.0-16 (2022-08-01)}{
  \itemize{
  \item{Fix a bug in plotSmoothsComparison that meant that attempting to change the x-axis breaks using the function \code{ggplot2::scale_x_continuous} were ignored.}
  \item{Fix a bug in plotSmoothsMedianDevn in plotting the envelope note when there are more than two variables in \code{facet.x.med}.}
  \item{Add \code{breaks.spacing.x}, argument to plotSmoothsComparison, plotMedianDeviations, probeSmooths and traitSmooth.}
  \item{Add \code{collapse.facets.x.pf} and \code{collapse.facets.y.pf} arguments and replace \code{plots.include.raw} with the \code{include.raw.pf} argument in plotSmoothsComparison and probeSmooths so that facetting can be based on an arbitrary combinations of smoothing-parameters and other factors.}
  \item{Change the argument \code{suffices.rates} in traitExtractFeatures to \code{suffices.growth.rates} abd add the arguments \code{water.trait.types}, \code{suffix.water.rate}, \code{suffix.water.index}.}
  \item{Add check to probeSmooths that there is no more than one observation for each combination of individuals and times.}
  \item{Add check to traitExtractFeatures that all of the columns speicified in a call to it are in \code{data}.}
  \item{Minor revisions to Tomato vignette.}
  \item{Change names of columns in RicePrepped.dat with \code{DAT} to \code{DAST}.}
  }}
\section{growthPheno Version 2.0-14 (2022-07-03)}{
  \itemize{
  \item{Add traitExtractFeatures as a wrapper function for producing traits that are single-valued features for each individual.}
  \item{Implement a revised naming convention for data column names. The names of the columns with \code{Area} and \code{Water.Loss} values are to be replaced with \code{PSA} and \code{WU}. The corresponding columns of smoothed data are to have an \code{s} prepended to the unsmoothed column names. A factor name is singular, the name of a numeric for a factor has an \code{x} prefixed to the factor name and the name of a centred numeric for a factor has a \code{c} prefixed to the factor name.} 
  \item{Modify the Rice and tomato vignettes and the associated data to follow the new naming convention.}
  \item{Replace plotLongitudinal with plotProfiles, the latter having the \code{xname} argument replaced by \code{times} that bevomes the primary argument for specifying the values for plotting on the x-axis. The argument \code{times} can be supplied as a factor, character or numeric and \code{x} is only required if an expression involving \code{times} is to be used for plotting on the x-axis. The naming of columns follows the new naming conventions.}
  \item{Add \code{times} argument to plotAnom. and change default for \code{x} to \code{NULL}.}
  \item{Replace longitudinalPrime with prepImageData, the latter implementing the new naming conventions. The column \code{xPosn} is no longer included by prepImagesData, but is added by designfactors; the order of the Id and times columns has been changed. The function designfactors has also been modified to produce columns that follow the new naming conventions.}
  \item{Add the function smoothSpline with the same functionality as fitSpline, but simplified specification of growth rates and derivatives.}
  \item{The defaults for the functions getTimesSubset, byIndv4Intvl_GRsAvg, byIndv4Intvl_GRsDiff, byIndv4Intvl_ValueCalc and byIndv4Intvl_WaterUse, byIndv4Times_SplinesGRs, byIndv4Times_GRsDiff, probeSmooths, plotSmoothsComparison, plotsSmoothsMedianDevns, plotProfiles have also been modified to reflect these naming conventions.}
  \item{Fix bugs when \code{extra.derivs} is used in byIndv4Times_SplinesGRs.}
  \item{Fix a big in byIndv_ValueCalc and byIndv4Intvl_ValueCalc that occurs when all the values for a response are \code{NA}.}
  \item{Ensure that unwanted growth rates are not retained by byIndv4Times_SplinesGRs when \code{rates.method}  is set to \code{derivatives}.}
  \item{Have growth rates for smoothed traits computed across the segments for segmented smooths when ntimes2span is 2 in bIndv4Times_SplinesGRs and probeSmooths.}
\item{Make \code{trimx} the default for \code{na.y.action} in fitSpline, byIndv4Times_SplinesGRs and probeSmooths.} 
  \item{Modify the Rice and tomato vignettes and the associated data to follow the naming convention that a factor name is singular, the name of a numeric for a factor has an \code{x} prefixed and the name of a centred numeric for a factor has a \code{c} prefixed.}
  }}
\section{growthPheno Version 2.0-13 (2022-06-19)}{
  \itemize{
  \item{The set of functions traitSmooth, probeSmooths, plotSmoothsComparison, plotSmoothMedianDeviations, splitSplines, splitContGR, splitValueCalculate, intervalGRaverage, intervalGRdiff, intervalValueCalculate, and intervalWUI have been overhauled to make them a consistent set of functions for smoothing and calculations performed on each individual. Where possible the arguments have been simplified. The changes that have been made include:
  \enumerate{
  \item{To make clear what the different functions do, all of The \code{split} and \code{interval} functions have been renamed to begin with \code{byIndv} to make clear that they all repeats operates on each individual from a set of individuals. In the renaming a distinction is made between those that produce, for each individual, (i) values for all of the times (4Times) and (ii) a value for an interval (4Intvl).}
     \item{The three \code{split} functions splitSplines, splitContGR and  splitValueCalculate have been renamed to byIndv4Times_SplinesGRs and byIndv4Times_GRsDiff and byIndv_ValueCalc,}
     \item{The new function byIndv4Times_SplinesGRs can smooth and calculate growth rates, either using derivative or differencing, in a single call using simplified arguments \code{rates.method} and \code{which.rates}.}
     \item{The \code{interval} functions have been renamed to byIndv4Intvl_GRsAvg, byIndv4Intvl_GRsDiff, byIndv4Intvl_ValueCalc and byIndv4Intvl_WaterUse.}
     \item{The new function byIndv4Intvl_WaterUse can produce the three water use traits \code{WU}, \code{WUR}, \code{WUI}, as well as the \code{AGR}, for a time interval. To achieve this \code{include.water.use} has been replaced by \code{trait.types}.}
     \item{The functions traitSmooth and probeSmooths have been revised to use the new \code{byIndv} functions.}
     \item{All functions now 
       \itemize{
         \item{use the \code{times} argument rather than \code{x} or \code{times.factor} and \code{times} accepts a \code{factor}, \code{character} or \code{numeric};}
        \item{use the \code{individuals} argument rather than \code{INDICES};}
        \item{have the initial arguments in the order \code{data}, \code{response(s)}, \code{individuals}, \code{times};}
        \item{have had their arguments reordered so that error may occur if arguments have been matched by position.}
        }
  }}
  \item{The \code{split} and \code{interval} functions have been soft-deprecated. They are no longer being maintained and  will eventually be removed.}
  \item{Changed the argument \code{times.factor} to \code{times} in getTimesSubset.}
  \item{Changed the arguments \code{start.times} and \code{end.times} to \code{start.time} and \code{end.time} in byIndv4Intvl_GRsDiff and byIndv4Intvl_WaterUse.}
  \item{Changed \code{Area} and \code{Water.Loss} to \code{PSA} and \code{WU} in the data set \code{tomato.dat}. Modified the vignette \code{Tomato} to reflect this and adapted it to version 2.x.x of \code{growthPheno}.}
}}}
\section{growthPheno Version 2.0-12 (2022-06-08)}{
  \itemize{
  \item{Generalize splitContGRdiff to calculate growth rates for any span of times using a new argument \code{ntimes2span} and add \code{ntimes2span} to probeSmooths.}
  \item{Add the argument \code{individuals} to splitContGRdiff and splitSplines as a pseudonym for \code{INDICES}, making them consistent with splitValueCalculate and probeSmooths.}
  \item{Changed the \code{pspline.segments} to \code{npspline.segments} in fitSpline, splitSplines and probeSmooths.}
}}
\section{growthPheno Version 2.0-11 (2022-05-25)}{
  \itemize{
  \item{Correct examples for probeSmooths and plotMediansDeviations.}
}}
\section{growthPheno Version 2.0-10 (2022-05-25)}{
  \itemize{
  \item{This version incorporates a major overhaul of the smoothing functionality of growthPheno with new and revised functions, new arguments in some functions and changes in the order of some arguments. The descriptiion of the function has been revised.}
  \item{Add \code{spline.type}, \code{lambda} and \code{npspline.segments} arguments to fitSpline, splitSplines, and  probeSmoothing to incorporate P-spline smoothing (\code{PS}) using the \code{JOPS} package.}
  \item{Replaced \code{AGR} and \code{RGR} with the argument \code{extra.rate} to simplify the specification.}
  \item{Allow the specification of \code{lambda}, as an alternative to \code{df}, for natural cubic smoothing splines (\code{NCSS}).}
  \item{Add the new functions probeSmooths and plotSmoothsMedianDevns to replace probeSmoothing and plotMedianDeviations and the additional function plotSmoothsComparison. The new functions are more flexible than those that they replace, their arguments and options having been substantially changed from those for the old functions. They allow for the comparison of the two \code{spline.types} of natural cubic smoothing splines (\code{NCSS}) and P-splines (\code{PS}) and organize the smoothed data in long format to produce a \code{data.frame} of class \code{smooths.frame}.}
  \item{Add the function traitSmooth as a wrapper function for probeSmooths.}
  \item{Add functions as.smooths.frame, is.smooths.frame and validSmoothsFrame to allow a \code{data.frame} to be assigned the class \code{smooths.frame} and to test whether a \code{data.frame} is of class \code{smooths.frame} and whether it is a valid object of class \code{smooths.frame}.}
  \item{Deprecate function probeDF.}
  \item{Implement segmented smoothing in splitSplines and probeSmooths.}
  \item{Fix bugs in plotLongitudinal and plotMedianDeviations when multiple factors are supplied to either \code{facet.x} or \code{facet.y}. Changed defaults for \code{facet.x} and \code{facet.y} in plotLongitudinal to \code{"."}.}
   \item{Add \code{PSAcolumn} argument to longitudinalPrime.}
   \item{Add \code{response.smoothed} argument to fitSpline, splitSplines, and  probeSmoothing to allow the user to specify arbitrary names for the smoothed response.}
   \item{Mofify fitSpline to return the \code{predictions} and the object returned by either \code{smooth.spline} or \code{JOPS::psNormal}.}
   \item{Fix an obscure bug in fitSpline that occurred when the \code{times.factor} differnces pre-exist in \code{data} and the growth rates are to be calculated for a subset of the times. In such cases, the first growth rate was not \code{NA}.}
   \item{Correct the documentation of WUI to reflect what is actually done in the code, namely that the response is divided by the water.}
}}
\section{growthPheno Version 1.0-36 (2022-02-17)}{
  \itemize{
   \item{Add \code{na.rm} and \code{...} arguments to cumulate.}
}}
\section{growthPheno Version 1.0-35 (2021-11-06)}{
  \itemize{
   \item{Miscellaneous code tidying.}
}}
\section{growthPheno Version 1.0-34 (2021-10-14)}{
  \itemize{
   \item{Deprecate anomPlot, corrPlot, imagetimesPlot and longiPlots.}
   \item{Modified a ggplot function call to deal with a deprecation option.}
   \item{Add printPlot argument to plotCorrmatrix, plotImagetimes and plotMedianDeviations.}
   \item{Have the \code{ggplot} object returned for plotCorrmatrix, plotDeviationsBoxes and plotMedianDeviations.}
}}
\section{growthPheno Version 1.0-33 (2021-05-07)}{
  \itemize{
   \item{Move legend and to the bottom and turn it into a colour bar in plotCorrmatrix.}
}}
\section{growthPheno Version 1.0-32 (2021-05-06)}{
  \itemize{
   \item{Fix bug labelling y-axis variables in plotCorrmatrix.}
}}
\section{growthPheno Version 1.0-31 (2021-05-06)}{
  \itemize{
   \item{Add option \code{axis.text.size} to plotCorrmatrix for the \code{heatmap} and made a legend box for it.}
   \item{Added \code{ggplotFuncs} option to plotCorrmatrix.}
}}
\section{growthPheno Version 1.0-30 (2020-12-12)}{
  \itemize{
   \item{Fix a bug in splitValueCalculate associated with the arguments \code{which.obs} and \code{which.values}.}
   \item{Fix a bug for R 4.1.0 in longitudinalPrime, importExcel and plotMedianDeviations.}
}}
\section{growthPheno Version 1.0-29 (2020-11-10)}{
  \itemize{
   \item{Convert PVA into an S3 method with methods with PVA.data.frame and PVA.matrix. The argument \code{data} has been renamed to obj and placed first.}
   \item{Convert rcontrib into an S3 method with methods with rcontrib.data.frame and rcontrib.matrix. The argument \code{data} has been renamed to obj and placed first.}
   \item{Make intervalPVA into a method for \code{data.frame}s. The order of the arguments \code{data} and \code{responses} has been swapped in plotCorrmatrix.}
   \item{The argument \code{data} has been renamed to obj and placed first.}
}}
\section{growthPheno Version 1.0-28 (2020-07-25)}{
  \itemize{
   \item{Modify splitValueCalculate, and hence intervalValueCaculate, so that \code{which.obs} gives the observation number for the first observation that is closest to the value of the function results; \code{which.values} will return the value corresponding to \code{which.obs}.}
}}
\section{growthPheno Version 1.0-27 (2020-07-21)}{
  \itemize{
   \item{Fix bug in splitValueCalculate, and hence intervalValueCaculate, where \code{which.levels} set to other than \code{NULL} was giving an error. Replace the \code{which.levels} argument with the argument \code{which.values}.}
}}
\section{growthPheno Version 1.0-26 (2020-07-09)}{
  \itemize{
   \item{Fix bug in intervalWUI such the WUI was incorrectly calculated for \code{include.total.water = FALSE}.}
}}
\section{growthPheno Version 1.0-25 (2020-04-26)}{
  \itemize{
   \item{Fix bug in importExcel when calculated \code{TimeAfterPlaning} resulting from a timezone difference between the dates read in by readxl and that specified by \code{timeStart}. As a result the dat could be in error by one day for some observations.}
}}
\section{growthPheno Version 1.0-24 (2020-04-16)}{
  \itemize{
   \item{Fix bug in plotMedianDeviations when the plot involves a single DF.}
}}
\section{growthPheno Version 1.0-23 (2020-04-03)}{
  \itemize{
   \item{Add \code{x} to the call to plotMedianDeviations in probeSmoothing.}
   \item{Update name of  Rice vignette.}
   \item{Change the Maintainer email address.}
}}
\section{growthPheno Version 1.0-22 (2020-03-17)}{
  \itemize{
   \item{Revise Rice vignette and add Tomato vignette.}
   \item{Fix a bug for \code{which.plots = "compare"} in probeSmoothing.}
}}
\section{growthPheno Version 1.0-21 (2020-01-11)}{
  \itemize{
   \item{Corrrect link in DESCRIPTION and resubmit to CRAN.}
}}
\section{growthPheno Version 1.0-20 (2020-01-08)}{
  \itemize{
   \item{Put in a warning for the \code{addMEdianWhiskers} argument of plotLongitudinal for when the \code{x} and |code{xname} arguments are not equal; this can occur when \code{x} is a \code{formula} that is based on \code{xname}.}
}}
\section{growthPheno Version 1.0-19 (2019-11-28)}{
  \itemize{
   \item{Use package readxl instead of XLConnect to read Excel worksheets to avoid dependency on Java.}
}}
\section{growthPheno Version 1.0-18 (2019-09-30)}{
  \itemize{
   \item{Fix bug in passing individuals to plotMedianDeviations from probeSmoothing.}
}}
\section{growthPheno Version 1.0-17 (2019-07-19)}{
  \itemize{
   \item{Added options \code{methods+rawcompare} and \code{df+rawcompare} that allows comparison with and without the raw data.}
}}
\section{growthPheno Version 1.0-16 (2019-07-16)}{
  \itemize{
   \item{Fix bug in plotDeviationsBoxes that stopped boxplots being plotted.}
   \item{Modifed facet labels when \code{which.plot = "dfcompare"} to include the \code{smoothing.method}.}
}}
\section{growthPheno Version 1.0-15 (2019-07-11)}{
  \itemize{
   \item{Change Days to DAP in \code{tomato.dat}.}
   \item{Fix two bugs in probeSmoothing.}
}}
\section{growthPheno Version 1.0-14 (2019-06-05)}{
  \itemize{
   \item{Add Treatments column to tomato.dat.}
   \item{Modify probeSmoothing so that growth rates are caluclated and saved even when they are not plotted.}
}}
\section{growthPheno Version 1.0-13 (2019-05-12)}{
  \itemize{
    \item{Modify the Title and Description of the package, as well as the Overview of the package.}
    \item{Make a vignette using the rice data from Al-Tamimi.}
    \item{Include data from a tomato experiment.}
    \item{This version includes a major revamp of the plotting functions so that all of them begin with plot. The functions anomPlot, corrPlot, imagetimesPlot have been renamed to plotAnom, plotCorrmatrix and plotImagetimes. While the original function continue to work, they will be removed in future versions.}
    \item{The new functions plotDeviationsBoxes and plotMedianDeviations have been added. These two functions are not restricted to comparing observed with smoothed values generated using growthPheno functions.}
    \item{The new function probeSmoothing was created by extensively revising probeDF to facilitate comparison of different \code{smoothing.methods} and \code{df}, including plots of the median of deviations, the latter calculated as observed values minus smoothed values. The function probeDF has been retained for compatability purposes, but with the addition of a warning that probeDF is soft-deprecated. It will be removed in future versions.}
    \item{Replaced getDates with getTimesSubsets, deprecating getDates.}
    \item{Change the argument \code{smoothing.scale} to \code{smoothing.method}, replacing the option \code{identity} with \code{direct} in both fitSpline and splitSplines.}
}}
