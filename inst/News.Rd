\name{news}
\title{News for Package growthPheno}
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
