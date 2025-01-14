\name{growthPheno-deprecated}
\alias{growthPheno-deprecated}
\alias{anomPlot}
\alias{corrPlot}
\alias{getDates}
\alias{imagetimesPlot}
\alias{longiPlot}
\alias{probeDF}
\alias{fitSpline}
\alias{intervalGRaverage}
\alias{intervalGRdiff}
\alias{intervalValueCalculate}
\alias{intervalWUI}
\alias{longitudinalPrime}
\alias{plotLongitudinal}
\alias{plotMedianDeviations}
\alias{probeSmoothing}
\alias{splitContGRdiff}
\alias{splitSplines}
\alias{splitValueCalculate}
\title{Deprecated Functions in the Package growthPheno}
\description{These functions have been renamed and deprecated in \code{growthPheno}: 
    \enumerate{
      \item getDates -> \code{\link{getTimesSubset}}
      \item anomPlot -> \code{\link{plotAnom}}
      \item corrPlot -> \code{\link{plotCorrmatrix}}
      \item fitspline -> \code{\link{smoothSpline}}
      \item imagetimesPlot -> \code{\link{plotImagetimes}}
      \item intervalGRaverage -> \code{\link{byIndv4Intvl_GRsAvg}}
      \item intervalGRdiff -> \code{\link{byIndv4Intvl_GRsDiff}}
      \item intervalValueCalculate -> \code{\link{byIndv4Intvl_ValueCalc}}
      \item intervalWUI -> \code{\link{byIndv4Intvl_WaterUse}}
      \item longiPlot -> \code{\link{plotProfiles}}
      \item longitudinalPrime -> \code{\link{prepImageData}}
      \item plotLongitudinal -> \code{\link{plotProfiles}}
      \item plotMedianDeviations -> \code{\link{plotSmoothsMedianDevns}}
      \item probeDF -> \code{\link{probeSmooths}}
      \item probeSmoothing -> \code{\link{probeSmooths}}
      \item splitContGRdiff -> \code{\link{byIndv4Times_GRsDiff}}
      \item splitSplines -> \code{\link{byIndv4Times_SplinesGRs}}
      \item splitValueCalculate -> \code{\link{byIndv4Intvl_ValueCalc}}
}} 
\usage{getDates(...)
anomPlot(...)
corrPlot(...)
fitSpline(...)
imagetimesPlot(...)
intervalGRdiff(...)
intervalGRaverage(...)
intervalValueCalculate(...)
intervalWUI(...)
longiPlot(...)
longitudinalPrime(...)
plotLongitudinal(...)
plotMedianDeviations(...)
probeDF(...)
probeSmoothing(...)
splitContGRdiff(...)
splitSplines(...)
splitValueCalculate(...)
}
\arguments{
  \item{...}{absorbs arguments passed from the old functions of the style foo.bar().}
}
\author{Chris Brien}
