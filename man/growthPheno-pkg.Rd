\newcommand{\packageVer}{\Sexpr[results=rd,stage=build]{utils::packageDescription("#1", fields="Version")}}
\newcommand{\packageDate}{\Sexpr[results=rd,stage=build]{utils::packageDescription("#1", fields="Date")}}
\name{growthPheno-pkg}
\alias{growthPheno-pkg}
\alias{growthPheno}
\docType{package}
\title{
\packageTitle{growthPheno}
}
\description{\packageDescription{growthPheno}

\bold{Version:} {\packageVer{growthPheno}}

\bold{Date:} \packageDate{growthPheno}
}


\section{Index}{
  The following list of functions does not include those that are soft-deprecated, i.e. those that have been available in previous versions of \code{growthPheno} but will be removed in futer versions. For an description of the use of the listed functions and vignettes that are available, see the Overview section below.
  \tabular{ll}{
    (i) Apex wrapper functions    \cr\cr
    \code{\link{traitSmooth}}
    \tab Obtain smooths for a trait by fitting spline \cr
    \tab functions and, having compared several smooths, \cr
    \tab allows one of them to be chosen and returned in \cr
    \tab a \code{data.frame}. \cr
    \code{\link{traitExtractFeatures}} 
    \tab Extract features, that are single-valued for each \cr
    \tab individual, from smoothed traits over time. \cr
\cr    
    (ii) Data\cr\cr
    
    \code{\link{exampleData}} 
    \tab A small data set to use in function examples. \cr
    \code{\link{RicePrepped.dat}} 
    \tab Prepped data from an experiment to investigate \cr
    \tab a rice germplasm panel. \cr
    \code{\link{RiceRaw.dat}} 
    \tab Data for an experiment to investigate a rice \cr
    \tab germplasm panel.\cr
    \code{\link{tomato.dat}} 
    \tab Longitudinal data for an experiment to investigate \cr
    \tab tomato response to mycorrhizal fungi and zinc. \cr
    \cr
    (iii) Plots\cr\cr
    
    \code{\link{plotAnom}} 
    \tab Identifies anomalous individuals and produces\cr
    \tab profile plots without them and with just them.\cr
    \code{\link{plotCorrmatrix}} 
    \tab Calculates and plots correlation matrices for a\cr
    \tab set of responses.\cr
    \code{\link{plotDeviationsBoxes}}
    \tab Produces boxplots of the deviations of the observed \cr 
    \tab values from the smoothed values over values of x. \cr
    \code{\link{plotImagetimes}} 
    \tab Plots the time within an interval versus the interval.\cr 
    \tab For example, the hour of the day carts are imaged\cr
    \tab against the days after planting (or some other\cr
    \tab number of days after an event).\cr
    \code{\link{plotMedianDeviations}} 
    \tab Calculates and plots the medians of the deviations \cr 
    \tab of the smoothed values from the observed values. \cr
    \code{\link{plotProfiles}} 
    \tab Produces profile plots of longitudinal data for a set \cr
    \tab of individuals. \cr
    \code{\link{probeSmooths}} 
    \tab Computes and compares, for a set of smoothing parameters, \cr
    \tab a response and the smooths of it, possibly along with \cr
    \tab growth rates calculated from the smooths.\cr
    \code{\link{plotSmoothsComparison}}
    \tab Plots  several sets of smoothed values for a response, \cr
    \tab possibly along with growth rates and optionally including \cr 
    \tab the unsmoothed values, as well as deviations boxplots. \cr
    \code{\link{plotSmoothsMedianDevns}} 
    \tab Calculates and plots the medians of the deviations \cr 
    \tab from the observed values of several sets for smoothed \cr
    \tab values stored in a \code{data.frame} in long format. \cr
    \cr
    (iv) Smoothing and calculation \cr 
    of growth rates and WUI \cr 
    for each individual (Indv) \cr\cr

    \code{\link{byIndv4Intvl_GRsAvg}}
    \tab Calculates the growth rates for a specified\cr
    \tab time interval for individuals in a \code{data.frame} in \cr
    \tab long format by taking weighted averages of growth \cr
    \tab rates for times within the interval.\cr
    \code{\link{byIndv4Intvl_GRsDiff}}
    \tab Calculates the growth rates for a specified \cr
    \tab time interval for individuals in a \code{data.frame} in \cr
    \tab long format by differencing the values for \cr
    \tab a response within the interval.\cr
    \code{\link{byIndv4Intvl_ValueCalc}}
    \tab Calculates a single value that is a \code{function} of \cr
    \tab the values of an individual for a response in a \cr 
    \tab \code{data.frame} in long format over a specified \cr
    \tab time interval. \cr
    \code{\link{byIndv4Intvl_WaterUse}}
    \tab Calculates, for a set of responses, water use traits \cr
    \tab (WU, WUR, WUI), and the AGR, over a specified time \cr
    \tab interval for each individual in a \code{data.frame} \cr
    \tab in long format.\cr
    \code{\link{byIndv4Times_GRsDiff}}
    \tab Adds, to a 'data.frame', the growth rates \cr 
    \tab calculated for consecutive times for individuals in \cr
    \tab a \code{data.frame} in long format by differencing \cr
    \tab response values.\cr
    \code{\link{byIndv4Times_SplinesGRs}}
    \tab For a response in a \code{data.frame} in long format, \cr 
    \tab computes, for a single set of smoothing parameters, \cr 
    \tab smooths of the response, possibly along with \cr 
    \tab growth rates calculated from the smooths. \cr 
    \code{\link{byIndv_ValueCalc}}
    \tab Applies a \code{function} to calculate a single value from \cr
    \tab an individual's values for a response in a \code{data.frame} in \cr
    \tab long format.\cr
    \code{\link{probeSmooths}} 
    \tab For a response in a \code{data.frame} in long format, computes and \cr
    \tab compares, for sets of smoothing parameters, smooths \cr 
    \tab of the response, possibly along with growth rates \cr
    \tab calculated from the smooths.\cr
    \code{\link{smoothSpline}} 
    \tab Fit a spline to smooth the relationship between a \cr
    \tab \code{response} and an \code{x} in a \code{data.frame}, \cr
    \tab optionally computing growth rates using derivatives.\cr
\cr
    (v) Data frame manipulation\cr\cr
    
    \code{\link{as.smooths.frame}} 
    \tab Forms a \code{\link{smooths.frame}} from a \code{data.frame}, \cr
    \tab ensuring that the correct columns are present. \cr
    \code{\link{designFactors}} 
    \tab Adds the factors and covariates for a blocked,\cr
    \tab split-unit design.\cr
    \code{\link{getTimesSubset}} 
    \tab Forms a subset of 'responses' in 'data' that\cr
    \tab contains their values for the nominated times.\cr
    \code{\link{importExcel}} 
    \tab Imports an Excel imaging file and allows some\cr
    \tab renaming of variables.\cr
    \code{\link{is.smooths.frame}} 
    \tab Tests whether an object is of class smooths.frame. \cr
    \code{\link{prepImageData}} 
    \tab Selects a set variables to be retained in a\cr
    \tab data frame of longitudinal data.\cr
    \code{\link{smooths.frame}} 
    \tab Description of a smooths.frame object, \cr
    \code{\link{twoLevelOpcreate}} 
    \tab Creates a data.frame formed by applying, for\cr
    \tab each response, a binary operation to the values of \cr
    \tab two different treatments.\cr
    \code{\link{validSmoothsFrame}} 
    \tab Checks that an object is a valid \code{\link{smooths.frame}}. \cr
    \cr

    (vi) General calculations \cr\cr
    
    \code{\link{anom}} 
    \tab Tests if any values in a vector are anomalous\cr
    \tab in being outside specified limits.\cr
    \code{\link{calcLagged}} 
    \tab Replaces the values in a vector with the result\cr
    \tab of applying an operation to it and a lagged value.\cr
    \code{\link{calcTimes}} 
    \tab Calculates for a set of times, the time intervals \cr
    \tab after an origin time and the position of each \cr 
    \tab within a time interval\cr
    \code{\link{cumulate}} 
    \tab Calculates the cumulative sum, ignoring the\cr
    \tab first element if exclude.1st is TRUE.\cr
    \code{\link{GrowthRates}} 
    \tab Calculates growth rates (AGR, PGR, RGRdiff) \cr
    \tab between a pair of values in a vector.\cr
    \code{\link{WUI}} 
    \tab Calculates the Water Use Index (WUI) for a value \cr 
    \tab of the response and of the water use.\cr
\cr
    (vii) Principal variates analysis (PVA)\cr\cr
    
    \code{\link{intervalPVA.data.frame}} 
    \tab Selects a subset of variables using PVA, based on \cr
    \tab the observed values within a specified time interval\cr
    \code{\link{PVA.data.frame}} 
    \tab Selects a subset of variables stored in a data.frame \cr
    \tab using PVA.\cr
    \code{\link{PVA.matrix}} 
    \tab Selects a subset of variables using PVA based on a \cr
    \tab correlation matrix. \cr
    \code{\link{rcontrib.data.frame}} 
    \tab Computes a measure of how correlated each\cr
    \tab variable in a set is with the other variable,\cr
    \tab conditional on a nominated subset of them.\cr
    \code{\link{rcontrib.matrix}} 
    \tab Computes a measure of how correlated each\cr
    \tab variable in a set is with the other variable,\cr
    \tab conditional on a nominated subset of them.\cr
  }}
\author{
\packageAuthor{growthPheno}

Maintainer: \packageMaintainer{growthPheno}
}
\keyword{ hplot }
\keyword{ manip }
\keyword{ package }
\seealso{\code{\link[dae:dae-package]{dae}}
}
\section{Overview}{
This package can be used to analyse growth data using splines to smooth the trend of individual plant traces over time and then to extract features or tertiarty traits for further analysis. This process is called smoothing and extraction 
of traits (SET) by Brien et al. (2020), who detail the use of `growthPheno` for carrying out the method. The two primary functions that implement this approach are \code{\link{traitSmooth}} and \code{\link{traitExtractFeatures}}. These may be the only functions that are used in that the complete SET process can be carried out with using only them. The \code{Tomato} vignette illustrates this for the example presented in Brien et al. (2020). 

The function \code{\link{traitSmooth}} utilizes the secondary functions \code{\link{probeSmooths}}, \code{\link{plotSmoothsComparison}} and code{\link{plotSmoothsMedianDevns}}. The function \code{\link{probeSmooths}} utilizes the tertiary functions \code{\link{byIndv4Times_SplinesGRs}} and \code{\link{byIndv4Times_GRsDiff}}, which in turn call the function \code{\link{smoothSpline}}. The function \code{\link{plotSmoothsComparison}} calls \code{\link{plotDeviationsBoxes}}. The primary functions \code{\link{traitExtractFeatures}} uses the secondary functions \code{\link{getTimesSubset}} and the set of \code{byIndv4Intvl_} functions. Recourse to these functions may be necessary for special cases. The \code{Rice} vignetter illustrates the use of the secondary and tertiary functions.

Use \code{vignette("Tomato", package = "growthPheno")} or \code{vignette("Rice", package = "growthPheno")} to access either of the vignettes.

In addition to functions that implement SET approach, `growthPheno`  also has functions for importing and orgainizing the data that are generally applicable, although they do have defaults that  make them particularly adapted to data from a high-throughput phenotyping facility based on a Lemna-Tec Scananalyzer 3D system. 

Data suitable for use with this package consists of columns of data obtained from a set of units (plant, pots, carts or plots) over time. There should be a unique identifier for each unit and a time variable, such as Days after Planting (DAP),  that contain no repeats for a unit. The combination of the identifier and a time for a unit should be unique to that unit. For imaging data, the carts/pots may be arranged in a grid of Lanes \eqn{\times}{x} Positions. 
}

\references{
Brien, C., Jewell, N., Garnett, T., Watts-Williams, S. J., & Berger, B. (2020). Smoothing and extraction of traits in the growth analysis of noninvasive phenotypic data. *Plant Methods*, **16**, 36. \doi{10.1186/s13007-020-00577-6}.
}
