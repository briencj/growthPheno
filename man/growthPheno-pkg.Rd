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
  For an overview of the use of these functions and an example see below.
  \tabular{ll}{
    (i) Data\cr\cr
    
    \code{\link{exampleData}} 
    \tab A small data set to use in function examples. \cr
    \code{\link{RiceRaw.dat}} 
    \tab Data for an experiment to investigate a rice \cr
    \tab germplasm panel.\cr
    \code{\link{tomato.dat}} 
    \tab Longitudinal data for an experiment to investigate \cr
    \tab tomato response to mycorrhizal fungi and zinc. \cr
    \cr
    (ii) Data frame manipulation\cr\cr
    
    \code{\link{designFactors}} 
    \tab Adds the factors and covariates for a blocked,\cr
    \tab split-plot design.\cr
    \code{\link{getTimesSubset}} 
    \tab Forms a subset of 'responses' in 'data' that\cr
    \tab contains their values for the nominated times.\cr
    \code{\link{importExcel}} 
    \tab Imports an Excel imaging file and allows some\cr
    \tab renaming of variables.\cr
    \code{\link{longitudinalPrime}} 
    \tab Selects a set variables to be retained in a\cr
    \tab data frame of longitudinal data.\cr
    \code{\link{twoLevelOpcreate}} 
    \tab Creates a data.frame formed by applying, for\cr
    \tab each response, abinary operation to the values of \cr
    \tab two different treatments.\cr
    \cr
    (iii) Plots\cr\cr
    
    \code{\link{plotAnom}} 
    \tab Identifies anomalous individuals and produces\cr
    \tab longitudinal plots without them and with just them.\cr
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
    \code{\link{plotLongitudinal}} 
    \tab Plots longitudinal data for a set of indiividuals \cr
    \code{\link{plotMedianDeviations}} 
    \tab Calculates and plots the median of the deviations \cr 
    \tab of the smoothed values from the observed values. \cr
    \code{\link{probeSmoothing}} 
    \tab Compares, for a set of specified values of df and \cr
    \tab different smoothing methods, a response and the smooths \cr
    \tab of it, possibly along with growth rates calculated \cr
    \tab from the smooths.\cr
    \cr
    (iii) Smoothing\cr\cr

    \code{\link{fitSpline}} 
    \tab Produce the fits from a natural cubic smoothing\cr
    \tab spline applied to a response in a 'data.frame', \cr
    \tab and growth rates can be computed using derivatives.\cr
    \code{\link{splitSplines}} 
    \tab Adds the fits, and optionally growth rates computed \cr
    \tab from derivatives, after fitting natural cubic\cr
    \tab smoothing splines to subsets of a response to a\cr
    \tab 'data.frame'.\cr
    \cr
    (iv) Growth rate and WUI calculation\cr\cr
    
    \code{\link{fitSpline}} 
    \tab Produce the fits from a natural cubic smoothing\cr
    \tab spline applied to a response in a 'data.frame', \cr
    \tab and growth rates can be computed using derivatives.\cr
    \code{\link{GrowthRates}} 
    \tab Calculates growth rates (AGR, PGR, RGRdiff) \cr
    \tab between pairs of values in a vector.\cr
    \code{\link{intervalGRaverage}} 
    \tab Calculates the growth rates for a specified\cr
    \tab time interval by taking weighted averages of\cr
    \tab growth rates for times within the interval.\cr
    \code{\link{intervalGRdiff}} 
    \tab Calculates the growth rates for a specified\cr
    \tab time interval.\cr
    \code{\link{splitContGRdiff}} 
    \tab Adds the growth rates calculated continuously\cr
    \tab over time for subsets of a response to a\cr
    \tab 'data.frame'.\cr
    \code{\link{splitSplines}} 
    \tab Adds the fits, and optionally growth rates computed \cr
    \tab from derivatives, after fitting natural cubic\cr
    \tab smoothing splines to subsets of a response to a\cr
    \tab 'data.frame'.\cr
    \code{\link{WUI}} 
    \tab Calculates the Water Use Index (WUI).\cr
    \code{\link{intervalWUI}} 
    \tab Calculates water use indices (WUI) over a\cr
    \tab specified time interval to a data.frame.\cr
    \cr
    (v) General calculations \cr\cr
    
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
    \code{\link{intervalValueCalculate}}
    \tab Calculates a single value that is a function of\cr
    \tab an individual's values for a response over a\cr
    \tab specified time interval.\cr
    \code{\link{splitValueCalculate}}   
    \tab Calculates a single value that is a function of\cr
    \tab an individual's values for a response.\cr
    \cr
    (vi) Principal variates analysis (PVA)\cr\cr
    
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
This package can be used to analyse growth data using splines to smooth the trend of individual plant traces over time and then to extract traits for further analysis. This process is called smoothing and extraction of traits (SET) by Brien et al. (2020), who detail the use of `growthPheno` for carrying out the method. 

The package `growthPheno` has tools that aid in choosing the degree of smoothing and the selection of traits. There are also functions for importing and orgainizing the data that are generally applicable, although they do have defaults that  make them particularly adapted to data from a high-throughput phenotyping facility based on a Lemna-Tec Scananalyzer 3D system. 

Data suitable for use with this package consists of columns of data obtained from a set of units (pots, carts or plots) over time. There should be a unique identifier for each unit, which by default is \code{Snapshot.ID.Tag}, and variable giving the Days after Planting for each measurement, by default {\code{Time.after.Planting..d.}}. In some cases, it is expected that there will be a column labelled \code{Snapshot.Time.Stamp}, which reflects the imaging time from which a particular data value was obtained. For imaging data, the carts/pots may be arranged in a grid of Lanes \eqn{\times}{x} Positions. 

The vignettes \code{Tomato} and \code{Rice} illustrate this process, the former being the example presented in Brien et al. (2020). Use \code{vignette("Tomato", package = "growthPheno")} or \code{vignette("Rice", package = "growthPheno")} to access either of the vignettes.

}

\references{
Brien, C., Jewell, N., Garnett, T., Watts-Williams, S. J., & Berger, B. (2020). Smoothing and extraction of traits in the growth analysis of noninvasive phenotypic data. *Plant Methods*, **16**, 36. <http://dx.doi.org/10.1186/s13007-020-00577-6>.
}
