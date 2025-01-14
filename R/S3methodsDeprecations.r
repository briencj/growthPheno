#S3 method generics

intervalPVA <- function(obj, ...) UseMethod("intervalPVA")
PVA <- function(obj, ...) UseMethod("PVA")
rcontrib <- function(obj, ...) UseMethod("rcontrib")

#Deprecations

anomPlot <- function(...)
{ .Deprecated(new = "plotAnom", package = "growthPheno")
  invisible()
}

corrPlot <- function(...)
{ .Deprecated(new = "plotCorrmatrix", package = "growthPheno")
  invisible()
}

fitSpline <- function(...)
{ .Deprecated(new = "smoothSpline", package = "growthPheno")
  invisible()
}

getDates <- function(...)
{ .Deprecated(new = "getTimesSubset", package = "growthPheno")
  invisible()
}

getImageData <- function(...)
{ .Deprecated(new = "prepImageData", package = "growthPheno")
  invisible()
}


imagetimesPlot <- function(...)
{ .Deprecated(new = "plotImagetimes", package = "growthPheno")
  invisible()
}

intervalGRaverage <- function(...)
{
  .Deprecated(new = "byIndv4Intvl_GRsAvg", package = "growthPheno")
  invisible()
}

intervalGRdiff <- function(...)
{
  .Deprecated(new = "byIndv4Intvl_GRsDiff", package = "growthPheno")
  invisible()
}

intervalValueCalculate <- function(...)
{
  .Deprecated(new = "byIndv4Intvl_ValueCalc", package = "growthPheno")  
  invisible()
}

intervalWUI <- function(...)
{
  .Deprecated(new = "byIndv4Intvl_WaterUse", package = "growthPheno")
  invisible()
}

longiPlot <- function(...)
{ .Deprecated(new = "plotLongitudinal", package = "growthPheno")
  invisible()
}

longitudinalPrime <- function(...)
{
  .Deprecated(new = "prepImageData", package = "growthPheno")
  invisible()
}

plotLongitudinal <- function(...)
{
  .Deprecated(new = "plotProfiles", package = "growthPheno")
  invisible()
}

plotMedianDeviations <- function(...)
{
  .Deprecated(new = "plotSmoothsMedianDevns", package = "growthPheno")
  invisible()
}

probeDF <- function(...)
{ .Deprecated(new = "probeSmooths", package = "growthPheno")
  invisible()
}

probeSmoothing <- function(...)
{
  .Deprecated(new = "probeSmooths", package = "growthPheno")
  invisible()
}

splitContGRdiff <- function(...)
{
  .Deprecated(new = "byIndv4Times_GRsDiff", package = "growthPheno")
  invisible()
}

splitSplines <- function(...)
{
  .Deprecated(new = "byIndv4Times_SplinesGRs", package = "growthPheno")
  invisible()
}

splitValueCalculate <- function(...)
{
  .Deprecated(new = "byIndv4Intvl_ValueCalc", package = "growthPheno")
  invisible()
}
