#S3 method generics

intervalPVA <- function(obj, ...) UseMethod("intervalPVA")
PVA <- function(obj, ...) UseMethod("PVA")
rcontrib <- function(obj, ...) UseMethod("rcontrib")

#Deprecations

getDates <- function(...)
{ .Deprecated(new = "getTimesSubset", package = "growthPheno")
  invisible()
}

anomPlot <- function(...)
{ .Deprecated(new = "plotAnom", package = "growthPheno")
  invisible()
}

corrPlot <- function(...)
{ .Deprecated(new = "plotCorrmatrix", package = "growthPheno")
  invisible()
}

imagetimesPlot <- function(...)
{ .Deprecated(new = "plotImagetimes", package = "growthPheno")
  invisible()
}

longiPlot <- function(...)
{ .Deprecated(new = "plotLongitudinal", package = "growthPheno")
  invisible()
}

probeDF <- function(...)
{ .Deprecated(new = "probeSmooths", package = "growthPheno")
  invisible()
}
