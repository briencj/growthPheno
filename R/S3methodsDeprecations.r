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

longiPlot <- function(...)
{ .Deprecated(new = "plotLongitudinal", package = "growthPheno")
  invisible()
}

probeDF <- function(...)
{ .Deprecated(new = "probeSmooths", package = "growthPheno")
  invisible()
}
