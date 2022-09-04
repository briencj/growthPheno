#Tesat for the trait wrapper functions

cat("#### Test traitSmooth with small example\n")
test_that("exampleData_traitSmooth", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(growthPheno)
  
  data(exampleData)
  vline <- list(ggplot2::geom_vline(xintercept=29, linetype="longdash", size=1),
                ggplot2::scale_x_continuous(breaks=seq(28, 42, by=2)))
  yfacets <- c("Smarthouse", "Treatment.1")
  testthat::expect_message(
    smth.dat <- traitSmooth(data = longi.dat, 
                            response = "PSA", response.smoothed = "sPSA",
                            individuals = "Snapshot.ID.Tag", times = "DAP", 
                            keep.columns = yfacets, 
                            facet.y.pf = yfacets, 
                            facet.y.chosen = yfacets, 
                            addMediansWhiskers = TRUE, #used  whenever plotLongitudinal is used
                            ggplotFuncsProfile = vline))
  testthat::expect_equal(nrow(smth.dat), 280)
  testthat::expect_equal(ncol(smth.dat), 37)
  testthat::expect_true(all(names(longi.dat) %in% names(smth.dat)))
  testthat::expect_true(all(longi.dat$Snapshot.ID.Tag == smth.dat$Snapshot.ID.Tag))
  testthat::expect_true(all(c("Smarthouse","Treatment.1","PSA","PSA.AGR","PSA.RGR",
                              "sPSA","sPSA.AGR","sPSA.RGR") %in% names(smth.dat)))
})

cat("#### Test traitExtractFeatures with tomato example\n")
test_that("tomato_traitExtractFeatures", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(dae)
  library(growthPheno)
  
  data(tomato.dat)
  DAP.endpts   <- c(18,22,27,33,39,43,51)
  nDAP.endpts <- length(DAP.endpts)
  DAP.starts <- DAP.endpts[-nDAP.endpts]
  DAP.stops   <- DAP.endpts[-1]
  DAP.mids <- (DAP.starts + DAP.stops)/2
  DAP.segs <- list(c(DAP.endpts[1]-1, 39), 
                   c(40, DAP.endpts[nDAP.endpts]))
  #Add PSA rates and smooth PSA, also producing sPSA rates
  tom.dat <- byIndv4Times_SplinesGRs(data = tomato.dat, 
                                 response = "PSA", response.smoothed = "sPSA", 
                                 times = "DAP", rates.method = "differences", 
                                 smoothing.method = "log", 
                                 spline.type = "PS", lambda = 1, 
                                 smoothing.segments = DAP.segs)
  
  #Smooth WU
  tom.dat <- byIndv4Times_SplinesGRs(data = tom.dat, 
                                 response = "WU", response.smoothed = "sWU",
                                 rates.method = "none", 
                                 times = "DAP", 
                                 smoothing.method = "direct", 
                                 spline.type = "PS", lambda = 10^(-0.5), 
                                 smoothing.segments = DAP.segs)
  testthat::expect_equal(nrow(tom.dat), 1120)
  testthat::expect_equal(ncol(tom.dat), 20)
  
  ### Omit responses for the outlier plant
  omit <- with(tom.dat, Zn==90 & AMF=="+" & Block ==4)
  responses.all <- names(tom.dat)[match("Weight.After", names(tom.dat)):length(tom.dat)]
  tom.dat[responses.all] <- lapply(tom.dat[responses.all], 
                                   function(kcol, omit) 
                                   {
                                     kcol[omit] <- NA
                                     return(kcol)
                                   }, omit = omit)
  #Set up for individual traits
  indv.cols <- c("Snapshot.ID.Tag", "Lane", "Position", "Block", "Cart", "AMF", "Zn")
  indv.ini <- subset(tom.dat, subset = DAP == DAP.endpts[1], 
                     select = indv.cols)

  #'## Extract single-valued smoothed traits for each individual
  indv.dat <- traitExtractFeatures(data = tom.dat, 
                                   starts.intvl = DAP.starts, stops.intvl = DAP.stops, 
                                   responses4singletimes = "sPSA", 
                                   responses4intvl.rates = "sPSA", growth.rates = c("AGR", "RGR"), 
                                   water.use4intvl.traits = "sWU", 
                                   responses4water = "sPSA", 
                                   responses4overall.total = "sWU",
                                   responses4overall.max = "sPSA.AGR",
                                   mergedata = indv.ini)
  testthat::expect_equal(nrow(indv.dat), 32)
  testthat::expect_equal(ncol(indv.dat), 47)

  #'## Extract single-valued unsmoothed and smoothed traits in parallel for each individual
  indv.dat <- traitExtractFeatures(data = tom.dat, times = "DAP", 
                                   starts.intvl = DAP.starts, stops.intvl = DAP.stops, 
                                   responses4singletimes = c("PSA", "sPSA"), 
                                   responses4intvl.rates = c("PSA", "sPSA"), growth.rates = c("AGR", "RGR"), 
                                   water.use4intvl.traits = c("WU","sWU"), 
                                   responses4water = c("PSA","sPSA"),
                                   responses4overall.rates = c("PSA", "sPSA"),
                                   water.use4overall.water = c("WU","sWU"), 
                                   responses4overall.water = c("PSA","sPSA"),
                                   times.overall = c(18,51),
                                   mergedata = indv.ini)
  testthat::expect_equal(nrow(indv.dat), 32)
  testthat::expect_equal(ncol(indv.dat), 7 + (2*7) + (4*6) + (6*6) + 4 + 6) #91
  suffs <- paste(DAP.starts, DAP.stops, sep = "to")
  testthat::expect_true(all(names(indv.dat)[-(1:7)] == c(as.vector(outer(c("PSA","sPSA"), DAP.endpts, paste, sep = ".")),
                                                         as.vector(outer(c("PSA.AGR","PSA.RGR"), suffs, paste, sep = ".")),
                                                         as.vector(outer(c("sPSA.AGR","sPSA.RGR"), suffs, paste, sep = ".")),
                                                         as.vector(outer(c("WU","WUR","PSA.WUI"), suffs, paste, sep = ".")),
                                                         as.vector(outer(c("sWU","sWUR","sPSA.sWUI"), suffs, paste, sep = ".")),
                                                         "PSA.AGR","PSA.RGR","sPSA.AGR","sPSA.RGR","WU","WUR","PSA.WUI",
                                                         "sWU","sWUR","sPSA.sWUI")))

  #'## Extract single-valued unsmoothed and smoothed traits in parallel for each individual with "_" separator
  indv.dat <- traitExtractFeatures(data = tom.dat, times = "DAP", 
                                   starts.intvl = DAP.starts, stops.intvl = DAP.stops, 
                                   responses4singletimes = c("PSA", "sPSA"), 
                                   responses4intvl.rates = c("PSA", "sPSA"), growth.rates = c("AGR", "RGR"), 
                                   water.use4intvl.traits = c("WU","sWU"), 
                                   responses4water = c("PSA","sPSA"), 
                                   responses4overall.rates = c("PSA", "sPSA"),
                                   water.use4overall.water = c("WU","sWU"), 
                                   responses4overall.water = c("PSA","sPSA"),
                                   times.overall = c(18,51),
                                   sep.growth.rates = "_", sep.water.traits = "_", 
                                   sep.suffix.times = "_", sep.times.intvl = "_", 
                                   mergedata = indv.ini)
  testthat::expect_equal(nrow(indv.dat), 32)
  testthat::expect_equal(ncol(indv.dat), 7 + (2*7) + (4*6) + (6*6) + 4 + 6) #91
  suffs <- paste(DAP.starts, DAP.stops, sep = "_")
  testthat::expect_true(all(names(indv.dat)[-(1:7)] == c(as.vector(outer(c("PSA","sPSA"), DAP.endpts, paste, sep = "_")),
                                                         as.vector(outer(c("PSA_AGR","PSA_RGR"), suffs, paste, sep = "_")),
                                                         as.vector(outer(c("sPSA_AGR","sPSA_RGR"), suffs, paste, sep = "_")),
                                                         as.vector(outer(c("WU","WU_R","PSA_WU_I"), suffs, paste, sep = "_")),
                                                         as.vector(outer(c("sWU","sWU_R","sPSA_sWU_I"), suffs, paste, sep = "_")),
                                                         "PSA_AGR","PSA_RGR","sPSA_AGR","sPSA_RGR","WU","WU_R","PSA_WU_I",
                                                         "sWU","sWU_R","sPSA_sWU_I")))
  #Check the overall values
  testthat::expect_true(all((indv.dat[1, c("PSA_AGR","PSA_RGR","sPSA_AGR","sPSA_RGR","WU","WU_R","PSA_WU_I",
                                           "sWU","sWU_R","sPSA_sWU_I")] - 
                               c( 4.899273,0.08852807,4.897457,0.08655332,932,28.24242,0.1734721,
                                  921.4677,27.92326,0.1753898)) < 1e-04))
  
  
  #'## Extract single-valued unsmoothed and smoothed traits in parallel for each individual with no separator
  indv.dat <- traitExtractFeatures(data = tom.dat, times = "DAP", 
                                   starts.intvl = DAP.starts, stops.intvl = DAP.stops, 
                                   responses4singletimes = c("PSA", "sPSA"), 
                                   responses4intvl.rates = c("PSA", "sPSA"), growth.rates = c("AGR", "RGR"), 
                                   water.use4intvl.traits = c("WU","sWU"), 
                                   responses4water = c("PSA","sPSA"), 
                                   responses4overall.rates = c("PSA", "sPSA"),
                                   water.use4overall.water = c("WU","sWU"), 
                                   responses4overall.water = c("PSA","sPSA"),
                                   times.overall = c(18,51),
                                   sep.growth.rates = "", sep.water.traits = "", 
                                   sep.suffix.times = "", sep.times.intvl = "",
                                   mergedata = indv.ini)
  testthat::expect_equal(nrow(indv.dat), 32)
  testthat::expect_equal(ncol(indv.dat), 7 + (2*7) + (4*6) + (6*6) + 4 + 6) #91
  suffs <- paste(DAP.starts, DAP.stops, sep = "")
  testthat::expect_true(all(names(indv.dat)[-(1:7)] == c(as.vector(outer(c("PSA","sPSA"), DAP.endpts, paste, sep = "")),
                                                         as.vector(outer(c("PSAAGR","PSARGR"), suffs, paste, sep = "")),
                                                         as.vector(outer(c("sPSAAGR","sPSARGR"), suffs, paste, sep = "")),
                                                         as.vector(outer(c("WU","WUR","PSAWUI"), suffs, paste, sep = "")),
                                                         as.vector(outer(c("sWU","sWUR","sPSAsWUI"), suffs, paste, sep = "")),
                                                         "PSAAGR","PSARGR","sPSAAGR","sPSARGR","WU","WUR","PSAWUI",
                                                         "sWU","sWUR","sPSAsWUI")))
  
  #one AGR for sPSA and its overall AGR
  indv.dat <- traitExtractFeatures(data = tom.dat, times = "DAP", 
                                   starts.intvl = DAP.starts, stops.intvl = DAP.stops, 
                                   responses4intvl.rates = "sPSA",
                                   growth.rates = "AGR", 
                                   responses4overall.rates = "sPSA",
                                   times.overall = c(18,51),
                                   mergedata = indv.ini)
  testthat::expect_equal(nrow(indv.dat), 32)
  testthat::expect_equal(ncol(indv.dat), 14)
  
  #Overall values only for both unsmoothed and smoothed traits in parallel
  indv.dat <- traitExtractFeatures(data = tom.dat, times = "DAP", 
                                   growth.rates = c("AGR", "RGR"), 
                                   responses4overall.rates = c("PSA", "sPSA"),
                                   water.use4overall.water = c("WU","sWU"), 
                                   responses4overall.water = c("PSA","sPSA"),
                                   times.overall = c(18,51),
                                   mergedata = indv.ini)
  #Check the overall values
  testthat::expect_true(all((indv.dat[1, c("PSA.AGR","PSA.RGR","sPSA.AGR","sPSA.RGR","WU","WUR","PSA.WUI",
                                           "sWU","sWUR","sPSA.sWUI")] - 
                               c( 4.899273,0.08852807,4.897457,0.08655332,932,28.24242,0.1734721,
                                  921.4677,27.92326,0.1753898)) < 1e-04))
  testthat::expect_equal(nrow(indv.dat), 32)
  testthat::expect_equal(ncol(indv.dat), 17)
  
  #Overall values only for smoothed traits
  testthat::expect_error(indv.diff.dat <- traitExtractFeatures(data = tom.dat, times = "DAP", 
                                                               responses4overall.rates = "sPSA",
                                                               water.use4overall.water = "sWU", 
                                                               responses4overall.water = "sPSA",
                                                               times.overall = c(18,51),
                                                               mergedata = indv.ini),
                         regexp = "growth.rates needs to be set for responses4overall.rates")
  
  indv.diff.dat <- traitExtractFeatures(data = tom.dat, times = "DAP", 
                                        growth.rates = "AGR", 
                                        responses4overall.rates = "sPSA",
                                        water.use4overall.water = "sWU", 
                                        responses4overall.water = "sPSA",
                                        times.overall = c(18,51),
                                        mergedata = indv.ini)
  testthat::expect_equal(nrow(indv.diff.dat), 32)
  testthat::expect_equal(ncol(indv.diff.dat), 11)
  
  #only overall water traits
  indv.diff.dat <- traitExtractFeatures(data = tom.dat, times = "DAP", 
                                        water.use4overall.water = "sWU", 
                                        responses4overall.water = "sPSA",
                                        times.overall = c(18,51),
                                        mergedata = indv.ini)
  testthat::expect_equal(nrow(indv.diff.dat), 32)
  testthat::expect_equal(ncol(indv.diff.dat), 10)
  
  
  #Overall values only for unsmoothed and smoothed traits in parallel using ratesaverage
  testthat::expect_silent(
    indv.dat <- traitExtractFeatures(data = tom.dat, times = "DAP", 
                                     growth.rates = c("AGR", "RGR"), rates.method = "ratesaverage",
                                     responses4overall.rates = c("PSA", "sPSA"),
                                     water.use4overall.water = c("WU","sWU"), 
                                     responses4overall.water = c("PSA","sPSA"),
                                     times.overall = c(18,51),
                                     mergedata = indv.ini))
  testthat::expect_equal(nrow(indv.dat), 32)
  testthat::expect_equal(ncol(indv.dat), 17)
  
  #Overall values only for smoothed traits using ratesaverage
  indv.dat <- traitExtractFeatures(data = tom.dat, times = "DAP", 
                                   starts.intvl = DAP.starts, stops.intvl = DAP.stops, 
                                   responses.rates = "sPSA",
                                   growth.rates = "AGR", rates.method = "ratesaverage",
                                   responses4overall.rates = "sPSA",
                                   water.use4overall.water = "sWU", 
                                   responses4overall.water = "sPSA",
                                   times.overall = c(18,51),
                                   mergedata = indv.ini)
  testthat::expect_equal(nrow(indv.dat), 32)
  testthat::expect_equal(ncol(indv.dat), 11)
  
  #Check the overall values
  indv.dat <- traitExtractFeatures(data = tom.dat, times = "DAP", 
                                   growth.rates = c("AGR", "RGR"), rates.method = "ratesaverage",
                                   responses4overall.rates = c("PSA","sPSA"),
                                   water.use4overall.water = c("WU","sWU"), 
                                   responses4overall.water = c("PSA","sPSA"),
                                   times.overall = c(18,51),
                                   mergedata = indv.ini)
  testthat::expect_true(all((indv.dat[1, c("PSA.AGR","PSA.RGR","sPSA.AGR","sPSA.RGR","WU","WUR","PSA.WUI",
                                           "sWU","sWUR","sPSA.sWUI")] - 
                               c( 4.899273,0.08852807,4.897457,0.08655332,932,28.24242,0.1734721,
                                  921.4677,27.92326,0.1753898)) < 1e-04))
  
})
