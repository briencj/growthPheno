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
  
  #'## Extract single-valued traits for each individual
  indv.cols <- c("Snapshot.ID.Tag", "Lane", "Position", "Block", "Cart", "AMF", "Zn")
  indv.dat <- subset(tom.dat, subset = DAP == DAP.endpts[1], 
                     select = indv.cols)
  indv.dat <- traitExtractFeatures(data = tom.dat, 
                                   starts.intvl = DAP.starts, stops.intvl = DAP.stops, 
                                   responses.singletimes = "sPSA", 
                                   responses.rates = "sPSA", growth.rates = c("AGR", "RGR"), 
                                   water.use = "sWU", responses.water = "sPSA", 
                                   responses.total = "sWU",
                                   responses.max = "sPSA.AGR",
                                   mergedata = indv.dat)
  testthat::expect_equal(nrow(indv.dat), 32)
  testthat::expect_equal(ncol(indv.dat), 47)
  
})