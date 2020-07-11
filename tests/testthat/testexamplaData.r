
cat("#### Test \\dontrun examples for exampleData\n")
test_that("exampleData_growthPheno", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(growthPheno)
  library(ggplot2)
  
  #Load the data
  data(exampleData)
  
  #plotImagetimes
  longi.dat <- calcTimes(longi.dat, imageTimes = "Snapshot.Time.Stamp",
                         timePositions = "Hour")
  testthat::expect_message(
    plotImagetimes(data = longi.dat, intervals = "Days", timePositions = "Hour",
                   ggplotFuncs=list(scale_colour_gradient(low="grey20", high="black"), 
                                    geom_line(aes(group=Snapshot.ID.Tag, colour=Lane)))))
  
  #longitudinalPrime
  longiPrime.dat <- longitudinalPrime(data=raw.dat, smarthouse.lev=1)
  testthat::expect_equal(ncol(longiPrime.dat), 36)
  
  longiPrime.dat <- longitudinalPrime(data=raw.dat, smarthouse.lev=1, 
                                      traits = list(a = "Area", c = "Compactness"),
                                      labsCamerasViews = list(all = c("SV1", "SV2", "TV"), 
                                                              t = "TV"))
  testthat::expect_equal(ncol(longiPrime.dat), 21)
  
  longiPrime.dat <- longitudinalPrime(data=raw.dat, smarthouse.lev=1, 
                                        traits = c("Area.SV1", "Area.SV2", "Area.TV", 
                                                   "Compactness.TV"),
                                        labsCamerasViews = NULL)
  testthat::expect_equal(ncol(longiPrime.dat), 21)
  
  longiPrime.dat <- longitudinalPrime(data=raw.dat, smarthouse.lev=1, 
                                      calcWaterLoss = FALSE, 
                                      traits = list(img = c("Area", "Compactness"), 
                                                    H20 = c("Weight.Before","Weight.After",
                                                            "Water.Amount")),
                                      labsCamerasViews = list(all = c("SV1", "SV2", "TV"), 
                                                              H2O = NULL))
  testthat::expect_equal(ncol(longiPrime.dat), 22)
  
  #plotLongitudinal
  testthat::expect_silent(  
    plotLongitudinal(data = longi.dat, response = "Area.smooth"))
  
  testthat::expect_silent(  
    plt <- plotLongitudinal(data = longi.dat, response = "Area.smooth", x.title = "DAP",  
                            y.title = "Area.smooth", x="xDays+35.42857143", printPlot=FALSE))
  testthat::expect_silent(  
    plt <- plt + ggplot2::geom_vline(xintercept=29, linetype="longdash", size=1) +
      ggplot2::scale_x_continuous(breaks=seq(28, 42, by=2)) + 
      ggplot2::scale_y_continuous(limits=c(0,750)))
  testthat::expect_silent(  
    print(plt))
    
  testthat::expect_silent(  
    plotLongitudinal(data = longi.dat, response = "Area.smooth", x.title = "DAP",  
                     y.title = "Area.smooth", x="xDays+35.42857143", 
                     ggplotFuncs = list(ggplot2::geom_vline(xintercept=29, linetype="longdash", 
                                                            size=1), 
                                        ggplot2::scale_x_continuous(breaks=seq(28, 42, by=2)), 
                                        ggplot2::scale_y_continuous(limits=c(0,750)))))
    
  #plotAnom
  testthat::expect_silent(  
    anomalous <- plotAnom(longi.dat, response="Area.smooth.AGR", 
                          lower=2.5, start.time=40, 
                          x = "xDays+35.42857143", vertical.line=29, 
                          breaks=seq(28, 42, by=2), 
                          whichPrint=c("innerPlot"), 
                          y.title="Area.smooth.AGR"))
  
  #plotDeviationsBoxes
  testthat::expect_silent(  
    plotDeviationsBoxes(longi.dat, observed = "Area", smoothed = "Area.smooth",
                        x.factor="Days", facet.x = ".", facet.y= ".", df =5))
  
  
  #plotMedianDeviations
  vline <- list(ggplot2::geom_vline(xintercept=20, linetype="longdash", size=1),
                ggplot2::scale_x_continuous(breaks=seq(12, 36, by=2)))
  testthat::expect_silent(  
    traits <- probeSmoothing(data = longi.dat, response = "Area", 
                             df = c(4:7), x="xDays+24.16666667", 
                             facet.x = ".", facet.y = ".",
                             which.plots = "none",
                             deviations.plots = "none", 
                             propn.types = NULL))
  testthat::expect_silent(  
    med <- plotMedianDeviations(data = traits, 
                                response = "Area", response.smoothed = "Area.smooth", 
                                x="xDays+24.16666667", xname = "xDays", 
                                df = c(4,7), x.title = "DAP", 
                                facet.x = ".", facet.y = ".",
                                trait.types = "response", propn.types = 0.05,
                                ggplotFuncsMedDevn = vline))
  
  #probeSmoothing
  testthat::expect_silent(  
    probeSmoothing(data = longi.dat, response = "Area", df = c(4,7), x="xDays+24.16666667", 
                   ggplotFuncs=vline))
  
  #fitSpline
  testthat::expect_silent(  
    fit <- fitSpline(longi.dat, response="Area", , x="xDays", df = 4,
                     deriv=c(1,2), suffices.deriv=c("AGRdv","Acc")))
  
  #splitContGRdiff
  testthat::expect_silent(  
    tmp <- splitContGRdiff(longi.dat, response="Area.smooth", 
                           INDICES = "Snapshot.ID.Tag", which.rates=c("AGR", "RGR")))

    #intervalGRaverage
  testthat::expect_silent(  
    tmp <- splitSplines(longi.dat, response="Area", x="xDays", 
                        INDICES = "Snapshot.ID.Tag", 
                        df = 4, deriv=1, suffices.deriv="AGRdv", RGR="RGRdv"))
  testthat::expect_silent(  
    Area.smooth.GR <- intervalGRaverage("Area.smooth", which.rates = c("AGR","RGR"), 
                                        suffices.rates = c("AGRdv","RGRdv"), 
                                        start.time = 31, end.time = 35, 
                                        suffix.interval = "31to35",
                                        data = tmp))
  
  #twoLevelOpcreate
  responses <- c("Area.smooth.AGR","Area.smooth.RGR")
  cols.retained <-  c("Snapshot.ID.Tag","Smarthouse","Lane","Position",
                      "Days","Snapshot.Time.Stamp", "Hour", "xDays",
                      "Zones","xZones","SHZones","ZLane","ZMainplots",
                      "xMainPosn", "Genotype.ID")
  longi.SIIT.dat <- 
    twoLevelOpcreate(responses, longi.dat, suffices.treatment=c("C","S"),
                     operations = c("-", "/"), 
                     suffices.results = c("diff", "SIIT"), 
                     columns.retained = cols.retained, 
                     by = c("Smarthouse","Zones","ZMainplots","Days"))
  longi.SIIT.dat <- with(longi.SIIT.dat, 
                         longi.SIIT.dat[order(Smarthouse,Zones,ZMainplots,Days),])
  testthat::expect_equal(ncol(longi.SIIT.dat), 21)  
  testthat::expect_true("Area.smooth.RGR.SIIT" %in% names(longi.SIIT.dat))
  testthat::expect_true(is.na(longi.SIIT.dat$Area.smooth.RGR.SIIT[1]) && 
                          abs(longi.SIIT.dat$Area.smooth.RGR.SIIT[2] - 0.854679) < 1e-03)
})