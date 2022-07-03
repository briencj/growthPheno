
cat("#### Test \\dontrun examples for exampleData\n")
test_that("exampleData_growthPheno", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(growthPheno)
  library(ggplot2)
  
  #Load the data
  data(exampleData)
  
  #longi.dat
  longit.dat <- prepImageData(data=raw.dat, smarthouse.lev=1)
  testthat::expect_equal(ncol(longit.dat), 35)
  
  longit.dat <- prepImageData(data=raw.dat, smarthouse.lev=1, 
                              traits = list(a = "Area", c = "Compactness"),
                              labsCamerasViews = list(all = c("SV1", "SV2", "TV"), 
                                                      t = "TV"))
  testthat::expect_equal(ncol(longit.dat), 20)
  
  longit.dat <- prepImageData(data=raw.dat, smarthouse.lev=1, 
                              traits = c("Area.SV1", "Area.SV2", "Area.TV", 
                                         "Compactness.TV"),
                              labsCamerasViews = NULL)
  testthat::expect_equal(ncol(longit.dat), 20)
  
  longit.dat <- prepImageData(data=raw.dat, smarthouse.lev=1, 
                              calcWaterUse = FALSE, 
                              traits = list(img = c("Area", "Compactness"), 
                                            H20 = c("Weight.Before","Weight.After",
                                                    "Water.Amount")),
                              labsCamerasViews = list(all = c("SV1", "SV2", "TV"), 
                                                      H2O = NULL))
  testthat::expect_equal(ncol(longit.dat), 21)

  #plotImagetimes
  longit.dat <- calcTimes(longit.dat, imageTimes = "Snapshot.Time.Stamp",
                         timePositions = "Hour")
  testthat::expect_equal(ncol(longit.dat), 21)
  
  testthat::expect_message(
    plotImagetimes(data = longi.dat, intervals = "DAP", timePositions = "Hour",
                   ggplotFuncs=list(scale_colour_gradient(low="grey20", high="black"), 
                                    geom_line(aes(group=Snapshot.ID.Tag, colour=Lane)))))
  
  #plotProfiles
  testthat::expect_true("sPSA" %in% names(longi.dat))
  testthat::expect_silent(  
    plotProfiles(data = longi.dat, response = "sPSA"))
  
  testthat::expect_silent(  
    plt <- plotProfiles(data = longi.dat, response = "sPSA", times = "DAP", 
                        y.title = "sPSA (kpixels)", 
                        facet.x = "Treatment.1", facet.y = "Smarthouse", 
                        printPlot=FALSE))
  testthat::expect_silent(  
    plt <- plt + ggplot2::geom_vline(xintercept=29, linetype="longdash", size=1) +
      ggplot2::scale_x_continuous(breaks=seq(28, 42, by=2)) + 
      ggplot2::scale_y_continuous(limits=c(0,750)))
  testthat::expect_silent(  
    print(plt))
    
  testthat::expect_silent(  
    plotProfiles(data = longi.dat, response = "sPSA", times = "DAP", 
                     y.title = "sPSA (kpixels)",  
                     facet.x = "Treatment.1", facet.y = "Smarthouse", 
                     ggplotFuncs = list(ggplot2::geom_vline(xintercept=29, linetype="longdash", 
                                                            size=1), 
                                        ggplot2::scale_x_continuous(breaks=seq(28, 42, by=2)), 
                                        ggplot2::scale_y_continuous(limits=c(0,750)))))
    
  #plotAnom
  testthat::expect_silent(  
    anomalous <- plotAnom(longi.dat, response="sPSA.AGR", times = "DAP", 
                          lower=2.5, start.time=40, 
                          vertical.line=29, 
                          breaks=seq(28, 42, by=2), 
                          whichPrint=c("innerPlot"), 
                          y.title="sPSA AGR (kpixels)"))
  
  #plotDeviationsBoxes
  testthat::expect_silent(  
    plotDeviationsBoxes(longi.dat, observed = "PSA", smoothed = "sPSA",
                        x.factor="DAP", df =5))

  #probeSmooths
  vline <- list(ggplot2::geom_vline(xintercept=29, linetype="longdash", size=1))
  testthat::expect_silent(tmp <- probeSmooths(data = longi.dat, 
                                              response = "PSA", times = "DAP", 
                                              df = c(4,7), 
                                              facet.x.pf = "Tuning", facet.y.pf = "Treatment.1",
                                              ggplotFuncsProfile = vline))
  testthat::expect_equal(nrow(tmp), 560)
  testthat::expect_equal(ncol(tmp), 15)
  testthat::expect_silent(  
    traits <- probeSmooths(data = longi.dat, response = "PSA", 
                           individuals = "Snapshot.ID.Tag", 
                           df = c(4:7), times = "DAP", keep.columns = "Treatment.1", 
                           which.plots = "none",
                           propn.types = NULL))
  testthat::expect_silent(  
    med <- plotSmoothsMedianDevns(data = traits, 
                                  response = "PSA", response.smoothed = "sPSA", 
                                  times = "DAP", 
                                  x.title = "DAP", 
                                  trait.types = "response", propn.types = 0.05,
                                  plots.group.med = "Tuning", facet.y.med = "Treatment.1",
                                  ggplotFuncsMedDevn = vline))
  
  #fitSpline
  testthat::expect_silent(  
    fit <- fitSpline(longi.dat, response = "PSA", response.smoothed = "sPSA", 
                     x="xDAP", df = 4,
                     deriv=c(1,2), suffices.deriv=c("AGRdv","Acc")))
  testthat::expect_equal(nrow(fit$predictions), 280)
  testthat::expect_equal(ncol(fit$predictions), 4)
  testthat::expect_true(all(c("xDAP","sPSA","sPSA.AGRdv","sPSA.Acc") 
                            %in% names(fit$predictions)))
  
  #byIndv4Times_GRsDiff
  testthat::expect_silent(  
     tmp <- byIndv4Times_GRsDiff(data = longi.dat, response="sPSA", 
                                 which.rates=c("AGR", "RGR")))

  #byIndv4Times_SplinesGRs
  testthat::expect_silent(  
    tmp <- byIndv4Times_SplinesGRs(data = longi.dat, response="PSA", times="DAP", 
                                   df = 4, rates.method = "deriv", 
                                   which.rates = c("AGR", "RGR"), 
                                   suffices.rates = c("AGRdv", "RGRdv")))
  testthat::expect_silent(  
    sPSA.GR <- byIndv4Intvl_GRsAvg(data = tmp, responses = "sPSA", 
                                          which.rates = c("AGR","RGR"), 
                                          suffices.rates = c("AGRdv","RGRdv"), 
                                          start.time = 31, end.time = 35, 
                                          suffix.interval = "31to35"))
  
  #twoLevelOpcreate
  responses <- c("sPSA.AGR","sPSA.RGR")
  cols.retained <-  c("Snapshot.ID.Tag","Smarthouse","Lane","Position",
                      "DAP","Snapshot.Time.Stamp", "Hour", "xDAP",
                      "Zone","cZone","SHZone","ZLane","ZMainunit",
                      "cMainPosn", "Genotype.ID")
  longi.SIIT.dat <- 
    twoLevelOpcreate(dat = longi.dat, responses = responses, 
                     suffices.treatment=c("C","S"),
                     operations = c("-", "/"), 
                     suffices.results = c("diff", "SIIT"), 
                     columns.retained = cols.retained, 
                     by = c("Smarthouse","Zone","ZMainunit","DAP"))
  longi.SIIT.dat <- with(longi.SIIT.dat, 
                         longi.SIIT.dat[order(Smarthouse,Zone,ZMainunit,DAP),])
  testthat::expect_equal(ncol(longi.SIIT.dat), 21)  
  testthat::expect_true("sPSA.RGR.SIIT" %in% names(longi.SIIT.dat))
  testthat::expect_true(is.na(longi.SIIT.dat$sPSA.RGR.SIIT[1]) && 
                          abs(longi.SIIT.dat$sPSA.RGR.SIIT[2] - 0.854679) < 1e-03)
})

cat("#### Test byIndv_ValueCalc for exampleData\n")
test_that("exampleData_byIndv_ValueCalc", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(growthPheno)
 
  #Load the data
  data(exampleData)
  
  #Test byIndv_ValueCalc
  tmp <- byIndv_ValueCalc(data=longi.dat, response = "sPSA.AGR", FUN="max")
  testthat::expect_equal(ncol(tmp), 2)  
  testthat::expect_equal(nrow(tmp), 20)  
  testthat::expect_true(abs(tmp[1,2] - 42.82888) < 1e-03)  
  testthat::expect_true("sPSA.AGR.max" %in% names(tmp))
  tmp <- byIndv_ValueCalc(data=longi.dat, response = "sPSA.AGR", FUN="max", which.values = "DAP")
  testthat::expect_equal(ncol(tmp), 3)  
  testthat::expect_equal(nrow(tmp), 20)  
  testthat::expect_true(tmp[1,3] ==  42)  
  testthat::expect_true("sPSA.AGR.max.DAP" %in% names(tmp))
  tmp <- byIndv_ValueCalc(data=longi.dat, response = "sPSA.AGR", FUN="max", which.obs = T)
  testthat::expect_equal(ncol(tmp), 3)  
  testthat::expect_equal(nrow(tmp), 20)  
  testthat::expect_equal(tmp[1,3], 14)  
  testthat::expect_true("sPSA.AGR.max.obs" %in% names(tmp))
  tmp <- byIndv_ValueCalc(data=longi.dat, response = "sPSA.AGR", FUN="max", 
                          which.values = "DAP", which.obs = T)
  testthat::expect_equal(ncol(tmp), 4)  
  testthat::expect_equal(nrow(tmp), 20)  
  testthat::expect_true(all(c("sPSA.AGR.max.DAP",
                              "sPSA.AGR.max.obs") %in% names(tmp)))
  #Test numeric which.values
  tmp <- byIndv_ValueCalc(data=longi.dat, response = "sPSA.AGR", FUN="max", 
                          which.values = "xDAP", which.obs = T)
  testthat::expect_equal(ncol(tmp), 4)  
  testthat::expect_equal(nrow(tmp), 20)  
  testthat::expect_true(all(c("sPSA.AGR.max.xDAP",
                              "sPSA.AGR.max.obs") %in% names(tmp)))
  testthat::expect_true(abs(tmp[1,2] - 42.82888) < 1e-03)  
  
  #Test quantile that involves the probs argument and is not an exact observed value
  tmp <- byIndv_ValueCalc(data=longi.dat, response = "sPSA", FUN="quantile", 
                          which.values = "DAP", probs = 0.1)
  tmp <- cbind(tmp, longi.dat$sPSA[longi.dat$DAP == 28], 
               longi.dat$sPSA[longi.dat$DAP == 30])
  names(tmp)[4:5] <- c("sPSA.28", "sPSA.30")
  testthat::expect_equal(nrow(tmp), 20)  
  #CHeck that sPSA.quantile is closer to sPSA.30 than to sPSA.28 
  testthat::expect_true(all(with(tmp, (sPSA.28 - sPSA.quantile) <
                                   (sPSA.30 - sPSA.quantile))))
  
  #Test byIndv4Intvl_ValueCalc
  tmp <- byIndv4Intvl_ValueCalc(data=longi.dat, response = "sPSA.AGR", FUN="max", which.obs=T)
  testthat::expect_equal(ncol(tmp), 3)  
  testthat::expect_equal(nrow(tmp), 20)  
  testthat::expect_true(tmp[1,3] ==  14)  
  testthat::expect_true("sPSA.AGR.max.obs" %in% names(tmp))
  tmp <- byIndv4Intvl_ValueCalc(data=longi.dat, response = "sPSA.AGR", FUN="max", 
                                which.obs = T, which.values = "DAP")
  testthat::expect_equal(ncol(tmp), 4)  
  testthat::expect_equal(nrow(tmp), 20)  
  testthat::expect_true(all(c("sPSA.AGR.max.DAP",
                              "sPSA.AGR.max.obs") %in% names(tmp)))
  
  AGR.max.dat <- byIndv4Intvl_ValueCalc(data=longi.dat, response = "sPSA.AGR", 
                                        FUN="max", 
                                        start.time = 31, end.time = 35, 
                                        suffix.interval = "31to35", 
                                        which.values = "DAP", which.obs = TRUE)
  testthat::expect_equal(ncol(AGR.max.dat), 4)  
  testthat::expect_equal(nrow(AGR.max.dat), 20)  
  testthat::expect_true(all(c("Snapshot.ID.Tag", "sPSA.AGR.max.31to35",
                              "sPSA.AGR.max.obs.31to35", 
                              "sPSA.AGR.max.DAP.31to35") %in% names(AGR.max.dat)))
  testthat::expect_true(abs(AGR.max.dat[1,2] - 29.24427)  < 1e-03)  
  testthat::expect_true(AGR.max.dat[1,3] ==  5)  
  testthat::expect_true(AGR.max.dat[1,4] ==  35)  
  
  #Test byIndv_ValueCalc handling of which.levels, a deprecaed argument
  testthat::expect_error(byIndv_ValueCalc(data=longi.dat, response = "sPSA.AGR", 
                                          FUN="max", which.levels = "DAP"))

  })