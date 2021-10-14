
cat("#### Test probeSmoothing with Judith Atieno 0278\n")
test_that("chickpea_growthPheno", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(growthPheno)
  library(ggplot2)
  
  data(dat1)
  
  #'## Values of df for which to obtain plots
  df <- c(4,7)
  
  
  #'## Obtain plots
  testthat::expect_warning(t <- probeSmoothing(dat1, df=df, xname="TimeAfterPlanting", 
                                               response="ShootArea1000",
                                               which.plots=c("bothseparately"), 
                                               ggplotFuncs=list(scale_x_continuous(breaks=
                                                                                     seq(30, 50, by=2)))))
  testthat::expect_equal(nrow(t), 8480)
  testthat::expect_equal(ncol(t), 15)
  
  testthat::expect_silent(t <- probeSmoothing(dat1, df=df, xname="TimeAfterPlanting", 
                                               response="ShootArea1000", 
                                               trait.types=c("AGR")))
  testthat::expect_equal(nrow(t), 8480)
  testthat::expect_equal(ncol(t), 12)
  
  testthat::expect_silent(t <- probeSmoothing(dat1, df=df, xname="TimeAfterPlanting", 
                                              response="ShootArea1000", 
                                              trait.types=c("response"), 
                                              which.plots="dfcompare", get.rates=FALSE))
  testthat::expect_equal(nrow(t), 8480)
  testthat::expect_true(all(c("Snapshot.ID.Tag", "Days", "TimeAfterPlanting", "ShootArea1000", 
                              "Treatment.1", "Smarthouse",  "ShootArea1000", 
                              "ShootArea1000.smooth.Direct.4", "ShootArea1000.smooth.Direct.7") 
                            %in% names(t)))
  
  #'## Include deviations boxplots
  testthat::expect_error(t <- probeSmoothing(dat1, df=df, xname="TimeAfterPlanting", 
                                             response="ShootArea1000",
                                             x.title = "DAP", 
                                             trait.types=c("response"), get.rates=FALSE,
                                             which.plots="dfcompare", 
                                             deviations.boxplots = "absolute")) 
  #deviations.boxplots deprecated
  testthat::expect_silent(t <- probeSmoothing(dat1, df=df, xname="TimeAfterPlanting", 
                                              response="ShootArea1000",
                                              x.title = "DAP", 
                                              trait.types=c("response"), get.rates=FALSE,
                                              which.plots="dfcompare", 
                                              deviations.plots = "absolute"))
  testthat::expect_equal(nrow(t), 8480)
  testthat::expect_true(all(c("Snapshot.ID.Tag", "Days", "TimeAfterPlanting", "ShootArea1000", 
                              "Treatment.1", "Smarthouse", "ShootArea1000", 
                              "ShootArea1000.smooth.Direct.4", "ShootArea1000.smooth.Direct.7") 
                            %in% names(t)))
  
  testthat::expect_warning(t <- probeSmoothing(dat1, df=df, xname="TimeAfterPlanting", 
                                               response="ShootArea1000",
                                               x.title = "DAP", 
                                               trait.types=c("response", "AGR"), 
                                               which.plots="methodscompare", 
                                               deviations.plots = c("absolute", "relative")))
  testthat::expect_equal(nrow(t), 8480)
  testthat::expect_true(all(c("Snapshot.ID.Tag", "Days", "TimeAfterPlanting", "ShootArea1000", 
                              "Treatment.1", "Smarthouse", "ShootArea1000", 
                              "Days.diffs", "ShootArea1000.AGR", 
                              "ShootArea1000.smooth.Direct.4", "ShootArea1000.smooth.AGR.Direct.4", 
                              "ShootArea1000.smooth.Direct.7", 
                              "ShootArea1000.smooth.AGR.Direct.7") 
                            %in% names(t)))
  testthat::expect_lt(max(abs(t["ShootArea1000"] - t["ShootArea1000.smooth.Direct.4"]), 
                          na.rm = TRUE) - 5.563407, 1e-07)
  testthat::expect_lt(max(abs(t["ShootArea1000"] - t["ShootArea1000.smooth.Direct.4"])/
                            t["ShootArea1000.smooth.Direct.4"]) - 370.7951, 1e-04)
  
})

cat("#### Test probeSmoothing with small example\n")
test_that("exampleData_growthPheno", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(growthPheno)
  
  data(exampleData)
  
  vline <- list(ggplot2::geom_vline(xintercept=20, linetype="longdash", size=1),
                ggplot2::scale_x_continuous(breaks=seq(12, 36, by=2)))
  
  plotDeviationsBoxes(longi.dat, observed = "Area", smoothed = "Area.smooth",
                    x.factor="Days", facet.x = ".", facet.y= ".", df =5)
  
  testthat::expect_silent(tmp <- probeSmoothing(data = longi.dat, response = "Area", 
                                                df = c(4,7), x="xDays+24.16666667", 
                                                ggplotFuncs=vline))
  testthat::expect_equal(nrow(tmp), 280)
  testthat::expect_equal(ncol(tmp), 15)
  
  testthat::expect_warning(tmp <- probeSmoothing(data = longi.dat, response = "Area", 
                                                 df = c(4,7), x="xDays+24.16666667", 
                                                 which.plots = "methodsc", 
                                                 facet.y = ".",
                                                 get.rates = FALSE,
                                                 ggplotFuncs=vline))
  testthat::expect_equal(nrow(tmp), 280)
  testthat::expect_equal(ncol(tmp), 7)
  testthat::expect_warning(tmp <- probeSmoothing(data = longi.dat, response = "Area", 
                                                 df = c(4,7), x="xDays+24.16666667", 
                                                 which.plots = "dfcompare", 
                                                 facet.y = ".",
                                                 get.rates = FALSE,
                                                 ggplotFuncs=vline))
  testthat::expect_equal(nrow(tmp), 280)
  testthat::expect_equal(ncol(tmp), 7)
  
  testthat::expect_warning(tmp <- probeSmoothing(data = longi.dat, response = "Area", 
                                                 df = c(4:7), x="xDays+24.16666667", 
                                                 which.plots = "dfcompare", 
                                                 facet.y = ".",
                                                 alpha = 0.6, get.rates = FALSE,
                                                 ggplotFuncs=vline))
  testthat::expect_equal(nrow(tmp), 280)
  testthat::expect_equal(ncol(tmp), 9)
  
  testthat::expect_silent(tmp <- probeSmoothing(data = longi.dat, response = "Area", 
                                                df = c(4:7), x="xDays+24.16666667", 
                                                which.plots = "dfcompare", facet.y = ".",
                                                alpha = 0.5, trait.types = "AGR",
                                                ggplotFuncs=vline))
  testthat::expect_equal(nrow(tmp), 280)
  testthat::expect_equal(ncol(tmp), 15)
  
  testthat::expect_warning(tmp <- probeSmoothing(data = longi.dat, response = "Area", 
                                                 df = c(4:7), x="xDays+24.16666667", 
                                                 facet.y = ".", 
                                                 deviations.plots = "compare",
                                                 propn.types = c(0.025,0.2, 0.25),
                                                 ggplotFuncs=vline))
  testthat::expect_equal(nrow(tmp), 280)
  testthat::expect_equal(ncol(tmp), 20)
  
  testthat::expect_warning(tmp <- probeSmoothing(data = longi.dat, response = "Area", 
                                                 df = c(4:7), x="xDays+24.16666667", 
                                                 facet.x = ".", facet.y = ".",
                                                 which.plots = "none",
                                                 deviations.plots = "compare", 
                                                 propn.types = NULL, 
                                                 ggplotFuncs=vline))
  testthat::expect_equal(nrow(tmp), 280)
  testthat::expect_equal(ncol(tmp), 19)
  
  testthat::expect_silent(traits <- probeSmoothing(data = longi.dat, response = "Area", 
                                                   df = c(4:7), x="xDays+24.16666667", 
                                                   facet.x = ".", facet.y = ".",
                                                   which.plots = "none",
                                                   deviations.plots = "none", 
                                                   propn.types = NULL))
  testthat::expect_silent(med <- plotMedianDeviations(data = traits, 
                                                      response = "Area", 
                                                      response.smoothed = "Area.smooth", 
                                                      x="xDays+24.16666667", xname = "xDays", 
                                                      df = c(4,7), x.title = "DAP", 
                                                      facet.x = ".", facet.y = ".",
                                                      trait.types = "response", 
                                                      propn.types = 0.05,
                                                      ggplotFuncsMedDevn = vline))
  testthat::expect_equal(length(med), 2)
  testthat::expect_true(all(names(med) == c("plots", "med.devn.dat")))
  testthat::expect_equal(nrow(med$med.devn.dat), 28)
  testthat::expect_equal(ncol(med$med.devn.dat), 6)
  testthat::expect_equal(length(med$plots), 1)
  testthat::expect_warning(traits <- probeSmoothing(data = longi.dat, response = "Area", 
                                                    df = c(4:7), x="xDays+24.16666667", 
                                                    facet.x = ".", facet.y = ".",
                                                    which.plots = "none",
                                                    deviations.plots = c("compare", 
                                                                         "absolute"),
                                                    propn.types = NULL))
  testthat::expect_equal(nrow(traits), 280)
  testthat::expect_equal(ncol(traits), 19)
  testthat::expect_warning(traits <- probeSmoothing(data = longi.dat, response = "Area", 
                                                   df = c(4:7), x="xDays+24.16666667", 
                                                   facet.x = ".", facet.y = ".",
                                                   which.plots = "methodsc",
                                                   deviations.plots = c("compare", 
                                                                        "absolute"),
                                                   propn.types = NULL))
  testthat::expect_equal(nrow(traits), 280)
  testthat::expect_equal(ncol(traits), 19)
  
})


cat("#### Test probeSmoothing with tomato example\n")
test_that("tomato_growthPheno", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(dae)
  library(growthPheno)
  
  data(tomato.dat)
  tomato.dat <- within(tomato.dat, xDAP <- as.numfac(DAP))
  
  df.vec        <- c(4:6,12)
  labelMyc <- as_labeller(function(lev) paste(lev, "AMF"))
  labelZn <- as_labeller(function(lev) paste("Zn:", lev, "ppm"))
  
  #'## Gives error that the Length of propn.types is not the same as the number of traits
  testthat::expect_warning(testthat::expect_error(
    tmp <- probeSmoothing(data = tomato.dat, 
                          response = "Area", 
                          times.factor = "DAP", xname = "xDAP", 
                          smoothing.methods = c("dir", "log"), 
                          facet.x = "Zn", facet.y = "AMF",
                          df = c(4,7), x="xDAP", get.rates = FALSE,
                          which.plots = "methodsc", 
                          deviations.plots = "compare",
                          labeller = labeller(Zn = labelZn, 
                                              AMF = labelMyc)),
    regexp = "Length of propn.types is not the same as the number of trait.types"), 
    regexp = "get.rates is FALSE and so trait.types changed to response")
  testthat::expect_warning(tom <- probeSmoothing(data = tomato.dat, 
                                                 response = "Area", 
                                                 times.factor = "DAP", xname = "xDAP", 
                                                 smoothing.methods = c("dir", "log"), 
                                                 facet.x = "Zn", facet.y = "AMF",
                                                 df = c(4,7), x="xDAP", get.rates = FALSE,
                                                 propn.types = 0.1,
                                                 which.plots = "methodsc", 
                                                 deviations.plots = "compare",
                                                 labeller = labeller(Zn = labelZn, 
                                                                     AMF = labelMyc)), 
                           regexp = "get.rates is FALSE and so trait.types changed to response")
  testthat::expect_equal(nrow(tom), 1120)
  testthat::expect_equal(ncol(tom), 10)
  
  testthat::expect_silent(med <- plotMedianDeviations(data = tom, 
                                                      response = "Area", 
                                                      response.smoothed = "Area.smooth", 
                                                      xname = "xDAP", 
                                                      smoothing.methods = c("dir", "log"), 
                                                      df = c(4,7), x.title = "DAP", 
                                                      y.titles = "PSA deviation (kpixels)",
                                                      facet.x = "Zn", facet.y = "AMF",
                                                      trait.types = "response", 
                                                      propn.types = 0.1,
                                                      labeller = labeller(Zn = labelZn, 
                                                                          AMF = labelMyc)))
  
  testthat::expect_equal(length(med$plots), 1)
  testthat::expect_equal(nrow(med$med.devn.dat), 1120)
  testthat::expect_equal(ncol(med$med.devn.dat), 8)
  
  #Multiple df, single methods
  testthat::expect_warning(tom <- probeSmoothing(data = tomato.dat, response = "Area", 
                                                 times.factor = "DAP", xname="xDAP", df=5:6,
                                                 smoothing.methods = c("logarithmic"),
                                                 facet.x = ".", facet.y = ".",
                                                 which.plots = "none",
                                                 propn.types = c(0.02, 0.2, 0.5),
                                                 deviations.plots = "compare.medians"))
  testthat::expect_equal(nrow(tom), 1120)
  testthat::expect_equal(ncol(tom), 13)
  
  #'Single `df`, multiple methods and trait.types
  testthat::expect_warning(tom <- probeSmoothing(data = tomato.dat, response = "Area", 
                                                 times.factor = "DAP", xname="xDAP", df=5,
                                                 smoothing.methods = c("direct","logarithmic"),
                                                 facet.x = ".", facet.y = ".",
                                                 which.plots = "none",
                                                 deviations.plots = "none"))
  testthat::expect_equal(nrow(tom), 1120)
  testthat::expect_equal(ncol(tom), 13)
  testthat::expect_warning(med <- plotMedianDeviations(data = tom, response = "Area", 
                                                       response.smoothed = "Area.smooth", 
                                                       xname = "xDAP", 
                                                       smoothing.methods = c("dir", "log"), 
                                                       df = 5, x.title = "DAP", 
                                                       facet.x = ".", facet.y = ".",
                                                       propn.types = c(0.02, 0.2, 0.5)))
  
  testthat::expect_equal(length(med), 2)
  testthat::expect_true(all(names(med) == c("plots", "med.devn.dat")))
  testthat::expect_equal(nrow(med$med.devn.dat), 70)
  testthat::expect_equal(ncol(med$med.devn.dat), 8)
  testthat::expect_equal(length(med$plots), 3)
  testthat::expect_true(all(names(med$plots) == c("Area","Area.AGR","Area.RGR")))
  testthat::expect_warning(print(med$plots$Area.AGR))
  
  #'Single `df`, single method
  testthat::expect_warning(tom <- probeSmoothing(data = tomato.dat, response = "Area", 
                                                 times.factor = "DAP", xname="xDAP", df=5,
                                                 smoothing.methods = c("direct"),
                                                 facet.x = ".", facet.y = ".",
                                                 which.plots = "none",
                                                 propn.types = c(0.02, 0.2, 0.5),
                                                 deviations.plots = "compare.medians"))
  testthat::expect_equal(nrow(tom), 1120)
  testthat::expect_equal(ncol(tom), 10)

})

