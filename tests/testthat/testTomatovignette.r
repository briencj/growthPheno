#devtools::test("growthPheno")

cat("#### Test using Tomato vignette\n")
test_that("Tomatovignette_growthPheno", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(dae)
  library(ggplot2)
  library(growthPheno)
  library(nlme)

  ## Set up
  # The responses
  responses <- c("PSA", paste("PSA", c("AGR", "RGR"), sep = "."))
  responses.smooth <- paste0("s", responses)
  
  #Specify time intervals of homogeneous growth dynamics
  DAP.endpts   <- c(18,22,27,33,39,43,51)
  nDAP.endpts <- length(DAP.endpts)
  DAP.starts <- DAP.endpts[-nDAP.endpts]
  DAP.stops   <- DAP.endpts[-1]
  DAP.segs <- list(c(DAP.endpts[1]-1, 39), 
                   c(40, DAP.endpts[nDAP.endpts]))
  tune.fac <- c("Method","Type","Tuning")
  #Functions to label the plot facets
  labelAMF <- as_labeller(function(lev) paste(lev, "AMF"))
  labelZn <- as_labeller(function(lev) paste("Zn:", lev, "mg/kg"))
  vline.water <- list(geom_vline(xintercept=39, linetype="longdash", 
                                 alpha = 0.3, size=1))
  vline.DAP.endpts <- list(geom_vline(xintercept=DAP.starts, linetype="longdash", 
                                      colour = "blue", alpha = 0.5, size=0.75))
  
  # Step I: Import the longitudinal data
  
  ## Load the pre-prepared data
  data(tomato.dat)
  
  ## Copy the data to preserve the original data.frame
  longi.dat <- tomato.dat
  
  # Step II: Investigate the smoothing of the PSA and obtain growth rates
  
  ## Fit three-parameter logistic curves logistic curves to compare with spline curves
  
  ### Organize non-missing data into a grouped object
  logist.dat<- na.omit(longi.dat)
  logist.grp <- nlme::groupedData(PSA ~ cDAP | Snapshot.ID.Tag, 
                                  data = logist.dat)
  ### Fit the logistics and obtain fitted values
  logist.lis <- nlme::nlsList(SSlogis, logist.grp)
  logist.dat$sPSA <- fitted(logist.lis)
  logist.dat <- cbind(Type = factor("Logistic"), 
                      Tuning = factor("Logistic"), 
                      logist.dat)
  
  ## Probe the smoothing methods, DF and lambdas
  longi.nam <- names(longi.dat)
  suppressWarnings(
    smth.dat <- traitSmooth(data = tomato.dat, 
                            response = "PSA", response.smoothed = "sPSA", 
                            individuals = "Snapshot.ID.Tag", times = "DAP", 
                            keep.columns = c("AMF","Zn"), 
                            get.rates = TRUE, trait.types = c("response", "AGR", "RGR"), 
                            smoothing.args = args4smoothing(df = c(4:6,12),  
                                                            smoothing.segments = DAP.segs, 
                                                            external.smooths = logist.dat),
                            profile.plot.args = 
                              args4profile_plot(facet.y = "AMF", 
                                                colour.column = "Zn", 
                                                facet.labeller = labeller(AMF = labelAMF)), 
                            meddevn.plot.args = 
                              args4meddevn_plot(facet.y = "AMF", 
                                                facet.labeller = labeller(AMF = labelAMF)), 
                            chosen.plot.args = 
                              args4chosen_plot(facet.y = "AMF",
                                               facet.labeller = labeller(AMF = labelAMF),
                                               colour.column = "Zn", 
                                               ggplotFuncs = vline.DAP.endpts), 
                            mergedata = tomato.dat))

  testthat::expect_equal(nrow(smth.dat), 1120)
  testthat::expect_equal(ncol(smth.dat), 21)
  testthat::expect_true(all(c(responses, responses.smooth) %in% names(smth.dat)))
  testthat::expect_true(all(longi.nam %in% names(smth.dat)))
  
  ## Plot the profile plots comparing log smoothing for NCSS with DF = 6 and PS with lambda = 1
  smth.dat <- traitSmooth(data = longi.dat, 
                          response = "PSA", response.smoothed = "sPSA", 
                          individuals = "Snapshot.ID.Tag", times = "DAP", 
                          keep.columns = c("AMF","Zn"),
                          smoothing.args = 
                            args4smoothing(smoothing.methods = c("log", "log"), 
                                           spline.types = c("N", "P"), 
                                           df = c(6, NA), lambdas = c(NA, 1),
                                           combinations = "parallel",
                                           smoothing.segments = DAP.segs),
                          chosen.smooth.args = NULL, 
                          profile.plot.args = 
                            args4profile_plot(plots.by = NULL, 
                                              facet.x = tune.fac, facet.y = "AMF", 
                                              facet.labeller = labeller(AMF = labelAMF), 
                                              colour.column = "AMF"), 
                          meddevn.plot.args = 
                            args4meddevn_plot(plots.by = NULL, plots.group = tune.fac, 
                                              facet.x = ".", facet.y = "AMF", 
                                              facet.labeller = labeller(AMF = labelAMF)))
  testthat::expect_equal(nrow(smth.dat), 2240)
  testthat::expect_equal(ncol(smth.dat), 16)
  
  ## Extract the chosen smooth, adding it to longi.dat
  
  longi.dat <- traitSmooth(data = smth.dat, 
                           response = "PSA", response.smoothed = "sPSA", 
                           individuals = "Snapshot.ID.Tag", times = "DAP", 
                           keep.columns = c("AMF","Zn"), 
                           smoothing.args = NULL, which.plots = "none", 
                           chosen.smooth.args = 
                             args4chosen_smooth(smoothing.methods = "log", 
                                                spline.types = "PS", lambdas = 1), 
                           chosen.plot.args = 
                             args4chosen_plot(facet.y = "AMF",
                                              facet.labeller = labeller(AMF = labelAMF),
                                              colour.column = "Zn", 
                                              ggplotFuncs = vline.DAP.endpts), 
                           mergedata = tomato.dat)
  
  testthat::expect_equal(nrow(longi.dat), 1120)
  testthat::expect_equal(ncol(longi.dat), 21)
  testthat::expect_true(all(c(responses, responses.smooth) %in% names(longi.dat)))
  testthat::expect_true(all(longi.nam %in% names(longi.dat)))
  
  
  # Step III: Investigate the smoothing of the WU
  
  ## Explore the smooths of WU for a range of smoothing parameters
  suppressWarnings(
    smth.dat <- traitSmooth(data = longi.dat, 
                            response = "WU", response.smoothed = "sWU", 
                            individuals = "Snapshot.ID.Tag", times = "DAP", 
                            keep.columns = c("AMF","Zn"), 
                            trait.types = "response", 
                            smoothing.args = 
                              args4smoothing(smoothing.methods = "direct", 
                                             smoothing.segments = DAP.segs),
                            chosen.smooth.args = NULL, 
                            profile.plot.args = 
                              args4profile_plot(plots.by = NULL,
                                                facet.y = "AMF", 
                                                colour.column = "Zn", 
                                                facet.labeller = labeller(AMF = labelAMF)),
                            meddevn.plot.args = 
                              args4meddevn_plot(plots.by = NULL,
                                                facet.y = "AMF",  
                                                facet.labeller = labeller(AMF = labelAMF))))
  testthat::expect_equal(nrow(smth.dat), 7840)
  testthat::expect_equal(ncol(smth.dat), 11)
  
  #Compare two
  suppressWarnings(
    traitSmooth(data = smth.dat, 
                response = "WU", response.smoothed = "sWU", 
                individuals = "Snapshot.ID.Tag", times = "DAP",
                trait.types = "response", 
                smoothing.args = args4smoothing(smoothing.methods = c("dir", "dir"), 
                                                spline.types = c("N", "P"), 
                                                df = c(6, NA), lambdas = c(NA, 0.316),
                                                smoothing.segments = DAP.segs, 
                                                combinations = "parallel"),
                chosen.smooth.args = NULL, 
                profile.plot.args = 
                  args4profile_plot(plots.by = NULL, 
                                    facet.x = tune.fac, facet.y = "AMF",  
                                    colour.column = "AMF", 
                                    facet.labeller = labeller(AMF = labelAMF)), 
                meddevn.plot.args = 
                  args4meddevn_plot(plots.by = NULL, plots.group = tune.fac, 
                                    facet.x = ".", facet.y = "AMF", 
                                    facet.labeller = labeller(AMF = labelAMF))))  
  
  #Extract chosen smooth from two compared
  longi.dat <- traitSmooth(data = smth.dat, 
                           response = "WU", response.smoothed = "sWU", 
                           individuals = "Snapshot.ID.Tag", times = "DAP", 
                           trait.types = "response",  
                           smoothing.args = NULL, which.plots = NULL, 
                           chosen.smooth.args = 
                             args4chosen_smooth(smoothing.method = "direct",
                                                spline.type = "PS", 
                                                lambdas = 0.316),  #tried 1 first
                           chosen.plot.args = 
                             args4chosen_plot(facet.y = "AMF", 
                                              facet.labeller = labeller(AMF = labelAMF),
                                              colour.column = "Zn",
                                              ggplotFuncs = vline.DAP.endpts),
                           mergedata = longi.dat)
  testthat::expect_equal(nrow(longi.dat), 1120)
  testthat::expect_equal(ncol(longi.dat), 22)
  testthat::expect_true(all(c(responses, responses.smooth, "sWU") %in% names(longi.dat)))
  testthat::expect_true(all(longi.nam %in% names(longi.dat)))
  
  # Step IV: Identify potential outliers and remove if justified

  ## Omit responses for the outlier plant
  omit <- with(longi.dat, Zn==90 & AMF=="+" & Block ==4)
  responses.all <- names(longi.dat)[match("Weight.After", names(longi.dat)):length(longi.dat)]
  longi.dat[responses.all] <- lapply(longi.dat[responses.all], 
                                     function(kcol, omit) 
                                     {
                                       kcol[omit] <- NA
                                       return(kcol)
                                     }, omit = omit)
  
  # Step V:  Extract single-valued traits for each individual
  
  ## Produce indv.dat that contains the extracted traits based on the time intervals of homogeneous growth dynamics
  indv.cols <- c("Snapshot.ID.Tag", "Lane", "Position", "Block", "Cart", "AMF", "Zn")
  indv.dat <- subset(longi.dat, subset = DAP == DAP.endpts[1], 
                     select = indv.cols)
  indv.dat <- traitExtractFeatures(data = longi.dat, 
                                   starts.intvl = DAP.starts, stops.intvl = DAP.stops, 
                                   responses4singletimes = "sPSA", 
                                   responses4intvl.rates = "sPSA", growth.rates = c("AGR", "RGR"), 
                                   water.use4intvl.traits = "sWU", responses4water = "sPSA", 
                                   responses4overall.total = "sWU",
                                   responses4overall.max = "sPSA.AGR",
                                   mergedata = indv.dat)
  
  ## Finalise
  indv.dat <- with(indv.dat, indv.dat[order(Snapshot.ID.Tag), ])
  summary(indv.dat)
  head(indv.dat)
  testthat::expect_equal(nrow(indv.dat), 32)
  testthat::expect_equal(ncol(indv.dat), 47)

})

