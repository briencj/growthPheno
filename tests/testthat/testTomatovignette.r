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
  responses <- c("PSA", paste("PSA", c("AGR", "RGR"), sep = "."))
  responses.smooth <- c("sPSA", paste("sPSA", c("AGR", "RGR"), sep = "."))
  responses.logis <- paste(responses.smooth, "Logistic", sep = ".")
  resptitles <- c("PSA", "PSA AGR", "PSA RGR")
  resptitles.smooth <- c("sPSA", "sPSA AGR", "sPSA RGR")
  respunits <- c("(kpixels)", "(kpixels / day)", "( / day)")
  y.titles <- c("PSA (kpixels)", "PSA AGR (kpixels / day)", "PSA RGR ( / day)")
  names(y.titles) <- responses.smooth
  ypred.titles <- paste0("Predicted s", y.titles)
  names(ypred.titles) <- responses.smooth
  pred.type <- c("Predicted", "Backtransformed predicted")
  devn.titles <- c("PSA deviation (kpixels)", "PSA AGR deviation (kpixels / day)", 
                   "PSA RGR deviation ( / day)")
  names(devn.titles) <- responses.smooth
  x.title <- "DAP"
  # Specify time intervals of homogeneous growth dynamics
  DAP.endpts   <- c(18,22,27,33,39,43,51)
  nDAP.endpts <- length(DAP.endpts)
  DAP.starts <- DAP.endpts[-nDAP.endpts]
  DAP.stops   <- DAP.endpts[-1]
  DAP.mids <- (DAP.starts + DAP.stops)/2
  DAP.segs <- list(c(DAP.endpts[1]-1, 39), 
                   c(40, DAP.endpts[nDAP.endpts]))
  #Functions to label the plot facets
  labelAMF <- as_labeller(function(lev) paste(lev, "AMF"))
  labelZn <- as_labeller(function(lev) paste("Zn:", lev, "mg/kg"))
  vline.water <- list(geom_vline(xintercept=39, linetype="longdash", 
                                 alpha = 0.5, size=1))
  x.axis <- list(theme(axis.text.x = element_text(angle = 90),
                       panel.grid.minor.x = element_blank()))
  vline.DAP.endpts <- list(geom_vline(xintercept=DAP.starts, linetype="longdash", 
                                     alpha = 0.5, size=0.75))
  theme.profile <- list(vline.DAP.endpts,x.axis)
  
  # Step 1: Import, select and derive longitudinal data
  
  ## Load the pre-prepared data
  data(tomato.dat)
  
  ## Copy the data to preserve the original data.frame
  longi.dat <- tomato.dat
  
  ## Add continuous growth rates for raw data
  longi.dat <- byIndv4Times_GRsDiff(data = longi.dat, response = responses[1], 
                                    individuals = "Snapshot.ID.Tag", times = "DAP", 
                                    which.rates = c("AGR", "RGR"))
  testthat::expect_equal(nrow(longi.dat), 1120)
  testthat::expect_equal(ncol(longi.dat), 18)
  
  # Steps 2 & 3: Explore PSA and its AGR and RGR; investigate the smoothing of the PSA 
  
  ## Fit three-parameter logistic curves logistic curves to compare with spline curves
  
  ### Organize non-missing data into a grouped object
  logist.dat<- na.omit(longi.dat)
  logist.grp <- nlme::groupedData(PSA ~ cDAP | Snapshot.ID.Tag, 
                                  data = logist.dat)
  ### Fit the logistics and obtain fitted values
  logist.lis <- nlme::nlsList(SSlogis, logist.grp)
  logist.dat$sPSA <- fitted(logist.lis)
  
  ### Calculate the growth rates from the logistic fits
  logist.dat <- byIndv4Times_GRsDiff(logist.dat, responses = responses.smooth[1], 
                                     individuals = "Snapshot.ID.Tag", times = "DAP", 
                                     which.rates = c("AGR", "RGR"))
  logist.dat <- cbind(Type = factor("Logistic"), 
                      Tuning = factor("Logistic"), 
                      logist.dat)
  
  ## Probe the smoothing methods, DF and lambdas
  lambdas <- round(10^c(-0.5, 0, 0.5, 1), digits = 3)
  df = c(4:6,12)
  traits <- c("PSA","PSA.AGR","PSA.RGR")
  
  longi.nam <- names(longi.dat)
  longi.dat <- traitSmooth(data = longi.dat, 
                           response = "PSA", response.smoothed = "sPSA", 
                           individuals = "Snapshot.ID.Tag", times = "DAP", 
                           keep.columns = c("AMF","Zn"),
                           external.smooths = logist.dat, 
                           smoothing.segments = DAP.segs,
                           df = df, smoothing.methods = "log", 
                           facet.y.pf = "AMF", facet.y.med = "AMF", 
                           facet.y.chosen = "AMF",
                           labeller.chosen = labeller(Zn = labelZn, 
                                                      AMF = labelAMF),
                           colour.column.pf = "Zn", colour.column.chosen = "Zn",
                           ggplotFuncsProfile = theme.profile, 
                           ggplotFuncsChosen = c(theme.profile, vline.DAP.endpts))

  testthat::expect_equal(nrow(longi.dat), 1120)
  testthat::expect_equal(ncol(longi.dat), 21)
  testthat::expect_true(all(c(responses, responses.smooth) %in% names(longi.dat)))
  testthat::expect_true(all(longi.nam == names(longi.dat)[1:length(longi.nam)]))
  
  ## Plot the profile plots comparing log smoothing for NCSS with DF = 6 and PS with lambda = 1
  spar.schemes <- data.frame(Type = c("N", "P"), 
                             TunePar = c("df", "lam"),
                             TuneVal = c(6, 1),
                             Method = c("log", "log"))
  tune.fac <- c("Method","Type","Tuning")
  traitSmooth(data = longi.dat, 
              response = "PSA", response.smoothed = "sPSA", 
              individuals = "Snapshot.ID.Tag", times = "DAP", 
              keep.columns = c("AMF","Zn"),
              smoothing.schemes= spar.schemes,
              smoothing.segments = DAP.segs,
              df = df, smoothing.methods = "log", 
              chosen.smooth = NULL, 
              plots.by.pf = NULL, facet.x.pf = tune.fac, 
              facet.y.pf = "AMF",  
              facet.x.med = ".", facet.y.med = "AMF", 
              plots.group.med = tune.fac, 
              colour.column.pf = "AMF", 
              labeller = labeller(Zn = labelZn, 
                                  AMF = labelAMF),
              ggplotFuncsProfile = theme.profile)

  ### Omit responses for the outlier plant
  omit <- with(longi.dat, Zn==90 & AMF=="+" & Block ==4)
  responses.all <- names(longi.dat)[match("Weight.After", names(longi.dat)):length(longi.dat)]
  longi.dat[responses.all] <- lapply(longi.dat[responses.all], 
                                     function(kcol, omit) 
                                     {
                                       kcol[omit] <- NA
                                       return(kcol)
                                     }, omit = omit)
  

  ## Compute smooths of  WU for a range of smoothing parameters
  ## Here we take a two-step approach
  suppressWarnings(
    traitSmooth(data = longi.dat, 
                response = "WU", response.smoothed = "sWU", 
                individuals = "Snapshot.ID.Tag", times = "DAP", 
                keep.columns = c("AMF","Zn"), 
                trait.types = "response", 
                smoothing.segments = DAP.segs,
                df = df, smoothing.methods = "direct",
                chosen.smooth = NULL, 
                facet.y.pf = "AMF", facet.y.med = "AMF", 
                colour.column.pf = "Zn", colour.column.chosen = "Zn",
                ggplotFuncsProfile = theme.profile))
  
  #Compare
  spar.schemes <- data.frame(Type = c("N", "P"), 
                             TunePar = c("df", "lam"),
                             TuneVal = c(6, 0.316),
                             Method = c("dir", "dir"))
  suppressWarnings(
    traitSmooth(data = longi.dat, 
                response = "WU", response.smoothed = "sWU", 
                individuals = "Snapshot.ID.Tag", times = "DAP",
                get.rates = FALSE,
                smoothing.schemes= spar.schemes,
                smoothing.segments = DAP.segs,
                chosen.smooth = NULL, 
                plots.by.pf = NULL, facet.x.pf = tune.fac, 
                facet.y.pf = "AMF",  
                facet.x.med = ".", facet.y.med = "AMF", 
                plots.group.med = tune.fac, 
                colour.column.pf = "AMF", 
                labeller = labeller(Zn = labelZn, 
                                    AMF = labelAMF),
                ggplotFuncsProfile = theme.profile))
  
  #Extract chosen smooth from two compared
  longi.dat <- traitSmooth(data = longi.dat, 
                           response = "WU", response.smoothed = "sWU", 
                           individuals = "Snapshot.ID.Tag", times = "DAP", 
                           keep.columns = c("AMF","Zn"), 
                           trait.types = "response", 
                           smoothing.schemes= spar.schemes,
                           smoothing.segments = DAP.segs,
                           chosen.smooth = list(spline.type = "PS", 
                                                df = NULL, 
                                                lambda = 0.316,  #tried 1 first
                                                smoothing.method = "direct"), 
                           which.plots = NULL, 
                           facet.y.chosen = "AMF",
                           labeller.chosen = labeller(Zn = labelZn, 
                                                      AMF = labelAMF),
                           colour.column.chosen = "Zn",
                           ggplotFuncsChosen = c(theme.profile, vline.DAP.endpts))
  testthat::expect_equal(nrow(longi.dat), 1120)
  testthat::expect_equal(ncol(longi.dat), 22)
  testthat::expect_true(all(c(responses, responses.smooth, "sWU") %in% names(longi.dat)))
  testthat::expect_true(all(longi.nam == names(longi.dat)[1:length(longi.nam)]))
  
  # Step 5:  Per-cart trait formulation and extraction
  
  ## The time intervals of homogeneous growth dynamics (set in global.r)
  (suffices   <- paste(DAP.starts, DAP.stops, sep="to"))

  ## Commence cart.dat
  nidcols <- match("Weight.After", names(longi.dat))-1
  idcols.cart <- names(longi.dat)[1:nidcols]
  cart.dat <- longi.dat[longi.dat$DAP == DAP.endpts[1], idcols.cart]
  cart.dat <- cart.dat[, -match(c("DAP","cDAP"), names(cart.dat))]
  cart.dat <- cart.dat[do.call(order, cart.dat), ]
  testthat::expect_true(all(names(cart.dat) == c("Snapshot.ID.Tag", "Lane", "Position", 
                                                 "cPosn", "Block", "Cart", 
                                                 "AMF", "Zn", "Treatments")))

  ## Populate cart.dat at prescribed time-points
  for (day in DAP.endpts)
    cart.dat <- cbind(cart.dat,
                      getTimesSubset(responses.smooth, data = longi.dat,
                                     times = "DAP", which.times = day,
                                     suffix = as.character(day)))

  ## Populate cart.dat and cart.allt.intervals at prescribed intervals
  for (r in responses.smooth[1])
  { 
    for (k in 1:length(suffices))
    { 
      cart.dat <- merge(cart.dat, 
                        byIndv4Intvl_GRsDiff(data = longi.dat, responses = r, 
                                             individuals = "Snapshot.ID.Tag", 
                                             times = "DAP", 
                                             which.rates = c("AGR","RGR"), 
                                             start.time = DAP.starts[k], 
                                             end.time = DAP.stops[k], 
                                             suffix.interval = suffices[k]),
                        by = "Snapshot.ID.Tag")
      
    }
  }

  ## Finalise
  cart.dat <- with(cart.dat, cart.dat[order(Snapshot.ID.Tag), ])
  testthat::expect_equal(nrow(cart.dat), 32)
  testthat::expect_equal(ncol(cart.dat), 42)

})
