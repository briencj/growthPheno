#devtools::test("growthPheno")

cat("#### Test using Rice vignette with no plots\n")
test_that("Rice2015_growthPheno", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(dae)
  library(ggplot2)
  library(growthPheno)
  
  ## A dummy example to illustrate the use of growthPheno
  #'# Step 1: Import the raw data
  data(RiceRaw.dat)
  raw.dat <- RiceRaw.dat[1:280, ]
  raw.dat$Smarthouse <- 1
  
  #'# Step 2: Select imaging variables and add covariates and factors (produces longi.dat)
  longi.dat <- prepImageData(data=raw.dat, smarthouse.lev=1, 
                             idcolumns = c("Genotype.ID", "Treatment.1", "Treatment.2"))
  
  longi.dat <- designFactors(longi.dat, insertName = "Reps",
                             nzones = 1, nlanesperzone = 1, nmainunitsperlane = 10, 
                             designfactorMethod="StandardOrder")
  testthat::expect_equal(nrow(longi.dat), 280)
  testthat::expect_equal(ncol(longi.dat), 44)
  
  
  #'## Particular edits to longi.dat - add Days after treatment (DAT)
  longi.dat <- within(longi.dat, 
                      DAT <- xDAP - 29)
  
  #'# Step 3: Form derived traits that result in a value for each observation
  #'### Set responses
  responses.image <- c("PSA")
  responses.smooth <- paste0("s", responses.image)
  
  # Form growth rates for each observation of a subset of responses by differencing
  longi.dat <- byIndv4Times_GRsDiff(longi.dat, responses = responses.image, 
                                    times = "DAP", 
                                    which.rates = c("AGR","RGR"))
  
  # Form PSA.WUI 
  longi.dat <- within(longi.dat, 
                      PSA.WUI <- WUI(PSA.AGR*DAP.diffs, WU))
  
  # Add cumulative responses 
  longi.dat <- within(longi.dat, 
                      { 
                        WU.cum <- unlist(by(WU, Snapshot.ID.Tag, 
                                            cumulate, exclude.1st=TRUE))
                        WUI.cum <- PSA / WU.cum 
                      })
  # Check longi.dat
  head(longi.dat)

  #'# Step 4: Fit splines to smooth the longitudinal trends in the primary traits and calculate their growth rates
  #'
  #'## Smooth responses and form growth rates by differences
  #+
  for (response in c(responses.image, "WU"))
    longi.dat <- byIndv4Times_SplinesGRs(data = longi.dat, response = response, 
                                         response.smoothed = paste0("s", response),
                                         individuals = "Snapshot.ID.Tag", times="DAP",  
                                         df = 4)
  testthat::expect_equal(nrow(longi.dat), 280)
  testthat::expect_equal(ncol(longi.dat), 57)
  
  #'## Finalize longi.dat
  longi.dat <- with(longi.dat, longi.dat[order(Snapshot.ID.Tag, xDAP), ])
  
  #'# Step 5: Do exploratory plots on unsmoothed and smoothed longitudinal data

  #'# Step 6: Form single-value plant responses in Snapshot.ID.Tag order.
  #'

  #'### Set up intervals
  #+
  DAP.endpts <- c(31,35,38,42)
  DAP.starts <- c(31,35,31,38)
  DAP.stops   <- c(35,38,38,42)
  DAP.mids <- (DAP.starts + DAP.stops)/2
  suffices <- paste(DAP.starts, DAP.stops, sep = "to")
  
  #'## 6a) Set up a data frame with factors only
  #+
  cart.dat <- longi.dat[longi.dat$DAP == DAP.endpts[1], 
                        c("Smarthouse","Lane","Position","Snapshot.ID.Tag",
                          "cPosn","cMainPosn",
                          "Zone","cZone","SHZone","ZLane","ZMainunit", "Subunit",
                          "Genotype.ID","Treatment.1")]
  cart.dat <- cart.dat[do.call(order, cart.dat), ]
  
  #'## 6b) Get responses based on first and last date.
  #'
  #'### Observation for first and last date
  cart.dat <- cbind(cart.dat, getTimesSubset(data = longi.dat, responses = responses.image, 
                                             times = "DAP", which.times = DAP.endpts[1], 
                                             suffix = "first"))
  cart.dat <- cbind(cart.dat, getTimesSubset(data = longi.dat, responses = responses.image, 
                                             times = "DAP", 
                                             which.times = DAP.endpts[length(DAP.endpts)], 
                                             suffix = "last"))
  cart.dat <- cbind(cart.dat, getTimesSubset(data = longi.dat, responses = "WUI.cum", 
                                             times = "DAP", 
                                             which.times = DAP.endpts[length(DAP.endpts)], 
                                             suffix = "last"))
  responses.smooth <- paste0("s", responses.image)
  cart.dat <- cbind(cart.dat, getTimesSubset(data = longi.dat, responses = responses.smooth, 
                                             times = "DAP", which.times = DAP.endpts[1], 
                                             suffix = "first"))
  cart.dat <- cbind(cart.dat, getTimesSubset(data = longi.dat, responses = responses.smooth, 
                                             times = "DAP", 
                                             which.times = DAP.endpts[length(DAP.endpts)], 
                                             suffix = "last"))
  testthat::expect_equal(nrow(cart.dat), 20)
  testthat::expect_equal(ncol(cart.dat), 19)
  testthat::expect_true(all(c( "PSA.first", "PSA.last", "WUI.cum.last", 
                               "sPSA.first", "sPSA.last") %in% names(cart.dat)))
  
  # Growth rates over whole period.
  (tottime <- DAP.endpts[length(DAP.endpts)] - DAP.endpts[1]) #= 11
  cart.dat <- within(cart.dat, 
                     { 
                       PSA.AGR.full <- (PSA.last - PSA.first)/tottime
                       PSA.RGR.full <- log(PSA.last / PSA.first)/tottime
                     })
  
  # Calculate water index over whole period
  cart.dat <- merge(cart.dat, 
                    byIndv4Intvl_WaterUse(data = longi.dat, 
                                          water.use = "WU", response = "PSA", 
                                          trait.types = c("WUI","WUR", "WU"), 
                                          times = "DAP", 
                                          start.time = DAP.endpts[1], 
                                          end.time = DAP.endpts[length(DAP.endpts)]),
                    by = c("Snapshot.ID.Tag"))
  

  #'## 6c) Add growth rates and water indices for intervals
  
  #'### Rates for specific intervals from the smoothed data by differencing
  #+
  for (r in responses.smooth)
  { 
    for (k in 1:length(suffices))
    { 
      cart.dat <- merge(cart.dat, 
                        byIndv4Intvl_GRsDiff(data = longi.dat, responses = r, 
                                             times = "DAP", 
                                             which.rates = c("AGR","RGR"), 
                                             start.time = DAP.starts[k], 
                                             end.time = DAP.stops[k], 
                                             suffix.interval = suffices[k]),
                        by = "Snapshot.ID.Tag")
    }
  }
  testthat::expect_true(all(c(paste("sPSA.AGR", suffices, sep = "."), 
                              paste("sPSA.RGR", suffices, sep = ".")) %in% names(cart.dat)))
  
  #'### Water indices for specific intervals from the unsmoothed and smoothed data
  #+
  for (k in 1:length(suffices))
  { 
    cart.dat <- merge(cart.dat, 
                      byIndv4Intvl_WaterUse(data = longi.dat, 
                                            water.use = "WU", responses = "PSA", 
                                            times = "DAP", 
                                            trait.types = c("WU","WUR","WUI"), 
                                            start.time = DAP.starts[k], 
                                            end.time = DAP.stops[k], 
                                            suffix.interval = suffices[k]),
                      by = "Snapshot.ID.Tag")
  }
  testthat::expect_true(all(c(paste("WU", suffices, sep = "."),
                              paste("PSA.WUI", suffices, sep = "."), 
                              paste("WU", suffices, sep = ".")) %in% names(cart.dat)))
  
  cart.dat <- with(cart.dat, cart.dat[order(Snapshot.ID.Tag), ])
  testthat::expect_equal(nrow(cart.dat), 20)
  testthat::expect_equal(ncol(cart.dat), 44)
  
  #'# Step 7: Form continuous and interval SIITs
  #'
  #'## 7a) Calculate continuous
  #+
  cols.retained <-  c("Snapshot.ID.Tag","Smarthouse","Lane","Position",
                      "DAP","Snapshot.Time.Stamp", "Hour", "xDAP",
                      "Zone","cZone","SHZone","ZLane","ZMainunit",
                      "cMainPosn", "Genotype.ID")
  responses.GR <- c("sPSA.AGR","sPSA.AGR","sPSA.RGR")
  suffices.results <- c("diff", "SIIT", "SIIT")
  responses.SIIT <- unlist(Map(paste, responses.GR, suffices.results,sep="."))
  
  longi.SIIT.dat <- 
    twoLevelOpcreate(data = longi.dat, responses = responses.GR, 
                     suffices.treatment=c("C","S"),
                     operations = c("-", "/", "/"), suffices.results = suffices.results, 
                     columns.retained = cols.retained, 
                     by = c("Smarthouse","Zone","ZMainunit","DAP"))
  longi.SIIT.dat <- with(longi.SIIT.dat, 
                         longi.SIIT.dat[order(Smarthouse,Zone,ZMainunit,DAP),])
  testthat::expect_equal(nrow(longi.SIIT.dat), 140)
  testthat::expect_equal(ncol(longi.SIIT.dat), 22)
  
  #' ### Plot SIIT profiles 
  #' 
  #+ "03-SIITProfiles"
  k <- 2
  nresp <- length(responses.SIIT)
  limits <- with(longi.SIIT.dat, list(c(min(sPSA.AGR.diff, na.rm=TRUE),
                                        max(sPSA.AGR.diff, na.rm=TRUE)),
                                      c(0,3),
                                      c(0,1.5)))
  #Plots

  #'## 7b) Calculate interval SIITs 
  #+ "01-SIITIntClean"
  response <- "sPSA.RGR.31to35"
  SIIT <- paste(response, "SIIT", sep=".")
  responses.SIITinterval <- as.vector(outer("sPSA.RGR", suffices, paste, sep="."))
  
  cart.SIIT.dat <- twoLevelOpcreate(data = cart.dat, responses = responses.SIITinterval, 
                                    suffices.treatment=c("C","S"), 
                                    suffices.results="SIIT", 
                                    columns.suffixed="Snapshot.ID.Tag")
  testthat::expect_equal(nrow(cart.SIIT.dat), 10)
  testthat::expect_equal(ncol(cart.SIIT.dat), 23)
 
})
