#devtools::test("asremlPlus")
context("model_selection")

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
  longi.prime.dat <- longitudinalPrime(data=raw.dat, smarthouse.lev=1, 
                                       idcolumns = c("Genotype.ID", "Treatment.1", "Treatment.2"))
  
  longi.dat <- designFactors(longi.prime.dat, insertName = "xDays",
                             nzones = 1, nlanesperzone = 1, nmainplotsperlane = 10, 
                             designfactorMethod="StandardOrder")
  testthat::expect_equal(nrow(longi.dat), 280)
  testthat::expect_equal(ncol(longi.dat), 44)
  
  
  #'## Particular edits to longi.dat
  longi.dat <- within(longi.dat, 
                      { 
                        Days.after.Salting <- as.numfac(Days) - 29
                      })
  
  #'# Step 3: Form derived traits that result in a value for each observation
  #'### Set responses
  responses.image <- c("Area")
  responses.smooth <- paste(responses.image, "smooth", sep=".")
  
  #'## Form growth rates for each observation of a subset of responses by differencing
  longi.dat <- splitContGRdiff(longi.dat, responses.image, 
                               INDICES="Snapshot.ID.Tag",
                               which.rates = c("AGR","RGR"))
  testthat::expect_equal(nrow(longi.dat), 280)
  testthat::expect_equal(ncol(longi.dat), 48)
  
  #'## Form Area.WUE 
  longi.dat <- within(longi.dat, 
                      { 
                        Area.WUE <- WUI(Area.AGR*Days.diffs, Water.Loss)
                      })
  
  #'## Add cumulative responses 
  longi.dat <- within(longi.dat, 
                      { 
                        Water.Loss.Cum <- unlist(by(Water.Loss, Snapshot.ID.Tag, 
                                                    cumulate, exclude.1st=TRUE))
                        WUE.cum <- Area / Water.Loss.Cum 
                      })
  testthat::expect_equal(nrow(longi.dat), 280)
  testthat::expect_equal(ncol(longi.dat), 51)
  
  #'# Step 4: Fit splines to smooth the longitudinal trends in the primary traits and calculate their growth rates
  #'
  #'## Smooth responses
  #+
  for (response in c(responses.image, "Water.Loss"))
    longi.dat <- splitSplines(longi.dat, response, x="xDays", INDICES = "Snapshot.ID.Tag", 
                              df = 4)
  longi.dat <- with(longi.dat, longi.dat[order(Snapshot.ID.Tag, xDays), ])
  testthat::expect_equal(nrow(longi.dat), 280)
  testthat::expect_equal(ncol(longi.dat), 53)
  
  #'## Loop over smoothed responses, forming growth rates by differences
  #+
  responses.GR <- paste(responses.smooth, "AGR", sep=".")
  longi.dat <- splitContGRdiff(longi.dat, responses.smooth, 
                               INDICES="Snapshot.ID.Tag",
                               which.rates = c("AGR","RGR"))
  testthat::expect_equal(nrow(longi.dat), 280)
  testthat::expect_equal(ncol(longi.dat), 55)
  
  #'## Finalize longi.dat
  longi.dat <- with(longi.dat, longi.dat[order(Snapshot.ID.Tag, xDays), ])
  
  #'# Step 5: Do exploratory plots on unsmoothed and smoothed longitudinal data

  #'# Step 6: Form single-value plant responses in Snapshot.ID.Tag order.
  #'
  #'## 6a) Set up a data frame with factors only
  #+
  cart.dat <- longi.dat[longi.dat$Days == 31, 
                        c("Smarthouse","Lane","Position","Snapshot.ID.Tag",
                          "xPosn","xMainPosn",
                          "Zones","xZones","SHZones","ZLane","ZMainplots", "Subplots",
                          "Genotype.ID","Treatment.1")]
  cart.dat <- cart.dat[do.call(order, cart.dat), ]
  
  #'## 6b) Get responses based on first and last date.
  #'
  #'### Observation for first and last date
  cart.dat <- cbind(cart.dat, getTimesSubset(responses.image, data = longi.dat, 
                                             which.times = c(31), suffix = "first"))
  cart.dat <- cbind(cart.dat, getTimesSubset(responses.image, data = longi.dat, 
                                             which.times = c(42), suffix = "last"))
  cart.dat <- cbind(cart.dat, getTimesSubset(c("WUE.cum"), 
                                             data = longi.dat, 
                                             which.times = c(42), suffix = "last"))
  responses.smooth <- paste(responses.image, "smooth", sep=".")
  cart.dat <- cbind(cart.dat, getTimesSubset(responses.smooth, data = longi.dat, 
                                             which.times = c(31), suffix = "first"))
  cart.dat <- cbind(cart.dat, getTimesSubset(responses.smooth, data = longi.dat, 
                                             which.times = c(42), suffix = "last"))
  testthat::expect_equal(nrow(cart.dat), 20)
  testthat::expect_equal(ncol(cart.dat), 19)
  
  #'### Growth rates over whole period.
  #+
  tottime <- 42 - 31
  cart.dat <- within(cart.dat, 
                     { 
                       Area.AGR <- (Area.last - Area.first)/tottime
                       Area.RGR <- log(Area.last / Area.first)/tottime
                     })
  
  #'### Calculate water index over whole period
  tmp.dat <- cart.dat
  cart.dat <- merge(cart.dat, 
                    intervalWUI("Area", water.use = "Water.Loss", 
                                start.times = c(31), 
                                end.times = c(42), 
                                suffix = NULL, 
                                data = longi.dat, include.total.water = TRUE),
                    by = c("Snapshot.ID.Tag"))
  names(cart.dat)[match(c("Area.WUI","Water.Loss.Total"),names(cart.dat))] <- c("Area.Overall.WUE", 
                                                                                "Water.Loss.Overall")
  cart.dat$Water.Loss.rate.Overall <- cart.dat$Water.Loss.Overall / (42 - 31)
  testthat::expect_equal(nrow(cart.dat), 20)
  testthat::expect_equal(ncol(cart.dat), 25)
  #Check same value of the WUE for include.total.water both TRUE and FALSE
  tmp.dat <- merge(tmp.dat, 
                   intervalWUI("Area", water.use = "Water.Loss", 
                               start.times = c(31), 
                               end.times = c(42), 
                               suffix = NULL, 
                               data = longi.dat, include.total.water = FALSE),
                   by = c("Snapshot.ID.Tag"))
  names(tmp.dat)[match(c("Area.WUI"),names(tmp.dat))] <- "Area.Overall.WUE"
  testthat::expect_true(all(abs(cart.dat$Area.Overall.WUE - tmp.dat$Area.Overall.WUE) < 1e-05))
  
  #'## 6c) Add growth rates and water indices for intervals
  #'### Set up intervals
  #+
  start.days <- list(31,35,31,38)
  end.days <- list(35,38,38,42)
  suffices <- list("31to35","35to38","31to38","38to42")
  
  #'### Rates for specific intervals from the smoothed data by differencing
  #+
  for (r in responses.smooth)
  { for (k in 1:length(suffices))
  { 
    cart.dat <- merge(cart.dat, 
                      intervalGRdiff(r, 
                                     which.rates = c("AGR","RGR"), 
                                     start.times = start.days[k][[1]], 
                                     end.times = end.days[k][[1]], 
                                     suffix.interval = suffices[k][[1]], 
                                     data = longi.dat),
                      by = "Snapshot.ID.Tag")
  }
  }
  
  #'### Water indices for specific intervals from the unsmoothed and smoothed data
  #+
  for (k in 1:length(suffices))
  { 
    cart.dat <- merge(cart.dat, 
                      intervalWUI("Area", water.use = "Water.Loss", 
                                  start.times = start.days[k][[1]], 
                                  end.times = end.days[k][[1]], 
                                  suffix = suffices[k][[1]], 
                                  data = longi.dat, include.total.water = TRUE),
                      by = "Snapshot.ID.Tag")
    names(cart.dat)[match(paste("Area.WUI", suffices[k][[1]], sep="."), 
                          names(cart.dat))] <- paste("Area.WUE", suffices[k][[1]], sep=".")
    cart.dat[paste("Water.Loss.rate", suffices[k][[1]], sep=".")] <- 
      cart.dat[[paste("Water.Loss.Total", suffices[k][[1]], sep=".")]] / 
      ( end.days[k][[1]] - start.days[k][[1]])
  }
  
  cart.dat <- with(cart.dat, cart.dat[order(Snapshot.ID.Tag), ])
  testthat::expect_equal(nrow(cart.dat), 20)
  testthat::expect_equal(ncol(cart.dat), 49)
  
  #'# Step 7: Form continuous and interval SIITs
  #'
  #'## 7a) Calculate continuous
  #+
  cols.retained <-  c("Snapshot.ID.Tag","Smarthouse","Lane","Position",
                      "Days","Snapshot.Time.Stamp", "Hour", "xDays",
                      "Zones","xZones","SHZones","ZLane","ZMainplots",
                      "xMainPosn", "Genotype.ID")
  responses.GR <- c("Area.smooth.AGR","Area.smooth.AGR","Area.smooth.RGR")
  suffices.results <- c("diff", "SIIT", "SIIT")
  responses.SIIT <- unlist(Map(paste, responses.GR, suffices.results,sep="."))
  
  longi.SIIT.dat <- 
    twoLevelOpcreate(responses.GR, longi.dat, suffices.treatment=c("C","S"),
                     operations = c("-", "/", "/"), suffices.results = suffices.results, 
                     columns.retained = cols.retained, 
                     by = c("Smarthouse","Zones","ZMainplots","Days"))
  longi.SIIT.dat <- with(longi.SIIT.dat, 
                         longi.SIIT.dat[order(Smarthouse,Zones,ZMainplots,Days),])
  testthat::expect_equal(nrow(longi.SIIT.dat), 140)
  testthat::expect_equal(ncol(longi.SIIT.dat), 22)
  
  #' ### Plot SIIT profiles 
  #' 
  #+ "03-SIITProfiles"
  k <- 2
  nresp <- length(responses.SIIT)
  limits <- with(longi.SIIT.dat, list(c(min(Area.smooth.AGR.diff, na.rm=TRUE),
                                        max(Area.smooth.AGR.diff, na.rm=TRUE)),
                                      c(0,3),
                                      c(0,1.5)))
  #Plots

  #'## 7b) Calculate interval SIITs 
  #+ "01-SIITIntClean"
  suffices <- list("31to35","35to38","31to38","38to42")
  response <- "Area.smooth.RGR.31to35"
  SIIT <- paste(response, "SIIT", sep=".")
  responses.SIITinterval <- as.vector(outer("Area.smooth.RGR", suffices, paste, sep="."))
  
  cart.SIIT.dat <- twoLevelOpcreate(responses.SIITinterval, cart.dat,
                                    suffices.treatment=c("C","S"), 
                                    suffices.results="SIIT", 
                                    columns.suffixed="Snapshot.ID.Tag")
  testthat::expect_equal(nrow(cart.SIIT.dat), 10)
  testthat::expect_equal(ncol(cart.SIIT.dat), 23)
 
})
