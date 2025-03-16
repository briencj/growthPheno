cat("#### Test prepImageData\n")
test_that("prepImageData", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(dae)
  library(ggplot2)
  library(growthPheno)
  
  # A set of RGB images with all names using defaults
  raw.RGB.dat <- suppressWarnings(importExcel(file = "./data/rawRGBdatarow.csv",
                                              timeAfterStart = "Time.after.Plantind..d.", 
                                              startTime = "09/01/2017 0:00 AM",
                                              timeFormat = "%d/%m/%Y %H:%M",
                                              plotImagetimes = FALSE))
  testthat::expect_true(all(names(raw.RGB.dat)[c(18,56,94)] == 
                              c("Area.SV1", "Area.SV2", "Area.TV")))
  raw.RGB.dat <- rbind(raw.RGB.dat, raw.RGB.dat)
  longi.RGB <- prepImageData(data = raw.RGB.dat, 
                             timeAfterStart = "Time.after.Plantind..d.", 
                             idcolumns = c("Genotype.ID", "Treatment"))
  testthat::expect_equal(nrow(longi.RGB), 2)
  testthat::expect_equal(ncol(longi.RGB), 35)
  testthat::expect_true("Reps" %in% names(longi.RGB))
  
  # Area RGB images with some imaging traits and copying Water traits
  raw.RGB.dat <- suppressWarnings(importExcel(file = "./data/rawRGBdatarow.csv",
                                              timeAfterStart = "Time.after.Plantind..d.", 
                                              startTime = "09/01/2017 0:00 AM",
                                              timeFormat = "%d/%m/%Y %H:%M",
                                              plotImagetimes = FALSE))
  testthat::expect_true(all(names(raw.RGB.dat)[c(18,56,94)] == 
                              c("Area.SV1", "Area.SV2", "Area.TV")))
  raw.RGB.dat <- rbind(raw.RGB.dat, raw.RGB.dat)
  longi.RGB <- prepImageData(data = raw.RGB.dat, 
                             timeAfterStart = "Time.after.Plantind..d.", 
                             idcolumns = c("Genotype.ID", "Treatment"),
                             calcWaterUse = FALSE,
                             traits = list(img = c("Area", "Compactness"),
                                           H2O = c("Weight.Before","Weight.After",
                                                   "Water.Amount")),
                             labsCamerasViews = list(img = c("SV1", "SV2", "TV"),
                                                     H2O = NULL))
  testthat::expect_equal(nrow(longi.RGB), 2)
  testthat::expect_equal(ncol(longi.RGB), 21)
  
  #Check cameraType = FLUO and keepCameraType
  raw.FLUO.dat <- suppressWarnings(importExcel(file = "./data/rawFLUOdatarow.csv",
                                               timeAfterStart = "Time.after.Plantind..d.", 
                                               startTime = "09/01/2017 0:00 AM", 
                                               timeFormat = "%d/%m/%Y %H:%M",
                                               plotImagetimes = FALSE,
                                               cameraType = "FLUO",
                                               keepCameraType = TRUE))
  testthat::expect_true(all(names(raw.FLUO.dat)[c(18,54)] == 
                              c("Area.FLUO_SV1", "Area.FLUO_SV2")))
  raw.FLUO.dat <- rbind(raw.FLUO.dat, raw.FLUO.dat)
  longi.FLUO <- prepImageData(data = raw.FLUO.dat, 
                              timeAfterStart = "Time.after.Plantind..d.", 
                              idcolumns = c("Genotype.ID", "Treatment"),
                              traits =    list(all = "Area"), 
                              labsCamerasViews = list(all = c("FLUO_SV1", "FLUO_SV2")))
  testthat::expect_equal(nrow(longi.FLUO), 2)
  testthat::expect_equal(ncol(longi.FLUO), 18)
  testthat::expect_lt(abs(longi.FLUO$PSA[1] - 13.43), 1e-02)
  testthat::expect_lt(abs(longi.FLUO$PSA.FLUO_SV1[1] - 6.41), 1e-02)
  
  
  camview.labels <- c("SF0", "SL0", "SU0", "TV0")
  names(camview.labels) <- c("RGB_Side_Far_0", "RGB_Side_Lower_0", 
                             "RGB_Side_Upper_0", "RGB_TV_0")
  
  #Test name change with move to suffix  and supply characters instead of lists
  raw.19.dat <- suppressWarnings(importExcel(file = "./data/raw19datarow.csv",
                                             individualId = "Snapshot.ID.Tags",
                                             startTime = "06/10/2017 0:00 AM",
                                             timeFormat = "%d/%m/%Y %H:%M",
                                             labsCamerasViews = camview.labels, 
                                             plotImagetimes = FALSE))
  testthat::expect_true(all(c("Area.SF0", "Area.SL0", "Area.SU0", "Area.TV0") %in% 
                              names(raw.19.dat)))
  raw.19.dat <- rbind(raw.19.dat, raw.19.dat)
  longi.19 <- prepImageData(data = raw.19.dat, 
                            individualId = "Snapshot.ID.Tags",
                            idcolumns = c("Plant.Species", "Mycorrhiza", "Zn"),
                            calcWaterUse = FALSE, 
                            traits = "Area", 
                            labsCamerasViews = camview.labels)
  testthat::expect_equal(nrow(longi.19), 2)
  testthat::expect_equal(ncol(longi.19), 17)
  testthat::expect_lt(abs(longi.19$PSA[1] - 6.73), 1e-02)
  testthat::expect_lt(abs(longi.19$PSA.SF0[1] - 1.18), 1e-02)
  
  #Test name change with move to suffix  and pass through just the names to retain
  raw.19.dat <- suppressWarnings(importExcel(file = "./data/raw19datarow.csv",
                                             individualId = "Snapshot.ID.Tags",
                                             startTime = "06/10/2017 0:00 AM",
                                             timeFormat = "%d/%m/%Y %H:%M",
                                             labsCamerasViews = camview.labels, 
                                             plotImagetimes = FALSE))
  testthat::expect_true(all(c("Area.SF0", "Area.SL0", "Area.SU0", "Area.TV0") %in% 
                              names(raw.19.dat)))
  raw.19.dat <- rbind(raw.19.dat, raw.19.dat)
  longi.19 <- prepImageData(data = raw.19.dat, 
                            individualId = "Snapshot.ID.Tags",
                            idcolumns = c("Plant.Species", "Mycorrhiza", "Zn"),
                            calcWaterUse = FALSE, 
                            traits = c("Area.SF0", "Area.SL0", "Area.SU0", "Area.TV0"), 
                            labsCamerasViews = NULL)
  testthat::expect_equal(nrow(longi.19), 2)
  testthat::expect_equal(ncol(longi.19), 17)
  testthat::expect_lt(abs(longi.19$PSA[1] - 6.73), 1e-02)
  testthat::expect_lt(abs(longi.19$PSA.SF0[1] - 1.18), 1e-02)
  
  #Test remove cameraType with move to suffix
  raw.19.dat <- suppressWarnings(importExcel(file = "./data/raw19datarow.csv",
                                             individualId = "Snapshot.ID.Tags",
                                             startTime = "06/10/2017 0:00 AM",
                                             timeFormat = "%d/%m/%Y %H:%M",
                                             cameraType = "RGB", 
                                             plotImagetimes = FALSE))
  testthat::expect_true(all(c("Area.Side_Far_0", "Area.Side_Lower_0", "Area.Side_Upper_0", 
                              "Area.TV_0") %in% names(raw.19.dat)))
  raw.19.dat <- rbind(raw.19.dat, raw.19.dat)
  longi.19 <- prepImageData(data = raw.19.dat, 
                            individualId = "Snapshot.ID.Tags",
                            idcolumns = c("Plant.Species", "Mycorrhiza", "Zn"),
                            calcWaterUse = FALSE, 
                            traits = list(t = "Area"), 
                            labsCamerasViews = list(c = c("Side_Far_0", "Side_Lower_0", 
                                                          "Side_Upper_0", "TV_0")))
  testthat::expect_equal(nrow(longi.19), 2)
  testthat::expect_equal(ncol(longi.19), 17)
  testthat::expect_lt(abs(longi.19$PSA[1] - 6.73), 1e-02)
  testthat::expect_lt(abs(longi.19$PSA.Side_Far_0[1] - 1.18), 1e-02)
  
})

cat("#### Test prepImageData Replicates\n")
test_that("prepImageData Replicates", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(dae)
  library(ggplot2)
  library(growthPheno)
  
  # A set of RGB images with all names using defaults
  raw.RGB.dat <- suppressWarnings(importExcel(file = "./data/rawRGBdatarow.csv",
                                              timeAfterStart = "Time.after.Plantind..d.", 
                                              startTime = "09/01/2017 0:00 AM",
                                              timeFormat = "%d/%m/%Y %H:%M",
                                              plotImagetimes = FALSE))
  testthat::expect_true(all(names(raw.RGB.dat)[c(18,56,94)] == 
                              c("Area.SV1", "Area.SV2", "Area.TV")))
  raw.RGB.dat <- rbind(raw.RGB.dat, raw.RGB.dat)
  longi.RGB.Rep <- prepImageData(data = raw.RGB.dat, 
                                 timeAfterStart = "Time.after.Plantind..d.", 
                                 potIDcolumns = c("Genotype.ID", "Treatment", "Replicate"))
  testthat::expect_equal(nrow(longi.RGB.Rep), 2)
  testthat::expect_equal(ncol(longi.RGB.Rep), 35)
  testthat::expect_true(!("Reps" %in% names(longi.RGB.Rep)))
  testthat::expect_true(is.factor(longi.RGB.Rep$Replicate))
  
  longi.RGB.Rep <- prepImageData(data = raw.RGB.dat, 
                                 timeAfterStart = "Time.after.Plantind..d.", 
                                 potIDcolumns = c("Genotype.ID", "Treatment"))
  testthat::expect_equal(nrow(longi.RGB.Rep), 2)
  testthat::expect_equal(ncol(longi.RGB.Rep), 34)
  testthat::expect_true(!("Reps" %in% names(longi.RGB.Rep)))
  testthat::expect_true(!("Replicate" %in% names(longi.RGB.Rep)))

  #Test when no Replicate in data.frame  
  tmp <- raw.RGB.dat[, -match("Replicate", names(raw.RGB.dat))]
  testthat::expect_error(
    longi.RGB.Rep <- prepImageData(data = tmp, 
                                   timeAfterStart = "Time.after.Plantind..d.", 
                                   potIDcolumns = c("Genotype.ID", "Treatment", "Replicate")),
    regexp = "The following variables are not present in data:  Replicate")
  
})


cat("#### Test prepImageData exampleData\n")
test_that("prepImageData exampleData", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(growthPheno)

  data(exampleData)  
  
  longi.dat <- prepImageData(data=raw.dat, 
                           traits= list(all = "Area"), 
                           labsCamerasViews = list(all = NULL), 
                           potIDcolumns = c("Genotype.ID", "Treatment.1", "Replicate" ), 
                           smarthouse.lev=1)
  testthat::expect_true(all(c("Genotype.ID", "Treatment.1", "Replicate" ) %in% names(longi.dat)))
  testthat::expect_true(!("Reps" %in% names(longi.dat)))
  testthat::expect_true(!("Treatment.2" %in% names(longi.dat)))
  
})