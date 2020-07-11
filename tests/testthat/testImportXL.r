cat("#### Test importExcel arguments\n")
test_that("importXLRename", {
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
  
  #Check keepCameraType
  raw.RGB.dat <- importExcel(file = "./data/rawRGBdatarow.csv",
                             timeAfterStart = "Time.after.Plantind..d.", 
                             startTime = "09/01/2017 0:00 AM", 
                             timeFormat = "%d/%m/%Y %H:%M",
                             plotImagetimes = FALSE,
                             keepCameraType = TRUE)
  testthat::expect_true(all(names(raw.RGB.dat)[c(18,56,94)] == 
                              c("Area.RGB_SV1", "Area.RGB_SV2", "Area.RGB_TV")))
  #Check cameraType = FLUO
  raw.FLUO.dat <- suppressWarnings(importExcel(file = "./data/rawFLUOdatarow.csv",
                                               timeAfterStart = "Time.after.Plantind..d.", 
                                               startTime = "09/01/2017 0:00 AM", 
                                               timeFormat = "%d/%m/%Y %H:%M",
                                               plotImagetimes = FALSE,
                                               cameraType = "FLUO"))
  testthat::expect_true(all(names(raw.FLUO.dat)[c(18,54)] == 
                              c("Area.SV1", "Area.SV2")))
  
  #Check cameraType = FLUO and keepCameraType
  raw.FLUO.dat <- suppressWarnings(importExcel(file = "./data/rawFLUOdatarow.csv",
                                               timeAfterStart = "Time.after.Plantind..d.", 
                                               startTime = "09/01/2017 0:00 AM", 
                                               timeFormat = "%d/%m/%Y %H:%M",
                                               plotImagetimes = FALSE,
                                               cameraType = "FLUO", keepCameraType = TRUE))
  testthat::expect_true(all(names(raw.FLUO.dat)[c(18,54)] == 
                              c("Area.FLUO_SV1", "Area.FLUO_SV2")))
  
  #Check no cameraType and keepCameraType (can only move if can recognize camera type)
  testthat::expect_warning(raw.FLUO.dat <- importExcel(file = "./data/rawFLUOdatarow.csv",
                                                       timeAfterStart = "Time.after.Plantind..d.", 
                                                       startTime = "09/01/2017 0:00 AM", 
                                                       timeFormat = "%d/%m/%Y %H:%M",
                                                       plotImagetimes = FALSE, 
                                                       keepCameraType = TRUE))
  testthat::expect_true(all(names(raw.FLUO.dat)[c(18,54)] == 
                              c("FLUO_SV1.Area", "FLUO_SV2.Area")))
  
  #Test do nothing
  raw.19.dat <- suppressWarnings(importExcel(file = "./data/raw19datarow.csv",
                                             cartId = "Snapshot.ID.Tags",
                                             startTime = "06/10/2017 0:00 AM",
                                             timeFormat = "%d/%m/%Y %H:%M",
                                             plotImagetimes = FALSE, 
                                             prefix2suffix = FALSE, keepCameraType = TRUE))
  testthat::expect_true(all(c("RGB_Side_Far_0.Area", "RGB_Side_Lower_0.Area", 
                              "RGB_Side_Upper_0.Area", "RGB_TV_0.Area") %in% names(raw.19.dat)))
  
  #Test do nothing except remove cameraType
  raw.19.dat <- suppressWarnings(importExcel(file = "./data/raw19datarow.csv",
                                             cartId = "Snapshot.ID.Tags",
                                             startTime = "06/10/2017 0:00 AM",
                                             timeFormat = "%d/%m/%Y %H:%M",
                                             prefix2suffix = FALSE,
                                             plotImagetimes = FALSE))
  testthat::expect_true(all(c("Side_Far_0.Area", "Side_Lower_0.Area", 
                              "Side_Upper_0.Area", "TV_0.Area") %in% names(raw.19.dat)))
  
  camview.labels <- c("SF0", "SL0", "SU0", "TV0")
  names(camview.labels) <- c("RGB_Side_Far_0", "RGB_Side_Lower_0", 
                             "RGB_Side_Upper_0", "RGB_TV_0")
  
  #Test name change only
  raw.19.dat <- suppressWarnings(importExcel(file = "./data/raw19datarow.csv",
                                             cartId = "Snapshot.ID.Tags",
                                             startTime = "06/10/2017 0:00 AM",
                                             timeFormat = "%d/%m/%Y %H:%M",
                                             prefix2suffix = FALSE, 
                                             labsCamerasViews = camview.labels, 
                                             plotImagetimes = FALSE))
  testthat::expect_true(all(c("SF0.Area", "SL0.Area", "SU0.Area", "TV0.Area") %in% 
                              names(raw.19.dat)))
  
  #Test name change with move to suffix
  raw.19.dat <- suppressWarnings(importExcel(file = "./data/raw19datarow.csv",
                                             cartId = "Snapshot.ID.Tags",
                                             startTime = "06/10/2017 0:00 AM",
                                             timeFormat = "%d/%m/%Y %H:%M",
                                             labsCamerasViews = camview.labels, 
                                             plotImagetimes = FALSE))
  testthat::expect_true(all(c("Area.SF0", "Area.SL0", "Area.SU0", "Area.TV0") %in% 
                              names(raw.19.dat)))
  
  #Test remove cameraType with move to suffix
  raw.19.dat <- suppressWarnings(importExcel(file = "./data/raw19datarow.csv",
                                             cartId = "Snapshot.ID.Tags",
                                             startTime = "06/10/2017 0:00 AM",
                                             timeFormat = "%d/%m/%Y %H:%M",
                                             cameraType = "RGB", 
                                             plotImagetimes = FALSE))
  testthat::expect_true(all(c("Area.Side_Far_0", "Area.Side_Lower_0", "Area.Side_Upper_0", 
                              "Area.TV_0") %in% names(raw.19.dat)))
  
  #Test name change with retain cameraType and  move to suffix
  camview.labels <- c("RGB_SF0", "RGB_SL0", "RGB_SU0", "RGB_TV0")
  names(camview.labels) <- c("RGB_Side_Far_0", "RGB_Side_Lower_0", "RGB_Side_Upper_0", 
                             "RGB_TV_0")
  raw.19.dat <- suppressWarnings(importExcel(file = "./data/raw19datarow.csv",
                                             cartId = "Snapshot.ID.Tags",
                                             startTime = "06/10/2017 0:00 AM",
                                             timeFormat = "%d/%m/%Y %H:%M",
                                             labsCamerasViews = camview.labels, 
                                             keepCameraType = TRUE, 
                                             plotImagetimes = FALSE))
  testthat::expect_true(all(c("Area.RGB_SF0", "Area.RGB_SL0", "Area.RGB_SU0", 
                              "Area.RGB_TV0") %in% names(raw.19.dat)))
  
  #Test name change with remove cameraType and  move to suffix
  camview.labels <- c("RGB_SF0", "RGB_SL0", "RGB_SU0", "RGB_TV0")
  names(camview.labels) <- c("RGB_Side_Far_0", "RGB_Side_Lower_0", "RGB_Side_Upper_0", 
                             "RGB_TV_0")
  raw.19.dat <- suppressWarnings(importExcel(file = "./data/raw19datarow.csv",
                                             cartId = "Snapshot.ID.Tags",
                                             startTime = "06/10/2017 0:00 AM",
                                             timeFormat = "%d/%m/%Y %H:%M",
                                             labsCamerasViews = camview.labels, 
                                             plotImagetimes = FALSE))
  testthat::expect_true(all(c("Area.SF0", "Area.SL0", "Area.SU0", 
                              "Area.TV0") %in% names(raw.19.dat)))
})


cat("#### Test importExcel times\n")
test_that("importXLRename", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(dae)
  library(ggplot2)
  library(growthPheno)
  
  camera.names <- c("SV1", "SV2", "TV")
  names(camera.names) <- paste(c("RGB_Side_Far_0", "RGB_Side_Far_4", "RGB_TV"), "Result", sep=".")
  
  #Cannot test xlsx files because these cannot be included in a package
  #Test for csv file  
  #The warning occurs because the first Snaopshot.ID.Tag does not have a DAP15
  testthat::expect_warning(raw.dat <- importExcel(file = "./data/0498_Rice_Plett_test.csv",
                                                  startTime = "27/02/2020 0:00 AM",
                                                  timeFormat = "%d/%m/%Y %H:%M",
                                                  labsCamerasViews = camera.names))
  testthat::expect_true(all(dim(raw.dat) == c(863, 26)))
  testthat::expect_true(all(c("area.SV1", "area.SV2", "area.TV") %in% names(raw.dat)))
  testthat::expect_true(all(table(raw.dat$Time.after.Planting..d.) == c(215,216,216,216)))
  
  #Non-matching time formats is a problem for csv files
  testthat::expect_warning(testthat::expect_error(
    raw.dat <- importExcel(file = "./data/0498_Rice_Plett test.csv",
                           startTime = "2020-02-27 0:00 AM",
                           timeFormat = "%Y-%m-%d %H:%M",
                           labsCamerasViews = camera.names)))
})