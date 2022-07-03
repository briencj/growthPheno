cat("#### Test smoothSpline using NCSS with leaf data when there are missing values\n")
test_that("leaf_smoothSpline", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(dae)
  library(ggplot2)
  library(growthPheno)
  
  
  # A small subset of Exp 270 leaf tracking data
  data(testSpline)
  responses <- names(test)[5:ncol(test)]
  
  ##Test omit in fitSpline - Length 3
  leaf.dat <- test
  carts <- levels(leaf.dat$Snapshot.ID.Tag)
  nrows <- list(6,4,5,3,2,6,1,1)
  names(nrows) <- carts
  fit <- list()
  for (cart in carts)
  {
    fit[[cart]] <- smoothSpline(subset(leaf.dat, Snapshot.ID.Tag == cart), 
                                response = "Length.3", response.smoothed = "sLength.3", 
                                x="xDays", 
                                df = 4, na.x.action = "omit", na.y.action = "omit", 
                                rates = c("AGR", "RGR"), 
                                suffices.rates = c("AGRdv", "RGRdv"))
    testthat::expect_equal(nrows[[cart]], nrow(fit[[cart]]$predictions))
    testthat::expect_true(all(unlist(c("xDays", "sLength.3", "sLength.3.AGRdv", "sLength.3.RGRdv") %in% 
                                       names(fit[[cart]]$predictions))))
  }

  ##Test omit in snoothSpline - Length 2 with a 0 length data.frame
  nrows <- list(11,12,12,12,9,12,0,9)
  names(nrows) <- carts
  fit <- list()
  for (cart in carts)
  {
    fit[[cart]] <- smoothSpline(subset(leaf.dat, Snapshot.ID.Tag == cart), 
                                response = "Length.2", response.smoothed = "sLength.2", 
                                x="xDays", 
                                df = 4, na.x.action = "omi", na.y.action = "omit", 
                                rates = c("AGR", "RGR"), 
                                suffices.rates = c("AGRdv", "RGRdv"))
    testthat::expect_equal(nrows[[cart]], nrow(fit[[cart]]$predictions))
    testthat::expect_equal(ncol(fit[[cart]]$predictions), 4)
    testthat::expect_true(all(unlist(c("xDays", "sLength.2", "sLength.2.AGRdv", "sLength.2.RGRdv") %in% 
                                       names(fit[[cart]]$predictions))))
  }
  
  ##Test omit in smoothSpline - Length 2 with a 0 length data.frame
  leaf.dat <- test
  carts <- levels(leaf.dat$Snapshot.ID.Tag)
  nrows <- list(6,4,5,3,2,6,1,1)
  names(nrows) <- carts
  fit <- list()
  for (cart in carts)
  {
    fit[[cart]] <- smoothSpline(subset(leaf.dat, Snapshot.ID.Tag == cart), 
                                response = "Length.3", response.smoothed = "sLength.3", 
                                x="xDays", correctBoundaries = FALSE,
                                df = 4, na.x.action = "omit", na.y.action = "omit")
    testthat::expect_equal(ncol(fit[[cart]]$predictions), 2)
    testthat::expect_equal(nrow(fit[[cart]]$predictions), nrows[[cart]])
    testthat::expect_true(all(unlist(c("xDays", "sLength.3") %in% names(fit[[cart]]$predictions))))
  }  
  nrows <- list(6,4,5,3,2,6,1,1)
  names(nrows) <- carts
  fitC <- list()
  for (cart in carts)
  {
    fitC[[cart]] <- smoothSpline(subset(leaf.dat, Snapshot.ID.Tag == cart), 
                                 response = "Length.3", 
                                 x="xDays", correctBoundaries = TRUE,
                                 df = 4, na.x.action = "omit", na.y.action = "omit")
    testthat::expect_equal(ncol(fitC[[cart]]$predictions), 2)
    testthat::expect_equal(nrow(fitC[[cart]]$predictions), nrows[[cart]])
    testthat::expect_true(all(unlist(c("xDays", "sLength.3") %in% names(fit[[cart]]$predictions))))
  }  
  testthat::expect_true(all(abs(fit[["047162-C"]]$sLength.3 - 
                                  fitC[["047162-C"]]$sLength.3) > 0.01))
  testthat::expect_true(all(abs(fit[["047164-S"]]$sLength.3 - 
                                  fitC[["047164-S"]]$sLength.3) < 1e-05))
})
