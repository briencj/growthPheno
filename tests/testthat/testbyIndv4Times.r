
cat("#### Test spline fitting using NCSS with leaf data when there are missing values\n")
test_that("leaf_growthPheno", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(dae)
  library(ggplot2)
  library(growthPheno)
  
  
  # A small subset of Exp 270 leaf tracking data
  data(testSpline)
  #Set Days to NA to match xDays
  test$Days[is.na(test$xDays)] <- NA
  responses <- names(test)[5:ncol(test)]
  
  ##Test the fail options for NCSS
  testthat::expect_error(byIndv4Times_SplinesGRs(data = test, response = "Area", times="Days", 
                                                 individuals = "Snapshot.ID.Tag", 
                                                 df = 4, na.x.action = "fail"))
  testthat::expect_error(byIndv4Times_SplinesGRs(data = test, response = "Length.1", 
                                                 individuals = "Snapshot.ID.Tag", times="Days", 
                                                 df = 4, na.y.action = "fail"))
  
  ##Fit some splines - exclude y
  leaf.dat <- test
  testthat::expect_warning(
    for (resp in responses)
      leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, response = resp, 
                                          times="Days", df = 4),
    regexp = "Need at least 4 distinct x values to fit a spline - all fitted values set to NA")
  testthat::expect_equal(sum(!is.na(leaf.dat$Length.1)), 85)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.1)), 86)
  testthat::expect_equal(sum(is.na(leaf.dat$sArea)), 3)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.3[13:24])), 4)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.3[49:60])), 0)  
  
  ##trim
  leaf.dat <- test
  for (resp in responses)
    leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, response = resp, times = "Days", 
                                        df = 4, na.y.action = "trim")
  testthat::expect_equal(sum(!is.na(leaf.dat$Length.1)), 85)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.1)), 86)
  testthat::expect_lt(abs(leaf.dat$sLength.1[57] - 11.067455), 1e-05)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.3[13:24])), 4)
  
  ##Test ltrimx
  leaf.dat <- test
  for (resp in responses)
    leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, response = resp, times = "Days", 
                                        df = 4, na.y.action = "ltrim", 
                                        rates.method = "deriv", 
                                        suffices.rates = c("AGRdv", "RGRdv"))
  testthat::expect_equal(length(lapply(responses, 
                                       function(resp, data)
                                       {
                                         resp <- paste0("s", resp)
                                         testthat::expect_equal(sum(is.na(data[[resp]])), 
                                                                sum(is.na(data[[paste(resp,"AGRdv",sep = ".")]])),
                                                                sum(is.na(data[[paste(resp,"RGRdv",sep = ".")]])))
                                       }, data = leaf.dat)), 4)
  testthat::expect_equal(nrow(leaf.dat), 96)
  testthat::expect_equal(ncol(leaf.dat), 20)
  testthat::expect_equal(sum(!is.na(leaf.dat$Length.1)), 85)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.1)), 87)
  testthat::expect_lt(abs(leaf.dat$sLength.1[57] - 11.067455), 1e-05)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.3[13:24])), 7)
  leaf.dat.noC <- leaf.dat
  
  ##Test ltrimx with correctBoundaries 
  # - need to use xDays as Days causes correctBoundaries to go singular???
  leaf.dat <- test
  for (resp in responses)
    leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, response = resp, times = "xDays", 
                                        df = 4, rates.method = "none", 
                                        na.y.action = "ltrim", 
                                        correctBoundaries = TRUE)
  testthat::expect_equal(nrow(leaf.dat), 96)
  testthat::expect_equal(ncol(leaf.dat), 12)
  testthat::expect_equal(sum(!is.na(leaf.dat$Length.1)), 85)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.1)), 87)
  testthat::expect_lt(abs(leaf.dat$sLength.1[57] - 11.05895), 1e-05)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.3[13:24])), 7)
  leaf.dat$sLength.1.noC <- leaf.dat.noC$sLength.1
  leaf.dat$sLength.2.noC <- leaf.dat.noC$sLength.2
  leaf.dat$sLength.3.noC <- leaf.dat.noC$sLength.3
  ggplot(leaf.dat, aes(x = xDays, y = Length.1)) + 
    geom_line() + geom_line(aes(x = xDays, sLength.1), colour = "blue") +
    geom_line(aes(x = xDays, sLength.1.noC), colour = "red") +
    facet_wrap(~ Snapshot.ID.Tag)
  
  ##Test utrimx
  leaf.dat <- test
  for (resp in responses)
    leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, response = resp, times = "Days", 
                                        df = 4, rates.method = "none", 
                                        na.y.action = "utrim")
  testthat::expect_equal(nrow(leaf.dat), 96)
  testthat::expect_equal(ncol(leaf.dat), 12)
  testthat::expect_equal(sum(!is.na(leaf.dat$Length.1)), 85)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.1)), 92)
  testthat::expect_lt(abs(leaf.dat$sLength.1[57] - 11.067455), 1e-05)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.3[13:24])), 9)
  
  ##Test utrimx
  leaf.dat <- test
  for (resp in responses)
    leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, response = resp, times = "Days", 
                                        df = 4, rates.method = "none", 
                                        na.y.action = "allx")
  testthat::expect_equal(nrow(leaf.dat), 96)
  testthat::expect_equal(ncol(leaf.dat), 12)
  testthat::expect_equal(sum(!is.na(leaf.dat$Length.1)), 85)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.1)), 93)
  testthat::expect_lt(abs(leaf.dat$sLength.1[57] - 11.067455), 1e-05)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.3[13:24])), 12)
  testthat::expect_equal(sum(is.na(leaf.dat$sLength.3[94:96])), 3)
  
  
  ##Test omit in fitSpline - Length 3
  leaf.dat <- test
  carts <- levels(leaf.dat$Snapshot.ID.Tag)
  fit <- vector(mode = "list", length = 0)
  nrows <- vector(mode = "list", length = length(carts))
  nrows <- list(6,4,5,3,2,6,1,1)
  names(nrows) <- carts
  for (cart in carts)
  {
    fit[[cart]] <- fitSpline(subset(leaf.dat, Snapshot.ID.Tag == cart), 
                             response = "Length.3", response.smoothed = "sLength.3", 
                             x = "xDays", 
                             df = 4, na.x.action = "omi", na.y.action = "omit", 
                             deriv=1, suffices.deriv="AGRdv",  
                             extra.rate = c(RGRdv = "RGR"))
    testthat::expect_equal(nrows[[cart]], nrow(fit[[cart]]$predictions))
  }
  
  ##Test omit in fitSpline - Length 2 with a 0 length data.frame
  fit <- vector(mode = "list", length = 0)
  nrows <- list(11,12,12,12,9,12,0,9)
  names(nrows) <- carts
  for (cart in carts)
  {
    fit[[cart]] <- fitSpline(subset(leaf.dat, Snapshot.ID.Tag == cart), 
                             response = "Length.2", response.smoothed = "sLength.2", 
                             x = "xDays", 
                             df = 4, na.x.action = "omi", na.y.action = "omit", 
                             deriv=1, suffices.deriv="AGRdv",  
                             extra.rate = c(RGRdv = "RGR"))
    testthat::expect_equal(nrows[[cart]], nrow(fit[[cart]]$predictions))
  }
  testthat::expect_equal(ncol(fit[[cart]]$predictions), 4)
  
  ##Test omit in byIndv4Times_SplinesGRs - get full data.frame because of merge in byIndv4Times_SplinesGRs
  leaf.dat <- test
  for (resp in responses)
    leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, response = resp, times = "Days", 
                                        df = 4, rates.method = "none", 
                                        na.x.action = "omi", na.y.action = "omit")
  testthat::expect_equal(nrow(leaf.dat), 96)
  testthat::expect_equal(ncol(leaf.dat), 12)
  testthat::expect_equal(sum(!is.na(leaf.dat$Length.1)), 85)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.1)), 82)
  
  leaf.dat <- test
  for (resp in responses)
    leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, response = resp, times = "Days", 
                                        df = 4, rates.method = "none", na.y.action = "omit")
  
  
  ##Test omit in fitSpline - Length 2 with a 0 length data.frame
  leaf.dat <- test
  carts <- levels(leaf.dat$Snapshot.ID.Tag)
  fit <- vector(mode = "list", length = 0)
  nrows <- list(6,4,5,3,2,6,1,1)
  names(nrows) <- carts
  for (cart in carts)
  {
    fit[[cart]] <- fitSpline(subset(leaf.dat, Snapshot.ID.Tag == cart), 
                             response = "Length.3", response.smoothed = "sLength.3", 
                             x="xDays", correctBoundaries = FALSE,
                             df = 4, na.x.action = "omit", na.y.action = "omit")
    testthat::expect_equal(ncol(fit[[cart]]$predictions), 2)
    testthat::expect_equal(nrow(fit[[cart]]$predictions), nrows[[cart]])
  }  
  fitC <- vector(mode = "list", length = 0)
  nrows <- list(6,4,5,3,2,6,1,1)
  names(nrows) <- carts
  for (cart in carts)
  {
    fitC[[cart]] <- fitSpline(subset(leaf.dat, Snapshot.ID.Tag == cart), 
                              response = "Length.3", response.smoothed = "sLength.3", 
                              x="xDays", correctBoundaries = TRUE,
                              df = 4, na.x.action = "omit", na.y.action = "omit")
    testthat::expect_equal(ncol(fitC[[cart]]$predictions), 2)
    testthat::expect_equal(nrow(fitC[[cart]]$predictions), nrows[[cart]])
  }  
  testthat::expect_true(all(abs(fit[["047162-C"]]$sLength.3 - 
                                  fitC[["047162-C"]]$sLength.3) > 0.01))
  testthat::expect_true(all(abs(fit[["047164-S"]]$sLength.3 - 
                                  fitC[["047164-S"]]$sLength.3) < 1e-05))
})


cat("#### Test byIndv4Times_SplinesGRs using NCSS with leaf data for log-smoothing\n")
test_that("leaf_growthPheno", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(dae)
  library(ggplot2)
  library(growthPheno)
  
  # A small subset of Exp 270 leaf tracking data
  data(testSpline)
  test$Days[is.na(test$xDays)] <- NA
  
  ##Smooth splines using identity and logarithm transformations - exclude y
  leaf.dat <- test
  leaf.dat$Area.log <- log(leaf.dat$Area)
  #Investigate AGR and RGR calculations
  leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, response="Area", times = "Days", 
                                      df = 4, smoothing.method = "log",
                                      rates.method = "deriv", 
                                      suffices.rates = c("AGRdv", "RGRdv"),
                                      na.y.action = "excl")
  names(leaf.dat)[(match(c("sArea", "sArea.AGRdv","sArea.RGRdv"), 
                         names(leaf.dat)))] <- paste(c("sArea", "sArea.AGRdv",
                                                       "sArea.RGRdv"), "log", sep = ".")
  testthat::expect_equal(sum(is.na(leaf.dat$sArea.log)), 3)
  testthat::expect_false(any(abs(leaf.dat$sArea.AGRdv.log[1:3] - 
                                   c(2.993617, 3.682534, 4.568932)) > 1e-03, 
                             na.rm = TRUE))
  testthat::expect_false(any(abs(leaf.dat$sArea.RGRdv.log[1:3] - 
                                   c(0.1878992, 0.1913271, 0.1955440)) > 1e-03, 
                             na.rm = TRUE))
  #Manual calculation of log smooth
  leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, response = "Area.log", times = "Days", 
                                      df = 4, rates.method = "none")
  leaf.dat$sArea.log <- exp(leaf.dat$sArea.log)
  testthat::expect_false(any(abs(leaf.dat$sArea.log - leaf.dat$sArea.log) > 0.1, 
                             na.rm = TRUE))
  testthat::expect_equal(sum(is.na(leaf.dat$sArea.log)), 3)
  
  #identity smoothing scale calculation of AGR and RGR
  leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, response="Area", times = "Days", 
                                      df = 4, rates.method = "deriv", 
                                      suffices.rates = c("AGRdv", "RGRdv"))
  testthat::expect_false(any(abs(leaf.dat$sArea[1:3] - 
                                   c(14.48536, 18.89667, 23.69953)) > 1e-03, 
                             na.rm = TRUE))
  testthat::expect_false(any(abs(leaf.dat$sArea.AGRdv[1:3] - 
                                   c(4.341943, 4.550298, 5.099377)) > 1e-03, 
                             na.rm = TRUE))
  testthat::expect_false(any(abs(leaf.dat$sArea.RGRdv[1:3] - 
                                   c(0.2997469, 0.2407990, 0.2151678)) > 1e-03, 
                             na.rm = TRUE))
  
  ##Smooth splines using log-smoothing and rates.method = "deriv", with second deriv
  leaf.dat <- test
  #Investigate AGR and RGR calculations
  leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, response="Area", times = "Days", 
                                      df = 4, smoothing.method = "log",
                                      rates.method = "deriv", 
                                      suffices.rates = c("AGRdv", "RGRdv"),
                                      extra.derivs = 2, suffices.extra.derivs = "Acc")
  testthat::expect_equal(sum(is.na(leaf.dat$sArea)), 3)
  testthat::expect_false(any(abs(leaf.dat$sArea.AGRdv[1:3] - 
                                   c(2.993617, 3.682534, 4.568932)) > 1e-03, 
                             na.rm = TRUE))
  testthat::expect_false(any(abs(leaf.dat$sArea.RGRdv[1:3] - 
                                   c(0.1878992, 0.1913271, 0.1955440)) > 1e-03, 
                             na.rm = TRUE))

})


cat("#### Test correctBoundaries for NCSS in fitSpline using a single plant from Rice germplasm\n")
test_that("area_correctBoundaries", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(dae)
  library(ggplot2)
  library(growthPheno)
  data(area.dat)
  
  fity <- smooth.spline(area.dat$xDays,area.dat$Area1,cv=FALSE)
  fit <- as.data.frame(fity[1:2])
  fit$AGR <- c(NA, diff(fit$y)/diff(fit$x))
  fit$RGR <- c(NA, diff(log(fit$y))/diff(fit$x))
  
  fit$yC <- fitSpline(area.dat, x = "xDays", 
                      response = "Area1", response.smoothed = "Area1.smooth", 
                      correctBoundaries = TRUE)$predictions$Area1.smooth
  ggplot(fit) + geom_line(aes(x=x, y=y)) + geom_line(aes(x=x, y=yC), colour = "red")
  fit$AGRC <- c(NA, diff(fit$yC)/diff(fit$x))
  fit$RGRC <- c(NA, diff(log(fit$yC))/diff(fit$x))
  ggplot(fit) + geom_line(aes(x=x, y=AGR)) + geom_line(aes(x=x, y=AGRC), colour = "red")
  ggplot(fit) + geom_line(aes(x=x, y=RGR)) + geom_line(aes(x=x, y=RGRC), colour = "red")
  testthat::expect_true(abs(var(fit$yC) - 66391.75) < 1e-02)
  testthat::expect_true(abs(var(fit$AGRC, na.rm = TRUE) - 153.9782) < 1e-03)
  testthat::expect_true(abs(var(fit$RGRC, na.rm = TRUE) - 0.0009603766) < 1e-05)
  
  #specify df
  fit <- area.dat
  fit$Area1.smooth <- fitSpline(area.dat, x = "xDays", 
                                response = "Area1", response.smoothed = "Area1.smooth", 
                                df = 4, 
                                correctBoundaries = FALSE)$predictions$Area1.smooth
  ggplot(fit) + geom_line(aes(x=xDays, y=Area1)) + geom_line(aes(x=xDays, y=Area1.smooth),
                                                             colour = "blue")
  fit$AGR <- c(NA, diff(fit$Area1.smooth)/diff(fit$xDays))
  fit$RGR <- c(NA, diff(log(fit$Area1.smooth))/diff(fit$xDays))
  ggplot(fit) + geom_line(aes(x=xDays, y=Area1.AGR)) + geom_line(aes(x=xDays, y=AGR),colour = "blue")
  ggplot(fit) + geom_line(aes(x=xDays, y=Area1.RGR)) + geom_line(aes(x=xDays, y=RGR),colour = "blue")
  testthat::expect_true(abs(var(fit$Area1.smooth) - 65264.71) < 1e-02)
  testthat::expect_true(abs(var(fit$AGR, na.rm = TRUE) - 123.8462) < 1e-03)
  testthat::expect_true(abs(var(fit$RGR, na.rm = TRUE) - 0.005474027) < 1e-03)
  
  
  #Correct the boundaries
  fit$Area1.smooth.C <- fitSpline(area.dat, x = "xDays", 
                                  response = "Area1", response.smoothed = "Area1.smooth", 
                                  df = 4, 
                                  correctBoundaries = TRUE)$predictions$Area1.smooth
  ggplot(fit) + geom_line(aes(x=xDays, y=Area1)) + 
    geom_line(aes(x=xDays, y=Area1.smooth),colour = "blue") + 
    geom_line(aes(x=xDays, y=Area1.smooth.C), colour = "red")
  fit$AGRC <- c(NA, diff(fit$Area1.smooth.C)/diff(fit$xDays))
  fit$RGRC <- c(NA, diff(log(fit$Area1.smooth.C))/diff(fit$xDays))
  ggplot(fit) + geom_line(aes(x=xDays, y=Area1.AGR)) + 
    geom_line(aes(x=xDays, y=AGR), colour = "blue") + 
    geom_line(aes(x=xDays, y=AGRC), colour = "red")
  ggplot(fit) + geom_line(aes(x=xDays, y=Area1.RGR)) + 
    geom_line(aes(x=xDays, y=RGR), colour = "blue") + 
    geom_line(aes(x=xDays, y=RGRC), colour = "red")
  testthat::expect_true(abs(var(fit$Area1.smooth.C) - 66372.73) < 1e-02)
  testthat::expect_true(abs(var(fit$AGRC, na.rm = TRUE) - 145.0717) < 1e-03)
  testthat::expect_true(abs(var(fit$RGRC, na.rm = TRUE) - 0.000893429) < 1e-04)
  
})

cat("#### Test byIndv4Times_SplinesGRs using PS with leaf data when there are missing values\n")
test_that("leaf_growthPheno", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(dae)
  library(ggplot2)
  library(growthPheno)
  
  # A small subset of Exp 270 leaf tracking data
  data(testSpline)
  test$Days[is.na(test$xDays)] <- NA
  responses <- names(test)[5:ncol(test)]
  
  ##Test the fail options for PS
  testthat::expect_error(byIndv4Times_SplinesGRs(data = test, response = "Area", times = "Days", 
                                                 spline.type = "PS", lambda = 0.1, npspline.segments = 4, 
                                                 rates.method = "none", na.x.action = "fail"))
  testthat::expect_error(byIndv4Times_SplinesGRs(data = test, response = "Length.1", times = "Days", 
                                                 spline.type = "PS", lambda = 0.1, npspline.segments = 4, 
                                                 rates.method = "none", na.y.action = "fail"))
  
  ##Fit some splines - exclude y
  leaf.dat <- test
  resp <- responses[2]
  testthat::expect_warning(
    for (resp in responses)
      leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, 
                                          response = resp, times = "Days", 
                                          spline.type = "PS", lambda = 0.1, npspline.segments = 4, 
                                          rates.method = "none", na.y.action = "exclude"),
    regexp = "Need at least 4 distinct x values to fit a spline - all fitted values set to NA")
  testthat::expect_equal(sum(!is.na(leaf.dat$Length.1)), 85)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.1)), 82)
  testthat::expect_equal(sum(is.na(leaf.dat$sArea)), 3)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.3[13:24])), 4)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.3[49:60])), 0)  
  
  ##trim
  leaf.dat <- test
  testthat::expect_warning(
    for (resp in responses)
      leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, response = resp, times = "Days", 
                                          spline.type = "PS", lambda = 0.1, npspline.segments = 4, 
                                          rates.method = "none", na.y.action = "trim"),
    regexp = "Need at least 4 distinct x values to fit a spline - all fitted values set to NA")
  testthat::expect_equal(sum(!is.na(leaf.dat$Length.1)), 85)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.1)), 86)
  testthat::expect_lt(abs(leaf.dat$sLength.1[57] - 11.07072), 1e-05)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.3[13:24])), 4)
  
  ##Test ltrimx
  leaf.dat <- test
  for (resp in responses)
    leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, response = resp, times = "Days", 
                                        spline.type = "PS", lambda = 0.1, npspline.segments = 4, 
                                        na.y.action = "ltrim", 
                                        rates.method = "deriv", 
                                        suffices.rates = c("AGRdv", "RGRdv"))
  testthat::expect_equal(length(lapply(responses, 
                                       function(resp, data)
                                       {
                                         resp <- paste(resp,"smooth", sep = ".")
                                         testthat::expect_equal(sum(is.na(data[[resp]])), 
                                                                sum(is.na(data[[paste(resp,"AGRdv",sep = ".")]])),
                                                                sum(is.na(data[[paste(resp,"RGRdv",sep = ".")]])))
                                       }, data = leaf.dat)), 4)
  testthat::expect_equal(nrow(leaf.dat), 96)
  testthat::expect_equal(ncol(leaf.dat), 20)
  testthat::expect_equal(sum(!is.na(leaf.dat$Length.1)), 85)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.1)), 87)
  testthat::expect_lt(abs(leaf.dat$sLength.1[57] - 11.07072), 1e-05)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.3[13:24])), 7)
  leaf.dat.noC <- leaf.dat
  
  ##Test ltrimx with correctBoundaries (which is ignored)
  leaf.dat <- test
  for (resp in responses)
    leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, response = resp, times = "Days", 
                                        spline.type = "PS", lambda = 0.1, npspline.segments = 4, 
                                        rates.method = "none", na.y.action = "ltrim", 
                                        correctBoundaries = TRUE)
  testthat::expect_equal(nrow(leaf.dat), 96)
  testthat::expect_equal(ncol(leaf.dat), 12)
  testthat::expect_equal(sum(!is.na(leaf.dat$Length.1)), 85)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.1)), 87)
  testthat::expect_lt(abs(leaf.dat$sLength.1[57] - 11.07072), 1e-05)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.3[13:24])), 7)
  
  ##Test utrimx
  leaf.dat <- test
  for (resp in responses)
    leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, response = resp, times = "Days", 
                                        spline.type = "PS", lambda = 0.1, npspline.segments = 4, 
                                        rates.method = "none", na.y.action = "utrim")
  testthat::expect_equal(nrow(leaf.dat), 96)
  testthat::expect_equal(ncol(leaf.dat), 12)
  testthat::expect_equal(sum(!is.na(leaf.dat$Length.1)), 85)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.1)), 92)
  testthat::expect_lt(abs(leaf.dat$sLength.1[57] - 11.02593), 1e-05)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.3[13:24])), 9)
  
  ##Test utrimx
  leaf.dat <- test
  for (resp in responses)
    leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, response = resp, times = "Days", 
                                        spline.type = "PS", lambda = 0.1, npspline.segments = 4, 
                                        rates.method = "none", na.y.action = "allx")
  testthat::expect_equal(nrow(leaf.dat), 96)
  testthat::expect_equal(ncol(leaf.dat), 12)
  testthat::expect_equal(sum(!is.na(leaf.dat$Length.1)), 85)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.1)), 93)
  testthat::expect_lt(abs(leaf.dat$sLength.1[57] - 11.02593), 1e-05)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.3[13:24])), 12)
  testthat::expect_equal(sum(is.na(leaf.dat$sLength.3[94:96])), 3)
  
  
  ##Test omit in fitSpline - Length 3
  leaf.dat <- test
  carts <- levels(leaf.dat$Snapshot.ID.Tag)
  fit <- vector(mode = "list", length = 0)
  nrows <- vector(mode = "list", length = length(carts))
  nrows <- list(6,4,5,3,2,6,1,1)
  names(nrows) <- carts
  for (cart in carts)
  {
    fit[[cart]] <- fitSpline(subset(leaf.dat, Snapshot.ID.Tag == cart), 
                             response = "Length.3", response.smoothed = "sLength.3", 
                             x="xDays", 
                             spline.type = "PS", lambda = 0.1, npspline.segments = 4, 
                             na.x.action = "omi", na.y.action = "omit", 
                             deriv=1, suffices.deriv="AGRdv",  
                             extra.rate = c(RGRdv = "RGR"))
    testthat::expect_equal(nrows[[cart]], nrow(fit[[cart]]$predictions))
  }
  
  ##Test omit in fitSpline - Length 2 with a 0 length data.frame
  fit <- list()
  nrows <- list(11,12,12,12,9,12,0,9)
  names(nrows) <- carts
  for (cart in carts)
  {
    fit[[cart]] <- fitSpline(subset(leaf.dat, Snapshot.ID.Tag == cart), 
                             response = "Length.2", response.smoothed = "sLength.2", 
                             x="xDays", 
                             spline.type = "PS", lambda = 0.1, npspline.segments = 4, 
                             na.x.action = "omi", na.y.action = "omit", 
                             deriv=1, suffices.deriv="AGRdv",  
                             extra.rate = c(RGRdv = "RGR"))
    testthat::expect_equal(nrows[[cart]], nrow(fit[[cart]]$predictions))
  }
  testthat::expect_equal(ncol(fit[[cart]]$predictions), 4)
  
  ##Test omit in byIndv4Times_SplinesGRs - get full data.frame because of merge in byIndv4Times_SplinesGRs
  leaf.dat <- test
  for (resp in responses)
    leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, response = resp, times = "Days", 
                                        spline.type = "PS", lambda = 0.1, npspline.segments = 4, 
                                        rates.method = "none", 
                                        na.x.action = "omi", na.y.action = "omit")
  testthat::expect_equal(nrow(leaf.dat), 96)
  testthat::expect_equal(ncol(leaf.dat), 12)
  testthat::expect_equal(sum(!is.na(leaf.dat$Length.1)), 85)
  testthat::expect_equal(sum(!is.na(leaf.dat$sLength.1)), 82)
  
  leaf.dat <- test
  for (resp in responses)
    leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, response = resp, times = "Days", 
                                        spline.type = "PS", lambda = 0.1, npspline.segments = 4, 
                                        rates.method = "none", na.y.action = "omit")
  
  
  ##Test omit in fitSpline - Length 2 with a 0 length data.frame
  leaf.dat <- test
  carts <- levels(leaf.dat$Snapshot.ID.Tag)
  fit <- vector(mode = "list", length = 0)
  nrows <- list(6,4,5,3,2,6,1,1)
  names(nrows) <- carts
  for (cart in carts)
  {
    fit[[cart]] <- fitSpline(subset(leaf.dat, Snapshot.ID.Tag == cart), 
                             response = "Length.3", response.smoothed = "sLength.3", 
                             x="xDays", correctBoundaries = FALSE,
                             spline.type = "PS", lambda = 0.1, npspline.segments = 4, 
                             na.x.action = "omit", na.y.action = "omit")
    testthat::expect_equal(ncol(fit[[cart]]$predictions), 2)
    testthat::expect_equal(nrow(fit[[cart]]$predictions), nrows[[cart]])
  }  
  fitC <- vector(mode = "list", length = 0)
  nrows <- list(6,4,5,3,2,6,1,1)
  names(nrows) <- carts
  for (cart in carts)
  {
    fitC[[cart]] <- fitSpline(subset(leaf.dat, Snapshot.ID.Tag == cart), 
                              response = "Length.3", response.smoothed = "sLength.3", 
                              x="xDays",
                              spline.type = "PS", lambda = 0.1, npspline.segments = 4, 
                              na.x.action = "omit", na.y.action = "omit")
    testthat::expect_equal(ncol(fitC[[cart]]$predictions), 2)
    testthat::expect_equal(nrow(fitC[[cart]]$predictions), nrows[[cart]])
  }  
})


cat("#### Test byIndv4Times_SplinesGRs using PS with leaf data, including log-smoothing\n")
test_that("leaf_growthPheno", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(dae)
  library(ggplot2)
  library(growthPheno)
  
  # A small subset of Exp 270 leaf tracking data
  data(testSpline)
  test$Days[is.na(test$xDays)] <- NA
  
  ##Smooth splines using identity and logarithm transformations - exclude y
  leaf.dat <- test
  leaf.dat$Area.log <- log(leaf.dat$Area)
  #Investigate AGR and RGR calculations
  leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, 
                                      response="Area", response.smoothed = "Area.smooth", 
                                      times = "Days", 
                                      spline.type = "PS", lambda = 0.1, npspline.segments = 4, 
                                      smoothing.method = "log",
                                      rates.method = "deriv", 
                                      suffices.rates = c("AGRdv", "RGRdv"))
  names(leaf.dat)[(match(c("Area.smooth", "Area.smooth.AGRdv","Area.smooth.RGRdv"), 
                         names(leaf.dat)))] <- paste(c("Area.smooth", "Area.smooth.AGRdv", 
                                                       "Area.smooth.RGRdv"), "log", sep = ".")
  testthat::expect_equal(sum(is.na(leaf.dat$Area.smooth.log)), 3)
  testthat::expect_false(any(abs(leaf.dat$Area.smooth.AGRdv.log[1:3] - 
                                   c(3.009742, 3.704214, 4.556166)) > 1e-03, 
                             na.rm = TRUE))
  testthat::expect_false(any(abs(leaf.dat$Area.smooth.RGRdv.log[1:3] - 
                                   c(0.1893088, 0.1924901, 0.1950459)) > 1e-03, 
                             na.rm = TRUE))
  #Manual calculation of log smooth
  leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, 
                                      response = "Area.log", response.smoothed = "Area.log.smooth", 
                                      times = "Days", 
                                      spline.type = "PS", lambda = 0.1, npspline.segments = 4,
                                      rates.method = "none")
  leaf.dat$Area.log.smooth <- exp(leaf.dat$Area.log.smooth)
  testthat::expect_false(any(abs(leaf.dat$Area.log.smooth - leaf.dat$Area.smooth.log) > 0.1, 
                             na.rm = TRUE))
  testthat::expect_equal(sum(is.na(leaf.dat$Area.log.smooth)), 3)
  
  #identity smoothing scale calculation of AGR and RGR
  leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, 
                                      response="Area", response.smoothed = "Area.smooth", 
                                      times = "Days", 
                                      spline.type = "PS", lambda = 0.1, npspline.segments = 4, 
                                      rates.method = "deriv", 
                                      suffices.rates = c("AGRdv", "RGRdv"))
  testthat::expect_false(any(abs(leaf.dat$Area.smooth[1:3] - 
                                   c(14.89518, 18.97489, 23.57340)) > 1e-03, 
                             na.rm = TRUE))
  testthat::expect_false(any(abs(leaf.dat$Area.smooth.AGRdv[1:3] - 
                                   c(3.910074, 4.294226, 4.947678)) > 1e-03, 
                             na.rm = TRUE))
  testthat::expect_false(any(abs(leaf.dat$Area.smooth.RGRdv[1:3] - 
                                   c(0.2625059, 0.2263110, 0.2098839)) > 1e-03, 
                             na.rm = TRUE))
  
  #identity smoothing scale calculation of AGR and RGR - default calculation of npspline.segments
  leaf.dat <- byIndv4Times_SplinesGRs(data = leaf.dat, 
                                      response="Area", response.smoothed = "Area.smooth", 
                                      times = "Days", 
                                      spline.type = "PS", lambda = 0.1, 
                                      rates.method = "deriv", 
                                      suffices.rates = c("AGRdv", "RGRdv"))
  testthat::expect_false(any(abs(leaf.dat$Area.smooth[1:3] - 
                                   c(16.07077, 18.99996, 22.99637)) > 1e-03, 
                             na.rm = TRUE))
  testthat::expect_false(any(abs(leaf.dat$Area.smooth.AGRdv[1:3] - 
                                   c(2.599351, 3.367556, 4.687257)) > 1e-03, 
                             na.rm = TRUE))
  testthat::expect_false(any(abs(leaf.dat$Area.smooth.RGRdv[1:3] - 
                                   c(0.1617440, 0.1772402, 0.2038260)) > 1e-03, 
                             na.rm = TRUE))
})


cat("#### Test fitSpline using PS with leaf data when there are missing values\n")
test_that("leaf_growthPheno", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(dae)
  library(ggplot2)
  library(growthPheno)
  
  # A small subset of Exp 270 leaf tracking data
  data(testSpline)
  test$Days[is.na(test$xDays)] <- NA
  responses <- names(test)[5:ncol(test)]
  leaf.dat <- test
  
  ##Check the contents of the fit.spline object for PS
  carts <- levels(leaf.dat$Snapshot.ID.Tag)
  fit <- fitSpline(subset(leaf.dat, Snapshot.ID.Tag == carts), 
                   response = "Length.2", response.smoothed = "Length.2.smooth", 
                   x="xDays", 
                   spline.type = "PS", lambda = 0.1, npspline.segments = 4, 
                   na.x.action = "exclude", na.y.action = "allx", 
                   rates.method = "deriv", 
                   suffices.rates = c("AGRdv", "RGRdv"))
  testthat::expect_true(all(c(fit$fit.spline$lambda, fit$fit.spline$uncorrected.fit$lambda) == 0.1))
  testthat::expect_true(all(c(fit$fit.spline$npspline.segments, fit$fit.spline$uncorrected.fit$npspline.segments) == 4))
  testthat::expect_true(all((c(fit$fit.spline$df, fit$fit.spline$uncorrected.fit$effdim) -3.816562) < 1e-04))
  testthat::expect_equal(length(subset(leaf.dat, Snapshot.ID.Tag == carts)$Length.2) - 
                           length(fit$fit.spline$x) , 1)
  testthat::expect_equal(length(subset(leaf.dat, Snapshot.ID.Tag == carts)$Length.2) - 
                           length(fit$fit.spline$uncorrected.fit$xgrid) , 1)
})

cat("#### Test byIndv4Times_SplinesGRs with small example\n")
test_that("exampleData_byIndv4Times_SplinesGRs", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(growthPheno)
  
  data(exampleData)
  
  #Segmented smoothing
  t <- byIndv4Times_SplinesGRs(data = longi.dat, 
                               response = "PSA", response.smoothed = "sPSA", 
                               times = "DAP", 
                               smoothing.segments = list(c(28,34), c(35,42)), 
                               df = 5, rates.method = "none")
  testthat::expect_equal(nrow(t), 280)
  testthat::expect_equal(ncol(t), 37)
  testthat::expect_true(all(unique(t$DAP) == c(28, 30:42)))
  testthat::expect_false(all(abs(t$sPSA-longi.dat$sPSA) < 1e-04))
  
  
  testthat::expect_silent(
    plotProfiles(t, x = "DAP", response = "sPSA", 
                 facet.y = "Treatment.1", alpha = 0.75, 
                 y.title = "sPSA (kpixels)", 
                 ggplotFuncs = list(ggplot2::geom_vline(xintercept=35, linetype="longdash", 
                                               alpha = 0.5, size=0.75, colour = "blue"))))

  #smooth = direct, rates.methods = "none
  dat <- longi.dat[,-match(c("sPSA.AGR", "sPSA.RGR"), names(longi.dat))]
  testthat::expect_silent(  
    t <- byIndv4Times_SplinesGRs(data = dat, response="sPSA", 
                                 times = "DAP", 
                                 rates.method = "none", df = 4))
  testthat::expect_equal(nrow(t), 280)
  testthat::expect_equal(ncol(t), 36)

  #smooth = direct, rates.methods = "deriv", which.rates = AGR & RGR
  testthat::expect_silent(  
    all.rates <- byIndv4Times_SplinesGRs(data = dat, response="PSA", 
                                         times = "DAP", 
                                         rates.method = "deriv", which.rates = c("AGR","RGR"),
                                         df = 4))
  testthat::expect_equal(nrow(all.rates), 280)
  testthat::expect_equal(ncol(all.rates), 37)
  testthat::expect_true(any(abs(all.rates$sPSA.AGR[1:3] - 
                                  c(17.78304, 19.16926, 20.99404)) < 1e-03))
  testthat::expect_true(any(abs(all.rates$sPSA.RGR[1:3] - 
                                  c(0.3474298, 0.2186439, 0.1949627)) < 1e-05))
  
  #smooth = direct, rates.methods = "deriv", which.rates = AGR
  testthat::expect_silent(  
    t <- byIndv4Times_SplinesGRs(data = dat, response="PSA", 
                                 times = "DAP", 
                                 rates.method = "deriv", which.rates = "AGR",
                                 df = 4))
  testthat::expect_equal(nrow(t), 280)
  testthat::expect_equal(ncol(t), 36)
  testthat::expect_true(any(abs(t$sPSA.AGR - all.rates$sPSA.AGR) < 1e-05))
  
  #smooth = direct, rates.methods = "deriv", which.rates = RGR
  testthat::expect_silent(  
    t <- byIndv4Times_SplinesGRs(data = dat, response="PSA", 
                                 times = "DAP", 
                                 rates.method = "deriv", which.rates = "RGR",
                                 df = 4))
  testthat::expect_equal(nrow(t), 280)
  testthat::expect_equal(ncol(t), 36)
  testthat::expect_true(any(abs(t$sPSA.RGR - all.rates$sPSA.RGR) < 1e-05))
  
  #smooth = log, rates.methods = "none
  dat <- longi.dat[,-match(c("sPSA.AGR", "sPSA.RGR"), names(longi.dat))]
  testthat::expect_silent(  
    t <- byIndv4Times_SplinesGRs(data = dat, response="PSA", 
                                 times = "DAP", 
                                 rates.method = "none", 
                                 smoothing.method = "log", df = 4))
  testthat::expect_equal(nrow(t), 280)
  testthat::expect_equal(ncol(t), 35)
  
  #smooth = log, rates.methods = "deriv", which.rates = AGR & RGR
  dat <- longi.dat[,-match(c("sPSA.AGR", "sPSA.RGR"), names(longi.dat))]
  testthat::expect_silent(  
    all.rates <- byIndv4Times_SplinesGRs(data = dat, response="PSA", 
                                         times = "DAP", 
                                         rates.method = "deriv", which.rates = c("AGR","RGR"),
                                         smoothing.method = "log", df = 4))
  testthat::expect_equal(nrow(all.rates), 280)
  testthat::expect_equal(ncol(all.rates), 37)
  testthat::expect_true(any(abs(all.rates$sPSA.AGR[1:3] - 
                                  c(11.65445, 16.87781, 19.97875)) < 1e-03))
  testthat::expect_true(any(abs(all.rates$sPSA.RGR[1:3] - 
                                  c(0.1986342, 0.1939544, 0.1895384)) < 1e-05))
  
  #smooth = log, rates.methods = "deriv", which.rates = AGR
  testthat::expect_silent(  
    t <- byIndv4Times_SplinesGRs(data = dat, response="PSA", 
                                 times = "DAP", 
                                 rates.method = "deriv", which.rates = "AGR",
                                 smoothing.method = "log", df = 4))
  testthat::expect_equal(nrow(t), 280)
  testthat::expect_equal(ncol(t), 36)
  testthat::expect_true(any(abs(t$sPSA.AGR - all.rates$sPSA.AGR) < 1e-05))
  
  #smooth = log, rates.methods = "deriv", which.rates = RGR
  testthat::expect_silent(  
    t <- byIndv4Times_SplinesGRs(data = dat, response="PSA", 
                                 times = "DAP", 
                                 rates.method = "deriv", which.rates = "RGR",
                                 smoothing.method = "log", df = 4))
  testthat::expect_equal(nrow(t), 280)
  testthat::expect_equal(ncol(t), 36)
  testthat::expect_true(any(abs(t$sPSA.RGR - all.rates$sPSA.RGR) < 1e-05))
  
  #smooth = direct, rates.methods = "deriv", which.rates = AGR & RGR, P-splines
  testthat::expect_silent(  
    all.rates <- byIndv4Times_SplinesGRs(data = dat, response="PSA", 
                                         times = "DAP", 
                                         rates.method = "deriv", which.rates = c("AGR","RGR"),
                                         spline.type = "PS", lambda = 1))
  testthat::expect_equal(nrow(all.rates), 280)
  testthat::expect_equal(ncol(all.rates), 37)
  testthat::expect_true(any(abs(all.rates$sPSA.AGR[1:3] - 
                                  c(15.77705, 18.01763, 20.59344)) < 1e-03))
  testthat::expect_true(any(abs(all.rates$sPSA.RGR[1:3] - 
                                  c(0.2886963, 0.2051182, 0.1923348)) < 1e-05))
  
  #smooth = direct, rates.methods = "deriv", which.rates = AGR, P-splines
  testthat::expect_silent(  
    t <- byIndv4Times_SplinesGRs(data = dat, response="PSA", 
                                 times = "DAP", 
                                 rates.method = "deriv", which.rates = "AGR",
                                 spline.type = "PS", lambda = 1))
  testthat::expect_equal(nrow(t), 280)
  testthat::expect_equal(ncol(t), 36)
  testthat::expect_true(any(abs(t$sPSA.AGR - all.rates$sPSA.AGR) < 1e-05))
  
  #smooth = direct, rates.methods = "deriv", which.rates = RGR
  testthat::expect_silent(  
    t <- byIndv4Times_SplinesGRs(data = dat, response="PSA", 
                                 times = "DAP", 
                                 rates.method = "deriv", which.rates = "RGR",
                                 spline.type = "PS", lambda = 1))
  testthat::expect_equal(nrow(t), 280)
  testthat::expect_equal(ncol(t), 36)
  testthat::expect_true(any(abs(t$sPSA.RGR - all.rates$sPSA.RGR) < 1e-05))
  
  #Fitting an extra.derivs
  testthat::expect_silent(  
    t <- byIndv4Times_SplinesGRs(data = longi.dat, response="PSA", 
                                   times = "DAP", 
                                   rates.method = "deriv", df = 4, 
                                   extra.derivs = 2, suffices.extra.derivs = "Accel"))
  testthat::expect_equal(nrow(t), 280)
  testthat::expect_equal(ncol(t), 38)
  testthat::expect_true(all(c("sPSA","sPSA.AGR","sPSA.RGR",
                              "sPSA.Accel") %in% names(t)))
  testthat::expect_true(any(abs(t$sPSA.AGR[1:3] - 
                                   c(17.78304, 19.16926, 20.99404)) < 1e-03))
  testthat::expect_true(any(abs(t$sPSA.RGR[1:3] - 
                                  c(0.3474298, 0.2186439, 0.1949627)) < 1e-05))
  testthat::expect_true(any(abs(t$sPSA.Accel[1:3] - 
                                  c(-0.002044201,  1.388259633,  2.261297448)) < 1e-05))
  
  #Calculate 2nd deriv by differencing the AGR and compare with deriv =2 using r
  t <- byIndv4Times_GRsDiff(data = t, response="sPSA.AGR", times = "DAP", 
                            which.rates = "AGR", suffices.rates = "Acc", 
                            ntimes2span = 3)
  testthat::expect_true(abs(cor(na.omit(t[c("sPSA.Accel", "sPSA.AGR.Acc")]))[1,2] - 
                              0.9965408) < 1e-05)
  
  #Calculate 2nd deriv using PS to compare with NCSS
  #Fitting an extra.derivs
  testthat::expect_silent(  
    PS <- byIndv4Times_SplinesGRs(data = longi.dat, response="PSA", 
                                 times = "DAP", 
                                 spline.type = "PS", lambda = 1, 
                                 rates.method = "deriv",
                                 extra.derivs = 2, suffices.extra.derivs = "Accel"))
  testthat::expect_equal(nrow(PS), 280)
  testthat::expect_equal(ncol(PS), 38)
  testthat::expect_true(all(c("sPSA","sPSA.AGR","sPSA.RGR",
                              "sPSA.Accel") %in% names(t)))
  testthat::expect_true(any(abs(PS$sPSA.AGR[1:3] - 
                                  c(15.77705, 18.01763, 20.59344)) < 1e-03))
  testthat::expect_true(any(abs(PS$sPSA.RGR[1:3] - 
                                  c(0.2886963, 0.2051182, 0.1923348)) < 1e-05))
  testthat::expect_true(any(abs(PS$sPSA.Accel[1:3] - 
                                  c(0.2378153, 2.0798768, 2.8410223)) < 1e-05))
  testthat::expect_true(abs(
    cor(na.omit(data.frame(NCSS = t$sPSA.Accel, PS = PS$sPSA.Accel)))[1,2] - 0.9721318) < 1e-05)
})
