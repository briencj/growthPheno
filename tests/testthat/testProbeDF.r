
cat("#### Test probeDF with Judith Atieno 0278\n")
test_that("chickpea_growthPheno", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(growthPheno)
  library(ggplot2)
  
  data(dat1)
  
  #'## Values of df for which to obtain plots
  df <- c(4,7)
  
  
  #'## Obtain plots
  testthat::expect_warning(t <- probeDF(dat1, df=df, xname="TimeAfterPlanting", 
                                       response="ShootArea1000",
                                       which.plots=c("bothseparately"), 
                                       ggplotFuncs=list(scale_x_continuous(breaks=
                                                                             seq(30, 50, by=2)))))
  testthat::expect_equal(nrow(t), 8480)
  testthat::expect_equal(ncol(t), 15)
  
  testthat::expect_warning(t <- probeDF(dat1, df=df, xname="TimeAfterPlanting", 
                                        response="ShootArea1000",which.traits=c("AGR")))
  testthat::expect_equal(nrow(t), 8480)
  testthat::expect_equal(ncol(t), 13)
  
  testthat::expect_warning(t <- probeDF(dat1, df=df, xname="TimeAfterPlanting", 
                                        response="ShootArea1000", which.traits=c("response"), 
                                        which.plots="compare", get.rates=FALSE))
  testthat::expect_equal(nrow(t), 8480)
  testthat::expect_true(all(c("Snapshot.ID.Tag", "Days", "TimeAfterPlanting", "ShootArea1000", 
                              "Treatment.1", "Smarthouse", "ShootArea1000.smooth.4", 
                              "ShootArea1000.smooth.7") %in% names(t)))
  
  #'## Include deviations boxplots
  testthat::expect_warning(t <- probeDF(dat1, df=df, xname="TimeAfterPlanting", 
                                        response="ShootArea1000",
                                        x.title = "DAP", 
                                        which.traits=c("response"), get.rates=FALSE,
                                        which.plots="compare", deviations.boxplots = "absolute"))
  testthat::expect_equal(nrow(t), 8480)
  testthat::expect_true(all(c("Snapshot.ID.Tag", "Days", "TimeAfterPlanting", "ShootArea1000", 
                              "Treatment.1", "Smarthouse", "ShootArea1000.smooth.4", 
                              "ShootArea1000.smooth.7") %in% names(t)))
  
  testthat::expect_warning(t <- probeDF(dat1, df=df, xname="TimeAfterPlanting", 
                                        response="ShootArea1000",
                                        x.title = "DAP", 
                                        which.traits=c("response", "AGR"), 
                                        which.plots="compare", 
                                        deviations.boxplots = c("absolute", "relative")))
  testthat::expect_equal(nrow(t), 8480)
  testthat::expect_true(all(c("Snapshot.ID.Tag", "Days", "TimeAfterPlanting", "ShootArea1000", 
                              "Treatment.1", "Smarthouse", "ShootArea1000.smooth.4", 
                              "ShootArea1000.smooth.7") %in% names(t)))
  testthat::expect_lt(max(abs(t["ShootArea1000"] - t["ShootArea1000.smooth.4"]), na.rm = TRUE) - 
                        5.563407, 1e-07)
  testthat::expect_lt(max(abs(t["ShootArea1000"] - t["ShootArea1000.smooth.4"])/
                      t["ShootArea1000.smooth.4"]) - 370.7951, 1e-04)
  
})

cat("#### Test probeDF with small example\n")
test_that("exampleData_growthPheno", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(growthPheno)
  library(ggplot2)
  
  data(exampleData)
  
  vline <- list(ggplot2::geom_vline(xintercept=20, linetype="longdash", size=1),
                ggplot2::scale_x_continuous(breaks=seq(12, 36, by=2)))
  testthat::expect_warning(tmp <- probeDF(data = longi.dat, response = "Area", 
                                          df = c(4,7), x="xDays+24.16666667", 
                 ggplotFuncs=vline))
  testthat::expect_equal(nrow(tmp), 280)
  testthat::expect_equal(ncol(tmp), 13)
  
})


