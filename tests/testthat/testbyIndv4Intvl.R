

cat("#### Test byIndv4Intvl_WaterUse with small example\n")
test_that("exampleData_waterUse", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(growthPheno)
  
  data(exampleData)
  
  

  test <- byIndv4Times_SplinesGRs(data = longi.dat, response="PSA", times = "DAP", 
                                  df = 5, rates.method = "none")
  testthat::expect_equal(nrow(test), 280)
  testthat::expect_equal(ncol(test), 37)
  testthat::expect_true(all(unique(test$Day) == c(28, 30:42)))
  testthat::expect_false(all(abs(test$sPSA-longi.dat$sPSA) < 1e-04))
  
  test <- merge(test, 
                byIndv4Intvl_WaterUse(data = test, 
                                      water.use = "WU", 
                                      responses = "PSA", 
                                      trait.types = c("WUR","WUI"), 
                                      suffix.rate = ".Rate", 
                                      suffix.index = ".Index",
                                      start.time = 31, end.time = 35, 
                                      suffix.interval = "31to35"),
                by = "Snapshot.ID.Tag")
  testthat::expect_equal(nrow(test), 280)
  testthat::expect_equal(ncol(test), 39)
  testthat::expect_true(all(c("WU.Rate.31to35", "PSA.WU.Index.31to35") %in% 
                              names(test)))
  
  WU.WUI_31_35 <- byIndv4Intvl_WaterUse(data = longi.dat, 
                                        water.use = "WU", 
                                        responses = "PSA", 
                                        trait.types = "WU", 
                                        start.time = 31, end.time = 35, 
                                        suffix.interval = "31to35")
  testthat::expect_equal(nrow(WU.WUI_31_35), 20)
  testthat::expect_equal(ncol(WU.WUI_31_35), 2)
  testthat::expect_true(all(names(WU.WUI_31_35) == c("Snapshot.ID.Tag", "WU.31to35")))
  
  #Test values for all
  Tag1 <- subset(longi.dat, subset = Snapshot.ID.Tag == "045451-C" & xDAP >= 31 & xDAP <= 35)
  WU.WUI_31_35 <- byIndv4Intvl_WaterUse(data = longi.dat, 
                                        water.use = "WU", 
                                        responses = "PSA", 
                                        trait.types = "all", 
                                        start.time = 31, end.time = 35, 
                                        suffix.interval = "31to35")
  testthat::expect_equal(nrow(WU.WUI_31_35), 20)
  testthat::expect_equal(ncol(WU.WUI_31_35), 5)
  testthat::expect_true(all(names(WU.WUI_31_35) == c("Snapshot.ID.Tag", "WU.31to35", 
                                                     "WUR.31to35", "PSA.AGR.31to35", 
                                                     "PSA.WUI.31to35")))
  testthat::expect_equal(WU.WUI_31_35$WU.31to35[1], sum(Tag1$WU)-Tag1$WU[1])
  testthat::expect_equal(WU.WUI_31_35$WUR.31to35[1], (sum(Tag1$WU)-Tag1$WU[1])/(35-31))
  testthat::expect_equal(WU.WUI_31_35$PSA.AGR.31to35[1], Tag1$PSA[5]-Tag1$PSA[1])
  testthat::expect_equal(WU.WUI_31_35$PSA.WUI.31to35[1], 
                         (Tag1$PSA[5]-Tag1$PSA[1])/(sum(Tag1$WU)-Tag1$WU[1]))

  #Test NULL suffix.interval  
  WU.WUI_31_35 <- byIndv4Intvl_WaterUse(data = longi.dat, 
                                        water.use = "WU", 
                                        responses = "PSA", 
                                        trait.types = "all", 
                                        start.time = 31, end.time = 35, 
                                        suffix.interval = NULL)
  testthat::expect_equal(nrow(WU.WUI_31_35), 20)
  testthat::expect_equal(ncol(WU.WUI_31_35), 5)
  testthat::expect_true(all(names(WU.WUI_31_35) == c("Snapshot.ID.Tag", "WU", "WUR",
                                                     "PSA.AGR", "PSA.WUI")))
  #Test for WU and all
  names(test)[match("WU", names(test))] <- "Water.Use"
  WU.WUI_31_35 <- byIndv4Intvl_WaterUse(data = test, 
                                        water.use = "Water.Use", 
                                        responses = "PSA", 
                                        trait.types = "all", 
                                        suffix.rate = ".Rate", 
                                        suffix.index = ".Index",
                                        start.time = 31, end.time = 35, 
                                        suffix.interval = "31to35")
  testthat::expect_equal(nrow(WU.WUI_31_35), 20)
  testthat::expect_equal(ncol(WU.WUI_31_35), 5)
  testthat::expect_true(all(names(WU.WUI_31_35) == c("Snapshot.ID.Tag", "Water.Use.31to35", 
                                                     "Water.Use.Rate.31to35", "PSA.AGR.31to35", 
                                                     "PSA.Water.Use.Index.31to35")))
  
  #Test for WU and default
  WU.WUI_31_35 <- byIndv4Intvl_WaterUse(data = test, 
                                        water.use = "Water.Use", 
                                        responses = "PSA", 
                                        suffix.rate = ".Rate", 
                                        suffix.index = ".Index",
                                        start.time = 31, end.time = 35, 
                                        suffix.interval = "31to35")
  testthat::expect_equal(nrow(WU.WUI_31_35), 20)
  testthat::expect_equal(ncol(WU.WUI_31_35), 4)
  testthat::expect_true(all(names(WU.WUI_31_35) == c("Snapshot.ID.Tag", "Water.Use.31to35", 
                                                     "Water.Use.Rate.31to35", 
                                                     "PSA.Water.Use.Index.31to35")))
  #Test for sWU and all
  names(test)[match("sWU", names(test))] <- "sWater.Use"
  WU.WUI_31_35 <- byIndv4Intvl_WaterUse(data = test, 
                                        water.use = "sWater.Use", 
                                        responses = "sPSA", 
                                        trait.types = "all", 
                                        suffix.rate = ".Rate", 
                                        suffix.index = ".Index",
                                        start.time = 31, end.time = 35, 
                                        suffix.interval = "31to35")
  testthat::expect_equal(nrow(WU.WUI_31_35), 20)
  testthat::expect_equal(ncol(WU.WUI_31_35), 5)
  testthat::expect_true(all(names(WU.WUI_31_35) == c("Snapshot.ID.Tag", "sWater.Use.31to35", 
                                                     "sWater.Use.Rate.31to35", "sPSA.AGR.31to35", 
                                                     "sPSA.sWater.Use.Index.31to35")))
  
  #Test for AGR only
  testthat::expect_warning(
    WU.WUI_31_35 <- byIndv4Intvl_WaterUse(data = longi.dat, 
                                          water.use = "WU", 
                                          responses = "PSA", 
                                          trait.types = "AGR", 
                                          start.time = 31, end.time = 35, 
                                          suffix.interval = "31to35"),
    regexp = "no water use traits have specified - only AGR")
  testthat::expect_equal(nrow(WU.WUI_31_35), 20)
  testthat::expect_equal(ncol(WU.WUI_31_35), 2)
  testthat::expect_true(all(names(WU.WUI_31_35) == c("Snapshot.ID.Tag", "PSA.AGR.31to35")))

  #Test for WUI only with multiple responses
  WU.WUI_31_35 <- byIndv4Intvl_WaterUse(data = longi.dat, 
                                        water.use = "WU", 
                                        responses = c("PSA", "sPSA"),
                                        trait.types = "WUI", 
                                        start.time = 31, end.time = 35, 
                                        suffix.interval = "31to35")
  testthat::expect_equal(nrow(WU.WUI_31_35), 20)
  testthat::expect_equal(ncol(WU.WUI_31_35), 3)
  testthat::expect_true(all(names(WU.WUI_31_35) == c("Snapshot.ID.Tag", "PSA.WUI.31to35", 
                                                     "sPSA.WUI.31to35")))
  
  #Test for NULL responses and WUR only
  WU.WUI_31_35 <- byIndv4Intvl_WaterUse(data = longi.dat, 
                                        water.use = "WU", 
                                        trait.types = "WUR", 
                                        start.time = 31, end.time = 35, 
                                        suffix.interval = "31to35")
  testthat::expect_equal(nrow(WU.WUI_31_35), 20)
  testthat::expect_equal(ncol(WU.WUI_31_35), 2)
  testthat::expect_true(all(names(WU.WUI_31_35) == c("Snapshot.ID.Tag", "WUR.31to35")))
})
