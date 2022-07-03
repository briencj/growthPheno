cat("#### Test splitContGRdiff with small example\n")
test_that("exampleData_splitContGRdiff", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(growthPheno)
  
  data(exampleData)
  
  #ntimes2span == 2, avail.times.diffs = TRUE
  t <- splitContGRdiff(longi.dat, response="PSA", 
                       times.factor = "DAP", 
                       which.rates=c("AGR", "RGR"), avail.times.diffs = TRUE) 
  testthat::expect_equal(nrow(t), 280)
  testthat::expect_equal(ncol(t), 37)
  testthat::expect_true(all(unique(t$DAP) == c(28, 30:42)))
  testthat::expect_true(all(is.na(t[t$DAP == "28", "PSA.AGR"])))
  testthat::expect_true(all(is.na(t[t$DAP == "28", "PSA.RGR"])))
  testthat::expect_true(all(t$PSA.AGR[2:4] - (t$PSA[2:4] - t$PSA[1:3])/c(2,1,1) < 1e-05))
  testthat::expect_true(all(t$PSA.RGR[2:4] - (log(t$PSA[2:4]) - log(t$PSA[1:3]))/c(2,1,1) < 1e-05))
  testthat::expect_true(all(abs(t$sPSA-longi.dat$sPSA) < 1e-04))

  #ntimes2span == 2, avail.times.diffs = FALSE
  t <- splitContGRdiff(longi.dat, response="PSA", times.factor = "DAP", 
                       which.rates=c("AGR", "RGR")) 
  testthat::expect_equal(nrow(t), 280)
  testthat::expect_equal(ncol(t), 37)
  testthat::expect_true(all(unique(t$Day) == c(28, 30:42)))
  testthat::expect_true(all(is.na(t[t$DAP == "28", "PSA.AGR"])))
  testthat::expect_true(all(is.na(t[t$DAP == "28", "PSA.RGR"])))
  testthat::expect_true(all(t$PSA.AGR[2:4] - (t$PSA[2:4] - t$PSA[1:3])/c(2,1,1) < 1e-05))
  testthat::expect_true(all(t$PSA.RGR[2:4] - (log(t$PSA[2:4]) - log(t$PSA[1:3]))/c(2,1,1) < 1e-05))
  testthat::expect_true(all(abs(t$sPSA-longi.dat$sPSA) < 1e-04))

  #ntimes2span = 3  
  t <- splitContGRdiff(longi.dat, response="PSA", times.factor = "DAP", 
                       which.rates=c("AGR", "RGR"), 
                       ntimes2span = 3, avail.times.diffs = FALSE) 
  testthat::expect_equal(nrow(t), 280)
  testthat::expect_equal(ncol(t), 37)
  testthat::expect_true(all(unique(t$Day) == c(28, 30:42)))
  testthat::expect_true(all(is.na(t[t$DAP %in% c("28","42"), "PSA.AGR"])))
  testthat::expect_true(all(is.na(t[t$DAP %in% c("28","42"), "PSA.RGR"])))
  testthat::expect_true(all(t$DAP.diffs[2:4] == c(3,2,2)))
  testthat::expect_true(all(t$PSA.AGR[2:4] - (t$PSA[3:5] - t$PSA[1:3])/c(3,2,2) < 1e-05))
  testthat::expect_true(all(t$PSA.RGR[2:4] - (log(t$PSA[3:5]) - log(t$PSA[1:3]))/c(3,2,2) < 1e-05))
  testthat::expect_true(all(abs(t$sPSA-longi.dat$sPSA) < 1e-04))
  
  #ntimes2span = 4
  t <- splitContGRdiff(longi.dat, response="PSA", times.factor = "DAP", 
                       which.rates=c("AGR", "RGR"), 
                       ntimes2span = 4, avail.times.diffs = FALSE) 
  testthat::expect_equal(nrow(t), 280)
  testthat::expect_equal(ncol(t), 37)
  testthat::expect_true(all(unique(t$Day) == c(28, 30:42)))
  testthat::expect_true(all(is.na(t[t$DAP %in% c("28","29","42"), "PSA.AGR"])))
  testthat::expect_true(all(is.na(t[t$DAP %in% c("28","29","42"), "PSA.RGR"])))
  testthat::expect_true(all(t$DAP.diffs[3:5] == c(4,3,3)))
  testthat::expect_true(all(t$PSA.AGR[3:5] - (t$PSA[4:6] - t$PSA[1:3])/c(4,3,3) < 1e-05))
  testthat::expect_true(all(t$PSA.RGR[3:5] - (log(t$PSA[4:6]) - log(t$PSA[1:3]))/c(4,3,3) < 1e-05))
  testthat::expect_true(all(abs(t$sPSA-longi.dat$sPSA) < 1e-04))
  
  #ntimes2span = 5
  t <- splitContGRdiff(longi.dat, response="PSA", times.factor = "DAP", 
                       which.rates=c("AGR", "RGR"), 
                       ntimes2span = 5, avail.times.diffs = FALSE) 
  testthat::expect_equal(nrow(t), 280)
  testthat::expect_equal(ncol(t), 37)
  testthat::expect_true(all(unique(t$Day) == c(28, 30:42)))
  testthat::expect_true(all(is.na(t[t$DAP %in% c("28","29","41","42"), "PSA.AGR"])))
  testthat::expect_true(all(is.na(t[t$DAP %in% c("28","29","41","42"), "PSA.RGR"])))
  testthat::expect_true(all(t$DAP.diffs[3:5] == c(5,4,4)))
  testthat::expect_true(all(t$PSA.AGR[3:5] - (t$PSA[5:7] - t$PSA[1:3])/c(5,4,4) < 1e-05))
  testthat::expect_true(all(t$PSA.RGR[3:5] - (log(t$PSA[5:7]) - log(t$PSA[1:3]))/c(5,4,4) < 1e-05))
  testthat::expect_true(all(abs(t$sPSA-longi.dat$sPSA) < 1e-04))
  
  #  smoothing.segments = list(c(28,34), c(35,42)), df = 5)
  #ntimes2span = 3, DAP = 30:34
  t <- splitContGRdiff(subset(longi.dat, DAP %in% as.character((30:34))), 
                       response="PSA", times.factor = "DAP", 
                       which.rates=c("AGR", "RGR"), 
                       ntimes2span = 3, avail.times.diffs = FALSE) 
  testthat::expect_equal(nrow(t), 100)
  testthat::expect_equal(ncol(t), 37)
  testthat::expect_true(all(unique(t$Day) == c(30:34)))
  testthat::expect_true(all(is.na(t[t$DAP %in% c("30","34"), "PSA.AGR"])))
  testthat::expect_true(all(is.na(t[t$DAP %in% c("30","34"), "PSA.RGR"])))
  testthat::expect_true(all(t$DAP.diffs[2:4] == 2))
  testthat::expect_true(all(t$PSA.AGR[2:4] - (t$PSA[3:5] - t$PSA[1:3])/2 < 1e-05))
  testthat::expect_true(all(t$PSA.RGR[2:4] - (log(t$PSA[3:5]) - log(t$PSA[1:3]))/2 < 1e-05))
})
