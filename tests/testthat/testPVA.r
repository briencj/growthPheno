cat("#### Test PVA\n")
test_that("PVA", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(growthPheno)

  data(exampleData)
  responses <- c("Area","Area.SV","Area.TV", "Image.Biomass", "Max.Height","Centre.Mass",
                 "Density", "Compactness.TV", "Compactness.SV")
  selected.responses.95 <- PVA(responses = responses, 
                               data = longi.dat, 
                               p.variance = 0.95, plot = FALSE)
  testthat::expect_equal(nrow(selected.responses.95), 4)
  testthat::expect_lt(abs(selected.responses.95[1, "Cumulative.Propn"] - 0.6715883), 1e-04)
  
  #test when name in include is not in responses
  testthat::expect_error(PVA(responses = responses, 
                             data = longi.dat, 
                             p.variance = 0.999,
                             include = c("Area.SF"), plot = FALSE))
  
  #Test when response not in data
  testthat::expect_error(PVA(responses = c(responses, "Area.SF"), 
                             data = longi.dat, 
                             p.variance = 0.999,
                             include = c("Area.TV"), plot = FALSE))
  
  #Test that when only the variable in include is seleted
  testthat::expect_equal(nrow(PVA(responses = responses, data = longi.dat, p.variance = 0.60,
                                  include = "Area", plot = FALSE)), 1)
})