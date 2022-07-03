cat("#### Test PVA\n")
test_that("PVA", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(growthPheno)

  data(exampleData)
  longi.dat <- prepImageData(data=raw.dat, smarthouse.lev=1)
  longi.dat <- within(longi.dat, 
                      {
                        Max.Height <- pmax(Max.Dist.Above.Horizon.Line.SV1,  
                                           Max.Dist.Above.Horizon.Line.SV2)
                        Density <- PSA/Max.Height
                        PSA.SV = (PSA.SV1 + PSA.SV2) / 2
                        Image.Biomass = PSA.SV * (PSA.TV^0.5)
                        Centre.Mass <- (Center.Of.Mass.Y.SV1 + Center.Of.Mass.Y.SV2) / 2
                        Compactness.SV = (Compactness.SV1 + Compactness.SV2) / 2
                      })
  responses <- c("PSA","PSA.SV","PSA.TV", "Image.Biomass", "Max.Height","Centre.Mass",
                 "Density", "Compactness.TV", "Compactness.SV")
  selected.responses.95 <- PVA(obj = longi.dat, 
                               responses = responses, 
                               p.variance = 0.95, plot = FALSE)
  testthat::expect_equal(nrow(selected.responses.95), 4)
  testthat::expect_lt(abs(selected.responses.95[1, "Cumulative.Propn"] - 0.6834477), 1e-04)
  
  #test when name in include is not in responses
  testthat::expect_error(PVA(obj = longi.dat, 
                             responses = responses, 
                             p.variance = 0.999,
                             include = c("PSA.SF"), plot = FALSE))
  
  #Test when response not in data
  testthat::expect_error(PVA(obj = longi.dat, 
                             responses = c(responses, "Area.SF"), 
                             p.variance = 0.999,
                             include = c("Area.TV"), plot = FALSE))
  
  #Test that when only the variable in include is selected
  testthat::expect_equal(nrow(PVA(obj = longi.dat, responses = responses, p.variance = 0.60,
                                  include = "PSA", plot = FALSE)), 1)
})