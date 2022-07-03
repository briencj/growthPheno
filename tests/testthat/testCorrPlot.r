cat("#### Test corrPlot\n")
test_that("corrPlot", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(growthPheno)
  library(ggplot2)
  library(GGally)
  
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
  testthat::expect_invisible(plt <- plotCorrmatrix(longi.dat, responses, 
                                                   pairs.sets=list(c(1:4),c(5:7))))
  testthat::expect_true("ggplot" %in% class(plt))
})