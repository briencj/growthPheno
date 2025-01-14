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
  
  #Both types of plots and printing
  testthat::expect_invisible(plt <- plotCorrmatrix(longi.dat, responses, 
                                                   pairs.sets=list(c(1:4),c(5:7))))
  testthat::expect_true("list" %in% class(plt))
  testthat::expect_equal(names(plt), c("heatmap", "matrixplots"))
  testthat::expect_true("ggplot" %in% class(plt$heatmap))
  testthat::expect_true("list" %in% class(plt$matrixplots))
  testthat::expect_equal(length(plt$matrixplots), 2)
  testthat::expect_true("ggmatrix" %in% class(plt$matrixplots[[1]]))

  #Both types of plots with no printing
  testthat::expect_invisible(plt <- plotCorrmatrix(longi.dat, responses, 
                                                   printPlot = FALSE,
                                                   pairs.sets=list(c(1:4),c(5:7))))
  testthat::expect_true("list" %in% class(plt))
  testthat::expect_equal(names(plt), c("heatmap", "matrixplots"))
  testthat::expect_true("ggplot" %in% class(plt$heatmap))
  testthat::expect_true("list" %in% class(plt$matrixplots))
  testthat::expect_equal(length(plt$matrixplots), 2)
  testthat::expect_true("ggmatrix" %in% class(plt$matrixplots[[1]]))
  
  #Heatmap only with printing
  testthat::expect_invisible(plt <- plotCorrmatrix(longi.dat, responses, 
                                                   which.plots = "heatmap",
                                                   pairs.sets=list(c(1:4),c(5:7))))
  testthat::expect_true("list" %in% class(plt))
  testthat::expect_equal(names(plt), c("heatmap", "matrixplots"))
  testthat::expect_true("ggplot" %in% class(plt$heatmap))
  testthat::expect_true(is.null(plt$matrixplots))
  

  #Matrixplots only
  testthat::expect_invisible(plt <- plotCorrmatrix(longi.dat, responses, 
                                                   which.plots = "matrixplot",
                                                   pairs.sets=list(c(1:4),c(5:7))))
  testthat::expect_true("list" %in% class(plt))
  testthat::expect_equal(names(plt), c("heatmap", "matrixplots"))
  testthat::expect_true(is.null(plt$heatmap))
  testthat::expect_equal(length(plt$matrixplots), 2)
  testthat::expect_true("ggmatrix" %in% class(plt$matrixplots[[1]]))
  
  #Heatmap only with printing
  testthat::expect_invisible(plt <- plotCorrmatrix(longi.dat, responses, 
                                                   which.plots = "heatmap", 
                                                   printPlot = FALSE,
                                                   pairs.sets=list(c(1:4),c(5:7))))
  testthat::expect_true("list" %in% class(plt))
  testthat::expect_equal(names(plt), c("heatmap", "matrixplots"))
  testthat::expect_true("ggplot" %in% class(plt$heatmap))
  testthat::expect_true(is.null(plt$matrixplots))
})
