library(datasets)
library(scPipe)
library(scuttle)
library(iDA)
library(dplyr)
data(iris)
data("sc_sample_data")
data("sc_sample_qc")

sce <- SingleCellExperiment(assays = list(counts = as.matrix(sc_sample_data)))
logcounts(sce) <- normalizeCounts(sce,  size.factors = sizeFactors(sce))

#test VariableGenesGeneric
var.genes <- VariableGenesGeneric(logcounts(sce), 
                        mean.low.cutoff = 0.1,
                        mean.high.cutoff = 8,
                        dispersion.cutoff = 1)
test_that("dispersion selection option returns expected values", {
  expect_equal(var.genes[["var.features"]][1:3], c("ENSMUSG00000053044", "ENSMUSG00000071561", "ENSMUSG00000037944"))
  expect_equal(length(x = var.genes$var.features), 141)
  expect_equal(var.genes[["dispersions"]][["scaled.dispersions"]][1:3], c(2.243072656, 0.126951097, 2.978339717), tolerance = 1e-6)
  expect_equal(var.genes[["dispersions"]][["ExpMeans"]][1:3], c(0.93090698, 0.08538151, 0.10246037), tolerance = 1e-6)
})

# test WCS
test_that("calculating within class scatter matrix works correctly", {
  iris_input <- select(iris, -c("Species"))
  sp_dfs <- split(iris_input, f = iris$Species)
  wcs_mat <- WCS(sp_dfs, diag = TRUE)
  expect_equal(nrow(wcs_mat), ncol(iris_input))
  expect_equal(ncol(wcs_mat), ncol(iris_input))
 })

# test BCS
test_that("calculating between class scatter matrix works correctly", {
  iris_input <- select(iris, -c("Species"))
  sp_dfs <- split(iris_input, f = iris$Species)
  bcs_mat <- BCS(sp_dfs)
  expect_equal(nrow(bcs_mat), ncol(iris_input))
  expect_equal(ncol(bcs_mat), ncol(iris_input))
})

# test decomposesvd
test_that("decomposeSVD works correctly", {
    iris_input <- select(iris, -c("Species"))
    sp_dfs <- split(iris_input, f = iris$Species)
    bcs_mat <- BCS(sp_dfs)
    wcs_mat <- WCS(sp_dfs)
    svd <- decomposeSVD(wcs_mat, bcs_mat, nu = 4)
    expect_equal(ncol(svd$eigenvecs), 4)
    expect_equal(nrow(svd$eigenvecs), 4)
    expect_equal(length(svd$eigenvalues), 4)
})

# test getSNN




# test getLouvain









