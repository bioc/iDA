library(datasets)
library(scPipe)
library(scuttle)
library(iDA)
library(dplyr)
data(iris)


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









