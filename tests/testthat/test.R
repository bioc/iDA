library(datasets)
library(scPipe)
library(scuttle)
library(iDA)
data(iris)
data("sc_sample_data")
data("sc_sample_qc")

sce <- SingleCellExperiment(assays = list(counts = as.matrix(sc_sample_data)))
logcounts(sce) <- normalizeCounts(sce,  size.factors = sizeFactors(sce))

#test VariableGenes
var.genes <- VariableGenes(logcounts(sce), 
                        mean.low.cutoff = 0.1,
                        mean.high.cutoff = 8,
                        dispersion.cutoff = 1)
test_that("dispersion selection option returns expected values", {
  expect_equal(var.genes[["var.features"]][1:3], c("ENSMUSG00000053044", "ENSMUSG00000071561", "ENSMUSG00000037944"))
  expect_equal(length(x = var.genes$var.features), 141)
  expect_equal(var.genes[["dispersions"]][["scaled.dispersions"]][1:3], c(2.243072656, 0.126951097, 2.978339717), tolerance = 1e-6)
  expect_equal(var.genes[["dispersions"]][["ExpMeans"]][1:3], c(0.93090698, 0.08538151, 0.10246037), tolerance = 1e-6)
  #expect_true(all(var.genes[["var.features"]] == rownames(var.genes[["dispersions"]][scaled.dispersions > 1 & ExpMeans < 8 & ExpMeans > 0.1, ])))
})


#test split_clusters
test_that("split_clusters has correct dataframes in list ", {
  expect_equal(length(split_clusters(data = iris, clusterIDcol = iris$Species)), length(unique(iris$Species)))
  expect_equal(dim(split_clusters(data = iris, clusterIDcol = iris$Species)[[1]]), dim(iris[which(iris$Species == "setosa"),]))
  expect_equal(dim(split_clusters(data = iris, clusterIDcol = iris$Species)[[2]]), dim(iris[which(iris$Species == "versicolor"),]))
  expect_equal(dim(split_clusters(data = iris, clusterIDcol = iris$Species)[[3]]), dim(iris[which(iris$Species == "virginica"),]))
  expect_equal(dim(split_clusters(data = iris, clusterIDcol = iris$Species)[[1]])[1] + dim(split_clusters(data = iris, clusterIDcol = iris$Species)[[2]])[1] + dim(split_clusters(data = iris, clusterIDcol = iris$Species)[[3]])[1], dim(iris)[1])
})


# test withinclass_scattermatrix_LDA
test_that("calculating within class scatter matrix works correctly", {
  sp_dfs <- split_clusters(data = iris, clusterIDcol = iris$Species)
  # plan: rework main iDA function 
  #expect_equal(length(split_clusters(data = iris, clusterIDcol = iris$Species)), length(unique(iris$Species)))
  #expect_equal(dim(split_clusters(data = iris, clusterIDcol = iris$Species)[[1]]), dim(iris[which(iris$Species == "setosa"),]))
  #expect_equal(dim(split_clusters(data = iris, clusterIDcol = iris$Species)[[2]]), dim(iris[which(iris$Species == "versicolor"),]))
  #expect_equal(dim(split_clusters(data = iris, clusterIDcol = iris$Species)[[3]]), dim(iris[which(iris$Species == "virginica"),]))
  #expect_equal(iDA:::withinclass_scattermatrix_LDA(sp_dfs, diag = TRUE), dim(iris)[1])
})



# test betweenclass_scatter_matrix



# test decomposesvd


# test decomposeirlba



# test getSNN




# test getLouvain









