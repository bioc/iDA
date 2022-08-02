#' Compute each cluster's within class scatter matrix
#'
#' Takes in the output from split_clusters() and computes the within class 
#' scatter matrix
#'
#' @param splitclusters A list of dataframes with scaled data from each cluster 
#' (output from split_clusters())
#' @param diag if off diagonal entries in within class scatter matrix should be
#'  zeroed
#' @importFrom dplyr select
#' @return returns the within class scatter matrix
WCS <- function(splitclusters, diag = FALSE) {
    ## calculate means vector for each cluster
    clustermeans <- c()
    k <- 1
    for (i in splitclusters) {
        clustermeans[[k]] <- colMeans(i)
        k = k + 1
    }
    
    ## calculate within class scatter matrix for each cluster
    wcsm <- c()
    k <- 1
    for (i in splitclusters) {
        dataMatrix <- t(i)
        wcsm[[k]] <- (t(i - clustermeans[[k]])) %*% (i - clustermeans[[k]])
        k <- k + 1
    }
    
    ## add all within class scatter matrices together
    Sw <- array(0L, dim(wcsm[[1]]))
    k <- 1
    
    list <- vapply(splitclusters, 
                   function(l) l[1], FUN.VALUE = list(numeric(1)))
    n_obs <- sum(lengths(list))
    
    for (i in wcsm) {
        Sw <- Sw + ((dim(splitclusters[[k]])[1]) / n_obs) * i
        k <- k + 1
    }
    
    if (diag == TRUE) {
        ## set off-diagonal entries to 0
        Sw <- diag(diag(Sw))
    }
    return(Sw)
}

#' Compute the between class scatter matrix
#'
#' Takes in a list of dataframes with scaled data (output from split()) and 
#' returns the between class scatter matrix
#' @param splitclusters A list of dataframes (from the output of split()) with 
#' scaled data from each cluster
#' @importFrom dplyr select
#' @return returns the between class scatter matrix
BCS <- function(splitclusters) {
    ## calculate means vector for each cluster
    clustermeans <- c()
    k <- 1
    for (i in splitclusters) {
        clustermeans[[k]] <- colMeans(i)
        k <- k + 1
    }
    ## calculate overallMeans for each feature
    overallMeanVector <- c()
    for (i in seq_along(clustermeans[[1]])) {
        overallMeanVector[[i]] <- mean(vapply(clustermeans, 
                                              function(l) l[[i]], 
                                              FUN.VALUE = 0))
    }
    ## calculate each btsc matrix per cluster
    btsc <- c()
    for (i in seq_along(clustermeans)) {
        btsc[[i]] <- ((clustermeans[[i]] - unlist(overallMeanVector)) %*%
                          t(clustermeans[[i]] - unlist(overallMeanVector)))
    }
    
    ## add all btsc's together
    Sb <- array(0L, dim(btsc[[1]]))
    k <- 1
    for (i in btsc) {
        Sb <- Sb + i
        k <- k + 1
    }
    return(Sb)
}

#' Compute LDA using within and between cluster covariance matrices
#'
#' Takes in within and between cluster covariance matrices (output from WCS
#' and BSC) and returns the eigenvectors and eigenvalues for the decomp.
#' @param WCSmat Output from WCS (matrix) of within class scatter 
#' @param BCSmat Output from BCS (matrix) of between class scatter 
#' @param nu The number of columns (eigenvectors) to keep 
#' @importFrom dplyr select
#' @return returns the between class scatter matrix
decomposeSVD <- function(WCSmat,
                         BCSmat,
                         nu = 10) {
    svd <- svd(solve(WCSmat) %*% BCSmat, nu)
    top_eigenvectors <- svd$u[, seq(nu)]
    top_eigenvalues <- svd$d[seq(nu)]
    return(list(eigenvecs = top_eigenvectors, eigenvalues = top_eigenvalues))
}


#' Cluster Determination
#'
#' Calculate k-nearest neighbors and construct a shared nearest neighbor 
#' (SNN) graph.
#'
#' @param data.use (matrix) Matrix with scaled data to find nearest neighbors
#' @param k.param (numeric) Defines k for the k-nearest neighbor algorithm
#' @param prune.SNN (numeric) Sets the cutoff for acceptable Jaccard index when
#'  computing the neighborhood overlap for the SNN construction. Any edges with
#'  values less than or equal to this will be set to 0 and removed from the SNN
#'  graph. Essentially sets the strigency of pruning (0 --- no pruning, 1 ---
#'  prune everything).
#' building KNN graph
#' @import igraph
#' @import scran
#' @return The SNN graph (igraph object)
getSNN <- function(data.use,
                   k.param = 10,
                   prune.SNN = 1/15) {
    data.use <- as.matrix(data.use)
    n.obs <- nrow(x = data.use)
    if (n.obs < k.param) {
        warning("k.param set larger than number of cells. 
            Setting k.param to number of cells - 1.",
                call. = FALSE)
        k.param <- n.obs - 1
    }
    # Smaller 'k' usually yields finer clusters
    SNN_igraph <- scran::buildKNNGraph(data.use, k = k.param, transposed = TRUE)
    snn.matrix <- similarity(SNN_igraph, method = "jaccard")
    snn.matrix[snn.matrix < prune.SNN] <- 0
    rownames(x = snn.matrix) <- rownames(x = data.use)
    colnames(x = snn.matrix) <- rownames(x = data.use)
    snn.graph <- graph_from_adjacency_matrix(snn.matrix, 
                                             weighted = TRUE, 
                                             mode = "undirected")
    return(snn.graph)
}

#' Louvain Cluster Determination
#'
#' Identify clusters of cells by a shared nearest neighbor (SNN) modularity
#' optimization based clustering algorithm. Optimize the modularity function to
#' determine clusters. For a full description of the algorithms, see Waltman and
#' van Eck (2013) \emph{The European Physical Journal B}.
#'
#' @param data.use (matrix) Matrix with scaled data to find nearest neighbors
#' @param k.param (numeric) Defines k for the k-nearest neighbor algorithm
#' @param prune.SNN (numeric) Sets the cutoff for acceptable Jaccard index when
#'  computing the neighborhood overlap for the SNN construction. Any edges with
#'  values less than or equal to this will be set to 0 and removed from the SNN
#'  graph. Essentially sets the strigency of pruning (0 --- no pruning, 1 ---
#'  prune everything).
#' @importFrom NetworkToolbox louvain
#' @return a list of identities for clustering
getLouvain <- function(data.use, k.param, prune.SNN){
    snn <- getSNN(data.use = data.use, 
                  k.param = k.param, 
                  prune.SNN = prune.SNN)
    louvain_clusters <- cluster_louvain(snn)
    idents <- louvain_clusters$membership
    return(idents)
}

#' Set Feature Weights from iDA
#'
#' 
#' @param object (SE or SCE) an SE or SCE object in which to fill in feature 
#' weights
#' @param feature.weights (matrix) matrix of features x LDs which is output 
#' from iDA 
#' @return an SE or SCE object with rowData() for the iDA feature weights
setFeatureWeights <- function(object, feature.weights){
    #add column indicating if a gene is an iDA variable feature or not
    rowData(object)[rownames(feature.weights), "iDA_var.feature"] <- TRUE
    false_rows <- is.na(rowData(object)[["iDA_var.feature"]])
    rowData(object)[false_rows, "iDA_var.feature"] <- FALSE
    for (cn in colnames(feature.weights)) { 
        rowData(object)[[cn]] <- 0
        keep_rows <- rownames(feature.weights)
        rowData(object)[keep_rows, ][[cn]] <- feature.weights[,cn]
        }
    return(object)
}


