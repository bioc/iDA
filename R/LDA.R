#' Find Variable Genes
#'
#' Identify rows in a data frame with high sclaed dispersion. This rule was
#' taken from the Seurat package.
#'
#' @param NormCounts (data.frame) features are rows, samples are columns. 
#' Rows must be named.
#' @param dispersion.cutoff (numeric) rows returned will have scaled dispersion 
#' higher provided cutoff
#' @param mean.low.cutoff (numeric) rows returned will have average higher than 
#' this cutoff
#' @param mean.high.cutoff (numeric) rows returned will have average lower than
#'  this cutoff
#' @importFrom stats sd var
#' @return (character) a list of row names with high dispersion rows
VariableGenesGeneric <- function(NormCounts,
                                 dispersion.cutoff = 1,
                                 mean.low.cutoff = 0.1,
                                 mean.high.cutoff = 8) {
    ## calculate logged means and VMR
    ExpMeans <- apply(NormCounts, 1, FUN = function(x) log(mean(exp(x) - 1) + 1))
    finite_idx <- is.finite(ExpMeans)
    data.use <- NormCounts[finite_idx, ]
    ExpMeans <- ExpMeans[finite_idx]
    dispersions <- apply(data.use, 1, FUN = function(x) {log(var(exp(x) - 1) / mean( exp(x) - 1))})
    dispersions <- dispersions[finite_idx]
    dispersions[is.na(x = dispersions)] <- 0
    ExpMeans[is.na(x = ExpMeans)] <- 0
    num.bin <- 20
    data.x.bin <- cut(x = ExpMeans, breaks = num.bin)
    names(x = data.x.bin) <- names(x = ExpMeans)
    mean.y <- tapply(X = dispersions, INDEX = data.x.bin, FUN = mean)
    sd.y <- tapply(X = dispersions, INDEX = data.x.bin, FUN = sd)
    
    ## scale dispersions
    scaled.dispersions <- (dispersions - mean.y[as.numeric(x = data.x.bin)]) /
        sd.y[as.numeric(x = data.x.bin)]
    names(x = scaled.dispersions) <- names(x = ExpMeans)
    ## find variable features
    variable_idx <- scaled.dispersions > dispersion.cutoff &
        ExpMeans > mean.low.cutoff &
        ExpMeans < mean.high.cutoff
    var.features <- names(dispersions[variable_idx])
    retlist <- list(
        "dispersions" = data.frame(scaled.dispersions, ExpMeans),
        "use.data" = data.use,
        "var.features" = var.features)
    return(retlist)
}

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
        wcsm[[k]] <- (t(t(dataMatrix) - clustermeans[[k]])) %*% (t(dataMatrix) - clustermeans[[k]])
        k <- k + 1
    }
    
    ## add all within class scatter matrices together
    Sw <- array(0L, dim(wcsm[[1]]))
    k <- 1
    
    list <- vapply(splitclusters, function(l) l[1], FUN.VALUE = list(numeric(1)))
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

#' Cluster Determination
#'
#' Identify clusters of cells by a shared nearest neighbor (SNN) modularity
#' optimization based clustering algorithm. Optimize the modularity function to
#' determine clusters. For a full description of the algorithms, see Waltman and
#' van Eck (2013) \emph{The European Physical Journal B}.
#'
#' @param SNN a matrix of shared nearest neighbors (output from getSNN)
#' @importFrom NetworkToolbox louvain
#' @return a list of identities for clustering
getLouvain <- function(SNN){
    louvain_clusters <- cluster_louvain(SNN)
    idents <- louvain_clusters$membership
    return(idents)
}

