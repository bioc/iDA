#'  Find embedding space using iterative LDA
#'
#'  Takes scaled data and iterates between clustering using the Louvain 
#'  community detection method and embedding in LDA space, then recluster
#'  in the LDA transformed data space.
#'
#' @param var.data (matrix) A matrix of scaled count data for variable genes 
#' to find embedding for. (sample x variable features)
#' @param k.param (numeric) Defines k for the k-nearest neighbor algorithm 
#' (passed to [`getSNN`])
#' @param prune.SNN (numeric) Sets the cutoff for acceptable Jaccard index when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the strigency of pruning (0 --- no pruning, 1 ---
#' prune everything). Passed to [`getSNN`]
#' @param dims.use (numeric) A vector of the dimensions to use in construction 
#' of the SNN graph (e.g. To use the first 10 PCs, pass 10) 
#' (passed to [`getSNN`])
#' @param diag Diagonalize the within class scatter matrix (assume the features
#' are independent within each cluster)
#' @param c.param (numeric) Defines the number of desired clusters to be found
#' in the embedding
#' @param cluster.method What clustering method to use 
#'
#' @import irlba
#' @import igraph
#' @import plyr
#' @importFrom stats kmeans
#' @importFrom mclust adjustedRandIndex
#' @return n number of dataframes for each cluster's data

.iDA_core <- function(var.data,
                     k.param = 10,
                     prune.SNN = 1/15,
                     dims.use = 10,
                     diag = TRUE,
                     c.param = NULL,
                     cluster.method = "walktrap") {
    #calculate svd for covariance matrix of variable_features
    if (ncol(var.data) < dims.use){
        dims.use <- ncol(var.data)
    }
    svd <- svdr(as.matrix(var.data), k = dims.use)
    # transform data
    transformed <- svd$v
    rownames(transformed) <- colnames(var.data)
    #cluster
    if(cluster.method == "louvain") {
        snn <- getSNN(data.use = transformed, 
                      k.param = k.param, 
                      prune.SNN = prune.SNN)
        louvainClusters <- getLouvain(SNN = snn)
        . <- cbind(start = rep(1,dim(transformed)[1]), 
                          currentclust = louvainClusters)
    } else if (cluster.method == "kmeans"){
        kmeansclusters <- kmeans(transformed, centers = c.param)
        clusters <- cbind(start = rep(1,dim(transformed)[1]), 
                          currentclust = kmeansclusters$cluster)
    } else if (cluster.method == "walktrap"){
        snn <- getSNN(data.use = transformed, 
                      k.param = k.param, 
                      prune.SNN = prune.SNN)
        walktrapClusters <- cluster_walktrap(snn)
    }
    #pick highest modularity
    if (is.null(c.param)){
        modularity <- c() 
        for (i in 1:min(ncol(transformed) - 1 ,15)){
            # at max 15 clusters
            modularity <- c(modularity, modularity(snn, 
                                                   cut_at(walktrapClusters, 
                                                          no = i)))
        }
        maxmodclust <- cut_at(walktrapClusters, no = which.max(modularity))
        clusters <- as.data.frame(cbind(start = rep(1,nrow(transformed)), 
                                        currentclust = maxmodclust))
    } else if (is.numeric(c.param)) {
        maxmodclust <- cut_at(walktrapClusters, 
                              no = c.param)
        clusters <- as.data.frame(cbind(start = rep(1,nrow(transformed)), 
                                        currentclust = maxmodclust))
    } else {
        stop("Invalid c.param")
    }
    rownames(clusters) <- rownames(transformed)
    if (length(unique(clusters$currentclust)) == ncol(var.data)) {
        stop(paste("There are ", sum(table(clusters$currentclust) == 1), 
                   " clusters which only have one sample in them. 
                   Consider increasing k.param."))
    }
    #calculate concordance between last and current iteration's clustering 
    concordance <- adjustedRandIndex(clusters[,(ncol(clusters)-1)], 
                                     clusters[,(ncol(clusters))])
    #start iterations
    i = 1        
    while(concordance < .98) {
        if(i > 1){
            message("iteration ", i-1)
            message("concordance: ", concordance)
        }
        # split by cluster
        splitclusters <- split(x =  as.data.frame(t(var.data)),
                               f = factor(clusters[,ncol(clusters)]))
        # calculate within cluster scatter matrix
        Sw <- WCS(splitclusters = splitclusters, diag = diag)
        # calculate between cluster scatter matrix
        Sb <- BCS(splitclusters = splitclusters)
        # LDA on these scatter matrices (this is the time consuming step)
        decomp <- decomposeSVD(Sw, Sb, nu = length(splitclusters) - 1)
        eigenvecs <- decomp[["eigenvecs"]]
        #transform data
        eigenvectransformed <- t(var.data) %*% eigenvecs
        #calculate SNN matrix for top LDs
        if (cluster.method == "louvain") {
            snn <- getSNN(data.use = eigenvectransformed, 
                          k.param = k.param, 
                          prune.SNN = prune.SNN)
            louvainClusters <- getLouvain(SNN = as.matrix(snn))
            clusters <- cbind(clusters, currentclust = louvainClusters)
        } else if (cluster.method == "kmeans"){
            kmeansclusters <- kmeans(eigenvectransformed, 
                                     centers = c.param)
            clusters <- cbind(clusters, 
                              currentclust = kmeansclusters$cluster)
        } else if (cluster.method == "walktrap"){
            snn_transformed <- getSNN(data.use = eigenvectransformed, 
                                      k.param = k.param, 
                                      prune.SNN = prune.SNN)
            #cluster
            walktrapClusters <- cluster_walktrap(snn_transformed)
            #pick highest modularity 
            if (is.null(c.param)){
                modularity <- c()
                for (i in 1:min(ncol(transformed) - 1 ,15)){
                    # at max 15 clusters
                    modularity <- c(modularity, modularity(snn, 
                                                           cut_at(walktrapClusters, 
                                                                  no = i)))
                }
                maxmodclust <- cut_at(walktrapClusters, 
                                      no = which.max(modularity))
                clusters <- cbind(clusters, 
                                  currentclust = maxmodclust)
            } else if (is.numeric(c.param)) {
                maxmodclust <- cut_at(walktrapClusters, 
                                      no = c.param)
                clusters <- cbind(clusters, 
                                  currentclust = maxmodclust)
            } else {
                stop("Invalid c.param")
            }
        }
        concordance <- adjustedRandIndex(clusters[,(ncol(clusters)-1)], 
                                         clusters[,(ncol(clusters))])
        i = i + 1
    }
    geneweights <- as.data.frame(eigenvecs)
    rownames(geneweights) <- rownames(var.data)
    colnames(geneweights) <- paste("LD", seq(ncol(geneweights)), sep = "")
    rownames(eigenvectransformed) <- rownames(transformed)
    colnames(eigenvectransformed) <- paste("LD", 
                                           seq(ncol(eigenvectransformed)), 
                                           sep = "")
    message("final concordance: ")
    message(concordance)
    retlist <- list(clusters = clusters[, dim(clusters)[2]],
                    LDs = eigenvectransformed,
                    feature_weights = geneweights)
    return(retlist)
}
