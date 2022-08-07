#' Generic method to input data to iDA
#'
#' @param object The object to run iDA on
#' @param ... Additional arguments passed to object constructors
#' @return iDA output with clustering, gene weights, and cell weights
#' @examples 
#' 
#' @export
setGeneric("iDA", signature=c("object"),
           function(object, ...) standardGeneric("iDA"))


#' Set method for matrix input to iDA
#'
#' @param object The normalized matrix (feature x sample) to run iDA on
#' @param nFeatures The number of HVG features to use in reduction
#' @param ... Additional arguments passed to object constructors
#' @importFrom genefilter rowVars
#' @importFrom utils head
#' @return iDA output with clustering, gene weights, and cell weights as a list
#' @examples 
#' exp <- matrix(rpois(20000, 5), ncol=20)
#' colnames(exp) <- paste0("donor", seq_len(ncol(exp)))
#' rownames(exp) <- paste0("gene", seq_len(nrow(exp)))
#' logexp <- logexp <- log2(exp + 1)
#' set.seed(11)
#' ida <- iDA(logexp)
#' 
#' @export
setMethod("iDA", "matrix",
          function(object, nFeatures = 2000, ...) {
              #center and scale
              scale.data <- scale(object, center = TRUE, scale = TRUE)
              
              rowvars <- rowVars(scale.data)
              names(rowvars) <- rownames(scale.data)
              topVarGenes <- names(head(rowvars[order(-rowvars)], n = nFeatures))
              message(length(topVarGenes), 
                      " variable features found using rowVars(). \n")
              var.data <- scale.data[topVarGenes,]
              iDAoutput <- .iDA_core(var.data, ...)

              return(iDAoutput)
          })


#' Set method for SummarizedExperiment to input data to iDA
#'
#' @param object The object to run iDA on
#' @param nFeatures The number of HVG features to use in reduction
#' @param ... Additional arguments passed to object constructors
#' @import airway
#' @import SummarizedExperiment
#' @import DESeq2
#' @importFrom genefilter rowVars
#' @importFrom utils head
#' @return iDA output with clustering, gene weights, and cell weights
#' @examples 
#' exp <- matrix(rpois(20000, 5), ncol=20)
#' colnames(exp) <- paste0("donor", seq_len(ncol(exp)))
#' rownames(exp) <- paste0("gene", seq_len(nrow(exp)))
#' logexp <- logexp <- log2(exp + 1)
#' conditions <- factor(rep(1:4, 5))
#' se <- SummarizedExperiment::SummarizedExperiment(
#' assays = list(counts = exp, logcounts = logexp),
#' colData = data.frame(conditions = conditions))
#' set.seed(11)
#' se <- iDA(se)
#' 
#' @export
setMethod("iDA", "SummarizedExperiment",
          function(object, nFeatures = 2000, ...) {
              # Filtering counts < 10
              if (!is.null(assays(object)[["logcounts"]])){
                  scale.data <- assays(object)[["logcounts"]]
              } else if(!is.null(assays(object)[["counts"]])) {
                  keep <- rowSums(assays(object)[["counts"]]) >= 10
                  object <- object[keep,]
    
                  #normalize counts
                  dds <- DESeqDataSet(object, design = ~ 1)
                  dds <- estimateSizeFactors(dds)
                  #variance stabilizing transformation
                  message("Transforming counts with varianceStabilizingTransformation().")
                  scale.data <- assay(varianceStabilizingTransformation(dds, 
                                                                  blind = TRUE))
              } else {
                  stop("Did not find 'counts' or 'logcounts' in assays().")
              }
              rowvars <- rowVars(scale.data)
              names(rowvars) <- rownames(scale.data)
              topVarGenes <- names(head(rowvars[order(-rowvars)], n = nFeatures))
              message(length(topVarGenes), 
                      " variable features found using rowVars(). \n")
              var.data <- scale.data[topVarGenes,]
              iDAoutput <- .iDA_core(var.data, ...)
              #add metadata back to object
              if(is.null(iDAoutput)){
                  message("Only one cluster found and no LDA to compute. \nIf this is an error, please try adjusting the clustering parameters.")
              } else {
              if (all(rownames(colData(object)) == rownames(iDAoutput$LDs))) {
                  colData(object) <- cbind(colData(object), 
                                           iDAoutput$LDs, 
                                           "iDA_clusters" = iDAoutput$clusters)
                  }
              }
              object <- setFeatureWeights(object, iDAoutput[["feature_weights"]])
              return(object)
          })

#' Set method for DESeqDataSet to input data to iDA
#'
#' @param object The object to run iDA on
#' @param nFeatures The number of HVG features to use in reduction
#' @param ... Additional arguments passed to object constructors
#' @return iDA output with clustering, gene weights, and cell weights
#' @import DESeq2
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom genefilter rowVars
#' @importFrom utils head
#' @importFrom S4Vectors DataFrame
#' 
#' @export
setMethod("iDA", "DESeqDataSet",
          function(object, nFeatures = 2000, ...) {
              # Filtering counts < 10
              keep <- rowSums(counts(object)) >= 10
              object <- object[keep,]
              #normalize counts
              object <- estimateSizeFactors(object)
              #variance stabilizing transformation
              message("Transforming counts with vst().")
              scale.data <- varianceStabilizingTransformation(object)
              rowvars <- rowVars(assay(scale.data))
              names(rowvars) <- rownames(scale.data)
              if(is.null(names(head(rowvars[order(-rowvars)], n = nFeatures)))) {
                  stop("Please set rownames() to gene IDs or other identifier.")
              }
              topVarGenes <- names(head(rowvars[order(-rowvars)], n = nFeatures))
              message(length(topVarGenes), 
                      " variable features found using rowVars(). \n")
              var.data <- assay(scale.data)[topVarGenes,]
              iDAoutput <- .iDA_core(var.data, ...)
              #add metadata back to object
              if (is.null(iDAoutput)){
                message("Only one cluster found and no LDA to compute. \nIf this is an error, please try adjusting the clustering parameters.")
              } else if (all(rownames(colData(object)) == rownames(iDAoutput$LDs))) {
                  colData(object) <- cbind(colData(object), 
                                           iDAoutput$LDs, 
                                           "iDA_clusters" = as.factor(iDAoutput$clusters))
              }
              object <- setFeatureWeights(object, iDAoutput[["feature_weights"]])
              return(object)
          })

#' Method for SingleCellExperiment object to input data to iDA
#'
#' @param object The single cell experiment object to run iDA on
#' @param nFeatures The number of HVG features to use in reduction
#' @param ... Additional arguments passed to object constructors
#' @import  SingleCellExperiment 
#' @importFrom SummarizedExperiment assays
#' @importFrom scuttle normalizeCounts
#' @return SingleCellExperiment object with iDA cell weights and gene weights 
#' stored in reducedDims and cluster assignments stored in rowLabels
#' @examples 
#' library(SingleCellExperiment)
#' data(sc_sample_data, package = "scPipe")
#' sce <- SingleCellExperiment(assays = list(counts = as.matrix(sc_sample_data)))
#' set.seed(11)
#' sce <- iDA(sce)
#' 
#' @export
setMethod("iDA", "SingleCellExperiment",
          function(object, nFeatures = 2000, ...) {
              if (!('logcounts' %in% names(assays(object)))){
                  logcounts(object) <- normalizeCounts(object,  
                                            size.factors = sizeFactors(object))
              }
              normcounts <-  logcounts(object)
              if (nrow(normcounts) < 2000) {
                  var.features <- rownames(normcounts)
              } else {
                  stats <- modelGeneVar(normcounts)
                  var.features <- getTopHVGs(stats, n = nFeatures)
              }
              message(length(var.features), 
                      " variable features found using scran::getTopHVGs. \n")
              var.data <- normcounts[var.features, ]
              iDA_sce <- .iDA_core(var.data, ...)
              if (is.null(iDA_sce)){
                  message("Only one cluster found and no LDA to compute. \nIf this is an error, please try adjusting the clustering parameters.") 
              } else {
              reducedDims(object) <- list(iDAcellweights = iDA_sce[["LDs"]])
              #reducedDims(object) <- list(iDAgeneweights = iDA_sce[["feature_weights"]])
              colData(object)[["iDAclusters"]] <- iDA_sce[["clusters"]]
              object <- setFeatureWeights(object, iDA_sce[["feature_weights"]])
              return(object)
              }
          })

