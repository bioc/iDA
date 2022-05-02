#' Generic method to input data to iDA
#'
#' @param object The object to run iDA on
#' @param ... Additional arguments passed to object constructors
#' @return iDA output with clustering, gene weights, and cell weights
#' @export
setGeneric("iDA", signature=c("object"),
           function(object, ...) standardGeneric("iDA"))


#' Set method for SummarizedExperiment to input data to iDA
#'
#' @param object The object to run iDA on
#' @param ... Additional arguments passed to object constructors
#' @import airway
#' @import SummarizedExperiment
#' @import DESeq2
#' @importFrom genefilter rowVars
#' @importFrom utils head
#' @return iDA output with clustering, gene weights, and cell weights
#' @examples 
#' data(airway, package="airway")
#' se <- airway
#' set.seed(11)
#' se <- iDA(se)
#' 
#' @export
setMethod("iDA", "SummarizedExperiment",
          function(object, ...) {
              # Filtering counts < 10
              keep <- rowSums(assays(object)[["counts"]]) >= 10
              object <- object[keep,]
              #normalize counts
              dds <- DESeqDataSet(object, design = ~ 1)
              dds <- estimateSizeFactors(dds)
              #variance stabilizing transformation
              message("Transforming counts with vst().")
              scale.data <- assay(vst(dds, blind = TRUE))
              rowvars <- rowVars(scale.data)
              names(rowvars) <- rownames(scale.data)
              topVarGenes <- names(head(rowvars[order(-rowvars)], n = 2000))
              message(length(topVarGenes), 
                      " variable features found using rowVars(). \n")
              var.data <- logcounts[topVarGenes,]
              iDAoutput <- .iDA_core(var.data, ...)
              #add metadata back to object
              if (all(rownames(colData(object)) == rownames(iDAoutput$LDs))) {
                  colData(object) <- cbind(colData(object), 
                                           iDAoutput$LDs, 
                                           iDAoutput$clusters)
              }
              return(object)
          })

#' Set method for DESeqDataSet to input data to iDA
#'
#' @param object The object to run iDA on
#' @param ... Additional arguments passed to object constructors
#' @return iDA output with clustering, gene weights, and cell weights
#' @import DESeq2
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom genefilter rowVars
#' @importFrom utils head
#' @examples 
#' countData <- matrix(1:100,ncol=4)
#' condition <- factor(c("A","A","B","B"))
#' dds <- DESeq2::DESeqDataSetFromMatrix(countData, S4Vectors::DataFrame(condition), ~ condition)
#' set.seed(11)
#' dds <- iDA(dds)
#' 
#' @export
setMethod("iDA", "DESeqDataSet",
          function(object, ...) {
              # Filtering counts < 10
              keep <- rowSums(counts(object)) >= 10
              object <- object[keep,]
              #normalize counts
              object <- estimateSizeFactors(object)
              #variance stabilizing transformation
              message("Transforming counts with vst().")
              scale.data <- assay(vst(object, blind = TRUE))
              rowvars <- rowVars(scale.data)
              names(rowvars) <- rownames(scale.data)
              topVarGenes <- names(head(rowvars[order(-rowvars)], n = 2000))
              message(length(topVarGenes), 
                      " variable features found using rowVars(). \n")
              var.data <- scale.data[topVarGenes,]
              iDAoutput <- .iDA_core(var.data, ...)
              #add metadata back to object
              if (all(rownames(colData(object)) == rownames(iDAoutput$LDs))) {
                  colData(object) <- cbind(colData(object), 
                                           iDAoutput$LDs, 
                                           iDAoutput$clusters)
              }
              return(object)
          })

#' Method for SingleCellExperiment object to input data to iDA
#'
#' @param object The single cell experiment object to run iDA on
#' @param ... Additonal arguments passed to object constructors
#' @import  SingleCellExperiment 
#' @importFrom SummarizedExperiment assays
#' @importFrom scuttle normalizeCounts
#' @return SingleCellExperiment object with iDA cell weights and gene weights 
#' stored in reducedDims and cluster assignments stored in rowLabels
#' @examples 
#' data("sc_sample_data")
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(sc_sample_data)))
#' set.seed(11)
#' sce <- iDA(sce)
#' 
#' @export
setMethod("iDA", "SingleCellExperiment",
          function(object, ...) {
              if (!('logcounts' %in% names(assays(object)))){
                  logcounts(object) <- normalizeCounts(object,  
                                            size.factors = sizeFactors(object))
              }
              normcounts <-  logcounts(object)
              if (nrow(normcounts) < 2000) {
                  var.features <- rownames(normcounts)
              } else {
                  stats <- modelGeneVar(normcounts)
                  var.features <- getTopHVGs(stats, n = 2000)
              }
              message(length(var.features), 
                      " variable features found using scran::getTopHVGs. \n")
              var.data <- normcounts[var.features, ]
              iDA_sce <- .iDA_core(var.data, ...)
              reducedDims(object) <- list(iDAcellweights = iDA_sce[["LDs"]])
              colLabels(object) <- list(iDAclusters = iDA_sce[["clusters"]])
              return(object)
          })

#' Method for Seurat object to input data to iDA
#'
#' @param object The single cell experiment object to run iDA on
#' @param assay The assay to take counts from
#' @param selection.method The Seurat method to use for variable feature 
#' selection. 
#' @param ... Additional arguments passed to object constructors
#' @import Seurat 
#' @return Seurat object with iDA cell weights and gene weights stored in 
#' object[["iDA"]] and cluster assignments stored in rowLabels
#' 
#' @examples
#' data("pbmc_small")
#' force(pbmc_small)
#' set.seed(11)
#' pbmc_small <- iDA(pbmc_small, assay = "RNA")
#' @export
setMethod("iDA", "Seurat",
          function(object, assay, selection.method = "dispersion", ...) {
              if (length(GetAssayData(object, slot = "scale.data")) == 0){
                  object <- NormalizeData(
                      object, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000
                  ) #normalize data
                  object <- ScaleData(object) #center and scale
                  scale.data <- GetAssayData(object, slot = "scale.data") 
              } else {
                  scale.data <- GetAssayData(object, slot = "scale.data") 
              }
              if (length(VariableFeatures(object)) == 0) {
                  object <- FindVariableFeatures(object, 
                                                 selection.method = selection.method)
              }
              var.features <- VariableFeatures(object)
              message(length(var.features), 
                      " variable features found using FindVariableFeatures(). \n")
              var.data <- scale.data[var.features,]
              iDA_seurat <- .iDA_core(var.data, ...)
              object[["iDA"]] <- CreateDimReducObject(embeddings = as.matrix(iDA_seurat[["LDs"]]),
                                                      key = "LD",
                                                      loadings = as.matrix(iDA_seurat[["feature_weights"]]),
                                                      assay = assay)
              object <- AddMetaData(object = object, 
                                    metadata = iDA_seurat[["clusters"]], 
                                    col.name = "iDA_clust")
              return(object)
          })
