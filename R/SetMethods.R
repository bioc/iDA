#' Generic method to input data to iDA
#'
#' @param object The object to run iDA on
#' @param ... Additonal arguments passed to object constructors
#' @return iDA output with clustering, gene weights, and cell weights
#' @export
setGeneric("iDA", signature=c("object"),
           function(object, ...) standardGeneric("iDA"))

#' Set method for matrix to input data to iDA
#'
#' @param object The object to run iDA on
#' @param scaled Is the input matrix normalized and scaled already?
#' @param ... Additonal arguments passed to object constructors
#' @return iDA output with clustering, gene weights, and cell weights
#' @export
setMethod("iDA", "matrix",
          function(object, scaled, ...) {
            if (scaled == FALSE){
              message("scaled = FALSE. Normalizing and scaling matrix now.")
              data.use.norm <- Matrix::t((Matrix::t(object)/ Matrix::colSums(object))* 10000)
              data.use.norm <- log1p(data.use.norm)
              data.use <- scale(data.use.norm)
            }
            iDAoutput <- iDA_core(data.use, scaled = TRUE, ...)
            return(iDAoutput)
          })

#' Method for SingleCellExperiment object to input data to iDA
#'
#' @param object The single cell experiment object to run iDA on
#' @param selection.method The method to use for finding variable features
#' @param ... Additonal arguments passed to object constructors
#' @import  SingleCellExperiment 
#' @importFrom SummarizedExperiment assays
#' @importFrom scuttle normalizeCounts
#' @return SingleCellExperiment object with iDA cell weights and gene weights stored in reducedDims and cluster assignemts
#' stored in rowLabels
#' @export
setMethod("iDA", "SingleCellExperiment",
          function(object, selection.method = "scran", ...) {
              if (!('logcounts' %in% names(assays(object)))){
                logcounts(object) <- normalizeCounts(object,  size.factors = sizeFactors(object))
              }
              normcounts <-  logcounts(object)
              
              if (selection.method == "scran") {
                if (nrow(normcounts) < 3000) {
                  var.features <- rownames(normcounts)
                } else {
                  stats <- scran::modelGeneVar(normcounts)
                  var.features <- scran::getTopHVGs(stats, n = 3000)
                }
              }
              
              iDA_sce <- iDA_core(normcounts, NormCounts = normcounts, var.features = var.features, scaled = TRUE, ...)
              reducedDims(object) <- list(iDAcellweights = iDA_sce[[2]])
              colLabels(object) <- list(iDAclusters = iDA_sce[[1]])
              #rowData(object[iDA_sce[[4]],]) <- list(iDAgeneweights = iDA_sce[[3]])
              return(object)
          })

#' Method for Seurat object to input data to iDA
#'
#' @param object The single cell experiment object to run iDA on
#' @param assay The assay to take counts from
#' @param selection.method The Seurat method to use for variable feature selection. 
#' @param ... Additional arguments passed to object constructors
#' @import Seurat Seurat
#' @return Seurat object with iDA cell weights and gene weights stored in object[["iDA"]] and cluster assignments stored in rowLabels
#' @export

setMethod("iDA", "Seurat",
          function(object, assay, selection.method = "vst", ...) {
            if (length(GetAssayData(object, slot = "scale.data")) == 0){
              object <- NormalizeData(object, 
                                      normalization.method = "LogNormalize", 
                                      scale.factor = 10000) #normalize data
              object <- ScaleData(object) #center and scale
              counts <- GetAssayData(object, slot = "scale.data") 
              normCounts <- GetAssayData(object, slot = "data") 
            } else {
              counts <- GetAssayData(object, slot = "scale.data") 
              normCounts <- GetAssayData(object, slot = "data") 
            }
            
            if (length(VariableFeatures(object)) == 0) {
              object <- FindVariableFeatures(object, selection.method = selection.method)
            }
            var.features <- VariableFeatures(object)
            iDA_seurat <- iDA_core(counts, NormCounts = normCounts, scaled = TRUE, var.features = var.features, ...)

            object[["iDA"]] <- CreateDimReducObject(embeddings = as.matrix(iDA_seurat[[2]]),
                                                    key = "LD_",
                                                    loadings = as.matrix(iDA_seurat[[3]]),
                                                    assay = assay)
            
            object <- AddMetaData(object = object, metadata = iDA_seurat[[1]], col.name = "iDA_clust")
            object[[assay]]@var.features <- iDA_seurat[[4]]
            return(object)
          })
