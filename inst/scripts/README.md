
# Overview

The iDA package contains functions to create an embedding to reduce the
dimensionality of high dimensional data as well as a clustering
assignment which optimally separates latent clusters. This method
leverages linear discriminant analysis (LDA) in order to find
transformations which both separate classes from each other while also
minimizing variance within each class. These classes can arise from both
known (measured) and unknown (unmeasured) sources. This embedding can be
used as an interpretable reduction to characterize and investigate
features which drive class separation or to correct for unwanted sources
of variation in downstream analysis.

iDA can take in either a scRNA-seq dataset (as a SingleCellExperiment
object) or a bulk-RNAseq dataset (in either a DESeqDataSet or
SummarizedExperiment). iDA will check for needed normalization, do a
feature selection step, and then perform the iterative LDA. The iDA
method is performs the following steps:

------------------------------------------------------------------------

**Algorithm 1**: iDA Algorithm

------------------------------------------------------------------------

![](./iDAalgo.jpg) \# Installation

``` r
# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install iDA
BiocManager::install("iDA")
```

# Libraries

``` r
library(iDA)
library(scPipe)
library(celldex)
library(scater)
library(ggplot2)
library(Rtsne)
library(dplyr)
library(gridExtra)
```

# Applying iDA to scRNA-seq data to estimate cell types and top discriminative markers

A common motivation of studies using scRNA-seq is to detect novel
cell-types and identify marker genes that define these cell types. The
current standard scRNA approach is to use PCA to reduce the
computational cost of downstream analyses like unsupervised clustering
for putative cell-type identification or non-linear transformations like
tSNE or UMAP to visualize global cell-type differences.

Using PCA objective to detect latent discrete structure assumes that
latent classes will in fact account for most of the dataset variance.
However, since PCA is unsupervised, it does not necessarily separate
classes. A more useful approach would be to find an embedding which
defines the latent structure of the separation of cell types. iDA uses
linear discriminant analysis (LDA) in an iterative approach to maximize
the within-to-between cluster variance. The result is 1. vectors for
each LD with gene weights which are class-discriminative and 2. a vector
of cluster assignments for each cell. We can then evaluate the top
weighted genes in each LD to identify putative cell types.

Since iDA uses a single value decomposition (SVD) of the gene x gene
scatter matrices in each iterative step, this method scales relatively
well for increasing cell number counts, but not as well for large gene
numbers. To help with computational time, we perform a features
selection step at the beginning of the method to to limit the size of
the scatter matrices.

This dataset is an experimental dataset in the scPipe software generated
by Dr Christine Biben. This experiment profiles gene expression in 383
blood cells for a subset of 1000 genes. The rows names are Ensembl gene
ids and the columns are cell names, which is the wall position in the
384 plates.

## QC data

Since these are raw counts, we will first go through some QC as
recommended by OSCA for scRNAseq pre-processing.

``` r
data("sc_sample_data")
data("sc_sample_qc")
sce <- SingleCellExperiment(assays = list(counts = as.matrix(sc_sample_data)))
organism(sce) <- "mmusculus_gene_ensembl"
gene_id_type(sce) <- "ensembl_gene_id"
QC_metrics(sce) <- sc_sample_qc
demultiplex_info(sce) <- cell_barcode_matching
UMI_dup_info(sce) <- UMI_duplication
sce <- detect_outlier(sce)
sce <- remove_outliers(sce)
```

## iDA decomposition

To run iDA on a SingleCellExperiment object, we simply use the
following:

``` r
start <- Sys.time()
sce <- iDA(sce)
#> 1000 variable features found using scran::getTopHVGs.
#> final concordance:
#> 0.996689642042347
end <- Sys.time()
```

For each iteration, iDA will display a message of the concordance
between the previous iteration’s clustering and the current iteration’s
clustering. iDA will stop iterating when concordance is greater than
98%. For this dataset, iDA only iterates one time for a final
concordance of 99.7%.

    #> iDA ran in 20.96 seconds for the SCE object with 361 cells and 1000 variable features.

We can then evaluate the clustering and LDs produced from iDA. The
cluster labels are stored as a column `iDAclusters` in the colData() of
the SingleCellExperiment object, and the LDs for the cell weights are
stored as a reducedDims() under the `iDAcellweights` slot.

``` r
#cluster assignments
sce_iDAident <- colData(sce)[["iDAclusters"]]
#iDA cell embedding
sce_iDAembedding <- as.data.frame(reducedDims(sce)[["iDAcellweights"]])
```

There are three clustering methods available to use in the iDA
algorithm: louvain community detection, walktrap community detection, or
k-means clustering. These can be set in the cluster.method argument:

``` r
sce <- iDA(sce, cluster.method = c("louvain", "walktrap", "kmeans"))
```

You can also change the `nFeatures` parameter to use a different number
of HVGs for the embedding. The computation time does increase as
nFeatures increases since the embedding is computed with a single value
decomposition (SVD) of the gene x gene scatter matrix in each iteration.
Because of this though, the computation time scales relatively well for
increasing numbers of samples. the `nFeatures` parameter can be from the
default number (2000 HVGs) with the following:

``` r
sce <- iDA(sce, nFeatures = 3000)
```

### tSNE of iDA decomposition

Since iDA is still a linear dimensionality reduction method and the
dimension will likely be larger than 2 (the dimension will always be k -
1), we may still want to pair this with a non-linear method such as tSNE
to visualize the clustering in a 2D plot. We can now treat these LDs
just as we would PCs as input to a non-linear dimensionality reduction
to reduce this to 2 dimensions for cluster visualization.

``` r
set.seed(11)
tsne_iDA <- as.data.frame(Rtsne(sce_iDAembedding, pca = FALSE)$Y)

piDA_sce <- ggplot(tsne_iDA, aes(x = V1, y = V2, col = as.factor(sce_iDAident))) + 
  geom_point() + 
  xlab("tSNE 1") +
  ylab("tSNE 2") +
  labs(color = "iDA Cluster") +
  ggtitle("tSNE of the top LDs, colored by iDA Cluster")

piDA_sce
```

<img src="README_files/figure-gfm/unnamed-chunk-11-1.png" width="100%" />

## Comparing to PCA decomposition

To compare what information we are gaining from iDA over a PCA
reduction, we’ll look at the tSNE of the PCA reduction to compare
clusters. We will run scater’s `calculatePCA()` and then `Rtsne()` on
the top 10 PCs.

``` r
sce_pca <- calculatePCA(sce)
```

``` r
set.seed(11)
tsne_PCAsce <- as.data.frame(Rtsne(sce_pca[,1:10], pca = FALSE)$Y)

pPCA_sce <- ggplot(tsne_PCAsce, aes(x = V1, y = V2, col = as.factor(sce_iDAident))) + 
  geom_point() + 
  xlab("tSNE 1") +
  ylab("tSNE 2") +
  labs(color = "iDA Cluster") +
  ggtitle("tSNE of the top PCs, colored by iDA Cluster")

pPCA_sce
```

<img src="README_files/figure-gfm/unnamed-chunk-13-1.png" width="100%" />

Given this dataset with only highly variable genes, the PCA does a
pretty good job of finding the separation between the discrete types,
but we note a few key differences.

First, the two putative cell types in clusters 8 (purple) and 9 (pink)
have an ill-defined boundary in the PCA embedding. The iDA embedding
offers a much improved within-cluster distance in each cluster and also
describes a much better boundary between these two putative cell types.

Similarly, the separation between clusters 4 (kelly green) and 6 (light
blue) is much better defined in the iDA embedding compared to the PCA
one.

Lastly, we also see a much better boundary defined between clusters 1
(salmon) and 2 (orange) in the iDA embedding compared to the PCA one.

We can evaluate the iDA embedding to see what dimension is defining this
cluster split and what genes are highly weighted for this split.

## iDA embedding cell weights

For example, let’s evaluate the LD’s for which we see separation between
cluster 4 and 6. We can evaluate the highly weighted genes for those
dimensions to see the expression differences between the clusters.

``` r
ggplot(sce_iDAembedding, aes(x = LD7, y = seq(1, nrow(sce_iDAembedding)), color = as.factor(sce_iDAident))) + 
    geom_point()  + 
    labs(colour="iDA Cluster") + 
    xlab( "LD7 Weight") + 
    ylab("Cell Index")
```

<img src="README_files/figure-gfm/unnamed-chunk-14-1.png" width="100%" />

We see LD 7 is the dimension which separates clusters 4 (kelly green)
and 6 (light blue). We could then investigate the gene weights for LD 7
to identify genes which help in defining this cluster separation.

# Applying iDA to SummarizedExperiment object

This dataset is from the celldex package and contains normalized
expression values for 259 bulk RNA-seq samples generated by Blueprint
and ENCODE from pure populations of stroma and immune cells (Martens and
Stunnenberg, 2013; The ENCODE Consortium, 2012). The samples were
processed and normalized as described in Aran, Looney and Liu et
al. (2019).

``` r
BED <- BlueprintEncodeData(rm.NA = "none",
                           ensembl = FALSE,
                           cell.ont = c("nonna"))
#> snapshotDate(): 2021-10-19
#> see ?celldex and browseVignettes('celldex') for documentation
#> loading from cache
#> see ?celldex and browseVignettes('celldex') for documentation
#> loading from cache
```

## iDA decomposition

iDA can be run on the SummarizedExperiment object. iDA will normalize
(if not normalized already), identify highly variable genes (HVGs), and
then perform the iterative LDA on the normalized counts of the HVGs as
described in the algorithm above We can run iDA with the following:

``` r
set.seed(11)
start <- Sys.time()
BED <- iDA(BED, cluster.method = "louvain")
#> 2000 variable features found using rowVars().
#> iteration 1
#> concordance: 0.90962800628796
#> final concordance:
#> 0.99207246201568
end <- Sys.time()
```

    #> This SE dataset with 257 samples and 2000 HVGs ran in 2.82 minutes.

iDA produces both the LD’s (linear discriminants) and the cluster
assignments. The iDA embedding as well as cluster assignment labels will
be added as column data as seen in the colData().

``` r
head(colData(BED))
#> DataFrame with 6 rows and 14 columns
#>                                                    label.main  label.fine   label.ont       LD1       LD2       LD3       LD4
#>                                                   <character> <character> <character> <numeric> <numeric> <numeric> <numeric>
#> mature.neutrophil                                 Neutrophils Neutrophils  CL:0000775 -110.7855  -33.5680  -86.6750  29.02787
#> CD14.positive..CD16.negative.classical.monocyte     Monocytes   Monocytes  CL:0000576 -100.3858  -77.7118  -34.3153  47.63619
#> mature.neutrophil.1                               Neutrophils Neutrophils  CL:0000775 -107.6417  -37.0194  -82.9049  21.58429
#> megakaryocyte.erythroid.progenitor.cell                   HSC         MEP  CL:0000050  -16.6099  -23.5768   54.6583   3.54454
#> mature.neutrophil.2                               Neutrophils Neutrophils  CL:0000775 -110.4916  -36.9927  -88.6993  31.44123
#> CD14.positive..CD16.negative.classical.monocyte.1   Monocytes   Monocytes  CL:0000576 -100.3457  -80.0760  -37.2426  45.72353
#>                                                         LD5        LD6       LD7       LD8       LD9      LD10 iDA_clusters
#>                                                   <numeric>  <numeric> <numeric> <numeric> <numeric> <numeric>    <numeric>
#> mature.neutrophil                                   22.3722  -0.781729  -5.78720   3.22005  -22.2556   1.12689            1
#> CD14.positive..CD16.negative.classical.monocyte     23.1097 -23.204542 -56.92605 -32.85769  -13.3753 -23.64416            2
#> mature.neutrophil.1                                 25.0522  -3.126805  -7.84354   6.39983  -22.7654   3.85572            1
#> megakaryocyte.erythroid.progenitor.cell             60.2838  -7.024889 -36.61820   5.91998  -22.8179   2.64654            3
#> mature.neutrophil.2                                 27.3205  -0.768227 -14.93209  -2.67469  -18.1748  -4.40435            1
#> CD14.positive..CD16.negative.classical.monocyte.1   26.1455 -22.264428 -54.68060 -33.28467  -12.2926 -23.27457            2
```

There are three clustering methods available to use in the iDA
algorithm: louvain community detection, walktrap community detection, or
k-means clustering. These can be set in the cluster.method argument:

``` r
BED <- iDA(BED, cluster.method = c("louvain", "walktrap", "kmeans"))
```

You can also change the `nFeatures` parameter to use a different number
of HVGs for the embedding. The computation time does increase as
nFeatures increases since the embedding is computed with a single value
decomposition (SVD) of the gene x gene scatter matrix in each iteration.
Because of this though, the computation time scales relatively well for
increasing numbers of samples. the `nFeatures` parameter can be from the
default number (2000 HVGs) with the following:

``` r
BED <- iDA(BED, nFeatures = 3000)
```

We can visualize the cluster how many samples were assigned to each
cluster.

``` r
count <- table(colData(BED)[["iDA_clusters"]])
count.df <- data.frame(count)

plot <- ggplot(count.df, aes(x = Var1, y = Freq, fill = Var1))

plot + geom_bar(stat="identity") + 
        labs(title="Cluster Counts",
                     y="Count", x="Cluster") + 
        theme(legend.position="none")
```

<img src="README_files/figure-gfm/unnamed-chunk-21-1.png" width="100%" />

When using louvain clustering, iDA identifies 11 clusters, and therefore
10 LD dimensions to best separate these clusters.

### tSNE of the iDA dimensions

Since iDA is still a linear dimensionality reduction method and the
dimension will likely be larger than 2 (the dimension will always be k -
1), we may still want to pair this with a non-linear method such as tSNE
to visualize the clustering in a 2D plot. We can now treat these LDs
just as we would PCs as input to a non-linear dimensionality reduction
to reduce this to 2 dimensions for cluster visualization.

``` r
set.seed(11)
tsne_ida <- as.data.frame(colData(BED)) %>%
  select(contains("LD")) %>%
  Rtsne(pca = FALSE)

pBED <- tsne_ida$Y %>%
  as.data.frame() %>%
  ggplot(aes(x = V1, y = V2, color = as.factor(as.data.frame(colData(BED))$iDA_clusters))) + 
  geom_point() + 
  xlab("tSNE 1") +
  ylab("tSNE 2") +
  labs(color = "iDA Cluster")

pBED
```

<img src="README_files/figure-gfm/unnamed-chunk-22-1.png" width="100%" />

iDA produces an embedding which maximizes the separation between the
clusters found. We can then evaluate if these clusters are indicative of
known biological variation in the dataset. This bulk RNAseq dataset is
comprised of several cell types which were cell sorted and then
sequenced to produce homogenous cell type populations to be sequenced
together.

``` r
pBED_iDA <- tsne_ida$Y %>%
    as.data.frame() %>%
    ggplot(aes(x = V1, y = V2, color = as.data.frame(colData(BED))$label.main)) + 
    geom_point() + 
    xlab("tSNE 1") +
    ylab("tSNE 2") +
    labs(color = "iDA Cluster")

pBED_iDA
```

<img src="README_files/figure-gfm/unnamed-chunk-23-1.png" width="100%" />

We see some cell population’s expression is much more unique than
others. For example, the dendritic cells (DC - pea green) are extremely
unique and cluster together around tSNE coordinates (15,3). iDA also
captures the expression differences between keratinocytes, monocytes,
and Neutrophils. Interestingly, iDA finds trajectory-like differences in
the Hematopoietic stem cell (HSC) populations (sea-foam green) around
tSNE coordinates (-5, 2.5). These differences could be investigated for
potential sequencing batch effects or other sources of variation (cell
cycle, cell differentiation, maturity).

## Comparing to PCA decomposition

We may want to compare the iDA dimensions to what we would get from PCA.
Let’s compute the PCA dimension reduction for comparison.

``` r
BED_PCA <- scater::calculatePCA(BED)

set.seed(11)
tsne_pca <- BED_PCA %>%
    Rtsne(pca = FALSE)

pBED_pca <- tsne_pca$Y %>%
    as.data.frame() %>%
    ggplot(aes(x = V1, y = V2, 
            color = colData(BED)[["label.main"]])) + 
    geom_point() + 
    labs(color = "label.main")

pBED_pca
```

<img src="README_files/figure-gfm/unnamed-chunk-24-1.png" width="100%" />

We can see generally that the PCA embedding roughly has the same
architecture, but the cluster separation and discreteness from the other
clusters is not as defined as it is in the iDA projections. Therefore,
the top highly weighted genes in each iDA dimension are better
indicators of cluster (and therefore cell type) uniqueness. These genes
are likely to be candidates which have high fold changes in a
differential expression (DE) analysis.

# What to cite

The iDA package was published in Biostatistics in 2021. The full article
can be found here: <https://doi.org/10.1093/biostatistics/kxab030>

To cite the iDA package, please include the citation: Theresa A
Alexander, Rafael A Irizarry, Héctor Corrada Bravo, Capturing discrete
latent structures: choose LDs over PCs, Biostatistics, 2021;, kxab030,
<https://doi.org/10.1093/biostatistics/kxab030>

``` r
sessionInfo()
#> R version 4.1.3 (2022-03-10)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Big Sur 11.6.5
#> 
#> Matrix products: default
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] iDA_0.99.0                  gridExtra_2.3               scater_1.22.0               mclust_5.4.9               
#>  [5] igraph_1.2.11               irlba_2.3.5                 Matrix_1.4-1                scuttle_1.4.0              
#>  [9] scPipe_1.16.1               dplyr_1.0.8                 Rtsne_0.15                  ggplot2_3.3.5              
#> [13] celldex_1.4.0               scRNAseq_2.8.0              SingleCellExperiment_1.16.0 DESeq2_1.34.0              
#> [17] SummarizedExperiment_1.24.0 Biobase_2.54.0              GenomicRanges_1.46.1        GenomeInfoDb_1.30.1        
#> [21] IRanges_2.28.0              S4Vectors_0.32.4            BiocGenerics_0.40.0         MatrixGenerics_1.6.0       
#> [25] matrixStats_0.61.0         
#> 
#> loaded via a namespace (and not attached):
#>   [1] utf8_1.2.2                    RUnit_0.4.32                  tidyselect_1.1.2              RSQLite_2.2.11               
#>   [5] AnnotationDbi_1.56.2          grid_4.1.3                    BiocParallel_1.28.3           munsell_0.5.0                
#>   [9] ScaledMatrix_1.2.0            codetools_0.2-18              statmod_1.4.36                scran_1.22.1                 
#>  [13] withr_2.5.0                   colorspace_2.0-3              NetworkToolbox_1.4.2          filelock_1.0.2               
#>  [17] highr_0.9                     knitr_1.38                    rstudioapi_0.13               robustbase_0.93-9            
#>  [21] airway_1.14.0                 labeling_0.4.2                optparse_1.7.1                GenomeInfoDbData_1.2.7       
#>  [25] mnormt_2.0.2                  farver_2.1.0                  bit64_4.0.5                   rprojroot_2.0.2              
#>  [29] vctrs_0.3.8                   generics_0.1.2                xfun_0.30                     BiocFileCache_2.2.1          
#>  [33] R6_2.5.1                      ggbeeswarm_0.6.0              rsvd_1.0.5                    locfit_1.5-9.5               
#>  [37] AnnotationFilter_1.18.0       bitops_1.0-7                  cachem_1.0.6                  reshape_0.8.8                
#>  [41] DelayedArray_0.20.0           assertthat_0.2.1              promises_1.2.0.1              BiocIO_1.4.0                 
#>  [45] scales_1.1.1                  beeswarm_0.4.0                gtable_0.3.0                  beachmat_2.10.0              
#>  [49] org.Mm.eg.db_3.14.0           biocViews_1.62.1              ensembldb_2.18.4              rlang_1.0.2                  
#>  [53] genefilter_1.76.0             splines_4.1.3                 rtracklayer_1.54.0            lazyeval_0.2.2               
#>  [57] BiocManager_1.30.16           yaml_2.3.5                    GenomicFeatures_1.46.5        httpuv_1.6.5                 
#>  [61] RBGL_1.70.0                   tools_4.1.3                   psych_2.2.3                   ellipsis_0.3.2               
#>  [65] RColorBrewer_1.1-2            Rcpp_1.0.8.3                  plyr_1.8.7                    sparseMatrixStats_1.6.0      
#>  [69] progress_1.2.2                zlibbioc_1.40.0               purrr_0.3.4                   RCurl_1.98-1.6               
#>  [73] prettyunits_1.1.1             viridis_0.6.2                 ggrepel_0.9.1                 cluster_2.1.2                
#>  [77] magrittr_2.0.2                tmvnsim_1.0-2                 ProtGenerics_1.26.0           pkgload_1.2.4                
#>  [81] evaluate_0.15                 hms_1.1.1                     mime_0.12                     xtable_1.8-4                 
#>  [85] XML_3.99-0.9                  testthat_3.1.2                compiler_4.1.3                biomaRt_2.50.3               
#>  [89] tibble_3.1.6                  crayon_1.5.1                  htmltools_0.5.2               later_1.3.0                  
#>  [93] geneplotter_1.72.0            DBI_1.1.2                     ExperimentHub_2.2.1           dbplyr_2.1.1                 
#>  [97] MASS_7.3-56                   rappdirs_0.3.3                BiocStyle_2.22.0              getopt_1.20.3                
#> [101] brio_1.1.3                    cli_3.2.0                     parallel_4.1.3                metapod_1.2.0                
#> [105] pkgconfig_2.0.3               GenomicAlignments_1.30.0      xml2_1.3.3                    foreach_1.5.2                
#> [109] annotate_1.72.0               vipor_0.4.5                   dqrng_0.3.0                   stringdist_0.9.8             
#> [113] XVector_0.34.0                BiocCheck_1.30.0              stringr_1.4.0                 digest_0.6.29                
#> [117] graph_1.72.0                  Biostrings_2.62.0             rmarkdown_2.13                edgeR_3.36.0                 
#> [121] DelayedMatrixStats_1.16.0     restfulr_0.0.13               curl_4.3.2                    commonmark_1.8.0             
#> [125] shiny_1.7.1                   Rsamtools_2.10.0              rjson_0.2.21                  jsonlite_1.8.0               
#> [129] lifecycle_1.0.1               nlme_3.1-157                  BiocNeighbors_1.12.0          desc_1.4.1                   
#> [133] viridisLite_0.4.0             limma_3.50.1                  fansi_1.0.3                   pillar_1.7.0                 
#> [137] lattice_0.20-45               GGally_2.1.2                  Rhtslib_1.26.0                KEGGREST_1.34.0              
#> [141] fastmap_1.1.0                 httr_1.4.2                    DEoptimR_1.0-10               survival_3.3-1               
#> [145] interactiveDisplayBase_1.32.0 glue_1.6.2                    png_0.1-7                     iterators_1.0.14             
#> [149] bluster_1.4.0                 BiocVersion_3.14.0            bit_4.0.4                     stringi_1.7.6                
#> [153] blob_1.2.2                    org.Hs.eg.db_3.14.0           BiocSingular_1.10.0           AnnotationHub_3.2.2          
#> [157] memoise_2.0.1
```
