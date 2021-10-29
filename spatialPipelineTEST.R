###Test spatial pipeline###
#  following: https://lmweber.org/OSTA-book/analysis-steps.html#load-data

## TEST TEST


## install and load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version='devel')

if (!requireNamespace("scater", quietly = TRUE))
  BiocManager::install("scater")
if (!requireNamespace("SpatialExperiment", quietly = TRUE))
  BiocManager::install("SpatialExperiment")
if (!requireNamespace("spatialLIBD", quietly = TRUE))
  BiocManager::install("spatialLIBD")
if (!requireNamespace("STexampleData", quietly = TRUE))
  BiocManager::install("STexampleData")
if (!requireNamespace("ggspavis", quietly = TRUE))
  BiocManager::install("ggspavis")
if (!requireNamespace("STexampleData", quietly = TRUE))
  BiocManager::install("STexampleData", version = "devel")
if (!requireNamespace("scater", quietly = TRUE))
  BiocManager::install("scater")
if (!requireNamespace("scran", quietly = TRUE))
  BiocManager::install("scran")
if (!requireNamespace("pheatmap", quietly = TRUE))
  BiocManager::install("pheatmap")

library(scater)
library(spatialLIBD)
library(SpatialExperiment) #specialized object class for spatial data
library(STexampleData)
library(ggspavis)
library(scran)
library(pheatmap)


## load data
spe <- STexampleData::Visium_humanDLPFC_3_13()

## plot data 
#  Error: 'data' must be uniquely names but has duplicate columns
plotSpots(spe)

## quality control - calculate
#  keep only spots over tissue (we dont care about background) - using variable in_tissue
spe <- spe[, spatialData(spe)$in_tissue == 1]
dim(spe)
#  identify mitochrondrial genes - indicated by mt
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)
#  calculate QC per spot, store in colData
spe <- addPerCellQC(spe, subsets = list(mito=is_mito))
head(colData(spe))

## quality control - threshold

## library size
#  visualize library size
#    for example, look for spike at low library size
hist(colData(spe)$sum, breaks = 20)
#  visualize library size vs num cells per spot
#    look at arbitrary threshold of 500 for library size
#      does this select for biologically consistent group of spots??
plotQC(spe, type="scatter",metric_x="cell_count",metric_y="sum",threshold_y=500)
#  threshold
qc_lib_size <- colData(spe)$sum < 500
colData(spe)$qc_lib_size <- qc_lib_size
#   visualize number of spots below threshold, check for spatial patter
#      we don't want to see any pattern/biologically informative
#   Error: 'data' must be uniquely named but has duplicate columns
plotQC(spe, type = "spots", discard = "qc_lib_size")
#   visualize ground truth (manually annotated)
#   Error: 'data' must be uniquely named but has duplicate columns
plotSpots(spe, annotate="ground_truth",palette="libd_layer_colors")

##  number of expressed genes
#   visualize number of expressed genes
hist(colData(spe)$detected, breaks=20)
#   visualize num genes vs num cells
#   use this to select threshold
plotQC(spe, type = "scatter", metric_x = "cell_count", metric_y = "detected", threshold_y = 250)
#   threshold
qc_detected <- colData(spe)$detected < 250
colData(spe)$qc_detected <- qc_detected
#   visualize discarded spots
plotQc(spe, type="spots",discard="qc_detected")

## proportion of mitochondrial reads
## **indicates cell damage**
#  visualize distribution
hist(colData(spe)$subsets_mito_percent, breaks=20)
plotQC(spe,type="scatter",metric_x = "cell_count",
      metric_y = "subsets_mito_percent", threshold_y = 30)
#  threshold
qc_mito <- colData(spe)$subsets_mito_percent > 30
table(qc_mito)
colData(spe)$qc_mito <- qc_mito
# check discarded spots
plotQC(spe, type="spots", discard="qc_mito")

## number of cells per spot
#  **outliers indicate problems with cell segmentation
#  **high cell counts and low genes indicates failed experiments
#  visualize distribution
hist(colData(spe)$cell_count, breaks =20)
tbl_cells_per_spot <- table(colData(spe)$cell_count)
plotQC(spe, type="scatter", metric_x = "cell_count",
      metric_y = "detected", threshold_x = 12)
#  threshold
qc_cell_count <- colData(spe)$cell_count > 12
table(qc_cell_count)
colData(spe)$qc_cell_count <- qc_cell_count

## remove combined sets of low quality 
apply(cbind(qc_lib_size, qc_detected, qc_mito, qc_cell_count), 2, sum)
discard <- qc_lib_size | qc_detected | qc_mito | qc_cell_count
table(discard)
colData(spe)$discard <- discard
plotQC(spe,type="spots",discard="discard")
spe <- spe[, !colData(spe)$discard]
dim(spe)

## Normalize
#  treat each spot as one cell
#  log transform normalized counts
set.seed(123)
qclus <- quickCluster(spe)
table(qclus)
spe <- computeSumFactors(spe,cluster=qclus)
summary(sizeFactors(spe))
hist(sizeFactors(spe),breaks=20)

spe <- logNormCounts(spe)
assayNames(spe)
dim(counts(spe))
dim(logcounts(spe))


## Feature Selection

#  HVGs = highly variable genes
#  first mitochondrial since highly expressed and not of interest
spe <- spe[!is_mito, ]
dim(spe)
#  use scran for list of HVGS
dec <- modelGeneVar(spe) # fit mean var 
fit <- metadata(dec)
plot(fit$mean, fit$var, xlab = "mean of log-expression",
     ylab = "variance of log-expression")
curve(fit$trend(x),col="dodgerblue",add=TRUE,lwd=2)
top_hvgs <- getTopHVGs(dec, prop=0.1)
length(top_hvgs)

#  SVG = spatially variable genes
#  **not in tutorial


## Dimensionality Reduction
#  PCA
set.seed(123)
spe <- runPCA(spe,subset_row=top_hvgs)
plotDimRed(spe, type = "PCA")

# UMAP
set.seed(123)
spe <- runUMAP(spe,dimred="PCA")
reducedDimNames(spe)
dim(reducedDim(spe,"UMAP"))
colnames(reducedDim(spe,"UMAP")) <- paste0("UMAP", 1:2)
plotDimRed(spe, type = "UMAP")

## Clustering
#  cluster on HVG
set.seed(123)
k <- 10
g <- buildSNNGraph(spe,k=k, use.dimred="PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
table(clus)
colLabels(spe) <- factor(clus)

#  visualize
plotSpots(spe, annotate = "label", 
          palette = "libd_layer_colors") #clusters
plotSpots(spe, annotate = "ground_truth", 
          palette = "libd_layer_colors") #ground truth
plotDimRed(spe, type = "PCA", 
           annotate = "label", palette = "libd_layer_colors") #PCA
plotDimRed(spe, type = "UMAP", 
           annotate = "label", palette = "libd_layer_colors") #UMAP

## Marker Genes
#  test for marker genes
rownames(spe) <- rowData(spe)$gene_name
markers <- findMarkers(spe, test = "binom", direction = "up")
markers
#  plot log FC 
interesting <- markers[[1]]
best_set <- interesting[interesting$Top <= 5, ]
logFCs <- getMarkerEffects(best_set)
pheatmap(logFCs, breaks = seq(-5, 5, length.out = 101))

top_genes <- head(rownames(interesting))
plotExpression(spe, x = "label", features = top_genes)

## Spot level deconvolution
#  number of cells per spot
plotQC(spe, type = "bar", metric_x = "cell_count") + 
  xlab("number of cells") + ggtitle("Number of cells per spot")
# then what??





