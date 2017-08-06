library(Seurat)
library(dplyr)
library(Matrix)

# Load the PBMC dataset (33694 x 4287)
pbmc.data <- Read10X(data.dir = "./data/cr-team-2")

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes (16640 genes across 4287 samples)
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, project = "10X_PBMC")

# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.  NOTE: You must have the Matrix package loaded to
# calculate the percent.mito values.

# AddMetaData adds columns to object@data.info, and is a great place to
# stash QC stats

# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@data.info, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")

# We filter out cells that have unique gene counts over 2,500 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate' -Inf and Inf should be used if you don't want a lower or upper
# threshold.
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene"), 
                    low.thresholds = 0, high.thresholds = 4000)

GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")
 

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 4, y.cutoff = 0.4)
length(x = pbmc@var.genes)

pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))

# PCA
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:7, 
               genes.print = 5, pcs.compute = 7)

# PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
# VizPCA(object = pbmc, pcs.use = 1:2)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)

PCElbowPlot(object = pbmc)

#library("factoextra")
#get_clust_tendency(df.stand, n = 10,
                   #gradient = list(low = "steelblue", high = "white"))

# save.SNN = T saves the SNN so that the clustering algorithm can be rerun
# using the same graph but with a different resolution value (see docs for
# full details)
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:7, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = pbmc)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:7, do.fast = TRUE)

# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = pbmc)
save(pbmc, file = "~/R/single_cell_rna_seq")

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE)
pbmc.markers %>% group_by(cluster) %>% top_n(5, avg_diff)

pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_diff)

