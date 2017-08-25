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
# the non-normalized values within a cell

# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@data.info, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")

# We filter out cells that have unique gene counts over 4000 or less than
# 0 Note that low.thresholds and high.thresholds are used to define a
# 'gate' -Inf and Inf should be used if you don't want a lower or upper
# threshold.
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene"), 
                    low.thresholds = 0, high.thresholds = 4000)

GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")

# After removing unwanted cells from the dataset, the next step is to normalize the data.
# By default, we employ a global-scaling normalization method “LogNormalize” that normalizes
# the gene expression measurements for each cell by the total expression, multiplies this by a
# scale factor (10,000 by default), and log-transforms the result.

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# Seurat calculates highly variable genes and focuses on these for downstream analysis.
# FindVariableGenes calculates the average expression and dispersion for each gene,
# places these genes into bins, and then calculates a z-score for dispersion within each bin.
# This helps control for the relationship between variability and average expression.
# This function is unchanged from (Macosko et al.), but new methods for variable gene expression 
# identification are coming soon. We suggest that users set these parameters to mark visual outliers 
# on the dispersion plot, but the exact parameter settings may vary based on the data type, heterogeneity 
# in the sample, and normalization strategy.
# Change cuttoffs looking at the plot :)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 4, y.cutoff = 0.4)
# ~2,000 variable genes

length(x = pbmc@var.genes)
# Your single cell dataset likely contains ‘uninteresting’ sources of variation.
# This could include not only technical noise, but batch effects, or even biological sources
# of variation (cell cycle stage). As suggested in Buettner et al, NBT, 2015, regressing these signals
# out of the analysis can improve downstream dimensionality reduction and clustering.

pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))

# You can use either PCA or ICA for dimensional reduction
#########################################################
# PCA
# pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:20,
#              genes.print = 5, pcs.compute = 20)
# 
# Using this plot find optimal number of components
#
# PCElbowPlot(object = pbmc)
# 
# pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:7,
#                              genes.print = 5, pcs.compute = 7)
# 
# Find clusters using SNN-Cliq alg with PCA reduction (resolution param - find optimal yourself :) 
#
# pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:7,
#                     resolution = 0.6, print.output = 0, save.SNN = TRUE)
# 
# PrintFindClustersParams(object = pbmc)
#
# Running TSNE to reduce dimension (for visualizing clusters)
#
# pbmc <- RunTSNE(object = pbmc, dims.use = 1:7, do.fast = TRUE)
#########################################################
#ICA
pbmc <- RunICA(object = pbmc, ic.genes = pbmc@var.genes, ics.compute = 8)

# Number of components - find optimal yourself :)
pbmc <- FindClusters(object = pbmc, reduction.type = "ica", dims.use = 1:8,
                     resolution = 0.5, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = pbmc)

pbmc <- RunTSNE(object = pbmc, reduction.use = "ica", dims.use = 1:8, do.fast = TRUE)
#########################################################

# Vizualize clusters. Note that you can set do.label=T to help label individual clusters
TSNEPlot(object = pbmc)

save(pbmc, file = "~/R/single_cell_rna_seq")

# Find markers for every cluster compared to all remaining cells, report
# only the positive ones

pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE)

# Find top 5 genes of each cluster and save

top_genes <- pbmc.markers %>% group_by(cluster) %>% top_n(5, avg_diff)
write.csv(top_genes, file = "./gen.csv")

# Vizualize top-1 genes of each cluster on whole dataset

FeaturePlot(object = pbmc, features.plot = c("LEF1", "S100A8", "IGKC", "IL32", 
 "NKG7", "HLA-DQA1", "NKG7", "FCGR3A", "IRF7", "HIST1H2AC"), cols.use = c("grey", "blue"),
 reduction.use = "tsne")
