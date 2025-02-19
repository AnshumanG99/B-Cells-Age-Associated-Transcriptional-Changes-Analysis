R.version.string
install.packages("pheatmap")
library(installr)
updateR()

#Load and install any necessary packages. You need to individually install each of them the first time
#All subsequent uses just require you to load each of the packages to the workspace  
install.packages("readxl") #use this to install the packages
devtools::install_github('immunogenomics/presto') #some packages require a slightly different installation

setwd("C://Users/nitin/Downloads/b_cells")
library(readxl) #this laods the installed package
library(Matrix)
library(dplyr)
library(devtools)
library(SeuratObject)
library(Seurat)
library(SeuratData)
library(patchwork)
library(cowplot)
library(RCurl)
library(purrr)
library(tibble)
library(ggplot2)
library(hdf5r)
library(devtools)
library(presto)
library(KernSmooth)
library(survival)
library(Nebulosa)
library(multtest)
library(metap)
library(pheatmap)

#DE Analysis loading
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)

#Trajectory Mapping Loading

library(BiocParallel)
library("slingshot")
library("Polychrome")
library("ggbeeswarm")
library("ggthemes")
library("SingleCellExperiment")
library("SummarizedExperiment")
library(S4Vectors) 
library(grDevices)
library(RColorBrewer)
library("DelayedMatrixStats")

library(BiocParallel)
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library("mgcv")
library("magrittr")
library(Biobase)
library("pbapply")
library("igraph")
library("princurve")
library("TrajectoryUtils")
library("viridis")
library("matrixStats")
library("MASS")


#B Cells Creation

count_matrix <- readRDS("/Users/nitin/downloads/b_cells/b_cells_rna.rds")
b_cells <- CreateSeuratObject(counts = count_matrix)

umap_data <- read.csv("/Users/nitin/Downloads/b_cells/b_cells_umap_label.csv")
rownames(umap_data) <- umap_data$cell_id
b_cells[["umap"]] <- CreateDimReducObject(
  embeddings = as.matrix(umap_data[, c("UMAP_1", "UMAP_2")]),
  key = "UMAP_"
) 

metadata <- read.csv("/Users/nitin/Downloads/b_cells/b_cells_metadata.csv", row.names = 1)  
b_cells <- AddMetaData(b_cells, metadata = metadata) 

b_cells <- NormalizeData(b_cells)
b_cells <- ScaleData(b_cells)


Idents(b_cells) <- b_cells$Cluster_names 
DefaultAssay(b_cells) <- "RNA"
b_cells[["RNA"]] <- JoinLayers(b_cells[["RNA"]])




#Trajectory Mapping Script


Idents(b_cells) <- b_cells$Age_group
b_cells_A <- subset(b_cells, ident = "A")

# Create a single cell experiment object for the input for the slingshot analysis
sce_A <- as.SingleCellExperiment(b_cells_A)
sce_A$ident <- sce_A$Cluster_names

# Generate your slingshot object
slingshot_A <- slingshot(sce_A,  clusterLabels = sce_A$ident, reducedDim = "UMAP",
                         start.clus = "Transitional")

# To view the number of lineages mapped you can either open the resulting slingshot_x object from the environment
# or type slingshot_A@colData$ and view all the lineages (slingPseudotime_x) in the dropdown menu

# Generate a summary of a created pseudotime lineage

summary(slingshot_A$slingPseudotime_1)

# print and/or export the list of lineages and curves
SlingshotDataSet(slingshot_A)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)


# Calculate the average pseudotime across all trajectories (slingPseudotime_1, slingPseudotime_2 ... slingPseudotime_5)
pseudotimes_A <- slingPseudotime(slingshot_A)

pseudtimes_df <- as.data.frame(pseudotimes_A)
average_pseudotime <- rowMeans(pseudtimes_df, na.rm=TRUE)

# Replace the slingPseudotime_1 in plotcol (below) with the average pseudotime values
# Dictate the colors to be used for plotting
plotcol <- colors[cut(average_pseudotime, breaks=100)]

# Show different pathways
# plotcol[is.na(slingshot_A$slingPseudotime_5)] <- "purple"

# Plot the UMAP, color coded by the pseudotime data, with lineages traced in black lines 
# I don't have it included here but I would recommend adding a title to each plot
plot(reducedDims(slingshot_A)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(slingshot_A), lwd=2, col='black')
title(main = "Age Group A Trajectory Pathways")

# Label the centroids of clusters to add the cluster names to the plot
cluster_centroids <- aggregate(reducedDims(slingshot_A)$UMAP, 
                               by = list(cluster = colData(sce_A)$Cluster_names), 
                               FUN = mean)

text(cluster_centroids[, 2], cluster_centroids[, 3], 
     labels = cluster_centroids$cluster, 
     pos = 3, col = 'black', cex = 0.8)

#Lineage information
lineage_list <- slingLineages(slingshot_A)
lineage_info <- sapply(seq_along(lineage_list), function(i) {
  paste0("Lineage ", i, ": ", paste(lineage_list[[i]], collapse = " -> "))
})
legend("bottomright", legend = lineage_info, bty = "n", cex = 0.45, title = "Lineage Data")




#Subset age group of interest
b_cells_A <- subset(b_cells, ident = "A")

#Create single cell experiment
sce_A <- as.SingleCellExperiment(b_cells_A)

# Set seed for reproducibility, the code chooses a random starting point within the data
# To make sure it chooses the same starting point, and therefore arrives at the same output, we use set.seed
# set.seed(5)
# library(BiocParallel)
# BPPARAM <- BiocParallel::SnowParam(workers = 16)

# This is the step that is computationally heavy, try reducing the nGenes variable to only 100 genes and see if that helps
# icMat <- evaluateK(counts = b_cells_A@assays$RNA@layers$counts, sds = SlingshotDataSet(slingshot_A), k = 3:10, nGenes = 100, verbose = T, parallel = TRUE, BPPARAM = BPPARAM)
# write.csv(icMat, "C:/Users/nitin/Downloads/b_cells/evaluateKoutput.csv")

# Assign your pseudotime and cellWeight variables from the previous slingshot analysis
pseudotime <- slingPseudotime(slingshot_A, na = FALSE)
cellWeights <- slingCurveWeights(slingshot_A)

# Assign subset_genes as the key genes you identified in FindAllMarkers analysis
subset_genes <- read.csv("/Users/nitin/Downloads/b_cells/FindAllMarkers/SigGenes.csv")
subset_genes <- subset_genes$x
subset_genes <- subset_genes[!grepl("^[A-Z]+[0-9]+\\.[0-9]+$", subset_genes)]


gene_names <- b_cells_A@assays$RNA@features@.Data[, 0]
gene_names <- rownames(gene_names)
rownames(b_cells_A@assays$RNA@layers$counts) <- gene_names

cell_names <- b_cells_A@assays$RNA@cells@.Data[,0]
cell_names <- rownames(cell_names)
colnames(b_cells_A@assays$RNA@layers$counts) <- cell_names

set.seed(5)
library(BiocParallel)
BPPARAM <- BiocParallel::SnowParam(workers = 8)
bpstart(BPPARAM)

filtered_countsA <- b_cells_A@assays$RNA@layers$counts[subset_genes, ]

sce_A <- fitGAM(counts = filtered_counts, pseudotime = pseudotime, cellWeights = cellWeights,
                 nknots = 6, genes = subset_genes, verbose = TRUE, parallel = TRUE, BPPARAM = BPPARAM)




# Code for subsetting single lineage and running fitgam
# Not needed since we can run it on the full object as long as we subset genes before
# 
# pseudotime_Lineage1 <- pseudotime[1:23652, 0:1, drop = FALSE]
# cellWeights_Lineage1 <- cellWeights[1:23652, 0:1, drop = FALSE]
# 
# cellWeights_Lineage1 <- cellWeights[cellWeights[, "Lineage1"] > 0, "Lineage1", drop = FALSE]
# 
# valid_cells <- rownames(cellWeights_Lineage1)
# pseudotime_Lineage1 <- pseudotime_Lineage1[valid_cells, , drop = FALSE]
# 
# # Subset counts matrix to include only the valid genes
# filtered_counts <- b_cells_A@assays$RNA@layers$counts[subset_genes, valid_cells]
# 
# # Verify dimensions
# dim(filtered_counts)
# 
# # Fit a model to each age group
# sce_A <- fitGAM(counts = filtered_counts, pseudotime = pseudotime_Lineage1, cellWeights = cellWeights_Lineage1,
#                 nknots = 6, genes = subset_genes, verbose = TRUE, parallel = TRUE, BPPARAM = BPPARAM)


# Determine genes that vary over pseudotime
assoResA <- associationTest(sce_A)
assoResB <- associationTest(sce_B)

head(assoResA)
head(assoResB)


# Test for differences in expression patterns across groups and lineages
# Not currently working
condition_test_A <- conditionTest(sce_A)
condition_test_B <- conditionTest(sce_B)

# condition_test <- conditionTest(list(sce_A, sce_B))

# Run lineage-specific pattern test
pattern_testA <- patternTest(sce_A)
head(pattern_testA)

pattern_testB <- patternTest(sce_B)
head(pattern_testB)


# Analyze the expession of significant genes throughout lineages

startRes <- startVsEndTest(sce_A)
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce_A)[oStart[1]]
plotSmoothers(sce_A, counts = filtered_countsA, gene = sigGeneStart)

plotGeneCount(sce_A$crv, counts, gene = sigGeneStart)
