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
b_cells <- CreateSeuratObject(counts = count_matrix) #Create your Seurat object and populate the counts with the loaded count matrix data

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


sigGenesTradeSeq <- read.csv("tradeseq_filteredgenes.csv")[,2]
sigGenesfindAllMarkers <- read.csv("C:/Users/nitin/Downloads/b_cells/FindAllMarkers/FinalResults.csv")$gene

sigGenesFinal <- unique(sigGenesTradeSeq[sigGenesTradeSeq %in% sigGenesfindAllMarkers])


VlnPlot(object = b_cells, features = "FOS", log = TRUE, group.by = "Age_group", pt.size = 0)
