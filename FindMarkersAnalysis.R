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

library(BiocParallel)

#Gene list Creation

findAllMarkersA <- read.csv("/Users/nitin/Downloads/b_cells/FindAllMarkers/Group_A/FAMAnalysis_A.csv")
findAllMarkersB <- read.csv("/Users/nitin/Downloads/b_cells/FindAllMarkers/Group_B/FAMAnalysis_B.csv")
findAllMarkersC <- read.csv("/Users/nitin/Downloads/b_cells/FindAllMarkers/Group_C/FAMAnalysis_C.csv")
findAllMarkersD <- read.csv("/Users/nitin/Downloads/b_cells/FindAllMarkers/Group_D/FAMAnalysis_D.csv")
findAllMarkersE <- read.csv("/Users/nitin/Downloads/b_cells/FindAllMarkers/Group_E/FAMAnalysis_E.csv")

findAllMarkersA <- filter(findAllMarkersA, abs(findAllMarkersA$avg_log2FC) > 1.5)
findAllMarkersB <- filter(findAllMarkersB, abs(findAllMarkersB$avg_log2FC) > 1.5)
findAllMarkersC <- filter(findAllMarkersC, abs(findAllMarkersC$avg_log2FC) > 1.5)
findAllMarkersD <- filter(findAllMarkersD, abs(findAllMarkersD$avg_log2FC) > 1.5)
findAllMarkersE <- filter(findAllMarkersE, abs(findAllMarkersE$avg_log2FC) > 1.5)


findAllMarkersSigGenes <- rbind(findAllMarkersA, findAllMarkersB, findAllMarkersC, findAllMarkersD, findAllMarkersE)
findAllMarkersSigGenes <- findAllMarkersSigGenes$gene

findAllMarkersSigGenes <- unique(findAllMarkersSigGenes)
findAllMarkersSigGenes <- findAllMarkersSigGenes[!grepl("^[A-Z]+[0-9]+\\.[0-9]+$", findAllMarkersSigGenes)]

write.csv(findAllMarkersSigGenes, "/Users/nitin/Downloads/b_cells/FindAllMarkers/SigGenes.csv")



#b_cells Seurat Object Creation

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

#FindMarkers Legit Process

Idents(b_cells) <- b_cells$Age_group
FindAllMarkersResults <- FindAllMarkers(b_cells)
head(FindAllMarkersResults)

FMgeneNames <- FindAllMarkersResults$gene

#FindAllmarkers

Idents(b_cells) <- b_cells$Age_group
b_cells_A <- subset(b_cells, idents = "A")
Idents(b_cells_A) <- b_cells_A$Cluster_names

b_cells_B <- subset(b_cells, idents = "B")
Idents(b_cells_B) <- b_cells_B$Cluster_names

b_cells_C <- subset(b_cells, idents = "C")
Idents(b_cells_C) <- b_cells_C$Cluster_names

b_cells_D <- subset(b_cells, idents = "D")
Idents(b_cells_D) <- b_cells_D$Cluster_names

b_cells_E <- subset(b_cells, idents = "E")
Idents(b_cells_E) <- b_cells_E$Cluster_names
FM_group_A <- FindAllMarkers(b_cells_A)
FM_group_B <- FindAllMarkers(b_cells_B)
FM_group_C <- FindAllMarkers(b_cells_C)
FM_group_D <- FindAllMarkers(b_cells_D)
FM_group_E <- FindAllMarkers(b_cells_E)

FM_group_A <- FM_group_A %>%
  as.data.frame() %>% 
  arrange(desc(abs(avg_log2FC))) %>%
  filter(pct.2 > 0) %>%
  filter(pct.1 > 0)

FM_group_B <- FM_group_B %>%
  as.data.frame() %>% 
  arrange(desc(abs(avg_log2FC))) %>%
  filter(pct.2 > 0) %>%
  filter(pct.1 > 0)

FM_group_C <- FM_group_C %>%
  as.data.frame() %>% 
  arrange(desc(abs(avg_log2FC))) %>%
  filter(pct.2 > 0) %>%
  filter(pct.1 > 0)

FM_group_D <- FM_group_D %>%
  as.data.frame() %>% 
  arrange(desc(abs(avg_log2FC))) %>%
  filter(pct.2 > 0) %>%
  filter(pct.1 > 0)

FM_group_E <- FM_group_E %>%
  as.data.frame() %>% 
  arrange(desc(abs(avg_log2FC))) %>%
  filter(pct.2 > 0) %>%
  filter(pct.1 > 0)

head(FM_group_E)

write.csv(FM_group_A, "C:/Users/nitin/Downloads/b_cells/FindAllMarkers/Group_A/FAMAnalysis_A.csv", quote = F)
write.csv(FM_group_B, "C:/Users/nitin/Downloads/b_cells/FindAllMarkers/Group_B/FAMAnalysis_B.csv", quote = F)
write.csv(FM_group_C, "C:/Users/nitin/Downloads/b_cells/FindAllMarkers/Group_C/FAMAnalysis_C.csv", quote = F)
write.csv(FM_group_D, "C:/Users/nitin/Downloads/b_cells/FindAllMarkers/Group_D/FAMAnalysis_D.csv", quote = F)
write.csv(FM_group_E, "C:/Users/nitin/Downloads/b_cells/FindAllMarkers/Group_E/FAMAnalysis_E.csv", quote = F)

FM_sig_group_E <- FM_group_E %>%
  filter(abs(avg_log2FC) > 4)

write.csv(FM_sig_group_E, "C:/Users/nitin/Downloads/b_cells/FindAllMarkers/Group_E/FAM_sig_genes_E.csv", quote = F)

#Switched Memory Cells

FM_Switched_b_cells <- subset(b_cells, Cluster_names == "Switched memory")

Idents(FM_Switched_b_cells) <- "Age_group"
FM_Switched_data <- FindMarkers(FM_Switched_b_cells, ident.1 = "A", ident.2 = "E")

FM_Switched_sorted_data <- FM_Switched_data %>%
  as.data.frame() %>% 
  arrange(desc(abs(avg_log2FC)))

head(FM_Switched_sorted_data)

FM_Switched_t20_genes <- FM_Switched_sorted_data %>%
  head(n = 20)

write.csv(FM_Switched_t20_genes, "FM_Switched_sig_genes.csv")


#Non_Switched Memory Cells

FM_Non_Switched_b_cells <- subset(b_cells, Cluster_names == "Non-switched memory")

Idents(FM_Non_Switched_b_cells) <- "Age_group"
FM_Non_Switched_b_cells <- FindMarkers(FM_Non_Switched_b_cells, ident.1 = "A", ident.2 = "E")

FM_Non_Switched_sorted_data <- FM_Non_Switched_b_cells %>%
  as.data.frame() %>% 
  arrange(desc(abs(avg_log2FC)))

head(FM_Non_Switched_sorted_data)

FM_Non_Switched_t20_genes <- FM_Non_Switched_sorted_data %>%
  head(n = 20)

write.csv(FM_Non_Switched_t20_genes, "FM_Non_Switched_sig_genes.csv")

#Plasma

FM_Plasma_b_cells <- subset(FM_b_cells, Cluster_names == "Plasma cells")

Idents(FM_Plasma_b_cells) <- "Age_group"
FM_Plasma_b_cells <- FindMarkers(FM_Plasma_b_cells, ident.1 = "A", ident.2 = "E")

FM_Plasma_sorted_data <- FM_Plasma_b_cells %>%
  as.data.frame() %>% 
  arrange(desc(abs(avg_log2FC)))

head(FM_Plasma_sorted_data)

FM_Plasma_t20_genes <- FM_Plasma_sorted_data %>%
  head(n = 20)

write.csv(FM_Plasma_t20_genes, "FM_Plasma_sig_genes.csv")


