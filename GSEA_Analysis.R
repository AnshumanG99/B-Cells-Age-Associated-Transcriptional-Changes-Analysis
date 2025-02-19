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

#GSEA Analysis Loading

library(clusterProfiler)
library(org.Hs.eg.db)


count_matrix <- readRDS("/Users/nitin/downloads/b_cells/b_cells_rna.rds")
b_cells <- CreateSeuratObject(counts = count_matrix)

umap_data <- read.csv("/Users/nitin/Downloads/b_cells/b_cells_umap_label.csv")
head(umap_data)
rownames(umap_data) <- umap_data$cell_id

b_cells[["umap"]] <- CreateDimReducObject(
  embeddings = as.matrix(umap_data[, c("UMAP_1", "UMAP_2")]),
  key = "UMAP_"
)

metadata <- read.csv("/Users/nitin/Downloads/b_cells/b_cells_metadata.csv", row.names = 1)
b_cells <- AddMetaData(b_cells, metadata = metadata)

b_cells <- NormalizeData(b_cells)
b_cells <- ScaleData(b_cells)

#GSEA

#Set the active ident as the age group
Idents(b_cells) <- b_cells$Age_group

#Choose the age groups you want to compare between and save the resulting output as a .csv
b_cells_AvsB <- FindMarkers(b_cells, ident.1 = "A", ident.2 = "B")
write.csv(b_cells_AvsB, "/Users/nitin/Downloads/b_cells/GSEA_Analysis/b_cells_AvsB.csv", quote = F) 


#OPTION2: If you are going to filter the output in R then you will probably (?) need to first convert the output to a dataframe: 
#Then you can filter etc as needed
b_cells_AvsB <- as.data.frame(b_cells_AvsB)

b_cells_AvsB_filtered <- b_cells_AvsB %>%
  filter(!grepl("^[A-Z]+[0-9]+\\.[0-9]+$", rownames(b_cells_AvsB)))

#write.csv(b_cells_AvsB_filtered, "/Users/nitin/Downloads/b_cells/GSEA_Analysis/b_cells_AvsB_filteringtest.csv", quote = F)

# Convert rownames to a column
b_cells_AvsB <- b_cells_AvsB_filtered
b_cells_AvsB <- as.data.frame(b_cells_AvsB)
b_cells_AvsB$gene <- rownames(b_cells_AvsB)

# Move 'gene' to the first column without using select
b_cells_AvsB <- b_cells_AvsB[, c("gene", setdiff(names(b_cells_AvsB), "gene"))]

#write.csv(b_cells_AvsB, "/Users/nitin/Downloads/b_cells/GSEA_Analysis/b_cells_AvsB.csv")


#Change variable to res for convenience
res <- read.csv("/Users/nitin/Downloads/b_cells/GSEA_Analysis/b_cells_AvsB.csv")

#For this analysis you need to input an ordered set of genes. This orders them in descending order based on log2FC
res <- res[order(-res$avg_log2FC), ]

#Create your list of genes and expression values that you will input into gseGO
genelist <- res$avg_log2FC
names(genelist) <- res$gene  #Remember you will have needed to rename the gene column in the FindMarkers output as "gene"
genelist

#Run gseGO
gse <- gseGO(genelist, ont = "BP", key = "SYMBOL", OrgDb = "org.Hs.eg.db", eps = 1e-300, minGSSize = 1)
as.data.frame(gse)

#Save resulting GSEA analysis
#write.csv(gse, "/Users/nitin/Downloads/b_cells/GSEA_Analysis/GSEA_allcells_AvsB.csv", quote = F)

#Reorder the gse output based on the NES score
gse_ordered  <- gse[order(-gse$NES), ]
gse_ordered

#Plot the NES scores, be sure to change the title of the plot
ggplot(gse_ordered, aes(reorder(Description, NES), NES)) +
  geom_col(aes(fill=p.adjust<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="All cells, Age group A vs B") + 
  theme_minimal()

#Here is where you can input the genes from whichever pathways we decide to highlight
immuneresponse <- c('IGHV5-51',
                    'TNFSF13',
                    'IGHV3-48',
                    'IGLC2',
                    'IGHV1-24',
                    'GAPT',
                    'IGHD',
                    'C1RL',
                    'IGHV3-21',
                    'IGHV3-20',
                    'IGHV3-11',
                    'IGHV2-70',
                    'IGHV4-28',
                    'IGHV1-69D',
                    'IGHV4-31',
                    'TREX1',
                    'IGHV2-26',
                    'IGHV6-1',
                    'IGHV1-3',
                    'IGHV3-30',
                    'IGLC1',
                    'IGLC7',
                    'IGHV2-5',
                    'IGHV5-10-1',
                    'IGHV4-4')

#Generate the first heatmap showing the gene expression across the different age groups
avgExp = AggregateExpression(object = b_cells, return.seurat = TRUE, group.by = "Age_group")
DoHeatmap(avgExp, features = immuneresponse, size = 3)

#Subset the desired age group 
#AgeGroupE <- subset(b_cells, ident = "E")

#Generate a heatmap showing the gene expression across the different cell types within the desired age group
avgExp = AggregateExpression(object = AgeGroupE, return.seurat = TRUE, group.by = "Cluster_names")
DoHeatmap(avgExp, features = immuneresponse, size = 3)