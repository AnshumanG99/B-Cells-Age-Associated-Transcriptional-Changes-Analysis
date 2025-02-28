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

#Create a folder on your desktop to serve as the working directory. Ensure all the needed data is within this folder.
count_matrix <- readRDS("/Users/nitin/downloads/b_cells/b_cells_rna.rds") #This loads the sequencing data from the working directory, change the directory to match yours. The final .rds is the file name
b_cells <- CreateSeuratObject(counts = count_matrix) #Create your Seurat object and populate the counts with the loaded count matrix data

#Here, we will load coordinate data (x,y) that maps the Seurat object onto 2 dimensional space, mapping cells with similar gene expression profiles closer together and cells with dissimilar profiles farther away
umap_data <- read.csv("/Users/nitin/Downloads/b_cells/b_cells_umap_label.csv") #Load the UMAP data. This was provided as part of the public data set, typically you would have to generate your own UMAP. 
head(umap_data) #Check that the umap_data is in the correct format. The first column should be the cell identifier, the second column should be UMAP_1 (the x coordinates), and the third column should be UMAP_2 (the y coordinates)
rownames(umap_data) <- umap_data$cell_id  #Set the row names of the matrix to cell_id
#Create the dimensionally reduced object that will be used to plot the UMAP
b_cells[["umap"]] <- CreateDimReducObject(
  embeddings = as.matrix(umap_data[, c("UMAP_1", "UMAP_2")]),
  key = "UMAP_"
) 
DimPlot(b_cells, reduction = "umap")  #This plots the UMAP and a plot should appear in the "plots" window to the right

metadata <- read.csv("/Users/nitin/Downloads/b_cells/b_cells_metadata.csv", row.names = 1)  #load the meta data for the seurat object
head(metadata) #check that it contains the correct data in the correct form (I already checked, it does)
b_cells <- AddMetaData(b_cells, metadata = metadata) #Add the meta data to the seurat object

#View the meta data of the seurat object
b_cells@meta.data %>%
  View()

DimPlot(b_cells, reduction = "umap", group.by = "Age_group") #you can now plot the UMAP for the seurat object and color the UMAP by any of the variables listed in the meta data, e.g. Donor_id, Age_group, Sex, Age, etc

#Normalize the data
b_cells <- NormalizeData(b_cells)

#Identify the most variable genes in the merged object
b_cells <- FindVariableFeatures(b_cells,
                                selection.method = "vst",
                                nfeatures = 2000,
                                verbose = FALSE)

# Scale the counts
b_cells <- ScaleData(b_cells)

#Sanity check -- FeaturePlot shows the per-cell expression levels of selected genes. Below, PTPRC should be ubiquitously highly expressed whereas CD19 should have a more varied expression.
FeaturePlot(b_cells, features = c("PTPRC"), reduction = "umap")
FeaturePlot(b_cells, features = c("CD19"), reduction = "umap")

#Currently, the active identity being used is just the project name. We want to change this such that we are querying cell types instead
b_cells@active.ident #Check the active ident. It should return a list of the cell barcodes with the active identity (SeuratProject)
Idents(b_cells) <- b_cells$Cluster_names #This sets the active ident. We could choose any of the column headers from the metadata
b_cells@active.ident #Check that it is reset, it should now return the cell barcodes with the names of the clusters (eg. cell types)
DimPlot(b_cells, reduction = "umap")  #The UMAP is now color coded by cell type (Cluster_names)

# Set the dafault assay (either sct or rna, usually sct is better)
DefaultAssay(b_cells) <- "RNA"
b_cells[["RNA"]] <- JoinLayers(b_cells[["RNA"]])

# Find all gene markers for each cluster and save as a .csv
#gene.markers <- FindAllMarkers(b_cells, only.pos = TRUE, min.pct = 0.0, logfc.threshold = 0.0) #This will take a decent amount of time
# Save as a csv to open in excel
#write.csv(gene.markers, "/Users/nitin/Downloads/b_cells/FindAllMarkers.csv", quote = F) #Change the path to your working directory and set the file name to whatever you want



#Single-cell RNA-seq analysis - Pseudo bulk DE analysis with DESeq2
# Create single cell experiment object
sce <- as.SingleCellExperiment(b_cells)

## Check the counts matrix 
dim(counts(sce))
counts(sce)[1:6, 1:6]

# Explore the cellular metadata for the dataset
dim(colData(sce))
head(colData(sce))

# Extract the unique names of the variables
age_group <- unique(sce$Age_group)
cluster_names <- unique(sce$Cluster_names)
donor_id <- unique(sce$Donor_id)

# Check that the length makes sense and that it's returning the correct variables
length(age_group)
age_group
length(cluster_names)
cluster_names
length(donor_id)
donor_id

# Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
groups <- colData(sce)[, c("Cluster_names", "Age_group", "Donor_id")]
head(groups)

# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                                groupings = groups, fun = "sum")

# Explore output matrix
class(aggr_counts)
dim(aggr_counts)
aggr_counts[1:6, 1:6]


# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts <- t(aggr_counts)
aggr_counts[1:6, 1:6]

# Understanding tstrsplit()

## Exploring structure of function output (list)
tstrsplit(colnames(aggr_counts), "_") %>% str()

## Comparing the first 10 elements of our input and output strings

head(colnames(aggr_counts), n = 45)
head(tstrsplit(colnames(aggr_counts), "_")[[1]], n = 45)


# Using which() to look up tstrsplit() output
b_cell_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == "Non-switched memory")
b_cell_idx

colnames(aggr_counts)[b_cell_idx]
aggr_counts[1:10, b_cell_idx]

# Loop over all cell types to extract corresponding counts, and store information in a list

## Initiate empty list
counts_ls <- list()

for (i in 1:length(cluster_names)) {
  
  # Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
  
  # Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
  
}

# Explore the different components of the list
str(counts_ls)



# Review issue where metadata doesn't merge with the complete counts_ls matrix
# Current issue: Donor_id has REPEATING TERMS, with the same donor id through multiple samples of different cell groups
# This prevents the row names from becoming the Donor IDs since there are repeating row names
# In addition, changing the rows selected messes with the dimensions of the matrix, when we only want 5 rows for each sample ID



# Extract sample-level variables
metadata <- colData(sce) %>% 
  as.data.frame() %>% 
  dplyr::select(Age_group, Cluster_names, Donor_id) 

dim(metadata)
head(metadata)

# Exclude duplicated rows
metadata <- metadata[!duplicated(metadata), ]

# Add cell id for downstream dds
metadata$cell_id <- paste(metadata$Age_group, 
                          metadata$Cluster_names,
                          sep = "_")

# Add sample id to keep track of replicates per cell_id group
metadata$sample_id <- paste(metadata$Age_group, 
                            metadata$Cluster_names,
                            metadata$Donor_id,
                            sep = "_")

dim(metadata)
head(metadata)


# Rename rows
rownames(metadata) <- metadata$sample_id
head(metadata)


# Number of cells per sample and cluster
t <- table(colData(sce)$Age_group,
           colData(sce)$Cluster_names)
t[1:5, 1:6]
t


df1 <- as.data.frame(as.table(t))
colnames(df1) <- c("AgeGroup", "CellType", "Count")

ggplot(df1, aes(x = AgeGroup, y = Count, fill = AgeGroup)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Cell Count Comparison by Age Group", x = "Age Group", y = "Count") +
  scale_fill_manual(values = c("A" = "lightblue", "B" = "lightgreen", "C" = "lightcoral", "D" = "lightpink", "E" = "lightyellow")) +
  facet_wrap(~ CellType, scales = "free_y")

ggsave(paste0("C:/Users/nitin/Downloads/b_cells/results/", "CellCountDistributionGraph.png"))


# Creating metadata list

## Initiate empty list
metadata_ls <- list()

for (i in 1:length(counts_ls)) {
  ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(sample_id = colnames(counts_ls[[i]]))
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$Cluster_names <- tstrsplit(df$sample_id, "_")[[1]]
  df$Age_group  <- tstrsplit(df$sample_id, "_")[[2]]
  df$Donor_id <- tstrsplit(df$sample_id, "_")[[3]]
  df$cell_id <- paste(df$Age_group, 
                      df$Cluster_names,
                      sep = "_")
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t) == unique(df$Cluster_names))
  cell_counts <- t[, idx]
  
  ## Remove samples with zero cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 0]
  
  ## Match order of cell_counts and sample_ids
  sample_order <- match(df$Age_group, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts
  
  
  #REVIEW ISSUE: Cannot merge metadata and df data to add the metadata information to DF
  #continued past this because metadata not really required for the DESeq2 analysis, however good to have
  
  ## Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- plyr::join(df, metadata, by = intersect(names(df), names(metadata)))
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$sample_id
  
  ## Store complete metadata for cluster i in list
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$Cluster_names)
  
}
names(df)
names(metadata)
# Explore the different components of the list
str(metadata_ls)

head(df)




#Begin Multi Cell Type DE Analysis

clustx <- "Non-switched memory"

idx <- which(names(counts_ls) == "Non-switched memory")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]

# Print error message if sample names do not match
if ( all(colnames(cluster_counts) != rownames(cluster_metadata)) ) {
  print("ERROR: sample names in counts matrix columns and metadata rows do not match!")
}

cluster_specific_counts_integers <- round(cluster_counts)

dds <- DESeqDataSetFromMatrix(cluster_specific_counts_integers, 
                              colData = cluster_metadata, 
                              design = ~ Age_group)

# Transform counts for data visualization
param <- SnowParam(workers = 16, type = "SOCK") # Adjust the number of workers as needed


rld <- rlog(dds, blind = TRUE)


# Generate QC plots

## Plot and save PCA plot
DESeq2::plotPCA(rld, intgroup = "Age_group")
ggsave(paste0("C:/Users/nitin/Downloads/b_cells/results/", clustx, "_specific_PCAplot.png"))

## Extract rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

## Plot and save heatmap
png(paste0("C:/Users/nitin/Downloads/b_cells/results/", clustx, "_specific_heatmap.png"),
    height = 6, width = 7.5, units = "in", res = 300)
pheatmap(rld_cor, annotation = cluster_metadata[, c("Age_group"), drop = FALSE])
dev.off()


# Run DESeq2 differential expression analysis
dds <- DESeq(dds, parallel = TRUE, BPPARAM = param)


## Plot dispersion estimates
png(paste0("C:/Users/nitin/Downloads/b_cells/results/", clustx, "_dispersion_plot.png"),
    height = 5, width = 6, units = "in", res = 300)
plotDispEsts(dds)
dev.off()

## Output and shrink results of Wald test for contrast A vs B
contrast <- paste(c("Age_group", "E", "vs", "A"), collapse = "_")
resultsNames(dds)

res <- results(dds, name = contrast, alpha = 0.05, parallel = TRUE, BPPARAM = param)
res <- lfcShrink(dds, coef = contrast, res = res, parallel = TRUE, BPPARAM = param)

## Turn the results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()

write.csv(res_tbl,
          paste0("C:/Users/nitin/Downloads/b_cells/results/", clustx, "_", contrast, "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)

## Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

write.csv(sig_res,
          paste0("C:/Users/nitin/Downloads/b_cells/results/", clustx, "_", contrast, "_signif_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)


# Generate results visualization plots

## Extract normalized counts from dds object
normalized_counts <- counts(dds, normalized = TRUE)

## Extract top 20 DEG from resLFC (make sure to order by padj)
top20_sig_genes <- sig_res %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n = 20)

## Extract matching normalized count values from matrix
top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes,]

numrowsgenes <- nrow(top20_sig_counts)

if (is.null(numrowsgenes)) {
  numrowsgenes <- 1
}

top20_sig_counts <- matrix(top20_sig_counts, nrow = numrowsgenes, dimnames = list(top20_sig_genes, colnames(normalized_counts)))


## Convert wide matrix to long data frame for ggplot2
top20_sig_df <- data.frame(top20_sig_counts)
top20_sig_df$gene <- rownames(top20_sig_counts)

top20_sig_df <- melt(setDT(top20_sig_df), 
                     id.vars = c("gene"),
                     variable.name = "sample_id") %>% 
  data.frame()

## Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
top20_sig_df$sample_id <- gsub("\\.", " ", top20_sig_df$sample_id)
top20_sig_df$sample_id <- gsub("\\  ", "+ ", top20_sig_df$sample_id)
top20_sig_df$sample_id <- gsub("Non switched memory", "Non-switched memory", top20_sig_df$sample_id)

## Join counts data frame with metadata
top20_sig_df <- plyr::join(top20_sig_df, as.data.frame(colData(dds)),
                           by = "sample_id")

## Generate plot
ggplot(top20_sig_df, aes(y = value, x = Age_group, col = Age_group)) +
  geom_jitter(height = 0, width = 0.15) +
  scale_y_continuous(trans = 'log10') +
  ylab("log10 of normalized expression level") +
  xlab("condition") +
  ggtitle("Top 20 Significant DE Genes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ gene)

ggsave(paste0("C:/Users/nitin/Downloads/b_cells/results/", clustx, "_", contrast, "_top20_DE_genes.png"))


b_cells_Non_Switched <- subset(b_cells, subset = Cluster_names == "Non-switched memory")

express_data <- FetchData(b_cells_Non_Switched, vars = c("LILRA4", "Age_group")) %>%
  filter(LILRA4 > 0)

b_cells_Non_Switched <- subset(b_cells_Non_Switched, cells = rownames(express_data))

VlnPlot(
  object = b_cells_Non_Switched,
  features = "log2_nCount",
  log = TRUE,
  group.by = "Age_group",
  pt.size = .25,
) +
  labs(title = "LILRA4 log2_nCounts in Non-Switched Memory B Cells by Age Group", x = "Age Group", y = "LILRA4 Expression Level") +
  theme_minimal()