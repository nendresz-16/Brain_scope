setwd("P1_Brain_scope")

#Load packages 
library(readxl)
library(readr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)

#Load in Reference
Counts_ref <- readRDS("/project/rrg-shreejoy/nendresz/raw_counts_ref.rds")

counts_hodge_matrix <- as.matrix(Counts_ref)
counts_hodge <- Matrix(counts_hodge_matrix, sparse=TRUE)

meta_hodge <- readRDS("/project/rrg-shreejoy/nendresz/Neurotypical_ref_metadata.rds")

all(rownames(counts_hodge) %in% rownames(meta_hodge))  # Ensure they align

#Load the testing count matrix
counts_sn <- readRDS("/project/rrg-shreejoy/nendresz/Brain_scope/UCLA_counts.rds")
counts_sn <- t(counts_sn)
meta_sn <- read_tsv("/project/rrg-shreejoy/nendresz/Brain_scope/Cell_metadata/UCLA-ASD_cell_metadata_mismatches_removed.tsv")
# convert tibble → data.frame
meta_sn <- as.data.frame(meta_sn)

# set rownames to barcodes
rownames(meta_sn) <- meta_sn$barcodekey

meta_sn <- meta_sn[rownames(counts_sn), ]
all(rownames(counts_sn) == rownames(meta_sn))

# filter counts matrices
## Get the common gene names between the two datasets
common_genes <- intersect(colnames(counts_hodge), colnames(counts_sn)) #32550 common genes

## Filter each counts matrix to include only the common genes
counts_hodge <- counts_hodge[, common_genes]
counts_sn <- counts_sn[, common_genes]

### SEURAT INTEGRATION

Seu_hodge_for_int <- CreateSeuratObject(counts = t(counts_hodge), 
                                       meta.data = meta_hodge) 

Seu_sn_for_int <- CreateSeuratObject(counts = t(counts_sn), meta.data = meta_sn) 



Seu.list <- c(Seu_hodge_for_int, Seu_sn_for_int)

rm(meta_sn, common_genes, counts_hodge,      
counts_hodge_matrix, counts_sn, meta_hodge)


# normalize and identify variable features for each dataset independently
Seu.list <- lapply(X = Seu.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 1000000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  #print(x)
})



# Extract the reference and query datasets
Seu_hodge <- Seu.list[[1]]  # Reference dataset
Seu_sn <- Seu.list[[2]]     # Query dataset

# Scale, PCA, and clustering for the reference dataset
Seu_hodge <- ScaleData(Seu_hodge)
Seu_hodge <- RunPCA(Seu_hodge)
Seu_hodge <- FindNeighbors(Seu_hodge, dims = 1:30)
Seu_hodge <- FindClusters(Seu_hodge)

# Find transfer anchors using the reference and query datasets
anchors <- FindTransferAnchors(reference = Seu_hodge, query = Seu_sn, dims = 1:30, 
                               reference.reduction = "pca")

# Transfer data (e.g., cell types) from reference to query
predictions <- TransferData(anchorset = anchors, refdata = Seu_hodge$Supertype, dims = 1:30)

# Add predictions as metadata to the query dataset
Seu_sn <- AddMetaData(Seu_sn, metadata = predictions)


#Read old counts
counts_sn <- readRDS("/project/rrg-shreejoy/nendresz/Brain_scope/UCLA_counts.rds")
counts_sn <- t(counts_sn)

#Extract new metadata 
new_meta <- Seu_sn@meta.data

New_seu <- CreateSeuratObject(counts = t(counts_sn), meta.data = new_meta) 

library(ggalluvial)
library(ggplot2)
library(dplyr)


alluvial_df <- new_meta %>%
  count(subclass, predicted.id)

# Plot
p <- ggplot(alluvial_df,
       aes(axis1 = subclass, axis2 = predicted.id, y = n)) +
  geom_alluvium(aes(fill = subclass), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Subclass", "Predicted"), expand = c(.1, .05)) +
  theme_classic() +
  labs(y = "Cell count")

ggsave("Figures/Riverplot_UCLA.png",  plot = p,  width = 12, height = 20, dpi = 300)

saveRDS(New_seu, "Files/LT_UCLA_Supertypes.rds")

