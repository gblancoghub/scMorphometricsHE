library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
# install.packages("harmony")
library(harmony) # !! NEW: Required for Harmony integration !!


# 

setwd("D:/EVOSdata/Ceci/Desarrollo2025/CPro(2)/CCA/IntegracionHarmony")
getwd()


# 1. Load Individual Objects and Add Metadata
OBJ00 <- readRDS(file = "CProSeuratOBJNEG_00.rds")
OBJ00$Patient <- "00_Rosetto" # Add patient ID
OBJ01 <- readRDS(file = "CProSeuratOBJNEG_01.rds")
OBJ01$Patient <- "01_Macchia"
OBJ02 <- readRDS(file = "CProSeuratOBJNEG_02.rds")
OBJ02$Patient <- "02_Thomas"
OBJ04 <- readRDS(file = "CProSeuratOBJNEG_04.rds")
OBJ04$Patient <- "04_Paciente4"


# Create a list of Seurat objects
obj.list <- list(OBJ00, OBJ01, OBJ02, OBJ04)


setToExclude <- read.csv("Exclude.csv")
setToExclude  <- gsub("_", "-", setToExclude[,1])

features_for_pca_per_obj <- setdiff(rownames(obj.list[[1]]), setToExclude )

# 2. Pre-process Each Object Independently (Normalize/FindVariableFeatures)
# Harmony, like CCA, benefits from finding common features first.
for (i in 1:length(obj.list)) {
  obj.list[[i]] <- NormalizeData(obj.list[[i]])
  # Finding variable features is less critical for Harmony but good practice
  obj.list[[i]] <- FindVariableFeatures(obj.list[[i]]) 
}

# 3. Merge All Objects into a Single Seurat Object
# Harmony operates on a single object where the 'batch' is a metadata column.
seurat.merged <- merge(obj.list[[1]], 
                       y = obj.list[2:length(obj.list)], 
                       add.cell.ids = c("Pac1", "Pac2", "Pac3", "Pac4"), 
                       project = "CLL_Morpho")

# Ensure the 'Patient' metadata column is present in the merged object
# It was added during the individual object loading, so this should be fine.
# We will use this column as the 'group.by.vars' for Harmony.

# 4. Scale Data on the merged object (required before PCA)
# Scale using only the features selected for PCA
seurat.merged <- ScaleData(seurat.merged, 
                           features = features_for_pca_per_obj)

# 5. Run PCA on the merged, scaled object
# Determine the number of dimensions (e.g., 50) based on an elbow plot 
# or a reasonable maximum for morphometric data.
seurat.merged <- RunPCA(seurat.merged, 
                        features = features_for_pca_per_obj, 
                        npcs = 50, # Arbitrary number; check elbow plot later
                        verbose = FALSE)


# Visualize the unintegrated PCA/UMAP for comparison (optional)
# seurat.merged <- RunUMAP(seurat.merged, dims = 1:50)
 seurat.merged  <- RunUMAP(seurat.merged, min.dist = 0.9, n.neighbors = 200L, spread=1, dims = 1:20)


DimPlot(seurat.merged, group.by = "Patient") # This will show batch effects
DimPlot(seurat.merged, split.by = "Patient") # This will show batch effects


# 6. Harmony Integration 
# RunHarmony corrects the PCA embeddings stored in the 'pca' reduction.
# 'group.by.vars' tells Harmony which metadata column defines the batches.
# 'assay' specifies which assay to use (your 'RNA' assay containing morphometrics)

seurat.integrated.harmony <- RunHarmony(seurat.merged, 
                                        group.by.vars = "Patient", 
                                        assay.use = "RNA", 
                                        reduction = "pca", # Apply to the 'pca' reduction
                                        reduction.save = "harmony", # Name the new, integrated reduction
                                        verbose = FALSE)

# The integrated data is now stored in the 'harmony' reduction slot.
# Unlike CCA, no new 'integrated' assay is created in Seurat 5 by Harmony.
# The downstream steps will use the 'harmony' reduction instead of the 'integrated' assay.


# 7. Clustering and UMAP/t-SNE on Harmony Embeddings
# The subsequent steps use the batch-corrected 'harmony' dimensions.

# Find Neighbors (using the 'harmony' reduction)
seurat.integrated.harmony <- FindNeighbors(seurat.integrated.harmony, 
                                             reduction = "harmony", 
                                             dims = 1:50)

# 

# Find Clusters (Louvain algorithm)
# seurat.integrated.harmony <- FindClusters(seurat.integrated.harmony, 
#                                            resolution = 0.8) # Adjust resolution as needed

seurat.integrated.harmony <- FindClusters(seurat.integrated.harmony, resolution = c(0.1, 0.3, 0.5, 0.7, 1))


# Run UMAP (using the 'harmony' reduction)
# seurat.integrated.harmony <- RunUMAP(seurat.integrated.harmony, 
#                                        reduction = "harmony", 
#                                        dims = 1:50)

seurat.integrated.harmony  <- RunUMAP(seurat.integrated.harmony, reduction = "harmony", min.dist = 0.9, n.neighbors = 200L, spread=1, dims = 1:20)


# 8. Save the Integrated Object
saveRDS(seurat.integrated.harmony, "seurat_integrated_harmony.rds")
# seurat.integrated.harmony <- readRDS("seurat_integrated_harmony.rds")

# To visually compare the integration:
# DimPlot(seurat.integrated.harmony, group.by = "Patient", reduction = "umap")
DimPlot(seurat.integrated.harmony, group.by = "Patient") # This will show batch effects
DimPlot(seurat.integrated.harmony, split.by = "Patient", label=TRUE) # This will show batch effects
DimPlot(seurat.integrated.harmony, group.by = "seurat_clusters",split.by = "Patient", label=TRUE)

view(seurat.integrated.harmony@meta.data)

seurat.integrated.harmony@meta.data$seurat_clusters
seurat.integrated.harmony@meta.data$RNA_snn_res.0.1
seurat.integrated.harmony@meta.data$RNA_snn_res.0.3
seurat.integrated.harmony@meta.data$RNA_snn_res.0.5
seurat.integrated.harmony@meta.data$RNA_snn_res.0.7
seurat.integrated.harmony@meta.data$RNA_snn_res.1


metadata <- seurat.integrated.harmony@meta.data
Proportions <- table(metadata[,c("Patient","seurat_clusters")])
Proportions <- table(metadata[,c("Patient","seurat_clusters")])


write.csv( Proportions, "ProportionsHarmony.csv")

