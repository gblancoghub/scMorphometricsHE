
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)



getwd()


# 1. Load Individual Objects and Add Metadata
OBJ00 <- readRDS(file = "CProSeuratOBJNEG_00.rds")
OBJ00$Patient <- "00_Rosetto" # Add patient ID
OBJ01 <- readRDS(file = "CProSeuratOBJNEG_01.rds")
OBJ01$Patient <- "01_PAC02"
OBJ02 <- readRDS(file = "CProSeuratOBJNEG_02.rds")
OBJ02$Patient <- "02_PAC03"
OBJ04 <- readRDS(file = "CProSeuratOBJNEG_04.rds")
OBJ04$Patient <- "04_Paciente4"

#-------------------------------------------------
#Este prefijo se une al numero de celula o barcode

merged_seurat <- merge(OBJ00, OBJ01, add.cell.ids = c("00_PAC01", "01_PAC02"),   project ="CLL")

merged_seurat <- merge(merged_seurat , OBJ02, add.cell.ids = c("","_02_PAC03"), project ="CLL")

merged_seurat <- merge(merged_seurat , OBJ04, add.cell.ids = c("","__04_PAC04"), project ="CLL")

view(merged_seurat@meta.data)
str(merged_seurat)

#-------------------------------------------------
# creo la columna sample con rownmaes que tiene el prefijo  + el barcode original o numero de cell

merged_seurat$sample <- rownames(merged_seurat@meta.data)

#-------------------------------------------------
# LA FUNCION SEPARATE ELIMINA SAMPLE Y LA REEMPLZA POR 3 NUEVAS COLUMNAS
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = "sample", into = c("num1","num2","num3","Patient","Barcode"), sep = "_")


setToExclude <- read.csv("Exclude.csv")
setToExclude  <- gsub("_", "-", setToExclude[,1])

features_for_pca <- setdiff(rownames(merged_seurat), setToExclude )


# 5. Scaling -------------

merged_seurat <- ScaleData(merged_seurat, features = features_for_pca)



# 6. Perform Linear dimensionality reduction --------------
merged_seurat <- RunPCA(merged_seurat, features = features_for_pca)

# Find Neighbors and Clusters on the integrated data
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:20)
merged_seurat <- FindClusters(merged_seurat, resolution = c(0.1, 0.3, 0.5, 0.7, 1))



# ESTOY USANDO ESTE
merged_seurat <- RunUMAP(merged_seurat, min.dist = 0.9, n.neighbors = 200L, spread=1, dims = 1:20)

#-------------- Optiona tSNE plot for control and comparison to CPro Analyst
merged_seurat <- RunTSNE(merged_seurat, dims = 1:20)


saveRDS(merged_seurat,"MergeSinIntegrar.rds")


