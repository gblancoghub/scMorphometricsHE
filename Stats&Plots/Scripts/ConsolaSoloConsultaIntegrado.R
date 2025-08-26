
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)



# 5. Integrate Data
seurat.integrated <- readRDS("seurat_integrated(RAW).rds")

seurat.integrated <- readRDS("seurat_integrated.rds")

# Load the dplyr library if you haven't already
library(dplyr)
seurat.integrated@meta.data <- seurat.integrated@meta.data %>%
  mutate(
    CellType = case_when(
      seurat_clusters %in% c("0", "1") ~ "M1",
      seurat_clusters %in% c("40", "43", "45", "46", "19", "37") ~ "M2",
      seurat_clusters %in% c("14", "47") ~ "M3",
      seurat_clusters %in% c("13", "11", "18", "28", "5") ~ "R0",
      seurat_clusters %in% c("35", "24") ~ "D1",
      seurat_clusters %in% c("6", "42", "30") ~ "D2",
      seurat_clusters %in% c("34", "39", "12", "3") ~ "D3",
      seurat_clusters %in% c("20", "25", "44") ~ "ADB1",
      seurat_clusters %in% c("7", "10", "23") ~ "ADB2",
      seurat_clusters %in% c("41", "8") ~ "ADBR1",
      seurat_clusters %in% c("38", "33", "9", "32") ~ "ADMR1",
      seurat_clusters %in% c("16", "17", "27", "21", "22", "31","36", "26") ~ "ADTR1",
      seurat_clusters %in% c("4", "29", "15", "2") ~ "ADT1",
      TRUE ~ NA_character_ # Fallback for any unassigned clusters (shouldn't happen with full coverage)
    )
  )

# Optional: Convert CellType to a factor with a specific order if desired for plotting
# (e.g., to control the order in legends or facets)
seurat.integrated@meta.data$CellType <- factor(
  seurat.integrated@meta.data$CellType,
  levels = c("M1", "M2", "M3", "R0", "D1", "D2", "D3", "ADB1", "ADB2", "ADBR1", "ADMR1", "ADTR1", "ADT1")
)

table(seurat.integrated@meta.data$CellType)

# visualize it with DimPlot

DimPlot(seurat.integrated, reduction = 'umap', group.by = "CellType", label = TRUE)
DimPlot(seurat.integrated, reduction = 'umap', group.by = "CellType", split.by = "Patient", label = TRUE)
DimPlot(seurat.integrated, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)

table(seurat.integrated@meta.data$CellType,seurat.integrated@meta.data$Patient )
table(seurat.integrated@meta.data$Patient, seurat.integrated@meta.data$CellType)


seurat.integrated@meta.data <- seurat.integrated@meta.data %>%
  mutate(
    BroadCellType = case_when(
      CellType %in% c("M1", "M2", "M3") ~ "Mature",
      CellType %in% c("D1", "D2", "D3") ~ "Differentiation",
      CellType %in% c("ADB1", "ADB2") ~ "AlternativeDiff_B",
      CellType %in% c("ADBR1", "ADMR1") ~ "AlternativeDiff_MR",
       CellType %in% c("ADTR1","ADT1") ~ "AlternativeDiff_TR",
     CellType %in% c("R0") ~ "RootStem",
      TRUE ~ NA_character_ # Fallback for any unassigned clusters (shouldn't happen with full coverage)
    )
  )


table(seurat.integrated@meta.data$BroadCellType)
# Optional: Convert CellType to a factor with a specific order if desired for plotting
# (e.g., to control the order in legends or facets)
seurat.integrated@meta.data$BroadCellType <- factor(
  seurat.integrated@meta.data$BroadCellType,
  levels = c("RootStem", "Differentiation", "Mature", "AlternativeDiff_B", "AlternativeDiff_MR", "AlternativeDiff_TR")
)



DimPlot(seurat.integrated, reduction = 'umap', group.by = "BroadCellType", label = TRUE)
DimPlot(seurat.integrated, reduction = 'umap', group.by = "BroadCellType", split.by = "Patient", label = TRUE)
DimPlot(seurat.integrated, reduction = 'tsne', group.by = "BroadCellType", split.by = "Patient", label = TRUE)
DimPlot(seurat.integrated, reduction = 'tsne', group.by = "BroadCellType", label = TRUE)

table(seurat.integrated@meta.data$CellType,seurat.integrated@meta.data$Patient )
table(seurat.integrated@meta.data$Patient, seurat.integrated@meta.data$CellType)




# visualize and cluster the integrated data with batch effects removed
DimPlot(seurat.integrated, reduction = "umap", group.by = "Patient") # To see if batch effects are removed
# Se ven todos los pacientes apilados perfecto
DimPlot(seurat.integrated, reduction = "tsne", group.by = "Patient") # To see if batch effects are removed
# Se ven todos los pacientes apilados perfecto
DimPlot(seurat.integrated, reduction = "pca", group.by = "Patient") # To see if batch effects are removed
# Se ven todos los pacientes apilados perfecto




DimPlot(seurat.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
# se ven todos lso clusters como un mapa politico perfecto

DimPlot( seurat.integrated  , reduction="umap", split.by = "Patient", group.by = "RNA_snn_res.0.1", label = TRUE)
# Se ve el mapa politico de cada paciente


DimPlot( seurat.integrated  , reduction="tsne", split.by = "Patient", group.by = "RNA_snn_res.0.1", label = TRUE)
# Se ve el mapa politico de cada paciente


DimPlot( seurat.integrated  , reduction="pca", split.by = "Patient", group.by = "RNA_snn_res.0.1", label = TRUE)
# Se ve el mapa politico de cada paciente

DimPlot( seurat.integrated  , reduction="umap", split.by = "Patient", group.by = "RNA_snn_res.0.7", label = TRUE)
# Se ve el mapa politico de cada paciente
DimPlot( seurat.integrated  , reduction="umap", split.by = "Patient", group.by = "seurat_clusters", label = TRUE)
# Se ve el mapa politico de cada paciente




view(seurat.integrated@meta.data)
colnames(seurat.integrated@meta.data)

# Idents(seurat.integrated) <- "RNA_snn_res.0.7"
Idents(seurat.integrated) <- "seurat_clusters"

# EL UNICO QUE FUNCIONA BIEN ES EL SLOT SEURAT CLUSTERS
# NUNCA PISARLO

# seurat.integrated@meta.data$seurat_clusters <- Idents(seurat.integrated)
DimPlot( seurat.integrated  , reduction="umap", split.by = "Patient", label = TRUE)
DimPlot( seurat.integrated  , reduction="tsne", split.by = "Patient", label = TRUE)



DimPlot(seurat.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot( seurat.integrated  , reduction="umap", split.by = "Patient", label = TRUE)
DimPlot( seurat.integrated  , reduction="tsne", split.by = "Patient", label = TRUE)
DimPlot(seurat.integrated, reduction = "tsne", group.by = "seurat_clusters", label = TRUE)






# Visualize with t-SNE, coloring by your integrated clusters
DimPlot(seurat.integrated, reduction = "tsne", group.by = "seurat_clusters", label = TRUE)

# Visualize with t-SNE, splitting by patient to see batch mixing
DimPlot(seurat.integrated, reduction = "tsne", split.by = "Patient", label = TRUE)



