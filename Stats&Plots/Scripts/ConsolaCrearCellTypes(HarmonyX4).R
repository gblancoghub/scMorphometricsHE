
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)



setwd("D:/EVOSdata/Ceci/Desarrollo2025/temp/A-HarmonyML/1-CrearCellTypes")
getwd()


# 5. Integrate Data

#-------------------------------------------
# LEO DATASET INTEGRADO
#-------------------------------------------
seurat.integrated <- readRDS("seurat_integrated_harmony.rds")


#-------------------------------------------
# CREO LEL CELL TYPE DE PRIMER NIVEL SEGUN TRAYECTORIA
#-------------------------------------------

# Load the dplyr library if you haven't already
library(dplyr)
seurat.integrated@meta.data <- seurat.integrated@meta.data %>%
  mutate(
    CellType = case_when(
      seurat_clusters %in% c("3", "0", "9") ~ "Root",
      seurat_clusters %in% c("17", "11","4","1") ~ "Stem",
      seurat_clusters %in% c("14", "15", "5", "16", "20") ~ "AltDiff",
      seurat_clusters %in% c("2", "10", "19","8") ~ "Diff",
      seurat_clusters %in% c("6", "22", "13", "7", "12", "21","23", "18") ~ "Mat",
      TRUE ~ NA_character_ # Fallback for any unassigned clusters (shouldn't happen with full coverage)
    )
  )

#-------------------------------------------
# CONVIERTO CELL TYPES A FACTOR
#-------------------------------------------

# Optional: Convert CellType to a factor with a specific order if desired for plotting
# (e.g., to control the order in legends or facets)
seurat.integrated@meta.data$CellType <- factor(
  seurat.integrated@meta.data$CellType,
  levels = c("Root", "Stem", "AltDiff", "Diff", "Mat")
)

#-------------------------------------------
# MIRO LOS % POR CELL TYPE
#-------------------------------------------

table(seurat.integrated@meta.data$CellType)

# You can also visualize it with DimPlot

#-------------------------------------------
# MIRO LAS DISTRIBUCION POR CELL TYPE Y POR PACIENTE
#-------------------------------------------


DimPlot(seurat.integrated, reduction = 'umap', group.by = "CellType", label = TRUE)
DimPlot(seurat.integrated, reduction = 'umap', group.by = "CellType", split.by = "Patient", label = TRUE)
DimPlot(seurat.integrated, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)

#-------------------------------------------
# CREO TABLAS CELL TYPE X PACIENTE Y PACIENTE X CELL TYPE
#-------------------------------------------

table(seurat.integrated@meta.data$CellType,seurat.integrated@meta.data$Patient )
x <- table(seurat.integrated@meta.data$Patient, seurat.integrated@meta.data$CellType)
write.csv (x,"TablaCellTypesHarmony.csv")


# saveRDS(seurat.integrated, "seurat_integrated_Harmony_wCellType.rds")
seurat.integrated <- readRDS("seurat_integrated_Harmony_wCellType.rds")

