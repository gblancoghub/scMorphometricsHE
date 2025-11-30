
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)



setwd("D:/EVOSdata/Ceci/Desarrollo2025/temp/A-HarmonyML/2-CrearCellTypesHarmonyFull")
getwd()


# 5. Integrate Data

#-------------------------------------------
# LEO DATASET INTEGRADO
#-------------------------------------------
seurat.integrated <- readRDS("seurat_integrated_harmonyX6RTFull.rds")


#-------------------------------------------
# CREO LEL CELL TYPE DE PRIMER NIVEL SEGUN TRAYECTORIA
#-------------------------------------------

# Load the dplyr library if you haven't already
library(dplyr)
seurat.integrated@meta.data <- seurat.integrated@meta.data %>%
  mutate(
    CellType = case_when(
      seurat_clusters %in% c("23", "0", "1", "2") ~ "Root",
      seurat_clusters %in% c("17", "18","19","20","21","22","9","24","25","26","27") ~ "Stem1",
      seurat_clusters %in% c("3","5","8","10","11","12") ~ "Stem2",
      seurat_clusters %in% c("15","28","31") ~ "Diff",
      seurat_clusters %in% c("29","30","7","13","14","31","4","35","16","32","33","34") ~ "Mat",
      seurat_clusters %in% c("6","36") ~ "MatRe",

      TRUE ~ NA_character_ # Fallback for any unassigned clusters (shouldn't happen with full coverage)
    )
  )

# Root	 23-0-1-2
# Stem1	"17", "18","19","20","21","22","9","24","25","26","27"
# Stem2	"3","5","8","10","11","12"
# Diff	"15","28","31"
# Mat	 "29","30","7","14","13","31","4","35","16","32","33","34"
# MatRe	6-36

#-------------------------------------------
# CONVIERTO CELL TYPES A FACTOR
#-------------------------------------------

# Optional: Convert CellType to a factor with a specific order if desired for plotting
# (e.g., to control the order in legends or facets)
seurat.integrated@meta.data$CellType <- factor(
  seurat.integrated@meta.data$CellType,
  levels = c("Root", "Stem1", "Stem2", "Diff", "Mat", "MatRe")
)

#-------------------------------------------
# MIRO LOS % POR CELL TYPE
#-------------------------------------------

table(seurat.integrated@meta.data$CellType)

# You can also visualize it with DimPlot

#-------------------------------------------
# MIRO LAS DISTRIBUCION POR CELL TYPE Y POR PACIENTE
#-------------------------------------------


DimPlot(seurat.integrated, reduction = 'umap', group.by = "CellType", label = TRUE, pt.size=2)
DimPlot(seurat.integrated, reduction = 'umap', group.by = "CellType", split.by = "Patient", label = TRUE,pt.size=2)
DimPlot(seurat.integrated, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)


plot <- DimPlot(seurat.integrated, reduction = 'umap', group.by = "CellType", split.by = "Patient", pt.size=2)
plot_repelled_clusters <- LabelClusters(plot = plot, id = "CellType", repel = TRUE)
print(plot_repelled_clusters)


#-------------------------------------------
# CREO TABLAS CELL TYPE X PACIENTE Y PACIENTE X CELL TYPE
#-------------------------------------------

table(seurat.integrated@meta.data$CellType,seurat.integrated@meta.data$Patient )
x <- table(seurat.integrated@meta.data$Patient, seurat.integrated@meta.data$CellType)
write.csv (x,"TablaCellTypesHarmonyFull.csv")


# saveRDS(seurat.integrated, "seurat_integrated_HarmonyFull_wCellType.rds")


