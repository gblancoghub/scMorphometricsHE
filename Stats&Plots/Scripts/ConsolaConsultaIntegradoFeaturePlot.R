
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)




seurat.integrated <- readRDS("seurat_integrated_wCellType.rds")


view(seurat.integrated@meta.data)

DimPlot(seurat.integrated, reduction = 'umap', group.by = "CellType", label = TRUE)
DimPlot(seurat.integrated, reduction = 'umap', group.by = "CellType", split.by = "Patient", label = TRUE)


table(seurat.integrated@meta.data$CellType,seurat.integrated@meta.data$Patient )
table(seurat.integrated@meta.data$Patient, seurat.integrated@meta.data$CellType)


DimPlot(seurat.integrated, reduction = 'umap', group.by = "BroadCellType", label = TRUE)
DimPlot(seurat.integrated, reduction = 'umap', group.by = "BroadCellType", split.by = "Patient", label = TRUE)

DimPlot(seurat.integrated, reduction = 'tsne', group.by = "BroadCellType", split.by = "Patient", label = TRUE)
DimPlot(seurat.integrated, reduction = 'tsne', group.by = "BroadCellType", label = TRUE)

table(seurat.integrated@meta.data$BroadCellType,seurat.integrated@meta.data$Patient )
table(seurat.integrated@meta.data$Patient, seurat.integrated@meta.data$BroadCellType)

str(seurat.integrated@assays$RNA)
str(seurat.integrated@assays$RNA$data) # 516
str(seurat.integrated@assays$RNA$scale.data) # 99


str(seurat.integrated@assays$integrated@data)  # 99
str(seurat.integrated@assays$integrated@scale.data)  # 87

# Get all current features
all_morph_features <- rownames(seurat.integrated@assays$integrated@scale.data) 
all_morph_features <- rownames(seurat.integrated@assays$RNA$data) 

AP <- all_morph_features 

DefaultAssay(seurat.integrated) <- "integrated"
DefaultAssay(seurat.integrated) <- "RNA"


FeaturePlot(seurat.integrated, features = AP[184], min.cutoff = 'q50', split.by = "Patient",
pt.size=1.8, order =TRUE, cols= c("green","red"))+
theme(legend.position = c(0.04,0.85))

selec <- read.csv("0-Seleccion.txt")
selec[,1]
AP[selec[,1]]
AP <- AP[selec[,1]]

for (i in 1:length(AP)) {
sp <- FeaturePlot(seurat.integrated, features = AP[i], split.by = "Patient",
min.cutoff = 'q40', pt.size=1.8, order =TRUE, cols= c("green","red"))+
theme(legend.position = c(0.04,0.85))

name <- paste(i,"_", AP[i],sep="")
plotPNG(name,sp)
}



#-----------------------------------------------------------
# Funcion ploteo individual en png
plotPNG <- function(name,p){
WithExt <- paste(name,".png",sep="")
# png(WithExt,width = 780, height = 780,)
png(WithExt,width = 1800, height = 480,)
 print(p)
dev.off()
}
#-----------------------------------------------------------
# FIN Funcion ploteo individual en png


