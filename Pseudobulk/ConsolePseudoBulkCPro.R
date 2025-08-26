


# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)


backup <- readRDS(file = "seurat_integrated.rds")

CProAll <- backup
View(CProAll@meta.data)
colnames(CProAll@meta.data)
Idents(CProAll) 

CProAll@meta.data$"RNA_snn_res.0.1"

DimPlot(CProAll, reduction = 'umap', group.by = "Patient")

# str(CProAll)
# ? AggregateExpression()

Bulk.CProAll.NotInteg <- AggregateExpression(
  CProAll,
  return.seurat = FALSE,
  group.by = "Patient",
  slot = "counts",
  verbose = TRUE
)
# Separo la lista RNA - es la suma de las mediciones de todas las células
Bulk.CPro.Ori <- as.data.frame(Bulk.CProAll.NotInteg$RNA)

Bulk.CProAll.Integ <- AggregateExpression(
  CProAll,
  return.seurat = FALSE,
  group.by = "Patient",
  slot = "data",
  verbose = TRUE
)

# str(Bulk.CProAll.Integ)
Bulk.CPro.Integ.RNA <- as.data.frame(Bulk.CProAll.Integ$RNA)
Bulk.CPro.Integ.integrated <- as.data.frame(Bulk.CProAll.Integ$integrated)

dim(Bulk.CPro.Integ.RNA)
dim(Bulk.CPro.Integ.integrated)



table(CProAll@meta.data$Patient)
for (i in 1:nrow(Bulk.CPro.Integ.RNA)) {
Bulk.CPro.Integ.RNA[i,] <- Bulk.CPro.Integ.RNA[i,] /table(CProAll@meta.data$Patient)
}

for (i in 1:nrow(Bulk.CPro.Integ.integrated)) {
Bulk.CPro.Integ.integrated[i,] <- Bulk.CPro.Integ.integrated[i,] /table(CProAll@meta.data$Patient)
}

for (i in 1:nrow(Bulk.CPro.Ori)) {
Bulk.CPro.Ori[i,] <- Bulk.CPro.Ori[i,] /table(CProAll@meta.data$Patient)
}
# Bulk.CPro.Ori[1:25,]
# head(Bulk.CPro.Ori)


write.csv(Bulk.CPro.Ori ,"BulkCProOri.csv")
write.csv(Bulk.CPro.Integ.RNA,"BulkCProIntegratedRNA.csv")
write.csv(Bulk.CPro.Integ.integrated ,"BulkCProIntegratedCCA.csv")





