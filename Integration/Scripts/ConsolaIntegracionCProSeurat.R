
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)


getwd()

# 1. Load Individual Objects and Add Metadata
OBJ00 <- readRDS(file = "CProSeuratOBJNEG_00.rds")
OBJ00$Patient <- "00_PAC01" # Add patient ID
OBJ01 <- readRDS(file = "CProSeuratOBJNEG_01.rds")
OBJ01$Patient <- "01_Macchia"
OBJ02 <- readRDS(file = "CProSeuratOBJNEG_02.rds")
OBJ02$Patient <- "02_PAC03"
OBJ04 <- readRDS(file = "CProSeuratOBJNEG_04.rds")
OBJ04$Patient <- "04_PAC04"


# Create a list of Seurat objects
obj.list <- list(OBJ00, OBJ01, OBJ02, OBJ04)

# 2. Pre-process - Each Object Independently
# 

setToExclude <- read.csv("Exclude.csv")
setToExclude  <- gsub("_", "-", setToExclude[,1])

features_for_pca_per_obj <- setdiff(rownames(obj.list[[1]]), setToExclude )


# 3. Select Integration Features
features <- SelectIntegrationFeatures(object.list = obj.list)



# 4. Find Integration Anchors (CCA)

anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

saveRDS(anchors, "anchors.rds")

# anchors <- readRDS("anchors.rds")
str(anchors)

# 5. Integrate Data
seurat.integrated <- IntegrateData(anchorset = anchors)

saveRDS(seurat.integrated,"seurat_integrated(RAW).rds")

seurat.integrated <- readRDS("seurat_integrated(RAW).rds")

# 6. Post-Integration Processing (on the integrated assay)
# Set the default assay to "integrated" for downstream analysis
DefaultAssay(seurat.integrated) <- "integrated"

# seurat.integrated <- readRDS("seurat_integrated(RAW).rds")

features_for_pca <- features_for_pca_per_obj
# seurat.integrated <- FindVariableFeatures(seurat.integrated, selection.method = "vst",
                                         nfeatures = 50, features = features_for_pca) # Apply to subset


# 5. Scaling -------------

seurat.integrated  <- ScaleData(seurat.integrated , features = features_for_pca)



# 6. Perform Linear dimensionality reduction --------------
seurat.integrated <- RunPCA(seurat.integrated, features = features_for_pca)

# Find Neighbors and Clusters on the integrated data
seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:20)
seurat.integrated <- FindClusters(seurat.integrated, resolution = c(0.1, 0.3, 0.5, 0.7, 1))



# ESTOY USANDO ESTE
seurat.integrated <- RunUMAP(seurat.integrated, min.dist = 0.9, n.neighbors = 200L, spread=1, dims = 1:20)

#-------------- Optiona tSNE plot for control and comparison to CPro Analyst
seurat.integrated <- RunTSNE(seurat.integrated, dims = 1:20)


DimPlot(seurat.integrated, reduction = "umap", group.by = "Patient") # To see if batch effects are removed
# Se ven todos los pacientes apilados perfecto
DimPlot(seurat.integrated, reduction = "tsne", group.by = "Patient") # To see if batch effects are removed
# Se ven todos los pacientes apilados perfecto
DimPlot(seurat.integrated, reduction = "pca", group.by = "Patient") # To see if batch effects are removed





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




Idents(seurat.integrated) <- "seurat_clusters"



# seurat.integrated@meta.data$seurat_clusters <- Idents(seurat.integrated)
DimPlot( seurat.integrated  , reduction="umap", split.by = "Patient", label = TRUE)

# Corremos tambien un tSNE

DefaultAssay(seurat.integrated) <- "integrated"# Run t-SNE on the integrated data (using the integrated PCA dimensions)

seurat.integrated <- RunTSNE(object = seurat.integrated, dims = 1:20, perplexity = 30)

DimPlot( seurat.integrated  , reduction="tsne", split.by = "Patient", label = TRUE)



DimPlot(seurat.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot( seurat.integrated  , reduction="umap", split.by = "Patient", label = TRUE)
DimPlot( seurat.integrated  , reduction="tsne", split.by = "Patient", label = TRUE)
DimPlot(seurat.integrated, reduction = "tsne", group.by = "seurat_clusters", label = TRUE)





#----------------------------------------------------------------------------------------


saveRDS(seurat.integrated,"seurat_integrated.rds")





# Visualize with t-SNE, coloring by your integrated clusters
DimPlot(seurat.integrated, reduction = "tsne", group.by = "seurat_clusters", label = TRUE)

# Visualize with t-SNE, splitting by patient to see batch mixing
DimPlot(seurat.integrated, reduction = "tsne", split.by = "Patient", label = TRUE)




#-----------------------------------------------------------------
# Comienza una seccion para obtner las coordenadas UMAP 
# tanto originales como luego de la integracion por CCA
#-----------------------------------------------------------------


seurat.integrated@assays$RNA@data # for the normalised values or 
seurat.integrated$RNA@counts # for the raw counts.

str(seurat.integrated$RNA@counts )

Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
  ..@ i       : int [1:37674] 1836 2297 108 472 1735 1899 2129 2131 5 67 ...
  ..@ p       : int [1:3417] 0 2 8 22 25 35 44 45 45 47 ...
  ..@ Dim     : int [1:2] 3416 3416
  ..@ Dimnames:List of 2
  .. ..$ : chr [1:3416] "AAA" "AAE" "AAL" "AAN" ...
  .. ..$ : chr [1:3416] "AAA" "AAE" "AAL" "AAN" ...
  ..@ x       : num [1:37674] 1 1 1 1 1 1 1 1 1 1 ...
  ..@ factors : list()

#-----------------------------------------------------------------
# paso la matriz esparsa de Seurat a matriz comun

Sparemat <- seurat.integrated$RNA@counts
mat <- as.matrix(Sparemat)


mat <- as.matrix(seurat.integrated@assays$RNA@layers$counts)
Sparemat_v5 <- seurat.integrated@assays$RNA@layers$counts

dim(mat)
colnames(mat)
rownames(mat)
mat <- t(mat)
mat <- as.data.frame(mat)

#-----------------------------------------------------------------
# Esta es la tabla de todo como la de CPro pero sin coordenadas UMAP

write.csv(mat,"TableSeurat.csv")
getwd()

colnames(seurat.integrated@meta.data)
metadata <- seurat.integrated@meta.data
write.csv(metadata ,"Metadata.csv")

seurat.integrated@assays$RNA@layers@counts.1.CLL.CLL.1 # for the raw counts.
str(seurat.integrated@assays$RNA@layers)
seurat.integrated@assays$RNA@layers$counts.SeuratProject.4
str(seurat.integrated@assays$RNA@layers$counts.SeuratProject.4)
seurat.integrated@assays$RNA@layers$data.SeuratProject.4 
seurat.integrated@assays$RNA@layers$scale.data.SeuratProject.4  


#-----------------------------------------------------------------
# coordenadas UMAP post integracion

umap.integrated <- as.data.frame(seurat.integrated$umap@cell.embeddings) 
write.csv(umap.integrated ,"UMAPintegrated.csv")

#-----------------------------------------------------------------
# coordenadas UMAP originales

umap.original <- merged_seurat_filtered$umap@cell.embeddings 
colnames(umap.original) <- c("UMAP_1ori", "UMAP_2ori")


write.csv(umap.original ,"UMAPoriginal.csv")

#-----------------------------------------------------------------
# Agrego las  coordenadas UMAP integradas y originales a la tabla de datos 

rownames(umap.integrated)
dat.long <-   merge(mat, umap.integrated, by = "row.names") 
rownames(dat.long) <- dat.long[,1]
dat.long <-   dat.long[,-1]




dat.long2 <-   merge(dat.long, umap.original, by = "row.names") 
# colnames(dat.long2)
# dim(dat.long2)
# colnames(dat.long)

rownames(dat.long2) <- dat.long2[,"Row.names"]
dat.long <-   dat.long[,-1]

write.csv(dat.long2  ,"TablaCPro4PacUMAPs.csv")

#-----------------------------------------------------------------
# Divido por pacientes el objeto integrado
# grabo objetos separado ya con las nueva coordenadas

obj.list <- SplitObject(seurat.integrated,split.by = "Patient")

for (i in 1:length(obj.list)){
saveRDS(obj.list[[i]], file = paste0("OBJpostInt_",i,".rds"))
}


1=PAC01
2=PAC02
3=PAC03
4-PAC04
#----------- UMAP del integrado PAC04
i=1

umap.integrated <- as.data.frame(obj.list[[i]]$umap@cell.embeddings) 
IntRes01 <- as.numeric(obj.list[[i]]@meta.data$"RNA_snn_res.0.1")
IntRes03 <- as.numeric(obj.list[[i]]@meta.data$"RNA_snn_res.0.3")
IntRes05 <- as.numeric(obj.list[[i]]@meta.data$"RNA_snn_res.0.5")

IntRes01 <- as.data.frame(IntRes01)
IntRes03 <- as.data.frame(IntRes03)
IntRes05 <- as.data.frame(IntRes05)

colnames(IntRes01) <- c("IntRes01")
colnames(IntRes03) <- c("IntRes03")
colnames(IntRes05) <- c("IntRes05")

table(IntRes01)
table(IntRes03)
table(IntRes05)


library("RSQLite")



conn <- dbConnect(RSQLite::SQLite(), "ASC_classDB.db")


# List all the tables available in the database

dbListTables(conn)

> dbListTables(conn)
 [1] "ASC_ASC_dot"                                
 [2] "ASC_Per_Experiment"                         
 [3] "ASC_Per_Image"                              
 [4] "ASC_Per_Object"                             
 [5] "ASC_Per_RelationshipTypes"                  
 [6] "ASC_Per_Relationships"                      
 [7] "ASC_Per_RelationshipsView"                  
 [8] "Experiment"                                 
 [9] "Experiment_Properties"                      
[10] "_link_columns_ASC_Per_Image_ASC_Per_Object_"
[11] "_link_tables_ASC_Per_Image_ASC_Per_Object_" 
[12] "pca_table"                                  
[13] "sqlite_sequence"                            
> 


# Get the car names and horsepower of the cars with 8 cylinders
dbGetQuery(conn,"SELECT * FROM Experiment")
dbGetQuery(conn,"SELECT field, value FROM Experiment_Properties")

dbGetQuery(conn,"SELECT * FROM ASC_Per_Experiment")
dbGetQuery(conn,"SELECT * FROM ASC_Per_Experiment")  # aca esta todo lo de las iamgenes y thumbnails

dbGetQuery(conn,"SELECT * FROM ASC_Per_Object LIMIT 10")  # aca esta todo lo de objetos - lo que uso siempre

dbGetQuery(conn,"SELECT * FROM ASC_Per_RelationshipTypes")  # 
dbGetQuery(conn,"SELECT * FROM ASC_Per_Relationships LIMIT 10")  # relacion por pares entre objetos ID_obj+ID_img

dbGetQuery(conn,"SELECT * FROM ASC_Per_RelationshipsView LIMIT 10")  # 

dbGetQuery(conn,"SELECT * FROM _link_columns_ASC_Per_Image_ASC_Per_Object_")  # 

dbGetQuery(conn,"SELECT * FROM _link_tables_ASC_Per_Image_ASC_Per_Object_")  # 

dbGetQuery(conn,"SELECT * FROM sqlite_sequence")  # 




my_old_table <- dbGetQuery(conn,"SELECT * FROM ASC_Per_Object")  
my_tableA01 <- cbind(my_old_table,IntRes01,IntRes03,IntRes05,umap.integrated)

# colnames(my_tableA01)


dbRemoveTable(conn, "ASC_Per_Object")

dbWriteTable(conn, "ASC_Per_Object", my_tableA01)

# dbGetQuery(conn,"SELECT * FROM MyExpt_K562Per_Object_UMAP LIMIT 10")  # 
# dbRemoveTable(conn, "MyExpt_K562Per_Object_UMAP")

# Close the database connection to CarsDB
dbDisconnect(conn)

