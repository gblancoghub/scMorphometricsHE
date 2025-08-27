# Construccion de un objeto Seurat
# Clustering
# Actualizacion de la base de datos sqlite de Cell profiler con los resultados de clutering y coordenadas UMAP

library(Matrix)
############################################
# Leo CPro y creo una matriz comun
#
myfile <- read.csv("my_table00.csv")

dim(myfile)
str(myfile)
colnames(myfile)

table(is.na(myfile))
myfile[is.na(myfile)] <- 0


############################################
# Modificacion de el input para desplazar los valores positivos
# Hay que correr una vez y salen todos los negativos
# En ese paso se corrigen todos los valores negativos
# Luego vuelvo a correr y tiene que dar cero columna con valores negativos

############################################
# El caso particular de Euler


vec <- myfile[,"CellObjects_Mean_ASCspecs_AreaShape_EulerNumber"]* (-1) + 2
myfile[,"CellObjects_Mean_ASCspecs_AreaShape_EulerNumber"] <- vec

df <- myfile
# Initialize a list to store information about columns with negative values
columns_with_negatives <- list()

# Loop through each column of the dataframe
for (col_name in colnames(df)) {
  current_column <- df[[col_name]] # Use double brackets to extract as a vector

  # Check if the column contains any negative values
  if (any(current_column < 0, na.rm = TRUE)) { # na.rm=TRUE handles potential NA values
    min_val <- min(current_column, na.rm = TRUE)
    max_val <- max(current_column, na.rm = TRUE)
    df[[col_name]] <- current_column + abs(min_val)
    
    # Store the information
    columns_with_negatives[[col_name]] <- list(
      min_value = min_val,
      max_value = max_val
    )

    # Print the information to the console
    cat(paste0("Column '", col_name, "' contains negative values.\n"))
    cat(paste0("  Min value: ", round(min_val, 4), "\n"))
    cat(paste0("  Max value: ", round(max_val, 4), "\n"))
    cat("  Suggested transformation to make values non-negative (min becomes 0): \n")
    cat(paste0("    `df$", col_name, "_shifted <- df$", col_name, " + ", round(abs(min_val), 4), "`\n"))
    cat(paste0("    (This shifts the range from [", round(min_val, 4), ", ", round(max_val, 4), "] to [0, ", round(max_val + abs(min_val), 4), "])\n\n"))
  }
}

myfile <- df


write.csv(myfile, "my_tableNEG00.csv")

#
# Acá terminó la correción de valores negativos
############################################





myStd.matrix <- as.matrix(myfile)

############################################
# preparo los "genes" = parametros de CPro
#

dim(myStd.matrix)
colnames(myStd.matrix)
parameters <-  as.data.frame(colnames(myStd.matrix))
# write.csv(parameters,"featuresR.csv")

# Find duplicate rows (it will return TRUE for the second occurrence of a duplicate)
duplicate_rows_logical <- duplicated(parameters)
duplicate_values <- parameters[duplicate_rows_logical, , drop = FALSE]
print(duplicate_values)

parameters_modified <- parameters
parameters_modified[[1]] <- gsub("_", "-", parameters[[1]])
parameters <- parameters_modified 


############################################
# Preparo los barcodes o cells que es una ID de celu obtenido de image+object
#

# rownames(myStd.matrix) <- paste(myStd.matrix[,"ImageNumber"],myStd.matrix[,"ObjectNumber"],sep="-")
rownames(myStd.matrix) <- sprintf("%02d-%03d", myStd.matrix[,"ImageNumber"], myStd.matrix[,"ObjectNumber"])
cells <-  as.data.frame(rownames(myStd.matrix))
# write.csv(cells ,"barcodesR.csv")


# Find duplicate rows (it will return TRUE for the second occurrence of a duplicate)
duplicate_rows_logical <- duplicated(cells)
duplicate_values <- cells[duplicate_rows_logical, , drop = FALSE]
print(duplicate_values)

############################################
# hago la transpuesta de la matriz comun
#
myStd.matrixT <- t(myStd.matrix)
head(myStd.matrixT)

rownames(myStd.matrixT)
colnames(myStd.matrixT)

# creo la matriz esparsa
sparse.gbm <- Matrix(myStd.matrixT, sparse = T )
head(sparse.gbm)

writeMM(obj = sparse.gbm, file="matrix.mtx")


# save genes and cells names
write(x = rownames(myStd.matrixT), file = "genes.tsv")
write(x = colnames(myStd.matrixT), file = "barcodes.tsv")

############################################
# Me dejï¿½ en la capeta matrix.mtx, features y barcodes
#

# Leo con Seurat



# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
# library(DoubletFinder)

# create counts matrix
cts <- ReadMtx(mtx = 'matrix.mtx',
        features = 'genes.tsv',
        feature.column = 1,
        cells = 'barcodes.tsv')


# create Seurat object
pbmc.seurat <- CreateSeuratObject(counts = cts)
# str(pbmc.seurat)
nsclc.seurat.obj <- pbmc.seurat



# 3. Normalize data ----------
#nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# OR
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)
# str(nsclc.seurat.obj)






# 4. Identify highly variable features --------------
# nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "mvp", nfeatures = 60)
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 50)

# ? FindVariableFeatures
# selection.method = "vst" - es el default

# Identify the 10 most highly variable genes
top30 <- head(VariableFeatures(nsclc.seurat.obj), 30)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top30, repel = TRUE)



# Identify the features to exclude
# features_to_exclude <- c("ImageNumber", "ObjectNumber","CellObjects-Number-Object-Number") # Add any other features you don't want to use

features_to_exclude <- read.csv("Exclude.csv")
features_to_exclude  <- gsub("_", "-", features_to_exclude[,1])

# Get all current features
all_morph_features <- rownames(nsclc.seurat.obj[["RNA"]]@features) # This assumes you're using 'RNA' assay for your morphological data

# Select features to include in PCA by removing the unwanted ones
features_for_pca <- setdiff(all_morph_features, features_to_exclude)



# 5. Scaling -------------

nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = features_for_pca)


# 6. Perform Linear dimensionality reduction --------------
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = features_for_pca)


# visualize PCA results
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(nsclc.seurat.obj, dims = 1:2, cells = 500, balanced = TRUE)

# determine dimensionality of the data
ElbowPlot(nsclc.seurat.obj)



# 7. Clustering ------------
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:20)

# several resolutions
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(nsclc.seurat.obj@meta.data)

DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.1", label = TRUE)

# setting identity of clusters
Idents(nsclc.seurat.obj)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"
Idents(nsclc.seurat.obj)




# UMAP
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, min.dist = 0.9, n.neighbors = 200L, spread=1, dims = 1:20)

#-------------- Optiona tSNE plot for control and comparison to CPro Analyst
 nsclc.seurat.obj <- RunTSNE(nsclc.seurat.obj, dims = 1:20)
# 5. Plot the tSNE
DimPlot(nsclc.seurat.obj, reduction = "tsne")
#-----------------------------------------------------------------------------------------

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

DimPlot(nsclc.seurat.obj, reduction = "umap")

DimPlot(nsclc.seurat.obj, reduction = "pca")
Idents(nsclc.seurat.obj) <- "seurat_clusters"
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.5"
DimPlot(nsclc.seurat.obj, reduction = "umap", label = TRUE)


Idents(nsclc.seurat.obj) <- "RNA_snn_res.1"
DimPlot(nsclc.seurat.obj, reduction = "umap")
DimPlot(nsclc.seurat.obj, reduction = "tsne")

# Save an object to a file - Explicacion genï¿½rcia
saveRDS(nsclc.seurat.obj, file = "CProSeuratOBJNEG_00.rds")
# Restore the object
# readRDS(file = "CProSeuratOBJ_00.rds")

DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.1", label = TRUE)

View(nsclc.seurat.obj@meta.data)
colnames(nsclc.seurat.obj@meta.data)
str(nsclc.seurat.obj)
str(nsclc.seurat.obj@assays$RNA@data)


CellObjects-Mean-ASCspecs-AreaShape-EulerNumber
                 
FeaturePlot(nsclc.seurat.obj, features = c('CellObjects-AreaShape-Area'), min.cutoff = 'q10')
FeaturePlot(nsclc.seurat.obj, features = c('CellObjects-Mean-ASCspecs-Intensity-IntegratedIntensity-HEGray'), min.cutoff = 'q50')
FeaturePlot(nsclc.seurat.obj, features = c('CellObjects-Mean-Nucelolos-AreaShape-MeanRadius'), min.cutoff = 'q50')
FeaturePlot(nsclc.seurat.obj, features = c('CellObjects-Mean-Nucelolos-AreaShape-MajorAxisLength'), min.cutoff = 'q50')
FeaturePlot(nsclc.seurat.obj, features = c('CellObjects-Mean-ASCspecs-AreaShape-EulerNumber'), min.cutoff = 'q50')
FeaturePlot(nsclc.seurat.obj, features = c('CellObjects-Mean-ASCspecs-AreaShape-EulerNumber'), min.cutoff = 'q50', slot="scale.data")
FeaturePlot(nsclc.seurat.obj, features = c('CellObjects-Mean-ASCspecs-AreaShape-EulerNumber'), min.cutoff = 'q50', slot="data")
FeaturePlot(nsclc.seurat.obj, features = c('CellObjects-Mean-ASCspecs-AreaShape-EulerNumber'), order=TRUE, min.cutoff = 'q50', slot="counts", pt.size=1)



###------------------------------------------

library("RSQLite")



conn <- dbConnect(RSQLite::SQLite(), "ASC_classDB.db")


# List all the tables available in the database

dbListTables(conn)
> dbListTables(conn)
 [1] "ASC_ASC_dot"                                 "ASC_Per_Experiment"                         
 [3] "ASC_Per_Image"                               "ASC_Per_Object"                             
 [5] "ASC_Per_RelationshipTypes"                   "ASC_Per_Relationships"                      
 [7] "ASC_Per_RelationshipsView"                   "Experiment"                                 
 [9] "Experiment_Properties"                       "_link_columns_ASC_Per_Image_ASC_Per_Object_"
[11] "_link_tables_ASC_Per_Image_ASC_Per_Object_"  "pca_table"                                  
[13] "sqlite_sequence"                            
>

# controls exploring DBase
dbGetQuery(conn,"SELECT * FROM Experiment")
dbGetQuery(conn,"SELECT field, value FROM Experiment_Properties")

dbGetQuery(conn,"SELECT * FROM ASC_Per_Experiment")
dbGetQuery(conn,"SELECT * FROM ASC_Per_Experiment")  # aca esta todo lo de las iamgenes y thumbnails

dbGetQuery(conn,"SELECT * FROM ASC_Per_Object LIMIT 10")  # aca esta todo lo de objetos - lo que uso siempre

dbGetQuery(conn,"SELECT * FROM ASC_Per_RelationshipTypes")  # 
dbGetQuery(conn,"SELECT * FROM ASC_Per_Relationships LIMIT 10")  # relacion por pares entre objetos ID_obj+ID_img

dbGetQuery(conn,"SELECT * FROM ASC_Per_RelationshipsView LIMIT 10")  # 

dbGetQuery(conn,"SELECT * FROM _link_tables_ASC_Per_Image_ASC_Per_Object_")  # 

dbGetQuery(conn,"SELECT * FROM _link_tables_ASC_Per_Image_ASC_Per_Object_")  # 

dbGetQuery(conn,"SELECT * FROM sqlite_sequence")  # 




my_tableA02 <- read.csv("my_table00.csv")
head(my_tableA02)


umap.original <- nsclc.seurat.obj$umap@cell.embeddings 
str(umap.original)
dim(umap.original)


colnames(umap.original) <- c("UMAP_1ori", "UMAP_2ori")

clusterNum <- as.numeric(nsclc.seurat.obj@meta.data$RNA_snn_res.0.1)
clusterNum <- as.data.frame(clusterNum)
table(clusterNum)
clusterNum1 <- clusterNum-1
table(clusterNum1)
str(clusterNum)


# str(nsclc.seurat.obj)
dim(clusterNum)

colnames(clusterNum) <- c("cluster_RES01")

# umap@cell.embeddings 


my_tableA01 <- cbind(my_tableA02,clusterNum,umap.original)
head(my_tableA01)


my_old_table <- dbGetQuery(conn,"SELECT * FROM ASC_Per_Object")  
head(my_old_table)
head(my_tableA01)

dim(my_old_table)
dim(my_tableA01)


colnames(my_old_table)
colnames(my_tableA01)

dbRemoveTable(conn, "ASC_Per_Object")

dbWriteTable(conn, "ASC_Per_Object", my_tableA01)


# Close the database connection 
dbDisconnect(conn)





