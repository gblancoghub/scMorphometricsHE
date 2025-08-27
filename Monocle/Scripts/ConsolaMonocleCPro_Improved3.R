
# A basic script for building trajectories for each patient
# with plotting and exploring effects of selection of initital clusters
# There is also an intersting box plot of clusters all along trajectories

library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)




# Restore the object
# 

# DEL11q.integrated <- readRDS(file = "merged_seurat_integrated.rds")
  DEL11q.integrated <- readRDS(file = "seurat_integrated.rds")

View(DEL11q.integrated@meta.data)
colnames(DEL11q.integrated@meta.data)
DEL11q.integrated@meta.data$seurat_clusters
DEL11q.integrated@meta.data$RNA_snn_res.0.1 


# Process pateints one by one by subsetting

# DEL11q.integrated <- subset(DEL11q.integrated, subset = (Patient == "PAC02"))

# DEL11q.integrated <- subset(DEL11q.integrated, subset = (Patient == "PAC01"))

# DEL11q.integrated <- subset(DEL11q.integrated, subset = (Patient == "PAC03"))

# DEL11q.integrated <- subset(DEL11q.integrated, subset = (Patient == "PAC04"))

cds <- as.cell_data_set(DEL11q.integrated)


### Exploración del objeto cds

colData(cds)
rownames(colData(cds)) # los barcodes
colnames(colData(cds))


dff <- colData(cds)

fData(cds) #metadata


rownames(fData(cds))[1:10]
# tiene los symbol en rownames y no tiene ninguna columna
> rownames(fData(cds))[1:10]
 [1] "CellObjects-Mean-Nucelolos-AreaShape-MajorAxisLength"       
 [2] "CellObjects-Mean-Nucelolos-AreaShape-MaxFeretDiameter"      
 [3] "CellObjects-Mean-Nucelolos-Number-Object-Number"            
 [4] "CellObjects-Mean-Nucelolos-Distance-Centroid-CellObjects"   
 [5] "CellObjects-Mean-Nucelolos-AreaShape-EquivalentDiameter"    
 [6] "CellObjects-Mean-Nucelolos-Texture-Variance-HEGray-3-02-256"
 [7] "CellObjects-Mean-Nucelolos-AreaShape-Area"                  
 [8] "CellObjects-Mean-Nucelolos-Texture-Variance-HEGray-3-00-256"
 [9] "CellObjects-Mean-Nucelolos-Texture-Variance-HEGray-3-01-256"
[10] "CellObjects-Mean-Nucelolos-AreaShape-ConvexArea"            
> 

# la crea con los rownames
# ESTA COLUMNA  gene_short_name LA NCESITA MONOCLE
# since it misses the gene_short_name column, let's add it

fData(cds)$gene_short_name <- rownames(fData(cds))


# to get counts
counts(cds)

# PARA PASAR EL CLUSTERING DE SEURAT A MONOCLE HAY QUE REPRODUCIR EL CLSUTERING DE MONOCLE QUE NO ES IGUAL A SEURAT
# MONOCLE USA METACLUSTERING QUE LA LLAMA PARTITION
# POR ESO CREAMOS UNA PARTITION UNICA QUE CONTIEN A TODOS LOS CLUSTERS


# LE TENGO QUE PASAR DE CADA CELULA:
# 1- PARTITION
# 2- CLUSTERS
# 3-COORDENADAS UMAP (embeddings)


# 1- PARTITION
# A assign paritions
cds@colData@rownames # todas las cells

# B  armo una lista llena de unos, tan larga como la cantidad de cells
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames

# C ES UNA LISTA NOMINADA O NOMBRADA
reacreate.partition <- as.factor(reacreate.partition)

# Con eso tengo la partición

cds@clusters$UMAP$partitions <- reacreate.partition


2- CLUSTERS
# Assign the cluster info 

list_cluster <- DEL11q.integrated@active.ident
cds@clusters$UMAP$clusters <- list_cluster


# 3-COORDENADAS UMAP (embeddings)
# Assign UMAP coordinate - cell embeddings

cds@int_colData@listData$reducedDims$UMAP <- DEL11q.integrated@reductions$umap@cell.embeddings


# PLOTEO CON MONOCLE


cluster.before.trajectory <- plot_cells(cds,
           color_cells_by = 'seurat_clusters',
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  theme(legend.position = "right")



# ...3. Learn trajectory graph ------------------------
cds <- learn_graph(cds, use_partition = FALSE)

# MONOCLE HACE UNA TRAJECTORIA RAMIFICADA SEPARADA PARA CADA PARTITION
# PERO COMO ACÁ ES UNA SOLA DIRECTAMENTE LE PONEMOS FALSE


# DESPUES DEL LEARN HAGO UN PLOT
plot_cells(cds,
           color_cells_by = 'seurat_clusters',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)+
  theme(legend.position = "right")






# ...4. Order the cells in pseudotime -------------------

cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == c(9)]))


plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)

# root+branch+leaves
plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = TRUE,
           label_roots = TRUE,
           label_leaves = TRUE,
           group_label_size = 5) +
		theme(legend.position = "right")

# solo root
plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = TRUE,
           label_leaves = FALSE,
           group_label_size = 5) +
		theme(legend.position = "right")
# solo leaves
plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = TRUE,
           group_label_size = 5) +
		theme(legend.position = "right")

# solo branches
plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = TRUE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5) +
		theme(legend.position = "right")



# cells ordered by monocle3 pseudotime

# pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = seurat_clusters)) +
  geom_boxplot()










