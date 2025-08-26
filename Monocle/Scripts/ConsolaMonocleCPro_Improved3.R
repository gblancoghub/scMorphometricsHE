

set.seed(1234)

library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)



setwd("D:/EVOSdata/Ceci/Desarrollo2025/CPro(2)/A-Monocle")
getwd()

# Restore the object
# El que tiene ya los clusters generado en carpeta ScType

# DEL11q.integrated <- readRDS(file = "merged_seurat_integrated.rds")
  DEL11q.integrated <- readRDS(file = "seurat_integrated(AIplus).rds")

View(DEL11q.integrated@meta.data)
colnames(DEL11q.integrated@meta.data)
DEL11q.integrated@meta.data$seurat_clusters
DEL11q.integrated@meta.data$RNA_snn_res.0.1 


Idents(DEL11q.integrated) 
# Idents(DEL11q.integrated) <- "RNA_snn_res.0.1"

DEL11q.integrated@meta.data$seurat_clusters 
DEL11q.integrated@meta.data$RNA_snn_res.0.1


 Idents(DEL11q.integrated) <- "RNA_snn_res.0.1"
DimPlot(DEL11q.integrated , reduction = 'umap', group.by = "RNA_snn_res.0.3", label = TRUE)
subOBJ <- subset(DEL11q.integrated, idents = c("0"))
DimPlot(DEL11q.integrated , reduction = 'umap', label = TRUE)

DimPlot(subOBJ  , reduction = 'umap', split.by = "Patient", label = TRUE)


# backup <- DEL11q.integrated
DEL11q.integrated <- backup 

# DEL11q.integrated <- subset(DEL11q.integrated, subset = (tissue=="LN"))
# DEL11q.integrated <- subset(DEL11q.integrated, idents = c("0","1","2","3","6"))

# DEL11q.integrated <- subset(DEL11q.integrated, subset = (Patient == "Macchia"))
# DEL11q.integrated <- subset(DEL11q.integrated, idents = c("0","2","5"))

# DEL11q.integrated <- subset(DEL11q.integrated, subset = (Patient == "Rosetto"))
# DEL11q.integrated <- subset(DEL11q.integrated, idents = c("0","2","3","5"))

# DEL11q.integrated <- subset(DEL11q.integrated, subset = (Patient == "Thomas"))
# DEL11q.integrated <- subset(DEL11q.integrated, idents = c("0","2","5"))

# DEL11q.integrated <- subset(DEL11q.integrated, subset = (Patient == "Paciente4"))
# DEL11q.integrated <- subset(DEL11q.integrated, idents = c("1","4","5","6"))

cds <- as.cell_data_set(DEL11q.integrated)


### Exploración del objeto cds

colData(cds)
rownames(colData(cds)) # los barcodes
colnames(colData(cds))
> colnames(colData(cds))
 [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "RNA_snn_res.0.1"
 [5] "RNA_snn_res.0.3" "RNA_snn_res.0.5" "RNA_snn_res.0.7" "RNA_snn_res.1"  
 [9] "seurat_clusters" "num1"            "num2"            "Patient"        
[13] "Barcode"         "ident"          
>
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

# plot  -plot_cells parece de Monocle3
# ACA EL COLOR ES POR EL SEURAT CLUSTER

cluster.before.trajectory <- plot_cells(cds,
           color_cells_by = 'seurat_clusters',
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  theme(legend.position = "right")


# ACA EL COLOR ES POR EL CUSTOM CLASSIFICATION
cluster.names <- plot_cells(cds,
           color_cells_by = "seurat_clusters",
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  scale_color_manual(values = c('red', 'blue', 'green', 'maroon', 'yellow', 'grey', 'cyan','black','pink','orange','magenta')) +
  theme(legend.position = "right")

cluster.before.trajectory | cluster.names

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


# ...5. Finding genes that change as a function of pseudotime --------------------
deg_bcells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

head(deg_bcells)
> head(deg_bcells)
           status   p_value morans_test_statistic      morans_I
OR4F5        FAIL        NA                    NA            NA
FO538757.3     OK       NaN                   NaN -0.0005900174
FO538757.2     OK 0.8115361            -0.8835711 -0.0077445452
OR4F29       FAIL        NA                    NA            NA
OR4F16         OK       NaN                   NaN -0.0008526141
SAMD11         OK 0.6718432            -0.4450086 -0.0043132580
           gene_short_name   q_value
OR4F5                OR4F5 1.0000000
FO538757.3      FO538757.3       NaN
FO538757.2      FO538757.2 0.8337307
OR4F29              OR4F29 1.0000000
OR4F16              OR4F16       NaN
SAMD11              SAMD11 0.7403549



Los_genes <- deg_bcells %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') 

write.csv(Los_genes,"DEgenes_Monocle.csv")

b.seu <- DEL11q.integrated
FeaturePlot(b.seu, features = c('ENO1', 'STMN1', 'CDCA8'))
FeaturePlot(b.seu, features = c('YBX1', 'CDC20', 'PRDX1'))

FeaturePlot(b.seu, features = c('CXCR4', 'CD69', 'B2M', 'HMGB2'))
FeaturePlot(b.seu, features = c('AURKB', 'BIRC5', 'PCNA', 'UBE2C'))
FeaturePlot(b.seu, features = c('MYC', 'CD37', 'CD2', 'KLF2'))
FeaturePlot(b.seu, features = c('CD19', 'CBX2', 'CBX4', 'CCL4'))

FeaturePlot(b.seu, features = c('MKI67', 'TOP2A', 'GAPDH','YBX1'))
HIST1H4C
LDHA
JUNB
RPL29
FeaturePlot(b.seu, features = c('HIST1H4C', 'LDHA', 'JUNB','CDC20'))




FeaturePlot(b.seu, features = c('ENO1', 'CXCR4'),split.by ="tissue", label = T)

1 3 5 6 
7 - 2 4 - 8


# visualizing pseudotime in seurat

DEL11q.integrated$pseudotime <- pseudotime(cds)
# Idents(DEL11q.integrated) <- DEL11q.integrated$seurat_clusters
FeaturePlot(DEL11q.integrated, features = "pseudotime", label = T)


COMIENZA SECCION QUE NO ESTA EN EL TP DE K PATTEL - ES DEL TUTORIAL
PASO 1 - PREPROCESO DEL OBJETO CDS

cds2 <- preprocess_cds(cds, num_dim = 50)
# cds <- align_cds(cds, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")

cds2 <- reduce_dimension(cds2)
plot_cells(cds2, label_groups_by_cluster=FALSE,  color_cells_by = "seurat_clusters")

cds2 <- cluster_cells(cds2)
plot_cells(cds2, color_cells_by = "seurat_clusters",group_label_size=5,min_expr =500 )
?plot_cells
?cluster_cells
cds2 <- cluster_cells(cds2)

PASO 2 LEARNING DEL NUEVO CDS2
cds2 <- learn_graph(cds2)

plot_cells(cds2,
           color_cells_by = "seurat_clusters",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

# RECORDAR QUE ES UNA APLICACION SHINY Y HAY QUE IR AL BROWSER!!!
cds2 <- order_cells(cds2)
# order_cells(
#   cds2,
#   reduction_method = "UMAP",
#   root_pr_nodes = NULL,
#   root_cells = NULL,
#   verbose = FALSE
# )
?order_cells
plot_cells(cds2,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

3- GENES DEG CON GRAPH_TEST

# deg_bcells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)
pr_graph_test_res <- graph_test(cds2, neighbor_graph="knn", cores=8)

pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
head(pr_deg_ids)




4-FIND MODULES - AGRUPA GENES EN MODULOS - GRABO CSV - SOLO GENE Y MODULO - MISMA PARTICION

gene_module_df <- find_gene_modules(cds2[pr_deg_ids,], resolution=1e-2)

write.csv(gene_module_df,"gene_module.csv" )

cds2
dff[,"seurat_clusters"] 

5- HAGO AGRUPAMIENTO POR CLUSTER SEURAT
cell_group_df <- tibble::tibble(cell=row.names(colData(cds2)), 
                                cell_group=dff[,"seurat_clusters"] )

6- USO LA FUNCION AGGREGATE GENE EXPRESION QUE COMBIAN CLUSTER CON MODULE - GRABO CSV
TENGO 24 MODULES Y VEO EN QUE CLUSTER VA CADA MODULO

agg_mat <- aggregate_gene_expression(cds2, gene_module_df, cell_group_df)
head(agg_mat)
write.csv(agg_mat,"cluster_module.csv")

row.names(agg_mat) <- stringr::str_c("", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Cluster ", colnames(agg_mat))

write.csv(agg_mat,"cluster_moduleWithLabels.csv")
pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)


plot_cells(cds2, 
           genes=gene_module_df %>% filter(module %in% c(10, 4, 1, 11)),
           group_cells_by="cluster",
           color_cells_by="cluster",
           show_trajectory_graph=FALSE)

plot_cells(cds2, 
           genes=gene_module_df %>% filter(module %in% c(16, 3, 1, 9)),
           group_cells_by="cluster",
           color_cells_by="cluster",
           show_trajectory_graph=FALSE)



7-PLOTEO UN FUNCION EXPRESION VS TIEMPO PARA VER COMO SE MUEVE CADA GEN 
-EJ ENO1 AL COMIENZO CD69 Y CXCR4/CD69 LEVANTAN AL FINAL 
ANDA FEO - PERO LOS DATOS ESTAN

AFD_genes <- c("F3", "CD55", "ICAM1")
AFD_genes <- c("PKM", "TGFB1", "BNIP3")
AFD_lineage_cds2 <- cds2[rowData(cds2)$gene_short_name %in% AFD_genes,
                       colData(cds2)$seurat_clusters %in% c("0","1","2","3","4")]
# c("0","1","2","3","4","6","8")

# AFD_lineage_cds2 <- cds2[rowData(cds2)$gene_short_name %in% AFD_genes,]


The function plot_genes_in_pseudotime() takes a small set of genes and shows 
you their dynamics as a function of pseudotime: 

plot_genes_in_pseudotime(AFD_lineage_cds2,
                         color_cells_by="monocle3_pseudotime",vertical_jitter = NULL)+
ylim(0,1)


?plot_genes_in_pseudotime

trend_formula = "~ splines::ns(pseudotime, df=3)",
min_expr = NULL,
  cell_size = 0.75,
plot_genes_in_pseudotime(AFD_lineage_cds2,
                         color_cells_by="pseudotime",trend_formula = "~ splines::ns(pseudotime, df=3)",
min_expr =1e-5, cell_size = 0.75,  ncol = 1 )



colData(AFD_lineage_cds2)

? plot_genes_in_pseudotime()

getwd()














cds_sub <- choose_graph_segments(cds) # va a un shiny en browser Mozilla

plot_cells(cds_sub ,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = TRUE,
           label_leaves = FALSE,
           group_label_size = 5) +
		theme(legend.position = "right")

# gene_fits <- fit_models(cds_subset, model_formula_str = "~pseudotime") # no termina mas esto




