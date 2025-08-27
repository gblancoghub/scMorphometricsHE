# A simplified code for building and plotting  trjacetories for each patients
# Includes options to customize plotting of trajectories and dots representing cells (transparency)
# Includes reporting nummber of roots (expansion appears relevant to a-CLL cases)
# may be linked to genomic instability

set.seed(1234)

library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)




getwd()

# Restore the object
# El que tiene ya los clusters generado en carpeta ScType

# DEL11q.integrated <- readRDS(file = "merged_seurat_integrated.rds")
  DEL11q.integrated <- readRDS(file = "seurat_integrated.rds")

View(DEL11q.integrated@meta.data)
colnames(DEL11q.integrated@meta.data)
DEL11q.integrated@meta.data$seurat_clusters


# backup <- DEL11q.integrated
DEL11q.integrated <- backup 



# DEL11q.integrated <- subset(DEL11q.integrated, subset = (Patient == "01_PAC02"))

 DEL11q.integrated <- subset(DEL11q.integrated, subset = (Patient == "00_PAC01"))

#  DEL11q.integrated <- subset(DEL11q.integrated, subset = (Patient == "02_PAC03"))
#
# DEL11q.integrated <- subset(DEL11q.integrated, subset = (Patient == "04_PAC04"))

cds <- as.cell_data_set(DEL11q.integrated)
dff <- colData(cds)
fData(cds)$gene_short_name <- rownames(fData(cds))
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)
cds@clusters$UMAP$partitions <- reacreate.partition
list_cluster <- DEL11q.integrated@active.ident
cds@clusters$UMAP$clusters <- list_cluster
cds@int_colData@listData$reducedDims$UMAP <- DEL11q.integrated@reductions$umap@cell.embeddings

cluster.before.trajectory <- plot_cells(cds,
           color_cells_by = 'seurat_clusters',
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  theme(legend.position = "right")
cluster.before.trajectory



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

cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == c(13)]))


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


ROOTS
# str(cds@principal_graph_aux)
# cds@principal_graph_aux$UMAP$"root_pr_nodes"
str(cds@principal_graph_aux$UMAP$"root_pr_nodes")



library(ggrepel)
plot_cells(cds,
      #     color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_cell_groups = FALSE, # added label_cell_groups
           label_branch_points = TRUE,
           label_roots = FALSE,
           label_leaves = FALSE,
           cell_size = 0)


           cell_size = 0.00000000000001)




