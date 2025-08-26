


matrixReady <- read.csv(file="PseudoBulk86P.csv")
# matrixReady <- read.csv(file="PseudoBulk19P.csv")
# matrixReady <- read.csv(file="PseudoBulk41P.csv")


rownames(matrixReady)<- matrixReady[,1]
matrixReady <- matrixReady[,-1]



TrasposeReady <- t(matrixReady)
MatrixFinal <- t(TrasposeReady)

head(TrasposeReady)



tabele <-  scale(TrasposeReady)

#########################################################
# ACA EL CODIGO PARA HEATMAPS UNO AMPLIADO -con dendextend()

library(dendextend)
# ESTO TERMINA EN UN DENDROGRAMA ROWV
# order for rows
Rowv <- tabele %>% scale %>% dist %>% hclust %>% as.dendrogram %>%
set("branches_k_color", k = 3) %>% set("branches_lwd", 1.2) %>%
ladderize
heatmap(tabele, Rowv = Rowv, scale = "none")

str(Rowv)

# ESTO TERMINA EN UN DENDROGRAMA COLV
# Order for columns: We must transpose the data

Colv <- tabele %>% scale %>% t %>% dist %>% hclust %>% as.dendrogram %>%
set("branches_k_color", k = 3) %>%
set("branches_lwd", 1.2)

str(Colv )
leaf_labels <- labels(Colv )
print(leaf_labels)

# COMENTARIOS DE LOS PARAMETROS
# set("branches_k_color", k = 4) - ES PARA DEFINIR 3 COLORES - EL K DEFINE 3 CLUSTERS -  SI PONGO 4 PINTA 4 CLUSTERS
# set("branches_lwd", 1.2)  - DEFINE GROSOR DE LAS RAMAS -  CON 5.2 QUEDA GRUESO Y SE VE BIEN
# ladderize - ORDENA DE MENOR A MAYOR EN ESCALERA LAS RAMAS
# set("branches_k_color", k = 2, value = c("orange", "blue")) - DOS CLUSTERS PARA LAS COLUMNAS DOS COLORES PARA LAS RAMAS


heatmap(tabele, Rowv = Rowv, Colv = Colv,scale = "none",cexCol = 1.4 )


/////////////////////////
