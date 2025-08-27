
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)


setwd("D:/EVOSdata/Ceci/Desarrollo2025/temp/RepasoPIBparaML")
getwd()


OBJ <- readRDS("CProSeuratOBJNEG_04.rds")
# str(OBJ@assays$RNA@layers)
# str(OBJ@assays$RNA@layers$counts)


Sparemat <- OBJ@assays$RNA@layers$counts
mat <- as.matrix(Sparemat)
mat <- t(mat)
mat <- as.data.frame(mat)

dim(mat)

rownames(mat)<-  colnames(OBJ@assays$RNA)
colnames(mat)<- rownames(OBJ@assays$RNA)

colnames(mat)

ML <- read.csv("my_tableClass_04Paciente4.csv")
colnames(ML)




###-----------------------------------------------------------------------------
### PLOTEO CLUSTER CON LEYENDA  - USAR ESTE PARA PLOTEO
###-----------------------------------------------------------------------------

ML1 <- as.data.frame(ML[,"class"])
colnames(ML1)<- c("CellType")
dim(ML1)
table(ML1[,1])

dat1 <- cbind(mat,ML1)
colnames(dat1)

table(dat1$"ImageNumber")


dat2 <- dat1[dat1$"ImageNumber"== 10, ]
dim(dat1)
dim(dat2)



# PLOTEO CON TAMAÑO Y COLOR VARIABLE SEGUN TIPO CELULAR
sp <- ggplot(dat2, aes(x=dat2$"NucObjects-AreaShape-Center-X" , 
y= -(dat2$"NucObjects-AreaShape-Center-Y"), color=dat2$"CellType",size=dat2$"CellType" ))+
scale_x_continuous(limits = c(0, 2040))+
geom_point()+
# scale_colour_brewer(type = 'qual', palette = "Dark2")
scale_size_manual(values = c("PIB" = 8, "PL" =6, "SL" = 3))+
scale_color_manual(values = c("PIB" = "red", "PL" = "orange", "SL"="blue"))

sp

# plotPNG("Thomas",sp)
#-----------------------------------------------------------
# Funcion ploteo individual en png
plotPNG <- function(name,p){
WithExt <- paste(name,".png",sep="")
png(WithExt,width = 780, height = 480,)
 print(p)
dev.off()
}
#-----------------------------------------------------------
# FIN Funcion ploteo individual en png




#-------- Seccion imprimir a png en bloque de un paciente
# usa la función plotPNG

nombrePac <- "Paciente4"
numberOfimages <- max(dat1$"ImageNumber")
dim(dat1)
colnames(dat1)

# i=1
for (i in 1:numberOfimages) {
name <- paste(nombrePac ,"_",i,sep="")
# print(name)
dat2 <- dat1[dat1$"ImageNumber"== i, ]
dim(dat2)

sp <- ggplot(dat2, aes(x=dat2$"NucObjects-AreaShape-Center-X" , 
y= -(dat2$"NucObjects-AreaShape-Center-Y"), color=dat2$"CellType",size=dat2$"CellType" ))+
scale_x_continuous(limits = c(0, 2040))+
geom_point()+
# scale_colour_brewer(type = 'qual', palette = "Dark2")
scale_size_manual(values = c("PIB" = 8, "PL" =6, "SL" = 3))+
scale_color_manual(values = c("PIB" = "red", "PL" = "orange", "SL"="blue"))

sp

plotPNG(name,sp)
}

#--------Fin Seccion imprimir a png en bloque de un paciente