


CproINum <- read.csv("my_IMAGES_tableP01.csv")
colnames(CproINum)
# [1] "ImageNumber" 
# [2] "Image_FileName_HE"  



CproINum <- CproINum[,1:2]
colnames(CproINum) 
  
colnames(CproINum) <- c("ImageNumber","full_filename")
# > colnames(CproINum) 
# [1] "ImageNumber"   "full_filename"
# > 

write.csv(CproINum ,"CproINum_To_FileName.csv" ) 

#-------------------------------------------------
Atención: hay que modificar en el csv final 3 cosas:

d4.tif a d0.TIF

_R_  a _M_

eliminar clara y oscura que estan repetidas. Los numero  comienzan en 3 hasta 28
son 26 totales