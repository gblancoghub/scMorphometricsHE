
# install.packages("sf")
library(sf)
library(dplyr)




#------ Actualizar la altura del Tile para futuro flipping
#

h <- 30720 

#--------------------------------------------------------------
cells_global <- read.csv("CLL_Global_Spatial_Data.csv")
colnames(cells_global)


# 3. Convert to a Spatial Object (GIS Point Layer)
cells_sf <- st_as_sf(cells_global, coords = c("xGlobUn", "yGlobUn"))


# xGlob and yGlob are your 'Ground Truth'
# In your H&E plot, you used (xGlob) and (-yGlob)
# To match the Polygons (which are 0 to 30720), we need:
# NewX = xGlob
# NewY = (yGlob * -1) + 30720  <-- This flips the origin to the top

coords_matrix <- matrix(c(
  cells_sf$xGlob, 
  (cells_sf$yGlob * -1) + h), ncol = 2)

# 2. Rebuild the geometry column from scratch
st_geometry(cells_sf) <- st_sfc(lapply(1:nrow(coords_matrix), function(i) {
  st_point(coords_matrix[i, ])
}))

# 3. Final Bounding Box Check
print(st_bbox(cells_sf))

> print(st_bbox(cells_sf))
       xmin        ymin        xmax        ymax 
10267.02677    18.92203 26602.25274 15342.99580 
>

global_h <- 30720 # Your known max height

saveRDS(cells_sf, "cells_spatial_object.rds")


colnames(cells_sf)

cells_sf <- readRDS( "cells_spatial_object.rds")
write.csv(cells_sf , "control.csv")

