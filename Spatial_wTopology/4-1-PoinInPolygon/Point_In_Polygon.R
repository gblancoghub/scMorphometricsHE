



# input1 - "cells_spatial_object.rds"
# input2 - "pc_polygonsFiltered100.rds" o mejor el disuelto pc_polygons_dissolved_cleanFiltered100.rds

# output - "cells_analysis_sf.rds"  tiene el resultado de la operacion punto en poligono


library(sf)
library(dplyr)

cells_analysis_sf <- readRDS("cells_spatial_object.rds")
pc_polygons_dissolved_clean <- readRDS("pc_polygons_dissolved_cleanFiltered100.rds")

dim(cells_analysis_sf )

#-------- Este no lo use  por que no esta disuelto y cada cell hace overlay con mas de un poligono
# tipo en un tiro al blanco con poligonos concentricos
# 1. Point-in-Polygon: Assign each cell to a PC ID
# left = TRUE keeps all cells; cells not in a PC will have NA in PC_ID
# --------------------------------------------------------------------------------------------

# Use the dissolved version for the True/False flag
cells_analysis_sf <- cells_analysis_sf %>%
  mutate(In_PC = lengths(st_intersects(., pc_polygons_dissolved_clean)) > 0)
dim(cells_analysis_sf)


cells_analysis_sf <- st_join(cells_analysis_sf, 
                             pc_polygons_dissolved_clean %>% mutate(PC_ID = row_number()), 
                             join = st_intersects)

# 2. Update the 'In_PC' logical flag
cells_analysis_sf <- cells_analysis_sf %>%
  mutate(In_PC = !is.na(PC_ID))

# 3. Re-calculate the Distance to the nearest *Clean* PC Centroid
pc_centroids_clean <- st_centroid(pc_polygons_dissolved_clean)
dist_matrix_final <- st_distance(cells_analysis_sf, pc_centroids_clean)
cells_analysis_sf$dist_to_pc_center <- apply(dist_matrix_final, 1, min)

# 4. Convert to microns
cells_analysis_sf$dist_um <- as.numeric(cells_analysis_sf$dist_to_pc_center) * 0.25

saveRDS(cells_analysis_sf, "cells_analysis_sf.rds")
write.csv(cells_analysis_sf, "cells con overlay sobre PC layer.csv")
dim(cells_analysis_sf)
