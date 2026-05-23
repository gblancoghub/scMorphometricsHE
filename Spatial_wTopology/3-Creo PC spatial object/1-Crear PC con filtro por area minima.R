



library(terra)
library(sf)
library(dplyr)


# 1. Load the small probability map
# prob_map <- rast("scan_Bottom Slide_TM1_p00_0_A01f00d0-25_Probabilities_1.tif")

prob_map <- rast("scan_Bottom Slide_TM1_p00_0_A01f00d0-25_Probabilities_1.tif")

# prob_map <- rast("scan_Bottom Slide_TM1_ProbabilitiesMODIFIED.tif")

# 2. Vectorize the SMALL raster (High speed, Low RAM)

pc_mask <- prob_map > 0.6
names(pc_mask) <- "class_value"

# This creates polygons in the small coordinate space (0 to 2335)
pc_polygons_small <- as.polygons(pc_mask, values = TRUE) %>% 
                     st_as_sf() %>%
                     filter(class_value == 1)

# 3. Calculate scaling factors
global_width  <- 26624 
global_height <- 30720 
scale_x <- global_width / ncol(prob_map)
scale_y <- global_height / nrow(prob_map)

# 4. Scale the Geometry directly (Mathematically, not graphically)
# We treat the geometry as a simple list of coordinates to avoid sf errors
scaled_geometry <- st_geometry(pc_polygons_small) * c(scale_x, scale_y)

# 5. Create the final object by combining the data and the new geometry
pc_polygons_global <- st_sf(
  class_value = pc_polygons_small$class_value, 
  geometry = scaled_geometry
)

# 6. Verify with a quick plot of just the polygons
plot(st_geometry(pc_polygons_global), main="Global PC Polygons")

# 1. Take your global polygons (before the union)
# Or if you only have the dissolved one, use st_cast to break it apart:
pc_islands <- st_cast(pc_polygons_global, "POLYGON")

# 2. Calculate the area of each individual island
# Note: Since CRS is NA, this will be in 'square pixels'
pc_islands <- pc_islands %>%
  mutate(area_px = as.numeric(st_area(.)))

dim(pc_islands)


# threshold_area <- 10000 
threshold_area <-     # 160000 


pc_clean <- pc_islands %>%
  filter(area_px > threshold_area)


saveRDS( pc_clean, "pc_polygonsFiltered100.rds" )

pc_clean <- readRDS("pc_polygonsFiltered100.rds") 

# 4. (Optional) Re-dissolve the remaining large polygons
 pc_polygons_dissolved_clean <- pc_clean %>%
   st_union() %>%
   st_as_sf() %>%
   mutate(In_PC = TRUE)

saveRDS( pc_polygons_dissolved_clean, "pc_polygons_dissolved_cleanFiltered100.rds" )

# 5. Check the difference
cat("Original island count:", nrow(pc_islands), "\n")
cat("Filtered island count:", nrow(pc_clean), "\n")




