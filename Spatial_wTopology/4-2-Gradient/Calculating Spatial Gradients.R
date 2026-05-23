



# input1 cells_analysis_sf.rds
# input2 pc_polygons_dissolved_cleanFiltered100.rds
# output cells_analysis_sf_WithDistances.rds



library(sf)
library(dplyr)
pc_polygons_dissolved <- readRDS("pc_polygons_dissolved_cleanFiltered100.rds")

str(pc_polygons_dissolved)

#---------------------------------------------------------
# R Code: Calculating Spatial Gradients
# 
#----------------------------------------------------------
# 1. Extract the Boundary (the perimeter) of the PC Niche
pc_boundary <- st_cast(pc_polygons_dissolved, "MULTILINESTRING")

cells_final <-  readRDS("cells_analysis_sf.rds")

# 2. Calculate the distance from every cell to the boundary
# This returns a 'units' object (pixels)
distances <- st_distance(cells_final, pc_boundary)


# 3. Add to our dataframe
# We convert to numeric to make it easier for plotting
cells_final <- cells_final %>%
  mutate(dist_to_PC_edge = as.numeric(distances))

colnames(cells_final)
# 4. Optional: Create 'Inside vs Outside' distance logic
# (Inside = negative, Outside = positive)
cells_final <- cells_final %>%
  mutate(dist_signed = ifelse(In_PC == TRUE, 
                              -dist_to_PC_edge, 
                               dist_to_PC_edge))

# 5. Check the distribution
summary(cells_final$dist_signed)


saveRDS(cells_final,"cells_analysis_sf_WithDistances.rds")

# _______________________________________
# The Logic
# We will use st_distance() to calculate the shortest distance 
# from every cell centroid to the boundary (perimeter) of your Proliferation Center polygons.

# distances <- st_distance(cells_final, pc_boundary)

# •	Cells outside the PC will have a positive distance (distance to reach the PC).
# •	Cells inside the PC will have a distance of 0 (by default).
# •	Advanced Tip: We can calculate "negative" distances for cells inside the PC 
# to represent how deep they are from the border.



library(ggplot2)

ggplot(cells_final, aes(x = dist_signed, fill = In_PC)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "CLL Cell Distribution relative to PC Border",
       x = "Distance to Border (pixels) [Negative = Inside PC]",
       y = "Cell Density")


#--------- In microns

# Convert the signed distance to microns for biological relevance
cells_final <- cells_final %>%
  mutate(dist_signed_um = dist_signed * 0.25)

# Plotting the Density with Microns
ggplot(cells_final, aes(x = dist_signed_um, fill = In_PC)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "CLL Spatial Compartmentalization",
       subtitle = "Negative = Inside PC Niche | Positive = Interfollicular Space",
       x = "Distance to PC Border (µm)",
       y = "Density of Cells")

#------------------------------
