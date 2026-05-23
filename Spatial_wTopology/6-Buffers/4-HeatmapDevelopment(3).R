

library(sf)
library(dplyr)
library(ggplot2)

# Load the data
cells_final <- readRDS("cells_analysis_sf_WithDistances.rds")
pc_polygons_dissolved <- readRDS("pc_polygons_dissolved_cleanFiltered100.rds")


# 1. Separate populations and calculate distances (Optimized)
active_cells <- cells_final %>% filter(class %in% c("PIB", "PL"))
sl_cells <- cells_final %>% filter(class == "SL")

# Find Nearest Neighbor Index
nearest_idx <- st_nearest_feature(sl_cells, active_cells)

# Calculate distances (Convert to microns: * 0.25)
sl_cells <- sl_cells %>%
  mutate(
    dist_to_active_um = as.numeric(st_distance(geometry, active_cells[nearest_idx, ], by_element = TRUE)) * 0.25
  )

# 2. CRITICAL STEP: Define the "Cell-Only" Window
# We grab the bounding box of the SL cells to crop the polygon layer
cell_bbox <- st_bbox(sl_cells)

#-----------------------------------------------
# 3. Improved Heatmap Plot
#-----------------------------------------------

interaction_map <- ggplot() +
  # Layer 1: PC Boundaries (Cropped to cell area)
  # We use coord_sf later to restrict the view
  geom_sf(data = pc_polygons_dissolved, 
          fill = NA, 
          color = "cyan", 
          linewidth = 0.3, 
          alpha = 0.4) +
  
  # Layer 2: The "Heat" - SL cells
  geom_sf(data = sl_cells, aes(color = dist_to_active_um), 
          size = 0.4, 
          alpha = 0.8) +
  
  # THE FIX: Restrict the coordinate system to the cell bounding box
  coord_sf(xlim = c(cell_bbox["xmin"], cell_bbox["xmax"]),
           ylim = c(cell_bbox["ymin"], cell_bbox["ymax"]),
           datum = NA) + 
  
  # Styling
  scale_color_viridis_c(
    option = "inferno", 
    direction = -1, 
    name = "Distance to \nActive Cell (µm)",
    limits = c(0, 30), # Focus on the immediate niche
    oob = scales::squish # Ensures values > 30 don't disappear (stay dark purple)
  ) +
  theme_void() + 
  theme(
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    legend.position = "right",
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    plot.title = element_text(color = "white", hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray80", hjust = 0.5)
  ) +
  labs(
    title = "CLL Interaction Heatmap",
    subtitle = "Restricted to Scanned Area (CellProfiler Coverage)"
  )

# Display
print(interaction_map)



#-----------------------------------------------
# 4. Custom Heatmap Plot - with centroids - fixing the dissolve problem
#-----------------------------------------------

# --- USER SETTINGS ---
Tcon <- 10    # Contact threshold (Red)
Tbio <- 12    # Biological/Paracrine threshold (Green)
# Note: Anything > Tbio will be categorized as Resting (Blue)

# Size settings for easy adjustment
size_contact  <-  0.8 # 1.4
size_paracrine <- 0.4 # 0.8
size_resting   <- 0.4 # 0.8
# ---------------------

# 1. Prepare Data and Classify
# We classify EVERY cell into one of the three tiers
niche_data <- sl_cells %>%
  mutate(
    dist_to_active_um = as.numeric(st_distance(geometry, active_cells[nearest_idx, ], by_element = TRUE)) * 0.25,
    niche_class = case_when(
      dist_to_active_um <= Tcon ~ "Contact",
      dist_to_active_um > Tcon & dist_to_active_um <= Tbio ~ "Paracrine",
      dist_to_active_um > Tbio ~ "Resting"
    )
  )

# Set the factor levels to control the plotting order (Resting on bottom, Contact on top)
niche_data$niche_class <- factor(niche_data$niche_class, levels = c("Resting", "Paracrine", "Contact"))

cell_bbox <- st_bbox(sl_cells)

# 2. Generate the Categorical Map
interaction_map_categories <- ggplot() +
  # Layer 1: PC Polygons (Light Plum Fill)
  geom_sf(data = pc_polygons_dissolved, 
          fill = "#FFBBFF",      
          color = "#551A8B",     
          linewidth = 0.5, 
          alpha = 0.4) +         
  
  # Layer 2: Plot the cells with sizes mapped to their class
  geom_sf(data = niche_data, aes(color = niche_class, size = niche_class), alpha = 0.8) +
  
  # Coordinate restriction to scanned area
  coord_sf(xlim = c(cell_bbox["xmin"], cell_bbox["xmax"]),
           ylim = c(cell_bbox["ymin"], cell_bbox["ymax"]),
           datum = NA) + 
  
  # Manual Color Scale: Contact=Red, Paracrine=Green, Resting=Blue
  scale_color_manual(
    values = c("Contact" = "#FF0000", "Paracrine" = "#7FFF00", "Resting" = "#FFD700"),
    name = "Niche Tier"
  ) +

  # Manual Size Scale: Using your user-defined variables
  scale_size_manual(
    values = c("Contact" = size_contact, "Paracrine" = size_paracrine, "Resting" = size_resting),
    name = "Niche Tier"
  ) +
  
  # Theme and Styling
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "#fdfdfd", color = NA),
    legend.position = "right",
    plot.title = element_text(color = "black", hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40", hjust = 0.5, size = 10)
  ) +
  labs(
    title = "A-CLL Spatial Niche Stratification",
    subtitle = paste("Contact (≤", Tcon, "µm) | Paracrine (", Tcon, "-", Tbio, "µm) | Resting (>", Tbio, "µm)")
  )

print(interaction_map_categories)



///////////////////////
Variantes para ver como marcar Polygons y cells superpuestas


# Size settings for easy adjustment
size_contact  <-  1.45 # 1.4
size_paracrine <- 0.85 # 0.8
size_resting   <- 0.7 # 0.8

# https://r-charts.com/colors/
# 2. Generate the Categorical Map
interaction_map_categories <- ggplot() +
  # Layer 1: PC Polygons (Light Plum Fill)
  geom_sf(data = pc_polygons_dissolved, 
          fill = "#5C5C5C",      
          color = "#030303",     
          linewidth = 0.8, 
          alpha = 0.9) +    # 0.4 
  
  # Layer 2: Plot the cells with sizes mapped to their class
  geom_sf(data = niche_data, aes(color = niche_class, size = niche_class), alpha = 0.15) +
  
  # Coordinate restriction to scanned area
  coord_sf(xlim = c(cell_bbox["xmin"], cell_bbox["xmax"]),
           ylim = c(cell_bbox["ymin"], cell_bbox["ymax"]),
           datum = NA) + 
  
  # Manual Color Scale: Contact=Red, Paracrine=Green, Resting=Blue
  scale_color_manual(
    values = c("Contact" = "#FF0000", "Paracrine" = "#00FF00", "Resting" = "#0000FF"),
    name = "Niche Tier"
  ) +

  # Manual Size Scale: Using your user-defined variables
  scale_size_manual(
    values = c("Contact" = size_contact, "Paracrine" = size_paracrine, "Resting" = size_resting),
    name = "Niche Tier"
  ) +
  
  # Theme and Styling
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "#fdfdfd", color = NA),
    legend.position = "right",
    plot.title = element_text(color = "black", hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40", hjust = 0.5, size = 10)
  ) +
  labs(
    title = "A-CLL Spatial Niche Stratification",
    subtitle = paste("Contact (≤", Tcon, "µm) | Paracrine (", Tcon, "-", Tbio, "µm) | Resting (>", Tbio, "µm)")
  )

print(interaction_map_categories)


