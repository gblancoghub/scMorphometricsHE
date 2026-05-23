
# Version con shap factor 0.744 en lugar de 0.8

library(sf)
library(dplyr)
library(ggplot2)

# Load the data
cells_final <- readRDS("cells_analysis_sf_WithDistances.rds")

# 1. Separate populations
sl_cells <- cells_final %>% filter(class == "SL")
active_cells <- cells_final %>% filter(class %in% c("PIB", "PL"))

# 2. Optimized Nearest Neighbor Calculation (Avoids 25GB Matrix)
# Find the index of the closest active cell
nearest_idx <- st_nearest_feature(sl_cells, active_cells)

# Calculate distances only for those specific pairs
sl_cells$dist_to_nearest_active <- st_distance(
  sl_cells, 
  active_cells[nearest_idx, ], 
  by_element = TRUE
)

# 3. Apply the 3-Tier Stratification
sl_cells <- sl_cells %>%
  mutate(
    dist_to_active_um = as.numeric(dist_to_nearest_active) * 0.25,
    # New stratification logic
    sl_subclass = case_when(
      dist_to_active_um < 10.0              ~ "SL_contact",
      dist_to_active_um >= 10.0 & 
        dist_to_active_um <= 12.0156          ~ "SL_paracrine",
      dist_to_active_um > 12.0156             ~ "SL_resting"
    )
  )





# View the breakdown for the meeting report
print(table(sl_cells$sl_subclass))

# 4. Join status back to cells_final
sl_lookup <- sl_cells %>%
  st_drop_geometry() %>%
  select(ImageNumber, ObjectNumber, sl_subclass)

cells_final <- cells_final %>%
  left_join(sl_lookup, by = c("ImageNumber", "ObjectNumber")) %>%
  mutate(final_class = ifelse(class == "SL", sl_subclass, class))

# 5. Quick Summary Table for Hussein
summary_table <- cells_final %>%
  st_drop_geometry() %>%
  group_by(final_class) %>%
  summarise(
    Count = n(),
    Percentage = (n() / nrow(cells_final)) * 100
  )

print(summary_table)

# Save the updated object
saveRDS(cells_final, "cells_final_WithProximityStratification_3Tier744.rds")

# cells_final <- readRDS("cells_final_WithProximityStratification_3Tier.rds")

#-----------------------------------------------------------------------------------------
# 2. Plot the Fraction of SL Contact, Paracrine, and Resting along the Spatial Gradient
#-----------------------------------------------------------------------------------------

# 1. Prepare binned data
cells_binned <- cells_final %>%
  st_drop_geometry() %>%
  mutate(dist_signed_um = dist_signed * 0.25) %>%
  # Focus on a window of 150um around the PC boundary
  filter(dist_signed_um > -150 & dist_signed_um < 150) %>%
  mutate(bin = cut(dist_signed_um, breaks = seq(-150, 150, by = 10)))

# 2. Calculate proportions specifically for the 3 SL tiers
composition_sl <- cells_binned %>%
  # We filter for our three new subclasses created in the previous step
  filter(final_class %in% c("SL_contact", "SL_paracrine", "SL_resting")) %>% 
  group_by(bin, final_class) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(
    midpoint = (as.numeric(sub("\\((.+),(.+)]", "\\1", bin)) + 
                as.numeric(sub("\\((.+),(.+)]", "\\2", bin))) / 2,
    proportion = (count / sum(count)) * 100
  ) %>%
  ungroup()

# 3. Plotting the 3-Tier Distribution
ggplot(composition_sl, aes(x = midpoint, y = proportion, color = final_class)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  # Boundary line: Positive distances are outside the PC, negative are inside
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  # Using a high-contrast palette for the 3 tiers
  scale_color_manual(values = c(
    "SL_contact" = "#E41A1C",   # Red for high-intensity contact
    "SL_paracrine" = "#4DAF4A",  # Green for paracrine influence
    "SL_resting" = "#377EB8"     # Blue for resting/distant cells
  )) +
  theme_minimal() +
  labs(
    title = "SL Sub-population Spatial Distribution",
    subtitle = "Contact (<10µm), Paracrine (10-12µm), and Resting (>12µm) vs PC boundary",
    x = "Distance to PC Border (µm) [Inside < 0 > Outside]",
    y = "Percentage of SL Population (%)",
    color = "SL Tier"
  ) +
  theme(legend.position = "bottom")

#-----------------------------------------------------------------------------------------
# Area Under the Curve (AUC)
#-----------------------------------------------------------------------------------------
# Calculating the Area Under the Curve (AUC) is a brilliant way to quantify the "Halo." 
# Instead of just saying the halo is "bigger," you can tell Hussein: 
# "The total volume of paracrine-influenced space in A-CLL is X times larger than in C-CLL."
# In this context, the AUC represents the cumulative spatial influence. 
# A larger AUC for the Paracrine line indicates that the influence of the Proliferation Center 
# is not just strong, but geographically expansive.
#-----------------------------------------------------------------------------------------

# We will use the Trapezoidal Rule to integrate the area under the proportion lines across the distance bins.

library(dplyr)

# 1. Function to calculate AUC using the trapezoidal rule
calculate_auc <- function(x, y) {
  # Ensure data is sorted by distance
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  # Calculate area: sum of ((y_i + y_i+1)/2 * (x_i+1 - x_i))
  sum(diff(x) * (head(y, -1) + tail(y, -1))) / 2
}

# 2. Calculate AUC for each SL Tier
# We filter for the 'Outside' area (midpoint > 0) to measure the "Halo" expansion
auc_stats <- composition_sl %>%
  filter(midpoint >= 0) %>% # Focus on the signal extending into the tissue
  group_by(final_class) %>%
  summarise(
    Total_Spatial_Influence_AUC = calculate_auc(midpoint, proportion)
  )

print("--- Spatial Influence AUC (Outside PC) ---")
print(auc_stats)

# 3. Calculate the 'Halo Ratio' (Paracrine / Contact)
# This shows if the influence is spreading more through distance than direct contact
if(all(c("SL_contact", "SL_paracrine") %in% auc_stats$final_class)) {
  contact_auc <- auc_stats$Total_Spatial_Influence_AUC[auc_stats$final_class == "SL_contact"]
  paracrine_auc <- auc_stats$Total_Spatial_Influence_AUC[auc_stats$final_class == "SL_paracrine"]
  
  cat("\nHalo Expansion Ratio (Paracrine AUC / Contact AUC):", round(paracrine_auc / contact_auc, 2))
}




library(dplyr)
library(tidyr)

# 1. Define the Zones (Inside vs Outside)
niche_leak_summary <- cells_final %>%
  st_drop_geometry() %>%
  # Filter only for our SL sub-populations
  filter(final_class %in% c("SL_contact", "SL_paracrine", "SL_resting")) %>%
  mutate(zone = if_else(dist_signed < 0, "Inside_PC", "Outside_PC")) %>%
  
  # 2. Group by Zone and Sub-class
  group_by(zone, final_class) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  
  # 3. Calculate Percentages within each Zone
  mutate(Percentage = (n / sum(n)) * 100) %>%
  ungroup()

# 4. Format for a "Presentation-Ready" Table
niche_intensity_table <- niche_leak_summary %>%
  select(zone, final_class, Percentage) %>%
  pivot_wider(names_from = final_class, values_from = Percentage) %>%
  # Calculate Total Interaction Intensity (Contact + Paracrine)
  mutate(Total_Interaction_Intensity = SL_contact + SL_paracrine)

print("--- Niche Interaction Intensity: Inside vs. Outside PC ---")
print(niche_intensity_table)

# 5. Calculate the "Leak Coefficient"
# This is the ratio of Outside Intensity to Inside Intensity
outside_val <- niche_intensity_table$Total_Interaction_Intensity[niche_intensity_table$zone == "Outside_PC"]
inside_val <- niche_intensity_table$Total_Interaction_Intensity[niche_intensity_table$zone == "Inside_PC"]

cat("\nNiche Leak Coefficient:", round(outside_val / inside_val, 3))


# How to interpret the "Leak Coefficient" :
# 
# •	Low Leak (< 0.5): The activation is well-contained. 
# 
# Most SLs inside the PC are interacting, but they quickly become "Resting" once they are outside. 
# This is typical for stable C-CLL.
# 
# •	High Leak (> 0.8): The boundary is blurred. 
# 
# The SL cells outside the PC are being activated almost as intensely as those inside. 
# This is the "Activation Halo" Hussein is looking for in A-CLL.



#---------------------------------------------------------------------
# 
# The "Sub-Class" Mass Profiler Code
#
#---------------------------------------------------------------------
library(dplyr)
library(ggplot2)

# 1. Prepare binned data focusing ONLY on SL cells
sl_prepared <- cells_final %>%
  st_drop_geometry() %>%
  # Filter for our 3 specific SL subclasses
  filter(final_class %in% c("SL_contact", "SL_paracrine", "SL_resting")) %>% 
  mutate(dist_signed_um = dist_signed * 0.25) %>%
  filter(dist_signed_um > -150 & dist_signed_um < 150) %>%
  mutate(bin = cut(dist_signed_um, breaks = seq(-150, 150, by = 10))) %>%
  mutate(midpoint = (as.numeric(sub("\\((.+),(.+)]", "\\1", bin)) + 
                     as.numeric(sub("\\((.+),(.+)]", "\\2", bin))) / 2)

# 2. Identify features (Handling the exclusion list)
exclude_cols <- read.csv("ParametersExclude.csv")
all_features <- colnames(sl_prepared)

# Ensure exclude_cols is a character vector for setdiff
exclude_vec <- as.character(exclude_cols[[1]])

features_for_plot <- setdiff(all_features, 
                             c(exclude_vec, "geometry", "dist_signed", "dist_signed_um",
                               "class", "final_class", "bin", "midpoint", 
                               "ObjectNumber", "ImageNumber", "dist_to_nearest_active"))

# 3. Updated Mass Plotting Function
plot_sl_subclasses <- function(df, param_name, index) {
  
  profile_data <- df %>%
    group_by(midpoint, final_class) %>%
    summarise(
      mean_val = mean(.data[[param_name]], na.rm = TRUE),
      sd_val   = sd(.data[[param_name]], na.rm = TRUE),
      n = n(),
      se_val   = sd_val / sqrt(n),
      .groups = "drop"
    )
  
  # Standardize colors for the 3-Tier Strategy
  tier_colors <- c(
    "SL_contact" = "#E41A1C",   # Red
    "SL_paracrine" = "#4DAF4A",  # Green
    "SL_resting" = "#377EB8"     # Blue
  )
  
  p <- ggplot(profile_data, aes(x = midpoint, y = mean_val, color = final_class, group = final_class)) +
    geom_line(size = 1.2) +
    # Error ribbons help Hussein see if differences are statistically significant
    geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val, fill = final_class), 
                alpha = 0.15, color = NA) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    scale_color_manual(values = tier_colors) +
    scale_fill_manual(values = tier_colors) +
    theme_minimal() +
    labs(
      title = paste("Spatial Morphometric Gradient:", param_name),
      subtitle = "Red (Contact) | Green (Paracrine) | Blue (Resting)",
      x = "Distance to PC Border (µm) [Inside < 0 > Outside]",
      y = paste("Mean", param_name),
      color = "SL Tier",
      fill = "SL Tier"
    ) +
    theme(legend.position = "bottom")
  
  file_name <- paste0("sl_subclass_plots/", sprintf("%03d", index), "_", param_name, ".png")
  ggsave(file_name, plot = p, width = 10, height = 5, dpi = 150)
}

# 4. Create directory and run loop
if(!dir.exists("sl_subclass_plots")) dir.create("sl_subclass_plots")

for (i in seq_along(features_for_plot)) {
  param <- features_for_plot[i]
  # try() prevents one bad parameter from crashing the whole loop
  try(plot_sl_subclasses(sl_prepared, param, i))
}






