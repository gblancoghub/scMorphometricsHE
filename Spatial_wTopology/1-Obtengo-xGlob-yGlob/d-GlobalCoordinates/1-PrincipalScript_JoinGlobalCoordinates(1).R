

library(dplyr)



# 1. Load all three tables
cells     <- read.csv("my_table_WithClasses.csv")  # no tiene filename
bridge    <- read.csv("CproINum_To_FileName.csv") # tiene Image_FileName_HE d4.TIF (lo paso a d0.TIF
offsets   <- read.csv("2-evos_global_coordinates.csv") # tiene full_filename d0.TIF

dim(cells)
dim(bridge)
dim(offsets)

# 2. Step One: Link Cells to Filenames
# This adds the "full_filename" column to your 16k cell records
cells_with_names <- cells %>%
  left_join(bridge, by = "ImageNumber")

colnames(cells_with_names)
cells_with_names[,"full_filename"]

colnames(offsets)



# 3. Step Two: Link Filenames to Spatial Offsets
# Note: We use 'full_filename' or 'Image_FileName_HE' as the key.
# Adjust the 'by' argument if the column name in offsets is different.
final_table <- cells_with_names %>%
  left_join(offsets %>% select(full_filename, offX, offY), 
            by = join_by(full_filename)) # Este es el que anda


head(final_table)
dim(final_table)

colnames(final_table)
final_table[19000:19005,1:4]
final_table[19000:19005,516:524]

# 4. Step Three: Calculate Global Coordinates
final_table <- final_table %>%
  mutate(
    xGlob = offX + NucObjects_AreaShape_Center_X,
    yGlob = offY + NucObjects_AreaShape_Center_Y
  )

dim(final_table)
colnames(final_table)
final_table[19000:19005,516:526]

final_table <- final_table %>%
  mutate(
    xGlobUn = xGlob,
    yGlobUn = yGlob
  )

dim(final_table)
colnames(final_table)
final_table[19000:19005,516:528]

# Las GlobUn son un backup de las coordenadas globales
# para controlar y comparar con las geometrias cuando se crea el spatial object
# Porque al crear la geometria xGlob e yGlob van a desaparecer

# 5. Export the Final "Intelligent Map"
write.csv(final_table, "CLL_Global_Spatial_Data.csv", row.names = FALSE)

# Check for any NAs (which would indicate a mismatch in filenames)
if(any(is.na(final_table$xGlob))) {
  warning("Some cells did not receive global coordinates. Check filename consistency.")
} else {
  cat("Success! All cells registered to global coordinate system.")
}

