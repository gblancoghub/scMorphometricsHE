



library(dplyr)
library(readr)

# 1. Load your inventory CSV
df <- read_csv("1-scan_file_list.csv")

# HACER UPDATE DE K=14
# 2. Define grid parameters
cols_n <- 14
width_px <- 2048
height_px <- 1536

# 3. Calculate Serpentine Offsets
df_spatial <- df %>%
  mutate(
    # Row index (0 to 19)
    i_offY = field_index %/% cols_n,
    
    # Column index with serpentine logic
    i_offX = ifelse(i_offY %% 2 == 0, 
                    field_index %% cols_n,               # Even rows: left to right
                    (cols_n - 1) - (field_index %% cols_n) # Odd rows: right to left
             ),
    
    # Calculate pixel offsets
    offX = i_offX * width_px,
    offY = i_offY * height_px
  )

# 4. Save the "Intelligent" Coordinate Table
write_csv(df_spatial, "2-evos_global_coordinates.csv")

# Diagnostic: View the first 2 rows of the first 3 "lines" of the scan
print(df_spatial %>% select(field_index, i_offX, i_offY, offX, offY) %>% slice(c(1:2, 13:14, 26:27)))
