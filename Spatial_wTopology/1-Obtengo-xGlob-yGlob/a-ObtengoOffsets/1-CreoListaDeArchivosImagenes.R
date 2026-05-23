


# Load libraries
library(dplyr)
library(readr)
library(stringr)

# 1. Define the folder path (use "." for current folder)
folder_path <- getwd()

# 2. Get list of .tif files
file_names <- list.files(path = folder_path, pattern = "\\.TIF$", full.names = FALSE)



# 3. Create a clean data frame with metadata extraction
file_inventory <- data.frame(full_filename = file_names) %>%
  mutate(
    # Extracts the number following 'f' (e.g., 0 from f00 or 259 from f259)
    field_index = as.numeric(str_extract(full_filename, "(?<=f)\\d+")),
    # Captures the file extension
    extension = "tif",
    # Provides the full relative path for future R processing
    file_path = paste0(folder_path, "/", full_filename)
  ) %>%
  # Sort by field_index so the CSV follows the scanning order
  arrange(field_index)

# 4. Save to CSV
write_csv(file_inventory, "1-scan_file_list.csv")

# Print summary to console
cat("Success: xx files processed and saved to scan_file_list.csv\n")
