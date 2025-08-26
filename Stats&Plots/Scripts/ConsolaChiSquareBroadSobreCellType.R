
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)




# 5. Integrate Data

seurat.integrated <- readRDS("seurat_integrated.rds")

# Load the dplyr library if you haven't already
library(dplyr)
seurat.integrated@meta.data <- seurat.integrated@meta.data %>%
  mutate(
    CellType = case_when(
      seurat_clusters %in% c("0", "1") ~ "M1",
      seurat_clusters %in% c("40", "43", "45", "46", "19", "37") ~ "M2",
      seurat_clusters %in% c("14", "47") ~ "M3",
      seurat_clusters %in% c("13", "11", "18", "28", "5") ~ "R0",
      seurat_clusters %in% c("35", "24") ~ "D1",
      seurat_clusters %in% c("6", "42", "30") ~ "D2",
      seurat_clusters %in% c("34", "39", "12", "3") ~ "D3",
      seurat_clusters %in% c("20", "25", "44") ~ "ADB1",
      seurat_clusters %in% c("7", "10", "23") ~ "ADB2",
      seurat_clusters %in% c("41", "8") ~ "ADBR1",
      seurat_clusters %in% c("38", "33", "9", "32") ~ "ADMR1",
      seurat_clusters %in% c("16", "17", "27", "21", "22", "31","36", "26") ~ "ADTR1",
      seurat_clusters %in% c("4", "29", "15", "2") ~ "ADT1",
      TRUE ~ NA_character_ # Fallback for any unassigned clusters (shouldn't happen with full coverage)
    )
  )

# Optional: Convert CellType to a factor with a specific order if desired for plotting
# (e.g., to control the order in legends or facets)
seurat.integrated@meta.data$CellType <- factor(
  seurat.integrated@meta.data$CellType,
  levels = c("M1", "M2", "M3", "R0", "D1", "D2", "D3", "ADB1", "ADB2", "ADBR1", "ADMR1", "ADTR1", "ADT1")
)

table(seurat.integrated@meta.data$CellType)

# You can also visualize it with DimPlot

DimPlot(seurat.integrated, reduction = 'umap', group.by = "CellType", label = TRUE)
DimPlot(seurat.integrated, reduction = 'umap', group.by = "CellType", split.by = "Patient", label = TRUE)
DimPlot(seurat.integrated, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)

table(seurat.integrated@meta.data$CellType,seurat.integrated@meta.data$Patient )
table(seurat.integrated@meta.data$Patient, seurat.integrated@meta.data$CellType)


seurat.integrated@meta.data <- seurat.integrated@meta.data %>%
  mutate(
    BroadCellType = case_when(
      CellType %in% c("M1", "M2", "M3") ~ "Mature",
      CellType %in% c("D1", "D2", "D3") ~ "Differentiation",
      CellType %in% c("ADB1", "ADB2") ~ "AlternativeDiff_B",
      CellType %in% c("ADBR1", "ADMR1") ~ "AlternativeDiff_MR",
       CellType %in% c("ADTR1","ADT1") ~ "AlternativeDiff_TR",
     CellType %in% c("R0") ~ "RootStem",
      TRUE ~ NA_character_ # Fallback for any unassigned clusters (shouldn't happen with full coverage)
    )
  )


table(seurat.integrated@meta.data$BroadCellType)
# Optional: Convert CellType to a factor with a specific order if desired for plotting
# (e.g., to control the order in legends or facets)
seurat.integrated@meta.data$BroadCellType <- factor(
  seurat.integrated@meta.data$BroadCellType,
  levels = c("RootStem", "Differentiation", "Mature", "AlternativeDiff_B", "AlternativeDiff_MR", "AlternativeDiff_TR")
)



DimPlot(seurat.integrated, reduction = 'umap', group.by = "BroadCellType", label = TRUE)
DimPlot(seurat.integrated, reduction = 'umap', group.by = "BroadCellType", split.by = "Patient", label = TRUE)
DimPlot(seurat.integrated, reduction = 'tsne', group.by = "BroadCellType", split.by = "Patient", label = TRUE)
DimPlot(seurat.integrated, reduction = 'tsne', group.by = "BroadCellType", label = TRUE)


table(seurat.integrated@meta.data$BroadCellType,seurat.integrated@meta.data$Patient )
table(seurat.integrated@meta.data$Patient, seurat.integrated@meta.data$BroadCellType)

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
# 
# STATS
# 
###########################################################################################
 
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(reshape2) # For melting data for ggplot2
library(ggpubr) # For adding p-value to plot (optional)
library(corrplot) # For visualizing residuals (alternative to ggplot heatmap)


seurat.integrated <- readRDS("seurat_integrated_wCellType.rds")


cell_type_patient_table <- table(seurat.integrated@meta.data$BroadCellType, seurat.integrated@meta.data$Patient)

# 1. Proportions 
# Row proportions (proportion of each cell type within a patient)
row_props <- prop.table(cell_type_patient_table, margin = 1)
print("Row Proportions (Cell Type Distribution within Patients):")
print(round(row_props, 3))

# Column proportions (proportion of each patient within a cell type)
col_props <- prop.table(cell_type_patient_table, margin = 2)
print("Column Proportions (Patient Distribution within Cell Types):")
print(round(col_props, 3))

# Overall proportions
total_props <- prop.table(cell_type_patient_table)
print("Overall Proportions:")
print(round(total_props, 3))


# 2. Chi-squared Test
chi_sq_result <- chisq.test(cell_type_patient_table)
print("\nChi-squared Test Result:")
print(chi_sq_result)

# 3. Standardized Residuals

standardized_residuals <- chi_sq_result$stdres
print("\nStandardized Residuals:")
print(round(standardized_residuals, 2))


# 4. Visualize Standardized Residuals (Heatmap)

# Convert residuals matrix to a long format for ggplot2
x <- as.data.frame(standardized_residuals)
colnames(x) <- c("BroadCellType","Patient","Residual")
residuals_df <- x

# Create the heatmap
ggplot(residuals_df, aes(x = Patient, y = BroadCellType, fill = Residual)) +
  geom_tile(color = "white", size = 0.5) + # Add white borders for clarity
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, # do not use this paraneter - limit = max(abs(residuals_df$Residual)), 
                       space = "Lab", name = "Std. Residual") +
  geom_text(aes(label = round(Residual, 1)), color = "black", size = 3) + # Add residual values as text
  labs(
    title = "Standardized Residuals of CellType by Patient",
    subtitle = paste0("Chi-squared p-value: ", format.pval(chi_sq_result$p.value, digits = 3)),
    x = "Patient",
    y = "Broad Cell Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

# 5. Calculate Effect Size (Cramer's V)
# Cramer's V ranges from 0 to 1, where 0 means no association and 1 means perfect association.

n <- sum(cell_type_patient_table)
min_dim <- min(nrow(cell_type_patient_table), ncol(cell_type_patient_table))
cramers_v <- sqrt(chi_sq_result$statistic / (n * (min_dim - 1)))
print(paste0("\nCramer's V (Effect Size): ", round(cramers_v, 3)))


# 	Cramer's V quantifies the strength of the association regardless of sample size.
# 	Interpretation:
# 	0 to 0.1: Very weak association
# 	0.1 to 0.3: Weak to moderate association
# 	0.3 to 0.5: Moderate to strong association
# 	> 0.5: Strong association

> print(paste0("\nCramer's V (Effect Size): ", round(cramers_v, 3)))
[1] "\nCramer's V (Effect Size): 0.33"
>


# Column proportions (proportion of each patient within a cell type)
col_props <- prop.table(cell_type_patient_table, margin = 2)
print("Column Proportions (Patient Distribution within Cell Types):")
print(round(col_props, 3))

x1 <- as.data.frame(col_props)
colnames(x1) <- c("BroadCellType","Patient","PropOfCell")
Props <- x1

# Create the heatmap FOR PROPORTIONS OF CELLS
ggplot(Props , aes(x = Patient, y = BroadCellType, fill = PropOfCell)) +
  geom_tile(color = "white", size = 0.5) + # Add white borders for clarity
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       # midpoint = 0, # do not use this paraneter - limit = max(abs(residuals_df$Residual)), 
                       space = "Lab", name = "Std. Residual") +
  geom_text(aes(label = round(PropOfCell*100, 1)), color = "black", size = 3) + # Add residual values as text
  labs(
    title = "Proportion of Cells - CellType by Patient",
    subtitle = paste0("Chi-squared p-value: ", format.pval(chi_sq_result$p.value, digits = 3)),
    x = "Patient",
    y = "Broad Cell Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )


###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
# 
# END STATS
# 
###########################################################################################
 

