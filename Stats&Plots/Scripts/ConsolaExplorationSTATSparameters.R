
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)




seurat.integrated <- readRDS("seurat_integrated_wCellType.rds")


view(seurat.integrated@meta.data)
colnames(seurat.integrated@meta.data)
rownames(seurat.integrated@meta.data)

PIB <- read.csv("my_tableClass_All.csv")
seurat.integrated@meta.data <- cbind(seurat.integrated@meta.data, PIB)


DimPlot(seurat.integrated, reduction = 'umap', group.by = "class", 
label = FALSE, shuffle = TRUE, order ="PIB",pt.size=NULL,  alpha=0.3)

str(seurat.integrated@assays$RNA@layers$counts.1)




Sparemat <- seurat.integrated@assays$RNA@layers$counts.1
mat <- as.matrix(Sparemat)
mat <- t(mat)
mat <- as.data.frame(mat)
mat00 <- mat


Sparemat <- seurat.integrated@assays$RNA@layers$counts.2
mat <- as.matrix(Sparemat)
mat <- t(mat)
mat <- as.data.frame(mat)
mat01 <- mat

Sparemat <- seurat.integrated@assays$RNA@layers$counts.3
mat <- as.matrix(Sparemat)
mat <- t(mat)
mat <- as.data.frame(mat)
mat02 <- mat

Sparemat <- seurat.integrated@assays$RNA@layers$counts.4
mat <- as.matrix(Sparemat)
mat <- t(mat)
mat <- as.data.frame(mat)
mat04 <- mat

mat <- rbind(mat00,mat01,mat02,mat04)

dim(mat)


rownames(mat)<-  colnames(seurat.integrated@assays$RNA)
colnames(mat)<- rownames(seurat.integrated@assays$RNA)

rownames(mat)
colnames(mat)

mat[, c(1,2)]

colnames(seurat.integrated@meta.data)
head()
Mdata <- as.data.frame(seurat.integrated@meta.data[,c("Patient", "CellType", "BroadCellType","class")])
colnames(Mdata)

dat1 <- cbind(mat,Mdata)
colnames(dat1)
dim(dat1)
table(dat1$Patient, dat1$class)
write.csv(dat1,"tablaMainALllPatients77K.csv")
# dat1 <- read.csv("tablaMainALllPatients77K.csv")

 
dat1[,"CellObjects-Mean-ASCspecs-AreaShape-EulerNumber"]
max(dat1[,"CellObjects-AreaShape-Eccentricity"])

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
colnames(dat1)

df <- as.data.frame(dat1[,c("Patient", "CellType", "BroadCellType","class","NucObjects-AreaShape-EquivalentDiameter")])
colnames(df)[5] <- "size"

df <- as.data.frame(dat1[,c("Patient", "CellType", "BroadCellType","class","CellObjects-AreaShape-EquivalentDiameter")])
colnames(df)[5] <- "size"

df <- as.data.frame(dat1[,c("Patient", "CellType", "BroadCellType","class","CellObjects-Mean-ASCspecs-AreaShape-EulerNumber"  )])
colnames(df)[5] <- "size"

  
df <- as.data.frame(dat1[,c("Patient", "CellType", "BroadCellType","class","CellObjects-Children-ASCspecs-Count"  )])
colnames(df)[5] <- "size"

###------------
df <- as.data.frame(dat1[,c("Patient", "CellType", "BroadCellType","class","CellObjects-Mean-ASCspecs-Granularity-11-HEGray")])
colnames(df)[5] <- "size"


df <- as.data.frame(dat1[,c("Patient", "CellType", "BroadCellType","class","NucObjects-Granularity-2-HEGray")])
colnames(df)[5] <- "size"


df <- as.data.frame(dat1[,c("Patient", "CellType", "BroadCellType","class","Cytoplasm-Granularity-9-HEGray")])
colnames(df)[5] <- "size"


###------------
df <- as.data.frame(dat1[,c("Patient", "CellType", "BroadCellType","class","CellObjects-Mean-ASCspecs-Texture-Variance-HEGray-3-02-256")])
colnames(df)[5] <- "size"


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
First, let's get the summary statistics you need. dplyr is the best tool for this.
R
# Load necessary libraries
library(dplyr)
library(ggplot2) # For plotting



# 'class' and 'Patient' are factors for proper grouping and plotting
df$class <- factor(df$class, levels = c("SL", "PL", "PIB")) # Set order for plotting
df$Patient <- factor(df$Patient, levels = c("00_PAC01", "04_PAC04", "01_PAC02", "02_PAC03"))
# --- End of Sample Data ---


# Calculate mean and SD of 'size' by 'class' and 'Patient'
summary_df <- df %>%
  group_by(BroadCellType, Patient) %>%
  summarise(
    Mean_Size = mean(size, na.rm = TRUE),
    SD_Size = sd(size, na.rm = TRUE),
    N_Cells = n(), # Also good to know the number of cells contributing to each mean
    SE_Size = SD_Size / sqrt(N_Cells) # Standard Error of the Mean, often preferred for error bars
  ) %>%
  ungroup() # ungroup after summarise for subsequent operations

# Print the summary table
print(summary_df)


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

R
# Box plot: 'size' distribution by 'class', split by 'Patient'
ggplot(df, aes(x = BroadCellType, y = size, fill = Patient)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.alpha = 0.3) + # position_dodge separates boxes
  labs(
    title = "Distribution of Size by Cell Class and Patient",
    x = "Cell Class",
    y = "Size (Arbitrary Units)",
    fill = "Patient"
  ) +
  scale_fill_brewer(palette = "Set2") + # Use a colorblind-friendly palette
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels if needed
    legend.title = element_text(face = "bold")
  )+
  coord_cartesian(ylim = c(0, 5))

# Alternative Box Plot: Faceting by Patient (good if many patients or classes)
ggplot(df, aes(x = class, y = size, fill = class)) + # Fill by class for within-patient comparison
  geom_boxplot(outlier.alpha = 0.3) +
  facet_wrap(~ Patient, scales = "free_y", ncol = 2) + # Separate plots for each patient
  labs(
    title = "Distribution of Size by Cell Class (Faceted by Patient)",
    x = "Cell Class",
    y = "Size (Arbitrary Units)",
    fill = "Cell Class"
  ) +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"), # Facet titles
    legend.position = "none" # Hide redundant legend if filling by facet variable
  )

///////////
# Violin plot: 'size' distribution by 'class', split by 'Patient'
ggplot(df, aes(x = class, y = size, fill = Patient)) +
  geom_violin(position = position_dodge(width = 0.8), trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), outlier.alpha = 0.3) + # Overlay boxplot for median/quartiles
  labs(
    title = "Distribution of Size by Cell Class and Patient (Violin Plot)",
    x = "Cell Class",
    y = "Size (Arbitrary Units)",
    fill = "Patient"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(face = "bold")
  )

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


custom_class_colors <- c(
  "PIB" = "red",
  "PL" = "yellow",
  "SL" = "blue"
)

custom_class_sizes <-  c(
  "PIB" = 1.4,
  "PL" = 1,
  "SL" = 0.6
)



seurat.integrated@meta.data$class == "PIB"
var <- seurat.integrated@meta.data$class
varPIB <- var == "PIB"
varPL <- var == "PL"
varSL <- var == "SL"

custom_class_sizes <- c(
  varPIB = 1.8,
  varPL = 1.0,
  varSL = 0.6
)

seurat.integrated@meta.data$class[var]


 
#---------- 1
DimPlot(
  object = seurat.integrated,
  reduction = 'umap',
  group.by = "class", # Map 'class' to color
  label = FALSE,
  shuffle = TRUE,
  order = "PIB",
  pt.size= NULL 
  ) +
  # Apply the custom color scale
  scale_color_manual(values = custom_class_colors) +
  # Apply the custom size scale
  scale_size_manual(values = custom_class_sizes) +
  # Optional: Customize legend and titles for clarity
  labs(
    title = "Title",
    color = "Cell Class", # Legend title for color
    size = "Point Size" # Legend title for size
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10)
  )


#---------- 2

DimPlot(seurat.integrated, reduction = 'umap', group.by = "class", split.by = "Patient", repel=TRUE,
label = FALSE, shuffle = TRUE, order ="PIB", alpha=0.2,pt.size= y )+
  scale_color_manual(values = custom_class_colors) +
  labs(
    title = "Title",
    color = "Cell Class", # Legend title for color
    size = "Point Size" # Legend title for size
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10)
  )


#---------- 3

x <- as.numeric(as.factor(seurat.integrated@meta.data$class))
# y <- x/3


DimPlot(seurat.integrated, reduction = 'umap', group.by = "class", split.by = "Patient", repel=TRUE,
label = FALSE, shuffle = TRUE, order ="PIB", alpha=0.2, pt.size= y )+
 # scale_size_manual(values =  c("PIB" = "1.4","PL" = "1","SL" = "0.6"), breaks = c("PIB", "PL", "SL"),limits = c("PIB", "PL", "SL")) +
  labs(
    title = "Title",
    color = "Cell Class", # Legend title for color
    size = "Point Size" # Legend title for size
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10)
  )


DimPlot(seurat.integrated, reduction = 'umap', group.by = "class", split.by = "Patient", repel=TRUE,
label = TRUE, shuffle = TRUE, order ="PIB", alpha=0.2, cells.highlight=list(PIBlist, PL_list),
 pt.size=NULL, cols.highlight = c("darkblue", "darkred"), ncol=2)


DimPlot(seurat.integrated, reduction = 'umap', group.by = "Patient", split.by = "class", label = TRUE)
Idents(seurat.integrated)<- "BroadCellType"
DimPlot(seurat.integrated, reduction = 'umap', split.by = "class", label = TRUE, order =TRUE)

DimPlot(seurat.integrated, reduction = 'umap', group.by = "class", split.by = "BroadCellType", 
repel=TRUE,
label = TRUE, shuffle = TRUE, order ="PIB", alpha=0.2, cells.highlight=list(PIBlist, PL_list),
 pt.size=NULL, cols.highlight = c("darkblue", "darkred", "green"), ncol=2)



Idents(seurat.integrated)<- "class"
PIBlist <- WhichCells(seurat.integrated, idents = c("PIB"))
PL_list <- WhichCells(seurat.integrated, idents = c( "PL"))
SL_list <- WhichCells(seurat.integrated, idents = c( "SL"))


table(seurat.integrated@meta.data$CellType,seurat.integrated@meta.data$class)
table(seurat.integrated@meta.data$class, seurat.integrated@meta.data$CellType)


DimPlot(seurat.integrated, reduction = 'umap', group.by = "BroadCellType", label = TRUE)
DimPlot(seurat.integrated, reduction = 'umap', group.by = "BroadCellType", split.by = "Patient", label = TRUE)

DimPlot(seurat.integrated, reduction = 'tsne', group.by = "BroadCellType", split.by = "Patient", label = TRUE)
DimPlot(seurat.integrated, reduction = 'tsne', group.by = "BroadCellType", label = TRUE)

table(seurat.integrated@meta.data$BroadCellType,seurat.integrated@meta.data$class )
table(seurat.integrated@meta.data$class, seurat.integrated@meta.data$BroadCellType)


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

table(seurat.integrated@meta.data$BroadCellType,seurat.integrated@meta.data$class )
# table(seurat.integrated@meta.data$class, seurat.integrated@meta.data$BroadCellType)


cell_type_class_table <- table(seurat.integrated@meta.data$BroadCellType,seurat.integrated@meta.data$class )


# 1. Calculate Proportions (for better understanding)
# Row proportions (proportion of each cell type within a class)
row_props <- prop.table(cell_type_class_table, margin = 1)
print("Row Proportions (Cell Type Distribution within classes):")
print(round(row_props, 3))

# Column proportions (proportion of each class within a cell type)
col_props <- prop.table(cell_type_class_table, margin = 2)
print("Column Proportions (class Distribution within Cell Types):")
print(round(col_props, 3))

# Overall proportions
total_props <- prop.table(cell_type_class_table)
print("Overall Proportions:")
print(round(total_props, 3))


# 2. Chi-squared Test
chi_sq_result <- chisq.test(cell_type_class_table)
print("\nChi-squared Test Result:")
print(chi_sq_result)

# 3. Standardized Residuals
# how much each cell contributes to the chi-squared statistic.
# Values > 2 or < -2 (approx) are often considered significant deviations from expectation.
# Positive residuals indicate enrichment (more cells than expected)
# Negative residuals indicate depletion (fewer cells than expected)
standardized_residuals <- chi_sq_result$stdres
print("\nStandardized Residuals:")
print(round(standardized_residuals, 2))


# 4. Visualize Standardized Residuals (Heatmap)

# Convert residuals matrix to a long format for ggplot2
x <- as.data.frame(standardized_residuals)
colnames(x) <- c("BroadCellType","class","Residual")
residuals_df <- x

# Create the heatmap
ggplot(residuals_df, aes(x = class, y = BroadCellType, fill = Residual)) +
  geom_tile(color = "white", size = 0.5) + # Add white borders for clarity
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, # do not use this paraneter - limit = max(abs(residuals_df$Residual)), 
                       space = "Lab", name = "Std. Residual") +
  geom_text(aes(label = round(Residual, 1)), color = "black", size = 3) + # Add residual values as text
  labs(
    title = "Standardized Residuals of CellType by Class",
    subtitle = paste0("Chi-squared p-value: ", format.pval(chi_sq_result$p.value, digits = 3)),
    x = "class",
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
[1] "\nCramer's V (Effect Size): 0.159"
> 
>

# Column proportions (proportion of each class within a cell type)
col_props <- prop.table(cell_type_class_table, margin = 2)
print("Column Proportions (Class Distribution within Cell Types):")
print(round(col_props, 3))

x1 <- as.data.frame(col_props)
colnames(x1) <- c("BroadCellType","Class","PropOfCell")
Props <- x1

# Create the heatmap FOR PROPORTIONS OF CELLS
ggplot(Props , aes(x = Class, y = BroadCellType, fill = PropOfCell)) +
  geom_tile(color = "white", size = 0.5) + # Add white borders for clarity
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       # midpoint = 0, # do not use this paraneter - limit = max(abs(residuals_df$Residual)), 
                       space = "Lab", name = "Std. Residual") +
  geom_text(aes(label = round(PropOfCell*100, 1)), color = "black", size = 3) + # Add residual values as text
  labs(
    title = "Proportion of Cells - CellType by Class",
    subtitle = paste0("Chi-squared p-value: ", format.pval(chi_sq_result$p.value, digits = 3)),
    x = "Class",
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
 


