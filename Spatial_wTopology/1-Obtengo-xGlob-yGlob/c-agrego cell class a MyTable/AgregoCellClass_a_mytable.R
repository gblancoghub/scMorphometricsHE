




# input1:  "my_table.csv" - la tabla estandar de CPro results


# 
# Usar directo


# output: "my_table_WithClasses.csv"


library(dplyr)

# 1. Load your table from CPro results without classes (MyTable.csv)

#Table has more than 130.000 cells so I save an RDS object before

# TableWOclasses <- read.csv("my_table_00Full-Rosetto.csv")
colnames(TableWOclasses)
dim(TableWOclasses)

# saveRDS(TableWOclasses, "TableWOclasses.rds")

TableWOclasses <- readRDS("TableWOclasses.rds")

# 2. Load your class reference CSV
classesRef <- read.csv("my_table_Classes00Full_Rosetto.csv")
colnames(classesRef)
dim(classesRef)

# 3. Use base R merge (very robust for sf objects)
# We merge by both ID columns to ensure a perfect match
table_with_classes <- merge(
  TableWOclasses, 
  classesRef[, c("ImageNumber", "ObjectNumber", "class", "class_number")], 
  by = c("ImageNumber", "ObjectNumber"), 
  all.x = TRUE
)

colnames(table_with_classes)
dim(table_with_classes)

write.csv(table_with_classes, "my_table_WithClasses.csv")

x <- read.csv("my_table_WithClasses.csv")
dim(x)

control <- table_with_classes[,c(1,2,517,518,519)]
write.csv(control, "controlClasses.csv")

saveRDS(table_with_classes, "table_with_classes.rds")
