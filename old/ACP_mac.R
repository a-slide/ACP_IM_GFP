#! /usr/bin/Rscript

library(calibrate) # Pour la fonction textxy
library(FactoMineR)
x11()

# Inport data and transpose matrix
file <- commandArgs(TRUE)[1]
raw_data <- read.table(file, dec=".", sep="\t")
raw_data <- t(raw_data)

# Extract numeric data from raw data and modify the class to numeric
data <- raw_data[-1,-1]
class(data) <- "numeric"

# Extract columns and row names and affect them to the data matrix
genes <- raw_data[1,-1]
animals <- raw_data[-1,1]
rownames(data) <- animals
colnames(data) <- genes

# TODO = create a dataframe instead and implement a for loop

# Create sub group for separate analyses
J7 <- data[c(1,5,9,13,17,21),]
J30 <- data[c(2,6,10,14,18,22),]
J90 <- data[c(3,7,11,15,19),]
JF <- data[c(4,8,12,16,20,23),]

# Perform ACP with FactoMineR pakage
ACP_J7 <- PCA(J7, scale.unit=TRUE, graph = FALSE)
ACP_J30 <- PCA(J30, scale.unit=TRUE, graph = FALSE)
ACP_J90 <- PCA(J90, scale.unit=TRUE, graph = FALSE)
ACP_JF <- PCA(JF, scale.unit=TRUE, graph = FALSE)

# Create and export plots for all ACP results
plot(ACP_J7, axes = c(1, 2), choix = "ind", title = "ACP individus à J7")
dev.print(file="ACP ind J7.svg", device=svg)
plot(ACP_J7, axes = c(1, 2), choix = "var", title = "ACP variables à J7", cex = 0.5)
dev.print(file="ACP var J7.svg", device=svg)

plot(ACP_J30, axes = c(1, 2), choix = "ind", title = "ACP individus à J30")
dev.print(file="ACP ind J30.svg", device=svg)
plot(ACP_J30, axes = c(1, 2), choix = "var", title = "ACP variables à J30", cex = 0.5)
dev.print(file="ACP var J30.svg", device=svg)

plot(ACP_J90, axes = c(1, 2), choix = "ind", title = "ACP individus à J90")
dev.print(file="ACP ind J90.svg", device=svg)
plot(ACP_J90, axes = c(1, 2), choix = "var", title = "ACP variables à J90", cex = 0.5)
dev.print(file="ACP var J90.svg", device=svg)

plot(ACP_JF, axes = c(1, 2), choix = "ind", title = "ACP individus à JF")
dev.print(file="ACP ind JF.svg", device=svg)
plot(ACP_JF, axes = c(1, 2), choix = "var", title = "ACP variables à JF", cex = 0.5)
dev.print(file="ACP var JF.svg", device=svg)

# Create summary of highly corelated variables with dim1 and dim2 and export to CSV files

Dim_J7 <- dimdesc(ACP_J7, axes = c(1, 2, 3))

outputJ_7 <- rbind (
	c("Dim1.Correlation", "Dim1.p-Value"), Dim_J7$Dim.1$quanti,
	c("Dim2.Correlation", "Dim2.p-Value"), Dim_J7$Dim.2$quanti,
	c("Dim3.Correlation", "Dim3.p-Value"), Dim_J7$Dim.3$quanti)

write.table(outputJ_7, "ACP_J7.csv", sep="\t", col.names = FALSE)

Dim_J30 <- dimdesc(ACP_J30, axes = c(1, 2))

write.table(Dim_J30$Dim.1$quanti, "ACP_J30_Dim1.csv", sep="\t")
write.table(Dim_J30$Dim.2$quanti, "ACP_J30_Dim2.csv", sep="\t")

Dim_J90 <- dimdesc(ACP_J90, axes = c(1, 2))
write.table(Dim_J90$Dim.1$quanti, "ACP_J90_Dim1.csv", sep="\t")
write.table(Dim_J90$Dim.2$quanti, "ACP_J90_Dim2.csv", sep="\t")

Dim_JF <- dimdesc(ACP_JF, axes = c(1, 2))
write.table(Dim_JF$Dim.1$quanti, "ACP_JF_Dim1.csv", sep="\t")
write.table(Dim_JF$Dim.2$quanti, "ACP_JF_Dim2.csv", sep="\t")
