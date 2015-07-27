#! /usr/bin/Rscript

### VERIFY THE NUMBER OF ARG ###

args <- commandArgs(trailingOnly = TRUE)

if (length (args) != 2)
    stop("2 arguments are required:
    Filename of csv file containing the data
    Max percentage of values with NA allows per column
    Example : ACP_IM_GFP.R ./data/Sample.csv 50")

### LOAD PACKAGES ###

options(width = 160) # wider terminal size
library(FactoMineR)
library(RColorBrewer)
library(gplots)
x11()

### IMPORT DATA ###

# Inport data in a dataframe
filename = args[1]
sample_data = read.table(filename, dec=".", sep="\t", header = TRUE, row.names=1)

# Extract basename of the filename
sample_name = sub("^([^.]*).*", "\\1", basename(filename))
print (paste("ANALYSING SAMPLE  ", sample_name))

# Calculate the maximal number of NA allowed
NA_max_percent = as.numeric(args[2])
NA_max = round ((NA_max_percent*nrow(sample_data)/100), digits=0)
print (paste("MAX NA ALLOWED BY GROUP = " , NA_max))

### PREPROCESS DATA ###

# Open a file for writing a report of eliminated and retained mir
report_name = paste("Report_", sample_name,".csv", collapse="")

# List genes with less NA value than max NA and export to file
write ("ELIMINATED GENES", report_name)
slice = paste((names(sample_data[,colSums(is.na(sample_data)) > NA_max])), collapse = "\t")
write (slice, report_name, append=T)

# Remove column with more NA that the maximum allowed
sample_data = sample_data[,colSums(is.na(sample_data)) <= NA_max]

# List retained genes with NA
write ("RETAINED GENES WITH NA", report_name, append=T)
slice = paste((names(sample_data[,colSums(is.na(sample_data)) > 0] )), collapse = "\t")
write (slice, report_name, append=T)

# Retained Mir for analysis
write ("RETAINED GENES WITHOUT NA", report_name, append=T)
slice = paste((names(sample_data[,colSums(is.na(sample_data)) == 0])), collapse = "\t")
write (slice, report_name, append=T)

# replace NA by 0
sample_data[is.na(sample_data)] = 0

# Scale and center data
sample_data = scale(sample_data)
#print (sample_data)

### PCA ###

# Perform PCA with FactoMineR pakage
PCA.res = PCA (sample_data, graph = FALSE)

# Create and export plots for all ACP results
cos2_PC1_PC2 = rowSums (PCA.res$ind$cos2[,c("Dim.1","Dim.2") ], dims = 1)
plot(PCA.res, axes = c(1, 2), choix = "ind", title = sample_name, cex = cos2_PC1_PC2/2+0.6)
dev.print(file = paste(sample_name, "_PCA_samples_PC1_PC2.svg", collapse=""), device=svg)

cos2_PC1_PC3 = rowSums (PCA.res$ind$cos2[,c("Dim.1","Dim.3") ], dims = 1)
plot(PCA.res, axes = c(1, 3), choix = "ind", title = sample_name, cex = cos2_PC1_PC3/2+0.6)
dev.print(file = paste(sample_name, "_PCA_samples_PC1_PC3.svg", collapse=""), device=svg)

plot(PCA.res, axes = c(1, 2), choix = "var", title = sample_name, cex = 0.5, lim.cos2.var = 0.6)
dev.print(file = paste(sample_name, "_PCA_variables_PC1_PC2.svg", collapse=""), device=svg)

plot(PCA.res, axes = c(1, 3), choix = "var", title = sample_name, cex = 0.5, lim.cos2.var = 0.6)
dev.print(file = paste(sample_name, "_PCA_variables_PC1_PC3.svg", collapse=""), device=svg)

# Plot the variance explained
eigen_plot = barplot(PCA.res$eig[,2], names=paste("Dim",1:nrow(PCA.res$eig)), main = sample_name)
dev.print(file = paste(sample_name, "_PCA_variance.svg", collapse=""), device=svg)

# Create summary of highly correlated variables with dim1, dim2, dim3 and export to CSV files
Dim = dimdesc(PCA.res, axes = c(1, 2, 3))
output = rbind (
    c("Dim1.Correlation", "Dim1.p-Value"), Dim$Dim.1$quanti,
    c("Dim2.Correlation", "Dim2.p-Value"), Dim$Dim.2$quanti,
    c("Dim3.Correlation", "Dim3.p-Value"), Dim$Dim.3$quanti)
write.table(output, paste(sample_name, "_PCA_var_dim.csv", collapse=""), sep="\t", col.names = FALSE)

### HEATMAP ###

## Row clustering
dist_mat = as.dist(1-cor(sample_data, method="pearson"))
hr <- hclust(dist_mat, method="complete")
sample_data = t(sample_data)

heatmap.2(
    sample_data,
    Rowv=as.dendrogram(hr),
    Colv= FALSE,
    dendrogram = "row",
    scale="row",
    density.info="density",
    trace = "column",
    tracecol = "black",
    col = colorRampPalette(c("springgreen4","palegreen1","white","coral1","red4")),
    main = sample_name,
    cexCol = 0.7,
    cexRow = 0.5)

dev.print(file = paste(sample_name, "_Heatmap.svg", collapse=""), device=svg)

warnings()
