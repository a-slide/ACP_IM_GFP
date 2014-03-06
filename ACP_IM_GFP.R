#! /usr/bin/Rscript

options(width = 160) # wider terminal size
library(FactoMineR)
library(RColorBrewer)
library(gplots)
x11()

# Inport data in a dataframe
file = commandArgs(TRUE)[1]
data = read.table(file, dec=".", sep="\t", header = TRUE, row.names=1)

# Find and remove variables whith NA values
NA.mat = is.na(data)                        # Create a logical matrix T/F indicating if a value is NA
NA.col= apply(NA.mat, MARGIN=2, FUN = any)  # Check if any value in the column is NA
NA.index = which(NA.col)                    # Return the position of columns whith at least one NA
data = data[,-NA.index]                     # Remove NA containing columns from data

print (paste(length(NA.index), "genes had been removed due to NA values"))
write(names(NA.col[NA.col==TRUE]), "Eliminated_genes.txt")

# Create a list of row number to analyze data separatly
groups = c("J7", "J30", "J90", "Euthanasia", "Wilson", "Sheperd", "Monk", "Walcott","Ruben","Woody")

l = list(   c(1,9,17,25,33,41),     # Line range for J7
            c(2,10,18,26,34,42),    # Line range for J30
            c(3,11,19,27,35),       # Line range for J90
            c(4,12,20,28,36,43),    # Line range for eutha
            1:8,        # Line range for Wil
            9:16,       # Line range for Shep
            17:24,      # Line range for Monk
            25:32,      # Line range for Walc
            33:40,      # Line range for Rub
            41:46)      # Line range for Woo

names(l) = groups

#Loop foreach element in l

for (i in seq_along(l)){

    #Extract informations for the current group
    group.name = groups[i]
    group.index = l[[i]]
    group.data = data[group.index,]

    print (paste("ANALYSING GROUP ", group.name))
    print (row.names(group.data))

    # Perform ACP with FactoMineR pakage
    PCA.res = PCA(group.data, scale.unit=TRUE, graph = FALSE)

    # Create and export plots for all ACP results
    plot(PCA.res, axes = c(1, 2), choix = "ind", title = group.name)
    dev.print(file = paste("ACP_individus_", group.name,".svg", collapse=""), device=svg)

    plot(PCA.res, axes = c(1, 2), choix = "var", title = group.name, cex = 0.5)
    dev.print(file = paste("ACP_variables_", group.name,".svg", collapse=""), device=svg)

    # Create summary of highly corelated variables with dim1, dim2, dim3 and export to CSV files

    Dim = dimdesc(PCA.res, axes = c(1, 2, 3))

    output = rbind (
        c("Dim1.Correlation", "Dim1.p-Value"), Dim$Dim.1$quanti,
        c("Dim2.Correlation", "Dim2.p-Value"), Dim$Dim.2$quanti,
        c("Dim3.Correlation", "Dim3.p-Value"), Dim$Dim.3$quanti)

    write.table(output, paste("ACP_correlated_variables_", group.name,".csv", collapse=""), sep="\t", col.names = FALSE)


    ## Row clustering (adjust here distance/linkage methods to what you need!)
    hr <- hclust(as.dist(1-cor(group.data, method="pearson")), method="complete")
    group.data = t(group.data)

    ## Plot heatmap
    heatmap.2(
        group.data,
        Rowv=as.dendrogram(hr),
        Colv= FALSE,
        dendrogram = "row",
        scale="row",
        density.info="density",
        trace = "column",
        tracecol = "black",
        col = colorRampPalette(c("springgreen4","palegreen1","white","coral1","red4")),
        main = group.name,
        cexCol = 0.8,
        srtCol = 45)

    dev.print(file = paste("Heatmap_", group.name,".svg", collapse=""), device=svg)

}
