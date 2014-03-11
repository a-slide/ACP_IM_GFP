#! /usr/bin/Rscript

options(width = 160) # wider terminal size
library(FactoMineR)
library(RColorBrewer)
library(gplots)
x11()

########################################################################################################################
#   2 arguments are required
#       1 = filename containing the data
#       2 = Max Percent of NA allowed in data
########################################################################################################################

# Inport data in a dataframe
filename = commandArgs(TRUE)[1]
data = read.table(filename, dec=".", sep="\t", header = TRUE, row.names=1)

# Create a list of group names to analyze data separatly
groups = c("J7", "J30", "J90", "Eutha", "Wilson", "Sheperd", "Monk", "Walcott","Ruben","Woody")

eliminated_genes = list()
retained_genes = list()

write ("RETAINED GENES", file = "Retained_genes.txt")
write ("ELIMINATED GENES", file = "Eliminated_genes.txt")

for (group in groups){
    
    print (paste("ANALYSING GROUP ", group))
    
    #Extract informations for the current group
    group.data = data[grep(group, rownames(data), ignore.case = TRUE),]

    print ("INDIVIDUAL ANALYSED")
    print (rownames(group.data))
    
    # Calculate the maximal number of NA allowed
    NA_max_percent = as.numeric(commandArgs(TRUE)[2])
    NA_max = round ((NA_max_percent*nrow(group.data)/100), digits=0)
 
    print (paste("MAX NA ALLOWED BY GROUP = " , NA_max))
    
    # List genes with more NA value than max NA and export to file
    eliminated = paste(names(group.data[,colSums(is.na(group.data)) > NA_max]), collapse = " ")
    eliminated_genes = c(eliminated_genes, group, eliminated)
    
    # List genes with less NA value than max NA and export to file
    retained = paste((names(group.data[,colSums(is.na(group.data)) <= NA_max])), collapse = " ") 
    retained_genes = c(retained_genes, group, retained)
    
    # Remove columns with more NA than max_N and replacing the remaining values by 0
    group.data = group.data[,colSums(is.na(group.data)) <= NA_max]
    group.data[is.na(group.data)] = 0
    
    print ("VARIABLE ANALYSED ")
    print(colnames(group.data))
    
    # Perform ACP with FactoMineR pakage
    PCA.res = PCA(
        group.data,
        scale.unit=TRUE,
        graph = FALSE)

    # Create and export plots for all ACP results
    plot(PCA.res, axes = c(1, 2), choix = "ind", title = group)
    dev.print(file = paste("ACP_individus_", group,".svg", collapse=""), device=svg)

    plot(PCA.res, axes = c(1, 2), choix = "var", title = group, cex = 0.5)
    dev.print(file = paste("ACP_variables_", group,".svg", collapse=""), device=svg)

    # Create summary of highly corelated variables with dim1, dim2 and export to CSV files

    Dim = dimdesc(PCA.res, axes = c(1, 2))
    
    output = rbind (
        c("Dim1.Correlation", "Dim1.p-Value"), Dim$Dim.1$quanti,
        c("Dim2.Correlation", "Dim2.p-Value"), Dim$Dim.2$quanti)
        
    write.table(output, paste("ACP_correlated_variables_", group,".csv", collapse=""), sep="\t", col.names = FALSE)
    
    ## Row clustering
    dist_mat = as.dist(1-cor(group.data, method="pearson"))
    hr <- hclust(dist_mat, method="complete")
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
        main = group,
        cexCol = 0.7,
        cexRow = 0.5)
        
    dev.print(file = paste("Heatmap_", group,".svg", collapse=""), device=svg)
}

# Write the list of eliminated or retained genes in files
lapply(eliminated_genes, write, "Eliminated_genes.txt", append=TRUE)
lapply(retained_genes, write, "Retained_genes.txt", append=TRUE)

warnings()
