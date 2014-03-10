#! /usr/bin/Rscript

options(width = 160) # wider terminal size
library(FactoMineR)
library(RColorBrewer)
library(gplots)
x11()

########################################################################################################################
#   3 arguments are required
#       1 = filename containing the data
#       2 = Filtering NA contining genes for global data (TRUE) or group data (FALSE)
#       3 = Max Percent of NA allowed in data
########################################################################################################################

# Inport data in a dataframe
filename = commandArgs(TRUE)[1]
data = read.table(filename, dec=".", sep="\t", header = TRUE, row.names=1)

# Parse the second parameter to set up global_NA_filtering variable
global_NA_filtering = as.logical(commandArgs(TRUE)[2])

# if TRUE NA filtering is done on the global dataset
 if (global_NA_filtering == TRUE){
    
    # Calculate the maximal number of NA allowed
    NA_max_percent = as.numeric(commandArgs(TRUE)[3])
    NA_max = round ((NA_max_percent*nrow(data)/100), digits=0)
    
    print ("GLOBAL DATA FILTERING")
    print (paste("MAX NA ALLOWED = " , NA_max))
    
    # List genes with more NA value than max NA and export to file
    eliminated_genes = names(data[,colSums(is.na(data)) > NA_max])
    write(eliminated_genes, "Eliminated_genes.txt")
    print (paste(length(eliminated_genes), " GENES HAVE BEEN FILTERED OUT"))
    
    # List genes with less NA value than max NA and export to file
    retained_genes = names(data[,colSums(is.na(data)) <= NA_max])
    write(retained_genes, "Retained_genes.txt")
    
    # Remove columns with more NA than max_N and replacing the remaining values by 0
    data = data[,colSums(is.na(data)) <= NA_max]
    data[is.na(data)] = 0
}

# Create a list of row number to analyze data separatly
groups = c("J7", "J30", "J90", "Euthanasia", "Wilson", "Sheperd", "Monk", "Walcott","Ruben","Woody")

l = list(   1:6,                          # Line range for J7
            13:18,                        # Line range for J30
            25:29,                        # Line range for J90
            35:40,                        # Line range for eutha
            c(1,7,13,19,25,30,35,41),     # Line range for Wil
            c(2,8,14,20,26,31,36,42),     # Line range for Shep
            c(3,9,15,21,27,32,37,43),     # Line range for Monk
            c(4,10,16,22,28,33,38,44),    # Line range for Walc
            c(5,11,17,23,29,34,39,45),    # Line range for Rub
            c(6,12,18,24,40,46))          # Line range for Woo

names(l) = groups

#Loop foreach element in l

for (i in seq_along(l)){
    
    #Extract informations for the current group
    group.name = groups[i]
    group.index = l[[i]]
    group.data = data[group.index,]

    # if FALSE NA filtering is done on the independant group dataset
     if (global_NA_filtering == FALSE){
         
        # Calculate the maximal number of NA allowed
        NA_max_percent = as.numeric(commandArgs(TRUE)[3])
        NA_max = round ((NA_max_percent*nrow(group.data)/100), digits=0)
     
        print (paste ("DATA FILTERING FOR GROUP", group.name))
        print (paste("MAX NA ALLOWED BY GROUP = " , NA_max))
        
        # List genes with more NA value than max NA and export to file
        eliminated_genes = names(group.data[,colSums(is.na(group.data)) > NA_max])
        write(eliminated_genes, paste("Eliminated_genes_", group.name,".txt", collapse=""))
        
        # List genes with less NA value than max NA and export to file
        retained_genes = names(group.data[,colSums(is.na(group.data)) <= NA_max])
        write(retained_genes, paste("Retained_genes_", group.name,".txt", collapse=""))
        
        # Remove columns with more NA than max_N and replacing the remaining values by 0
        group.data = group.data[,colSums(is.na(group.data)) <= NA_max]
        group.data[is.na(group.data)] = 0
    }
    
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
        main = group.name,
        cexCol = 0.8,
        srtCol = 45)

    dev.print(file = paste("Heatmap_", group.name,".svg", collapse=""), device=svg)

}
