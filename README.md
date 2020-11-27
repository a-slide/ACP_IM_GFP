# TaqMDA

**Multi-dimensional analyse for RT QPCR data obtained from Taqman QPCR arrays**

---

**Creation : 2015/08/28**

**Last update : 2020/11/27**

---

## Motivation

TaqMDA is an R script developed to analyze large scale data obtained by RT-QPCR arrays. It generates a principal component analysis for the first 3 dimensions and a heatmap.

Developed for the project LR-GFP-NHP

## Principle

TaqMDA.R is a simple R script which analyzes QPCR values from a CSV file (see example in *data/* directory).

For each qPCR target, if more than N % (defined by the user) of the values are undetermined, this will trigger its exclusion from the analysis. 

The lists of genes containing NA, eliminated and retained are written in a text file.

The remaining qPCR data are reduced, centered and classified by hierarchical cluster analysis.

A principal component analyze is performed with FactoMineR package and summarized by graphs and tables in both individual and variable spaces.

Finally, a distance matrix is generated by the Pearson method and subsequently used to create a heatmap scaled by individual sample thanks to gplots package.

## Get TaqMDA

* Download the program

With *git clone* from the repository

```	git clone https://github.com/a-slide/TaqMDA```

From the [zip archive](https://github.com/a-slide/TaqMDA/archive/master.zip) 

* Set the script to be executable

	``	sudo chmod u+x TaqMDA.R```

* Add the executable to your system PATH

## Usage

**Format your tab separated csv data sheet as follow:**

```  
            qPCR1   qPCR2   qPCR3   ...
sample1-1   0.23    0       12334   ...
sample1-2   12      0       8755    ...
sample1-3   65      0       87545   ...
...         ...     ...     ...     ...
```
See example file in the repository data folder

* Replace non detected qPCR values by NA (if not already done)

* Define a percentage of autorized NA values per gene to be considered as valid

* Run the script with the data sheet and the percent of allowed NA values

```
    USAGE: TaqMDA.R  <Sample_sheet.csv>  <max_percent_of_na_authorized>
        * Filename of csv file containing the data
        * Max percentage of values with NA allows per column

    EXAMPLE : TaqMDA.R ./data/Sample.csv 50
```

All files will be generated in the current folder

## Dependencies:

The program was developed under Linux Mint 16 "petra" but is compatible with other LINUX debian based distributions.

#### Program

* R 3.2.2 +

#### Third party R packages

* FactoMineR
* RColorBrewer
* gplots

# Data and Results

* [Biceps Femoris Data](https://raw.githubusercontent.com/a-slide/TaqMDA/master/Paper_Data/Gernoux_RawDatas_84_genes_BicepsFemoris.csv)
* [Tibialis Data](https://raw.githubusercontent.com/a-slide/TaqMDA/master/Paper_Data/Gernoux_RawDatas_84_genes_Tibialis.csv)
* [All results](https://raw.githubusercontent.com/a-slide/TaqMDA/master/Paper_Results/Whole analysis of inflammation in muscle.pdf)

## Authors and Contact

Adrien Leger - 2020
* Adrien Leger - <aleg {at} ebi.ac.uk>
* [Github](https://github.com/a-slide)

[Translational gene therapy for genetic diseases - UMR 1089](https://umr1089.univ-nantes.fr/)
