## ----load-packages-silently, cache = FALSE, message = FALSE, echo = FALSE----
library(ggplot2)
library(plyr)
library(reshape2)
library(gridExtra)
library(ape)
library(scales)
library(matrixStats)
library(data.table)
library(vegan)


# To install phyloseq :
options(repos = c(CRAN = "http://cran.univ-lyon1.fr/"))
if(!requireNamespace("phyloseq", quietly = TRUE)){  #use 'install_github("joey711/phyloseq")' which requires the devtools package
    print("Package installation in progress, this may take a couple of minutes")
    if(!requireNamespace("phyloseq", quietly = TRUE)){ 
        options(repos = c(igraph = 'https://igraph.r-universe.dev',CRAN = 'https://cloud.r-project.org' ))
        install.packages('igraph')}
    options(repos = c(CRAN = "http://cran.univ-lyon1.fr/"))
    devtools::install_github("joey711/phyloseq",quiet=TRUE)}
library(phyloseq)

## ----load filtering fonctions----
source("./Tccfilter.R")
source("./Tfafilter.R")
source("./Repfilter2.R")
