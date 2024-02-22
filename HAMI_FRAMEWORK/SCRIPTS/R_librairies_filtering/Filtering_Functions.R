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

library(phyloseq)# to install phyloseq source, use 'install_github("joey711/phyloseq")' which requires the devtools package

## ----load filtering fonctions----
source("./Tccfilter.R")
source("./Tfafilter.R")
source("./Repfilter2.R")
