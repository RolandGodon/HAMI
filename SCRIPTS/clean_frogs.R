# date : 2021/12/08
# version : 0.01
# authors : Sylvain Piry (CBGP)
# licence : GPL

# Update_date : 2024/02/20
# By : Benoit Penel (CBGP)
# Update_name : HAMI_FRAMEWORK 

rm(list=ls())# Reset

#####################
##                 ##
##      SHORT      ##
##   DESCRIPTION   ##
##                 ##
#####################

## Work on raw abundance and multihits text files from Frogs.

## Correct the taxonomy of each cluster of the abundance tsv file tagged 'multi-affiliation', from decision rules applied on the cleaned multihits tsv file:

# Clean the taxonomy (7 ranks: Kingdom,Phylum,Class,Order,Family,Genus,Species) of each hit of the two files for:
# - the special characters in all 7 taxonomic ranks.
# - 'unidentified', and 'sp.'  which are replaced by 'unknown species' in the Species rank .
# - unexpected number of words in all 7 taxonomic ranks.
# - for the species and genus ranks only: if the number of 'unknown species' hits < % threshold then ignore them.
# - for all taxonomic ranks: 
  # - if there is a unique taxonomic name, then rename by this taxonomic name.
  # - if there is a unique taxonomic name and an 'unknown', then rename by 'unknown'.
  # - if there are two or more taxonomic names, then name 'multi-affiliation'.

## Identify chimeric sequences using the isBimeraDenovo(...) function in the dada2 library

## Export an abundance table tsv file cleaned for taxonomy and filtered for chimeras.

###################
##               ##
##    INPUT      ##
## REQUIREMENTS  ##
##               ##
###################

## In the set folder: a multihit tsv file and an abundance tsv file.
## Dependencies: packages to install are listed below.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
options(repos = c(CRAN = "http://cran.univ-lyon1.fr/"))
if (!requireNamespace("BiocManager", quietly = TRUE))  {
  print("Package installation in progress, this may take a couple of minutes")
  install.packages("BiocManager",quiet=TRUE)
  }
if (!require("dada2", quietly = TRUE)){ 
  print("Package installation in progress, this may take a couple of minutes")
  BiocManager::install(version = '3.18',quiet=TRUE)
  BiocManager::install("dada2", version = "3.18",quiet=TRUE)
}
suppressPackageStartupMessages(library(dada2)) # BioConductor


###################
##               ##
##  DATA IMPORT  ##
##               ##
###################

## Load input files from arguments
args <- commandArgs(trailingOnly = TRUE)
pathdata<-args[1] # data directory
setwd(pathdata) 


multi_namefile<-args[2] #path to multihits csv file
abundance_namefile<-args[3] #path to abundance csv file

name<-args[4] #name of the project, it will be used to import files and write names before extension of the different output tables.
fragment<-args[5] # name of the DNA fragment sequenced,  it will be used to import files and write names before extension of the different output tables.
threads<-args[6] #threads used for chimera deletion

mt<-fread(multi_namefile,sep="\t",header=TRUE) #open multihits file
mh <- data.table(newtaxaname=mt$blast_taxonomy, blast_taxonomy=mt$blast_taxonomy, observation_name=mt$`#observation_name`)
rm(mt)

at<-fread(abundance_namefile,sep="\t",header=TRUE)# open the abundance file
ab <- data.table(newtaxaname=at$`blast_taxonomy`, blast_taxonomy=at$`blast_taxonomy`, observation_name=at$observation_name)


###### CREATE CHIMERAS DIRECTORY : 
RPRODUCT<-paste(pathdata,"/chimeras/",sep="")
if (!file.exists(RPRODUCT)) {
  # Create the directory if it does not exist
  dir.create(RPRODUCT, recursive = TRUE)}

###################
##               ##
##  CLEANING OF  ##
##      THE      ##
##   MULTIHITS   ##
##     FILE      ##
##               ##
###################

print("cleaning of the taxonomic nomenclature")
## Cleaning for the special characters []()_"0123456789.
mh$newtaxaname<-gsub("[\\[\\]()\"0-9]","",mh$newtaxaname,perl=TRUE)
mh$newtaxaname<-gsub(" +;",";",mh$newtaxaname,perl=TRUE)

## Deconstruction of the taxonomic column ('blast_taxonomy') in 7 taxonomic fields
taxa_names<- mh[ , "newtaxaname"] %>% separate(newtaxaname,c("tx1","tx2","tx3","tx4","tx5","tx6","tx7"),";")
mh<-cbind(taxa_names,mh[, c("observation_name" ,"blast_taxonomy")]) # observation name is the cluster ID and useful to link with the abundance tsv file
rm(taxa_names)



## Replacement in the species field of 'unidentified' anywhere in the field by 'unknown species' 
## Replacement in the species field of 'sp' and 'sp.' at the end of the field by 'unknown species'
## Replacement in the species field pattern similar to 'Genus_xx._specie' by  'unknown species'
## Replacement of the species field by 'unknown species' if empty

mh$tx7<-gsub("^.*(unidentified).*$","unknown species",mh$tx7,perl=TRUE)
mh$tx7<-gsub("^.*(sp)\\..*$","unknown species",mh$tx7,perl=TRUE)
mh$tx7<-gsub("^.*(sp)$","unknown species",mh$tx7,perl=TRUE)
mh$tx7<-gsub("^[A-Z][a-z]+_[a-z\\.]+_[a-z\\-]+$","unknown species",mh$tx7,perl=TRUE)
mh$tx7[grepl("^ *$",mh$tx7)]<-"unknown species"

## Delete in the species field upper letter at the end of species name

mh$tx7<-gsub("AS$","",mh$tx7)
mh$tx7<-gsub("DHJ$","",mh$tx7)
mh$tx7<-gsub("-GLA$","",mh$tx7)
mh$tx7<-gsub("_L.$","",mh$tx7)
mh$tx7<-gsub("_CHU$","",mh$tx7)

### This line allow to found other wrong pattern that should be add to the previous list 
# erreur=mh[!grepl("^[A-Z][a-z]+_[a-z\\-]+$", mh$tx7) & !(mh$tx7=="unknown species"),]



## Ignoring of 'unknown species' or 'unknown genus' under a threshold
thresh=0.5 # % threshold below which 'unknown' are ignored 
mh$clus <- as.integer(gsub("Cluster_","",mh$observation_name))
mh$blast_taxonomy<-NULL

my.genus.result <- as.data.table(unique(mh %>% add_count(clus, tx1, tx2, tx3, tx4, tx5, tx6, observation_name)))
my.genus.result$ntot<-as.data.table(unique(mh %>% add_count(observation_name)))$n
setkey(my.genus.result, clus, n)
my.genus.result$prop.genus <-my.genus.result$n/my.genus.result$ntot
my.genus.result$ignore.genus <-FALSE
my.genus.result$ignore.genus <- my.genus.result$ignore.genus | (grepl("unknown",my.genus.result$tx6,perl=TRUE,ignore.case=TRUE) & my.genus.result$prop.genus <thresh)
#nrow(my.genus.result)
#names(my.genus.result)

my.species.result <- as.data.table(unique(mh %>% add_count(clus, tx1, tx2, tx3, tx4, tx5, tx6, tx7, observation_name)))
my.species.result$ntot<-as.data.table(unique(mh %>% add_count(observation_name)))$n
setkey(my.species.result, clus, n)
my.species.result$prop.species<-my.species.result$n/my.species.result$ntot
my.species.result$ignore.species<-FALSE
my.species.result$ignore.species <- my.species.result$ignore.species | (grepl("unknown",my.species.result$tx7,perl=TRUE,ignore.case=TRUE) & my.species.result$prop.species<thresh)
#nrow(my.species.result)
#names(my.species.result)

my.first.result<-merge(my.genus.result, my.species.result, by = names(my.species.result)[1:8])[,c(1:8,13,18)]

my.second.result <- data.table(observation_name=unique(my.first.result$observation_name))
for(j in 1:7) {
  tx <- as.character(as.data.frame(my.first.result)[,j])
  tmp <- apply(my.second.result, 1, function(r){
    v <- unique(tx[(j==6 | !my.first.result$ignore.genus) & my.first.result$observation_name==r[1]])
    v <- unique(tx[(j==7 | !my.first.result$ignore.species) & my.first.result$observation_name==r[1]])
    if (length(v)==1) {v[1]} 
    else if ((length(v)==2)) {
      unk <- grepl("unknown",v)
      unkV <- unique(v[unk])
      if (length(unkV) == 0) {
        "Multi-affiliation"
      } else {
        unkV[1]
      }
    } else {
      "Multi-affiliation"
    }
  })
  my.second.result <- cbind(my.second.result, tmp)
}

#nrow(my.second.result)

## Reconstruction of the new full taxonomy from each taxonomic field for each cluster
my.mh.result <- data.table(observation_name=my.second.result$observation_name,  blast_taxonomy= apply(my.second.result, 1, function(r) {
  paste(r[2:8], sep=";", collapse=";")
}))
#nrow(my.mh.result)

###################
##               ##
##  CLEANING OF  ##
##      THE      ##
##   ABUNDANCE   ##
##     FILE      ##
##               ##
###################

## Cleaning for some special characters
ab$newtaxaname<-gsub("[\\[\\]()\"0-9]","",ab$newtaxaname,perl=TRUE)
ab$newtaxaname<-gsub(" +;",";",ab$newtaxaname,perl=TRUE)

## Deconstruction of the taxonomic column in 7 taxonomic fields
ab  <- ab[ab$newtaxaname!="no data"]
taxa_names<- ab[ , "newtaxaname"] %>% separate(newtaxaname,c("tx1","tx2","tx3","tx4","tx5","tx6","tx7"),";")
ab<-cbind(taxa_names,ab[, c("observation_name" ,"blast_taxonomy")])
rm(taxa_names)




## Replacement in the species field of 'unidentified'  anywhere in the field by 'unknown species' 
## Replacement in the species field of 'sp' and 'sp.' at the end of the field by 'unknown species'
## Replacement in the species field pattern similar to 'Genus_xx._specie' by  'unknown species'
## Replacement of the species field by 'unknown species' if empty
ab$tx7<-gsub("^.*(unidentified).*$","unknown species",ab$tx7,perl=TRUE)
ab$tx7<-gsub("^.*(sp)\\..*$","unknown species",ab$tx7,perl=TRUE)
ab$tx7<-gsub("^.*(sp)$","unknown species",ab$tx7,perl=TRUE)
ab$tx7<-gsub("^[A-Z][a-z]+_[a-z\\.]+_[a-z\\-]+$","unknown species",ab$tx7,perl=TRUE)
ab$tx7[grepl("^ *$",ab$tx7)]<-"unknown species"

## Delete in the species field upper letter at the end of species name

ab$tx7<-gsub("AS$","",ab$tx7)
ab$tx7<-gsub("DHJ$","",ab$tx7)
ab$tx7<-gsub("-GLA$","",ab$tx7)
ab$tx7<-gsub("_L.$","",ab$tx7)
ab$tx7<-gsub("_CHU$","",ab$tx7)


### This line allow to found other wrong pattern that should be add to the previous list 
#erreur=ab[!grepl("^[A-Z][a-z]+_[a-z\\-]+$", ab$tx7) & !(ab$tx7=="unknown species") & !(ab$tx7=="Multi-affiliation"),]


# reconstruction of the new full taxonomy from each taxonomic field
my.ab.result <- data.table(observation_name=ab$observation_name,  blast_taxonomy= apply(ab, 1, function(r) {
  paste(r[1:7], sep=";", collapse=";")
}))

########################
##                    ##
##  WRITING OF A NEW  ##
##  TAXONOMY IN A NEW ##
##   ABUNDANCE FILE   ##
##                    ##
########################

my.result <- merge(my.mh.result, my.ab.result,by="observation_name",all.y=TRUE)
final <- merge(my.result, at,by="observation_name",all.y=TRUE)
final$blast_taxonomy <- apply(final, 1, function(r){
  if (! is.na(r["blast_taxonomy.x"])) {r["blast_taxonomy.x"]} else {r["blast_taxonomy.y"]}
})
final$blast_taxonomy.x<-NULL
final$blast_taxonomy.y<-NULL
setcolorder(final,c(2,3,4,5,6,7,8,9,10,1,11:(ncol(final))))
fordada<-final[order(as.integer(gsub("Cluster_","",final$observation_name)))]
fordada<-as.data.frame(fordada)

####################
##                ##
##     REMOVE     ##
##    CHIMERAS    ##
##                ##    
####################

print("Remove chimeras")
## Tune data structure
samples <- colnames(fordada)[seq(1+which(colnames(fordada)=="observation_sum"), ncol(fordada))]
fordada.m <- data.matrix(fordada[ , seq(1+which(colnames(fordada)=="observation_sum"), ncol(fordada))])
rownames(fordada.m) <- fordada$seed_sequence
fordada.t <- t(fordada.m)

##
options(mc.cores = threads)
bimeras.v <- isBimeraDenovo(getUniques(fordada.t), minFoldParentOverAbundance=4, verbose=TRUE, multithread = TRUE)
bimeras.df <- data.frame(x=as.logical(bimeras.v), seq=names(bimeras.v))
write.table(fordada[bimeras.df$x==FALSE,], file = paste(name,fragment,"_cleaned_abundance.txt",sep=""), sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)#writing of filtered abundance table (without identified chimeras)
write.table(fordada[bimeras.df$x==TRUE,], file = paste(RPRODUCT,name,fragment,"_chimeras_list.txt",sep=""), sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)#downloading of identified chimeras

## Summary file
number_of_total_clusters = nrow(fordada)
number_of_bimeras = nrow(fordada[bimeras.df$x==TRUE,])
number_of_nobimeras = nrow(fordada[bimeras.df$x==FALSE,])
summary.df = data.frame(number_of_total_clusters, number_of_bimeras, number_of_nobimeras)
write.table(summary.df,file=paste(RPRODUCT,name,fragment,"_chimeras_summary.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

