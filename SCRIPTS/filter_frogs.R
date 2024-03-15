# date : 2021/12/08
# version : 0.01
# author : Marie-Pierre Chapuis (CBGP)
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


## Filter data:
# - filter1: transform abundance data to null under a theshold estimated from alien controls (set with alien_in_samples = yes) or from an input value (set with alien_in_samples = no and rfa).
#   Controls for index contamination following correction 2 from Galan et al. 2016 (Tfa).
# - filter2: keep positive data only if congruent between replicates (2 or 3) and merge replicates
# - filter3: transform abundance data to null under a theshold estimated from negative controls.
#   Controls for extraction & pcr contamination following correction 1 from Galan et al. 2016 (Tcc).




###################
##               ##
##    INPUT      ##
## REQUIREMENTS  ##
##               ##
###################

## In the bin folder: the file 'Filtering_Functions.R' and the custom filtering functions.
## In the set folder: a metadata file and a abundance  tsv file.

## Specification of the metadata file:
# - a column 'observation_name' with sample names.
# - a column 'control' with 2 or 3 factors corresponding to non-control samples ("no"), negative samples ("negative") and to alien controls ("alien" ; optional).
# Below are fields not required for applying filters 1 and 2 but advised for next filters (e.g. filters on replicates - see file filter_project.R)
# - a column 'rep' with 2 or 3 levels of technical replicates. 
# - a column 'biological_unit', which corresponds to the names of the biological units, such as number of biological units x number of technical replicates = number of samples.
# - technical replicates of a same biological unit must follow each other.
# - a column 'project' to subset the metadata file by project


## Dependencies: packages to install are listed in the file 'Filtering_Functions'.

###################
##               ##
##  DATA IMPORT  ##
##               ##
###################

print("Data and params import")
## Load input files from arguments
args <- commandArgs(trailingOnly = TRUE)

## Load packages and custom filtering functions from bin folder
pathscript <-args[1] # SCRIPT directory
pathscript<-paste(pathscript,"/R_librairies_filtering/",sep="")
setwd(pathscript)
source("Filtering_Functions.R")
suppressPackageStartupMessages(library(tibble))

## Load input files from data folder
pathdata<-args[2] 
pathdata<-paste(pathdata,"/METABARCODING/",sep="") #DATA DIRECTORY
setwd(pathdata)

#the samples must be arranged by biological units (technical replicates of a same biological unit follow each other)
#if it is not the case, arrange them properly in your sample_data (metadata), then activate the following lines
# correct.order.samples<-rownames(sample)
# table<-table[,correct.order.samples]

METADATA<-args[3]
sample<-read.csv(METADATA,row.names=1,sep=";") 
sample<-sample_data(sample)


rawname<-args[4]
raw<-fread(rawname,sep="\t",header=TRUE)
name<-args[5] #name will be used to write names before extension of the different output tables. to modify if necessary (the character "_" is used to substract basename)
fragment<-args[6]
name<-paste(name,fragment,sep="")#name will be used to write names before extension of the different output tables. to modify if necessary (the character "_" is used to substract basename)
raw$observation_sum<-NULL
rownames(raw)<-raw$observation_name

seq<-raw[,9]
table<-raw[,(-(1:10)),with=FALSE]
all(colnames(table) %in% rownames(sample)) #check name consistency
all(rownames(sample) %in% colnames(table)) #check name consistency
setcolorder(table,rownames(sample))
table<-as.matrix(table)
rownames(table)<-raw$observation_name
table<-otu_table(table,taxa_are_rows = TRUE)

tax<-raw[,(1:9),with=FALSE]
tax<-as.matrix(tax)
rownames(tax)<-raw$observation_name
tax <- tax_table(tax)

data<-phyloseq(sample,table,tax)

# Transform all variables in sample_data to factors just in case...
df <- as.data.frame(lapply(sample_data(data),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
row.names(df) <- sample_names(data)
sample_data(data) <- sample_data(df)

data 
rm(raw)



###### CREATE R PRODUCTS DIRECTORY : 
RPRODUCT<-paste(pathdata,"Rproducts/",sep="")
if (!file.exists(RPRODUCT)) {
  # Create the directory if it does not exist
  dir.create(RPRODUCT, recursive = TRUE)}

###################
##               ##
##     FILTER    ##       
##    SETTINGS   ## 
##       &       ##
##   REDUCTION   ##
##               ##
###################


## Setting for filter1: transform abundance data to null under a theshold : index contamination estimated from alien controls
# If there are alien controls in the biom, then active the following line:
# alien_in_samples<-"yes" # to modify
# If there are no alien controls in the biom, then active the following line (and modify accordingly):
alien_in_samples<-"no" # to modify
rfa<-0.0002 # to modify

## Setting for filter2: keep positive data only if congruent between replicates
# 2 or 3 replicates must be included in data
nrep <-length(levels(factor(get_variable(data, varName = "rep")))) # to modify if the variable name of the column for replicates is not "rep"
#nrep
# here we consider that ALL technical replicates are congruent
if (nrep == 2) {ncon = 2} else if (nrep == 3) {ncon = 3} else {stop('nrep should be 2 or 3')}
# if nrep == 3 and you want at least 2 technical replicates to be congruent, please activate the next line
# ncon <- 2 # to modify
#ncon

## No setting required for filter3: transform abundance data to null under a theshold : extraction & pcr contamination estimated from negative controls
# but negative controls must be included in data



#########################
##                     ##
##  DATA PARTITIONING  ##
##                     ##
#########################

## Writing of a file for alien control samples and identification of the cluster_name of the alien control
if (alien_in_samples == "yes") {
  alien.data<-((subset_samples(data, control =="alien"))) # to modify if the variable name of the column for technical controls is not "control" and # to modify if the factor name of the column 'control' for alien samples is not "alien"
  
  condition <- function(x) { sum(x) > 0 }
  taxaToKeep <- filter_taxa(alien.data, condition)
  alien.data<-prune_taxa(taxaToKeep, alien.data)
  sort(taxa_sums(alien.data), TRUE)

  write.table(otu_table(alien.data),file=paste(RPRODUCT,name,"_abundance.alien.txt",sep=""),sep="\t")
  write.table(tax_table(alien.data),file=paste(RPRODUCT,name,"_tax.alien.txt",sep=""),sep="\t")
  
  alien2.data<-alien.data
  condition <- function(x) { sum(x) > 1000 }
  taxaToKeep <- filter_taxa(alien2.data, condition)
  taxa_alien_name<-taxa_names(prune_taxa(taxaToKeep, alien2.data))
  #taxa_alien_name
  
}

## Writing of a file for negative control samples
negative.data<-((subset_samples(data, control =="negative"))) # to modify if the variable name of the column for technical controls is not "control" and # to modify if the factor name of the column 'control' for negative samples is not "negative"
condition <- function(x) { sum(x) > 0 }
taxaToKeep <- filter_taxa(negative.data, condition)
negative.data<-prune_taxa(taxaToKeep, negative.data)
sort(taxa_sums(negative.data), TRUE)

write.table(otu_table(negative.data),file=paste(RPRODUCT,name,"_abundance.negative.txt",sep=""),sep="\t")
write.table(tax_table(negative.data),file=paste(RPRODUCT,name,"_tax.negative.txt",sep=""),sep="\t")

# Writing of a file for biological samples only (without control samples)
raw.data<-subset_samples(data,control=="no") # to modify if the variable name of the column for technical controls is not named "control" and and # to modify if the factor name of the column 'control' for non-control samples is not "no"
condition <- function(x) { sum(x) > 0 }
taxaToKeep <- filter_taxa(raw.data, condition)
raw.data <- prune_taxa(taxaToKeep, raw.data)
sort(taxa_sums(raw.data),decreasing = TRUE)
write.table(otu_table(raw.data),file=paste(RPRODUCT,name,"_abundance.raw.txt",sep=""),sep="\t")
write.table(tax_table(raw.data),file=paste(RPRODUCT,name,"_tax.raw.txt",sep=""),sep="\t")

###########################
##                       ##
##    DATA FILTERING     ##
##                       ##
###########################

print("Data Filtering")
## Filter1: transform data to null under a theshold : index contamination estimated from alien controls
# correction 2 from Galan et al. 2016 (Tfa)
print('Filter 1/3')
if (alien_in_samples == "yes") {
  otu.clone<-get_sample(data, i = taxa_alien_name) 
  no.sample.clone<-get_sample(subset_samples(tcc.data, control !="alien"), i = taxa_alien_name) # to modify if the variable name of the column for technical controls is not "control" and and # to modify if the factor name of the column 'control' for alien samples is not "alien"
  rfa<-max(no.sample.clone) / sum(otu.clone)
  max(no.sample.clone)
  sum(otu.clone)
}
#rfa
f1.data<-filter_tfa(physeq = data, ratio = rfa)

# eliminate alien control samples
f1.data<-subset_samples(f1.data,control !="alien") # to modify if the variable name of the column for technical controls is not "control" and and # to modify if the factor name of the column 'control' for alien samples is not "alien"

# update of files after filter1
condition <- function(x) { sum(x) > 0 }
taxaToKeep <- filter_taxa(f1.data, condition)
f1.data <- prune_taxa(taxaToKeep, f1.data)
sort(taxa_sums(f1.data),decreasing = TRUE)
write.table(otu_table(f1.data),file=paste(RPRODUCT,name,"_abundance.filter1.txt",sep=""),sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)
write.table(tax_table(f1.data),file=paste(RPRODUCT,name,"_tax.filter1.txt",sep=""),sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)


print('Filter 2/3')
tryCatch({
## Filter2: keep positive data only if congruent between replicates
condition <- function(x) { sum(x) > 0 }
taxaToKeep <- filter_taxa(f1.data, condition)
f1.data <- prune_taxa(taxaToKeep, f1.data)

# pre-requisite :  prune samples without replicate (in case one replicate failed)
filter<-(as.logical(plyr::count(sample_data(f1.data)$biological_unit)[,2]-(nrep-1))) # to modify if the variable name of the column for biological units is not "biological_unit"
nbrep<-(plyr::count(sample_data(f1.data)$biological_unit)[,2]) # to modify if the variable name of the column for biological units is not "biological_unit"
samplesToKeep<-rep(filter,times = nbrep, length.out = nrow(sample_data(f1.data)))#dim = number of samples
f2.data<-prune_samples(samplesToKeep,f1.data)

# keep congruent data
# to be applicable, the variable name of the column for replicates should be "rep"
# to be applicable, the number of replicates should be either 2 or 3
# to be applicable, the number of congruent replicates should be either 2 (for 2 or 3 replicates) or 3 (for 3 replicates) 
f2.data<-filter_rep(physeq = f2.data, nb_rep = nrep, nb_con = ncon)

# update of files after filter2
condition <- function(x) { sum(x) > 0 }
taxaToKeep <- filter_taxa(f2.data, condition)
f2.data <- prune_taxa(taxaToKeep, f2.data)
sort(taxa_sums(f2.data),decreasing = TRUE)

write.table(otu_table(f2.data),file=paste(RPRODUCT,name,".abundance.filter2.txt",sep=""),sep="\t")
write.table(tax_table(f2.data),file=paste(RPRODUCT,name,".tax.filter2.txt",sep=""),sep="\t")
},error = function(e) {
  message("ERROR IN DUPLICATED SAMPLES: One or several samples are not duplicated. The R filtering step cannot be finalised. Please check your samples ","\n",
  "When the problems linked to the unduplicated samples have been resolved, before re-executing the HAMI pipeline, make sure you have deleted the old files produced by the HAMI pipeline")
})


## update of files (including metadata), by merging of replicates (by sum)
merged.data <- merge_samples(f2.data, "biological_unit") # to modify if the variable name of the column for biological units is not named "biological_unit"
otu_table(merged.data) <- t(otu_table(merged.data))
x<-sample_data(merged.data)
dx<-data.frame(x)
y<-sample_data(f2.data)
y<-y[y$rep=="1",]
dy<-data.frame(y)
for (i in 1:length(x)) {
  dx[,i]<-as.factor(dx[,i]) 
  dy[,i]<-as.factor(dy[,i]) 
  levels(dx[,i])<-levels(dy[,i])
}

sample_data(merged.data)<-dx
write.table(dx,file=paste(RPRODUCT,name,".merge_metadata.txt",sep=""),sep="\t")
condition <- function(x) { sum(x) > 0 }
taxaToKeep <- filter_taxa(merged.data, condition)
merged.data <- prune_taxa(taxaToKeep, merged.data)
sort(taxa_sums(merged.data),decreasing = TRUE)
write.table(otu_table(merged.data),file=paste(RPRODUCT,name,".abundance.filter2.merge.txt",sep=""),sep="\t")
write.table(tax_table(merged.data),file=paste(RPRODUCT,name,".tax.filter2.merge.txt",sep=""),sep="\t")

## Filter3: transform abundance data to null under a theshold : extraction & pcr contamination estimated from negative controls
# correction 1 from Galan et al. 2016 (Tcc)
print('Filter 3/3')
neg <- subset_samples(merged.data, control == "negative") # to modify if the factor name of the column 'control' for negative samples is not "negative"
# show contaminants
sort(taxa_sums(neg), TRUE)[1:15]
mean(taxa_sums(neg))
min(taxa_sums(neg))
max(taxa_sums(neg))
sum(taxa_sums(neg))
quantile(taxa_sums(neg),probs=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99,0.999,0.9999))

f3.data<-filter_tcc(physeq = merged.data, phyneg = neg) # the variable name of the column for technical controls should be "control"

# eliminate negative control samples 
f3.data<-subset_samples(f3.data, control !="negative") # to modify if the variable name of the column for technical controls is not "control" and # to modify if the factor name of the column 'control' for negative samples is not "negative"
dataframe <-as.data.frame(otu_table <- otu_table(f3.data))
dataframe <- rownames_to_column(dataframe, var = "cluster")
write.table(dataframe,file=paste(name,"_abundance.filter3.csv",sep=""),sep=";",row.names = FALSE)
dataframe <-as.data.frame(tax_table <- tax_table(f3.data))
dataframe <- rownames_to_column(dataframe, var = "cluster")
write.table(dataframe,file=paste(name,"_tax.filter3.csv",sep=""),sep=";",row.names = FALSE) 

print('End Data Filtering')


###################
##               ##
##    SUMMARY    ##
##     FILES      ##
##               ##
###################

## Summary file for data filtering
if (alien_in_samples == "yes")  {
  nsamples<-c(nsamples(data),nsamples(raw.data),nsamples(alien.data),nsamples(negative.data),nsamples(f1.data),nsamples(merged.data),nsamples(f3.data))
  ncluster<-c(ntaxa(data),ntaxa(raw.data),ntaxa(alien.data),ntaxa(negative.data),ntaxa(f1.data),ntaxa(merged.data),ntaxa(f3.data))
  nread<-c(sum(taxa_sums(data)),sum(taxa_sums(raw.data)),sum(taxa_sums(alien.data)),sum(taxa_sums(negative.data)),sum(taxa_sums(f1.data)),sum(taxa_sums(merged.data)),sum(taxa_sums(f3.data)))
  df<-data.frame(nsamples, ncluster, nread)
  rownames(df)<-c("all","samples","alien","negative","samples.filter1","samples.filter2","samples.filter3")
  df
  write.table(df,file=paste(RPRODUCT,name,".filters.summary.txt",sep=""),sep="\t")
  cat(paste("Log = file:",name,",number_aliens_filter1:",nsamples(alien.data),",rfa_filter1:",rfa,sep="",",number_negatives_filter2:",nsamples(negative.data)),file=paste(RPRODUCT,name,".filters.summary.txt",sep=""),append=TRUE)
}
if (alien_in_samples == "no")  {
  nsamples<-c(nsamples(data),nsamples(raw.data),nsamples(negative.data),nsamples(f1.data),nsamples(merged.data),nsamples(f3.data))
  ncluster<-c(ntaxa(data),ntaxa(raw.data),ntaxa(negative.data),ntaxa(f1.data),ntaxa(merged.data),ntaxa(f3.data))
  nread<-c(sum(taxa_sums(data)),sum(taxa_sums(raw.data)),sum(taxa_sums(negative.data)),sum(taxa_sums(f1.data)),sum(taxa_sums(merged.data)),sum(taxa_sums(f3.data)))
  df<-data.frame(nsamples, ncluster, nread)
  rownames(df)<-c("all","samples","negative","samples.filter1","samples.filter2","samples.filter3")
  
  df
  write.table(df,file=paste(RPRODUCT,name,".filters.summary.txt",sep=""),sep="\t")
  cat(paste("Log = file:",name,",rfa_filter1:",rfa,sep="",",number_negatives_filter2:",nsamples(negative.data)),file=paste(RPRODUCT,name,".filters.summary.txt",sep=""),append=TRUE)
}

print("End Summary")


###################
##               ##
##    QUALITY    ##
##    CHECKS     ##
##               ##
###################

## Check quality of replicates
# computation of R2 on the number of sequences for all OTUs = take a long time...
datachecks<- f2.data
rep1<-subset_samples(datachecks, rep %in% c("1"))
rep2<-subset_samples(datachecks, rep %in% c("2"))

allrep<-matrix(NA,nrow=(ntaxa(datachecks)*0.5*(nsamples(datachecks))),ncol=2)
for (i in 1:(ntaxa(datachecks))){
  nreads.rep1<-apply((otu_table(rep1)[i]),1,identity)
  nreads.rep2<-apply((otu_table(rep2)[i]),1,identity)
  rep<-cbind(nreads.rep1,nreads.rep2)
  allrep<-rbind(allrep,rep)
}
allrep<-data.frame(allrep)
names(allrep)<-c("nreads.rep1","nreads.rep2")
model<-(lm(allrep$nreads.rep1~allrep$nreads.rep2))
anova(model)
Rcarre<-anova(model)[1,2]/(anova(model)[1,2]+anova(model)[2,2])
p<-ggplot(allrep, aes(x=nreads.rep1,y=nreads.rep2))
p<-p+expand_limits(y=0)+expand_limits(x=0)
p<-p+ scale_x_log10()+scale_y_log10()
p<-p + geom_point()
p<-p + geom_abline(slope=1, intercept=0)
p<-p+ggtitle(paste("R2 = ",round(Rcarre,3)))


ggsave(filename = paste(RPRODUCT, name, ".nreads_by_replicate.jpeg", sep = ""),
       plot = p,
       width = 6,
       height = 6,
       units = "in",
       dpi = 600)


print('End')