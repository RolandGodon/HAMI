#!/usr/bin/env python
# -*-coding: utf-8 -*


'''

Author : Benoit Penel
Date : 28/08/2023
version : 1.0
PYTHON3

The aim of this script is to discriminate Metabarcoding, Barcoding and control samples which were sequenced at the same time in a context of HAMI FRAMEWORK.
In practical, we separate  Metabarcoding, Barcoding and control  data from the occurence file previously produce with frogs, using the sample name nomenclature .

'''

###### PACKAGES :

import sys
import pandas as pd
import os 
import re
import csv

###### PARAM YOUR PROJECT :
#Directory path
outputdir= sys.argv[2]
#Name file
name=sys.argv[3]
fragment=sys.argv[4]
#Parameter necessary to code different type of samples
prefixMETA = sys.argv[5]
prefixBARC = sys.argv[6]
prefixCONTROL = sys.argv[7]
numberID = sys.argv[8]
duplicat =sys.argv[9]
duplicatformule=duplicat.replace("/","|")



###############
#### BEGIN ####
###############


###OPEN FILES
csvfile = open(sys.argv[1])
Occurencefile = pd.read_csv(csvfile, delimiter="\t", header=0)

colname = list(Occurencefile.columns.values)


first=colname[0:11] #KEEP the 11th first colum of the dataframe previously produce by FROGS. It contains information about taxonomic affiliations, identity percentages, sequences etc...
BARCODE = []
METABARCODE = []
CONTROL = []


### FILTERING SAMPLES ACCORDING TO THEIR ID_NAME, WHICH INDICATES WHETHER THEY ARE BARCODING, METABARCODING, CONTROL SAMPLES 

for i in colname[11:]:   # 11th first colum are already keep in "first" variable
    if re.search(prefixMETA+ r'[0-9]{'+numberID+'}-('+ duplicatformule + r')', i):
        METABARCODE.append(i)
    if re.search(prefixBARC+ r'[0-9]{'+numberID+'}-('+ duplicatformule + r')', i):
        BARCODE.append(i)
    if re.search(prefixCONTROL+r'(I|E|P)', i):
        CONTROL.append(i)

### PRODUCE SUBSUET OF THE INPUT FILE  :  

## METABARCODING
#Metabarcoding subfile : 
METABARCODE=first+METABARCODE+CONTROL
METAfile= Occurencefile[METABARCODE]
# Do the read sum to update "observation_sum" variable of the subfile :
sum_by_row = METAfile.iloc[:, 11:].sum(axis=1) 
METAfile.loc[:,"observation_sum"]=sum_by_row
METAfile= METAfile[METAfile["observation_sum"] != 0]
#Create metabarcoding file : 
if not os.path.exists(os.path.join(outputdir,"METABARCODING")):#CHECKS WHETHER THE "METABARCODING" DIRECTORY EXISTS IN THE FILE DIRECTORY TREE 
    os.makedirs(os.path.join(outputdir,"METABARCODING"))
METAfile.to_csv(os.path.join(outputdir,"METABARCODING/"+name+fragment+"_abundance_raw_data_METABARCODING.tsv"), index=False,sep="\t")


##BARCODING
#barcoding subfile : 
if len(BARCODE) != 0: #If barcoding data 

    BARCODE=first+BARCODE+CONTROL
    BARCOfile= Occurencefile[BARCODE]
    # Do the read sum to update "observation_sum" variable of the subfile :
    sum_by_row = BARCOfile.iloc[:, 11:].sum(axis=1)
    BARCOfile.loc[:,"observation_sum"]=sum_by_row
    BARCOfile = BARCOfile[BARCOfile["observation_sum"]!= 0]
    #Create barcoding file : 

    if not os.path.exists(os.path.join(outputdir,"BARCODE")): #CHECKS WHETHER THE "BARCODING" DIRECTORY EXISTS IN THE FILE DIRECTORY TREE 
        os.makedirs(os.path.join(outputdir,"BARCODE"))
    BARCOfile.to_csv(os.path.join(outputdir,"BARCODE/"+name+fragment+"_abundance_raw_data_BARCODING.tsv"), index=False,sep="\t")
else: # If no barcoding data 
    if not os.path.exists(os.path.join(outputdir,"BARCODE")): #CHECKS WHETHER THE "BARCODING" DIRECTORY EXISTS IN THE FILE DIRECTORY TREE 
        os.makedirs(os.path.join(outputdir,"BARCODE"))
    Null=os.path.join(outputdir,"BARCODE/"+name+fragment+"_abundance_raw_data_BARCODING.tsv")
    with open(Null, mode='w', newline='') as fichier_tsv:
        writer = csv.writer(fichier_tsv)
        writer.writerow("NO BARCODING")


### PRODUCE METADATA for metabarcoding subfile : 

with open(os.path.join(outputdir,"METABARCODING/"+name+fragment+"_Metadata_METABARCODING.csv"),'w',newline='') as csvfile:
    csvwriter=csv.writer(csvfile,delimiter=";")
    header = ["observation_name", "rep", "control", "biological_unit", "project"]
    csvwriter.writerow(header)

    for i in METABARCODE:
        if re.search(prefixMETA+ r'[0-9]{'+numberID+'}',i):
            if re.search(duplicat[0],i):
                row = [i, "1", "no", i[:-2], name]
                csvwriter.writerow(row)
            elif re.search(duplicat[-1],i):
                row = [i, "2", "no", i[:-2], name]
                csvwriter.writerow(row)
        elif re.search(r'^NC(I|E|P)', i) and  re.search(r'-{}$'.format(duplicat[0]), i):
            row= [i, "1", "negative", i[:-2], name]
            csvwriter.writerow(row)
        elif re.search(r'^NC(I|E|P)', i) and  re.search(r'-{}$'.format(duplicat[-1]), i):
            row= [i, "2", "negative", i[:-2], name]
            csvwriter.writerow(row)