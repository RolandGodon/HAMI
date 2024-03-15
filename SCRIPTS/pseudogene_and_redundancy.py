#!/usr/bin/env python
# -*-coding: utf-8 -*

"""
Author : Benoit Penel
Date : 11/24/2022
version : 1.0
PYTHON3


DESCRIPTION :
-------------------------------------------------------------------------------------------------------------------------------
1- Firstly, this script merges the output produced by the R script "filter_frogs.R". And lead to the production of a occurence matrice file as produce by FROGS but cleaned.

2- A pseudogene (NUMPTs) filtering step is then performed as recommended in the literature (see: SONG et al., 2018 - Proceedings of the National Academy of Sciences) on the occurence file.
It mean that affilaition associated to sequences whose reading frame has shifted -as a result of indel events - or in which stop codons have appeared are discarded. 
Data associated with pseudogenes is saved in a first file, while data exempt from pseudogenes is saved in a second file.


3- The occurrences file is processed to reduce data redundancy caused by multiple OTUs associated with the same species, reflecting intraspecific variability.
This complexity is irrelevant to us and complicates data interpretation. 
We therefore merge the occurrences of OTUs whose taxonomic affiliations are identical and whose match threshold with the reference barcode is similar (up to the arbitrary threshold or lower than the threshold) 


4- A second cleaning step is then carried out. Taxa represented twice, i.e. represented by OTUs with values below the threshold and above the reference match threshold, are re-examined. 
If the 2 OTUs of the same taxon are present in exactly the same samples, the group with the lowest match threshold with the reference barcode is deleted as it is considered to be noise or a residual pseudogene. 


The resultant occurrence file is then saved (final_files/XXXXX_Occurence_file.tsv) and reviewed by a parataxonomist in the context of the HAMI framework.

--------------------------------------------------------------------------------------------------------------------------------


"""

import sys
import csv
import os 
import pandas as pd
import numpy 


###################################################
############# 1- START PRE-TREATMENT ##############
###################################################

###### PARAM YOUR PROJECT :
abundpath=sys.argv[1] #path to abundance file
taxpath=sys.argv[2] #path to taxonomic file
prefix=sys.argv[3] #Prefixe commun au ID des échantillon de métabarcoding 
affiliation=sys.argv[4] #Target affiliation
directorypath=sys.argv[5] #path to the directory of the project
Name=sys.argv[6] #Name of the project
fragment=sys.argv[7] #Name of the DNA fragment
size=int(sys.argv[8]) #Size of the DNA fragment
reading_frame=int(sys.argv[9]) #Reading frame of the DNA sequence
threshold=int(sys.argv[10])
stop_codon=sys.argv[11:] # List of codon stop

#Summary file
Diclog={}
Diclog["id"]=[]
Diclog["Ntaxa"]=[]
Diclog["Nseq"]=[]


###OPEN FILES
csvfile=open(abundpath) #abundance file open
Abondance = pd.read_csv(csvfile,delimiter=";",header=0)


csvfile2=open(taxpath) #affiliation file open
Tax = pd.read_csv(csvfile2,delimiter=";",header=0)


### CALCUL OF THE TOTAL NUMBER OF READ PER OTU 
# or sums all the values for a given row, for all rows 
Abondance=Abondance.assign(observation_sum = Abondance.loc[:,Abondance.columns!="cluster"].sum(axis=1))



###MERGE THE TWO OPEN FILES WITH "CLUSTER" COLUMN 
df=pd.merge(Tax,Abondance,on=["cluster"]) 

###REINDEXATION OF THE COLUMN
headdf=["#comment","blast_taxonomy","blast_subject","blast_perc_identity","blast_perc_query_coverage","blast_evalue",
"blast_aln_length","seed_id","seed_sequence","cluster","observation_sum"] # PREPARE THE DESIRED ORDER

for col in df.columns: ##SEARCH FOR COLUMNS ASSOCIATED WITH SAMPLES OCCURENCE USING THEIR PREFIX NAME 
    if prefix in col:
        headdf.append(col)

df=df.reindex(columns=headdf)#REINDEXATION

### CHOOSE CERTAIN CLUSTERS ACCORDING TO THE AFFILIATION 
col=df[df['blast_taxonomy'].str.contains(affiliation,na=False) ] 


### THE OCCURRENCE FILE ASSOCIATED WITH THE TARGET TAXON IS CREATED
if not os.path.exists(os.path.join(directorypath+"/METABARCODING/intermediary_step")):#CHECKS WHETHER THE "intermediary_step" DIRECTORY EXISTS IN THE FILE DIRECTORY TREE 
    os.makedirs(os.path.join(directorypath+"/METABARCODING/intermediary_step"))
col.to_csv(os.path.join(directorypath+"/METABARCODING/intermediary_step/"+Name+fragment+'_abundance_'+affiliation+'_raw.csv'),index_label=False,index=False,sep=';')


#Set logfile of initial file 
Diclog["id"].append("After.filter.3")
Diclog["Nseq"].append(int(df.shape[0]))
Diclog["Ntaxa"].append(int(df["blast_taxonomy"].nunique()))


###############################################################
### 2- FILTRE PSEUDOGENE (Length and stop codon) ##############
###############################################################

seq=pd.Series(col.seed_sequence)


lenaccepted=[size-9,size-6,size-3,size,size+3,size+6,size+9] #ACCEPTABLE SEQUENCE SIZE VALUE 
accepted=[]
deleted=[]


# CHECK INDEL IN EACH SEQUENCES
for i in range(len(seq)):
   if len(seq.iloc[i]) not in lenaccepted: #REMOVE AFFILIATIONS ASSOCIATED WITH SEQUENCES WITH INDEL EVENTS THAT ARE NOT MULTIPLES OF 3
      deleted.append(i) #Save OTU with wrong sequence

   else:      #CHECK STOP CODON IN EACH SEQUENCES WITHOUT INDEL 
      catch = numpy.arange(reading_frame, len(seq.iloc[i]), 3) #CREATES CODONS OF THREE NUCLEOTIDES FROM THE SECOND NUCLEOTIDE (ARGUMENT "READING FRAME") IN THE SEQUENCE. PLEASE CHECK THE READING FRAME ACCORDING TO THE TARGETED DNA FRAGMENT
      bool=False
      for j in catch:
         codon = seq.iloc[i][j:j + 3]
         if codon in stop_codon:
            deleted.append(i) #Save OTU with wrong sequence
            bool=True
            break
      if bool==False:  
         accepted.append(i)  #Save OTU without erreur of reading frame, codon stop 


pseudogene=col.iloc[deleted]
realcol=col.iloc[accepted]

#Save an occurence file without pseudogenic OTU
realcol.to_csv(os.path.join(directorypath+"/METABARCODING/intermediary_step/"+Name+fragment+'_abundance_'+affiliation+'_pseudogene_filtred.csv'),index_label=False,index=False,sep=';')

#Save an occurence file with pseudogenic OTU
if not os.path.exists(os.path.join(directorypath+"/METABARCODING/intermediary_step/pseudogene")):#CHECKS WHETHER THE "pseudogene" DIRECTORY EXISTS IN THE FILE DIRECTORY TREE 
    os.makedirs(os.path.join(directorypath+"/METABARCODING/intermediary_step/pseudogene"))
pseudogene.to_csv(os.path.join(directorypath+"/METABARCODING/intermediary_step/pseudogene/"+Name+fragment+'_pseudogene_deleted.csv'),index_label=False,index=False,sep=';')

#Set logfile for first pseudogene filter 
Diclog["id"].append("After.filter.pseudogene1")
Diclog["Nseq"].append(int(realcol.shape[0]))
Diclog["Ntaxa"].append(int(realcol["blast_taxonomy"].nunique()))


#############################################
############# 3- MERGE DATA #################
#############################################



###################
##### Pre-treatment
###################

def pre_trait(data): 
   
    ###### function for separating data by percentage of identity (numeric versus multi-affiliation)
    identity=pd.Series(data.blast_perc_identity)
   
    multi=data[identity=="multi-identity"]
    multi=multi.reset_index(drop=True) #Reset indexage 
    multi.to_csv(os.path.join(directorypath+"/METABARCODING/intermediary_step/before_merge_multi.csv"),sep=";",index=False) #temporary file

    db=data[identity!="multi-identity"]
    db=db.reset_index(drop=True)# Reset indexage 
    db=db.astype({'blast_perc_identity':numpy.float64}) #ensures that identity percentages are read as numerical values
    db.to_csv(os.path.join(directorypath+"/METABARCODING/intermediary_step/before_merge_without_multi.csv"),sep=";",index=False) #temporary file

    return db

data=pre_trait(realcol)


###################
#### Merge ########
###################

print("Start of the automatic OTU merger...")

##################################################################
### Automated merging for data without multi affiliation

DICO={} # Initialize an empty dictionary to store data
OTU=0  # Initialize OTU counter
donelist=[] # Initialize a list to keep track of processed taxa
sites=data.columns
sites=sites[11:]


print("Processing....")
for cell in range(len(data["blast_taxonomy"])): # For each OTU identified 

   if data["observation_sum"][cell]>0: # If observation count is greater than 0

      if data["blast_taxonomy"][cell] not in donelist:  # and the taxonomic name assigned is encountered for the first time
      
         #### Create specific lists for the new encoutered OTU
         DICO[OTU]={}
         DICO[OTU]["Taxon_name"]=[]
         DICO[OTU]["Nbr_cluster"]=[]
         DICO[OTU]["Cluster"]=[]
         DICO[OTU]["multi"]=[]
         DICO[OTU]["max_identity"]=[]
         DICO[OTU]["min_identity"]=[]
         DICO[OTU]["Observation"]=[]
         DICO[OTU]["Seq"]=[] 
         for i in sites:
            DICO[OTU][i]=[]

         ### and fill the lists with the following information               
         DICO[OTU]["Taxon_name"].append(data["blast_taxonomy"][cell])
         DICO[OTU]["Nbr_cluster"].append(int(1))
         if "Multi-affiliation" in data["blast_taxonomy"][cell]:
            DICO[OTU]["multi"].append(int(1))
         else :
            DICO[OTU]["multi"].append(int(0))
         DICO[OTU]["Observation"].append(int(data["observation_sum"][cell]))
         DICO[OTU]["Seq"].append(data["seed_sequence"][cell])
         if (data["blast_perc_identity"][cell]>0 and data["blast_perc_identity"][cell]<=100):
            DICO[OTU]["Cluster"].append(data["cluster"][cell])
            DICO[OTU]["max_identity"].append(data["blast_perc_identity"][cell])
            DICO[OTU]["min_identity"].append(data["blast_perc_identity"][cell])
         ### and saved OTU read occurences in each samples 
         for i in sites: 
            if pd.isnull(data[i][cell]) == True :
               DICO[OTU][i].append(int(0))
            elif data[i][cell] >=0:
               DICO[OTU][i].append(int(data[i][cell]))

         donelist.append( data["blast_taxonomy"][cell]) #Save the name of the OTU
         OTU=OTU+1 #anticipate the discovery of a new identified OTU to create a new list if necessary 
         
      

      else: # If taxon name of a given OTU has been encountered before 
         liste=0
         while liste < len(DICO): #looking for the specific list associated to the taxonomic name
            if data["blast_taxonomy"][cell] == DICO[liste]["Taxon_name"][0]: # If the list it found
            
               ## OTU with pourcentage of identity > Threshold % :
               # IF both OTU with similar taxonomic name also have similar pourcentage of identity (here > to the threshold)
               if  data["blast_perc_identity"][cell] >= int(threshold): 
                  if DICO[liste]["min_identity"][0]>= int(threshold): 
                     #merge of the OTU information 
                     DICO[liste]["Nbr_cluster"][0]=DICO[liste]["Nbr_cluster"][0]+int(1) #Number of cluster merge
                     #Merge of occurences 
                     DICO[liste]["Observation"][0]=int(DICO[liste]["Observation"][0]+data["observation_sum"][cell]) 
                     for i in sites:
                        if pd.isnull(data[i][cell]) == False and data[i][cell] >=1:
                           DICO[liste][i][0]=int(DICO[liste][i][0]+data[i][cell])

                     #Refreshes the identity percentage range 
                     if data["blast_perc_identity"][cell]<DICO[liste]["min_identity"][0]:
                        DICO[liste]["min_identity"][0]=data["blast_perc_identity"][cell]
                     if data["blast_perc_identity"][cell]>DICO[liste]["max_identity"][0]:
                        DICO[liste]["max_identity"][0]=data["blast_perc_identity"][cell]
                     #Add cluster name and sequence
                     DICO[liste]["Cluster"].append(data["cluster"][cell])
                     DICO[liste]["Seq"].append(data["seed_sequence"][cell])
                     break;
                  
               ## OTU with pourcentage of identity < Threshold % :
               # IF both OTU with similar taxonomic name also have similar pourcentage of identity (here < to the threshold)
               if  data["blast_perc_identity"][cell] <= int(threshold): 
                  if DICO[liste]["max_identity"][0]<=int(threshold): 
                     #merge of the OTU information 
                     DICO[liste]["Nbr_cluster"][0]=DICO[liste]["Nbr_cluster"][0]+int(1) 
                     #merge of occurences 
                     DICO[liste]["Observation"][0]=int(DICO[liste]["Observation"][0]+data["observation_sum"][cell])
                     for i in sites:
                           if pd.isnull(data[i][cell]) == False and data[i][cell] >=1:
                              DICO[liste][i][0]=int(DICO[liste][i][0]+data[i][cell])

                     #Refreshes the identity percentage range 
                     if data["blast_perc_identity"][cell]<DICO[liste]["min_identity"][0]:
                        DICO[liste]["min_identity"][0]=data["blast_perc_identity"][cell]
                     if data["blast_perc_identity"][cell]>DICO[liste]["max_identity"][0]:
                        DICO[liste]["max_identity"][0]=data["blast_perc_identity"][cell]
                     #Add cluster name and sequence
                     DICO[liste]["Cluster"].append(data["cluster"][cell])
                     DICO[liste]["Seq"].append(data["seed_sequence"][cell])
                     break;

            liste=liste+1
            
         ## If we reach the limits of the loop, it's because the OTU encountered, 
         #although it has the same taxonomic affiliation as a previously one, it has a different percentage of identity 
         #We thus create another group. Next time, if a new OTU with a similar taxonomic affiliation is found,
         # it will be in one of the two previous case (>threshold or <threshold)

         if liste==len(DICO):
               ###CREATES THE LISTS 
               DICO[OTU]={}
               DICO[OTU]["Taxon_name"]=[]
               DICO[OTU]["Nbr_cluster"]=[]
               DICO[OTU]["Cluster"]=[]
               DICO[OTU]["multi"]=[]
               DICO[OTU]["max_identity"]=[]
               DICO[OTU]["min_identity"]=[]
               DICO[OTU]["Observation"]=[]
               DICO[OTU]["Seq"]=[]
               for i in sites:
                  DICO[OTU][i]=[]

               ### FILLED THE LISTS
               DICO[OTU]["Taxon_name"].append(data["blast_taxonomy"][cell])
               DICO[OTU]["Nbr_cluster"].append(int(1))
               if "Multi-affiliation" in data["blast_taxonomy"][cell]:
                  DICO[OTU]["multi"].append(int(1))
               else :
                  DICO[OTU]["multi"].append(int(0))
               DICO[OTU]["Cluster"].append(data["cluster"][cell])
               DICO[OTU]["max_identity"].append(data["blast_perc_identity"][cell])
               DICO[OTU]["min_identity"].append(data["blast_perc_identity"][cell])
               DICO[OTU]["Observation"].append(int(data["observation_sum"][cell]))
               DICO[OTU]["Seq"].append(data["seed_sequence"][cell])
               ## saved occurences 
               for i in sites:
                  if pd.isnull(data[i][cell]) == True :
                     DICO[OTU][i].append(int(0))
                  elif data[i][cell] >=0:
                     DICO[OTU][i].append(int(data[i][cell]))

               OTU=OTU+1

#################################################
###  Automated merging for multi-identity data

#similar process is done 

print("Processing....")

donelist=[]
DICO1={}
OTU=0


csvfile=open(os.path.join(directorypath+"/METABARCODING/intermediary_step/before_merge_multi.csv"))
multidata = pd.read_csv(csvfile,delimiter=";",header=0,index_col=None)


for cell in range(len(multidata["blast_taxonomy"])): 
   if multidata["observation_sum"][cell]>0:
      if multidata["blast_taxonomy"][cell] not in donelist : 
      

         DICO1[OTU]={}
         DICO1[OTU]["Taxon_name"]=[]
         DICO1[OTU]["Nbr_cluster"]=[]
         DICO1[OTU]["Cluster"]=[]
         DICO1[OTU]["multi"]=[]
         DICO1[OTU]["max_identity"]=[]
         DICO1[OTU]["min_identity"]=[]
         DICO1[OTU]["Observation"]=[]
         DICO1[OTU]["Seq"]=[]
         for i in sites:
            DICO1[OTU][i]=[]
              
         DICO1[OTU]["Taxon_name"].append(multidata["blast_taxonomy"][cell])
         DICO1[OTU]["Nbr_cluster"].append(int(1))
         DICO1[OTU]["multi"].append(int(1))
         DICO1[OTU]["Observation"].append(int(multidata["observation_sum"][cell]))
         DICO1[OTU]["Cluster"].append(multidata["cluster"][cell])
         DICO1[OTU]["max_identity"].append("multi_identity")
         DICO1[OTU]["min_identity"].append("multi_identity")
         DICO1[OTU]["Seq"].append(multidata["seed_sequence"][cell])
         for i in sites:
            if pd.isnull(multidata[i][cell]) == True :
               DICO1[OTU][i].append(int(0))
            elif data[i][cell] >=0:
               DICO1[OTU][i].append(int(multidata[i][cell]))

         OTU=OTU+1
         donelist.append( multidata["blast_taxonomy"][cell])

      else: 
         liste=0
         while liste < len(DICO1):
            if multidata["blast_taxonomy"][cell] == DICO1[liste]["Taxon_name"][0]:         
                        DICO1[liste]["Nbr_cluster"][0]=DICO1[liste]["Nbr_cluster"][0]+int(1) 
                        DICO1[liste]["Observation"][0]=int(DICO1[liste]["Observation"][0]+multidata["observation_sum"][cell])
                        for i in sites:
                           if pd.isnull(multidata[i][cell]) == False:
                              DICO1[liste][i][0]=int(DICO1[liste][i][0]+multidata[i][cell])

                        DICO1[liste]["Cluster"].append(multidata["cluster"][cell])  
                        DICO1[liste]["Seq"].append(multidata["seed_sequence"][cell]) 
                        break;
                  
                     
            liste=liste+1


#############################
### 2nd pseudogenes check 
## Remove duplicates merged OTU < threshold  when they have the same occurrences as merged OTU >threshold



print("Processing....")

donelist=[]
pose=[]
it=0
for i in DICO: ##LIST ALL DUPLICATE TAXA
    if DICO[i]["Taxon_name"][0] not in donelist:
        donelist.append(DICO[i]["Taxon_name"][0])
        pose.append(it)
    it+=1



it=0
nbrsite=0
delposition=[]

for i in DICO:
    for j in range(len(donelist)): 
        if DICO[i]["Taxon_name"][0] == donelist[j] and it!=pose[j] :
            for k in DICO[0]:
                if prefix in k: #CHECK REDONDANCE IN EACH SAMPLES IN DUPLICATE TAXA
                     if DICO[i][k][0] == 0 and DICO[pose[j]][k][0] == 0:
                        nbrsite+=1
                     elif DICO[i][k][0] > 0 and DICO[pose[j]][k][0] > 0:
                        nbrsite+=1
                     else:
                        break
            if nbrsite==len(DICO[0])-8: # IF FULLY REDONDANT, CHOOSE THE TAXA WITH THE MAXIMUM OF IDENTITY (THE SECOND ONE IS PSEUDOGENE OR NOISE )
               if DICO[i]["max_identity"][0] >  DICO[pose[j]]["max_identity"][0]  :      
                  delposition.append(pose[j])
               if DICO[i]["max_identity"][0] <  DICO[pose[j]]["max_identity"][0]  : 
                  delposition.append(i)    
            nbrsite=0
    it+=1




#####################################################
### Compiling the CSV  ##############################
#####################################################

nseq=0
listtaxa=[]
Ntaxa=0

print('Automation complete ')
print('compilation of files')

if not os.path.exists(os.path.join(directorypath+"/METABARCODING/final_files")): #CHECKS WHETHER THE "final_files" DIRECTORY EXISTS IN THE FILE DIRECTORY TREE 
   os.makedirs(os.path.join(directorypath+"/METABARCODING/final_files"))

with open(os.path.join(directorypath+"/METABARCODING/final_files/"+Name+fragment+'_'+affiliation+'_final_abundance_file.tsv'), 'w') as f:
   i=0
   for key in DICO[0]:
      if key != 'Cluster' and key != 'Seq':
         [f.write(key)]
         [f.write('\t')]
      i=i+1
   i=0
   [f.write('\n')]

   for i in range(len(DICO)):
      if i not in delposition:
         #Calcul number of Taxa AND Seq saved in the final files
         nseq=nseq+int(DICO[i]["Nbr_cluster"][0])
         if DICO[i]["Taxon_name"][0] not in listtaxa:
            listtaxa.append(DICO[i]["Taxon_name"][0])
            Ntaxa=Ntaxa+1
         #Write file
         [f.write(DICO[i]["Taxon_name"][0] )]
         [f.write('\t')]
         [f.write(str(DICO[i]["Nbr_cluster"][0]))]
         [f.write('\t')]
         [f.write(str(DICO[i]["multi"][0]))]
         [f.write('\t')]
         if len(DICO[i]["min_identity"]) >=1 : 
            [f.write(str(DICO[i]["max_identity"][0] ))]
            [f.write('\t')]
            [f.write(str(DICO[i]["min_identity"][0] ))]
            [f.write('\t')]
         else:
            [f.write('')]
            [f.write('\t')]
            [f.write('')]
            [f.write('\t')]
         [f.write(str(DICO[i]["Observation"][0]))]       
         for j in sites:
            [f.write('\t')]
            [f.write(str(DICO[i][j][0] ))]
         [f.write('\n')]

   for i in range(len(DICO1)):
      #Calcul number of Taxa AND Seq saved in the final files
      nseq=nseq+int(DICO1[i]["Nbr_cluster"][0])
      if DICO1[i]["Taxon_name"][0] not in listtaxa:
         listtaxa.append(DICO1[i]["Taxon_name"][0])
         Ntaxa=Ntaxa+1
      #Write file
      [f.write(DICO1[i]["Taxon_name"][0] )]
      [f.write('\t')]
      [f.write(str(DICO1[i]["Nbr_cluster"][0]))]
      [f.write('\t')]
      [f.write(str(DICO1[i]["multi"][0]))]
      [f.write('\t')]
      if len(DICO1[i]["min_identity"]) >=1 : 
        [f.write(str(DICO1[i]["max_identity"][0] ))]
        [f.write('\t')]
        [f.write(str(DICO1[i]["min_identity"][0] ))]
        [f.write('\t')]
      else:
         [f.write('')]
         [f.write('\t')]
         [f.write('')]
         [f.write('\t')]
      [f.write(str(DICO1[i]["Observation"][0]))]       
      for j in sites:
         [f.write('\t')]
         [f.write(str(DICO1[i][j][0] ))]
      [f.write('\n')]
f.close()

with open(os.path.join(directorypath+"/METABARCODING/final_files/"+Name+fragment+'_'+affiliation+'_cluster_merged.txt'), 'w') as f:
   for key in DICO:
      if key not in delposition :
         [f.write(">"+ DICO[key]["Taxon_name"][0]+ "\n")]
         if DICO[key]["Cluster"]:
            for x in DICO[key]["Cluster"]:
               if pd.isna(x)==False:
                  [f.write(x+"\n")]
            [f.write('\n')]

   for key in DICO1:
      [f.write(">"+ DICO1[key]["Taxon_name"][0]+ "\n")]
      if DICO1[key]["Cluster"]:
         for x in DICO1[key]["Cluster"]:
            if pd.isna(x)==False:
               [f.write(x+"\n")]
         [f.write('\n')]
f.close()

with open(os.path.join(directorypath+"/METABARCODING/intermediary_step/pseudogene/"+Name+fragment+'_pseudogene_deleted2.tsv'), 'w') as f:
   i=0
   for key in DICO[0]:
      if key != 'Cluster' and key != 'Seq':
         [f.write(key)]
         [f.write('\t')]
      i=i+1

   i=0
   [f.write('\n')]

   for i in range(len(DICO)):
      if i in delposition:
         [f.write(DICO[i]["Taxon_name"][0] )]
         [f.write('\t')]
         [f.write(str(DICO[i]["Nbr_cluster"][0]))]
         [f.write('\t')]
         [f.write(str(DICO[i]["multi"][0]))]
         [f.write('\t')]
         if len(DICO[i]["min_identity"]) >=1 : 
            [f.write(str(DICO[i]["max_identity"][0] ))]
            [f.write('\t')]
            [f.write(str(DICO[i]["min_identity"][0] ))]
            [f.write('\t')]
         else:
            [f.write('')]
            [f.write('\t')]
            [f.write('')]
            [f.write('\t')]
         [f.write(str(DICO[i]["Observation"][0]))]       
         for j in sites:
            [f.write('\t')]
            [f.write(str(DICO[i][j][0] ))]
         [f.write('\n')]
f.close()


with open(os.path.join(directorypath+"/METABARCODING/final_files/"+Name+fragment+'_'+affiliation+'_cluster_seq_keep.txt'), 'w') as f:
   for key in DICO:
      if key not in delposition :
         if DICO[key]["Cluster"]:
            it=0
            for x in DICO[key]["Cluster"]:
               if pd.isna(x)==False:
                  [f.write(">"+ DICO[key]["Taxon_name"][0]+ "_" + x + "\n")]
                  [f.write(DICO[key]["Seq"][it]+"\n")]
                  it+=1
               [f.write('\n')]

   for key in DICO1:
         if DICO1[key]["Cluster"]:
            it=0
            for x in DICO1[key]["Cluster"]:
               if pd.isna(x)==False:
                  [f.write(">"+ DICO1[key]["Taxon_name"][0]+ "_" + x + "\n")]
                  [f.write(DICO1[key]["Seq"][it]+"\n")]
                  it+=1
               [f.write('\n')]
f.close()


#Set logfile for second pseudogene filter 
Diclog["id"].append("After.filter.pseudogene2")
Diclog["Nseq"].append(int(nseq))
Diclog["Ntaxa"].append(int(Ntaxa))

#Write log file
with open(os.path.join(directorypath+"/METABARCODING/final_files/"+Name+fragment+'_'+affiliation+'_summary.txt'), 'w') as f:
   col = list(Diclog.keys())
   val = list(Diclog.values())
   f.write("\t".join(col) + "\n")
   for line in zip(*val):
      f.write("\t".join(map(str, line)) + "\n")
f.close()

print("Programme completed. Have a nice day :) ")




################
###### END #####
################


