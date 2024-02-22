#!/usr/bin/env python
# -*-coding: utf-8 -*

"""
Author : Benoit Penel
Date : 11/24/2022
version : 1.0
PYTHON3


DESCRIPTION :
-------------------------------------------------------------------------------------------------------------------------------
Ce script permet de merger les sorti produit par le script R "filter_frogs.R". 

Une fois fusionner, un pré-nettoyage de données est réalisé. L'ensemble des clusters associées a des séquences de tailles non multiples de trois et présentant des codon stop
sont écartés du jeux de donnée.

S'en suit une phase de simplification des données d'abondances à l'image de ce que fait le script Merge_raw_data2.py avec les données AgriB02.

Une seconde étape de nettoyage est ensuite réalisé. Les taxon représenté deux fois, c'est a dire ayant un groupe avec de identity faible et un second avec des identité >97,
sont ne nouveau ré-examiner. Si ce doublon de taxon possèdent des occurences dans exactement les mêmes échantillons, le groupe ayant l'identité la plus faible est supprimer car considéré comme pseudogene.

Finalement, une phase dernière phase création de fichier qui constituent la base de donnée à révisé par Axel est lancé.

--------------------------------------------------------------------------------------------------------------------------------

RUN COMMANDE !
--------------------------------------------------------------------------------------------------------------------------------
      Python3 script/Merge_after_filtre.py  Abundance.filter3_file  Tax.filter3_file
--------------------------------------------------------------------------------------------------------------------------------

OUTPUT : 
--------------------------------------------------------------------------------------------------------------------------------
- AgriB03_abundance_coleoptera_raw_META_filtred.csv            Concatenation des deux fichier issue du script R et extraction des données Coléoptère
- AgriB03_abundance_coleoptere_pseudogene_filtred.csv          Homologue au fichier précédent, mais avec suppression des cluster associée à des pseudogene
- AgriB03_pseudogene_deleted.csv                               Fichier contenant les clusters pseudogene supprimé dans le fichier précédent
- AGRIBO3_AUTO_DATA.csv                                        Fichier homologue a AUTO_DATA de AGRIB02 a révisé par axel
- AGRIB03_Cluster_merged.txt                                   Fichier annexe a AUTO_DATA expliquants quel clusters à été fusionné
- AGRIB03_pseudogene_deleted2.csv
-------------------------------------------------------------------------------------------------------------------------------




"""

import sys
import csv
import os 
import pandas as pd
import numpy 



nameprojet="AgriB04"
outputdir="/home/penelben/Documents/THESE/AGRIBIODIV/AGRIB04/DATA/METABARCODE/"
prefixsample="CMEY"
phylum="Coleoptera"


###################################################
############# DEBUT PRETRAITEMENT #################
###################################################

###OPEN FILES
csvfile=open(sys.argv[1])
Abondance = pd.read_csv(csvfile,delimiter=";",header=0)


csvfile2=open(sys.argv[2])
Tax = pd.read_csv(csvfile2,delimiter=";",header=0)


### CALCUL OF THE TOTAL NUMBER OF READ PER CLUSTER
Abondance=Abondance.assign(observation_sum = Abondance.loc[:,Abondance.columns!="cluster"].sum(axis=1))


### MERGE AND REINDEXATION OF THE COLUMNS

df=pd.merge(Tax,Abondance,on=["cluster"])

headdf=["#comment","blast_taxonomy","blast_subject","blast_perc_identity","blast_perc_query_coverage","blast_evalue",
"blast_aln_length","seed_id","seed_sequence","cluster","observation_sum"]


for col in df.columns:
    if prefixsample in col:
        headdf.append(col)

df=df.reindex(columns=headdf)

### CHOOSE OF COLEOPTERA DATA AND THE 4 LEVELS OF TAXONOMY 
##### ICI A ENLEVER POUR AUTRE PROJET

col=df[df['blast_taxonomy'].str.contains(phylum,na=False) ] ### Selectionne seulement les lignes avec coleoptera dedans


Taxo=col['blast_taxonomy'].str.split(';',expand=True)
Taxo=Taxo.iloc[:,[3,4,5,6]]
Taxo.iloc[:,0]=Taxo.iloc[:,0]+','+Taxo.iloc[:,1]+','+Taxo.iloc[:,2]+','+Taxo.iloc[:,3]

col["blast_taxonomy"]=Taxo.iloc[:,0]


### CREATE CSV FILE

col.to_csv(os.path.join(outputdir+nameprojet+'_abundance_coleoptera_raw_META_filtred.csv'),index_label=False,index=False,sep=';')





###############################################################
### FILTRE PSEUDOGENE (Length and stop codon)##################
###############################################################

seq=pd.Series(col.seed_sequence)


lenaccepted=[400,403,406,409,412,415,418,421,424,427,430]
accepted=[]
deleted=[]


# CHECK INDEL IN EACH SEQUENCES
for i in range(len(seq)):
    if len(seq.iloc[i]) not in lenaccepted: #DELETE LINES WITH INDEL
        deleted.append(i)
    else:

        #CHECK STOP CODON IN EACH SEQUENCES WITHOUT INDEL 
        catch = numpy.arange(1, len(seq.iloc[i]), 3)
        stopCodonPositions = []
        bool=False
        for j in catch:
            codon = seq.iloc[i][j:j + 3]
            if codon == 'TAA' or codon == 'TAG' : # DELETE LINES WITH STOP CODON
                deleted.append(i)
                bool=True
                break

        if bool==False:  
            accepted.append(i)           


pseudogene=col.iloc[deleted]
realcol=col.iloc[accepted]
realcol.to_csv(os.path.join(outputdir+nameprojet+'_abundance_coleoptere_pseudogene_filtred.csv'),index_label=False,index=False,sep=';')
pseudogene.to_csv(os.path.join(outputdir+nameprojet+'_pseudogene_deleted.csv'),index_label=False,index=False,sep=';')



#############################################
############# MERGE DATA ####################
#############################################



###################
##### Pretraitement

def pre_trait(data): 
   
    ###### Séparation des data en fonction identity
    identity=pd.Series(data.blast_perc_identity)
   
    multi=data[identity=="multi-identity"]
    multi=multi.reset_index(drop=True)
    multi.to_csv(os.path.join(outputdir+nameprojet+'_before_merge_multi.csv'),sep=";",index=False)

    db=data[identity!="multi-identity"]
    db=db.reset_index(drop=True)
    db=db.astype({'blast_perc_identity':numpy.float64})
    db.to_csv(os.path.join(outputdir+nameprojet+'_before_merge_without_multi.csv'),sep=";",index=False)
   
   

    return db

data=pre_trait(realcol)


###################
#### Merge ########
###################

print("Début de l'automatisation....")

##################################################################
### Automatisation de la fusion pour les  données  sans les multi

DICO={}
OTU=0
donelist=[]
sites=data.columns
sites=sites[11:]


print("Processing....")
for cell in range(len(data["blast_taxonomy"])): # Pour chaque nom de taxa

   if data["observation_sum"][cell]>0:

      if data["blast_taxonomy"][cell] not in donelist: #Si nom jamais rencontré
      
         ###CRÉE LES LISTES
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

         ### REMPLI LES LISTES                
         DICO[OTU]["Taxon_name"].append(data["blast_taxonomy"][cell])
         DICO[OTU]["Nbr_cluster"].append(int(1))
         DICO[OTU]["multi"].append(int(0))
         DICO[OTU]["Observation"].append(int(data["observation_sum"][cell]))
         DICO[OTU]["Seq"].append(data["seed_sequence"][cell])
         if (data["blast_perc_identity"][cell]>0 and data["blast_perc_identity"][cell]<=100):
            DICO[OTU]["Cluster"].append(data["cluster"][cell])
            DICO[OTU]["max_identity"].append(data["blast_perc_identity"][cell])
            DICO[OTU]["min_identity"].append(data["blast_perc_identity"][cell])
         for i in sites:
            if pd.isnull(data[i][cell]) == True :
               DICO[OTU][i].append(int(0))
            else:
               if data[i][cell] >=0:
                  DICO[OTU][i].append(int(data[i][cell]))
               else :
                  DICO[OTU][i].append(int(data[i][cell]))
                  DICO[OTU]["multi"][0]=int(1)

         

         OTU=OTU+1
         donelist.append( data["blast_taxonomy"][cell])
      

      else: # Si déjà rencontré
         liste=0
         while liste < len(DICO):
            if data["blast_taxonomy"][cell] == DICO[liste]["Taxon_name"][0]: # Je rentre dans le taxon concerner
            
               ## New cluster avec > 97 % :
               if  data["blast_perc_identity"][cell] >=97: # Si le nouveau cluster a fusionné à une identity >97

                  if len(DICO[liste]["min_identity"])>0:
                     ## Groupe avec > 97 %:
                     if DICO[liste]["min_identity"][0]>=97: # Si le groupe de cluster est aussi a >97
                        DICO[liste]["Nbr_cluster"][0]=DICO[liste]["Nbr_cluster"][0]+int(1) #ajout 1 au nombre de cluster fusionner
                        DICO[liste]["Observation"]=int(DICO[liste]["Observation"]+data["observation_sum"][cell])
                        for i in sites:
                           if pd.isnull(data[i][cell]) == False:

                              if data[i][cell] >=1 and DICO[liste][i][0] >=0 : 
                                 DICO[liste][i][0]=int(DICO[liste][i][0]+data[i][cell])

                              elif data[i][cell] < 0 and DICO[liste][i][0] <=0:
                                 DICO[liste][i][0]=int(DICO[liste][i][0]+data[i][cell])
                                 DICO[liste]["multi"][0]=int(1)

                              elif data[i][cell] < 0 and DICO[liste][i][0] >=1 :
                                 DICO[liste][i][0] =int(data[i][cell] - DICO[liste][i][0])
                                 DICO[liste]["multi"][0]=int(1)
                              
                              elif data[i][cell] >=1 and DICO[liste][i][0] <0 :
                                 DICO[liste][i][0] =int(DICO[liste][i][0]- data[i][cell])
                                 DICO[liste]["multi"][0]=int(1)

                              

                        #Réactulise les reference 
                        if data["blast_perc_identity"][cell]<DICO[liste]["min_identity"][0]:
                           DICO[liste]["min_identity"][0]=data["blast_perc_identity"][cell]
                        if data["blast_perc_identity"][cell]>DICO[liste]["max_identity"][0]:
                           DICO[liste]["max_identity"][0]=data["blast_perc_identity"][cell]
                        #Ajout nom cluster
                        DICO[liste]["Cluster"].append(data["cluster"][cell])
                        DICO[liste]["Seq"].append(data["seed_sequence"][cell])
                        break;
                  else:
                     ## Groupe avec identity vide :
                     if not(DICO[liste]["min_identity"]): #Si le groupe de cluster à été rentré a la main et n'a pas de value identity
                        DICO[liste]["Nbr_cluster"][0]=DICO[liste]["Nbr_cluster"][0]+int(1) #ajout 1 au nombre de cluster fusionner
                        DICO[liste]["Observation"]=int(DICO[liste]["Observation"]+data["observation_sum"][cell])
                        for i in sites:
                           if pd.isnull(data[i][cell]) == False:

                              if data[i][cell] >=1 and DICO[liste][i][0] >=0 : 
                                 DICO[liste][i][0]=int(DICO[liste][i][0]+data[i][cell])

                              elif data[i][cell] < 0 and DICO[liste][i][0] <=0:
                                 DICO[liste][i][0]=int(DICO[liste][i][0]+data[i][cell])
                                 DICO[liste]["multi"][0]=int(1)

                              elif data[i][cell] < 0 and DICO[liste][i][0] >=1 :
                                 DICO[liste][i][0] =int(data[i][cell] - DICO[liste][i][0])
                                 DICO[liste]["multi"][0]=int(1)
                              
                              elif data[i][cell] >=1 and DICO[liste][i][0] <0 :
                                 DICO[liste][i][0] =int(DICO[liste][i][0]- data[i][cell])
                                 DICO[liste]["multi"][0]=int(1)
                  
                        #Réactulise les reference 
                        DICO[liste]["min_identity"].append(data["blast_perc_identity"][cell])
                        DICO[liste]["max_identity"].append(data["blast_perc_identity"][cell])
                        #Ajout nom cluster
                        DICO[liste]["Cluster"].append(data["cluster"][cell])
                        DICO[liste]["Seq"].append(data["seed_sequence"][cell])
                        break;

               ## Cluster avec < 97 % :
               if  data["blast_perc_identity"][cell] <=97: # Si le groupe concerné à une identity <97
                  if DICO[liste]["max_identity"][0]<=97: # Si le nouveau cluster est aussi a <97
                     DICO[liste]["Nbr_cluster"][0]=DICO[liste]["Nbr_cluster"][0]+int(1) #ajout 1 au nombre de cluster fusionner
                     DICO[liste]["Observation"]=int(DICO[liste]["Observation"]+data["observation_sum"][cell])
                     for i in sites:
                           if pd.isnull(data[i][cell]) == False:

                              if data[i][cell] >=1 and DICO[liste][i][0] >=0 : 
                                 DICO[liste][i][0]=int(DICO[liste][i][0]+data[i][cell])

                              elif data[i][cell] < 0 and DICO[liste][i][0] <=0:
                                 DICO[liste][i][0]=int(DICO[liste][i][0]+data[i][cell])
                                 DICO[liste]["multi"][0]=int(1)

                              elif data[i][cell] < 0 and DICO[liste][i][0] >=1 :
                                 DICO[liste][i][0] =int(data[i][cell] - DICO[liste][i][0])
                                 DICO[liste]["multi"][0]=int(1)
                              
                              elif data[i][cell] >=1 and DICO[liste][i][0] <0 :
                                 DICO[liste][i][0] =int(DICO[liste][i][0]- data[i][cell])
                                 DICO[liste]["multi"][0]=int(1)

                     #Réactulise les reference 
                     if data["blast_perc_identity"][cell]<DICO[liste]["min_identity"][0]:
                        DICO[liste]["min_identity"][0]=data["blast_perc_identity"][cell]
                     if data["blast_perc_identity"][cell]>DICO[liste]["max_identity"][0]:
                        DICO[liste]["max_identity"][0]=data["blast_perc_identity"][cell]
                     #Ajout nom cluster
                     DICO[liste]["Cluster"].append(data["cluster"][cell])
                     DICO[liste]["Seq"].append(data["seed_sequence"][cell])
                     break;

               ## New cluster sans value dans identity :
               if pd.isnull(data["blast_perc_identity"][cell]) == True :
                  if len(DICO[liste]["min_identity"])>0:
                     ## Groupe avec > 97 %:
                     if DICO[liste]["min_identity"][0]>=97: # Si le groupe de cluster est a >97
                        DICO[liste]["Nbr_cluster"][0]=DICO[liste]["Nbr_cluster"][0]+int(1) #ajout 1 au nombre de cluster fusionner
                        DICO[liste]["Observation"]=int(DICO[liste]["Observation"]+data["observation_sum"][cell])
                        for i in sites:
                           if pd.isnull(data[i][cell]) == False:

                              if data[i][cell] >=1 and DICO[liste][i][0] >=0 : 
                                 DICO[liste][i][0]=int(DICO[liste][i][0]+data[i][cell])

                              elif data[i][cell] < 0 and DICO[liste][i][0] <=0:
                                 DICO[liste][i][0]=int(DICO[liste][i][0]+data[i][cell])
                                 DICO[liste]["multi"][0]=int(1)

                              elif data[i][cell] < 0 and DICO[liste][i][0] >=1 :
                                 DICO[liste][i][0] =int(data[i][cell] - DICO[liste][i][0])
                                 DICO[liste]["multi"][0]=int(1)
                              
                              elif data[i][cell] >=1 and DICO[liste][i][0] <0 :
                                 DICO[liste][i][0] =int(DICO[liste][i][0]- data[i][cell])
                                 DICO[liste]["multi"][0]=int(1)

                        #Ajout nom cluster
                        DICO[liste]["Cluster"].append(data["cluster"][cell])
                        DICO[liste]["Seq"].append(data["seed_sequence"][cell])
                        break;

                  ## Groupe sans value d'identity
                  if not( DICO[liste]["min_identity"]): #Si le groupe de cluster à été rentré a la main et n'a pas de value identity
                     DICO[liste]["Nbr_cluster"][0]=DICO[liste]["Nbr_cluster"][0]+int(1) #ajout 1 au nombre de cluster fusionner
                     DICO[liste]["Observation"]=int(DICO[liste]["Observation"]+data["observation_sum"][cell])
                     for i in sites:
                           if pd.isnull(data[i][cell]) == False:

                              if data[i][cell] >=1 and DICO[liste][i][0] >=0 : 
                                 DICO[liste][i][0]=int(DICO[liste][i][0]+data[i][cell])

                              elif data[i][cell] < 0 and DICO[liste][i][0] <=0:
                                 DICO[liste][i][0]=int(DICO[liste][i][0]+data[i][cell])
                                 DICO[liste]["multi"][0]=int(1)

                              elif data[i][cell] < 0 and DICO[liste][i][0] >=1 :
                                 DICO[liste][i][0] =int(data[i][cell] - DICO[liste][i][0])
                                 DICO[liste]["multi"][0]=int(1)
                              
                              elif data[i][cell] >=1 and DICO[liste][i][0] <0 :
                                 DICO[liste][i][0] =int(DICO[liste][i][0]- data[i][cell])
                                 DICO[liste]["multi"][0]=int(1)

                     #Ajout nom cluster
                     DICO[liste]["Cluster"].append(data["cluster"][cell])
                     DICO[liste]["Seq"].append(data["seed_sequence"][cell])
                     break;

            liste=liste+1
            
         ## Si j'atteint les limites de ma boucle c'est que je n'ai pas rencontré le car de figure qui m'interesse et je recrée des groupes
         if liste==len(DICO):
               ###CRÉE LA LISTES QUI ME MANQUE
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

               ### REMPLI LES LISTES                
               DICO[OTU]["Taxon_name"].append(data["blast_taxonomy"][cell])
               DICO[OTU]["Nbr_cluster"].append(int(1))
               DICO[OTU]["multi"].append(int(0))
               DICO[OTU]["Cluster"].append(data["cluster"][cell])
               DICO[OTU]["max_identity"].append(data["blast_perc_identity"][cell])
               DICO[OTU]["min_identity"].append(data["blast_perc_identity"][cell])
               DICO[OTU]["Observation"].append(int(data["observation_sum"][cell]))
               DICO[OTU]["Seq"].append(data["seed_sequence"][cell])
               for i in sites:
                  if pd.isnull(data[i][cell]) == True :
                     DICO[OTU][i].append(int(0))
                  else:
                     if data[i][cell] >=0:
                        DICO[OTU][i].append(int(data[i][cell]))
                     else :
                        DICO[OTU][i].append(int(data[i][cell]))
                        DICO[OTU]["multi"][0]=int(1)
               OTU=OTU+1

########################################################################
### Automatisation de la fusion pour les  données  avec  multi-identité


print("Processing....")

donelist=[]
DICO1={}
OTU=0

multi="AgriB04_before_merge_multi.csv"

csvfile=open(os.path.join(outputdir+multi))
multidata = pd.read_csv(csvfile,delimiter=";",header=0,index_col=None)


for cell in range(len(multidata["blast_taxonomy"])): # Pour chaque nom de taxa
   if multidata["observation_sum"][cell]>0:
      if multidata["blast_taxonomy"][cell] not in donelist : #Si nom jamais rencontré
      
         ###CRÉE LES LISTES
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

         ### REMPLI LES LISTES                
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
            else:
               if data[i][cell] >=0:
                  DICO1[OTU][i].append(int(multidata[i][cell]))
               else :
                  DICO1[OTU][i].append(int(multidata[i][cell]))


         OTU=OTU+1
         donelist.append( multidata["blast_taxonomy"][cell])

      else: # Si déjà rencontré
         liste=0
         while liste < len(DICO1):
            if multidata["blast_taxonomy"][cell] == DICO1[liste]["Taxon_name"][0]: # Je rentre dans le taxon concerner         
                        DICO1[liste]["Nbr_cluster"][0]=DICO1[liste]["Nbr_cluster"][0]+int(1) #ajout 1 au nombre de cluster fusionner
                        DICO1[liste]["Observation"]=int(DICO1[liste]["Observation"]+multidata["observation_sum"][cell])
                        for i in sites:
                           if pd.isnull(multidata[i][cell]) == False:
                              DICO1[liste][i][0]=int(DICO1[liste][i][0]+multidata[i][cell])

                        #Ajout nom cluster
                        DICO1[liste]["Cluster"].append(multidata["cluster"][cell])  
                        DICO1[liste]["Seq"].append(multidata["seed_sequence"][cell]) 
                        break;
                  
                     
            liste=liste+1


#############################
### 2nd check des pseudogenes Enleve les doublons <97% quand ils ont les mêmes occurences que les clusters >97
###VERIFIER CA MARCHE PAS DERNIEREMENT


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
                if prefixsample in k: #CHECK REDONDANCE IN EACH SAMPLES IN DUPLICATE TAXA
                     if DICO[i][k][0] == 0 and DICO[pose[j]][k][0] == 0:
                        nbrsite+=1
                     elif DICO[i][k][0] > 0 and DICO[pose[j]][k][0] > 0:
                        nbrsite+=1
                     else:
                        break
            if nbrsite==len(DICO[0])-8: # IF FULLY REDONDANT, CHOOSE THE TAXA WITH THE MAXIMUM OF IDENTITY (THE SECOND ONE PROBABLY PSEUDOGENE)
               if DICO[i]["max_identity"][0] >  DICO[pose[j]]["max_identity"][0]  :      
                  delposition.append(pose[j])
               if DICO[i]["max_identity"][0] <  DICO[pose[j]]["max_identity"][0]  : 
                  delposition.append(i)    
            nbrsite=0
    it+=1



#####################################################
### Compilation du CSV ##############################
#####################################################


print('Automatisation terminé ')
print('compilation des fichiers')

with open(os.path.join(outputdir+nameprojet+'_AUTO_DATA.csv'), 'w') as f:
   i=0
   for key in DICO[0]:
      if key != 'Cluster' and key != 'Seq':
         [f.write(key)]
         [f.write(';')]
      i=i+1
   i=0
   [f.write('\n')]

   for i in range(len(DICO)):
      if i not in delposition:
         [f.write(DICO[i]["Taxon_name"][0] )]
         [f.write(';')]
         [f.write(str(DICO[i]["Nbr_cluster"][0]))]
         [f.write(';')]
         [f.write(str(DICO[i]["multi"][0]))]
         [f.write(';')]
         if len(DICO[i]["min_identity"]) >=1 : 
            [f.write(str(DICO[i]["max_identity"][0] ))]
            [f.write(';')]
            [f.write(str(DICO[i]["min_identity"][0] ))]
            [f.write(';')]
         else:
            [f.write('')]
            [f.write(';')]
            [f.write('')]
            [f.write(';')]
         [f.write(str(DICO[i]["Observation"]))]       
         for j in sites:
            [f.write(';')]
            [f.write(str(DICO[i][j][0] ))]
         [f.write('\n')]

   for i in range(len(DICO1)):
      [f.write(DICO1[i]["Taxon_name"][0] )]
      [f.write(';')]
      [f.write(str(DICO1[i]["Nbr_cluster"][0]))]
      [f.write(';')]
      [f.write(str(DICO1[i]["multi"][0]))]
      [f.write(';')]
      if len(DICO1[i]["min_identity"]) >=1 : 
        [f.write(str(DICO1[i]["max_identity"][0] ))]
        [f.write(';')]
        [f.write(str(DICO1[i]["min_identity"][0] ))]
        [f.write(';')]
      else:
         [f.write('')]
         [f.write(';')]
         [f.write('')]
         [f.write(';')]
      [f.write(str(DICO1[i]["Observation"]))]       
      for j in sites:
         [f.write(';')]
         [f.write(str(DICO1[i][j][0] ))]
      [f.write('\n')]
f.close()

with open(os.path.join(outputdir+nameprojet+'_Cluster_merged.txt'), 'w') as f:
   for key in DICO:
      if key not in delposition :
         [f.write("> "+ DICO[key]["Taxon_name"][0]+ " : \n")]
         if DICO[key]["Cluster"]:
            for x in DICO[key]["Cluster"]:
               if pd.isna(x)==False:
                  [f.write(x+"\n")]
            [f.write('\n')]

   for key in DICO1:
      [f.write("> "+ DICO1[key]["Taxon_name"][0]+ " : \n")]
      if DICO1[key]["Cluster"]:
         for x in DICO1[key]["Cluster"]:
            if pd.isna(x)==False:
               [f.write(x+"\n")]
         [f.write('\n')]
f.close()

with open(os.path.join(outputdir+nameprojet+'_pseudogene_deleted2.csv'), 'w') as f:
   i=0
   for key in DICO[0]:
      if key != 'Cluster' and key != 'Seq':
         [f.write(key)]
         [f.write(';')]
      i=i+1

   i=0
   [f.write('\n')]

   for i in range(len(DICO)):
      if i in delposition:
         [f.write(DICO[i]["Taxon_name"][0] )]
         [f.write(';')]
         [f.write(str(DICO[i]["Nbr_cluster"][0]))]
         [f.write(';')]
         [f.write(str(DICO[i]["multi"][0]))]
         [f.write(';')]
         if len(DICO[i]["min_identity"]) >=1 : 
            [f.write(str(DICO[i]["max_identity"][0] ))]
            [f.write(';')]
            [f.write(str(DICO[i]["min_identity"][0] ))]
            [f.write(';')]
         else:
            [f.write('')]
            [f.write(';')]
            [f.write('')]
            [f.write(';')]
         [f.write(str(DICO[i]["Observation"]))]       
         for j in sites:
            [f.write(';')]
            [f.write(str(DICO[i][j][0] ))]
         [f.write('\n')]
f.close()


with open(os.path.join(outputdir+nameprojet+'_Cluster_seq_keep.txt'), 'w') as f:
   for key in DICO:
      if key not in delposition :
         if DICO[key]["Cluster"]:
            it=0
            for x in DICO[key]["Cluster"]:
               if pd.isna(x)==False:
                  [f.write("> "+ DICO[key]["Taxon_name"][0]+ "_" + x + " : \n")]
                  [f.write(DICO[key]["Seq"][it]+"\n")]
                  it+=1
               [f.write('\n')]

   for key in DICO1:
         if DICO1[key]["Cluster"]:
            it=0
            for x in DICO1[key]["Cluster"]:
               if pd.isna(x)==False:
                  [f.write("> "+ DICO1[key]["Taxon_name"][0]+ "_" + x + " : \n")]
                  [f.write(DICO1[key]["Seq"][it]+"\n")]
                  it+=1
               [f.write('\n')]
f.close()




print("Programme terminé. Bonne journée")




################
###### END #####
################


