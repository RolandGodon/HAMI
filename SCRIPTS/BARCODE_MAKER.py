#!/usr/bin/env python
# -*-coding: utf-8 -*

"""
Author : Benoit Penel
Date : 10/15/2024
version : 3.0



run command : python3 BARCODE_MAKER.py /path/to/prim1.csv /path/to/prim2.csv    

#### Stage 1: Recovery of majority sequences in each samples for both primers set:

This script parses the two abundance files resulting from the sequencing of BARCODE samples (which were first filtered with BARCODE_pseudogene_filter.py) with the two sets of primers :
BF3 and BR2 primer (CCHGAYATRGCHTTYCCHCG / TCDGGRTGNCCRAARAAYCA; Elbrecht and Leese 2017; Elbrecht et al. 2019)
and LCO1490 and Ill_C_R (GGTCAACAAATCATAAAGATATTGG/ GGIGGRTAIACIGTTCAICC - Shokralla et al. 2015).
Different primers can be used but change need to be done.

Using two independent loops (one for each set of primers), the information associated with the sequences with the highest number of reads are saved
in each sample. This information is stored in a Python dictionary. Note that the majority sequences selected may not be the majority sequences if :

- The majority sequence is taxonomically affiliated to an organisms other than the targeted phylum with a percentage of identity superiror to the threshold.
- The majority sequence has no additional information (no data)

### Step 2: Merging of complementary DNA sequences target with two set of primers :

CR and BB sequences that have been previously identified as the most abundant in a given samples are merged together. If the merging do not produce a full Barcode of 658 bp and with zero mismatch in the overlapping regions, this later is not kept. 


output : For each samples, 1 fasta file containing the most abundante BB and CR sequence is produces as well as a second file with the reconstructed barcode if no issue in the lenght and mismatch in the overlapping region.
"""

###### PARAM YOUR PROJECT :


outputdir1="/home/lbenoit/Documents/Programmes/HAMI/DATA/RUN_AGRIB07/AgriB07_bar/FASTA1"
outputdir2="/home/lbenoit/Documents/Programmes/HAMI/DATA/RUN_AGRIB07/AgriB07_bar/FASTA1/BARCODE"
nameprojet=['FAUN','JHAR','CMEY','m2-JHAR']
prim1="CR" # target the 5'-----> middle of the DNA framgnent
prim2="BB" #target the middle of the DNA fragment -----> 3'
phylum='Insecta'
barcode_length=int(658)
threshold=int(97)

###### PACKAGES :

import sys
import pandas as pd
import os 
from os import getcwd, chdir, mkdir
from Bio.Align import PairwiseAligner

###############
#### BEGIN ####
###############

DICO={}
DICO["SAMPLE"]=[]
DICO["BARCODE"]=[]
DICO["COMMENT"]=[]

### INFORMATION FOR 1 MAJORITARY SEQUENCES
DICO["Prim_SET1"]=[]
DICO["Prim_SET2"]=[]
DICO["nbr_read_Prim_SET1"]=[]
DICO["nbr_read_Prim_SET2"]=[]

### INFORMATION FOR 2ND MAJORITARY SEQUENCES
DICO["Prim_SET1_2"]=[]
DICO["Prim_SET2_2"]=[]
DICO["nbr_read_Prim_SET1_2"]=[]
DICO["nbr_read_Prim_SET2_2"]=[]



###OPEN FILES
csvfile=open(sys.argv[1])
BARCODE_Prim_SET1 = pd.read_csv(csvfile,delimiter=";",header=0)

csvfile=open(sys.argv[2])
BARCODE_Prim_SET2 = pd.read_csv(csvfile,delimiter=";",header=0)


#####################################
####### SEARCH FOR MAJORITY SEQUENCE

print("Research of the majority sequence in each samples for ", prim1, " primer")
######
## Prim_SET1

MAX=0
it=0
position=-1
oldposition=-1
oldMAX=0


seq=pd.Series(BARCODE_Prim_SET1.seed_sequence)
name=pd.Series(BARCODE_Prim_SET1.blast_taxonomy)
identity=pd.Series(BARCODE_Prim_SET1.blast_perc_identity)

for i in BARCODE_Prim_SET1 :
    for prefix in nameprojet :
        if prefix in i :                                                              #Recognition of sample columns
            DICO["SAMPLE"].append(i)
            site=pd.Series(BARCODE_Prim_SET1[i])
            for j in site:
                if j > 0 and j > MAX :                            #If a sequence has a higher read abundance than a previously recorded sequence
                    if phylum not in str(name.iloc[it]) :         #If the taxonomic affiliation associated with this sequence is other than the expected phylum 
                        try :
                            num=float(identity.iloc[it])          #And the affiation have a pourcentage of identity > to the threshold, this majority sequence is discarded.
                            if num >=threshold:
                                it+=1
                                continue
                            else:                               #If not, we converse this new majority sequence
                            
                                if position >0 :                #Keep the second most abundant sequence
                                    oldposition=position
                                    oldMAX=MAX                                                                     
                                    position=it
                                    MAX=j
                                else:                           #preserve the most abundant sequence
                                    position=it
                                    MAX=j

                    
                        except ValueError:
                        
                            if "data" in str(identity.iloc[it]): #If this sequence has no additional information, it is discarded.
                                it+=1
                                continue
                            else:                               #If not, we converse this new majority sequence
                                if position >0 :                #Keep the second most abundant sequence
                                    oldposition=position
                                    oldMAX=MAX                                                                     
                                    position=it
                                    MAX=j
                                else:                           #preserve the most abundant sequence
                                    position=it
                                    MAX=j

                    else:                                #If the taxonomic affiliation associated with this sequence is other the expected phylum 
                        if position >0 :                 #Keep the second most abundant sequence
                            oldposition=position
                            oldMAX=MAX                                                                     
                            position=it
                            MAX=j
                        else:                           #Keep the most abundant sequence
                            position=it
                            MAX=j

                elif j>oldMAX and j<MAX:                      #If the number of reads associated with a sequence is lower than for the most important sequence but higher than for the second most important sequence 
                    if phylum not in str(name.iloc[it]) :     # and the taxonomic affiliation associated with this sequence is other than the expected target taxon
                        try :
                            num=float(identity.iloc[it])       #And the affiation have a pourcentage of identity > to the threshold, this majority sequence is discarded.
                            if num >=threshold:
                                it+=1
                                continue
                            else:                                #If not, we converse this new secondary majority sequence
                                oldposition=it
                                oldMAX=j                                                                    

                        except ValueError:
                        
                            if "data" in str(identity.iloc[it]):         #If this sequence has no additional information, it is discarded.
                                it+=1
                                continue
                            else:                               #If not, we converse this new secondary majority sequence
                                oldposition=it
                                oldMAX=j

                    else:
                        oldposition=it
                        oldMAX=j


                it+=1
        
            #SAVED POSITION OF THE MAJORITARIES SEQUENCES 
            if position != -1:  
                DICO["Prim_SET1"].append(seq.iloc[position])
            else:
                DICO["Prim_SET1"].append("")
            if oldposition != -1:
                DICO["Prim_SET1_2"].append(seq.iloc[oldposition])
            else:
                DICO["Prim_SET1_2"].append("")

            DICO["nbr_read_Prim_SET1"].append(MAX)
            DICO["nbr_read_Prim_SET1_2"].append(oldMAX)


            it=0
            MAX=0
            position=-1
            oldposition=-1
            oldMAX=0


print("Research of the majority sequence in each samples for ", prim2, " primer")
###############
### Prim_SET2

#Similar process is applied for the second set of sequence 
#Produce with the second primer set 


seq=pd.Series(BARCODE_Prim_SET2.seed_sequence)
name=pd.Series(BARCODE_Prim_SET2.blast_taxonomy)
identity=pd.Series(BARCODE_Prim_SET2.blast_perc_identity)


for i in BARCODE_Prim_SET2 :
    for prefix in nameprojet:
        if prefix in i and i in DICO["SAMPLE"]:
            site=pd.Series(BARCODE_Prim_SET2[i])
            for j in site:
                if j > 0 and j > MAX :
                    if phylum not in str(name.iloc[it]) :
                        try :
                            num=float(identity.iloc[it])
                            if num >=threshold:
                                it+=1
                                continue
                            else:
                                if position>0:
                                    oldposition=position
                                    oldMAX=MAX
                                    position=it
                                    MAX=j 
                                else:
                                    position=it
                                    MAX=j

                    
                        except ValueError:
                            if "data" in str(identity.iloc[it]):
                                it+=1
                                continue
                            else:
                                if position>0:
                                    oldposition=position
                                    oldMAX=MAX
                                    position=it
                                    MAX=j 
                                else:
                                    position=it
                                    MAX=j
                    else:
                        if position>0:
                            oldposition=position
                            oldMAX=MAX
                            position=it
                            MAX=j 
                        else:
                            position=it
                            MAX=j

                elif j>oldMAX and j<MAX:                                               
                    if phylum not in str(name.iloc[it]) :                    
                        try :
                            num=float(identity.iloc[it])                             
                            if num >=threshold:
                                it+=1
                                continue
                            else:                               
                                oldposition=it
                                oldMAX=j                                                                    

                        except ValueError:
                        
                            if "data" in str(identity.iloc[it]):         
                                it+=1
                                continue
                            else:                              
                                oldposition=it
                                oldMAX=j

                    else:
                        oldposition=it
                        oldMAX=j

                it+=1

            if position != -1:
                DICO["Prim_SET2"].append(seq.iloc[position])
            else:
                DICO["Prim_SET2"].append("")


            if oldposition != -1:
                DICO["Prim_SET2_2"].append(seq.iloc[oldposition])
            else:
                DICO["Prim_SET2_2"].append("")

            DICO["nbr_read_Prim_SET2"].append(MAX)
            DICO["nbr_read_Prim_SET2_2"].append(oldMAX)


            it=0
            MAX=0
            position=-1
            oldposition=-1
            oldMAX=0


##################################
###### CREATING FASTA FILES


print("Write fasta files")
if not os.path.exists(outputdir1): #Check that there is a fasta directory in which the majority of sequences are stored.
    os.mkdir(outputdir1)
os.chdir(outputdir1)

for i in range(len(DICO["SAMPLE"])):
    fasta=open(DICO["SAMPLE"][i]+".fasta","w")
    fasta.write(">"+ DICO["SAMPLE"][i]+"_"+prim1+"_"+str(DICO["nbr_read_Prim_SET1"][i])+"\n"+DICO["Prim_SET1"][i]+"\n")
    fasta.write(">"+ DICO["SAMPLE"][i]+"_"+prim2+"_"+str(DICO["nbr_read_Prim_SET2"][i])+"\n"+DICO["Prim_SET2"][i]+"\n") 
    fasta.close()


#############################################################
####Fragments alignment obtained with the two set of primers and barcode creation

print("Sequence alignment")
### First try :
for i in range(len(DICO["SAMPLE"])):
    try:

        ###################################
        ###Aligndment

        seq1 = DICO["Prim_SET1"][i]
        seq2 = DICO["Prim_SET2"][i]
    
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -50                      #Big penalty for introducing a gap in the sequence alignment. 
        aligner.extend_gap_score = -0.1

        aligned_seqs = aligner.align(seq1, seq2)
        
    
        ########################
        ###Alignment merger

        fused_seq = ""
        COMMENT=""
        aligned_seq1, aligned_seq2 = aligned_seqs[0]     

        for j in range(len(aligned_seq1)):                  #For each nucleotide of the alignment 
            if aligned_seq1[j] != "-":                      #As long as the alignment produced does not have a gap, its nucleotides will be recovered.
                if aligned_seq2[j]=='-': 
                    fused_seq += aligned_seq1[j]
                else:                                       
                    if aligned_seq1[j] != aligned_seq2[j]:  #If there is a mismatch 
                        fused_seq = ""                      #BARCODE FAIL
                        COMMENT="MISMATCH"                  #REASON OF THE FAIL
                        break;
                    else:                                   #If no mismatch detected
                        fused_seq += aligned_seq1[j]        #Saved alignement sequence

            else:                                            
                if aligned_seq2[j]!='-' :                   #If no mismatch detected as well as no gap
                    fused_seq += aligned_seq2[j]            #Saved alignement sequence
                else:                                       #If gap detected 
                    fused_seq = ""                          #BARCODE FAIL
                    COMMENT="GAP"                           #REASON OF THE FAIL
                    break;

        if len(fused_seq) == barcode_length :                           #IF fused_seq reach the expected size = Full Barcode recovery
            DICO["BARCODE"].append(fused_seq)
            DICO["COMMENT"].append("none")
        else:   
            if len(fused_seq)<barcode_length:                                            #Otherwise =  Barcode is not functional -> Fail
                DICO["BARCODE"].append("MISSING_BARCODE")               
                DICO["COMMENT"].append(COMMENT)
            elif len(fused_seq)>barcode_length:
                DICO["BARCODE"].append(fused_seq)
                DICO["COMMENT"].append("Too long")

    except:
        DICO["BARCODE"].append("MISSING_BARCODE")           #One of the DNA fragments were not extracted -> MISSING BARCODE
        DICO["COMMENT"].append("MISSING FRAGMENT")





### Second try : 
print("Second alignment try")
# Use the second majority sequence
DICO2={}
DICO2["SAMPLE"]=[]
DICO2["BARCODE2-1"]=[]
DICO2["BARCODE2-2"]=[]
DICO2["BARCODE2-3"]=[]
DICO2["COMMENT2-1"]=[]
DICO2["COMMENT2-2"]=[]
DICO2["COMMENT2-3"]=[]

for i in range(len(DICO["SAMPLE"])) :
    if DICO["COMMENT"][i]!="none" and DICO["COMMENT"][i]!="MISSING FRAGMENT": #If first trial succed, skip the sample
        DICO2["SAMPLE"].append(DICO["SAMPLE"][i])
        # 2-1
        try:

            ###################################
            ###Alignement ADN of primer set 1 majoritary and the second majoritary of primer set 2

            seq1 = DICO["Prim_SET1"][i]
            seq2 = DICO["Prim_SET2_2"][i]
    
            aligner = PairwiseAligner()
            aligner.mode = 'global'
            aligner.match_score = 2
            aligner.mismatch_score = -1
            aligner.open_gap_score = -50                      
            aligner.extend_gap_score = -0.1

            aligned_seqs = aligner.align(seq1, seq2)

            ########################
            ###Alignment merger

            fused_seq = ""
            COMMENT=""
            aligned_seq1, aligned_seq2 = aligned_seqs[0]        
            for j in range(len(aligned_seq1)):
                if aligned_seq1[j] != "-":                      
                    if aligned_seq2[j]=='-': 
                        fused_seq += aligned_seq1[j]
                    else:               
                        if aligned_seq1[j] != aligned_seq2[j]:  
                            fused_seq = ""
                            COMMENT="MISMATCH"
                            break;
                        else:                                   
                            fused_seq += aligned_seq1[j]

                else:                                           
                    if aligned_seq2[j]!='-' :                   
                        fused_seq += aligned_seq2[j]
                    else:                                       
                        fused_seq = ""
                        COMMENT="GAP"
                        break;

            if len(fused_seq) == barcode_length:                           
                DICO2["BARCODE2-1"].append(fused_seq)
                DICO2["COMMENT2-1"].append("none")
            else: 
                if len(fused_seq)<barcode_length:                                              
                    DICO2["BARCODE2-1"].append("MISSING_BARCODE")
                    DICO2["COMMENT2-1"].append(COMMENT)
                elif len(fused_seq)>barcode_length:
                    DICO2["BARCODE2-1"].append(fused_seq)
                    DICO2["COMMENT2-1"].append("Too long")

        except:
            DICO2["BARCODE2-1"].append("MISSING_BARCODE")           
            DICO2["COMMENT2-1"].append("MISSING FRAGMENT")
    
    # 2-2
        try:

            ###################################
            ###Alignement ADN of primer set 2 majoritary and the second majoritary of primer set 1

            seq1 = DICO["Prim_SET1_2"][i]
            seq2 = DICO["Prim_SET2"][i]

            aligner = PairwiseAligner()
            aligner.mode = 'global'
            aligner.match_score = 2
            aligner.mismatch_score = -1
            aligner.open_gap_score = -50                       
            aligner.extend_gap_score = -0.1

            aligned_seqs = aligner.align(seq1, seq2)

            ########################
            ###Alignment merger

            fused_seq = ""
            COMMENT=""
            aligned_seq1, aligned_seq2 = aligned_seqs[0]        
            for j in range(len(aligned_seq1)):
                if aligned_seq1[j] != "-":                      
                    if aligned_seq2[j]=='-': 
                        fused_seq += aligned_seq1[j]
                    else:               
                        if aligned_seq1[j] != aligned_seq2[j]:  
                            fused_seq = ""
                            COMMENT="MISMATCH"
                            break;
                        else:                                   
                            fused_seq += aligned_seq1[j]

                else:                                            
                    if aligned_seq2[j]!='-' :                   
                        fused_seq += aligned_seq2[j]
                    else:                                       
                        fused_seq = ""
                        COMMENT="GAP"
                        break;

            if len(fused_seq) ==barcode_length :                           
                DICO2["BARCODE2-2"].append(fused_seq)
                DICO2["COMMENT2-2"].append("none")
            else:                                               
                if len(fused_seq)<barcode_length:                                              
                    DICO2["BARCODE2-2"].append("MISSING_BARCODE")
                    DICO2["COMMENT2-2"].append(COMMENT)
                elif len(fused_seq)>barcode_length:
                    DICO2["BARCODE2-2"].append(fused_seq)
                    DICO2["COMMENT2-2"].append("Too long")

        except:
            DICO2["BARCODE2-2"].append("MISSING_BARCODE")          
            DICO2["COMMENT2-2"].append("MISSING FRAGMENT")
      
     # 2-3
        try:

            ###################################
            ###Alignement of the second majoritary DNA sequence for both primer set 

            seq1 = DICO["Prim_SET1_2"][i]
            seq2 = DICO["Prim_SET2_2"][i]

            aligner = PairwiseAligner()
            aligner.mode = 'global'
            aligner.match_score = 2
            aligner.mismatch_score = -1
            aligner.open_gap_score = -50                       
            aligner.extend_gap_score = -0.1

            aligned_seqs = aligner.align(seq1, seq2)

            ########################
            ###Alignment merger

            fused_seq = ""
            COMMENT=""
            aligned_seq1, aligned_seq2 = aligned_seqs[0]        
            for j in range(len(aligned_seq1)):
                if aligned_seq1[j] != "-":                      
                    if aligned_seq2[j]=='-': 
                        fused_seq += aligned_seq1[j]
                    else:               
                        if aligned_seq1[j] != aligned_seq2[j]:  
                            fused_seq = ""
                            COMMENT="MISMATCH"
                            break;
                        else:                                   
                            fused_seq += aligned_seq1[j]

                else:                                            
                    if aligned_seq2[j]!='-' :                   
                        fused_seq += aligned_seq2[j]
                    else:                                       
                        fused_seq = ""
                        COMMENT="GAP"
                        break;

            if len(fused_seq) ==barcode_length :                           
                DICO2["BARCODE2-3"].append(fused_seq)
                DICO2["COMMENT2-3"].append("none")
            else:                                               
                if len(fused_seq)<barcode_length:                                              
                    DICO2["BARCODE2-3"].append("MISSING_BARCODE")
                    DICO2["COMMENT2-3"].append(COMMENT)
                elif len(fused_seq)>barcode_length:
                    DICO2["BARCODE2-3"].append(fused_seq)
                    DICO2["COMMENT2-3"].append("Too long")

        except:
            DICO2["BARCODE2-3"].append("MISSING_BARCODE")          
            DICO2["COMMENT2-3"].append("MISSING FRAGMENT") 

##########################################
###### CREATING FASTA BARCODE FILES

print("Write barcode files")
if not os.path.exists(outputdir2):  #Check that there is a barcode directory in which the majority of sequences are stored.
    os.mkdir(outputdir2)
os.chdir(outputdir2)

texte=open("1.Log_BARCODING_MISSING.txt","w") ## FILE WHO SAVE MISSING BARCODE SAMPLE NAMES
fullbarcode=open("1.FULL_BARCODE.fasta","w") ## FILE WHO SAVE ALL BARCODE IN ONE FILE

for i in range(len(DICO["SAMPLE"])):
    if DICO["COMMENT"][i] =="none":                                   #IF FIRST BARCODE OK
        fasta=open(DICO["SAMPLE"][i]+".fasta","w")                          #CREATE FASTA FILE TO SAVE BARCODE FOR EACH SAMPLE  
        fasta.write(">"+ DICO["SAMPLE"][i]+"\n"+DICO["BARCODE"][i]+"\n")
        fasta.close()
        fullbarcode.write(">"+ DICO["SAMPLE"][i]+"\n"+DICO["BARCODE"][i]+"\n")   #SAVE BARCODE IN A FULL SAMPLES FILES BARCODE  

    else:                                                                           # IF FIRST BARCODE NOT OK, CHECK THE SECOND BARCODE PRODUCE
        
        for j in range (len(DICO2["SAMPLE"])):
            if DICO["SAMPLE"][i]==DICO2["SAMPLE"][j]:
                if  DICO2["BARCODE2-1"][j]!= "MISSING_BARCODE" or  DICO2["BARCODE2-2"][j]!= "MISSING_BARCODE" or DICO2["BARCODE2-3"][j]!= "MISSING_BARCODE" :     ## IF at least 1 NEW BARCODE ARE NOW AVAILABLE:
                    
                    outputdir3=outputdir2+"/2nd_alignement"
                    if not os.path.exists(outputdir3):
                        os.mkdir(outputdir3)
                    os.chdir(outputdir3)

                    fasta=open(DICO2["SAMPLE"][j]+"_2nd_alignement.fasta","w")   

                    if DICO2["BARCODE2-1"][j]!= "MISSING_BARCODE" :
                        fasta.write(">"+ DICO2["SAMPLE"][j]+"__1"+"\n"+DICO2["BARCODE2-1"][j]+"\n")
                    elif DICO2["BARCODE2-2"][j]!= "MISSING_BARCODE" :
                        fasta.write(">"+ DICO2["SAMPLE"][j]+"__2"+"\n"+DICO2["BARCODE2-2"][j]+"\n")
                    elif DICO2["BARCODE2-3"][j]!= "MISSING_BARCODE" :
                        fasta.write(">"+ DICO2["SAMPLE"][j]+"__3"+"\n"+DICO2["BARCODE2-3"][j]+"\n")
                    fasta.close()
                    
                    os.chdir(outputdir2)

                    if DICO2["BARCODE2-1"][j]!= "MISSING_BARCODE" :
                        fullbarcode.write(">"+ DICO2["SAMPLE"][j]+"___1"+"\n"+DICO2["BARCODE2-1"][j]+"\n")       # SAVE BARCODES IN A FULL SAMPLES FILES BARCODE
                    elif DICO2["BARCODE2-2"][j]!= "MISSING_BARCODE" :
                        fullbarcode.write(">"+ DICO2["SAMPLE"][j]+"___2"+"\n"+DICO2["BARCODE2-2"][j]+"\n")
                    elif DICO2["BARCODE2-3"][j]!= "MISSING_BARCODE" :
                        fullbarcode.write(">"+ DICO2["SAMPLE"][j]+"___3"+"\n"+DICO2["BARCODE2-3"][j]+"\n")


                else :               # SAVE NAME OF SAMPLES WITH A MISSING BARCODE
                    os.chdir(outputdir2)
                    texte.write(">"+DICO["SAMPLE"][i]+" : "+ DICO["COMMENT"][i]+"\n")


texte.close()
fullbarcode.close()




##########################################
###### CREATING BARCODE EXCEL FILES
import csv

#Ouverture d'un fichier CSV en mode écriture
with open('1.BARCODE.csv', 'w', newline='') as fichier_csv:
    writer = csv.writer(fichier_csv)

    writer.writerow(DICO.keys())
    for ligne in zip(*DICO.values()):
        writer.writerow(ligne)





###########
## END ####
###########     


###############
####  END  ####
###############
