
"""
Snakemake is a workflow management system written in Python. It facilitates the creation and execution of 
complex data analysis pipelines, particularly in bioinformatics and computational biology. It is associated to a config.yaml file 
as well as a environment.yaml file. Please check snakemake Tutorial for futher explications
(https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html)

Here HAMI_PIPELINE.smk is a snakemake workflow associated to HAMI framework. It is associateds to config.yaml file and HAMI_environement.yaml file.
Please check all these file to understand how it works.


The role of this snakemake pipeline is to process the sequencing data at the output of the FROGS pipeline as part of the HAMI framework.
There are five successive rules:

- All : It represents the final outputs that need to be generated by the workflow. It triggers the execution of the entire workflow.

- Clean_and_Chimeres :  R process for removing chimeric sequences using dada2 packages (Callahan et al., 2016 - Nature Methods)

- Separate_BARCODIN_METABARCODING_data : Python process to discriminate between metabarcoding and barcoding samples

- Filter_METABARCODING_DATA: R process to filter data and eliminate noise and contamination (3 filter steps)

- Pseudogene_Filter_and_reduce_redundancy : Python process to delete pseudogene sequence and reduce data redundancy associates to intraspecific diversty


When the conda environment has been setup and the config file adapted to your projet, the classic way of run the snakemake file is the following :

snakemake --cores 1  -s HAMI/HAMI_PIPELINE.smk   

"""


configfile :"HAMI/config.yaml"

rule all:
    """
    the rule all serves as a target rule that specifies the final output files or targets that should be generated at the end of the workflow.
    The HAMI pipeline is finished when the following three files are created during the last rule (Pseudogene_Filter_and_reduce_redundancy)
    """

    input:
        Occurence_table= expand(os.path.join(config["datadir"],"METABARCODING/final_files","{Name}{DNA}_{target_group}_occurence_file.tsv"),Name=config["Name"],
                        DNA=config["DNA_fragment"], target_group=config["Target_group"]),
        Cluster_merge= expand(os.path.join(config["datadir"],"METABARCODING/final_files","{Name}{DNA}_{target_group}_cluster_merged.txt"),Name=config["Name"],
                        DNA=config["DNA_fragment"], target_group=config["Target_group"]),
        seq_kept=expand(os.path.join(config["datadir"],"METABARCODING/final_files","{Name}{DNA}_{target_group}_cluster_seq_keep.txt"),Name=config["Name"],
                        DNA=config["DNA_fragment"], target_group=config["Target_group"])



rule Clean_and_Chimeres:
    """
    This rule firstly clean the taxonomy (7 ranks: Kingdom,Phylum,Class,Order,Family,Genus,Species) of each hit of the two input files for:
    - the special characters ; 
    -'unidentified' and 'sp.' affiliation which are replaced by 'unknown species'
    - unexpected number of words in all 7 taxonomic ranks; 
    - and more 

    Then it identify chimeric sequences using the isBimeraDenovo(...) function in the dada2 library. By default thread is set at 2, further can be used.

    Finally, it exports three .txt files, a summary of the chimeric sequence deletion step, the list of deleted chimeric sequences and the chimeric sequence-free abundance file.
    """
    input:
        abundance_data=expand(os.path.join(config["datadir"],"{Name}{DNA}_abundance_raw.csv"),Name=config["Name"],DNA=config["DNA_fragment"]),
        multi_data=expand(os.path.join(config["datadir"],"{Name}{DNA}_multi-affiliations.csv"),Name=config["Name"],DNA=config["DNA_fragment"])
    output: 
        summary=expand(os.path.join(config["datadir"],"chimeras","{Name}{DNA}_chimeras_summary.txt"),Name=config["Name"],DNA=config["DNA_fragment"]),
        chimera_list=expand(os.path.join(config["datadir"],"chimeras","{Name}{DNA}_chimeras_list.txt"),Name=config["Name"],DNA=config["DNA_fragment"]),
        cleaned_abundance=expand(os.path.join(config["datadir"],"{Name}{DNA}_cleaned_abundance.txt"),Name=config["Name"],DNA=config["DNA_fragment"])
    
    params: 
        path = lambda wilcards : config["datadir"],
        name = lambda wilcards : config["Name"],
        fragment = lambda wilcards : config["DNA_fragment"],
        threads= lambda wilcards : config["Threads"]

    shell:
        """
        Rscript SCRIPTS/clean_frogs.R {params.path} {input.multi_data} {input.abundance_data} {params.name} {params.fragment} {params.threads}
        """

rule Separate_BARCODIN_METABARCODING_data:
    """
    The aim of this rule is to discriminate Metabarcoding, Barcoding and control samples which were sequenced at the same time in a context of HAMI FRAMEWORK.
    In practical, it separates  Metabarcoding, Barcoding and control datas from the input file which were - previously produce with frogs- by using the sample name nomenclature .
    Your samples thus need to be discrimated on the basis of their name. Please use alphabectic prefix and number for it e.g : CMEY0001 
    Also pipeline is developped to handle duplicate of samples. Discrimination between duplicate will be done using suffix : e.g : CMEY0001A / CMEY0001B
    
    According to the input file, the rule will produce 3 output files : two of them are subset  of input file (Metabarcoding samples and Barcodign samples). The latter is METADATA file associated to metabarcoding datas
    
    Note : If no barcoding sample, the Barcoding subset file will be empty.

    """
    input:
        cleaned_abundance=expand(os.path.join(config["datadir"],"{Name}{DNA}_cleaned_abundance.txt"),Name=config["Name"],DNA=config["DNA_fragment"])
    
    output:
        META_data=expand(os.path.join(config["datadir"],"METABARCODING","{Name}{DNA}_abundance_raw_data_METABARCODING.tsv"),Name=config["Name"],
                        DNA=config["DNA_fragment"]),
        META_METADATA=expand(os.path.join(config["datadir"],"METABARCODING","{Name}{DNA}_Metadata_METABARCODING.csv"),Name=config["Name"],
                        DNA=config["DNA_fragment"]),
        BAR_data=expand(os.path.join(config["datadir"],"BARCODE","{Name}{DNA}_abundance_raw_data_BARCODING.tsv"),Name=config["Name"],
                        DNA=config["DNA_fragment"])


    params: 
        path = lambda wilcards : config["datadir"],
        name = lambda wilcards : config["Name"],
        fragment = lambda wilcards : config["DNA_fragment"],
        prefixMETA = lambda wilcards : config["Samplesprefix"]["Metabarcoding"],
        prefixBARC = lambda wilcards : config["Samplesprefix"]["Barcoding"],
        prefixCONTROL = lambda wilcards : config["Samplesprefix"]["Control"],
        numberID = lambda wilcards : config["Numberofnumber"],
        duplicat = lambda wilcards : config["Samplessuffix"]
    
    shell:
        """
        python3 SCRIPTS/Separate_data.py {input.cleaned_abundance}  {params.path} {params.name} {params.fragment} {params.prefixMETA} {params.prefixBARC} {params.prefixCONTROL} {params.numberID} {params.duplicat}
        """
        

rule Filter_METABARCODING_DATA:
    """
    Filtering of the data:

    - filter1: transform abundance data to null under a theshold estimated from an input value (set with alien_in_samples = no and rfa).
    Controls for index contamination following correction 2 from Galan et al. 2016 (Tfa).

    - filter2: keep positive data only if congruent between duplicates and merge it.

    - filter3: transform abundance data to null under a theshold estimated from negative controls.
    Controls for extraction & pcr contamination following correction 1 from Galan et al. 2016 (Tcc).

    Note : Here, the samples must be arranged by biological units (technical replicates of a same biological unit follow each other)
    if it is not the case, arrange them properly in your sample. 

    """
    input: 
        data=expand(os.path.join(config["datadir"],"METABARCODING","{Name}{DNA}_abundance_raw_data_METABARCODING.tsv"),Name=config["Name"],DNA=config["DNA_fragment"]),
        METADATA=expand(os.path.join(config["datadir"],"METABARCODING","{Name}{DNA}_Metadata_METABARCODING.csv"),Name=config["Name"],DNA=config["DNA_fragment"])
    output: 
        TAX=expand(os.path.join(config["datadir"],"METABARCODING","{Name}{DNA}_tax.filter3.csv"),Name=config["Name"],
                        DNA=config["DNA_fragment"]),
        ABUNDANCE=expand(os.path.join(config["datadir"],"METABARCODING","{Name}{DNA}_abundance.filter3.csv"),Name=config["Name"],
                        DNA=config["DNA_fragment"])

    params:
        scriptpath = lambda wilcards : config["Scriptdir"],
        datapath= lambda wilcards : config["datadir"],
        name = lambda wilcards : config["Name"],
        fragment = lambda wilcards : config["DNA_fragment"]

    shell:
        """
        Rscript SCRIPTS/filter_frogs.R {params.scriptpath} {params.datapath} {input.METADATA} {input.data} {params.name} {params.fragment}
        """


rule Pseudogene_Filter_and_reduce_redundancy:
    """
    Firstly, this rule allows to merge the output (tax_filter3.csv and abundance_filter3.csv) produced by the previous rule and lead to the production of a occurence matrice file as produce by FROGS, but cleaned.

    A pseudogene (NUMPTs) filtering step is then performed as recommended in the literature (see: SONG et al., 2018 - Proceedings of the National Academy of Sciences) on the occurence file.
    It mean that affilaition associated to sequences whose reading frame has shifted -as a result of indel events - or in which stop codons have appeared are discarded. 
    Data associated with pseudogenes is then saved in a first file, while data exempt from pseudogenes is saved in a second file.

    The occurrences file, free of pseudogene, is subsequently processed to reduce data redundancy caused by multiple OTUs associated with the same species, reflecting intraspecific variability.
    This complexity is irrelevant to us and complicates data interpretation. 
    The rule therefore allows to merge the occurrences of OTUs whose taxonomic affiliations are identical and whose match threshold with the reference barcode is similar (up to the arbitrary threshold or lower than the threshold) 

    A second cleaning step is then carried out. Taxa represented twice, i.e. represented by OTUs with values below the threshold and above the reference match threshold, are re-examined. 
    If the 2 OTUs of the same taxon are present in exactly the same samples, the group with the lowest match threshold with the reference barcode is deleted as it is considered to be noise or a residual pseudogene. 
    
    The resultant occurrence file is then saved (final_files/XXXXX_Occurence_file.tsv) and reviewed by a parataxonomist in the context of the HAMI framework.
    """
    input: 
        TAX=expand(os.path.join(config["datadir"],"METABARCODING","{Name}{DNA}_tax.filter3.csv"),Name=config["Name"],
                        DNA=config["DNA_fragment"]),
        ABUNDANCE=expand(os.path.join(config["datadir"],"METABARCODING","{Name}{DNA}_abundance.filter3.csv"),Name=config["Name"],
                        DNA=config["DNA_fragment"])
    
    output:
        Occurence_table= expand(os.path.join(config["datadir"],"METABARCODING/final_files","{Name}{DNA}_{target_group}_occurence_file.tsv"),Name=config["Name"],
                        DNA=config["DNA_fragment"], target_group=config["Target_group"]),
        Cluster_merge= expand(os.path.join(config["datadir"],"METABARCODING/final_files","{Name}{DNA}_{target_group}_cluster_merged.txt"),Name=config["Name"],
                        DNA=config["DNA_fragment"], target_group=config["Target_group"]),
        seq_kept=expand(os.path.join(config["datadir"],"METABARCODING/final_files","{Name}{DNA}_{target_group}_cluster_seq_keep.txt"),Name=config["Name"],
                        DNA=config["DNA_fragment"], target_group=config["Target_group"])


    params:
        prefix = lambda wilcards : config["Samplesprefix"]["Metabarcoding"],
        target_taxon = lambda wilcards : config["Target_group"],
        outputdir= lambda wilcards : config["datadir"],
        name = lambda wilcards : config["Name"],
        fragment = lambda wilcards : config["DNA_fragment"],
        seq_length = lambda wilcards : config["ADN_length"],
        reading_frame = lambda wilcards : config["Reading_frame"],
        codon_stop = lambda wilcards : config["codon_stop"],
        threshold= lambda wilcards : config["Threshold"]

    shell:
        """
        python3 SCRIPTS/pseudogene_and_redudancy.py  {input.ABUNDANCE} {input.TAX} {params.prefix} {params.target_taxon} {params.outputdir} {params.name} {params.fragment} {params.seq_length} {params.reading_frame} {params.threshold} {params.codon_stop} 
        """

