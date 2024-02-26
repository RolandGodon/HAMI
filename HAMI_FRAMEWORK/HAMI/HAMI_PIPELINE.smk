configfile :"HAMI/config.yaml"


rule all:
    input:
        Occurence_table= expand(os.path.join(config["datadir"],"METABARCODING/final_files","{Name}{DNA}_{target_group}_occurence_file.tsv"),Name=config["Name"],
                        DNA=config["DNA_fragment"], target_group=config["Target_group"]),
        Cluster_merge= expand(os.path.join(config["datadir"],"METABARCODING/final_files","{Name}{DNA}_{target_group}_cluster_merged.txt"),Name=config["Name"],
                        DNA=config["DNA_fragment"], target_group=config["Target_group"]),
        seq_kept=expand(os.path.join(config["datadir"],"METABARCODING/final_files","{Name}{DNA}_{target_group}_cluster_seq_keep.txt"),Name=config["Name"],
                        DNA=config["DNA_fragment"], target_group=config["Target_group"])



rule Clean_and_Chimeres:
    """
    Clean the taxonomy (7 ranks: Kingdom,Phylum,Class,Order,Family,Genus,Species) of each hit of the two files for:
    - the special characters ; 
    -'mitochondria'  and 'chloroplast' affiliation;
    -'unidentified', 'sp.' and 'bacterium' affiliation which are replaced by 'unknown species'
    - unexpected number of words in all 7 taxonomic ranks; 
    and more 

    Identify chimeric sequences using the isBimeraDenovo(...) function in the dada2 library.

    Export an abundance table tsv file cleaned for taxonomy and filtered for chimeras.
    """
    input: 
        abundance_data=expand(os.path.join(config["datadir"],"{Name}{DNA}_abundance_raw.csv"),Name=config["Name"],DNA=config["DNA_fragment"]),
        multi_data=expand(os.path.join(config["datadir"],"{Name}{DNA}_multi-affiliations.csv"),Name=config["Name"],DNA=config["DNA_fragment"])
    output:
        summary=expand(os.path.join(config["datadir"],"{Name}{DNA}_chimeras_summary.txt"),Name=config["Name"],DNA=config["DNA_fragment"]),
        chimera_list=expand(os.path.join(config["datadir"],"{Name}{DNA}_chimeras_list.txt"),Name=config["Name"],DNA=config["DNA_fragment"]),
        cleaned_abundance=expand(os.path.join(config["datadir"],"{Name}{DNA}_cleaned_abundance.txt"),Name=config["Name"],DNA=config["DNA_fragment"])
    
    params: 
        path = lambda wilcards : config["datadir"],
        name = lambda wilcards : config["Name"],
        fragment = lambda wilcards : config["DNA_fragment"]

    shell:
        """
        Rscript SCRIPTS/clean_frogs.R {params.path} {input.multi_data} {input.abundance_data} {params.name} {params.fragment}
        """

rule separate_BARCODIN_METABARCODING_data:
    """
    The aim of this rule is to discriminate Metabarcoding, Barcoding and control samples which were sequenced at the same time in a context of HAMI FRAMEWORK.
    In practical, we separate  Metabarcoding, Barcoding and control  data from the input file which were - previously produce with frogs- by using the sample name nomenclature .

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

    - filter1: transform abundance data to null under a theshold estimated from alien controls (set with alien_in_samples = yes) or from an input value (set with alien_in_samples = no and rfa).
    Controls for index contamination following correction 2 from Galan et al. 2016 (Tfa).

    - filter2: keep positive data only if congruent between replicates (2 or 3) and merge replicates.

    - filter3: transform abundance data to null under a theshold estimated from negative controls.
    Controls for extraction & pcr contamination following correction 1 from Galan et al. 2016 (Tcc).

    Note : Here, the samples must be arranged by biological units (technical replicates of a same biological unit follow each other)
    if it is not the case, arrange them properly in your sample. 
    Note 2 : Explore the R script if you want to change any details (search "#to modify")
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
        codon_stop = lambda wilcards : config["codon_stop"]

    shell:
        """
        python3 SCRIPTS/Merge_after_filtre.py  {input.ABUNDANCE} {input.TAX} {params.prefix} {params.target_taxon} {params.outputdir} {params.name} {params.fragment} {params.seq_length} {params.reading_frame} {params.codon_stop}
        """

