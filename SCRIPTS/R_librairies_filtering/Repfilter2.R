filter_rep <- function(physeq, nb_rep, nb_con) {
  ## Requirements:
  ## - the samples muste be arranged by biological units (technical replicates of a same biological unit follow each other)
  ##   if it is not the case in your biom:
  ##   arrange them properly in your sample_data then,   
  ##   use the "Import each component separately" which sorts the otu_table using the order in sample_data
  ## - physeq: phyloseq class object, otu abundances are extracted from this object
  ## - group: a variable name in the corresponding sample_data of phyloseq which correspond to the level of technical replicates
  ## - the number of technical replicates should be either 2 or 3 (nb_rep)
  ## - the number of technical replicates to be congruent can be either 2 (for nb_rep = 2 or 3) or 3 (for nb_rep = 3) (nb_con)
  ## Returns:
  ## - a physeq object: with abundance data transformed to null for incongruent replicates. 
    # Functions
    presenceabsence<-function(x){
      x[x>0] <- 1
      return(x)
    }
    
    expandRows <- function(dataset, count) {
      if (length(count) == 1) {
        dataset[rep(rownames(dataset), each = count), ]
      } else {
        if (length(count) != nrow(dataset)) {
          stop("Expand vector does not match number of rows in data.frame")
        }
        dataset[rep(rownames(dataset), count), ]
      }
    }
    
    # Produce a TRUE/FALSE table that check replicate congruence
    pres<-transform_sample_counts(physeq, presenceabsence)
    
    if (nb_rep == 2) {
      rep1<-subset_samples(pres, rep %in% c(1))
      rep2<-subset_samples(pres, rep %in% c(2))
      if (taxa_are_rows(physeq)) { 
        equal<-t(otu_table(rep1) == otu_table(rep2))
        equal_expanded<-expandRows (equal,nb_rep)
        equal_expanded<-t(equal_expanded)
      }
    else {
       equal<-otu_table(rep1) == otu_table(rep2)
       equal_expanded<-expandRows (equal,nb_rep) 
    }
    }
    
    if (nb_rep == 3) {
      rep1<-subset_samples(pres, rep %in% c(1))
      rep2<-subset_samples(pres, rep %in% c(2))
      rep3<-subset_samples(pres, rep %in% c(3))
      if (taxa_are_rows(physeq)) { 
        equal12<-t(otu_table(rep1) == otu_table(rep2) & otu_table(rep1) == 1)
        equal13<-t(otu_table(rep1) == otu_table(rep3) & otu_table(rep1) == 1)
        equal23<-t(otu_table(rep2) == otu_table(rep3) & otu_table(rep2) == 1)
        if (nb_con == 3) {
        equal<-(equal12 * equal13) 
        equal_expanded<-expandRows (equal,nb_rep)
        equal_expanded<-t(equal_expanded) }
        if (nb_con == 2) {
          equal<-(equal12 | equal13 | equal23) 
          equal_expanded<-expandRows (equal,nb_rep)
          equal_expanded<-t(equal_expanded) }
      }
      else {
        equal12<-otu_table(rep1) == otu_table(rep2) 
        equal13<-otu_table(rep1) == otu_table(rep3)
        equal23<-otu_table(rep2) == otu_table(rep3)
        if (nb_con == 3) {
          equal<-(equal12 * equal13) 
          equal_expanded<-expandRows (equal,nb_rep)
          equal_expanded<-t(equal_expanded) }
        if (nb_con == 2) {
          equal<-(equal12 | equal13 | equal23) 
          equal_expanded<-expandRows (equal,nb_rep)
          equal_expanded<-t(equal_expanded) }
      }
    }
  # Check parameters and physeq orientation, transpose if-needed to make apply work properly.
    if (!taxa_are_rows(physeq)) { newphyseq <- t(as(otu_table(physeq), "matrix")) } else { newphyseq <- as(otu_table(physeq), "matrix") }
  # Filtering
    colnames(equal_expanded)<-colnames(newphyseq)
    newphyseq<-ifelse(equal_expanded,newphyseq,0)
    if (!taxa_are_rows(physeq)) { newphyseq <- t(newphyseq) }
  # Check that original and new dimensions agree. Error if not.
    if( !identical(dim(newphyseq), dim(otu_table(physeq))) ){
    stop("Dimensions of OTU table change after apply-ing function. \n",
         "       Please check both function and table")
    }
  otu_table(physeq) <- otu_table(newphyseq, taxa_are_rows=taxa_are_rows(physeq))
  return(physeq)
}