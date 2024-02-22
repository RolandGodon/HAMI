filter_tcc <- function(physeq, phyneg) {
  ## Requirements:
  ## - a variable name 'control' in the corresponding sample_data of phyloseq
  ## Args:
  ## - physeq: phyloseq class object, otu abundances are extracted from this object
  ## - group: a variable name in the corresponding sample_data of phyloseq which correspond to the levels for controls
  ## - control.level: a single character string matching a level factor of the variable name 'control' in the corresponding sample_data of phyloseq
  ##          should correspond to negative controls (if filter == Tcc) 
  ## Returns:
  ## - a physeq object: with abundance data transformed to null under a theshold estimated from negative controls. 
  ##                    controls for extraction & pcr contamination following correction 1 from Galan et al. 2016 (Tcc).
  
 # Check orientation, transpose if-needed to make apply work properly.
  if (!taxa_are_rows(physeq)) { 
    newphyseq <- t(as(otu_table(physeq), "matrix")) 
    neg<-t(as(otu_table(phyneg), "matrix")) 
  } 
  else {
       newphyseq <- as(otu_table(physeq), "matrix") 
       neg <-as(otu_table(phyneg), "matrix") 
  }
  # Filtering
  newphyseq[newphyseq <= (rowMaxs(neg))] <- 0
  if (!taxa_are_rows(physeq)) { newphyseq <- t(newphyseq) }
  # Check that original and new dimensions agree. Error if not.
  if( !identical(dim(newphyseq), dim(otu_table(physeq))) ){
    stop("Dimensions of OTU table change after applying function. \n",
         "       Please check both function and table")
  }
  otu_table(physeq) <- otu_table(newphyseq, taxa_are_rows=taxa_are_rows(physeq))
  return(physeq)
}