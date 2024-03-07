filter_tfa <- function(physeq, ratio) {
  ## Args:
  ## - physeq: phyloseq class object, otu abundances are extracted from this object
  ## - ratio: the rfa value from Galan et al. (2016)
  ## Returns:
  ## - a physeq object: with abundance data transformed to null under a theshold estimated from ratio/rfa and the sums of the otus. 
  ##                    controls for index contamination following correction 2 from Galan et al. 2016 (Tcc).
    
  # Check orientation, transpose if-needed to make apply work properly.
  if (!taxa_are_rows(physeq)) { newphyseq <- t(as(otu_table(physeq), "matrix")) } else { newphyseq <- as(otu_table(physeq), "matrix") }
  # Filtering
  if (class(ratio) != "numeric" && length(ratio) != 1 && ratio<=1) {
    stop("ratio should be a single numeric value inferior to 1.")
  }
  newphyseq[newphyseq <= (ratio*rowSums(newphyseq))] <- 0
  if (!taxa_are_rows(physeq)) { newphyseq <- t(newphyseq) }
  # Check that original and new dimensions agree. Error if not.
  if( !identical(dim(newphyseq), dim(otu_table(physeq))) ){
    stop("Dimensions of OTU table change after apply-ing function. \n",
         "       Please check both function and table")
  }
  otu_table(physeq) <- otu_table(newphyseq, taxa_are_rows=taxa_are_rows(physeq))
  return(physeq)
}